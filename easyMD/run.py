def pre_production_wrapper(experiment_dir, params_file, custom_func=None):
    # Load json file of prep_params:
    import json
    with open(params_file, 'r') as f:
        params = json.load(f)

    from pathlib import Path
    input_processed_dir = Path(experiment_dir) / 'input' / 'processed'
    input_processed_dir.mkdir(parents=True, exist_ok=True)
    simulation_dir = Path(experiment_dir) / 'simulation'
    simulation_dir.mkdir(parents=True, exist_ok=True)

    if 'SIMULATION' in params:
        params = params['SIMULATION']

    print('Beginning energy minimization and equilibration...')

    # Load custom prep function if provided:
    if custom_func is not None:
        from easyMD.utils import import_function_from_path
        pre_prod_function = import_function_from_path(custom_func, 'pre_production')
        print('Using custom pre_production function:', custom_func)
    # Use default prep function:
    else:
        pre_prod_function = pre_production

    pre_prod_function(experiment_dir, params)

    print('Finished Pre-Production!')
    # Save output
    # Save prepped system to input/prepped dir

def pre_production(experiment_dir, params):
    #def run_sim_local(parameters_file, continue_sim=True, extend_sim=False, skip_pre_production=False, force_no_gpu=False):
    '''
    Run a simulation using OpenMM.

    Inputs:
    - experiment_dir (str): The experiment directory.
    - parameters_file (str): The path to the parameters file, which contains:
        - minimization_length (int): Duration of minimization, in nanoseconds.
        - equilibration_length (int): Duration of equilibration, in nanoseconds.
        - reporting_interval (int): Reporting interval for simulation, in nanoseconds.
        - step_size (float): Step size for simulation, in femtoseconds.
        - temperature (float): Temperature for simulation, in Kelvin.
        - pressure (float): Pressure for simulation, in atmospheres.
        - forcefield_files (list of str): List of forcefield files to use.
    '''
    from openmm.app import PDBFile, ForceField, Simulation, DCDReporter, StateDataReporter, CheckpointReporter, PME, HBonds
    from openmm import LangevinIntegrator, MonteCarloBarostat
    from openmm import unit as openmm_unit
    from openff.toolkit import Molecule

    ############################################
    ############ GET SIM PARAMETERS ############
    ############################################

    minimization_length      = params['minimization_length_NS']    * openmm_unit.nanoseconds
    equilibration_length     = params['equilibration_length_NS']    * openmm_unit.nanoseconds

    step_size               = params['step_size_FS']           * openmm_unit.femtoseconds          # in femtoseconds
    temperature             = params['temperature_K']         * openmm_unit.kelvin                # in Kelvin
    pressure                = params['pressure_ATM']            * openmm_unit.atmosphere            # in atmospheres
    reporting_interval      = params['save_frame_interval_NS']  * openmm_unit.nanoseconds           # in nanoseconds 
    forcefield_files        = params['forcefield']                                              # list of forcefield files

    ############################################
    ############### SETUP SYSTEM ###############
    ############################################
    # Load the processed PDB:
    from pathlib import Path
    processed_dir = Path(experiment_dir) / 'input' / 'processed'
    # Get the file labeled _solvated.cif or _solvated.pdb
    pdb_path = None
    for file in processed_dir.iterdir():
        if file.name[-13:] in ['_solvated.cif', '_solvated.pdb']:
            pdb_path = file
            break
    if pdb_path is None:
        raise FileNotFoundError('No solvated structure found in the processed directory. Please run the prep function first. (Must be named *_solvated.cif or *_solvated.pdb)')
    
    from openmm.app import PDBFile
    pdb = PDBFile(str(pdb_path))

    # Set up the Forcefield:
    forcefield = ForceField(*forcefield_files)
    
    # Are we using ligands? If so, add their templates to the forcefield.
    from openmmforcefields.generators import GAFFTemplateGenerator
    ligand_dir = Path(experiment_dir) / 'input' / 'processed' / 'ligands'
    for ligand_file in ligand_dir.iterdir():
        if Path(ligand_file).suffix == '.sdf':
            print('Adding ligand {} to the forcefield.'.format(ligand_file))
            ligand = Molecule.from_file(ligand_file)
            forcefield.registerTemplateGenerator( GAFFTemplateGenerator(molecules=ligand).generator)

    #Create system:
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*openmm_unit.nanometer, constraints=HBonds)
    integrator = LangevinIntegrator(temperature, 1/openmm_unit.picosecond, step_size.in_units_of(openmm_unit.picoseconds))
    barostat = MonteCarloBarostat(pressure, temperature, 25)
    system.addForce(barostat)
    simulation = Simulation(pdb.topology, system, integrator)
    
    ############################################
    ############ RUN MINIMIZATION ##############
    ############################################
    print('Minimizing energy...')
    
    simulation.context.setPositions(pdb.positions)
    import numpy as np
    minimization_steps = int(np.round(minimization_length / step_size))
    simulation.minimizeEnergy(maxIterations=minimization_steps)

    # write out the minimized pdb
    print('Writing out the minimized pdb...')
    output_minim_file = Path(experiment_dir) / 'simulation' / (pdb_path.stem + '_minimized.pdb')
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), str(output_minim_file))

    ############################################
    ############ RUN EQUILIBRATION #############
    ############################################
    print('Running equilibration...')
    import sys
    output_eq_file = Path(experiment_dir) / 'simulation' / (pdb_path.stem + '_equilibrated.dcd')
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.reporters.append(DCDReporter(str(output_eq_file), int(np.round(reporting_interval / step_size))))
    simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))
    equilibration_steps = int(np.round(equilibration_length / step_size))
    simulation.step(equilibration_steps)

def production_wrapper(experiment_dir, params_file, custom_func=None):
    # Load json file of prep_params:
    import json
    with open(params_file, 'r') as f:
        params = json.load(f)

    from pathlib import Path
    input_processed_dir = Path(experiment_dir) / 'input' / 'processed'
    input_processed_dir.mkdir(parents=True, exist_ok=True)
    simulation_dir = Path(experiment_dir) / 'simulation'
    simulation_dir.mkdir(parents=True, exist_ok=True)

    if 'SIMULATION' in params:
        params = params['SIMULATION']

    print('Beginning production simulation...')

    # Load custom prep function if provided:
    if custom_func is not None:
        from easyMD.utils import import_function_from_path
        prod_function = import_function_from_path(custom_func, 'production')
        print('Using custom production function:', custom_func)
    # Use default prep function:
    else:
        prod_function = production

    prod_function(experiment_dir, params)

    print('Finished Production!')

def production(experiment_dir, params):
    '''
    Run a simulation using OpenMM.

    Inputs:
    - experiment_dir (str): The experiment directory.
    - parameters_file (str): The path to the parameters file, which contains:
        - simulation_length (int): Duration of simulation, in nanoseconds.
        - reporting_interval (int): Reporting interval for simulation, in nanoseconds.
        - checkpoint_interval (int): Checkpoint interval for simulation, in nanoseconds.
        - step_size (float): Step size for simulation, in femtoseconds.
        - temperature (float): Temperature for simulation, in Kelvin.
        - pressure (float): Pressure for simulation, in atmospheres.
        - forcefield_files (list of str): List of forcefield files to use.
    '''
    from openmm.app import PDBFile, ForceField, Simulation, DCDReporter, StateDataReporter, CheckpointReporter, PME, HBonds
    from openmm import LangevinIntegrator, MonteCarloBarostat
    from openmm import unit as openmm_unit
    from openff.toolkit import Molecule

    ############################################
    ############ GET SIM PARAMETERS ############
    ############################################

    simulation_length       = params['simulation_length_NS']    * openmm_unit.nanoseconds
    reporting_interval      = params['save_frame_interval_NS']  * openmm_unit.nanoseconds
    checkpoint_interval     = params['checkpoint_interval_NS'] * openmm_unit.nanoseconds

    step_size               = params['step_size_FS']           * openmm_unit.femtoseconds          # in femtoseconds
    temperature             = params['temperature_K']         * openmm_unit.kelvin                # in Kelvin
    pressure                = params['pressure_ATM']            * openmm_unit.atmosphere            # in atmospheres
    forcefield_files        = params['forcefield']                                        # list of forcefield files

    #double check: If checkpoint_interval is larger than production_steps, set it to production_steps
    if checkpoint_interval > simulation_length: 
        checkpoint_interval = simulation_length

    ############################################
    ############### SETUP SYSTEM ###############
    ############################################
    # Load the processed PDB:
    from pathlib import Path
    simulation_dir = Path(experiment_dir) / 'simulation'
    # Get the file labeled _solvated.cif or _solvated.pdb
    pdb_path = None
    for file in simulation_dir.iterdir():
        if file.name[-14:] in ['_minimized.cif', '_minimized.pdb']:
            pdb_path = file
            break
    if pdb_path is None:
        raise FileNotFoundError('No minimized structure found in the simulation directory. Please run the pre_production function first. (Must be named *_minim.cif or *_minim.pdb)')
    
    from openmm.app import PDBFile
    pdb = PDBFile(str(pdb_path))

    # Set up the Forcefield:
    forcefield = ForceField(*forcefield_files)
    
    # Are we using ligands? If so, add their templates to the forcefield.
    from openmmforcefields.generators import GAFFTemplateGenerator
    ligand_dir = Path(experiment_dir) / 'input' / 'processed' / 'ligands'
    for ligand_file in ligand_dir.iterdir():
        if Path(ligand_file).suffix == '.sdf':
            print('Adding ligand {} to the forcefield.'.format(ligand_file))
            ligand = Molecule.from_file(ligand_file)
            forcefield.registerTemplateGenerator( GAFFTemplateGenerator(molecules=ligand).generator)

    #Create system:
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*openmm_unit.nanometer, constraints=HBonds)
    integrator = LangevinIntegrator(temperature, 1/openmm_unit.picosecond, step_size.in_units_of(openmm_unit.picoseconds))
    simulation = Simulation(pdb.topology, system, integrator)

    # Set positions:
    simulation.context.setPositions(pdb.positions)
        
    # ############################################
    # ############ RUN PRODUCTION ################
    # ############################################
    print('Running production...')
    simulation.reporters.clear()
    output_prod_file = Path(experiment_dir) / 'simulation' / (pdb_path.stem + '_production.dcd')
    output_checkpoint_file = Path(experiment_dir) / 'simulation' / (pdb_path.stem + '_checkpoint.chk')
    import sys
    import numpy as np
    simulation.reporters.append(DCDReporter(  str(output_prod_file), int(np.round(reporting_interval / step_size)), append=False)) #this creates a new file only if we just did minimization and equilibration
    simulation.reporters.append(StateDataReporter(  sys.stdout, int(np.round(reporting_interval / step_size)), step=True, potentialEnergy=True, temperature=True))
    simulation.reporters.append(CheckpointReporter( str(output_checkpoint_file), int(np.round(checkpoint_interval / step_size))))
    print('Running simulation...')
    
    simulation.step(int(np.round(simulation_length / step_size)))