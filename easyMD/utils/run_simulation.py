from .io import read_sim_parameters
from pathlib import Path
from openmm.app import PDBFile, ForceField, Simulation, DCDReporter, StateDataReporter, CheckpointReporter, PME, HBonds
from openmm import LangevinIntegrator
from openmm import unit as openmm_unit
from openff.toolkit import Molecule
import sys

def run_sim_local(sim_dir, continue_from_previous_sim=False, continue_sim_steps=None):
    '''
    Run a simulation using OpenMM.

    Parameters
    ----------
    sim_dir : str
        The path to the simulation directory.
    continue_from_previous_sim : bool, optional
        Whether to continue from a previous simulation. Default is False.
    continue_sim_steps : int, optional
        The number of steps to continue from. Required if continue_from_previous_sim is True.
    '''

    ############################################
    ############ GET SIM PARAMETERS ############
    ############################################

    ### Are we continuing a previous simulation?
    if continue_from_previous_sim:
        #make sure they've specified the number of steps to continue
        if continue_sim_steps is None:
            raise ValueError("When continuing a simulation, you must manually specify the number of steps to run using the 'continue_sim_steps' parameter.")

    sim_dir = Path(sim_dir)
    sim_parameters_file = sim_dir / 'simulation_dict.json'

    minimization_steps, \
    equilibration_steps,\
    production_steps,   \
    reporting_interval, \
    checkpoint_interval,\
    input_pdb_file,     \
    sim_dir,            \
    processed_dir,      \
    experiment_name,    \
    forcefield_files,   \
    step_size,          \
    temperature,        \
    only_production      =  read_sim_parameters(sim_parameters_file)

    #double check: If checkpoint_interval is larger than production_steps, set it to production_steps
    if checkpoint_interval > production_steps: 
        checkpoint_interval = production_steps

    ############################################
    ############### SETUP SYSTEM ###############
    ############################################

    # Load the processed PDB:
    pdb_path = processed_dir / Path(input_pdb_file).name
    pdb = PDBFile(str(pdb_path))

    # Set up the Forcefield:
    forcefield = ForceField(*forcefield_files)
    
    # Are we using ligands? If so, add their templates to the forcefield.
    from openmmforcefields.generators import GAFFTemplateGenerator
    ligand_dir = processed_dir / 'ligands'
    if ligand_dir.exists():
        #iterate over files:
        for ligand_file in ligand_dir.iterdir():
            if ligand_file.suffix == '.sdf':
                ligand = Molecule.from_file(ligand_file)
                forcefield.registerTemplateGenerator( GAFFTemplateGenerator(molecules=ligand).generator)

    # Create simulation:
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*openmm_unit.nanometer, constraints=HBonds)
    integrator = LangevinIntegrator(temperature, 1/openmm_unit.picosecond, step_size)
    simulation = Simulation(pdb.topology, system, integrator)

    # Are we continuing a simulation from a previous run?
    if continue_from_previous_sim:
        # Load the checkpoint file
        checkpoint_file = sim_dir / (str(Path(input_pdb_file).stem) + '_prod.chk')
        if not checkpoint_file.exists():
            raise ValueError(f'Tried to continue a simulation, but the checkpoint file does not exist. Please check the path.\n\
                            File Not Found: {checkpoint_file}')
        simulation.loadCheckpoint(str(checkpoint_file))
        print('Continuing simulation from checkpoint file...')
        # make sure to skip minimization and equilibration:
        only_production = True
    else:
        simulation.context.setPositions(pdb.positions)
        

    if not only_production:
        ############################################
        ############ RUN MINIMIZATION ##############
        ############################################
        print('Minimizing energy...')
        simulation.minimizeEnergy(maxIterations=minimization_steps)

        # write out the minimized pdb
        print('Writing out the minimized pdb...')
        PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), str(sim_dir / (str(Path(input_pdb_file).stem) + '_minim.pdb')))

        ############################################
        ############ RUN EQUILIBRATION #############
        ############################################
        print('Running equilibration...')
        simulation.context.setVelocitiesToTemperature(temperature)
        simulation.reporters.append(DCDReporter(sim_dir / (str(Path(input_pdb_file).stem) + '_equil.dcd'), reporting_interval))
        simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))
        simulation.step(equilibration_steps)

    ############################################
    ############ RUN PRODUCTION ################
    ############################################
    print('Running production...')
    simulation.reporters.clear()
    simulation.reporters.append(DCDReporter(        sim_dir / (str(Path(input_pdb_file).stem) + '_prod.dcd'), reporting_interval, append=continue_from_previous_sim))
    simulation.reporters.append(StateDataReporter(  sys.stdout, reporting_interval, step=True, potentialEnergy=True, temperature=True, append=continue_from_previous_sim))
    simulation.reporters.append(CheckpointReporter( str(sim_dir / (str(Path(input_pdb_file).stem) + '_prod.chk')), checkpoint_interval))
    simulation.step(production_steps)

    print('Done!')


def run_sim_wynton(sim_dir, continue_from_previous_sim=False, continue_sim_steps=None, mem='4G', max_runtime='2:00:00'):
    '''
    Queue a simulation on UCSF Wynton using OpenMM.
    This method queues the job using the SGE scheduler.

    Parameters
    ----------
    sim_dir : str
        The path to the simulation directory.
    continue_from_previous_sim : bool, optional
        Whether to continue from a previous simulation. Default is False.
    continue_sim_steps : int, optional
        The number of steps to continue from. Required if continue_from_previous_sim is True.
    '''

    ############################################
    ############ GET SIM PARAMETERS ############
    ############################################

    ### Are we continuing a previous simulation?
    if continue_from_previous_sim:
        #make sure they've specified the number of steps to continue
        if continue_sim_steps is None:
            raise ValueError("When continuing a simulation, you must manually specify the number of steps to run using the 'continue_sim_steps' parameter.")

    sim_dir = Path(sim_dir)
    sim_parameters_file = sim_dir / 'simulation_dict.json'

    _,_,_,_,_,_,        \
    sim_dir,            \
    _,                  \
    experiment_name,    \
    _,_,_,_             =  read_sim_parameters(sim_parameters_file)


    ###################################################
    ################## RUN SCRIPT #####################
    ###################################################
    bash_script = f'''
#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -N sim_{experiment_name.replace('/', '_')} # set the name of the job
#$ -j y               # STDERR and STDOUT should be joined
#$ -o {sim_dir}/       # set the destination for the STDOUT file
#$ -l mem_free={mem}     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=2G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt={max_runtime}    # job requires up to this many hours of runtime
#$ -r y               # if job crashes, it should be restarted
#$ -q gpu.q           # use the gpu queue

## load the required modules
module load CBI
module load miniconda3/4.12.0-py39
module load Sali
module load cuda/10.0.130

## print start time:
date

## Run the simulation
conda activate easyMD
python3 -c \
"from easyMD.utils import run_sim_local; \
run_sim_local('{sim_dir}', continue_from_previous_sim={continue_from_previous_sim}, continue_sim_steps={continue_sim_steps})"

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                        # e.g. "did my job exceed its memory request?"

## print end time:
date
    ''' 

    with open(sim_dir / 'run_simulation.sh', 'w') as file:
        file.write(bash_script)

    #Finally, let's queue the job:
    import subprocess
    command = ['qsub', '-cwd', str(sim_dir / 'run_simulation.sh')] #subprocess takes commands as a list of strings
    #result = subprocess.run(command, shell=True, capture_output=True, text=True)
    #output = result.stdout.strip()
    result = subprocess.run(command, capture_output=True)
    output = result.stdout.strip()
    print(output.decode()) #turn from bytes to string and print
    print("Use the command 'qstat' in the console to check the status of your job.")
