from .io import read_sim_parameters
from pathlib import Path
from openmm.app import PDBFile, ForceField, Simulation, DCDReporter, StateDataReporter, CheckpointReporter, PME, HBonds
from openmm import LangevinIntegrator, MonteCarloBarostat
from openmm import unit as openmm_unit
from openff.toolkit import Molecule
import sys
import os
import json

def run_sim_local(parameters_file, continue_sim=True, extend_sim=False, skip_pre_production=False, force_no_gpu=False):
    '''
    Run a simulation using OpenMM.

    Parameters
    ----------
    parameters_file : str
        The path to the parameters file.
    continue_sim : bool, optional
        This will automatically continue a simulation from the last checkpoint file. This includes running all preproduction (unless specified otherwise) 
        and production steps. This enables several jobs to be queued with the same jobscript, and they will continue after the previous. Default is True.
    extend_sim: bool, optional
        This will extend the simulation from the last checkpoint file, beyond what is specified in the parameters file. 
        This will only run the production steps. Default is False. Not currently implemented
    skip_pre_production : bool, optional
        Whether to skip minimization and equilibration steps. Default is True.
    force_no_gpu : bool, optional
        Whether to force the simulation to run on the CPU. Default is False, and will cancel if GPU is not available.
    '''

    ############################################
    ############ GET SIM PARAMETERS ############
    ############################################

    parameters = json.load(open(parameters_file))

    minimization_steps      = parameters['simulation']['minimization_steps']
    equilibration_steps     = parameters['simulation']['equilibration_steps']
    production_steps        = parameters['simulation']['production_steps']
    reporting_interval      = parameters['simulation']['reporting_interval']
    checkpoint_interval     = parameters['simulation']['checkpoint_interval']
    step_size               = parameters['simulation']['step_size']
    temperature             = parameters['simulation']['temperature']
    only_production         = parameters['simulation']['only_production']
    forcefield_files        = parameters['simulation']['forcefield_files']
    
    experiment_dir          = parameters['paths']['experiment_dir']
    input_dir               = os.path.join(experiment_dir, parameters['paths']['inputs']['input_dir'])
    input_pdb_file          = os.path.join(input_dir, parameters['paths']['inputs']['input_pdb_file'])
    input_ligand_files      = [os.path.join(input_dir, ligand_file) for ligand_file in parameters['paths']['inputs']['input_ligand_files']]
    sim_dir                 = os.path.join(experiment_dir, parameters['paths']['outputs']['output_dir'])
    output_minim_file       = os.path.join(sim_dir, parameters['paths']['outputs']['output_minim_file'])
    output_eq_file          = os.path.join(sim_dir, parameters['paths']['outputs']['output_eq_file'])
    output_prod_file        = os.path.join(sim_dir, parameters['paths']['outputs']['output_prod_file'])
    output_checkpoint_file  = os.path.join(sim_dir, parameters['paths']['outputs']['output_checkpoint_file'])


    #double check: If checkpoint_interval is larger than production_steps, set it to production_steps
    if checkpoint_interval > production_steps: 
        checkpoint_interval = production_steps

    ###########################################
    ############### VERIFY GPU ################
    ###########################################
    import subprocess
    if not force_no_gpu:
        result = subprocess.run(['nvidia-smi'], capture_output=True)
        assert result.returncode == 0, 'No GPU detected. \n\
            Have you setup apptainer to use a GPU with the easyMD_setup.sh script?'

    ############################################
    ############### SETUP SYSTEM ###############
    ############################################
    # Load the processed PDB:
    pdb_path = input_pdb_file
    pdb = PDBFile(str(pdb_path))

    # Set up the Forcefield:
    forcefield = ForceField(*forcefield_files)
    
    # Are we using ligands? If so, add their templates to the forcefield.
    from openmmforcefields.generators import GAFFTemplateGenerator

    for ligand_file in input_ligand_files:
        if Path(ligand_file).suffix == '.sdf':
            print('Adding ligand {} to the forcefield.'.format(ligand_file))
            # verify file exists:
            if not Path(ligand_file).exists():
                raise ValueError(f'Tried to add a ligand to the forcefield, but the file does not exist. Please check the path.\n\
                                File Not Found: {ligand_file}')
            ligand = Molecule.from_file(ligand_file)
            forcefield.registerTemplateGenerator( GAFFTemplateGenerator(molecules=ligand).generator)

    #Create system:
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*openmm_unit.nanometer, constraints=HBonds)
    integrator = LangevinIntegrator(temperature, 1/openmm_unit.picosecond, step_size)
    barostat = MonteCarloBarostat(1*openmm_unit.atmosphere, temperature, 25)
    system.addForce(barostat)
    simulation = Simulation(pdb.topology, system, integrator)
    
    ############################################
    ############ CONTINUING SIM? ###############
    ############################################
    if continue_sim and Path(output_checkpoint_file).exists():
        # Load the checkpoint file
        print('Continuing simulation from production checkpoint file...')
        simulation.loadCheckpoint(output_checkpoint_file)

        # make sure to skip minimization and equilibration:
        if not skip_pre_production:
            print('Skipping minimization and equilibration steps...')
        skip_pre_production = True

        # How many steps have we already taken?
        current_step = simulation.currentStep

        if current_step >= production_steps:
            print(f'Simulation has already completed {current_step} steps. Exiting...')
            return
    
        print(f'Continuing simulation from step {current_step}...')
        production_steps -= current_step
    else:
        simulation.context.setPositions(pdb.positions)

    if not skip_pre_production:
        ############################################
        ############ RUN MINIMIZATION ##############
        ############################################
        print('Minimizing energy...')
        
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy(maxIterations=minimization_steps)

        # write out the minimized pdb
        print('Writing out the minimized pdb...')
        PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), output_minim_file)

        ############################################
        ############ RUN EQUILIBRATION #############
        ############################################
        print('Running equilibration...')
        simulation.context.setVelocitiesToTemperature(temperature)
        simulation.reporters.append(DCDReporter(output_eq_file, reporting_interval))
        simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))
        simulation.step(equilibration_steps)
        

    ############################################
    ############ RUN PRODUCTION ################
    ############################################
    print('Running production...')
    #Turn off barostat:
    barostat.setFrequency(0)
    simulation.reporters.clear()
    simulation.reporters.append(DCDReporter(  output_prod_file, reporting_interval, append=True))
    simulation.reporters.append(StateDataReporter(  sys.stdout, reporting_interval, step=True, potentialEnergy=True, temperature=True, append=True))
    simulation.reporters.append(CheckpointReporter( output_checkpoint_file, checkpoint_interval))
    simulation.step(production_steps)

    print('Done!')


def run_sim_wynton(apptainer_path, parameters_file, mem='4G', max_runtime='2:00:00', repeat=1):
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
    parameters = json.load(open(parameters_file))
    
    experiment_dir          = parameters['paths']['experiment_dir']
    sim_dir                 = os.path.join(experiment_dir, parameters['paths']['outputs']['output_dir'])

    experiment_name         = parameters['metadata']['experiment_name']

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

## make sure we only use the assigned GPU:
export CUDA_VISIBLE_DEVICES=$SGE_GPU
echo $CUDA_VISIBLE_DEVICES

## Run the simulation
/{apptainer_path} sim {parameters_file}

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                        # e.g. "did my job exceed its memory request?"

## print end time:
date
    ''' 

    sim_script = Path(sim_dir) / 'run_simulation.sh'
    with open(sim_script, 'w') as file:
        file.write(bash_script)

    #Finally, let's queue the job:
    import subprocess

    previous_jid=None
    for i in range(repeat):
        #if i > 0: continue_from_previous_sim = True
        if previous_jid is not None:
            command = ['qsub', '-cwd', '-hold_jid', previous_jid, str(sim_script)]
        else:
            command = ['qsub', '-cwd', str(sim_script)] #subprocess takes commands as a list of strings
        result = subprocess.run(command, capture_output=True)
        if result.returncode != 0:
            print(result.stderr)
            break
        output = result.stdout.strip()
        print(output.decode()) #turn from bytes to string and print

        previous_jid = output.decode().split()[2]
