import sys
from pathlib import Path
from openmm import *
from openmm.app import *
from openmm import unit as openmm_unit

input_json_file = sys.argv[1]
### Are we continuing a previous simulation?
continue_from_previous_sim = True if sys.argv[2] == 'True' else False
continue_sim_steps = int(sys.argv[3]) #500000 = 1ns. If we are continuing a previous simulation, this overrides the production_steps variable.

############################################
############ GET SIM PARAMETERS ############
############################################

#add .. to the path so we can import from the parent directory
sys.path.append(str(Path(__file__).resolve().parents[1]))
from utils.io import read_sim_parameters

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
only_production      =  read_sim_parameters(input_json_file)

#double check: If checkpoint_interval is larger than production_steps, set it to production_steps
if checkpoint_interval > production_steps: 
    checkpoint_interval = production_steps

############################################
############### SETUP SYSTEM ###############
############################################

# Here, we load the processed pdb and set up the system:
pdb_path = processed_dir / Path(input_pdb_file).name
pdb = PDBFile(str(pdb_path))
forcefield = ForceField(*forcefield_files)
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


from easyMD.utils import run_sim 









