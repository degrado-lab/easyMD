import sys
from pathlib import Path
from openmm import *
from openmm.app import *
from openmm.unit import *

############################################
############# INPUT PARAMETERS ##############
############################################

#the input parameters are stored in a json file, of the form:
# {'minimization_steps': minimization_steps,
# 'equilibration_steps': equilibration_steps,
# 'production_steps': production_steps,
# 'reporting_interval': reporting_interval,
# 'input_pdb_file': input_pdb_file,
# 'sim_dir': str(sim_dir),
# 'processed_dir': str(processed_dir),
# 'experiment_name': experiment_name,
# 'forcefield_files': forcefield_files,
# 'step_size': step_size,
# 'temperature': temperature,
# 'only_production': False}}

input_json_file = sys.argv[1]

#read in the json file
import json
with open(input_json_file) as f:
    input_data = json.load(f)

#extract the input parameters
minimization_steps = input_data['minimization_steps']
equilibration_steps = input_data['equilibration_steps']
production_steps = input_data['production_steps']
reporting_interval = input_data['reporting_interval']
input_pdb_file = input_data['input_pdb_file']
sim_dir = Path(input_data['sim_dir'])
processed_dir = Path(input_data['processed_dir'])
experiment_name = input_data['experiment_name']
forcefield_files = input_data['forcefield_files']
step_size = input_data['step_size']
temperature = input_data['temperature']
only_production = input_data['only_production']

    


############################################
############## CREATE SYSTEM ###############
############################################
#create a new simulation object from the processed pdb
pdb_path = processed_dir / Path(input_pdb_file).name
pdb = PDBFile(str(pdb_path))
forcefield = ForceField(*forcefield_files)
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(temperature, 1/picosecond, step_size)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)


############################################
############ RUN MINIMIZATION ##############
############################################
if not only_production:
    print('Minimizing energy...')
    simulation.minimizeEnergy(maxIterations=minimization_steps)

    # write out the minimized pdb
    print('Writing out the minimized pdb...')
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), str(sim_dir / (str(Path(input_pdb_file).stem) + '_minim.pdb')))

############################################
############ RUN EQUILIBRATION #############
############################################
if not only_production:
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
simulation.reporters.append(DCDReporter(sim_dir / (str(Path(input_pdb_file).stem) + '_prod.dcd'), reporting_interval))
simulation.reporters.append(StateDataReporter(sys.stdout, reporting_interval, step=True, potentialEnergy=True, temperature=True))
simulation.step(production_steps)

print('Done!')