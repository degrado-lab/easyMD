import argparse
import json
import os

parser = argparse.ArgumentParser(description='Run a batch of simulations.')
parser.add_argument('experiment_dir', type=str, help='The path to the batch directory.')

args = parser.parse_args()

parameters_file = args.experiment_dir + '/parameters.json'

with open(parameters_file, 'r') as f:
    parameters = json.load(f)

'''
If we haven't yet started a simulation, we'll need to run our minimization + equilibration.
If we've already finished that, we'll start production.
If production has already started, continue the simulation until the specified # of steps is met.
'''

### CHECK IF MINIMIZATION AND EQUIL HAVE BEEN RUN ###
output_dir = parameters['paths']['outputs']['output_dir']
output_eq_file = parameters['paths']['outputs']['output_eq_file']
# does the eq_file exist?
minim_eq_complete = os.path.exists(os.path.join(output_dir, output_eq_file))

from easyMD.utils import run_sim_local

### BEGIN MINIMIZATION AND EQUILIBRATION ###
print('parameters_file:', parameters_file)
run_sim_local(parameters_file, pre_production=not minim_eq_complete)

