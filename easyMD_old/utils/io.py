import json #package up as, read/write parameters
from pathlib import Path

def read_sim_parameters(input_json_file):
    '''
    Read in the input parameters from a json file.

    Parameters
    ----------
    input_json_file : str
        The path to the input json file.
        
    Returns
    -------
    minimization_steps : int
        The number of minimization steps.
    equilibration_steps : int
        The number of equilibration steps.
    production_steps : int
        The number of production steps.
    reporting_interval : int
        The reporting interval for the simulation (in steps).
    checkpoint_interval: int
        The checkpoint interval for the simulation (in steps).
    input_pdb_file : str
        The path to the input pdb file.
    sim_dir : Path
        The path to the simulation directory.
    processed_dir : Path
        The path to the processed directory.
    experiment_name : str
        The name of the experiment.
    forcefield_files : list
        The list of forcefield files.
    step_size : int
        The step size for the simulation (in fs).
    temperature : int
        The temperature for the simulation (in K).
    only_production : bool
        Whether to only run the production step.
    '''
    with open(input_json_file) as f:
        input_data = json.load(f)
    #extract the input parameters
    minimization_steps = input_data['minimization_steps']
    equilibration_steps = input_data['equilibration_steps']
    production_steps = input_data['production_steps']
    reporting_interval = input_data['reporting_interval']
    checkpoint_interval = input_data['checkpoint_interval']
    input_pdb_file = input_data['input_pdb_file']
    sim_dir = Path(input_data['sim_dir'])
    processed_dir = Path(input_data['processed_dir'])
    experiment_name = input_data['experiment_name']
    forcefield_files = input_data['forcefield_files']
    step_size = input_data['step_size']
    temperature = input_data['temperature']
    only_production = input_data['only_production']

    return minimization_steps,  \
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
            only_production
            
def write_sim_parameters(input_json_file, input_dict):
    '''
    Write the input parameters to a json file.

    Parameters
    ----------
    input_json_file : str
        The path to the input json file.
    input_dict: dict
        The input parameters.
        Must contain the following keys:

        minimization_steps : int
            The number of minimization steps.
        equilibration_steps : int
            The number of equilibration steps.
        production_steps : int
            The number of production steps.
        reporting_interval : int
            The reporting interval for the simulation (in steps).
        checkpoint_interval: int
            The checkpoint interval for the simulation (in steps).
        input_pdb_file : str
            The path to the input pdb file.
        sim_dir : Path
            The path to the simulation directory.
        processed_dir : Path
            The path to the processed directory.
        experiment_name : str
            The name of the experiment.
        forcefield_files : list
            The list of forcefield files.
        step_size : int
            The step size for the simulation (in fs).
        temperature : int
            The temperature for the simulation (in K).
        only_production : bool
            Whether to only run the production step.
    '''
    
    #verify that the input_dict contains the required keys
    required_keys = ['minimization_steps', 'equilibration_steps', 'production_steps', 'reporting_interval', 'checkpoint_interval', 'input_pdb_file', 'sim_dir', 'processed_dir', 'experiment_name', 'forcefield_files', 'step_size', 'temperature', 'only_production']
    for key in required_keys:
        if key not in input_dict:
            raise ValueError(f'input_dict must contain the key {key}')
        
    with open(input_json_file, 'w') as f:
        json.dump(input_dict, f, indent=4)

def setup_experiment_dir(experiments_dir, experiment_name):
    '''
    Create the experiments directory if it does not exist.

    Parameters
    ----------
    experiments_dir : Path
        The path to the experiments directory.
    experiment_name : str
        The name of the experiment.

    Returns
    -------
    current_experiment_dir : Path
        The path to the current experiment directory.
    inputs_dir : Path
        The path to the inputs directory.
    raw_dir : Path
        The path to the raw directory.
    processed_dir : Path
        The path to the processed directory.
    sim_dir : Path
        The path to the simulation directory.
    '''
    experiments_dir = Path('experiments')
    current_experiment_dir = experiments_dir / experiment_name
    inputs_dir = current_experiment_dir / 'inputs'
    raw_dir = inputs_dir / 'raw'
    processed_dir = inputs_dir / 'processed'
    sim_dir = current_experiment_dir / 'simulations'

    #Make the directories if they don't exist
    dir_list = [experiments_dir, current_experiment_dir, inputs_dir, raw_dir, processed_dir, sim_dir]
    for dir in dir_list:
        if not dir.exists():
            #make directory recursively:
            dir.mkdir(parents=True)

    return current_experiment_dir, inputs_dir, raw_dir, processed_dir, sim_dir

def create_dirs(dir_list):
    '''
    Create the directories if they do not exist.
    Input:
    dir_list: list of Path objects
    
    Returns:
    None
    '''
    #Make the directories if they don't exist
    for dir in dir_list:
        if not dir.exists():
            #make directory recursively:
            dir.mkdir(parents=True)

###############################
### DOWNLOAD EXTERNAL FILES ###
###############################

import requests
def download_pdb(pdb_id):
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(f'{pdb_id}.pdb', 'w') as file:
            file.write(response.text)
        print(f'Successfully downloaded {pdb_id}.pdb')
    else:
        print(f'Failed to download {pdb_id}.pdb. Please check the PDB ID.')


