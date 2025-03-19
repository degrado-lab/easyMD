from openmm import app, unit
from easyMD.app import prepare_simulation, run_simulation, process_trajectory
from easyMD.io import write_pdb_from_simulation
import logging

if __name__ == "__main__":

    # Set up logger:
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('easyMD').setLevel(logging.INFO)
    
    top_file = 'output1.pdb'
    traj_file = 'output1.dcd'

    simulation = prepare_simulation("downloads/7H6N_A.cif", 
                                output_pdb='output1.pdb', 
                                fix=True)
    
    try:
        # Print info to stdout, and save to DCD file:
        run_simulation(simulation, traj_file, duration=1*unit.nanoseconds, output_frequency=2000*unit.femtoseconds)
    except KeyboardInterrupt:
        print("Interrupted")
    except Exception as e:
        raise e
    
    # Then, rewrap + center the trajectory:
    aligned_traj_path = 'output_aligned1.dcd'
    process_trajectory(top_file, traj_file, aligned_traj_path)