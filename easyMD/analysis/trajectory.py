# my_md_package/simulation/trajectory.py
import logging
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np
from tqdm import tqdm

logger = logging.getLogger(__name__)

def reduce_frames(top_path, traj_path, reduced_traj_path, num_frames):
    """
    Reduce the number of frames in a trajectory file by taking every nth frame.
    Parameters:
    - top_path: path to the topology file (PDB)
    - traj_path: path to the trajectory file (DCD)
    - reduced_traj_path: path to save the reduced trajectory file (DCD)
    - num_frames: number of frames to keep
    """

    u = mda.Universe(top_path, traj_path)
    print("Number of frames in the trajectory:", len(u.trajectory))

    all_atoms = u.select_atoms('all')

    step_size = len(u.trajectory)//num_frames
    all_atoms.write(reduced_traj_path, frames=u.trajectory[::step_size])

def rewrap_positions(positions, center, box):
    """
    Re-wrap atom positions within the periodic boundaries using modulo operation.
    Parameters:
    - positions: 3D numpy array of atom positions
    - box: 3D numpy array of box dimensions
    - selection: selection string for the group of atoms to center on (default: 'protein')
    Returns:
    - rewrapped_positions: 3D numpy array of re-wrapped atom positions
    """
    #center on protein:
    positions -= center
    positions += box/2
    positions = np.mod(positions, box)
    return positions

def create_rewrapped_universe(u, selection='protein'):
    """
    Creates a new universe with rewrapped positions for all frames.
    Parameters:
    - u: MDAnalysis Universe object
    Returns:
    - new_u: MDAnalysis Universe object with rewrapped positions
    """
    new_positions = []  # To store rewrapped positions for all frames

    # Iterate through each frame and apply re-wrapping
    for ts in tqdm(u.trajectory):
        box = ts.dimensions[:3]  # Get current box dimensions
        center = u.select_atoms(selection).center_of_geometry()
        wrapped_positions = rewrap_positions(ts.positions, center, box)
        new_positions.append(wrapped_positions)

    # Convert list of positions to a 3D numpy array
    new_positions_array = np.array(new_positions)

    # Create a new Universe with the rewrapped positions using MemoryReader
    new_u = mda.Universe(u.filename, u.trajectory.filename)
    new_u.load_new(new_positions_array, format=MemoryReader)

    return new_u

def rewrap_trajectory(top_path, traj_path, rewrapped_traj_path):
    """
    Rewrap the trajectory within the periodic boundaries.
    Parameters:
    - top_path: path to the topology file (PDB)
    - traj_path: path to the trajectory file (DCD)
    - rewrapped_traj_path: path to save the rewrapped trajectory file (DCD)
    """
    try:
        # Add debug information
        logger.info(f"Loading trajectory from {traj_path}")
        logger.info(f"Using topology from {top_path}")
        
        # Try loading with topology_format explicitly set for OpenMM PDB
        u = mda.Universe(
            top_path, 
            traj_path,
            topology_format='PDB',
            format='DCD',
            permissive=True
        )
        
        logger.info(f"Successfully loaded universe with {len(u.atoms)} atoms")
        logger.info(f"Number of frames: {len(u.trajectory)}")
        
        # Create a new universe with rewrapped atoms
        rewrapped_u = create_rewritten_universe(u)
        
        # Save the rewrapped trajectory
        with mda.Writer(rewritten_traj_path, rewrapped_u.atoms.n_atoms) as W:
            for ts in rewrapped_u.trajectory:
                W.write(rewritten_u.atoms)
                
    except Exception as e:
        logger.error(f"Error loading trajectory: {str(e)}")
        # Try alternative loading method
        try:
            u = mda.Universe(
                top_path, 
                traj_path,
                in_memory=True,
                guess_bonds=True
            )
        except Exception as e:
            logger.error(f"Both loading attempts failed. Final error: {str(e)}")
            raise

def align_trajectory(top_path, traj_path, aligned_traj_path, selection='protein and name CA'):
    """
    Align the trajectory to a reference structure.
    Parameters:
    - top_path: path to the topology file (PDB)
    - traj_path: path to the trajectory file (DCD)
    - aligned_traj_path: path to save the aligned trajectory file (DCD)
    """
    u = mda.Universe(top_path, traj_path)

    # Step 1: Calculate the average structure of the protein, as an alignment reference.
    print("Calculating the average structure of the protein...")
    from MDAnalysis.analysis import align
    average = align.AverageStructure(u, u, select=selection,
                                    ref_frame=0,
                                    verbose=True)
    average.run()
    ref = average.universe

    # Step 2: Align all frames to the average structure
    #here we're using carbon-alphas
    aligner = align.AlignTraj(u, ref,
                            select=selection,
                            filename=str(aligned_traj_path),
                            in_memory=False,
                            verbose=True)
    aligner.run()

    u = mda.Universe(top_path, aligned_traj_path)
    print("Number of frames in the trajectory:", len(u.trajectory))
