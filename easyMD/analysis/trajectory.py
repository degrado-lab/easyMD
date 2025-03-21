# my_md_package/simulation/trajectory.py
import logging
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np
from tqdm import tqdm
import tempfile

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
    logger.debug(f"Number of frames in the trajectory: {len(u.trajectory)}")

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
    - selection: selection string for the group of atoms to center on (default: 'protein')
    Returns:
    - new_u: MDAnalysis Universe object with rewrapped positions
    """
    new_positions = []  # To store rewrapped positions for all frames
    new_dimensions = [] # To store dimensions for all frames

    # Iterate through each frame and apply re-wrapping
    for ts in tqdm(u.trajectory):
        box = ts.dimensions[:3]  # Get current box dimensions
        center = u.select_atoms(selection).center_of_geometry()
        wrapped_positions = rewrap_positions(ts.positions, center, box)
        new_positions.append(wrapped_positions)
        new_dimensions.append(ts.dimensions)
    # Convert list of positions to a 3D numpy array
    new_positions_array = np.array(new_positions)

    # Create a new Universe with the rewrapped positions using MemoryReader
    new_u = mda.Universe(u.filename, u.trajectory.filename)
    new_u.load_new(new_positions_array, format=MemoryReader)

    # # Set the dimensions for each frame, from previous universe
    for ts in new_u.trajectory:
        ts.dimensions = new_dimensions[ts.frame]

    return new_u

# def fix_waters_on_boundary(u, cutoff=10.0):
#     # Select the water residues
#     waters = u.select_atoms("resname WAT or resname TIP3 or resname HOH")
#     water_residues = waters.residues  # Each residue is presumably 3 atoms of water

#     for ts in tqdm(u.trajectory):
#         box = ts.dimensions[:3]  # (Lx, Ly, Lz)
#         all_coords = u.atoms.positions  # <--- entire Universe, shape = (N_universe, 3)

#         # For each water residue, gather the global indices
#         #for res in water_residues:
#         ridx = water_residues.atoms.indices  # global indices in u.atoms
#         coords_res = all_coords[ridx]  # shape = (n_atoms_in_res, 3)

#         # Check dimension by dimension
#         boundary_mask = coords_res > (box - cutoff)

#         #print(box[0], box[1], box[2])
#         print(boundary_mask.shape)
#         print(coords_res[:, 0][boundary_mask[:, 0]])
#         print(box[0])
#         coords_res[:, 0][boundary_mask[:, 0]] -= box[0]
#         coords_res[:, 1][boundary_mask[:, 1]] -= box[1]
#         coords_res[:, 2][boundary_mask[:, 2]] -= box[2]
        
#         print(coords_res[:, 0][boundary_mask[:, 0]])
#         print()
#         all_coords[ridx] = coords_res

#         # Put updated coords back
#         u.atoms.positions = all_coords

#     return u

def rewrap_trajectory(top_path, traj_path, rewrapped_traj_path, perform_initial_rewrap=True, selection='protein'):
    """
    Rewrap the trajectory within the periodic boundaries.
    Parameters:
    - top_path: path to the topology file (PDB)
    - traj_path: path to the trajectory file (DCD)
    - rewrapped_traj_path: path to save the rewrapped trajectory file (DCD)
    """
    # Add debug information
    logger.debug(f"Loading trajectory from {traj_path}")
    logger.debug(f"Using topology from {top_path}")

    # Reduce frames:
    #reduce_frames(top_path, traj_path, 'output_red.dcd', 5)
    
    # Try loading with topology_format explicitly set for OpenMM PDB
    u = mda.Universe(
        top_path, 
        traj_path,
        topology_format='PDB',
        format='DCD',
        permissive=True
    )

    logger.debug(f"Successfully loaded universe with {len(u.atoms)} atoms")
    logger.debug(f"Number of frames: {len(u.trajectory)}")
    
    # We'll initially rewrap on the first chain. This should push everything roughly into the box.
    if perform_initial_rewrap:
        logger.debug("Performing initial rewrap on the first chain...")
        initial_rewrapped_u = create_rewrapped_universe(u, selection='segid A')
        # Make a named temporary file:
        temp_dcd = tempfile.NamedTemporaryFile(delete=False, suffix='.dcd').name

        logger.debug(f"Saving initial rewrap to temporary file: {temp_dcd}")
        with mda.Writer(temp_dcd, initial_rewrapped_u.atoms.n_atoms) as W:
            for ts in initial_rewrapped_u.trajectory:
                W.write(initial_rewrapped_u.atoms)
        
        # Reload from temporary file to ensure clean state
        initial_rewrapped_u = mda.Universe(top_path, temp_dcd)
    else:
        logger.debug("Skipping initial rewrap.")
        initial_rewrapped_u = u
    
    # Then, we rewrap on the chosen selection.
    rewrapped_u = create_rewrapped_universe(initial_rewrapped_u, selection=selection)

    # Fix the waters:
    # rewrapped_u = fix_waters_on_boundary(rewrapped_u)

    # Save the rewrapped trajectory
    with mda.Writer(rewrapped_traj_path, rewrapped_u.atoms.n_atoms) as W:
        for ts in rewrapped_u.trajectory:
            W.write(rewrapped_u.atoms)
    
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
    logger.debug("Calculating the average structure of the protein...")
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
    logger.debug(f"Number of frames in the trajectory: {len(u.trajectory)}")
