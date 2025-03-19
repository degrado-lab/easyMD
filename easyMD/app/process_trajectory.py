import logging
import tempfile

logger = logging.getLogger(__name__)

def process_trajectory(top_file, trajectory_file, output_file):
    from easyMD.analysis.trajectory import rewrap_trajectory, align_trajectory

    # Create a temporary file for the rewrapped trajectory
    rewrapped_traj_path = tempfile.NamedTemporaryFile(delete=False, suffix='.dcd').name
    rewrap_trajectory(top_file, trajectory_file, rewrapped_traj_path, perform_initial_rewrap=False)
    align_trajectory(top_file, rewrapped_traj_path, output_file)