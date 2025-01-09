import argparse
import logging
import sys
from pathlib import Path
from openmm import unit
from easyMD.app import prepare_system, run_simulation, process_trajectory

def parse_args():
    parser = argparse.ArgumentParser(
        description='Run molecular dynamics simulation with easyMD',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        'input_structure',
        type=str,
        help='Input structure file (e.g., PDB, CIF)'
    )
    
    # Optional arguments
    parser.add_argument(
        '--output-pdb',
        type=str,
        default='output.pdb',
        help='Output PDB file name'
    )
    
    parser.add_argument(
        '--output-traj',
        type=str,
        default='output.dcd',
        help='Output trajectory file name'
    )
    
    parser.add_argument(
        '--output-aligned-traj',
        type=str,
        default='output_aligned.dcd',
        help='Output aligned trajectory file name'
    )
    
    parser.add_argument(
        '--duration',
        type=float,
        default=1.0,
        help='Simulation duration in nanoseconds'
    )
    
    parser.add_argument(
        '--output-frequency',
        type=float,
        default=0.1,
        help='Output frequency in nanoseconds'
    )
    
    parser.add_argument(
        '--fix-structure',
        action='store_true',
        help='Apply structure fixes during preparation'
    )
    
    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='Set logging level'
    )
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Set up logging
    logging.basicConfig(level=getattr(logging, args.log_level))
    logging.getLogger('easyMD').setLevel(getattr(logging, args.log_level))
    
    # Convert input paths to Path objects
    input_path = Path(args.input_structure)
    output_pdb = Path(args.output_pdb)
    output_traj = Path(args.output_traj)
    output_aligned = Path(args.output_aligned_traj)
    
    # Verify input file exists
    if not input_path.exists():
        logging.error(f"Input structure file not found: {input_path}")
        sys.exit(1)
    
    # Prepare system
    simulation = prepare_system(
        str(input_path),
        output_pdb=str(output_pdb),
        fix=args.fix_structure
    )

    try:
        # Run simulation
        run_simulation(
            simulation,
            str(output_traj),
            duration=args.duration * unit.nanoseconds,
            output_frequency=args.output_frequency * unit.nanoseconds
        )
    
    except KeyboardInterrupt:
        logging.info("Simulation interrupted by user")
    except Exception as e:
        logging.error(f"Error during simulation: {str(e)}")
        raise

    # Process trajectory
    process_trajectory(
        str(output_pdb),
        str(output_traj),
        str(output_aligned)
    )

if __name__ == "__main__":
    main()
