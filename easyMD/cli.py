import logging
import sys
from pathlib import Path

import typer
from typing import List
from openmm import unit

from easyMD.app import prepare_system, run_simulation, process_trajectory

app = typer.Typer(help="Run molecular dynamics simulations with easyMD", add_completion=False)
# Listing here so we can check against user inputs later:
DEFAULT_FORCEFIELD = ["amber14-all.xml", "amber14/tip3p.xml"]
DEFAULT_WATER_MODEL = 'tip3p'


@app.command()
def run(
    input_structure: str = typer.Argument(
        ...,
        help="Input structure file (e.g., PDB, CIF)",
        show_default=False
    ),
    output: str = typer.Option(
        "./output",
        help="Prefix for output files. EasyMD will produce [output].pdb, [output].dcd, and [output]_aligned.dcd."
    ),
    duration: float = typer.Option(
        1.0,
        help="Simulation duration in nanoseconds"
    ),
    output_frequency: float = typer.Option(
        0.1,
        help="Output frequency in nanoseconds"
    ),
    fix: bool = typer.Option(
        True,
        help="Fix the structure during preparation with PDBFixer"
    ),
    forcefield: List[str] = typer.Option(
        DEFAULT_FORCEFIELD,
        help="List of forcefield files (e.g., amber14-all.xml amber14/tip3p.xml).\n \
                To list available files, run `easymd forcefields`."
    ),
    water_model: str = typer.Option(
        DEFAULT_WATER_MODEL,
        help = 'Water model used for solvation. Must match forcefield parameters. (tip3p|spce|tip4pew|tip5p|swm4ndp)'
    ),
    log_level: str = typer.Option(
        "INFO",
        help="Set logging level (DEBUG|INFO|WARNING|ERROR)",
        case_sensitive=False
    )
):
    """
    Command to run a molecular dynamics simulation, specifying input files,
    output files, simulation parameters, and logging.
    """
    # Set up logging
    logging.basicConfig(level=getattr(logging, log_level.upper()))
    logging.getLogger('easyMD').setLevel(getattr(logging, log_level.upper()))

    # If we changed the default ForceField or Water model, remind the user to make sure they match.
    if forcefield != DEFAULT_FORCEFIELD or water_model != DEFAULT_WATER_MODEL:
        logging.warning('Default forcefield or water model has been changed. Ensure the forcefield parameters and water model for solvation match.')

    # Convert input/output paths to Path objects
    input_path = Path(input_structure)
    output_name = Path(output).with_suffix('') # remove extension in case they provided one.
    output_pdb_path =       output_name.with_suffix('.pdb')
    output_traj_path =      output_name.with_suffix('.dcd')
    output_aligned_path =   output_name.with_name(f"{output_name.stem}_aligned.dcd")

    # Verify input file exists
    if not input_path.exists():
        logging.error(f"Input structure file not found: {input_path}")
        sys.exit(1)

    # Prepare system
    simulation = prepare_system(
        str(input_path),
        output_pdb=str(output_pdb_path),
        fix=fix,
        water_model=water_model
    )

    try:
        # Run simulation
        run_simulation(
            simulation,
            str(output_traj_path),
            duration=duration * unit.nanoseconds,
            output_frequency=output_frequency * unit.nanoseconds
        )
    except KeyboardInterrupt:
        logging.info("Simulation interrupted by user")
    except Exception as e:
        logging.error(f"Error during simulation: {str(e)}")
        raise

    # Process trajectory
    process_trajectory(
        str(output_pdb_path),
        str(output_traj_path),
        str(output_aligned_path)
    )

    logging.info("Done!")

@app.command()
def get_forcefields():
    """
    Display the path to the OpenMM forcefields.
    """
    import openmm
    openmm_path = Path(openmm.__file__).parent
    forcefield_path = openmm_path / 'app' / 'data'
    print_directory_tree(forcefield_path)

def print_directory_tree(directory: Path, prefix: str = "") -> None:
    """
    Print a directory tree structure starting from `directory` and 
    recursing through all subdirectories.

    :param directory: Base directory to traverse.
    :param prefix: (Internal) String used for indentation in recursive calls.
    """
    # Sort entries so that directories appear before files
    entries = sorted(directory.iterdir(), key=lambda e: (not e.is_dir(), e.name))
    last_entry = len(entries) - 1

    for index, entry in enumerate(entries):
        # Determine connector
        connector = "└──" if index == last_entry else "├──"
        
        if entry.is_dir():
            # Mark directories with a slash for clarity
            typer.echo(f"{prefix}{connector} {entry.name}/")
            # For the next level, extend the prefix appropriately
            extension = "    " if index == last_entry else "│   "
            print_directory_tree(entry, prefix + extension)
        else:
            # Files are printed as-is
            if entry.suffix == '.xml':
                typer.echo(f"{prefix}{connector} {entry.name}")

def main():
    app()

if __name__ == "__main__":
    main()