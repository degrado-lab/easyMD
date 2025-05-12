import logging
import sys
from pathlib import Path

import typer
from typing import List
from click import Context
from typer.core import TyperGroup
from typing_extensions import Annotated

class OrderCommands(TyperGroup):
  def list_commands(self, ctx: Context):
    return list(self.commands)

app = typer.Typer(help="Run molecular dynamics simulations with easyMD.", add_completion=False, no_args_is_help=True, cls=OrderCommands, rich_markup_mode="rich")
# Listing here so we can check against user inputs later:
DEFAULT_FORCEFIELD = ["amber14-all.xml", "amber14/tip3p.xml"]
DEFAULT_WATER_MODEL = 'tip3p'
DEFAULT_IONIC_STRENGTH = 0.15 #* unit.molar
DEFAULT_BOX_PADDING = 1.0 #* unit.nanometer


@app.command(rich_help_panel='Simulation')
def run(
    input_structure: str = typer.Argument(
        ...,
        help="Input structure file (e.g., PDB, CIF)",
        show_default=False,
    ),
    # Optional list of ligand SDF files:
    ligands: List[str] = typer.Option(
        [],
        "--ligand",
        "-l",
        help="The structure file (e.g., SDF) of a ligand in your system. This file specifies the correct bonds and protonation states. For multiple, use the flag multiple times.",
        rich_help_panel="Files",
    ),
    output: str = typer.Option(
        "./output",
        help="Prefix for output files. EasyMD will produce [output].pdb, [output].dcd, and [output]_aligned.dcd.",
        rich_help_panel="Files",
    ),
    duration: float = typer.Option(
        10,
        help="Simulation duration in nanoseconds",
        rich_help_panel="Simulation Settings",
    ),
    relax_duration: float = typer.Option(
        1.0,
        help="Duration of initial equilibration/relaxation simulation, in nanoseconds",
        rich_help_panel="Simulation Settings",
    ),
    output_frequency: float = typer.Option(
        1,
        help="Output frequency in nanoseconds",
        rich_help_panel="Simulation Settings",
    ),
    fix: bool = typer.Option(
        True,
        help="Fix the structure during preparation with PDBFixer",
        rich_help_panel="System Preparation",
    ),
    add_h: bool = typer.Option(
        True,
        help="Add hydrogens during preparation.",
        rich_help_panel="System Preparation",
    ),
    skip_residue_h: List[str] = typer.Option(
        [],
        "--skip-residue-h",
        "-srh",
        help="List of residues to skip adding hydrogens to. Specify the chain and residue number, then the variant. E.g. A:13, or B:98. \n \
            Use multiple times for multiple residues.",
        rich_help_panel="System Preparation",
    ),
    forcefield: List[str] = typer.Option(
        DEFAULT_FORCEFIELD,
        "--forcefield",
        "-f",
        help="List of forcefield files (e.g., amber14-all.xml amber14/tip3p.xml).\n \
                To list available files, run `easymd forcefields`. \n \
                Use multiple times for multiple forcefield files.",
        rich_help_panel="Forcefield",
    ),
    water_model: str = typer.Option(
        DEFAULT_WATER_MODEL,
        help = 'Water model used for solvation. Must match forcefield parameters. (tip3p|spce|tip4pew|tip5p|swm4ndp|implicit).',
        rich_help_panel="System Preparation",
    ),
    pH: float = typer.Option(
        7.0,
        help = 'pH of the system. Used for protonation state calculations.',
        rich_help_panel="System Preparation",
    ),
    hydrogen_variants: List[str] = typer.Option(
        [],
        "--hydrogen-variant",
        "-hv",
        help="List of hydrogen variants to use. Specify the chain and residue number, then the variant. E.g. A:13=HIE, or B:98=ASH. \n \
            Use multiple times for multiple variants.",
        rich_help_panel="System Preparation",
    ),
    ionic_strength: float = typer.Option(
        DEFAULT_IONIC_STRENGTH,
        help = 'Ionic strength used in solvation of the system (Molar).',
        rich_help_panel="System Preparation",
    ),
    box_padding: float = typer.Option(
        DEFAULT_BOX_PADDING,
        help = 'Padding of solvent around the system (nanometers).',
        rich_help_panel="System Preparation",
    ),
    custom_bonds: List[str] = typer.Option(
        [],
        "--custom-bond",
        "-cb",
        help="List of custom harmonic bonds to add. Specify the two atoms in the form 'chain_id:residue_id:atom_name', then the k-value (units kcal/(mol*A^2)), then the target distance (units A). \n \
            E.g. A:99:CG,B:99:CG,10.0,2.0 \n \
            Use multiple times for multiple bonds.",
        rich_help_panel="Forcefield",
    ),
    custom_angles: List[str] = typer.Option(
        [],
        "--custom-angle",
        "-ca",
        help="List of custom harmonic angles to add. Specify the three atoms in the form 'chain_id:residue_id:atom_name', then the k-value (units kcal/mol), the periodicity, and the phase value (units DEGREES). \n \
            E.g. A:99:CG,B:99:CG,C:99:CG,10.0,1,180 \n \
            Use multiple times for multiple angles.",
        rich_help_panel="Forcefield",
    ),
    custom_torsions: List[str] = typer.Option(
        [],
        "--custom-torsion",
        "-ct",
        help="List of custom harmonic torsions to add. Specify the four atoms in the form 'chain_id:residue_id:atom_name', then the k-value (units kcal/(mol*RADIANs^2)), the periodicity, then the target angle (units DEGREES). \n \
            E.g. A:99:CG,B:99:CG,C:99:CG,D:99:CG,10.0,1,109.5 \n \
            Use multiple times for multiple torsions.",
        rich_help_panel="Forcefield",
    ),
    minimize_only: bool = typer.Option(
        False,
        help="Only minimize the system and output the minimized structure. No simulation will be run.",
        rich_help_panel="Simulation Settings",
    ),
    # minimize_only_with_traj: bool = typer.Option(
    #     False,
    #     help="Only minimize the system and output a minimized trajectory. No simulation will be run. (Primarily for visualization.)",
    #     rich_help_panel="Simulation Settings",
    # ),
    log_level: str = typer.Option(
        "INFO",
        help="Set logging level (DEBUG|INFO|WARNING|ERROR)",
        case_sensitive=False
    )
):
    """
    Prepare and run a molecular dynamics simulation.\n
    Example usage: easymd run protein_ligand.cif -l ligand.sdf --output ./output/out --duration 100 --output-frequency 1
    """
    # Set up logging
    logging.basicConfig(level=getattr(logging, log_level.upper()))
    logging.getLogger('easyMD').setLevel(getattr(logging, log_level.upper()))

    # If we changed the default ForceField or Water model, remind the user to make sure they match.
    if forcefield != DEFAULT_FORCEFIELD or water_model != DEFAULT_WATER_MODEL:
        logging.warning('Default forcefield or water model has been changed. Ensure the forcefield parameters and water model for solvation match.')
    if water_model == 'implicit' and ligands != []:
        raise ValueError('Currently, when using implicit water model, ligands cannot be added dynamically to the forcefield (e.g., via SDF). \n \
                         Please use forcefield XML files that include the ligand parameters.')

    # Convert input/output paths to Path objects
    input_path = Path(input_structure)
    output_name = Path(output).with_suffix('') # remove extension in case they provided one.
    output_name.parent.mkdir(parents=True, exist_ok=True)
    output_pdb_path =       output_name.with_suffix('.pdb')
    output_traj_path =      output_name.with_suffix('.dcd')
    output_aligned_path =   output_name.with_name(f"{output_name.stem}_aligned.dcd")
    output_EM_pdb_path =   output_name.with_name(f"{output_name.stem}_EM.pdb")
    output_EM_traj_path =  output_name.with_name(f"{output_name.stem}_EM.dcd")

    # Verify input file exists
    if not input_path.exists():
        logging.error(f"Input structure file not found: {input_path}")
        sys.exit(1)

    from openmm import unit

    from easyMD.app import prepare_system, run_simulation, process_trajectory, energy_minimize_with_trajectory
    # Prepare system
    system, modeller = prepare_system(
        protein_input_path = str(input_path),
        ligand_sdf_paths = ligands,
        ionic_strength = ionic_strength * unit.molar,
        box_padding = box_padding * unit.nanometer,
        water_model = water_model,
        pH = pH,
        output_pdb=str(output_pdb_path),
        forcefield_files = forcefield,
        fix = fix,
        add_h = add_h,
        skip_residue_h = skip_residue_h,
        hydrogen_variants = hydrogen_variants,
        custom_bonds = custom_bonds,
        custom_angles = custom_angles,
        custom_torsions = custom_torsions
    )
    
    ### TODO: Haven't finished minimize_only_with_traj, the EM visualization function.
    minimize_traj_only = False
    
    try:
        # Run simulation
        run_em = True
        run_relax = True
        run_sim = True

        if minimize_only:
            run_relax = False
            run_sim = False

        if not minimize_traj_only:
            run_simulation(
                system=system,
                modeller=modeller,
                output_file=str(output_traj_path),
                duration=duration * unit.nanoseconds,
                output_frequency=output_frequency * unit.nanoseconds,
                energy_minimize=run_em,
                relax=run_relax,
                simulate=run_sim,
                relax_duration = relax_duration * unit.nanoseconds,
                minimized_pdb= str(output_EM_pdb_path) if minimize_only else None,
            )
        else:
            energy_minimize_with_trajectory(
                system=system,
                modeller=modeller,
                output_pdb=str(output_EM_pdb_path),
                output_traj=str(output_EM_traj_path),
            )
    except KeyboardInterrupt:
        logging.info("Simulation interrupted by user")
    except Exception as e:
        logging.error(f"Error during simulation: {str(e)}")
        raise

    # Process trajectory
    if not minimize_only and not (water_model == 'implicit'): # This isn't implemented for implicit solvent
        process_trajectory(
            str(output_pdb_path),
            str(output_traj_path),
            str(output_aligned_path)
        )
    
    logging.info("Done!")

@app.command(rich_help_panel='Utility')
def reduce(
    input_structure: str = typer.Argument(
        ...,
        help="Input structure file (e.g., PDB, CIF)",
        show_default=False
    ),
    # Optional list of ligand SDF files:
    ligands: List[str] = typer.Option(
        [],
        "--ligand",
        "-l",
        help="The structure file (e.g., SDF) of a ligand in your system. This file specifies the correct bonds and protonation states. For multiple, use the flag multiple times.",
    ),
    output_file: str = typer.Option(
        None,
        help="Output file for the reduced structure. If not provided, the output will be named [input]_h.pdb."
    ),
    fix: bool = typer.Option(
        True,
        help="Fix the structure during preparation with PDBFixer"
    ),
    forcefield: List[str] = typer.Option(
        DEFAULT_FORCEFIELD,
        "--forcefield",
        "-f",
        help="List of forcefield files (e.g., amber14-all.xml amber14/tip3p.xml).\n \
                To list available files, run `easymd get-forcefields`."
    ),
    log_level: str = typer.Option(
        "INFO",
        help="Set logging level (DEBUG|INFO|WARNING|ERROR)",
        case_sensitive=False
    )
):
    """
    Add hydrogens to a structure and save as a PDB file.
    """
    # Set up logging
    logging.basicConfig(level=getattr(logging, log_level.upper()))
    logging.getLogger('easyMD').setLevel(getattr(logging, log_level.upper()))

    # Convert input/output paths to Path objects
    input_path = Path(input_structure)
    output_file = output_file or input_path.with_name(f"{input_path.stem}_h.pdb")
    output_file = Path(output_file)

    # Verify input file exists
    if not input_path.exists():
        logging.error(f"Input structure file not found: {input_path}")
        sys.exit(1)

    # Add hydrogens:
    from easyMD.app.prepare_system import reduce
    reduce(
        protein_input_path=str(input_path),
        ligand_sdf_paths=ligands,
        output_pdb=str(output_file),
        forcefield_files=forcefield,
        fix=fix
    )
    logging.info("Done!")

@app.command(rich_help_panel='Utility')
def process(
    top_path: str = typer.Argument(
        ..., help="Path to the topology file (e.g., PDB)", show_default=False
    ),
    traj_path: str = typer.Argument(
        ..., help="Path to the trajectory file (e.g., DCD)", show_default=False
    ),
    output_prefix: str = typer.Option(
        "./output_processed",
        help="Prefix for output files. EasyMD will produce [output_prefix]_rewrapped.dcd and [output_prefix]_aligned.dcd."
    ),
    log_level: str = typer.Option(
        "INFO",
        help="Set logging level (DEBUG|INFO|WARNING|ERROR)",
        case_sensitive=False
    )
):
    """
    Align and re-wrap a trajectory around a protein.
    """
    # Set up logging
    logging.basicConfig(level=getattr(logging, log_level.upper()))
    logging.getLogger('easyMD').setLevel(getattr(logging, log_level.upper()))

    # Convert input/output paths to Path objects
    top_path = Path(top_path)
    traj_path = Path(traj_path)
    output_prefix = Path(output_prefix).with_suffix('')
    rewrapped_traj_path = output_prefix.with_name(f"{output_prefix.stem}_rewrapped.dcd")
    aligned_traj_path = output_prefix.with_name(f"{output_prefix.stem}_aligned.dcd")

    # Verify input files exist
    if not top_path.exists() or not traj_path.exists():
        logging.error(f"Input files not found: {top_path}, {traj_path}")
        sys.exit(1)

    # Process trajectory
    from easyMD.analysis.trajectory import rewrap_trajectory, align_trajectory
    rewrap_trajectory(str(top_path), str(traj_path), str(rewrapped_traj_path))
    align_trajectory(str(top_path), str(rewrapped_traj_path), str(aligned_traj_path))

    logging.info("Processing complete!")

@app.command(rich_help_panel='Info')
def get_forcefields():
    """
    Display available OpenMM Forcefields.
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

@app.command(rich_help_panel='Info')
def get_hydrogen_variants():
    """
    Display the valid protonation states for amino acids in OpenMM.
    """
    print("From OpenMM documentation:")
    print("""
Aspartic acid:
    ASH: Neutral form with a hydrogen on one of the delta oxygens
    ASP: Negatively charged form without a hydrogen on either delta oxygen

Cysteine:
    CYS: Neutral form with a hydrogen on the sulfur
    CYX: No hydrogen on the sulfur (either negatively charged, or part of a disulfide bond)

Glutamic acid:
    GLH: Neutral form with a hydrogen on one of the epsilon oxygens
    GLU: Negatively charged form without a hydrogen on either epsilon oxygen

Histidine:
    HID: Neutral form with a hydrogen on the ND1 atom
    HIE: Neutral form with a hydrogen on the NE2 atom
    HIP: Positively charged form with hydrogens on both ND1 and NE2
    HIN: Negatively charged form without a hydrogen on either ND1 or NE2

Lysine:
    LYN: Neutral form with two hydrogens on the zeta nitrogen
    LYS: Positively charged form with three hydrogens on the zeta nitrogen
""")
    
def main():
    app()

if __name__ == "__main__":
    main()