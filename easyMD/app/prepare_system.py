from easyMD.prep import forcefield as prep_forcefield
from easyMD.prep import structure as prep_structure
from easyMD.prep import ligand as prep_ligand
from easyMD.prep import system as prep_system
import openmm.app
from openmm import unit
from pathlib import Path
from tempfile import mkdtemp
import logging

logger = logging.getLogger(__name__)

'''
I think I want to split this a bit more:
- prepare_model:
    Should take all my inputs and prepare the model (via adding Hs, solvating, etc.)
    Options to add Hs, solvate, etc. And where to write out the intermediates (if I want)
    This will let me use the same fx for JUST adding Hs, JUST solvating, etc.
    Returns the modeller object - this is all I need to create the simulation!
- prepare_simulation:
    Takes the modeller object and creates the simulation object.
    This is where I can set up the integrator, platform, etc.
    This needs to be flexible so I can set up different integrators with different parameters.
    Returns the simulation object.
- run_simulation:
    Should take in a Modeller object and parameters about the simulations needed.
    We should create then run the simulation(s) here.
    Does that make sense?
'''
def prepare_model(
    protein_input_path: str,
    ligand_sdf_paths: list = [],
    forcefield_files: list = ['amber14-all.xml', 'amber14/tip3p.xml'],
    fix: bool = False,
    temp_dir: str = None,
    show_non_standard_residues: bool = True
):
    """
    High-level wrapper to set up an OpenMM `Modeller` from a protein (PDB or CIF)
    and (optionally) a set of ligand SDFs. Returns a Simulation object ready for minimization,
    equilibration, or production.
    """

    # Create a temporary directory for storing temporary files:
    if temp_dir is None:
        temp_dir = mkdtemp()
    temp_dir = Path(temp_dir)
    temp_dir.mkdir(exist_ok=True, parents=True)

    # Step 1: Convert CIF -> PDB if needed, fix PDB if requested
    protein_pdb_fixed_path = prep_structure.prepare_protein_pdb(
        input_path=protein_input_path,
        fix=fix,
        temp_dir=temp_dir
    )

    # Step 2: Identify non-standard residues:
    non_standard_residues = prep_forcefield.get_non_standard_residues(forcefield_files, protein_pdb_fixed_path)
    
    # Step 3: Add CONECT records for any non-standard residues
    ### TODO: This step is very yeehaw. 
    # Under the hood, we're using RDKit to infer bond connectivity from the 3D coordinates in the PDB.
    # This is used to later match the PDB residues with the ligand SDFs, and won't work without it. However, we're kinda comparing apples to oranges.
    # This should be changed in the future.
    # Additionally, this overwrites any existing CONECT records.
    # In the future, we should check if the CONECT records are already present, before attempting to add.
    prep_structure.add_conect_for_nonstandard_residues(
        protein_pdb_fixed_path,
        non_standard_residues,
        protein_pdb_fixed_path # We're overwriting our fixed file here.
    )
    prep_ligand.show_non_standard_residues(protein_pdb_fixed_path, non_standard_residues)

    # Step 4: Prepare ligand files.
    # This means checking/downloading ligand SDFs,
    # And also extracting a copy of each non-standard residue from our input PDBs to a unique temporary PDB.
    # This creates a dictionary mapping residue names to (ligand SDF path, ligand PDB path) tuples.
    NS_residue_dict = prep_ligand.prepare_ligand_files(
        protein_pdb_fixed_path,
        non_standard_residues,
        ligand_sdf_paths,
        temp_dir)

    # Step 5: Create our OpenMM ForceField object, which includes all non-standard residues using GAFF.
    # Soon, we'll add options for other forcefields.
    forcefield = prep_forcefield.make_forcefield( forcefield_files, NS_residue_dict)

    # Step 6: Create hydrogen templates
    # This tells the Modeller object how to add hydrogens to all non-standard residues.
    hydrogen_templates_dict = prep_ligand.write_all_hydrogen_templates( NS_residue_dict, temp_dir )

    # Step 7: Create modeller, add hydrogens, (optionally) load ligand H definitions, add solvent
    # To create the hydrogen definitions, we need to pass in a residue coupled to an SDF.
    # Which means before we can do this, we need to pair our SDFs with our non-standard residues.
    # This will be a problem I fix in the future!
    modeller = prep_system.create_modeller( protein_pdb_fixed_path )

    return modeller, forcefield, hydrogen_templates_dict

def add_hydrogens_to_model( modeller, forcefield, hydrogen_templates_dict, output_pdb=None, pH=7.0 ):
    for residue_name, xml_temp_file_path in hydrogen_templates_dict.items():
        logger.debug("Loading hydrogen definitions for ligand %s." % residue_name)
        modeller.loadHydrogenDefinitions( xml_temp_file_path )

    # Add missing hydrogens for entire complex
    logger.debug("Adding hydrogens to the entire system.")
    modeller.addHydrogens(forcefield, pH=pH)

    # Write out the hydrogenated structure
    if output_pdb is not None:
        logger.debug("Writing out hydrogenated structure.")
        with open(output_pdb, 'w') as f:
            openmm.app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

    return modeller

def add_solvent_to_model( modeller, forcefield, water_model, ionic_strength, box_padding, output_pdb=None ):
    logger.debug("Adding solvent (model=%s, ionic_strength=%s, padding=%s).",
                water_model, ionic_strength, box_padding)
    modeller.addSolvent(forcefield, model=water_model,
                        ionicStrength=ionic_strength, padding=box_padding)

    # Write out the solvated structure
    if output_pdb is not None:
        with open(output_pdb, 'w') as f:
            openmm.app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

    return modeller

# def make_sim(modeller, forcefield):
#     # Step 9: Create the OpenMM System
#     system = prep_system.create_openmm_system(modeller, forcefield)
    
#     # Step 10: Create integrator + Simulation
#     simulation = prep_system.create_simulation(
#         system=system,
#         modeller=modeller,
#         platform_name=platform_name,
#         temperature=temperature,
#         friction=friction,
#         timestep=timestep
#     )

#     logger.info("System setup complete.")
#     return simulation

def prepare_system(
    protein_input_path: str,
    ligand_sdf_paths: list = [],
    ionic_strength: float = 0.15 * unit.molar,
    box_padding: float = 1.0 * unit.nanometer,
    pH: float = 7.0,
    # platform_name: str = "CUDA",
    # temperature: float = 300.0,
    # friction: float = 1.0,
    # timestep: float = 0.002,
    water_model: str = 'tip3p',
    output_pdb: str = "system_solvated.pdb",
    forcefield_files: list = ['amber14-all.xml', 'amber14/tip3p.xml'],
    fix: bool = False,
    keep_temp_files: bool = False
):
    """
    High-level wrapper to set up an OpenMM `Simulation` from a protein (PDB or CIF)
    and (optionally) a set of ligand SDFs. Returns a Simulation object ready for minimization,
    equilibration, or production.
    """
    
    # Step 0: 
    # Create a temporary directory for storing temporary files:
    if not keep_temp_files:
        temp_dir = mkdtemp()
    else:
        temp_dir = Path('prep')

    # Step 1: Prepare the model, forcefield, and hydrogen templates:
    modeller, forcefield, hydrogen_templates_dict = prepare_model(
        protein_input_path=protein_input_path,
        ligand_sdf_paths=ligand_sdf_paths,
        forcefield_files=forcefield_files,
        fix=fix,
        temp_dir=temp_dir
    )

    # Step 2: Add hydrogens to the model
    modeller = add_hydrogens_to_model(
        modeller=modeller,
        forcefield=forcefield,
        hydrogen_templates_dict=hydrogen_templates_dict,
        output_pdb=None,
        pH=pH
    )

    # Step 3: Add solvent to the model
    modeller = add_solvent_to_model(
        modeller=modeller,
        forcefield=forcefield,
        water_model=water_model,
        ionic_strength=ionic_strength,
        box_padding=box_padding,
        output_pdb=output_pdb
    )

    # Step 9: Create the OpenMM Systemv
    # TODO:
    # These parameters should be passed in as arguments, so users can change things like PME if they want.
    logger.debug("Creating OpenMM System.")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=openmm.app.PME,
        nonbondedCutoff=1.0 * openmm.unit.nanometer,
        constraints=None,
        rigidWater=True,
        removeCMMotion=False,
        hydrogenMass=None,
    )

    logger.info("System setup complete.")
    return system, modeller

def reduce(
    protein_input_path: str,
    ligand_sdf_paths: list = [],
    output_pdb: str = None,
    forcefield_files: list = ['amber14-all.xml', 'amber14/tip3p.xml'],
    pH: float = 7.0,
    fix: bool = False,
):
    """
    This is nearly the same as `prepare_system`, but it only adds hydrogens and outputs a new PDB.
    """

    if output_pdb is None:
        output_pdb = Path(protein_input_path).stem + '_h.pdb'

    # Step 1: Prepare the model, forcefield, and hydrogen templates:
    modeller, forcefield, hydrogen_templates_dict = prepare_model(
        protein_input_path=protein_input_path,
        ligand_sdf_paths=ligand_sdf_paths,
        forcefield_files=forcefield_files,
        fix=fix,
    )

    # Step 2: Add hydrogens to the model
    modeller = add_hydrogens_to_model(
        modeller=modeller,
        forcefield=forcefield,
        hydrogen_templates_dict=hydrogen_templates_dict,
        output_pdb=output_pdb,
        pH=pH
    )