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

def add_hydrogens_to_model( modeller, forcefield, hydrogen_templates_dict, output_pdb=None, pH=7.0, hydrogen_variants=[], skip_residue_h_names=[] ):
    for residue_name, xml_temp_file_path in hydrogen_templates_dict.items():
        ### TODO: Add a skip_h_add flag, so we can use user-generated Hs (and to preserver h-names).
        if residue_name in skip_residue_h_names:
            logger.debug("Skipping hydrogen addition for residue %s." % residue_name)
            continue
        logger.debug("Loading hydrogen definitions for ligand %s." % residue_name)
        modeller.loadHydrogenDefinitions( xml_temp_file_path )

    # Parse variants (from string to list of [None, None, VARIANT, ...])
    if hydrogen_variants:
        logger.debug("Parsing hydrogen variants: %s" % hydrogen_variants)
        variant_list = parse_hydrogen_variants(hydrogen_variants, modeller)
    else:
        logger.debug("Using default protonation states.")
        variant_list = None

    # Add missing hydrogens for entire complex
    logger.debug("Adding hydrogens to the entire system.")
    modeller.addHydrogens(forcefield, pH=pH, variants=variant_list)

    # Write out the hydrogenated structure
    if output_pdb is not None:
        logger.debug("Writing out hydrogenated structure.")
        with open(output_pdb, 'w') as f:
            openmm.app.PDBFile.writeFile(modeller.topology, modeller.positions, f, keepIds=True)

    return modeller

def add_solvent_to_model( modeller, forcefield, water_model, ionic_strength, box_padding, output_pdb=None ):
    logger.debug("Adding solvent (model=%s, ionic_strength=%s, padding=%s).",
                water_model, ionic_strength, box_padding)
    
    if water_model != 'implicit':
        modeller.addSolvent(forcefield, model=water_model,
                            ionicStrength=ionic_strength, padding=box_padding)
    else:
        # Add implicit solvent
        logger.debug("Using implicit solvent (No explicit solvent added).")

    # Write out the solvated structure
    if output_pdb is not None:
        with open(output_pdb, 'w') as f:
            openmm.app.PDBFile.writeFile(modeller.topology, modeller.positions, f, keepIds=True)

    return modeller

def parse_hydrogen_variants(hydrogen_variants, modeller):
    """
    Inputs:
    - hydrogen_variants: list of strings of form [CHAIN:RESID=VARIANT, ...]. Example: ['A:1=HIS', 'C:102=HIE']

    Outputs:
    - list of length n, where n is the number of residues in the modeller.
    - Each entry is a string of the form VARIANT, where VARIANT is the name of the variant, or None if no variant is specified.
    """

    # Create a list of None values for each residue in the modeller
    residue_variants = [None] * len(list(modeller.topology.residues()))

    # Iterate over the hydrogen variants
    for variant in hydrogen_variants:
        # Split the variant string into chain, residue number, and variant name
        chain_residue, variant_name = variant.split('=')
        chain_id, residue_number = chain_residue.split(':')
        residue_number = int(residue_number)
        # Check if the chain ID is valid
        if chain_id not in [chain.id for chain in modeller.topology.chains()]:  
            raise ValueError(f"Invalid chain ID: {chain_id}")
        # Check if the residue number is valid
        if residue_number not in [int(residue.id) for residue in modeller.topology.residues() if residue.chain.id == chain_id]:
            raise ValueError(f"Invalid residue number: {residue_number}")

        # Find the corresponding residue in the modeller
        for i, residue in enumerate(modeller.topology.residues()):
            if (residue.chain.id == chain_id) and (int(residue.id) == residue_number):
                residue_variants[i] = variant_name

    return residue_variants

def parse_custom_bonds(custom_bonds):
    """
    Parses a list of custom bonds from a string to a list of tuples.
    Inputs:
    - custom_bonds: list of strings of form ["atom1,atom2,k,r0" , ...]

    Outputs:
    - list of tuples.
    """


    # Create a list of tuples for each custom bond
    custom_bond_list = []
    for bond in custom_bonds:
        split_string = bond.split(',')
        if len(split_string) != 4:
            raise ValueError("Custom bond string must be of the form 'atom1,atom2,k,r0'")
        # Parse the bond string into its components
        atom1, atom2, k, r0 = split_string
        k = float(k)
        r0 = float(r0)
        custom_bond_list.append((atom1, atom2, k, r0))

    return custom_bond_list

def parse_custom_angles(custom_angles):
    """
    Parses a list of custom angles from a string to a list of tuples.
    Inputs:
    - custom_angles: list of strings of form ["atom1,atom2,atom3,k,theta0" , ...]

    Outputs:
    - list of tuples.
    """

    # Create a list of tuples for each custom angle
    custom_angle_list = []
    for angle in custom_angles:
        split_string = angle.split(',')
        if len(split_string) != 5:
            raise ValueError("Custom angle string must be of the form 'atom1,atom2,atom3,k,theta0'")
        # Parse the angle string into its components
        atom1, atom2, atom3, k, theta0 = split_string
        k = float(k)
        theta0 = float(theta0)
        custom_angle_list.append((atom1, atom2, atom3, k, theta0))

    return custom_angle_list

def parse_custom_torsions(custom_torsions):
    """
    Parses a list of custom torsions from a string to a list of tuples.
    Inputs:
    - custom_torsions: list of strings of form ["atom1,atom2,atom3,atom4,k,phi0" , ...]

    Outputs:
    - list of tuples.
    """

    # Create a list of tuples for each custom torsion
    custom_torsion_list = []
    for torsion in custom_torsions:
        split_string = torsion.split(',')
        if len(split_string) != 7:
            raise ValueError("Custom torsion string must be of the form 'atom1,atom2,atom3,atom4,k,phase'")
        # Parse the torsion string into its components
        atom1, atom2, atom3, atom4, k, periodicity, phase = split_string
        k = float(k)
        phase = float(phase)
        periodicity = float(periodicity)
        custom_torsion_list.append((atom1, atom2, atom3, atom4, k, periodicity, phase))

    return custom_torsion_list

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
    fix: bool = True,
    add_h: bool = True,
    skip_residue_h: list = [],
    hydrogen_variants: list = [],
    keep_temp_files: bool = False,
    custom_bonds: list = [],
    custom_angles: list = [],
    custom_torsions: list = [],
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

    # Parse custom bonds, angles, and torsions:
    custom_bonds = parse_custom_bonds(custom_bonds)
    custom_angles = parse_custom_angles(custom_angles)
    custom_torsions = parse_custom_torsions(custom_torsions)

    # Step 1: Prepare the model, forcefield, and hydrogen templates:
    modeller, forcefield, hydrogen_templates_dict = prepare_model(
        protein_input_path=protein_input_path,
        ligand_sdf_paths=ligand_sdf_paths,
        forcefield_files=forcefield_files,
        fix=fix,
        temp_dir=temp_dir
    )

    # Step 2: Add hydrogens to the model
    if add_h:

        # If there are residues to skip, convert the chain:ID format to resnames:
        skip_residue_h_names = []
        for skip_residue in skip_residue_h:
            # Split the residue into chain and ID
            chain_id, residue_id = skip_residue.split(':')
            # Find the corresponding residue in the modeller
            for residue in modeller.topology.residues():
                if (residue.chain.id == chain_id) and (int(residue.id) == int(residue_id)):
                    skip_residue_h_names.append(residue.name)
                    break
            else:
                raise ValueError(f"Residue {skip_residue} not found in the modeller.")
        if len(skip_residue_h_names) > 0:
            logger.info("Skipping hydrogen addition for residues: %s" % skip_residue_h_names)
        
        modeller = add_hydrogens_to_model(
            modeller=modeller,
            forcefield=forcefield,
            hydrogen_templates_dict=hydrogen_templates_dict,
            output_pdb=None,
            pH=pH,
            hydrogen_variants=hydrogen_variants,
            skip_residue_h_names=skip_residue_h_names
        )
    else:
        logger.info("--no-add-h flag was used. Will not add hydrogens to model.")

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
        nonbondedMethod=openmm.app.PME if water_model != 'implicit' else openmm.app.NoCutoff, # for implicit water, we don't use a periodic system.
        nonbondedCutoff=1.0 * openmm.unit.nanometer,
        constraints=None,
        rigidWater=True,
        removeCMMotion=False,
        hydrogenMass=None,
    )

    # Step 10: Add a custom force to the system
    for custom_bond in custom_bonds:
        system = prep_forcefield.add_harmonic_force(system, modeller, *custom_bond)
    for custom_angle in custom_angles:
        system = prep_forcefield.add_angle_force(system, modeller, *custom_angle)
    for custom_torsion in custom_torsions:
        system = prep_forcefield.add_torsion_force(system, modeller, *custom_torsion)

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