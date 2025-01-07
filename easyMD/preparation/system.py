import logging
from openmm import app, unit, Platform, openmm
from pathlib import Path
from tempfile import NamedTemporaryFile

from . import structure, ligand, forcefield_utils, rcsb
from .data import common_residues

logger = logging.getLogger(__name__)


def setup_system_from_files(
    protein_input_path: str,
    ligand_sdf_paths: list = [],
    ionic_strength: float = 0.15 * unit.molar,
    box_padding: float = 1.0 * unit.nanometer,
    platform_name: str = "CUDA",
    temperature: float = 300.0,
    friction: float = 1.0,
    timestep: float = 0.002,
    water_model: str = 'tip3p',
    output_pdb: str = "system_solvated.pdb",
    forcefield_files: list = ['amber14-all.xml', 'amber14/tip3p.xml'],
    fix: bool = False
):
    """
    High-level wrapper to set up an OpenMM `Simulation` from a protein (PDB or CIF)
    and (optionally) a ligand SDF. Returns a Simulation object ready for minimization,
    equilibration, or production.
    """
    # Step 1: Convert CIF -> PDB if needed, fix PDB if requested
    protein_pdb_fixed_path = structure.prepare_protein_pdb(
        input_path=protein_input_path,
        fix=fix
    )

    # Step 2: Identify non-standard residues:
    non_standard_residues = forcefield_utils.get_non_standard_residues(forcefield_files, protein_pdb_fixed_path)

    # Step 2: Add CONECT records for any non-standard residues
    # Currently, this overwrites any existing CONECT records.
    # In the future, we should check if the CONECT records are already present, before attempting to add.
    structure.add_conect_for_nonstandard_residues(
        protein_pdb_fixed_path,
        non_standard_residues,
        protein_pdb_fixed_path # We're overwriting our fixed file here.
    )

    # Step 3: Load the protein topology/positions
    pdb = app.PDBFile(protein_pdb_fixed_path)
    protein_topology = pdb.topology
    protein_positions = pdb.positions

    # Step 4: Decide how to handle each ligand.
    # Will we use an SDF, or download from PDB?
    # For now, we just assume each non-standard residue has a proper SDF.
    # for non_standard_residue in non_standard_residues:
    #     logger.warning(
    #         "Found non-standard residue %s matching ligand name %s. "
    #         "Consider providing an SDF for this ligand." % (non_standard_residue[2])
    #     )

    ligand.show_non_standard_residues(protein_topology, non_standard_residues)
    
    # Step 4: Create forcefield. 
    forcefield = forcefield_utils.create_forcefield(forcefield_files)

    ### This will be the eventual way.
    # For now, let's simply add all SDFs to the forcefield.
    # Step 4: Add ligands to forcefield!
    # for non_standard_residue in non_standard_residues:
    #     # - Display structure with TerMol (I need to add simple RDKit input)
    #     #termol.draw(..)
    #     # - Does it match a provided SDF?
    #     # if (matches sdf):
    #     ligand_molecule = ligand.load_openff_ligand_from_sdf(ligand_sdf_path, sanitize=False, removeHs=False)
    #     # - If so, attempt to add to forcefield.
    #     forcefield = add_molecule_to_forcefield(forcefield, ligand_molecule)
    #     # - On fail, clearly alert user and attempt to download _ add from PDB.
    #     # - Otherwise, attempt to download + add from PDB.

    # For each non-standard residue, make sure they have a matched SDF.
    # If not provided, we'll need to download from PDB.
    download_dir = Path("./downloads")
    for non_standard_residue in non_standard_residues:
        chain_id, residue_number, res_name = non_standard_residue
        # Extract the residue to a temporary PDB:
        temp_ligand_file = NamedTemporaryFile(suffix='.pdb', delete=False)
        structure.extract_residue_from_pdb(protein_pdb_fixed_path, chain_id, residue_number, temp_ligand_file.name)

        for ligand_sdf_path in ligand_sdf_paths:
            if ligand.matches_residue_to_sdf(temp_ligand_file.name, ligand_sdf_path):
                logger.info("Matched ligand %s to residue %s." % (ligand_sdf_path, res_name))
                break
        else:
            logger.warning("Could not find a match for residue %s in provided SDFs." % res_name)
            # Attempt to download from PDB:
            try:
                new_ligand_sdf = rcsb.download_ligand(res_name, download_dir)
            except ValueError as e:
                logger.error(e)
                raise
            ligand_sdf_paths.append(new_ligand_sdf)

    for ligand_sdf_path in ligand_sdf_paths:
        ligand_molecule = ligand.load_openff_ligand_from_sdf(ligand_sdf_path, sanitize=False, removeHs=False)
        forcefield = forcefield_utils.add_molecule_to_forcefield(forcefield, ligand_molecule)

    # Step 6: Create modeller, add hydrogens, (optionally) load ligand H definitions, add solvent
    # To create the hydrogen definitions, we need to pass in a residue coupled to an SDF.
    # Which means before we can do this, we need to pair our SDFs with our non-standard residues.
    # This will be a problem I fix in the future!
    modeller = create_modeller_and_solvate(
        pdb_topology=protein_topology,
        pdb_positions=protein_positions,
        non_standard_residues=non_standard_residues,
        ligand_sdf_paths=ligand_sdf_paths,
        water_model=water_model,
        forcefield=forcefield,
        pdb_file=protein_pdb_fixed_path,
        ionic_strength=ionic_strength,
        box_padding=box_padding
    )

    # Write out the solvated structure
    with open(output_pdb, 'w') as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

    # Step 7: Create the System
    system = create_openmm_system(modeller, forcefield)

    # Step 8: Create integrator + Simulation, minimize
    simulation = create_simulation(
        topology=modeller.topology,
        system=system,
        positions=modeller.positions,
        platform_name=platform_name,
        temperature=temperature,
        friction=friction,
        timestep=timestep
    )

    logger.info("Minimizing energy...")
    simulation.minimizeEnergy()
    logger.info("System setup complete.")
    return simulation

def create_modeller_and_solvate(pdb_topology, pdb_positions, forcefield, pdb_file, water_model,
                                non_standard_residues,
                                ligand_sdf_paths,
                                ionic_strength=0.15*unit.molar,
                                box_padding=1.0*unit.nanometer):
    """
    Builds an OpenMM `Modeller` from a PDB topology/positions, optionally
    loading custom hydrogen definitions for a ligand, then solvates.
    """
    from openmm import app
    import tempfile

    modeller = app.Modeller(pdb_topology, pdb_positions)

    # For each non-standard residue, we first extract it to a PDB:
    non_standard_residue_dict = {}
    for non_standard_residue in non_standard_residues:
        chain_id, residue_number, res_name = non_standard_residue

        if res_name in non_standard_residue_dict:
            # Skip if we've already extracted this residue.
            # We're assuming that each residue in the PDB file with the same name 
            # has the exact same atom names/order
            continue

        temp_file = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
        structure.extract_residue_from_pdb(
            input_pdb=pdb_file,
            chain_id=chain_id,
            residue_number=residue_number,
            output_pdb=temp_file.name
        )
        non_standard_residue_dict[res_name] = temp_file.name
    
    # Now, for each extracted ligand, we need to match it to a provided SDF.
    residue_pdb_sdf_dict = {} #contains (ligand_sdf, name) for each PDB
    for res_name, pdb_file in non_standard_residue_dict.items():
        for ligand_sdf_path in ligand_sdf_paths:
            if ligand.matches_residue_to_sdf(pdb_file, ligand_sdf_path):
                logger.info("Matched ligand %s to residue %s." % (ligand_sdf_path, res_name))
                residue_pdb_sdf_dict[pdb_file] = (ligand_sdf_path, res_name)
                break
        else:
            logger.warning("Could not find a match for residue %s in provided SDFs." % res_name)
    
    # Now we have a residue -> SDF mapping, we can create our hydrogen template and add it to the modeller.
    for i, (pdb_file, (ligand_sdf_path, res_name)) in enumerate(residue_pdb_sdf_dict.items()):
        # Add the residue to the modeller
        xml_temp_file = tempfile.NamedTemporaryFile(suffix='.xml', delete=False)
        ligand.generate_hydrogen_template_xml_from_sdf(
            ligand_sdf_path,
            pdb_file,
            xml_temp_file.name,
            residue_name=res_name
        )
        logger.info("Loading hydrogen definitions for ligand %s." % ligand_sdf_path)
        modeller.loadHydrogenDefinitions(xml_temp_file.name)
    
    # Add missing hydrogens for entire complex
    modeller.addHydrogens(forcefield)

    # Add solvent
    logger.info("Adding solvent (model=%s, ionic_strength=%s, padding=%s).",
                water_model, ionic_strength, box_padding)
    modeller.addSolvent(forcefield, model=water_model,
                        ionicStrength=ionic_strength, padding=box_padding)
    return modeller


def create_openmm_system(modeller, forcefield):
    """
    Calls forcefield.createSystem(...) with typical settings to produce an OpenMM System.
    """
    logger.info("Creating OpenMM System.")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        removeCMMotion=True
    )
    return system


def create_simulation(topology, system, positions,
                      platform_name="CUDA",
                      temperature=300.0,
                      friction=1.0,
                      timestep=0.002):
    """
    Creates a Simulation object with a Langevin integrator and sets initial positions.
    """
    integrator = openmm.LangevinIntegrator(
        temperature * unit.kelvin,
        friction / unit.picosecond,
        timestep * unit.picoseconds
    )
    platform = Platform.getPlatformByName(platform_name)
    simulation = app.Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)
    return simulation
