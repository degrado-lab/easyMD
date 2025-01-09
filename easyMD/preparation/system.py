import logging
from openmm import app, unit, Platform, openmm
from pathlib import Path
from tempfile import NamedTemporaryFile

from . import structure, ligand, forcefield_utils, rcsb
from .data import common_residues

logger = logging.getLogger(__name__)

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


def create_simulation(system, modeller,
                      platform_name="CUDA",
                      temperature=300.0 * unit.kelvin,
                      friction=1.0 / unit.picosecond,
                      timestep=0.002 * unit.picoseconds):
    """
    Creates a Simulation object with a Langevin integrator and sets initial positions.
    """
    topology, positions = modeller.topology, modeller.positions
    integrator = openmm.LangevinIntegrator(
        temperature,
        friction,
        timestep 
    )
    platform = Platform.getPlatformByName(platform_name)
    simulation = app.Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)
    return simulation
