import openmm.app
from .data import common_residues
from easyMD.prep import structure as prep_structure
from easyMD.prep import ligand as prep_ligand
from easyMD.prep import rcsb
from tempfile import NamedTemporaryFile
from pathlib import Path

import logging

logger = logging.getLogger(__name__)

def load_forcefield_residue_names(forcefield: openmm.app.ForceField):
    '''
    Given an OpenMM ForceField object, return a list of names of the residues which are templated by default.
    '''
    return forcefield._templates.keys()

def get_non_standard_residues(forcefield_files: list, pdb_path: str):
    '''
    Given a set of forcefield files, and path to a PDB, find all the residues whose names are not in the forcefield files.
    Returns list of residues in the form [(chain id, index, name), ...]
    '''
    # Create the openMM forcefield:
    forcefield = openmm.app.ForceField(*forcefield_files)

    template_names = load_forcefield_residue_names(forcefield)

    # Load the PDB, and make a list of the residue names:
    pdb = openmm.app.PDBFile(pdb_path)
    pdb_residue_list = []
    for residue in pdb.topology.residues():
        pdb_residue_list.append((residue.chain.id, int(residue.id), residue.name))

    unmatched_residues = []
    for pdb_residue in pdb_residue_list:
        if pdb_residue[2] not in template_names:
            unmatched_residues.append(pdb_residue)

    # Sometimes our forcefield will not include common AA names (like HIS) in favor of more specific names (e.g. HID, HIE).
    # We assume these will still be captured by the forcefield, so I'll exclude the 20 'regular' AA names:
    unmatched_residues_without_common = []
    for i, residue in enumerate(unmatched_residues):
        if residue[2] not in common_residues:
            unmatched_residues_without_common.append(residue)
    
    return unmatched_residues_without_common

def add_molecule_to_forcefield(forcefield, molecule, name=None, generator='GAFF', generator_forcefield=None):
    '''
    Adds an OpenFF Molecule to a ForceField file:
    '''
    from openmmforcefields.generators import GAFFTemplateGenerator, SMIRNOFFTemplateGenerator, EspalomaTemplateGenerator

    logger.info(f"Adding {generator} template for the ligand {name}")

    # Choose the appropriate template generator:
    available_template_generators = {
        'GAFF': GAFFTemplateGenerator,
        'SMIRNOFF': SMIRNOFFTemplateGenerator,
        'ESPALOMA': EspalomaTemplateGenerator
    }
    TemplateGenerator = available_template_generators[generator]

    # Create it, with the chosen forcefield if necessary (Each template generator has a default forcefield):
    if generator_forcefield is not None:
        template_generator = TemplateGenerator(molecules=[molecule], forcefield=generator_forcefield)
    else:
        template_generator = TemplateGenerator(molecules=[molecule])
    
    if name is None:
        name = molecule.to_smiles()
    
    # Add the template generator to the forcefield:
    forcefield.registerTemplateGenerator(template_generator.generator)
    return forcefield

def make_forcefield(forcefield_files, NS_residue_dict):
    ''''
    Given a set of forcefield files, a system PDB, and a list of ligand SDFs, return an OpenMM ForceField object.'
    Inputs:
    - forcefield_files: list of strings, paths to forcefield files.
    - NS_residue_dict: dictionary of non-standard residues, in the form {residue_name: (SDF_path, PDB_path)}
    '''
    # Step 1: Create forcefield. 
    logger.info("Setting up forcefield with files: %s." % forcefield_files)
    forcefield = openmm.app.ForceField(*forcefield_files)

    # Step 2: Add ligands to forcefield
    ligand_sdf_paths = [NS_residue_dict[residue][0] for residue in NS_residue_dict]
    for ligand_sdf_path in ligand_sdf_paths:
        ligand_molecule = prep_ligand.load_openff_ligand_from_sdf(ligand_sdf_path, sanitize=False, removeHs=False)
        forcefield = add_molecule_to_forcefield(forcefield, ligand_molecule, name=ligand_sdf_path)

    return forcefield

def add_harmonic_force(system, modeller, atom1, atom2, k_value=10.0, r0_value=2.0):
    """
    Adds a custom harmonic force to pull two specified atoms together.
    
    Parameters:
      system   : an OpenMM System object.
      modeller      : a modeller object from which the topology is used.
      atom1   : string, the first atom in the format "chain_id:residue_id:atom_name".
      atom2   : string, the second atom in the format "chain_id:residue_id:atom_name".
      k_value  : force constant in kilocalories_per_mole/angstrom^2 (default: 10.0).
      r0_value : target distance in angstrom (default: 2.0).
    
    Returns:
      The modified system with the custom force added.
    
    Raises:
      ValueError if one or both atoms cannot be located in the pdb topology.
    """
    logger.info(f"Adding custom Harmonic force between {atom1} and {atom2}, with force constant {k_value} kcal/(mol*A^2) and target distance {r0_value} A.")

    from openmm import CustomBondForce, unit
    from easyMD.prep.structure import get_atom_indices

    # Get the atom indices from the modeller
    atom1_index, atom2_index = get_atom_indices(modeller, [atom1, atom2])

    # Define a harmonic potential: E = 0.5 * k * (r - r0)^2
    customBondForce = CustomBondForce("0.5 * k * (r - r0)^2")
    customBondForce.addPerBondParameter("k")   # Force constant
    customBondForce.addPerBondParameter("r0")  # Target distance

    # Set the parameters with proper unit conversions
    k = k_value * unit.kilocalories_per_mole/unit.angstrom**2
    r0 = r0_value * unit.angstrom

    # Add the bond between the two atoms
    customBondForce.addBond(atom1_index, atom2_index, [k, r0])

    # Add the custom force to the system
    system.addForce(customBondForce)
    
    return system

def add_angle_force(system, modeller, atom1, atom2, atom3, k_value=10.0, theta0_value=109.5):
    """
    Adds a custom harmonic angle force to maintain a specified angle between three atoms.
    
    Parameters:
      system      : an OpenMM System object.
      modeller    : a modeller object from which the topology is used.
      atom1, atom2, atom3 : strings, each specifying an atom in the format "chain_id:residue_id:atom_name".
      k_value     : force constant in kilocalories_per_mole/radian^2 (default: 10.0).
      theta0_value: target angle in degrees (default: 109.5).
    
    Returns:
      The modified system with the custom angle force added.
    
    Raises:
      ValueError if one or more atoms cannot be located in the pdb topology.
    """
    logger.info(f"Adding custom Harmonic Angle force among {atom1}, {atom2}, and {atom3}, "
                f"with force constant {k_value} kcal/(mol*rad^2) and target angle {theta0_value} degrees.")

    from openmm import CustomAngleForce, unit
    from easyMD.prep.structure import get_atom_indices
    import math

    # Get the atom indices from the modeller
    atom1_index, atom2_index, atom3_index = get_atom_indices(modeller, [atom1, atom2, atom3])

    # Define a harmonic angle potential: E = 0.5 * k * (theta - theta0)^2
    customAngleForce = CustomAngleForce("0.5 * k * (theta - theta0)^2")
    customAngleForce.addPerAngleParameter("k")      # Force constant
    customAngleForce.addPerAngleParameter("theta0")   # Target angle (in radians)

    # Convert target angle from degrees to radians and set unit for k
    k = k_value * unit.kilocalories_per_mole / unit.radian**2
    theta0 = math.radians(theta0_value)

    # Add the angle defined by the three atoms
    customAngleForce.addAngle(atom1_index, atom2_index, atom3_index, [k, theta0])

    # Add the custom force to the system
    system.addForce(customAngleForce)

    return system

def add_torsion_force(system, modeller, atom1, atom2, atom3, atom4, k_value=10.0, periodicity=1, phase_value=0.0):
    """
    Adds a custom torsion force using a periodic potential.
    
    Parameters:
      system      : an OpenMM System object.
      modeller    : a modeller object from which the topology is used.
      atom1, atom2, atom3, atom4 : strings, each specifying an atom in the format "chain_id:residue_id:atom_name".
      k_value     : force constant in kilocalories_per_mole (default: 10.0).
      periodicity : periodicity of the torsion potential (default: 1).
      phase_value : phase offset in degrees (default: 0.0).
    
    Returns:
      The modified system with the custom torsion force added.
    
    Raises:
      ValueError if one or more atoms cannot be located in the pdb topology.
    """
    logger.info(f"Adding custom Torsion force among {atom1}, {atom2}, {atom3}, and {atom4}, "
                f"with force constant {k_value} kcal/mol, periodicity {periodicity}, and phase {phase_value} degrees.")

    from openmm import CustomTorsionForce, unit
    from easyMD.prep.structure import get_atom_indices
    import math

    # Get the atom indices from the modeller
    atom1_index, atom2_index, atom3_index, atom4_index = get_atom_indices(modeller, [atom1, atom2, atom3, atom4])

    # Define a periodic torsion potential: E = k * (1 + cos(n*phi - phase))
    customTorsionForce = CustomTorsionForce("k*(1 + cos(n*theta - phase))")
    customTorsionForce.addPerTorsionParameter("k")     # Force constant
    customTorsionForce.addPerTorsionParameter("n")       # Periodicity
    customTorsionForce.addPerTorsionParameter("phase")   # Phase (in radians)

    # Convert phase from degrees to radians and set unit for k
    k = k_value * unit.kilocalories_per_mole
    phase = math.radians(phase_value)

    # Add the torsion defined by the four atoms
    customTorsionForce.addTorsion(atom1_index, atom2_index, atom3_index, atom4_index, [k, periodicity, phase])

    # Add the custom force to the system
    system.addForce(customTorsionForce)

    return system
