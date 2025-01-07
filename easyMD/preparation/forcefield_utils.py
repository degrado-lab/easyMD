from openmm.app import ForceField
from openmm.app import PDBFile
from .data import common_residues
import logging

logger = logging.getLogger(__name__)

def load_forcefield_residue_names(forcefield: ForceField):
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
    forcefield = ForceField(*forcefield_files)

    template_names = load_forcefield_residue_names(forcefield)

    # Load the PDB, and make a list of the residue names:
    pdb = PDBFile(pdb_path)
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

def create_forcefield(forcefield_files):
    """
    Creates an OpenMM ForceField object with AMBER for proteins and optionally GAFF for ligand.
    """
    logger.info("Setting up forcefield with files: %s." % forcefield_files)
    forcefield = ForceField(*forcefield_files)
    return forcefield

def add_molecule_to_forcefield(forcefield, molecule):
    '''
    Adds an OpenFF Molecule to a ForceField file:
    '''
    from openmmforcefields.generators import GAFFTemplateGenerator

    logger.info(f"Adding GAFF template for the ligand {molecule.to_smiles()}")
    gaff = GAFFTemplateGenerator(molecules=[molecule])
    forcefield.registerTemplateGenerator(gaff.generator)
    return forcefield