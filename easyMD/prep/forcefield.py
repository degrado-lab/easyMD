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

def add_molecule_to_forcefield(forcefield, molecule, name=None):
    '''
    Adds an OpenFF Molecule to a ForceField file:
    '''
    from openmmforcefields.generators import GAFFTemplateGenerator
    
    if name is None:
        name = molecule.to_smiles()
    logger.info(f"Adding GAFF template for the ligand {name}")
    gaff = GAFFTemplateGenerator(molecules=[molecule])
    forcefield.registerTemplateGenerator(gaff.generator)
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