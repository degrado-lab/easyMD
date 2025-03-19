# my_md_package/preparation/ligand.py
import logging
from openff.toolkit.topology import Molecule
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from openmm.app import PDBFile
from easyMD.prep.data import element_symbols
from tempfile import NamedTemporaryFile
from easyMD.prep import structure as prep_structure
from pathlib import Path
from easyMD.prep import rcsb
import termol

logger = logging.getLogger(__name__)


def load_openff_ligand_from_sdf(sdf_path: str, sanitize: bool = True, removeHs: bool = True) -> Molecule:
    '''
    Load a ligand from an SDF file using the OpenFF Toolkit.
    This allows us to toggle sanitization and hydrogen removal, which can confound kekulization (i.e. representation of aromatic rings as single/double bonds).
    '''
    logger.debug(f"Loading ligand from SDF: {sdf_path}")
    # Here, creating the openFF molecule may fail during kekulization.
    # We toggle sanitization and hydrogen removal, however, in Molecule.from_rdkit() it does this anyway.
    rdkit_mol = Chem.SDMolSupplier(sdf_path, sanitize=sanitize, removeHs=removeHs)[0]
    return Molecule.from_rdkit(rdkit_mol, hydrogens_are_explicit= not removeHs)

def generate_hydrogen_template_xml(ligand_sdf_path, ligand_pdb_path, xml_file, residue_name):
    '''
    Using a system PDB file and a ligand SDF file, generate an XML template for adding hydrogens to the ligand.
    '''
    
    # Here again, we toggle sanitization and hydrogen removal to avoid kekulization issues.
    # We assume this will cause no issues...
    ligand = load_openff_ligand_from_sdf(ligand_sdf_path, sanitize=False, removeHs=False)
    
    if ligand.n_conformers == 0:
        raise ValueError("No 3D conformers found in the SDF. Please provide a 3D conformer.")
    # We assume one molecule, one residue

    def symbol(atom):
        return element_symbols[atom.atomic_number]
    
    # From the PDB file, get a dict mapping the index of each atom to its name
    ### TODO: I've added this on-the-fly PDB extraction. Does it work?
    pdb = PDBFile(ligand_pdb_path)
    atom_names = {atom.index: atom.name for atom in pdb.topology.atoms()}

    # Build a connectivity lookup
    # We'll identify which hydrogen is bound to which heavy atom.
    # In OpenFF Molecule, bonds are stored as pairs of atom indices.(0-indexed)
    parent_index_map = {}  # Maps hydrogen atom index to its parent heavy atom index
    for bond in ligand.bonds:
        a1 = bond.atom1_index
        a2 = bond.atom2_index
        # Identify if one of these is H and the other is heavy
        if symbol(ligand.atoms[a1]) == "H" and symbol(ligand.atoms[a2]) != "H":
            parent_index_map[a1] = a2
        elif symbol(ligand.atoms[a2]) == "H" and symbol(ligand.atoms[a1]) != "H":
            parent_index_map[a2] = a1

    pdb_sdf_map = get_atom_mapping(ligand_pdb_path, ligand_sdf_path, resname=residue_name)
    # Reverse the dict:
    sdf_pdb_map = {v: k for k, v in pdb_sdf_map.items()}

    parent_names_map = {k: atom_names[sdf_pdb_map[v]] for k, v in parent_index_map.items()}

    # Now generate the XML. We only need entries for hydrogens, since the template
    # tells OpenMM how to add hydrogens given the heavy atoms.
    # If hydrogens are already present in the SDF, this template can also serve to confirm their presence.
    xml_lines = []
    xml_lines.append('<Residues>')
    xml_lines.append(f'  <Residue name="{residue_name}">')
    # For each hydrogen atom, write a line
    hydrogen_count = 0
    for i, atom in enumerate(ligand.atoms):
        if symbol(atom) == "H":
            # The parent is the heavy atom to which it's bonded
            parent = parent_names_map[i] if i in parent_names_map else None
            if parent is None:
                raise ValueError(f"Could not find parent heavy atom for hydrogen {atom_names[i]}")
            xml_lines.append(f'    <H name="H{hydrogen_count}" parent="{parent}"/>')
            hydrogen_count += 1
    xml_lines.append('  </Residue>')
    xml_lines.append('</Residues>')

    with open(xml_file, "w") as f:
        f.write("\n".join(xml_lines))

def write_all_hydrogen_templates(NS_residue_dict, temp_dir):
    '''
    Given a system PDB file and a dictionary of non-standard residues mapped to ligand SDFs and PDBs,
    generate XML templates for adding hydrogens to the ligands.
    Inputs:
    - NS_residue_dict: Dictionary mapping residue names to (ligand SDF path, ligand PDB path) tuples.
    - temp_dir: Path to the temporary directory where ligand files are stored.
    Returns:
    - None
    '''
    temp_dir = Path(temp_dir) # Ensure it's a Path object, if not already.
    xml_dir = temp_dir / 'ligands' / 'hydrogen_templates'
    xml_dir.mkdir(parents=True, exist_ok=True)

    hydrogen_templates_dict = {}

    for res_name, (ligand_sdf_path, ligand_pdb_path) in NS_residue_dict.items():
        
        xml_temp_file = NamedTemporaryFile(suffix='.xml', dir=str(xml_dir), delete=False)

        ### TODO: Fix this to match new definition
        generate_hydrogen_template_xml(
            ligand_sdf_path,
            ligand_pdb_path,
            xml_temp_file.name,
            residue_name=res_name
        )

        hydrogen_templates_dict[res_name] = xml_temp_file.name
    
    return hydrogen_templates_dict

def matches_residue_to_sdf(pdb_file, ligand_sdf_path):
    ''' Given a PDB file and a ligand sdf, check if the ligand in the SDF matches the residue in the PDB.
    Returns True if the ligand in the SDF matches the residue in the PDB, False otherwise.
    Ignores hydrogens.
    '''

    # To do this, we'll attempt to create a mapping
    try:
        atom_mapping = get_atom_mapping(pdb_file, ligand_sdf_path, resname="LIG", silent=True)
        return True
    except ValueError as e:
        return False

def get_atom_mapping(pdb_file, sdf_file, resname, silent=False):
    '''
    Given a PDB file and an SDF file, find the Maximum Common Substructure (MCS) between the two molecules and align the SDF molecule to the PDB molecule.
    Then, return a mapping of atom indices between the PDB molecule and the SDF molecule. Indexing is 0-based.
    '''
    pdb_mol = Chem.MolFromPDBFile(str(pdb_file), removeHs=False)
    sdf_suppl = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
    sdf_mol = [m for m in sdf_suppl if m is not None][0]
    
    # Remove hydrogens to focus on heavy atoms for MCS
    pdb_mol_noH = Chem.RemoveHs(pdb_mol)
    sdf_mol_noH = Chem.RemoveHs(sdf_mol)

    # Find MCS with relaxed bond order comparison
    mcs = rdFMCS.FindMCS(
        [pdb_mol_noH, sdf_mol_noH], 
        bondCompare=rdFMCS.BondCompare.CompareAny, 
        ringMatchesRingOnly=True
    )
    if not mcs.smartsString:
        raise ValueError("No MCS found between PDB and SDF molecules.")

    # Create a query molecule from the MCS
    mcs_query = Chem.MolFromSmarts(mcs.smartsString)
    if mcs_query is None:
        raise ValueError("Could not generate a query molecule from MCS.")

    # Match the MCS query against both molecules
    pdb_match = pdb_mol_noH.GetSubstructMatch(mcs_query)
    sdf_match = sdf_mol_noH.GetSubstructMatch(mcs_query)

    if not pdb_match or not sdf_match:
        raise ValueError("MCS does not match on one of the molecules after all.")

    # pdb_match and sdf_match are tuples of atom indices in pdb_mol_noH and sdf_mol_noH respectively
    # We need to map back to the original indices in pdb_mol and sdf_mol.
    pdb_heavy_indices = [i for i, a in enumerate(pdb_mol.GetAtoms()) if a.GetAtomicNum() > 1]
    sdf_heavy_indices = [i for i, a in enumerate(sdf_mol.GetAtoms()) if a.GetAtomicNum() > 1]

    # pdb_match[i] is the index in pdb_mol_noH. The atom at pdb_heavy_indices[pdb_match[i]] is the original pdb atom index.
    # Similarly for sdf_match.
    if len(pdb_match) != len(sdf_match):
        raise ValueError("The length of matched MCS atoms does not align between PDB and SDF.")

    atom_mapping = {}
    for p_idx, s_idx in zip(pdb_match, sdf_match):
        pdb_atom_idx = pdb_heavy_indices[p_idx]
        sdf_atom_idx = sdf_heavy_indices[s_idx]
        atom_mapping[pdb_atom_idx] = sdf_atom_idx

    # Check if the mapping contains all heavy atoms
    num_pdb_heavy_atoms = pdb_mol_noH.GetNumAtoms()
    num_sdf_heavy_atoms = sdf_mol_noH.GetNumAtoms()

    if len(atom_mapping) != num_pdb_heavy_atoms or len(atom_mapping) != num_sdf_heavy_atoms:
        # Which atoms aren't mapped?
        pdb_mapped_atoms = set(atom_mapping.keys())
        pdb_unmapped_atoms = set(range(num_pdb_heavy_atoms)) - pdb_mapped_atoms
        sdf_mapped_atoms = set(atom_mapping.values())
        sdf_unmapped_atoms = set(range(num_sdf_heavy_atoms)) - sdf_mapped_atoms

        message = f'The SDF provided does not contain all heavy atoms from the residue {resname} in the PDB file. \n\
        Unmapped atoms in PDB: {pdb_unmapped_atoms}, Unmapped atoms in SDF: {sdf_unmapped_atoms}'
        
        if not silent:
            termol.draw(Chem.MolToSmiles(pdb_mol), name=pdb_file, three_d=False)
            termol.draw(Chem.MolToSmiles(sdf_mol), name=sdf_file, three_d=False)
        raise ValueError(message)

    # Align the SDF molecule to the PDB molecule using this atom mapping
    AllChem.AlignMol(sdf_mol, pdb_mol, atomMap=list(atom_mapping.items()))
    
    return atom_mapping

def show_non_standard_residues(system_pdb_path, non_standard_residues, width=80, height=40):
    '''
    Given a pdb file and a list of non-standard residues, print out the residue names and atom names.
    Residues are in the format (chain_id, residue_id, residue_name).
    '''
    
    pdb = PDBFile(str(system_pdb_path))
    topology = pdb.topology

    residues_list = []
    for non_standard_residue in non_standard_residues:
        chain_id, residue_id, residue_name = non_standard_residue
        matched_residue = None
        for chain in topology.chains():
            if chain.id == chain_id:
                for residue in chain.residues():
                    if int(residue.id) == residue_id:
                        matched_residue = residue
                        break
        if matched_residue is not None:
            residues_list.append(matched_residue)
    
    # draw with TerMOL:
    logger.info("The following non-standard residues were found:")
    for residue in residues_list:
        try:
            termol.draw(residue_to_smiles(residue, topology), name=str(residue), three_d=False, width=width, height=height)
        except Exception as e:
            print(f"Could not draw residue {residue} due to error: {e}")

def residue_to_smiles(residue, topology):
    """
    Convert a single OpenMM Residue object into a SMILES string using RDKit.

    Parameters
    ----------
    residue : simtk.openmm.app.topology.Residue
        The OpenMM residue object.
    topology : simtk.openmm.app.topology.Topology
        The full OpenMM Topology that contains bond information for the system.

    Returns
    -------
    str
        The (canonical) SMILES string of the residue.
    """

    # Create a writable RDKit molecule
    mol = Chem.RWMol()

    # Keep track of the mapping from OpenMM atom index -> RDKit atom index
    atom_map = {}

    # 1) Add RDKit atoms for each atom in this residue
    for atom in residue.atoms():
        # Convert OpenMM element to atomic number
        atomic_num = atom.element.atomic_number
        # Create a new RDKit atom and add it to the RWMol
        rd_atom = Chem.Atom(atomic_num)
        # Store the RDKit index so we can add bonds later
        rd_idx = mol.AddAtom(rd_atom)
        atom_map[atom.index] = rd_idx

    # 2) Add bonds within this residue
    # We can iterate over all bonds in the entire topology, but only add those
    # that connect atoms belonging to this residue.
    for bond in topology.bonds():
        # bond is a tuple (atom1, atom2)
        atom1, atom2 = bond
        if atom1.residue == residue and atom2.residue == residue:
            # We have a bond within the residue
            rd_idx1 = atom_map[atom1.index]
            rd_idx2 = atom_map[atom2.index]

            # You can also attempt to guess bond orders if needed, but usually
            # the OpenMM topology doesn't store bond orders, just connectivity.
            # For simplicity, we add single bonds:
            mol.AddBond(rd_idx1, rd_idx2, order=Chem.rdchem.BondType.SINGLE)

    # 3) Convert RWMol to Mol and sanitize
    mol = mol.GetMol()
    Chem.SanitizeMol(mol)

    # 4) Generate SMILES (canonical by default)
    smiles = Chem.MolToSmiles(mol, canonical=True)
    return smiles

def prepare_ligand_files(system_pdb_path, non_standard_residues, ligand_sdf_paths, temp_dir):
    ''''
    Verify that the ligand SDFs match the non-standard residues in the system PDB. Extract the ligands to PDB files for later.
    Inputs:
    - system_pdb_path: Path to the system PDB file.
    - non_standard_residues: List of non-standard residues in the format (chain_id, residue_id, residue_name).
    - ligand_sdf_paths: List of paths to ligand SDF files.
    - temp_dir: Path to the temporary directory where ligand files will be stored.
    Returns:
    - NS_residue_dict: Dictionary mapping residue names to (ligand SDF path, ligand PDB path) tuples.
    '''
    
    # Vital for adding ligands to forcefield!
    # Make sure we can find a matching SDF for each non-standard residue.
    # If we have a residue without a matching SDF, we'll need to download it.
    # If that doesn't work, we'll raise an error.

    temp_dir = Path(temp_dir) # Ensure it's a Path object, if not already.
    sdf_dir = temp_dir / 'ligands' / 'SDFs'
    pdb_dir = temp_dir / 'ligands' / 'PDBs'
    sdf_dir.mkdir(parents=True, exist_ok=True)
    pdb_dir.mkdir(parents=True, exist_ok=True)

    NS_residue_dict = {}

    for non_standard_residue in non_standard_residues:
        chain_id, residue_number, res_name = non_standard_residue

        if res_name in NS_residue_dict:
            logger.info("Residue %s already matched to ligand %s." % (res_name, NS_residue_dict[res_name][0]))
            continue
        
        # Extract the residue to a temporary PDB:
        ligand_pdb_path = NamedTemporaryFile(suffix='.pdb', dir=str(pdb_dir), delete=False)
        prep_structure.extract_residue_from_pdb(system_pdb_path, chain_id, residue_number, ligand_pdb_path.name)

        for ligand_sdf_path in ligand_sdf_paths:
            if matches_residue_to_sdf(ligand_pdb_path.name, ligand_sdf_path):
                logger.info("Matched ligand %s to residue %s." % (ligand_sdf_path, res_name))
                NS_residue_dict[res_name] = (ligand_sdf_path, ligand_pdb_path.name)
                break
        else:
            # Attempt to download from PDB:
            try:
                logger.debug("Could not find a match for residue %s in provided SDFs. Attempting to download from PDB." % res_name)
                new_ligand_sdf = rcsb.download_ligand(res_name, sdf_dir)
                # Double check it matches:
                if matches_residue_to_sdf(ligand_pdb_path.name, new_ligand_sdf):
                    logger.info("Downloaded ligand %s from PDB and matched to residue %s." % (new_ligand_sdf, res_name))
                    NS_residue_dict[res_name] = (new_ligand_sdf, ligand_pdb_path.name)
                else:
                    logger.error("Downloaded ligand %s from PDB did not match residue %s." % (new_ligand_sdf, res_name))
            except ValueError as e:
                logger.error(e)
                raise
            ligand_sdf_paths.append(new_ligand_sdf)

    return NS_residue_dict