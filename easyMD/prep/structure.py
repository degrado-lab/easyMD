# my_md_package/preparation/structure.py
import logging
import os
from openmm.app import PDBFile, PDBxFile
from pdbfixer import PDBFixer
from .data import common_residues
from rdkit import Chem
from Bio.PDB import PDBParser, PDBIO
from pathlib import Path
logger = logging.getLogger(__name__)


def prepare_protein_pdb(input_path: str, fix: bool, temp_dir: str = '.') -> str:
    """
    - Converts mmCIF to PDB if needed
    - Optionally uses PDBFixer to fix missing atoms
    - Returns final path to a PDB file
    """
    temp_dir = Path(temp_dir)
    if input_path.endswith(".cif"):
        logger.debug("Converting CIF to PDB.")
        pdb_path = temp_dir / "protein.pdb"  # or some temp filename
        cif_to_pdb(input_path, str(pdb_path))
    else:
        pdb_path = input_path

    if fix:
        logger.debug("Fixing PDB with PDBFixer.")
        protein_pdb_fixed_path = temp_dir / "protein_fixed.pdb"
        fix_pdb(str(pdb_path), str(protein_pdb_fixed_path))
        return str(protein_pdb_fixed_path)

    return str(pdb_path)

def cif_to_pdb(cif_file, pdb_file):
    pdbx = PDBxFile(cif_file)
    with open(pdb_file, 'w') as f:
        PDBFile.writeFile(pdbx.topology, pdbx.positions, f, keepIds=True)

def fix_pdb(pdb_file, output_file, add_hydrogens=False):
    fixer = PDBFixer(filename=pdb_file)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    if add_hydrogens:
        fixer.addMissingHydrogens(pH=7.0)
    with open(output_file, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

def add_conect_records(pdb_filename, residues_to_modify, output_filename):
    """
    Add CONECT records for specified residues in a PDB file.
    
    Parameters
    ----------
    pdb_filename : str
        Path to the input PDB file.
    residues_to_modify : list of tuples
        A list of (chain_id, residue_number) tuples specifying which residues
        should have their bonds re-inferred and CONECT records added.
    output_filename : str
        Path to the output PDB file.
    """
    
    # Read the PDB file lines
    with open(pdb_filename, 'r') as f:
        pdb_lines = f.readlines()
    
    # Extract ATOM and HETATM lines and store their parsed info
    # We'll store: index in file, original line, chain, resnum, atom serial number
    atom_records = []
    for i, line in enumerate(pdb_lines):
        record_type = line[0:6].strip()
        if record_type in ("ATOM", "HETATM"):
            atom_serial_str = line[6:11].strip()
            atom_serial = int(atom_serial_str)
            chain_id = line[21].strip()
            resnum_str = line[22:26].strip()
            resnum = int(resnum_str)
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()

            atom_records.append({
                'file_index': i,
                'line': line,
                'record_type': record_type,
                'chain_id': chain_id,
                'resnum': resnum,
                'atom_serial': atom_serial,
                'atom_name': atom_name,
                'res_name': res_name
            })
    
    # Group atoms by (chain, residue number)
    from collections import defaultdict
    residues_atoms = defaultdict(list)
    for atom in atom_records:
        key = (atom['chain_id'], atom['resnum'])
        residues_atoms[key].append(atom)
    
    # Function to create RDKit Mol from a single residue's PDB lines
    def residue_to_mol(res_atom_list):
        # Sort by atom_serial just to be consistent
        res_atom_list = sorted(res_atom_list, key=lambda x: x['atom_serial'])
        # Construct a small PDB block for the residue
        pdb_block = "".join(a['line'] for a in res_atom_list)
        # RDKit expects a trailing newline
        if not pdb_block.endswith('\n'):
            pdb_block += '\n'
        mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False)
        return mol, res_atom_list
    
    # We will store all bonds found: a dict from atom_serial to list of bonded atom_serials
    all_bonds = defaultdict(set)
    
    # Process each specified residue
    for (chain_id, resnum) in residues_to_modify:
        if (chain_id, resnum) not in residues_atoms:
            logger.warning(f"Residue {chain_id} {resnum} not found in PDB, skipping.\n \
                            This may indicate a ligand is not being handled properly.")
            # If the specified residue is not found in the PDB, skip
            continue
        res_atom_list = residues_atoms[(chain_id, resnum)]
        mol, sorted_atom_list = residue_to_mol(res_atom_list)
        if mol is None:
            # If RDKit failed to parse the residue, skip
            continue
        
        # Infer bonds (sometimes RDKit might need a sanitization)
        Chem.SanitizeMol(mol)
        
        # Map RDKit atom index -> original atom_serial
        # since sorted_atom_list is in the order passed to MolFromPDBBlock
        rdkit_to_serial = [a['atom_serial'] for a in sorted_atom_list]
        
        # Loop through bonds in RDKit molecule and record them
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            serial1 = rdkit_to_serial[a1]
            serial2 = rdkit_to_serial[a2]
            # Add to both directions to ensure completeness
            all_bonds[serial1].add(serial2)
            all_bonds[serial2].add(serial1)
    
    # Now we have all bonds for the specified residues.
    # We must append CONECT records at the end of the PDB file.
    # If the PDB contains existing CONECT lines, we can just append after them.
    
    # Find the last ATOM/HETATM line to know where to append
    last_atom_line_idx = -1
    for i, line in enumerate(pdb_lines):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            last_atom_line_idx = i
    
    # Collect all CONECT lines to add
    conect_lines = []
    # For each atom that has bonds, produce CONECT lines
    # PDB format: "CONECTatom_serial bonded_atom_1 bonded_atom_2 ..."
    # Each line can have up to 4 bonds. If more, we create multiple lines.
    for atom_serial, bonded_set in all_bonds.items():
        bonded_list = sorted(bonded_set)
        # Break into chunks of 4
        for i in range(0, len(bonded_list), 4):
            chunk = bonded_list[i:i+4]
            line = ("CONECT%5d" % atom_serial) + "".join("%5d" % b for b in chunk) + "\n"
            conect_lines.append(line)
    
    # Insert the CONECT lines at the end of the file but before END/ENDMDL
    # If there's an END or ENDMDL line, we place CONECT before it.
    end_line_idx = None
    for i, line in enumerate(pdb_lines):
        if line.startswith("END") or line.startswith("ENDMDL"):
            end_line_idx = i
            break
    
    if end_line_idx is None:
        # No END line, just append at the end
        final_lines = pdb_lines + conect_lines
    else:
        # Insert CONECT lines before END line
        final_lines = pdb_lines[:end_line_idx] + conect_lines + pdb_lines[end_line_idx:]
    
    # Write out the new PDB file
    with open(output_filename, 'w') as out:
        out.writelines(final_lines)

def add_conect_for_nonstandard_residues(pdb_file, non_standard_residues, output_file):
    """
    Add CONECT records for non-standard residues in a PDB file.
    
    Parameters
    ----------
    pdb_file : str
        Path to the input PDB file.
    non_standard_residues : list of tuples
        A list of (chain_id, residue_number, name) tuples specifying which residues are nonstandard.
    output_file : str
        Path to the output PDB file.
    """
    # Remove the residue name from the list of non-standard residues
    residues_to_CONECT = [(chain_id, residue_number) for chain_id, residue_number, _ in non_standard_residues]

    # For each, add CONECTs:
    if len(residues_to_CONECT) > 0:
        add_conect_records(pdb_file, residues_to_CONECT, output_file)

def extract_residue_from_pdb(input_pdb: str, chain_id: str, residue_number: int, output_pdb: str):
    """
    Extracts a single residue (specified by chain and residue number) from the given PDB file.
    Retains any CONECT records relevant to that residue and updates atom numbering accordingly.
    
    :param input_pdb: Path to the input PDB file.
    :param chain_id: The chain identifier (e.g. 'A').
    :param residue_number: The residue number (an integer).
    :param output_pdb: Path to the output PDB file.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("input", input_pdb)

    # Locate the residue by chain and residue number.
    # Residue IDs in Biopython are tuples: (hetero_flag, seq_number, insertion_code)
    # Typically (blank, residue_number, ' ') for standard residues.
    target_residue = None
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for res in chain:
                    # Residue ID looks like (" ", residue_number, " ") for standard residues
                    if res.get_id()[1] == residue_number:
                        target_residue = res
                        break
    if target_residue is None:
        raise ValueError(f"Residue {residue_number} not found in chain {chain_id} of {input_pdb}")

    # Extract just that residue into a new structure
    # We'll create a new structure -> model -> chain -> residue hierarchy
    from Bio.PDB.StructureBuilder import StructureBuilder
    sb = StructureBuilder()
    sb.init_structure("extracted")
    sb.init_model(0)
    sb.init_chain(chain_id)
    sb.init_seg(" ")

    # residue_id of the target residue
    hetfield, resseq, icode = target_residue.get_id()
    sb.init_residue(target_residue.resname, hetfield, resseq, icode)
    for atom in target_residue:
        sb.init_atom(atom.name, atom.coord, atom.bfactor, atom.occupancy, atom.altloc, atom.fullname, atom.serial_number, element=atom.element)
    extracted_structure = sb.get_structure()

    # Now we have the extracted structure with one residue.
    # Next step: parse the original PDB for CONECT records relevant to the extracted residue's atoms.
    # We know original atom serial numbers from `target_residue` atoms.

    # Collect atom serials in the target residue.
    residue_atom_map = {}
    for atom in target_residue.get_atoms():
        residue_atom_map[atom.serial_number] = atom.name

    original_conect_lines = []
    with open(input_pdb, 'r') as f:
        for line in f:
            if line.startswith("CONECT"):
                # Format: CONECT atom1 atom2 atom3 atom4 ... 
                # Each atom number is in columns 7-11, 12-16, etc. (5-char fields)
                # We'll parse with regex.
                # According to PDB specification: 
                # Columns: 
                #    1-6  Record name "CONECT"
                #    7-11 serial Atom serial number
                #   12-16 serial Atom serial number
                #   17-21 serial Atom serial number
                #   22-26 serial Atom serial number
                parts = line.strip().split()
                # parts[0] = "CONECT", followed by up to 4 atom serial numbers
                atom_nums = [int(x) for x in parts[1:]]
                # We only keep this line if all atom numbers are in residue_atom_map.
                # Actually, we should keep if at least the first atom belongs to our residue and any connected ones also do.
                # Typically, CONECT lines start with the "central" atom, and the others are connected to it.
                # We only retain fully internal connections (i.e., all atoms on this line are part of the same residue).
                if all(a in residue_atom_map for a in atom_nums):
                    original_conect_lines.append(atom_nums)

    # Now we need to renumber the atoms in the extracted residue structure starting from 1.
    # We'll create a map from old serial -> new serial
    old_serials = list(residue_atom_map.keys())
    old_serials_sorted = sorted(old_serials)
    old_to_new = {old: i+1 for i, old in enumerate(old_serials_sorted)}

    # Update the structure with new serial numbers
    # Biopython doesn't have a direct method to renumber atom serials in place easily,
    # so we may need to rebuild the atoms in a similar manner or just set .serial_number.
    for model in extracted_structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.serial_number = old_to_new[atom.serial_number]

    # Update the CONECT records with new numbering
    updated_conect_lines = []
    for line_atoms in original_conect_lines:
        new_atoms = [old_to_new[a] for a in line_atoms]
        # Create a formatted CONECT line
        # PDB format for CONECT:
        # COLUMNS   DATA  TYPE     FIELD         DEFINITION
        #  1 - 6    Record name   "CONECT"
        #  7 -11    Integer       serial        Atom  serial number
        #  12-16    Integer       serial        Serial number of bonded atom
        # Repeat this pattern for up to four bonded atoms.
        # We'll only write what's needed:
        conect_line = "CONECT" + "".join(f"{num:5d}" for num in new_atoms)
        updated_conect_lines.append(conect_line)

    # Write out the extracted structure
    io = PDBIO()
    io.set_structure(extracted_structure)
    io.save(output_pdb)

    # Now we append the updated CONECT lines to the output file
    with open(output_pdb, 'a') as out_f:
        for cl in updated_conect_lines:
            out_f.write(cl + "\n")

    logger.debug(f"Extracted residue {residue_number} from chain {chain_id} saved to {output_pdb}")

def find_ligand_chain_and_resid(pdb_topology, residue_name):
    """
    Given a PDB topology and residue name, returns (chain_id, residue_number).
    Raises ValueError if not found.
    """
    for chain in pdb_topology.chains():
        for residue in chain.residues():
            if residue.name == residue_name:
                return chain.id, int(residue.id)
    raise ValueError(f"Residue {residue_name} not found in topology.")

def get_atom_indices(modeller, atoms):
    """
    Given a Modeller object and an atom specification string, return the indices of the specified atoms.

    Parameters:
        modeller (openmm.app.Modeller): The Modeller object containing the topology.
        atoms (list): The list of atom specification strings in the format "chain_id:residue_id:atom_name".

    Returns:
        list: A list of indices of the specified atoms.
    """
    indices = []

    for atom_spec in atoms:
        chain_id, residue_id, atom_name = atom_spec.split(':')
        
        # We'll keep track of these to make the debug messages better:
        found_chain = False
        found_residue = False
        found_atom_name = False
        for atom in modeller.topology.atoms():
            if (atom.residue.chain.id == chain_id):
                found_chain = True
                if (atom.residue.id.strip() == residue_id):
                    found_residue=True
                    if (atom.name == atom_name):
                        found_atom_name=True
                        # Add the atom index to the list
                        indices.append(atom.index)
                        break
        else:
            error_string = f"Could not find atom {atom_spec} in the topology."
            if not found_chain:
                error_string += f"\nChain {chain_id} not found.\nFound chains: {[chain.id for chain in modeller.topology.chains()]}"
            elif not found_residue:
                error_string += f"\nResidue {residue_id} not found.\nFound residues: {[residue.id for residue in modeller.topology.residues() if residue.chain.id == chain_id]}"
            elif not found_atom_name:
                error_string += f"\nAtom {atom_name} not found.\nFound atoms: {[atom.name for atom in modeller.topology.atoms() if atom.residue.chain.id == chain_id and atom.residue.id.strip() == residue_id]}"
            raise ValueError(error_string)

    return indices