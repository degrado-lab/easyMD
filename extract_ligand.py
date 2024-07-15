def set_bond_order_data(input_file, pdb, ligand_residues = [], repair_broken_bonds = True):
    '''
    Get bond order data from a PDBx file. Add bonds to the PDBx topology.
    Parameters:
    input_file (str): Path to the input PDBx file.
    pdb (PDBxFile): PDBxFile object

    Returns:
    None. Modified the PDBx topology in place.
    '''
    from openmm.app.internal.pdbx.reader.PdbxReader import PdbxReader
    with open(str(input_file)) as f:
        reader = PdbxReader(f)
        data = []
        reader.read(data)
    block = data[0]

    bond_order_data = block.getObj('chem_comp_bond')
    comp_id = bond_order_data.getAttributeIndex('comp_id')
    atom_id_1 = bond_order_data.getAttributeIndex('atom_id_1')
    atom_id_2 = bond_order_data.getAttributeIndex('atom_id_2')
    value_order = bond_order_data.getAttributeIndex('value_order')
    pdbx_aromatic_flag = bond_order_data.getAttributeIndex('pdbx_aromatic_flag')
    pdbx_stereo_config = bond_order_data.getAttributeIndex('pdbx_stereo_config')
    pdbx_ordinal = bond_order_data.getAttributeIndex('pdbx_ordinal')

    bond_dict = {}
    for row in bond_order_data.getRowList():
        resname = row[comp_id]
        if resname not in bond_dict:
            bond_dict[resname] = []
        bond_dict[resname].append((row[comp_id], row[atom_id_1], row[atom_id_2], row[value_order]))


    ### Fixing Atom Mislabels:
    ### In the PDBx/cif file, the atom names in the topology (derived from the RCSB PDB file) are not always the same as the atom names in the bond list (derived from the RCSB SDF file).
    ### We'll attempt to repair this by matching the atom names in the bond list to the atom names in the topology.
    if repair_broken_bonds:
        for residue in pdb.topology.residues():
            if residue.name in ligand_residues:
                if residue.name not in bond_dict:
                    print(f'No bond order data for residue {residue.name}')
                    continue
                atom_dict = {atom.name: atom for atom in residue.atoms()}

                # Find the atoms in the bond list which are not in the PDBx topology:
                all_atoms_in_bonds = set([bond[1] for bond in bond_dict[residue.name]]).union(set([bond[2] for bond in bond_dict[residue.name]]))
                all_atoms_in_topology = set(atom_dict.keys())

                unmatched_bond_atoms = all_atoms_in_bonds - all_atoms_in_topology
                unmatched_topology_atoms = all_atoms_in_topology - all_atoms_in_bonds

                if len(unmatched_bond_atoms) > 0 or len(unmatched_topology_atoms) > 0:
                    print('Mismatch found between atom names in topology and bond list:')
                    print(f'Unmatched bond atoms: {unmatched_bond_atoms}')
                    print(f'Unmatched topology atoms: {unmatched_topology_atoms}')
                    print('Attempting to repair...')

                # For each unmatched atom in the bond list, examine the unmatched topology atoms.
                # Is there a corresponding unmatched topology atom that is the same element, and is within bond distance to the unmatched bond atom's bonding partner?
                # If so, rename the unmatched bond atom to the unmatched topology atom name.
                # If not, print a warning.
                import numpy as np
                from openmm.unit import angstroms
                newly_matched_bond_atoms = []
                
                ### This is a mess. I hope I will clean this up one day.
                for unmatched_bond_atom in unmatched_bond_atoms:
                    for unmatched_topology_atom in unmatched_topology_atoms:

                        #Get the letter-part of the atom names to compare (e.g. 'H011' -> 'H')
                        unmatched_topology_atom_element = ''.join([i for i in unmatched_topology_atom if not i.isdigit()])
                        unmatched_bond_atom_element = ''.join([i for i in unmatched_bond_atom if not i.isdigit()])

                        if unmatched_topology_atom_element == unmatched_bond_atom_element:
                            for bond in bond_dict[residue.name]: # the rare for-else loop!
                                partner = None
                                target_atom_index = None # (1 or 2)
                                if bond[1] == unmatched_bond_atom:
                                    target_atom_index = 1
                                    partner = bond[2]
                                elif bond[2] == unmatched_bond_atom:
                                    target_atom_index = 2
                                    partner = bond[1]
                                if partner:
                                    # Is the partner close to the proposed match?
                                    partner_obj = atom_dict[partner]
                                    partner_pos = pdb.positions[partner_obj.index].value_in_unit(angstroms)
                                    unmatched_topology_atom_obj = atom_dict[unmatched_topology_atom]
                                    unmatched_topology_atom_pos = pdb.positions[unmatched_topology_atom_obj.index].value_in_unit(angstroms)

                                    distance = np.linalg.norm(np.array(partner_pos) - np.array(unmatched_topology_atom_pos))

                                    if distance < 3:
                                        print(f'Matching {unmatched_bond_atom} to {unmatched_topology_atom}')
                                        # Rename the atom in the bond list to the atom in the topology
                                        fixed_bond = list(bond)
                                        fixed_bond[target_atom_index] = unmatched_topology_atom
                                        bond_dict[residue.name].remove(bond)
                                        bond_dict[residue.name].append(tuple(fixed_bond))
                                        newly_matched_bond_atoms.append(unmatched_bond_atom)
                                        break
                            else:
                                print(f'Could not find a bonding partner for {unmatched_bond_atom}')
                                break
                    else:
                        if unmatched_bond_atom not in newly_matched_bond_atoms:
                            print(f'Could not find a matching atom for {unmatched_bond_atom}')
    
    for residue in pdb.topology.residues():
        if residue.name in ligand_residues:
            if residue.name not in bond_dict:
                print(f'No bond order data for residue {residue.name}')
                continue
            atom_dict = {atom.name: atom for atom in residue.atoms()}
            bonds = bond_dict[residue.name]
            for bond in bonds:
                try:
                    pdb.topology.addBond(atom_dict[bond[1]], atom_dict[bond[2]], type=bond[3])
                except KeyError:
                    print(f'Could not add bond between {bond[1]} and {bond[2]}')

def create_molecule(atom_types, atom_names, coordinates, bonds):
    '''
    Create an RDKit molecule from atom types, atom names, coordinates, and bonds.
    Parameters:
    atom_types (list): List of atom types. (e.g. ['C', 'C', 'O', 'N'])
    atom_names (list): List of atom names. (e.g. ['C1', 'C2', 'O1', 'N1'])
    coordinates (list): List of atom coordinates in Angstroms. (e.g. [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (1.5, 1.5, 0.0), (0.0, 1.5, 0.0)])
    bonds (list): List of bonds. Each bond is a tuple of the form (atom1_idx, atom2_idx, bond_order). (e.g. [(0, 1, Chem.rdchem.BondType.SINGLE), (1, 2, Chem.rdchem.BondType.DOUBLE), (1, 3, Chem.rdchem.BondType.SINGLE)])
    '''
    from rdkit import Chem
    mol = Chem.RWMol()
    atom_indices = {}

    #coordinates is currently like: Quantity(value=Vec3(x=-2.0592680000000003, y=-1.2019890000000002, z=-2.565534), unit=nanometer)
    #need to convert to a tuple:
    if type(coordinates[0]) != tuple:
        from openmm.unit import angstroms
        coordinates = [tuple([coord.value_in_unit(angstroms).x, coord.value_in_unit(angstroms).y, coord.value_in_unit(angstroms).z]) for coord in coordinates]
    
    # Add atoms to the molecule
    for i, (atom_type, atom_name, coord) in enumerate(zip(atom_types, atom_names, coordinates)):
        atom = Chem.Atom(atom_type)
        atom.SetProp('atomName', atom_name)
        atom_idx = mol.AddAtom(atom)
        atom_indices[i] = atom_idx

    # Add a conformer to the molecule
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, coord in enumerate(coordinates):
        conf.SetAtomPosition(i, coord)
    mol.AddConformer(conf)

    mol = repair_broken_bonds(mol)

    # Add bonds to the molecule
    for bond in bonds:
        atom1_idx, atom2_idx, bond_order = bond
        mol.AddBond(atom_indices[atom1_idx], atom_indices[atom2_idx], bond_order)

    return mol

def repair_broken_bonds(mol):
    '''
    Repair broken bonds in an RDKit molecule.
    Parameters:
    mol (rdkit.Chem.rdchem.Mol): RDKit molecule.
    '''
    from rdkit import Chem
    Chem.SanitizeMol(mol)
    return mol

def save_molecule_as_sdf(mol, output_file):
    from rdkit import Chem
    w = Chem.SDWriter(str(output_file))
    w.write(mol)
    w.close()

def save_ligand(pdb, ligands, output_file):
    '''
    Save ligand as an SDF file.
    Parameters:
    pdb (PDBxFile): PDBxFile object
    ligands (list): List of ligand residue names (e.g. ['YSB', ...])
    output_file (str): Path to the output SDF file.
    
    Returns:
    None. Saves the ligand as an SDF file.
    '''
    
    from rdkit import Chem
    from pathlib import Path
    output_file = Path(output_file)

    for residue in pdb.topology.residues():
        if residue.name in ligands:
            atom_types = []
            atom_names = []
            coordinates = []
            bonds = []
            residue_start_index = None
            for atom in residue.atoms():
                atom_types.append(atom.element.symbol)
                atom_names.append(atom.name)
                coordinates.append(pdb.positions[atom.index])

                # We need to get the index of the first atom in the residue. 
                # When we split the complex into its components, the index of the first atom in the residue will be 0.
                if residue_start_index is None:
                    residue_start_index = atom.index
                else:
                    residue_start_index = min(residue_start_index, atom.index)
            
            for bond in residue.bonds():
                bond_type_map = {'sing': Chem.rdchem.BondType.SINGLE, 
                                'doub': Chem.rdchem.BondType.DOUBLE, 
                                'trip': Chem.rdchem.BondType.TRIPLE}
                bonds.append((bond[0].index - residue_start_index, bond[1].index - residue_start_index, bond_type_map[bond.type]))
            
            mol = create_molecule(atom_types, atom_names, coordinates, bonds)
            ligand_output_file = output_file.parent / f'{residue.name}.sdf'
            save_molecule_as_sdf(mol, ligand_output_file)