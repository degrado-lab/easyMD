def add_missing_ligand_atoms(input_file, pdb, ligand_residues = [], repair_broken_bonds = True):
    '''
    Get bond order data from a PDBx file. Add bonds to the PDBx topology.
    Parameters:
    input_file (str): Path to the input PDBx file.
    pdb (PDBxFile): PDBxFile object

    Returns:
    None. Modified the PDBx topology in place.
    '''

    ### EXTRACT BOND DATA ###
    from openmm.app.internal.pdbx.reader.PdbxReader import PdbxReader
    from openmm.unit.quantity import Quantity
    from openmm.unit import angstroms
    from openmm.app import element
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
        bond_dict[resname].append((row[atom_id_1], row[atom_id_2], pdbx_bond_to_num(row[value_order])))
    

    ### ADD MISSING BONDS ###
    for residue in pdb.topology.residues():
        if residue.name in ligand_residues:
            # For each ligand residue, we first identify the atoms that are in the bond list but not in the PDBx topology.

            if residue.name not in bond_dict:
                print(f'No bond order data for residue {residue.name}')
                continue
            atom_dict = {atom.name: atom for atom in residue.atoms()}

            # Find the atoms in the bond list which are not in the PDBx topology:
            all_atoms_in_bonds = set([bond[0] for bond in bond_dict[residue.name]]).union(set([bond[1] for bond in bond_dict[residue.name]]))
            all_atoms_in_topology = set([atom.name for atom in residue.atoms()])

            # This has all the atoms in the bond list that are not in the PDBx topology
            unmatched_bond_atoms = all_atoms_in_bonds - all_atoms_in_topology

            # This has all the atoms in the PDBx topology that are not in the bond list
            unmatched_topology_atoms = all_atoms_in_topology - all_atoms_in_bonds
            
            # If there are atoms in the topology that have no bonds, let's throw an error:
            if len(unmatched_topology_atoms) > 0:
                print(f'Atoms in the PDBx topology that have no bonds: {unmatched_topology_atoms}')
                raise ValueError('Atoms in the PDBx topology that have no bonds')
            
            # For each unmatched atom in the bond list, are they only H's? If so, we can add them here:
            # Go through them all first:
            trying_to_add_heavy_atoms = False
            for atom in unmatched_bond_atoms:
                if get_element_of_atom_name(atom)[0] != 'H':
                    print(f'Atom {atom} is not a hydrogen atom. Cannot add it to the topology')
                    trying_to_add_heavy_atoms = True
            if trying_to_add_heavy_atoms:
                raise ValueError('Trying to add heavy atoms to the topology. This is not supported.')
            
            print(f"Topology is missing {len(unmatched_bond_atoms)} atoms")
            for atom in unmatched_bond_atoms:
                if get_element_of_atom_name(atom) == 'H':
                    print(f'Adding hydrogen atom {atom} to residue {residue.name}')

                    #what is the location of the bonding partner? Save as Quantity in partner_pos
                    partner = None
                    for bond in bond_dict[residue.name]:
                        if bond[0] == atom:
                            partner = bond[1]
                        elif bond[1] == atom:
                            partner = bond[0]
                        if partner:
                            partner_obj = atom_dict[partner]
                            partner_pos = pdb.positions[partner_obj.index]
                            break

                    #Let's come up with a tentative location for the hydrogen atom:
                    location = partner_pos.value_in_unit(angstroms) + get_random_offset()
                    quantity = Quantity(location, angstroms)

                    # Now we can add the atom to the topology, and the positions
                    pdb.topology.addAtom(atom, element.hydrogen, residue)
                    pdb.positions.append(quantity)
                else:
                    raise ValueError(f'Atom {atom} is not a hydrogen atom. Cannot add it to the topology')

def pdbx_bond_to_num(bond_order):
    '''
    Convert a PDBx bond order to a number.
    Parameters:
    bond_order (str): PDBx bond order (e.g. 'sing', 'doub', 'trip')
    
    Returns:
    int: Bond order as a number (e.g. 1, 2, 3)
    '''
    bond_order_map = {'sing': 1, 'doub': 2, 'trip': 3}
    return bond_order_map[bond_order]

def get_element_of_atom_name(atom_name):
    '''
    Get the element of an atom name.
    Parameters:
    atom_name (str): Atom name (e.g. 'C1', 'O1', 'N1')
    
    Returns:
    str: Element of the atom name (e.g. 'C', 'O', 'N')
    '''
    element = ''
    for char in atom_name:
        if char.isalpha():
            element += char
        else:
            break
    return element

def set_bond_order_data(input_file, pdb, ligand_residues = []):
    '''
    Get bond order data from a PDBx file. Add bonds to the PDBx topology.
    Parameters:
    input_file (str): Path to the input PDBx file.
    pdb (PDBxFile): PDBxFile object

    Returns:
    None. Modified the PDBx topology in place.
    '''
    ### EXTRACT BOND DATA ###
    from openmm.app.internal.pdbx.reader.PdbxReader import PdbxReader
    from openmm.unit.quantity import Quantity
    from openmm.unit import angstroms
    from openmm.app import element
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
        bond_dict[resname].append((row[atom_id_1], row[atom_id_2], pdbx_bond_to_num(row[value_order])))

    for residue in pdb.topology.residues():
        if residue.name in ligand_residues:
            if residue.name not in bond_dict:
                print(f'No bond order data for residue {residue.name}')
                continue
            atom_dict = {atom.name: atom for atom in residue.atoms()}
            bonds = bond_dict[residue.name]

            for bond in bonds:
                try:
                    pdb.topology.addBond(atom_dict[bond[0]], atom_dict[bond[1]], type=bond[2])
                except KeyError:
                    print(f'Could not add bond between {bond[0]} and {bond[1]}')

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

    # Add bonds to the molecule
    for bond in bonds:
        atom1_idx, atom2_idx, bond_order = bond
        mol.AddBond(atom_indices[atom1_idx], atom_indices[atom2_idx], bond_order)

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

            print('BONDS:', list(residue.bonds()))

            # Before we continue, let's make sure these all have bond info. If not, we'll skip this and do the next ligand:
            all_bonds_have_order = True
            for bond in residue.bonds():
                if bond.type is None: 
                    all_bonds_have_order = False
                    break
            if not all_bonds_have_order:
                print(f'Ligand {residue.name} is missing bond order for bond between {bond[0].name} and {bond[1].name}.')
                print('SDF file must be manually added using the --ligand_files argument.')
                continue
            
            for bond in residue.bonds():
                bond_type_map = {1: Chem.rdchem.BondType.SINGLE, 
                                2: Chem.rdchem.BondType.DOUBLE, 
                                3: Chem.rdchem.BondType.TRIPLE}
                bonds.append((bond[0].index - residue_start_index, bond[1].index - residue_start_index, bond_type_map[bond.type]))
            
            mol = create_molecule(atom_types, atom_names, coordinates, bonds)
            ligand_output_file = output_file.parent / f'{residue.name}.sdf'
            save_molecule_as_sdf(mol, ligand_output_file)

def get_random_offset():
    '''
    Return a Vec3 object with random x, y, and z values in the ranges -0.2 to -0.1, and 0.1 to 0.2.
    '''
    import numpy as np
    from openmm import Vec3
    
    x = np.random.uniform(1, 2) * np.random.choice([-1, 1])
    y = np.random.uniform(1, 2) * np.random.choice([-1, 1])
    z = np.random.uniform(1, 2) * np.random.choice([-1, 1])

    return Vec3(x, y, z)
