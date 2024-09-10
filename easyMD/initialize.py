def init(input_file, experiment_dir, chains=[], ligands=[], ligand_files=[], model_termini=False):
    '''
    Initialize an easyMD run by splitting the input structure into its components.
    
    Parameters
    ----------
    input_file : str
        Path to the input structure file.
    chains : list of str, optional
        List of chain IDs to keep. Default is to keep all chains.
    ligands : list of str, optional
        List of ligand names to split out. Default is to keep all ligands.
    ligand_files : list of str, optional
        List of paths to ligand SDF files to use instead of splitting out ligands. Default is to split out ligands.
    model_termini : bool, optional
        Whether to model the termini. Default is False.
        
    Returns
    -------
    None. Saves the repaired structure as a .cif file in [experiment_dir]/input/processed/. Saves the ligands as .sdf files in [experiment_dir]/input/processed/ligands/.'''
    
    # Set up the input/raw directory:
    from pathlib import Path
    input_raw_dir = Path(experiment_dir) / 'input' / 'raw'
    input_raw_dir.mkdir(parents=True, exist_ok=True)
    
    # Is the input file a PDB code?
    is_pdb_code = False
    if len(input_file) == 4 and input_file.isalnum():
        is_pdb_code = True
        print('Downloading PDB file...')
        # Download the PDB file
        input_file = download_pdb(input_file, input_raw_dir)
    
    # Otherwise, it's a file:
    else:
        #Verify it's a .pdb, .pdbx, or .cif file:
        if not Path(input_file).suffix in ['.pdb', '.cif', '.pdbx']:
            raise ValueError('Input file must be a .pdb, .cif, or .pdbx file. Received input file: {input_file}')
        
        # Copy the PDB file to the input/raw directory:
        from shutil import copy
        input_raw_dir = Path(experiment_dir) / 'input' / 'raw'
        input_raw_dir.mkdir(parents=True, exist_ok=True)
        input_pdb_path = input_raw_dir / Path(input_file).name
        copy(input_file, input_pdb_path)

    #Load the pdb file
    from openmm.app import PDBxFile
    from pdbfixer import PDBFixer
    pdb = PDBFixer(str(input_file))

    ### Extract ligands and display:
    from easyMD.extract_ligand import save_ligand
    from termol import draw
    ligand_dir = Path(experiment_dir) / 'input' / 'processed' / 'ligands'
    ligand_dir.mkdir(parents=True, exist_ok=True)
    non_standard_residues = []
    for residue in pdb.topology.residues():
        # Is it nonstandard?
        if residue.name not in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'DA', 'DC', 'DG', 'DT']:
            if residue.name not in non_standard_residues:
                print(f'Nonstandard residue: {residue.name}')
                non_standard_residues.append(residue.name)
    print(f'Nonstandard residues: {non_standard_residues}')
    save_ligand(pdb, non_standard_residues, Path(experiment_dir) / 'input' / 'processed' / 'ligands' / f'butt.sdf')
    
    # iterate over all ligands in the ligand directory:
    for ligand_file in ligand_dir.iterdir():
        if ligand_file.suffix == '.sdf':
            print(str(ligand_file))
            draw(str(ligand_file), three_d=False)
    
    ### Fix messed up residues
    pdb.findMissingResidues()
    pdb.findNonstandardResidues()

    #Let's not add missing residues that are on the N or C terminus, as these are often inherently disordered regions
    if not model_termini:
        missingResidues_fixed = {}

        for chain in pdb.topology.chains():
            termini_indices = (0, len(list(chain.residues())))
            for (chain_index, residue_index), residues in pdb.missingResidues.items():
                if chain.index == chain_index and residue_index not in termini_indices:
                    missingResidues_fixed[(chain_index, residue_index)] = residues
        
        pdb.missingResidues = missingResidues_fixed
    pdb.replaceNonstandardResidues() 
    pdb.findMissingAtoms()
    pdb.addMissingAtoms()    #this adds both missing atoms and residues 

    ### Delete chains:
    chains_to_delete = []
    if chains:
        for chain in pdb.topology.chains():
            if chain.id not in chains:
                chains_to_delete.append(chain.index)
        pdb.removeChains(chains_to_delete)

    # We'll add missing atoms and bonds to our file, but only if it's from the Protein Data Bank.
    # This is because the PDB uses a well-defined format for bond orders, which we can use to infer the bond orders of missing atoms.
    if is_pdb_code:
        ### Add missing ligand atoms:
        from easyMD.extract_ligand import add_missing_ligand_atoms
        add_missing_ligand_atoms(input_file, pdb, ligands)

        ### Add bond order info for pdbx:
        from easyMD.extract_ligand import set_bond_order_data
        set_bond_order_data(input_file, pdb, ligands)
    else:
        print('Input file is a PDB file. Without bond order information, no attempt will be made to repair the ligands.')

    ### Write out the fixed pdb
    from openmm.app import PDBxFile
    input_processed_dir = Path(experiment_dir) / 'input' / 'processed'
    input_processed_dir.mkdir(parents=True, exist_ok=True)
    fixed_pdb_path = input_processed_dir / (Path(input_file).stem + '_fixed.cif')
    PDBxFile.writeFile(pdb.topology, pdb.positions, open(fixed_pdb_path, 'w'))

    # Write out as sdf:
    from easyMD.extract_ligand import save_ligand
    input_processed_ligands_dir = input_processed_dir / 'ligands'
    input_processed_ligands_dir.mkdir(parents=True, exist_ok=True)
    ligand_path = input_processed_ligands_dir / (Path(input_file).stem + '_ligand.sdf')
    save_ligand(pdb, ligands, ligand_path)

    # If the user has provided their own ligand SDFs, use those instead:
    from shutil import copy
    if ligand_files:
        for ligand_file in ligand_files:
            #Copy to the input/processed/ligands directory
            ligand_file_path = input_processed_ligands_dir / Path(ligand_file).name
            copy(ligand_file, ligand_file_path)


def download_pdb(pdb_id, download_dir):
    import requests
    from pathlib import Path
    url = f'https://files.rcsb.org/download/{pdb_id}.cif'
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(download_dir / f'{pdb_id}.cif', 'w') as file:
            file.write(response.text)
        print(f'Successfully downloaded {pdb_id}.cif')
    else:
        raise ValueError(f'Failed to download {pdb_id}.cif. Please check the PDB ID.')
    
    return download_dir / f'{pdb_id}.cif'
