def init(input_file, experiment_dir, chains=[], ligands=[], model_termini=False):
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
    if len(input_file) == 4 and input_file.isalnum():
        print('Downloading PDB file...')
        # Download the PDB file
        input_file = download_pdb(input_file, input_raw_dir)
    
    # Otherwise, it's a file:
    else:
        #Verify it's a .pdb, .pdbx, or .cif file:
        print(input_file)
        if not Path(input_file).suffix in ['.pdb', '.cif', '.pdbx']:
            raise ValueError('Input file must be a .pdb, .cif, or .pdbx file.')
        
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
        print(pdb.missingResidues)
    pdb.replaceNonstandardResidues() 
    pdb.findMissingAtoms()
    pdb.addMissingAtoms()    #this adds both missing atoms and residues 

    ### Delete chains:
    chains_to_delete = []
    if chains:
        for chain in pdb.topology.chains():
            print(chain.id, chains)
            if chain.id not in chains:
                print(f'Deleting chain {chain.id}')
                chains_to_delete.append(chain.index)
        pdb.removeChains(chains_to_delete)

    ### Add bond order info for pdbx:
    from extract_ligand import set_bond_order_data
    set_bond_order_data(input_file, pdb, ligands, repair_broken_bonds=True)

    ### Write out the fixed pdb
    from openmm.app import PDBxFile
    input_processed_dir = Path(experiment_dir) / 'input' / 'processed'
    input_processed_dir.mkdir(parents=True, exist_ok=True)
    fixed_pdb_path = input_processed_dir / (Path(input_file).stem + '_fixed.cif')
    PDBxFile.writeFile(pdb.topology, pdb.positions, open(fixed_pdb_path, 'w'))

    # Write out as sdf:
    from extract_ligand import save_ligand
    input_processed_ligands_dir = input_processed_dir / 'ligands'
    input_processed_ligands_dir.mkdir(parents=True, exist_ok=True)
    ligand_path = input_processed_ligands_dir / (Path(input_file).stem + '_ligand.sdf')
    save_ligand(pdb, ligands, ligand_path)



def download_pdb(pdb_id, download_dir):
    import requests
    from pathlib import Path
    url = f'https://files.rcsb.org/download/{pdb_id}.cif'
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(download_dir / f'{pdb_id}.cif', 'w') as file:
            file.write(response.text)
        print(f'Successfully downloaded {pdb_id}.pdb')
    else:
        raise ValueError(f'Failed to download {pdb_id}.pdb. Please check the PDB ID.')
    
    return download_dir / f'{pdb_id}.cif'


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Initialize an easyMD run.')
    parser.add_argument('input_file', help='Input structure file.')
    parser.add_argument('experiment_dir', help='Directory to save the output files.')
    parser.add_argument('--chains', nargs='+', help='Chains to keep.', default=[])
    parser.add_argument('--ligands', nargs='+', help='Ligands to split out.', default=[])
    parser.add_argument('--model_termini', action='store_true', help='Model the termini.', default=False)
    args = parser.parse_args()

    init(args.input_file, args.experiment_dir, args.chains, args.ligands, args.model_termini)
