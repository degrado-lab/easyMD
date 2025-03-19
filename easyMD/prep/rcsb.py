import logging

logger = logging.getLogger(__name__)

def download_pdb(pdb_id, download_dir):
    import requests
    from pathlib import Path
    url = f'https://files.rcsb.org/download/{pdb_id}.cif'
    response = requests.get(url)
    
    if response.status_code == 200:
        download_dir = Path(download_dir)
        download_dir.mkdir(parents=True, exist_ok=True)
        with open(download_dir / f'{pdb_id}.cif', 'w') as file:
            file.write(response.text)
        logger.info(f'Successfully downloaded {pdb_id}.cif')
    else:
        raise ValueError(f'Failed to download {pdb_id}.cif. Please check the PDB ID.')
    
    return download_dir / f'{pdb_id}.cif'

def download_ligand(ligand_id, download_dir):
    import requests
    from pathlib import Path
    url = f'https://files.rcsb.org/ligands/download/{ligand_id}_ideal.sdf'
    response = requests.get(url)
    
    if response.status_code == 200:
        download_dir = Path(download_dir)
        download_dir.mkdir(parents=True, exist_ok=True)
        with open(download_dir / f'{ligand_id}_ideal.sdf', 'w') as file:
            file.write(response.text)
        logger.info(f'Successfully downloaded {ligand_id}_ideal.sdf')
    else:
        raise ValueError(f'Failed to download {ligand_id}_ideal.sdf. Please check the ligand ID.')
    
    return str(download_dir / f'{ligand_id}_ideal.sdf')