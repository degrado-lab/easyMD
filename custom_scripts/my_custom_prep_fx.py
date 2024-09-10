def prep(experiment_dir, prep_params):
    '''
    Prepares the input PDB file for simulation. This includes adding hydrogesn and solvating the system.
    Inputs:
    - experiment_dir (str): The experiment directory. The processed PDB file will be saved to [experiment_dir]/input/processed.
    - prep_params (dict): A dictionary containing the following parameters:
        - pH (float): The pH of the system.
        - ion_concentration (float): The ionic concentration of the system.
        - padding (float): The padding to add around the system, in angstroms. Default is 10.
        - water_model (str): The water model to use.
    Outputs:
    - None. Saves the repaired structure as a .cif file in [experiment_dir]/input/processed/. Saves the ligands as .sdf files in [experiment_dir]/input/processed/ligands/.
    '''

    # Check for required parameters:
    REQUIRED_PARAMS = ['pH', 'ion_concentration', 'padding']
    from easyMD.utils import check_required_params
    check_required_params(prep_params, REQUIRED_PARAMS)

    # Find the file labeled *_fixed.cif or *_fixed.pdb
    from pathlib import Path
    experiment_dir = Path(experiment_dir)
    processed_dir = experiment_dir / 'input' / 'processed'
    pdb_path = None
    for file in processed_dir.iterdir():
        if file.name[-10:] in ['_fixed.cif', '_fixed.pdb']:
            pdb_path = file
            break
    if pdb_path is None:
        raise FileNotFoundError('No repaired structure found in the processed directory. Please run the repair function first. (Must be named *_fixed.cif or *_fixed.pdb)')
    
    from pdbfixer import PDBFixer
    pdb = PDBFixer(str(pdb_path))
    
    ### add hydrogens:
    pdb.addMissingHydrogens(prep_params['pH'])

    ### print the names of all non-standard residues:
    standard_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    nonstandard_residues = []
    for residue in pdb.topology.residues():
        if residue.name not in standard_residues:
            nonstandard_residues.append(residue.name)

    ### MAKE MODELLER OBJECT:
    from openmm.app import Modeller
    modeller = Modeller(pdb.topology, pdb.positions)

    ### Get a list of ligand resnames from the processed/ligands dir:
    ligand_dir = processed_dir / 'ligands'
    ligand_resnames = [ligand.stem for ligand in ligand_dir.iterdir()]

    ### Delete all nonstandard residues except ligand:
    to_delete = []
    for residue in modeller.topology.residues():
        if residue.name in nonstandard_residues and residue.name not in ligand_resnames:
            to_delete.append(residue)
    modeller.delete(to_delete)

    ### Write out the cleaned PDB file:
    from openmm.app import PDBFile
    from pathlib import Path
    cleaned_pdb_path = processed_dir / (str(Path(pdb_path).stem) + '_cleaned.pdb')
    PDBFile.writeFile(modeller.topology, modeller.positions, str(cleaned_pdb_path))

    ### Add Solvent:
    pdb = PDBFixer(str(cleaned_pdb_path))
    from openmm import unit as openmm_unit
    pdb.addSolvent(padding=prep_params['padding']*openmm_unit.angstroms, ionicStrength=prep_params['ion_concentration']*openmm_unit.molar) #TODO: Add more options here for solvent addition!

    ### Write out the Processed pdb
    print('Writing out the processed input pdb...')
    solvated_pdb_path = processed_dir / (str(Path(pdb_path).stem) + '_solvated.pdb')
    PDBFile.writeFile(pdb.topology, pdb.positions, str(solvated_pdb_path))

    ### Output:
    print('Done!')
    print('Number of atoms:', pdb.topology.getNumAtoms())