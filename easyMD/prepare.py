def prep_wrapper(experiment_dir, params_file, custom_func=None):
    
    # Load json file with all the parameters:
    import json
    with open(params_file, 'r') as f:
        params = json.load(f)

    if 'PREP' in params:
        params = params['PREP']

    print('Beginning preparation...')

    # Load custom prep function if provided:
    if custom_func is not None:
        from easyMD.utils import import_function_from_path
        prep_function = import_function_from_path(custom_func, 'prep')
        print('Using custom prep function:', custom_func)
    # Use default prep function:
    else:
        prep_function = prep

    prep_function(experiment_dir, params)

    print('Finished preparation!')
    # Save output
    # Save prepped system to input/prepped dir

def prep(experiment_dir, prep_params):
    '''
    Prepares the input PDB file for simulation. This includes adding hydrogesn and solvating the system.
    Inputs:
    - experiment_dir (str): The experiment directory. The processed PDB file will be saved to [experiment_dir]/input/processed.
    - prep_params (dict): A dictionary containing the following parameters:
        - pH (float): The pH of the system.
        - ion_concentration (float): The ionic concentration of the system.
        - padding (float): The padding to add around the system, in angstroms. Default is 10.
    Outputs:
    - None. Saves the repaired structure as a .cif file in [experiment_dir]/input/processed/. Saves the ligands as .sdf files in [experiment_dir]/input/processed/ligands/.
    '''

    ### Check for required parameters:
    REQUIRED_PARAMS = ['pH', 'ion_concentration_MOLAR', 'padding_ANGSTROM']
    from easyMD.utils import check_required_params
    check_required_params(prep_params, REQUIRED_PARAMS)

    ### Find the processed PDB file:
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

    ### Add hydrogens according to the pH:
    pdb.addMissingHydrogens(prep_params['pH'])

    ### Make list of residues to keep. Here, we're keeping protein, DNA, RNA, and ligand residues:
    from pdbfixer.pdbfixer import proteinResidues, dnaResidues, rnaResidues
    ligand_dir = processed_dir / 'ligands'
    ligand_resnames = [ligand.stem for ligand in ligand_dir.iterdir()]
    residues_to_keep = proteinResidues + dnaResidues + rnaResidues + ligand_resnames

    ### Delete all nonstandard residues except ligands:
    from openmm.app import Modeller
    # Convert to Modeller
    modeller = Modeller(pdb.topology, pdb.positions)
    # Delete all nonstandard residues except ligands:
    to_delete = []
    for residue in modeller.topology.residues():
        if residue.name not in residues_to_keep:
            to_delete.append(residue)
    modeller.delete(to_delete)
    # Update PDBFixer object:
    pdb.topology = modeller.topology
    pdb.positions = modeller.positions

    ### Add Solvent:
    from openmm import unit as openmm_unit
    pdb.addSolvent(padding=prep_params['padding_ANGSTROM']*openmm_unit.angstroms, ionicStrength=prep_params['ion_concentration_MOLAR']*openmm_unit.molar) #TODO: Add more options here for solvent addition!

    ### Write out the Processed pdb
    print('Writing out the processed input pdb...')
    solvated_pdb_path = processed_dir / (str(Path(pdb_path).stem) + '_solvated.pdb')
    from openmm.app import PDBFile
    PDBFile.writeFile(pdb.topology, pdb.positions, str(solvated_pdb_path))

    ### Output:
    print('Done!')
    print('Number of atoms:', pdb.topology.getNumAtoms())