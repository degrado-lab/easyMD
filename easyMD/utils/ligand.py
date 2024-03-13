import openmm
import openmm.unit as openmm_unit

def create_ligand_topology(pdb, ligand_residue_idx, ligand_chain_idx):
    '''
    Create a new openMM topology object for a ligand from a PDBFixer object.
    pdb: PDBFixer object
    ligand_residue_idx_chain: tuple of (residue index, chain index)
    
    '''
    ligand = None
    ligand_name = None
    for residue in pdb.topology.residues():
        if int(residue.id) == ligand_residue_idx and int(residue.chain.index) == ligand_chain_idx:   #note: residue.index is NOT the same as residue.id!
            ligand = residue
            ligand_name = residue.name
            break
    if ligand is None:
        raise ValueError(f'Could not find ligand with residue index {ligand_residue_idx} and chain index {ligand_chain_idx} in the PDBFixer object.')
    print(f'Found ligand: {ligand.name} with {len(list(ligand.atoms()))} atoms')

    #create a new topology object:
    ligand_topology = openmm.app.Topology()

    #add all atoms to a new residue and chain:
    ligand_chain = ligand_topology.addChain()
    ligand_residue = ligand_topology.addResidue(ligand.name, ligand_chain)
    [ligand_topology.addAtom(atom.name, atom.element, ligand_residue) for atom in ligand.atoms()]

    #get positions from ligand.atoms():
    ligand_positions = [pdb.positions[atom.index] for atom in ligand.atoms()]
    ligand_positions = [pos.value_in_unit(openmm_unit.angstroms) for pos in ligand_positions]

    return ligand_topology, ligand_positions, ligand_name

import requests
def download_ligand_sdf(three_letter_code, file_path):
    '''
    Download an SDF file for a ligand from the RCSB PDB website.
    three_letter_code: str, e.g. 'ATP' - the ligand code in the PDB
    file_path: str, path to save the file
    '''

    url = f'https://files.rcsb.org/ligands/download/{three_letter_code}_ideal.sdf'
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(file_path, 'w') as file:
            file.write(response.text)
        print(f"SDF file for {three_letter_code} downloaded successfully.")
    else:
        print(f"Failed to download SDF file. HTTP status code: {response.status_code}")

import openmm
from openmm.app import PDBFile
from openff.toolkit import Molecule
import tempfile
from rdkit import Chem
def extract_and_correct_ligand(pdb, ligand_name, output_path, template_path=None):
    '''
    Extract a ligand from a PDBFixer object and correct its stereochemistry by aligning it to a template molecule.
    This is necessary, because a PDB file doesn't specify bond orders. To use a ligand from a PDB, we need to assign them using some template.
    pdb: PDBFixer object
    ligand_name: str, e.g. 'ATP'
    output_path: str, path to save the corrected ligand SDF file
    template_path:  str, path to an SDF file of the template molecule with correct bond orders. 
                    If None, the template will be downloaded from the RCSB PDB website. The ligand_name will be used as the three-letter code.
    '''
    ### EXTRACT LIGAND FROM PDB FILE

    # get ligand residue idx and chain idx from residue name:
    ligand_residue_idx = None
    ligand_chain_idx = None
    for residue in pdb.topology.residues():
        if residue.name == ligand_name:
            ligand_residue_idx = int(residue.id)
            ligand_chain_idx =   int(residue.chain.index)
            break

    #get the ligand topology and positions:
    ligand_topology, ligand_positions, _ = create_ligand_topology(pdb, ligand_residue_idx, ligand_chain_idx)

    #write the ligand to a new PDB file:
    #ligand_path_pdb = 'ligand.pdb'
    ligand_path_pdb = tempfile.NamedTemporaryFile(suffix='.pdb').name
    PDBFile.writeFile(ligand_topology, ligand_positions, str(ligand_path_pdb))

    #Then, load to an RDKit molecule and save as SDF:
    ligand_mol = Chem.MolFromPDBFile(str(ligand_path_pdb), removeHs=False)

    #specify the stereochemistry:
    Chem.AssignStereochemistryFrom3D(ligand_mol)

    #load the template molecule:
    #template_path = 'ligand_template.sdf'
    if template_path is None:
        print('No ligand SDF template provided. Downloading from RCSB PDB...')
        template_path = tempfile.NamedTemporaryFile(suffix='.sdf').name
        download_ligand_sdf(ligand_name, template_path)
    template_mol = Chem.MolFromMolFile(template_path)

    #remove hydrogens:
    template_mol = Chem.RemoveHs(template_mol)
    ligand_mol = Chem.RemoveHs(ligand_mol)

    #align the ligand to the template:
    from rdkit.Chem import AllChem
    mol = AllChem.AssignBondOrdersFromTemplate(template_mol, ligand_mol)

    #save to file:
    #make sure file exists:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_path))
    writer.write(mol)
    writer.close()