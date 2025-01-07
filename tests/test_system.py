# tests/test_system.py
import pytest
from my_md_package.preparation.system import setup_system_from_files

def test_setup_system_minimal_protein(tmp_path):
    protein_pdb = "data/test_protein.pdb"
    # call your function
    sim = setup_system_from_files(protein_pdb, fix=False)
    # check something about sim
    assert sim is not None
    # Finish this...