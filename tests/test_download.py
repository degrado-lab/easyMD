from easyMD.preparation.rcsb import download_pdb

import random

with open('data/rcsb_pdb_ids_10K.txt', 'r') as f:
    pdb_ids = f.read().split(',')
    
random_pdb = random.choice(pdb_ids)
download_pdb(random_pdb, 'downloads')
