"""
Author: Jarek v Dijk, Daniel Rademaker

takes random pdb IDs our of the list of all entries and saves it to a file.

Usage:
    python create_random_pdbs.py <amount>
Example:
    python create_random_pdbs.py 150
    python create_random_pdbs.py 10
"""

import numpy as np
import sys


AMOUNT = sys.argv[1]
RK = np.random.randint(0,999,1)

allfiles = open("./preprocess_results/simpleDSSPv2.txt").read().split('\n')[1:-1]
allfiles = [f.split('\t') for f in allfiles]
IDs = np.array([entry[0] + entry[1] for entry in allfiles])

randoms = np.random.randint(0,len(IDs),AMOUNT)
randomIDs = IDs[randoms]

pdbs = open("./pdbs/pdbs_random{}_{}.txt".format(AMOUNT,RK), 'w')

for ID in randomIDs:
	pdbs.write(ID + '\n')
pdbs.close()


