"""
Author: Jarek v Dijk, Daniel Rademaker

Takes the scores, IDs, headers and description of all files containing a given string and add them to a new folder.
Note that the file to take it from is hard coded in line 18

Usage:
    python analSubset.py <term>
Example:
    python analSubset.py finger
    python analSubset.py hemoglobin
"""

import sys

select = sys.argv[1]

allprot = open ("./analysis_results/simpleDSSPv2_all_desc.txt").read().split('\n')
b = open("./resultSubset/{}.txt".format(select),'w') 
selectedprot = [b.write(prot+'\n') for prot in allprot if select in prot.lower()]
b.close()
