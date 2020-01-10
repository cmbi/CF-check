"""
Author: Jarek v Dijk, Daniel Rademaker

Scrapes the discription of the pdb entries and aadds them to analysis resuls.
Note that the analyssi results themselves are outdated, but the ID and discription match will still be usefull

Usage:
    python scapeDescription.py
"""

import urllib2
from tqdm import tqdm

new = open('described.txt', 'w')
pdbs = open("preprocess_results/simpleDSSPv2_all.txt").read().split('\n')[1:-1]
pdbs = [i.split('\t') for i in pdbs]
pdbids = [i[0] for i in pdbs]
pdbs_write = [[i[0]+i[1],i[2]] for i in pdbs]

for i, pdb in enumerate(tqdm(pdbids)):
   gelukt = False
   while not gelukt:
    try:
		data = urllib2.urlopen('https://www.rcsb.org/structure/%s' % pdb).read()
		data = data[data.find('id="structureTitle">')+20:]
		data = data[:data.find('</span>')]

		[new.write(ii+'\t')for ii in pdbs_write[i]]
		new.write(data + '\n')
		gelukt = True
    except:pass
new.close()
	
