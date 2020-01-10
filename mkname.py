"""
Author: Jarek v Dijk, Daniel Rademaker

Creates a file containing the IDs matching a given discription.

Usage:
    python mkname.py <term>
Example:
    python mkname.py finger
    python mkname.py hemoglobin
"""

import sys
import pylab as plt
import numpy as np
from scipy.stats import mannwhitneyu, f_oneway

def createPdbs(wat):
	"""
	:param wat: term to look for
	:return: creats a file with IDs matching the term
	"""
	a = open('described.txt').read().split('\n')[:-1]
	a = [i for i in a if wat.lower() in i.lower()]
	a = [i.split('\t') for i in a]
	b = open('./pdbs/pdbs_%s.txt' % (wat), 'w')
	scores = []
	for i in a:
		b.write(i[0] + '\n')
	b.close()
