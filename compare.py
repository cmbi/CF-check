"""
Author: Jarek v Dijk, Daniel Rademaker

Two functions are made to plot the scores of the proteines trained on all proteins vs the scores based on a specific list of proteins

Usage:
    python compare.py <term>
Example:
    python compare.py finger
    python compare.py hemoglobin
"""

import numpy as np
import pylab as plt
import sys

def procFile(inp):
	"""
	normalized the values for given file
	:param inp: filename with Scores, ID, Head, description
	:return:
	"""
	tmp = open(inp).read().split('\n')[:-1]
	tmp = [i.split('\t')for i in tmp]
	tmp = [[float(i[0])]+i[1:]for i in tmp]
	nmbrs = np.array([i[0]for i in tmp])
	# Normalize to zero mean and unit standard deviation
	nmbrs -= nmbrs.mean()
	nmbrs /= nmbrs.std()

	tmp = [[str(float(nmbrs[i]))]+tmp[i][1:]for i in range(len(tmp))]
	
	return tmp

def link(a, b):
	"""
	links the entries of file a and b together
	:param a: filename of specific protein analysis results with description
	:param b: filename of al proteins analysis results with description
	:return:
	"""
	dic = {}
	for ind in xrange(len(b)):
		dic[b[ind][1]] = ind
	c = [b[dic[a[i][1]]]for i in xrange(len(b))]
	
	return c

def compare(allp, specific, IDlist = None):
	"""
	Makes a plot to compare the scores of all proteins based on all proteins to the scores based on specific protein training.

	:param allp: file with scores of all proteins, trained on all proteins
	:param specific: file with scores, trained on specific proteins
	:param IDlist: filename of Ids to higlight in the comparison plot (usefull to show proteins that were trained on)
	:return:
	"""
	name = specific
	group = 'analysis_results/simpleDSSPv2_{}_desc.txt'.format(name)

	group = procFile(group)
	allp = procFile(allp)
	allp = link(group, allp)
	
	h1 = [float(i[0]) for i in group]
	h2 = [float(i[0]) for i in allp]

	plt.title("specific vs general trained normalized scores")
	plt.plot(h2, h1, 'y+')

	if 1:
		#[plt.plot(h2[i], h1[i], 'ro')for i in range(len(group)) if 'ubiquitin'.lower() in ''.join(group[i]).lower()]
		[plt.plot(h2[i], h1[i], 'ko', ms=3.8)for i in range(len(group)) if 'hemoglobin' in ''.join(group[i]).lower()]
		
		if IDlist is None:	
			[plt.plot(h2[i], h1[i], 'bo', ms=3.8)for i in range(len(group)) if name.lower() in ''.join(group[i]).lower()]
		else:
			IDs = open(IDlist,'r').read().split('\n')[:-1]
			[plt.plot(h2[i], h1[i], 'bo', ms=3.8)for i in range(len(group)) if ''.join(group[i][1]) in IDs]
		
	# [plt.plot(h2[i], h1[i], 'ro', ms=3.8)for i in range(len(group)) if 'antibody'.lower() in ''.join(group[i]).lower()]

	#plt.plot(h2,h2,'r')
	

	plt.xlabel('allprot trained')
	plt.xlabel('allprot trained', fontsize=14, horizontalalignment='right')
	plt.ylabel('{} trained'.format(name), fontsize=14)
	frame1 = plt.gca()
	frame1.axes.xaxis.set_ticklabels([])
	frame1.axes.yaxis.set_ticklabels([])

	plt.savefig('Compare_results1%s.png' % name, dpi=300)
	
	#plt.ylim(min(h1),max(h1))
	#plt.xlim(min(h1),max(h1))
	plt.xlim(-3,3)
	plt.ylim(-3,3)
	plt.plot([min(h1), max(h1)], [min(h1), max(h1)], color='white')
	#plt.xticks[]

	plt.savefig("Compare_result2%s.png" % name, dpi=300)
	plt.show()

if __name__ == "__main__":
	name = sys.argv[1]
	compare(name)
