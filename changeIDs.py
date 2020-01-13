"""
Author: Jarek v Dijk, Daniel Rademaker

Script with function to create a file containing only the preprocess results of specific entries
"""


def changeID(preProcessedFile, IDlist, newname):
	"""
	Creates a file with preprocess results of entries in IDlist

	:param preProcessedFile: File with all preprocessed results
	:param IDlist: file with IDs of entries of interest
	:param newname: file where only entries in the IDlist will be stored
	:return:
	"""
	files = open(preProcessedFile, 'r').read().split('\n')[:-1]
	header = files[0]
	files = files[1:]
	IDs = open(IDlist, 'r').read().split('\n')[:-1]
	newfiles = [f for f in files if ''.join(f.split('\t')[:2]) in IDs]
	
	new = open(newname, 'w')
	new.write(header +'\n')
	[new.write(i+'\n')for i in newfiles]
	new.close()
