"""
Author: Jarek v Dijk, Daniel Rademaker

Script with function to link the entries of scores trained on different FC parameters

"""

def link2Desc(analysis_file, newname):
	"""
	
	:param analysis_file: file with the analysis results
	:param newname: new filename with analysis results and descriptions
	"""
	newEntries = open(newname, 'w')
	entries = open(analysis_file, 'r').read().split('\n')[:-1]
	entries = [e.split('\t') for e in entries]

	descriptions = open('described.txt', 'r').read().split('\n')[:-1]
	descriptions = [d.split('\t') for d in descriptions]
	
	descDict = {}
	for d in descriptions:
		descDict[d[0]] = d[2]
		
	for e in entries:
		desc = descDict[e[1]]
		for i in e:
			newEntries.write(str(i) + '\t')
		newEntries.write(desc + '\n')