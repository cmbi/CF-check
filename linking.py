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
	print(descDict)
		
	for e in entries:
		desc = descDict[e[1]]
		for i in e:
			newEntries.write(str(i) + '\t')
		newEntries.write(desc + '\n')

def linkSizes(linkFile):
	entries = open("/mnt/D_disk/My_data/Jarek/Documents/Studie/2018-2019/Protein_stage/Project/BIOPS_BDP_simple/preprocess_results/simpleDSSPv2.txt").read().split("\n")[1:-1]
	entries = [e.split("\t") for e in entries]

	lenDict = {}

	for e in entries:
		lenDict[e[0] + e[1]] = len(e[3])

	entries2LinkTo = open(linkFile, 'r').read().split('\n')[:-1]
	entries2LinkTo = [e.split('\t') for e in entries2LinkTo]
	
	writeFile = open(linkFile, 'w')
	for e in entries2LinkTo:
		writeFile.write("\t".join(e) + "\t" + str(lenDict[e[1]]) + "\n")
	writeFile.close()
