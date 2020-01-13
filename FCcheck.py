from preProcessClass import preProcessClass 
from countClass import countClass
from scoreClass import scoreClass
from changeIDs import changeID
import pickle
from compare import compare
from mkname import createPdbs
import sys

def link2Desc(analysis_file, newname):
	newEntries = open(newname, 'w')
	entries = open(analysis_file, 'r').read().split('\n')[:-1]
	entries = [e.split('\t') for e in entries]

	descriptions = open('proc_finished.txt', 'r').read().split('\n')[:-1]
	descriptions = [d.split('\t') for d in descriptions]
	
	descDict = {}
	for d in descriptions:
		descDict[d[1]] = d[4]
		
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

# --- set constants ---
N = 2
allProteinsFile = "simpleDSSPv2" # file name of all proteins files
protGroup = "ubiquitin" # protein group of interest
protGroup = sys.argv[1].strip()
removeX = True


# --- create dataset with headers of wanted group ---
try:
	a = open('./pdbs/pdbs_%s.txt' % (protGroup), 'r')
	a.close()
	print("loaded pdbs")
except Exception as e:
	print(e)
	createPdbs(protGroup)
	
newIDlist = "pdbs/pdbs_{}.txt".format(protGroup)

print(protGroup)

# ---create preprocess_file---
preprocesser = preProcessClass(outputFileName = "./preprocess_results/"+allProteinsFile,
                              debugFilename = "./preprocess_results/"+allProteinsFile,
                              DatabaseName = "./data/PDBfinder2.txt",
                              changeLetters = True, #change DSSP entries to H,E,T
                              wantedFile = "data/wanted.txt"
							  )
#preprocesser.makeTextWithWantedInfo() # Extract the needed info from the pDBfinder


# --- Train on all the data ---
try:
	## Getting back the previous results:
	with open("./countresults/"+allProteinsFile+"N"+str(N)) as f:
		expectedVals, observedVals, scoreMatrix, scoreDict = pickle.load(f)
except: 
	## creating the results:
	c = countClass(N,
	"./preprocess_results/"+allProteinsFile+".txt",
	"countresults/"+allProteinsFile+".txt",
	skipSymbols = True)
	
	expectedVals, observedVals, scoreMatrix, scoreDict = c.createPrefParams() 
	with open("./countresults/"+allProteinsFile+"N"+str(N), 'w') as f: 
		pickle.dump([expectedVals, observedVals, scoreMatrix, scoreDict], f)




# --- score the proteins with trained data ---
## define the score class	
s = scoreClass(N,
			scoreDict, #note that scoredict was obtained from the training
			"./preprocess_results/"+allProteinsFile+".txt",
			"./analysis_results/"+allProteinsFile+".txt",
			removeDAA= removeX
			)
		
## show a histrogram and save all protein scores
s.makeHistogram('./pdbs/pdbs_{}.txt'.format(protGroup))
s.saveScores("analysis_results/"+allProteinsFile+'.txt')




# --- train on the specific group ---
## filter out specific group
changeID("./preprocess_results/"+allProteinsFile+".txt",
		 newIDlist,
		 "./preprocess_results/"+allProteinsFile+"_{}.txt".format(protGroup))
try:
	## Getting back the previous results:
	with open("./countresults/"+allProteinsFile+"_"+protGroup+"N"+str(N)) as f: 
		expectedVals, observedVals, scoreMatrix, scoreDict = pickle.load(f)
except: 
	## Saving the results:
	c = countClass(N, 
		"./preprocess_results/"+allProteinsFile+"_{}.txt".format(protGroup), 
		"countresults/"+allProteinsFile+"_{}.txt".format(protGroup,
		skipSymbols = True))
	expectedVals, observedVals, scoreMatrix, scoreDict = c.createPrefParams() 
	with open("./countresults/"+allProteinsFile+"_"+protGroup+"N"+str(N), 'w') as f:  
		pickle.dump([expectedVals, observedVals, scoreMatrix, scoreDict], f)


#scoreDict['M|E'] = -0.6 #alter scoredict if needed for comparisons

s = scoreClass(N, scoreDict, 
				"./preprocess_results/"+allProteinsFile+".txt",
				"./analysis_results/{}_{}.txt".format(allProteinsFile, protGroup),
				removeDAA= removeX
				)
				
#s.makeHistogram('./pdbs/pdbs_{}.txt'.format(protGroup))
s.saveScores("analysis_results/"+allProteinsFile+'_{}.txt'.format(protGroup))

# --- make a plot for result comparison---
## link reslts to their descriptions
link2Desc("./analysis_results/{}.txt".format(allProteinsFile),
 		  "./analysis_results/{}_desc.txt".format(allProteinsFile)
 	)
link2Desc("./analysis_results/{}_{}.txt".format(allProteinsFile, protGroup),
		  "./analysis_results/{}_{}_desc.txt".format(allProteinsFile, protGroup)
		 )

## make comparison plot

if 1:
	### compare and highlight only protein of interest
	compare("./analysis_results/{}_desc.txt".format(allProteinsFile),
	 		protGroup,
	 		IDlist = './pdbs/pdbs_{}.txt'.format(protGroup)
	 		)
	#compare("./analysis_results/{}_desc.txt".format(allProteinsFile), protGroup)

linkSizes("./analysis_results/{}_{}_desc.txt".format(allProteinsFile, protGroup))




	
	
	
	


