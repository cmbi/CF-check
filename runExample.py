from preProcessClass import preProcessClass as ppc
from countClass import countClass
from scoreClass import scoreClass
from mkname import createPdbs
from changeIDs import changeID
from compare import compare
from linking import link2Desc


####s
## codeblock 1
####
preproc = ppc(outputFileName = 'preprocess_results/simpleDSSPv2_all', debugFilename='preprocess_results/debug', DatabaseName = "data/PDBfinder2.txt",  changeLetters = True, wantedFile = "data/wanted.txt")

preproc.makeTextWithWantedInfo()

####
## codeblock 2
####
N=2
c = countClass(N,"./preprocess_results/simpleDSSPv2_all.txt", "countresults/simpleDSSPv2.txt", skipSymbols = True)

expectedVals, observedVals, scoreMatrix, scoreDict = c.createPrefParams()

####
## codeblock 3
####
s = scoreClass(N,scoreDict,"./preprocess_results/simpleDSSPv2_all.txt","./analysis_results/simpleDSSPv2.txt",removeDAA= False)
s.makeHistogram()
s.saveScores("analysis_results/simpleDSSPv2_all.txt")

####
## codeblock 4 + 5
####
protGroup = "hemoglobin"
createPdbs(protGroup)
newIDlist = "pdbs/pdbs_{}.txt".format(protGroup)

changeID("./preprocess_results/simpleDSSPv2_all.txt", newIDlist, "./preprocess_results/simpleDSSPv2_{}.txt".format(protGroup))

s.makeHistogram(highligh_IdsFile = newIDlist)

####
## codeblock 6 + 7
####

c = countClass(N, 
	"./preprocess_results/simpleDSSPv2_{}.txt".format(protGroup), 
	"countresults/simpleDSSPv2_{}.txt".format(protGroup,
	skipSymbols = True))
expectedVals, observedVals, scoreMatrix, scoreDict = c.createPrefParams()

s = scoreClass(N, scoreDict, "./preprocess_results/simpleDSSPv2_all.txt",
	"./analysis_results/{}_{}.txt".format("simpleDSSPv2", protGroup),
	removeDAA= True
)


s.saveScores("analysis_results/simpleDSSPv2_{}.txt".format(protGroup))
s.makeHistogram(newIDlist)

####
## codeblock 8 + 9
####


link2Desc("./analysis_results/simpleDSSPv2_{}.txt".format(protGroup), "./analysis_results/simpleDSSPv2_{}_desc.txt".format(protGroup)
)


compare("./analysis_results/simpleDSSPv2_all.txt",
	 		protGroup,
	 		IDlist = './pdbs/pdbs_{}.txt'.format(protGroup)
)


