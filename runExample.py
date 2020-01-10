from preProcessClass import preProcessClass as ppc
from countClass import countClass
from scoreClass import scoreClass
from linking import link2Desc
from compare import compare
from mkname import createPdbs
from changeIDs import changeID

N=2

preproc = ppc(outputFileName = 'preprocess_results/simpleDSSPv2_all', debugFilename='preprocess_results/debug', DatabaseName = "data/PDBfinder2.txt",  changeLetters = True, wantedFile = "data/wanted.txt")
preproc.makeTextWithWantedInfo()

c = countClass(N,"./preprocess_results/simpleDSSPv2_all.txt",
	"countresults/simpleDSSPv2.txt",
	skipSymbols = True)

expectedVals, observedVals, scoreMatrix, scoreDict = c.createPrefParams()




s = scoreClass(N,scoreDict,"./preprocess_results/simpleDSSPv2_all.txt","./analysis_results/simpleDSSPv2.txt",removeDAA= True)
s.makeHistogram()
s.saveScores("analysis_results/simpleDSSPv2_all.txt")



protGroup = "hemoglobin"
createPdbs(protGroup)
newIDlist = "pdbs/pdbs_{}.txt".format(protGroup)

changeID("./preprocess_results/simpleDSSPv2_all.txt", newIDlist, "./preprocess_results/simpleDSSPv2_{}.txt".format(protGroup))




c = countClass(N, 
	"./preprocess_results/simpleDSSPv2_{}.txt".format(protGroup), 
	"countresults/simpleDSSPv2_{}.txt".format(protGroup,
	skipSymbols = True))
expectedVals, observedVals, scoreMatrix, scoreDict = c.createPrefParams()



s = scoreClass(N, scoreDict, "./preprocess_results/simpleDSSPv2_all.txt",
	"./analysis_results/{}_{}.txt".format("simpleDSSPv2", protGroup),
	removeDAA= False
)
s.saveScores("analysis_results/simpleDSSPv2_{}.txt".format(protGroup))
s.makeHistogram(newIDlist)



link2Desc("./analysis_results/simpleDSSPv2_{}.txt".format(protGroup), "./analysis_results/simpleDSSPv2_{}_desc.txt".format(protGroup)
)




compare("./analysis_results/simpleDSSPv2_desc.txt",
	 		protGroup,
	 		IDlist = './pdbs/pdbs_{}.txt'.format(protGroup)
		)
