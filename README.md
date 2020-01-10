# dependency lisy
`
pandas
tqdm
numpy
matlplotlib
scipy
`

# Use this program:
To use this program, we will first create the Chou Fastman (FC) parameters based on all the proteins. We will then do so for specific proteins and compare the results.

## get needed data
First gather the needed data (see readme in the data folder)

## preprocess the PDBfinder to needed entries
The pdbfinder contains much irrelevant information for our purpose. We will only retain the relavant information

Example with N=1 and simplifying the secondary structures:

`
mkdir preprocess_results
python
from preProcessClass import preProcessClass as ppc

preproc = ppc(outputFileName = 'preprocess_results/simpleDSSPv2_all', debugFilename='preprocess_results/debug', DatabaseName = "data/PDBfinder2.txt",  changeLetters = True, wantedFile = "data/wanted.txt")

preproc.makeTextWithWantedInfo()
'

Now we have all the non-redundant entries from the PDBfinder. Next we want to make a file with the descriptions for all the entries


## Determine the chou-fastman parameters
We now determine the CF parameters for all protein. 

Here we will create the CF parameters for N=1, skip any non-wanted symbols (hardcoded as ['?', '-']) and include d-amino acids:


`
mkdir countresults
python
from countClass import countClass as cc
N=1
c = cc(N,"./preprocess_results/simpleDSSPv2_all.txt",
	"countresults/simpleDSSPv2.txt",
	skipSymbols = True)
	
expectedVals, observedVals, scoreMatrix, scoreDict = c.createPrefParams() 


## determine protein scores
We now have the preference parameters saved in the scoreDict. In the following example all proteins are scores based on the CF preference parameters:

`
mkdir
from scoreClass import scoreClass as sc
s = sc(N,scoreDict,"./preprocess_results/simpleDSSPv2_all.txt","./analysis_results/simpleDSSPv2.txt",removeDAA= False)
s.makeHistogram()
s.saveScores("analysis_results/simpleDSSPv2_all.txt")
`


## training on specific proteins.
We will now repeat the creating of CF parameters but trained on specific proteins.

### get proteins descriptions
To train on specific proteins, we will create a file to link protein IDs to their description. This will create the `description.txt` file (present by default).

*note: it is assumed all 

`
python2 scrapeDescription.py 

`
