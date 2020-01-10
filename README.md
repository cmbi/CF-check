# dependency list
all is made in ```Python 2.7.17```


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

Example with N=2 and simplifying the secondary structures:

```
mkdir preprocess_results
python
from preProcessClass import preProcessClass as ppc

preproc = ppc(outputFileName = 'preprocess_results/simpleDSSPv2_all', debugFilename='preprocess_results/debug', DatabaseName = "data/PDBfinder2.txt",  changeLetters = True, wantedFile = "data/wanted.txt")

preproc.makeTextWithWantedInfo()
```

Now we have all the non-redundant entries from the PDBfinder. Next we want to make a file with the descriptions for all the entries


## Determine the chou-fastman parameters
We now determine the CF parameters for all protein. 

Here we will create the CF parameters for N=1, skip any non-wanted symbols (hardcoded as ['?', '-']) and include d-amino acids:


```
mkdir countresults
python
from countClass import countClass
N=1
c = countClass(N,"./preprocess_results/simpleDSSPv2_all.txt",
	"countresults/simpleDSSPv2.txt",
	skipSymbols = True)
	
expectedVals, observedVals, scoreMatrix, scoreDict = c.createPrefParams() 
```

## determine protein scores
We now have the preference parameters saved in the scoreDict. In the following example all proteins are scores based on the CF preference parameters:

```
mkdir

python
from scoreClass import scoreClass

s = scoreClass(N,scoreDict,"./preprocess_results/simpleDSSPv2_all.txt","./analysis_results/simpleDSSPv2.txt",removeDAA= False)

s.makeHistogram()

s.saveScores("analysis_results/simpleDSSPv2_all.txt")
```


## training on specific proteins.
We will now repeat the creating of CF parameters but trained on specific proteins.

### get proteins descriptions
To train on specific proteins, we will create a file to link protein IDs to their description. This will create the `description.txt` file (present by default).

*note: it is assumed all proteins are in ```analysis_results/simpleDSSPv2_all.txt```

```
python2 scrapeDescription.py 
```

### use descriptions to find proteins with description
With the description file, the ```createPdbs``` function can make a list of all proteins with a term of interest in their discription. With this, the proteins with that description will be extracted

```
mkdir pdbs
python

from mkname import createPdbs
from changeIDs import changeID

protGroup = "hemoglobin"
createPdbs(protGroup)
newIDlist = "pdbs/pdbs_{}.txt".format(protGroup)

changeID("./preprocess_results/simpleDSSPv2_all.txt", newIDlist, "./preprocess_results/simpleDSSPv2_{}.txt".format(protGroup))
```

### Train on entries and make scores
First we train on the entries with the countclass:

```
python
c = countClass(N, 
	"./preprocess_results/simpleDSSPv2_{}.txt".format(protGroup), 
	"countresults/simpleDSSPv2_{}.txt".format(protGroup,
	skipSymbols = True))
	expectedVals, observedVals, scoreMatrix, scoreDict = c.createPrefParams() 

```

Then we use the obtained CF parameters to score all proteins. We choose to remove de d-amino acids in this example

```
python
s = scoreClass(N, scoreDict, "./preprocess_results/simpleDSSPv2_all.txt",
	"./analysis_results/{}_{}.txt".format("simpleDSSPv2", protGroup),
	removeDAA= True
)
s.saveScores("analysis_results/simpleDSSPv2_{}.txt".format(protGroup))
s.makeHistogram(newIDlist)
```

### analyze results
For convenience it can be usefull to link the results to their descriptions:
```
python

from linking import link2Desc

link2Desc("./analysis_results/simpleDSSPv2_{}.txt".format(protGroup), "./analysis_results/simpleDSSPv2_{}_desc.txt".format(protGroup)
)

```

The final part is to create plots to compare the results trained on all amino acids with the results the specifically trained amino acids.
```
from compare import compare

compare("./analysis_results/simpleDSSPv2_desc.txt",
	 		protGroup,
	 		IDlist = './pdbs/pdbs_{}.txt'.format(protGroup)
)

```
