# About this program
The Chou Fasman check is made to determine preference parameters of amino acid pairs for secondary structures. This is incorporated in the paper "The future of protein secondary structure prediction was invented
by Oleg Ptitsyn"

# Use this program:
To use this program, we will first create the Chou Fasman (CF) parameters based on all the proteins. We will then do so for specific proteins and compare the results.
_`runexample.py` contains the code given here in a script_

## create folders
For this example we will create some folders to organise the data
```
mkdir preprocess_results countresults pdbs analysis_results
```

## Dependency list
Programs and scripts are develloped in ```Python 2.7.17```

Required packages include:
```
pandas
tqdm
numpy
matlplotlib
scipy
```

## get needed data
First gather the needed data (see readme in the data folder). These are PDBFINDER and the non-redundant proteins list.

## preprocess the PDBfinder to needed entries
The pdbfinder contains much irrelevant information for our purpose. We will only retain the relavant information. All proteins under the size of 50 or not present in the non-redundant list (wanted.txt) are also filtered out by default.

Example with N=2 and simplifying the secondary structures:

```
from preProcessClass import preProcessClass as ppc

preproc = ppc(outputFileName = 'preprocess_results/simpleDSSPv2_all', debugFilename='preprocess_results/debug', DatabaseName = "data/PDBfinder2.txt",  changeLetters = True, wantedFile = "data/wanted.txt")

preproc.makeTextWithWantedInfo()
```

Now we have all the non-redundant entries from the PDBfinder. Next we want to make a file with the descriptions for all the entries


## Determine the chou-fasman parameters
We now determine the CF parameters for all protein. 

Here we will create the CF parameters for N=1, skip any non-wanted symbols (hardcoded as ['?', '-']) and include d-amino acids:


```
from countClass import countClass
N=2
c = countClass(N,"./preprocess_results/simpleDSSPv2_all.txt",
	"countresults/simpleDSSPv2.txt",
	skipSymbols = True)
	
expectedVals, observedVals, scoreMatrix, scoreDict = c.createPrefParams() 
```

## Determine protein scores
We now have the preference parameters saved in the scoreDict. In the following example all proteins are scores based on the CF preference parameters:

```
from scoreClass import scoreClass

s = scoreClass(N,scoreDict,"./preprocess_results/simpleDSSPv2_all.txt","./analysis_results/simpleDSSPv2.txt",removeDAA= False)

s.makeHistogram()
s.saveScores("analysis_results/simpleDSSPv2_all.txt")
```


## Training on specific proteins.
We will now repeat the creating of CF parameters but trained on specific proteins.

### Get proteins descriptions
**only do this if you want to create a new descriptions file, otherwise use the one given (should be the same)**
To train on specific proteins, we will create a file to link protein IDs to their description. This will create the `description.txt` file (present by default).

_note: it is assumed all proteins are in ```analysis_results/simpleDSSPv2_all.txt```_

```
python2 scrapeDescription.py 
```

### Use descriptions to find proteins with description
With the description file, the ```createPdbs``` function can make a list of all proteins with a term of interest in their discription. With this, the proteins with that description will be extracted

```
from mkname import createPdbs
from changeIDs import changeID

protGroup = "hemoglobin"
createPdbs(protGroup)
newIDlist = "pdbs/pdbs_{}.txt".format(protGroup)

changeID("./preprocess_results/simpleDSSPv2_all.txt", newIDlist, "./preprocess_results/simpleDSSPv2_{}.txt".format(protGroup))
```
### see how the selected group is divided over the all-protein trained data
Having the Ids of the protein group, it will be interesting to compare the distribution of these proteins to the distribution of all proteins, trained on all proteins

```
s.makeHistogram(highligh_IdsFile = newIDlist)
```


### Train on entries and make scores
First we train on the entries with the countclass:

```
c = countClass(N, 
	"./preprocess_results/simpleDSSPv2_{}.txt".format(protGroup), 
	"countresults/simpleDSSPv2_{}.txt".format(protGroup,
	skipSymbols = True))
	expectedVals, observedVals, scoreMatrix, scoreDict = c.createPrefParams() 

```

Then we use the obtained CF parameters to score all proteins. We choose to remove de d-amino acids in this example

```
s = scoreClass(N, scoreDict, "./preprocess_results/simpleDSSPv2_all.txt",
	"./analysis_results/{}_{}.txt".format("simpleDSSPv2", protGroup),
	removeDAA= True
)
s.saveScores("analysis_results/simpleDSSPv2_{}.txt".format(protGroup))
s.makeHistogram(newIDlist)
```

### Analyze results
For convenience it can be usefull to link the results to their descriptions:
```
from linking import link2Desc

link2Desc("./analysis_results/simpleDSSPv2_{}.txt".format(protGroup), "./analysis_results/simpleDSSPv2_{}_desc.txt".format(protGroup)
)

```

The final part is to create plots to compare the results trained on all amino acids with the results the specifically trained amino acids.
```
from compare import compare

compare("./analysis_results/simpleDSSPv2_all.txt",
	 		protGroup,
	 		IDlist = './pdbs/pdbs_{}.txt'.format(protGroup)
)

```

# Use program on random proteins
To determine the CF parameters on randomly selected proteins, the above steps can be repeated with one sligt alteration. A list of randomly chosen proteins is used rather than a list based on a search key.
To create a list of (50) random proteins run:
```
python2 create_random_pdbs.py "./preprocess_results/simpleDSSPv2_all.txt" 50
```
```
protGroup = "random50"
newIDlist = "pdbs/pdbs_{}.txt".format(protGroup)
changeID("./preprocess_results/simpleDSSPv2_all.txt", newIDlist, "./preprocess_results/simpleDSSPv2_{}.txt".format(protGroup))
```
Now continue the steps as normal from "Train on entries and make scores" onward.
