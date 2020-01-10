import pandas as pd
import numpy as np
from tqdm import tqdm

class countClass:
    def __init__(self,N, allProteinsFile, outputfileName, skipSymbols = True, skipDAA = False):
        """ The initiation method for countClass
        :param N: The number of aa used in counting.
        :param pickleFileName: The name of the file containing the proteins to base the counts on
        :param outputfileName: The name of the outputfiles that are written not followed by a data-type
        :param skipSymbols: Boolian that decides whether the symbols listed in symbolsToSKip are skipped.
        :param RK: random key to add to the output name
        """

        #set global constants
        self.N = N
        self.outputfileName = "{}_N{}_SKPD={}".format(outputfileName ,self.N, skipSymbols)
        self.symbolsToSKip = np.array(['?', '-'])

        if skipDAA:
            self.symbolsToSKip = self.symbolsToSKip.append('X')

        self.skipSymbols = skipSymbols

        # get data to work on
        allProteins = open(allProteinsFile).read().split('\n')[:-1]
        self.allProteins = [p.split('\t')[:-1] for p in allProteins]
        self.columns = self.allProteins[0]
        self.allProteins = self.allProteins[1:]

        # set the place for each field
        if not self.columns == ['ID', 'Chain', 'HEAD', 'Sequence', 'DSSP']:
            raise ValueError ("Column names not as expected")
        self.IDIdx = 0
        self.chainIdx = 1
        self.HEADIdx = 2
        self.SeqIdx = 3
        self.DSSPIdx = 4



    def determineProtCounts(self, Seq, DSSP, AADict, DSSPdict, Counts):
        """ Produces the AA-,DSSP-combi counts for one AA-,DSSP-sequence.
        :param Seq: AA-sequence that is iterated over.
        :param DSSP: DSSP-sequence that is iterated over
        :param AADict: dictionary AA-combis which have already been iterated over.
        :param DSSPDict: dictionary DSSP-combis which have already been iterated over.
        :param Counts: Dictionary with the counts of an AA|DSSP combination. Order is relavant.
        """

        for i in xrange(len(Seq) - self.N + 1):
            SEQpart = ''.join(Seq[i: i + self.N])
            DSSPpart = ''.join(DSSP[i: i + self.N])

            try:
                AADict[SEQpart] += 1
            except KeyError:
                AADict[SEQpart] = 1

            try:
                DSSPdict[DSSPpart] += 1
            except KeyError:
                DSSPdict[DSSPpart] = 1
            try:
                Counts[SEQpart+'|'+DSSPpart] +=1
            except KeyError:
                Counts[SEQpart + '|' + DSSPpart] = 1

    def determineListCounts(self, listOfProteins = None):
        """ Produces a pandas dataframe with the counts per AA-,DSSP-combi
        :param Sequences: List of proteins with [ID, Chain, Header, Seq, DSSP]
        :param: listOfProteins: List with the entries to base the counts on.
        :return: A pandas dataframe with the counts per AA-, DSSP-combi
        """

        if listOfProteins is None:
            listOfProteins = self.allProteins

        AAcounts = {}
        DSSPcounts = {}
        combiCounts = {}

        print("Creating counts for protein entries")
        for entry in tqdm(listOfProteins):
            seq = entry[self.SeqIdx]
            DSSP = entry[self.DSSPIdx]
            self.determineProtCounts(seq,DSSP,AAcounts,DSSPcounts,combiCounts)

        if self.skipSymbols:  # pass the symbol if it occurs in DSSP or AA
            # Note that this method will not delete a DSSP counted in an incorrect Amino acid. The counts will thus be nearly corerct, but deviate slightly. Speed is greatly improved
            print("Deleting unwanted counts")
            for dict in [AAcounts, combiCounts, DSSPcounts]:
                for key in dict.keys():
                    if any([part in self.symbolsToSKip for part in key]): #check if any part of the key appears in list of symbols to skip
                        del dict[key]

        return (AAcounts, combiCounts, DSSPcounts)


    def createPrefParams(self, AAcounts = None, combiCounts = None, DSSPcounts = None):
        """ Produces the expected matrix from the countmatrix obtained from the makeCounts function.
        :return: The expected matrix in pandas dataframe format.
        """

        if any([dic is None for dic in [AAcounts, combiCounts, DSSPcounts]]):
            AAcounts, combiCounts, DSSPcounts = self.determineListCounts()

        expectedVals = {}

        # if (not(sum(AAcounts.values()) == sum(DSSPcounts.values()))
        #     or not(sum(AAcounts.values()) == sum(combiCounts.values()))):
        #     raise ValueError("Counts are not the same!")

        AAsum = float(sum(AAcounts.values()))
        DSSPsum = float(sum(DSSPcounts.values()))
        combiSum = float(sum(combiCounts.values()))
        scoreDict = {}

        expectedVals = np.zeros((len(DSSPcounts.keys()), len(AAcounts.keys())))
        observedVals = np.zeros((len(DSSPcounts.keys()), len(AAcounts.keys())))
        scoreMatrix = np.zeros((len(DSSPcounts.keys()), len(AAcounts.keys())))

        for i in xrange(len(AAcounts.keys())):
            for k in xrange(len(DSSPcounts.keys())):
                AA = AAcounts.keys()[i]
                DSSP = DSSPcounts.keys()[k]

                expectedVals[k,i] = (AAcounts[AA] / AAsum) * (DSSPcounts[DSSP] / DSSPsum)

                try:
                    observedVals[k,i] = combiCounts[AA + "|" + DSSP] / combiSum
                except:
                    "nothing"

                scoreMatrix[k,i] = np.log(np.divide(observedVals[k,i], expectedVals[k,i])).astype(np.float16)

                if not np.isinf(scoreMatrix[k,i]):
                    scoreDict[AA + "|" + DSSP] = scoreMatrix[k,i]


        scoreMatrix[np.isinf(scoreMatrix)] = -999  # if a value does not occur in the observed, it will be NaN in the divide step. This is very rare.


        pd.DataFrame(expectedVals, columns=AAcounts.keys(), index= DSSPcounts.keys()).to_csv(self.outputfileName+"_expected.csv")
        pd.DataFrame(observedVals, columns=AAcounts.keys(), index=DSSPcounts.keys()).to_csv(
            self.outputfileName + "_observed.csv")
        pd.DataFrame(scoreMatrix, columns=AAcounts.keys(), index=DSSPcounts.keys()).to_csv(
            self.outputfileName + "_scores.csv")

        return (expectedVals, observedVals, scoreMatrix, scoreDict)

# ---- Run code -------------------------------------------------------------------

