import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.mlab as mlab
from scipy.stats import norm
import pickle
import heapq


class outlierClass:

    def __init__(self, outlierfile, PreProcessedInputFile, pklInputfile, N, outputfilename):
        '''
        :param outlierfile: file where outliers are saved
        :param PreProcessedInputFile: file with all preprocessed PDB entries
        :param pklInputfile: pickle file of the count results
        :param N: N
        :param outputfilename: name to give output
        '''
        self.outlierfile = outlierfile + ".txt"
        self.PreProcessedInputFile = PreProcessedInputFile + ".txt"
        self.pklInputfile = pklInputfile + ".pkl"
        self.N = N
        self.outputfilename = outputfilename + ".txt"

    def parser(self):
        '''
        goes through each outlier and gets the IDs
        :return: IDlist
        '''
        infile = open(self.outlierfile, 'r').read()
        IDlist = []
        for line in infile.split('\n'):
            x = 0
            for tab in line.split('\t'):
                x += 1
                if x == 2:
                    IDlist.append(tab)
                else:
                    continue
        return IDlist

    def getOutlierInfo(self):
        '''
        gets the wanted information for each entry
        :return: ID, sequence and DSSP
        '''
        IDlist = self.parser()
        infile = open(self.PreProcessedInputFile, 'r').read()
        entries = infile.split('\n')
        IDIndex = [i for i, s in enumerate(entries[0].split('\t')) if "ID" in s]
        chainIndex = [i for i, s in enumerate(entries[0].split('\t')) if "Chain" in s]
        seqIndex = [i for i, s in enumerate(entries[0].split('\t')) if "Sequence" in s]
        DSSPIndex = [i for i, s in enumerate(entries[0].split('\t')) if "DSSP" in s]
        AAlist = []
        DSSPlist = []
        IDchainList = []
        for entry in entries:
            if entry.split('\t') == [""]:
                continue
            elif entry.split('\t')[IDIndex[0]] + entry.split('\t')[chainIndex[0]] in IDlist:
                IDchainList.append(entry.split('\t')[IDIndex[0]] + entry.split('\t')[chainIndex[0]])
                AAlist.append(entry.split('\t')[seqIndex[0]])
                DSSPlist.append(entry.split('\t')[DSSPIndex[0]])
            else:
                continue
        return IDchainList, AAlist, DSSPlist

    def calcScores(self):
        '''
        calculates the local scores for each entry
        :return:
        '''
        df = pd.read_pickle(self.pklInputfile)
        diction = df.to_dict()
        IDchainList, AAlist, DSSPlist = self.getOutlierInfo()
        scoresArray = []
        for entry in zip(IDchainList, AAlist, DSSPlist):
            scoresArray.append(self.calcScoresPerID(entry, diction))
        return scoresArray, IDchainList, AAlist, DSSPlist

    def calcScoresPerID(self, entry, diction):
        '''
        :param entry: an entry with ID, sequence and DSSP
        :param diction: a dictionary with the preference parameter per AA-DSSP combination
        :return: the score per position
        '''
        AAseq, DSSPseq = list(entry[1]), list(entry[2])
        seqLength = len(AAseq) - self.N + 1
        scoreList = []
        residueList = []
        for i in range(seqLength):
            DSSPCombi = ''.join(DSSPseq[i: i + self.N])
            AACombi = ''.join(AAseq[i: i + self.N])
            scoreList.append(diction[DSSPCombi][AACombi])

        return scoreList

    def calcTreshold(self, scoresArray):
        '''
        :param scoresArray: the scores of all entries
        :return: calculates a treshold to see which positions to highlight
        '''
        highestLowScores = []
        tresholdPercentage = 0.1
        for IDScores in scoresArray:
            highestLowScores.append(max(heapq.nsmallest(int(tresholdPercentage * len(IDScores)), IDScores)))
        return np.mean(highestLowScores)

    def findLowScores(self):
        '''
        uses the treshold and scores to find the most low scores (the most interesting positions)
        :return:
        '''
        scoresArray, IDchainList, AAlist, DSSPlist = self.calcScores()
        treshold = self.calcTreshold(scoresArray)
        outfile = open(self.outputfilename, 'w')
        entry = 0
        outfile.write("Found treshold: {}\n".format(treshold))
        for scoreList in scoresArray:
            combi = 0
            outfile.write("{}\n".format(IDchainList[entry]))
            for score in scoreList:
                combi += 1
                if score <= treshold:
                    # Schrijf de score met het residue nummer op in de outfile
                    combiMax = combi + self.N - 1
                    seq = AAlist[entry][combi - 1:combiMax]
                    DSSP = DSSPlist[entry][combi - 1:combiMax]
                    outfile.write("Residues {} to {}\t {}\t{}\t{}\n".format(combi, combi + self.N, seq, DSSP, score))
                else:
                    continue
            entry += 1
        outfile.close()


#####testcode:
if __name__ == '__main__':
    outliers1 = outlierClass("analysis_results/outliers_N1", "preprocess_results/PreProPDBfind_changedDSSP_RK231",
                        "countresults/changedPDB_RK365_N1_SKPD=False", 1, "analysis_results/outlierLocations_N1")

    outliers2 = outlierClass("analysis_results/outliers_N2", "preprocess_results/PreProPDBfind_changedDSSP_RK231",
                        "countresults/changedPDB_RK887_N2_SKPD=False", 2, "analysis_results/outlierLocations_N2")

    outliers1.findLowScores()
    outliers2.findLowScores()