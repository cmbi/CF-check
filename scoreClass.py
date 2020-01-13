import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.mlab as mlab
from scipy.stats import norm
import pickle


class scoreClass:
    def __init__(self, N, scoreDiction, proteinsToScoreFile, outlierFileName, removeDAA = True):
        """ The initiation method of the scoreClass.
        :param pklInputFile: The pkl file produced by the countClass containing the preference parameters.
        :param PreProcessedInputFile: The txt file produced by preProcessClass.
        :param N: The multiplicity used for sequence segments.
        :param outlierFileName: The name of the file to which outliers, later given as raw input, are written.
        """
        # constants

        self.scoreDiction = scoreDiction
        self.N = N
        self.removeDAA = removeDAA
        self.outlierFileName = outlierFileName
        allProteins = open(proteinsToScoreFile).read().split('\n')[:-1]
        self.allProteins = [p.split('\t')[:-1] for p in allProteins]
        self.columns = self.allProteins[0]
        self.allProteins = self.allProteins[1:]

        if not self.columns == ['ID', 'Chain', 'HEAD', 'Sequence', 'DSSP']:
            raise ValueError ("Column names not as expected")
        self.IDIdx = 0
        self.chainIdx = 1
        self.HEADIdx = 2
        self.SeqIdx = 3
        self.DSSPIdx = 4


    def calcScores(self, diction = None):
        """ Calculates scores of all PDB entries and returns their names and scores in two different lists
        :return: A list with all entry names and a list with all entry scores
        """
        if diction is None:
            diction = self.scoreDiction

        entryNames = []
        entryScores = []
        entryHEADS = []
        unused_percentage = []

        for entry in self.allProteins:
            if self.removeDAA and entry[self.SeqIdx].upper().count('X')/float(len(entry[2])) < 0.1:
                entryName, entryScore, entryHEAD, unused_P = self.calcScorePerID(entry, diction)
                entryNames.append(entryName)
                entryScores.append(entryScore)
                entryHEADS.append(entryHEAD)
                unused_percentage.append(unused_P)
            else:
                entryName, entryScore, entryHEAD, unused_P = self.calcScorePerID(entry, diction)
                entryNames.append(entryName)
                entryScores.append(entryScore)
                entryHEADS.append(entryHEAD)
                unused_percentage.append(unused_percentage)

        return entryNames, entryScores, entryHEADS, unused_percentage

    def calcScorePerID(self, entry, diction):
        """ Calculates the score of one single entry and normalizes the score by the legit sequence length.
        :param entry: The entry of which the score is calculated.
        :param diction: The score dictionary created by pkl2Dictionary.
        :return: The score of the single entry.
        """

        entryName = str(entry[self.IDIdx]) + str(entry[self.chainIdx])
        entryHEAD = str(entry[self.HEADIdx])
        AAseq, DSSPseq = list(entry[self.SeqIdx]), list(entry[self.DSSPIdx])
        entryScore = 0
        miscount = 0
        seqLength = len(AAseq) - self.N + 1

        for i in range(seqLength):
                DSSPCombi = ''.join(DSSPseq[i: i + self.N])
                AACombi = ''.join(AAseq[i: i + self.N])

                if self.removeDAA and ('X' in DSSPCombi or 'X' in AACombi):
                    miscount+=1
                else:
                    try:
                        entryScore += diction[AACombi + "|" + DSSPCombi]
                    except:
                        miscount += 1
        unused_percentage = miscount / seqLength
        entryScore = entryScore/seqLength


        return entryName, entryScore, entryHEAD, unused_percentage

    def makeHistogram(self, highligh_IdsFile = None):
        """ Makes a histogram of the calculated scores for each chain and draws a normal distribution.
        :return: Lists of the: entry names, -scores, -headers and found strange percentage. Values for sigma and mu.
        """

        entryNames, entryScores, entryHEADS, unused_percentage = self.calcScores()
        bin2use = len(entryScores) / 10
        print("TOTAL PROTEINS: ", len(entryScores))


        weights = np.ones_like(entryScores) / float(len(entryScores))
        plt.hist(entryScores, bin2use, density=True, color='y', weights=weights)
        mu, sigma = norm.fit(entryScores)

        # plot hemoglobin on its own
        if not highligh_IdsFile is None:
            highligh_Ids = open(highligh_IdsFile, 'r').read().split('\n')[:-1]

            specificScores = [entryScores[i] for i in range(len(entryScores)) if entryNames[i] in highligh_Ids]
            weights = np.ones_like(specificScores) / float(len(specificScores))
            bin2use = len(specificScores) / 10

            plt.hist(specificScores, bin2use, density=True, lw=1, edgecolor='k', alpha = 0.4, color='b', weights=weights)

        fig = plt.gcf()
        plt.grid(False)
        fig.savefig('Histogram.png', dpi=300)
        plt.show()

        str(raw_input('Press enter key to continue'))
        plt.close()

        return entryNames, entryScores, entryHEADS, sigma, mu, unused_percentage

    def saveScores(self, fileName, sort = True):
        """
        :param fileName: filename to save scores to
        :param sort: boolean: sort scores or not
        :return:
        """
        f = open(fileName, 'w')
        entryNames, entryScores, entryHEADS, unused_percentage = self.calcScores()
        if sort:
            order = np.argsort(entryScores)

            entryNames = np.array(entryNames)[order]
            entryScores = np.array(entryScores)[order]
            entryHEADS = np.array(entryHEADS)[order]

        for i in range(len(entryNames)):
            f.write(str(entryScores[i]) + '\t' + str(entryNames[i]) + '\t' + entryHEADS[i] + '\n')
        f.close()


    def findOutliers(self):
        """ Finds the outliers based on given raw input and writes their ID, score and header to the outputfile
        """
        def makeOutputList(indexList, scoreList, nameList, headList, reverseBool=False):
            """ Makes and sorts the outlier lists. This function is only needed within the findOutliers function.
            :param indexList: list of indices which entries are outliers
            :param scoreList: list of scores for all entries
            :param nameList: list of all IDs
            :param headList: list of all headers
            :param reverseBool: boolean to reverse the list
            :return:
            """
            tempScoreList = []
            tempNameList = []
            tempHeadList = []
            for i in indexList:
                tempScoreList.append(scoreList[i])
                tempNameList.append(nameList[i])
                tempHeadList.append(headList[i])

            return [[x, y, z] for x, y, z in
                    sorted(zip(tempScoreList, tempNameList, tempHeadList), reverse=reverseBool)]

        entryNames, entryScores, entryHEADS, sigma, mu, strangeP = self.makeHistogram()

        finish = False
        while not finish:
            right_treshold = -1000 #float(raw_input('Give treshold value for the right-hand side: '))
            left_treshold = -1200 #float(raw_input('Give treshold value for the left-hand side: '))
            strange_treshold = 8 #float(raw_input('Give treshold value for strange percentage: '))/100.
            indexListRightOutliers = [i for i, s in enumerate(entryScores) if s >= right_treshold]
            indexListLeftOutliers = [i for i, s in enumerate(entryScores) if s <= left_treshold]
            indexListStrangeOutliers = [i for i, s in enumerate(strangeP) if s >= strange_treshold]
            noEntries = float(len(entryScores))
            print ("{} right-hand outliers, {} left-hand outliers and {} strange outliers were found".format(
                len(indexListRightOutliers), len(indexListLeftOutliers), len(indexListStrangeOutliers)))
            print("The right-hand represents {} %, left-hand {} % and strange outliers {} %".format
                  (round(float(len(indexListRightOutliers)/noEntries)*100.,2),
                   round(float(len(indexListLeftOutliers)/noEntries)*100.,2),
                   round(float(len(indexListStrangeOutliers)/noEntries)*100.,2)))
            b = 'y' #raw_input('Are these parameters to liking (y/n): ')
            if b == 'y':
                finish = True

        outfile = open(self.outlierFileName, 'w')
        sortedRightOutliers = makeOutputList(indexListRightOutliers, entryScores, entryNames, entryHEADS)
        sortedLeftOutliers = makeOutputList(indexListLeftOutliers, entryScores, entryNames, entryHEADS,
                                            reverseBool=True)
        sortedStrangeOutliers = makeOutputList(indexListStrangeOutliers, strangeP, entryNames, entryHEADS)
        tR = "Sorted outliers with high scores:"
        tL = "\n\nSorted outliers with low scores:"
        tS = "\n\nSorted outliers with high strange percentage:"
        for x, l in zip([sortedRightOutliers, sortedLeftOutliers, sortedStrangeOutliers], [tR, tL, tS]):
            outfile.write(l)
            for i in x:
                outfile.write('\n')
                for j in i:
                    outfile.write(str(j) + '\t')
        outfile.close()

#### testcode
