import random
import numpy as np
import pickle
from tqdm import tqdm
class preProcessClass:

    def __init__(self, outputFileName, debugFilename, DatabaseName = "PDBfinder2",  changeLetters = False, wantedFile = "wanted.txt"):
        """
        :param outputFileName: filename to store the preprocessed output. Will automatically convert to .txt
        :param debugFilename: file to store entries that encountered an error
        :param DatabaseName: name with PDB entries to preprocess
        :param changeLetters: Boolean that specifies whether DSSP symbols should be changed to only helices, strands and loops.
        :param wantedFile : File with the wanted, non-redundant entries
        :param RK: a random key added to the end of the filenames. Default is randoms
        """

        self.DatabaseName = DatabaseName
        self.wantedFile = wantedFile
        self.outputFileName = outputFileName
        self.debugFilename = debugFilename
        self.changeLetters = changeLetters
        self.changeWhat = ['B', 'G', 'I', 'S', 'C'] # DSSP letters to change
        self.changeTo = ['E', 'H', 'H', 'T', 'T']   # what to change them to
        self.possibleAAs = ""
        self.possibleDSSPs = ""
        self.posDictonary = {}

    def getEntries(self):

        """ gets the entries from the PDBFINDER2 file
        :return: a 2-D matrix with each entry and all elements of the entry
        """
        text = open(self.DatabaseName, 'r').read()
        entries = text.split('//')
        entries.pop(0)
        entries2 = []
        for entry in entries:
            entries2.append(entry.split('\n')[1:-1])
        return entries2

    def extractWantedInfo(self, IDentry):
        '''
        gets all the wanted info from an entry.
        :param IDentry: list with all information of a single ID in the PDB database
        :return: list with only the wanted information of a single ID. May contain lists as elements
        '''
        infoList = []
        features = []

        # add ID
        infoList.append("ID")
        features.append(self.getIDfromFile(IDentry).split(":")[1].replace(" ", ""))

        # add Head
        infoList.append("HEAD")
        features.append(self.getHeadfromFile(IDentry).split(":")[1].replace(" ", ""))

        sequences = []
        DSSPs = []
        validChains = []

        chainEntries = self.getEntriesPerChain(IDentry)
        for chainEntry in chainEntries:
            try:
                seq, DSSP, chain = self.getValidSeqDSSP(chainEntry)
                if self.changeLetters:
                    DSSP = self.changeLettersOfDSSp(DSSP)
                sequences.append(seq)
                DSSPs.append(DSSP)
                validChains.append(chain)
            except:
                pass
        infoList.append("Chains")
        features.append(validChains)

        infoList.append("Sequences")
        features.append(sequences)

        infoList.append("DSSPs")
        features.append(DSSPs)

        return infoList, features

    def getEntriesPerChain(self, IDentry):
        '''
        splits the entry based on the data for each chain
        :param IDentry: List containing all information for an ID
        :return: list of lists with of all information given, split per chain.
        '''

        chainIndices = [idx for idx, content in enumerate(IDentry) if 'Chain:' in content.replace(' ', '')]
        indexLen = len(chainIndices)
        chainEentryList = []
        if indexLen > 0:
            for i in range(indexLen - 1):
                currentChainIdx = chainIndices[i]
                nextChainIdx = chainIndices[i + 1]
                chainEentryList.append(IDentry[currentChainIdx:nextChainIdx])
            chainEentryList.append(IDentry[chainIndices[(indexLen - 1)]:])
        return chainEentryList

    def getValidSeqDSSP(self, chainEntry):
        """    Function that takes 1 chainEntry as input and returns its sequence, DSSP and chain# if a valid chainEntry is found. Assumes that each chainEntry contains only 1 DSSP and sequence.
        :param chainEntry: chainEntry in the form of a list that is splitted on enters
        :return: Returns the sequence, DSSP and chain if the chainEntry is found valid
        """
        if self.hasDSSP(chainEntry):
            chain = chainEntry[0].split(":")[1].replace(" ", "")
            indexSeq = [i for i, s in enumerate(chainEntry) if 'Sequence' in s][0]
            indexDSSP = [i for i, s in enumerate(chainEntry) if ' DSSP ' in s][0]
            seq = chainEntry[indexSeq].split(":")[1].replace(" ", "")
            DSSP = chainEntry[indexDSSP].split(":")[1].replace(" ", "")
            if self.isValid(seq):
                return seq, DSSP, chain
            else:
                pass
        else:
            pass

    def getInfoOfAllEntries(self):
        '''
        returns a list with all information that is wanted per ID+Chain combination.
        :param filename: file to extract information from
        :return: list with each element being all wanted information of a ID+Chain combination
        '''
        IDentries = self.getEntries()
        Featurenames = ["ID", "Chain", "HEAD", "Sequence", "DSSP"]

        for i,f in enumerate(Featurenames):
            self.posDictonary[f] = i



        Features = []
        for IDentry in IDentries:
            info, feats = self.extractWantedInfo(IDentry)
            ID = feats[info.index("ID")]
            Head = feats[info.index("HEAD")]
            Seqs = feats[info.index("Sequences")]
            DSSPs = feats[info.index("DSSPs")]
            Chains = feats[info.index("Chains")]

            for seq, DSSP, chain in zip(Seqs, DSSPs, Chains):
                Features.append([ID, chain, Head, seq, DSSP])

        return Featurenames, Features

    def makeTextWithWantedInfo(self):
        '''
        calls all other functions to change to PDB file.
        :return: none
        '''
        outfile = open(self.outputFileName+'.txt', 'w')
        debugfile = open(self.debugFilename+'.txt', 'w')
        featureNames, entries = self.getInfoOfAllEntries()

        [outfile.write(name + '\t') for name in featureNames]
        outfile.write('\n')

        text = open(self.wantedFile, 'r').read()
        wantedList = text.split('\n')[:-2]
        for entry in tqdm(entries):
            try:
                if self.isWanted(entry, wantedList):
                    self.writeEntryToOutfile(entry, outfile)
                    self.updateAAandDSSP(entry)
            except:
                for e in entry:
                    debugfile.write(e)
                    debugfile.write('\n')
                debugfile.write("\\")
                debugfile.write('\n')
        outfile.close()
        debugfile.close()
        pickleName = self.outputFileName+'.pkl'
        pickle.dump(self, open(pickleName, 'w'))


    def updateAAandDSSP(self, entry):
        """
        :param entry: adds the amino acids and DSSPs to the list of encountered amino acids and dssps
        :return:
        """
        Seq = self.getSeq(entry)
        DSSP = self.getDSSP(entry)
        self.possibleAAs = "".join(set(self.possibleAAs.join(Seq)))
        self.possibleDSSPs = "".join(set(self.possibleDSSPs.join(DSSP)))


    # ---- helper functions --------------------------------------------------------

    def getID(self, entry):
        '''
        :param entry: entry with an ID on position 0
        :return: the ID of an entry
        '''
        return entry[self.posDictonary["ID"]]

    def getSeq(self,entry):
        return entry[self.posDictonary["Sequence"]]

    def getDSSP(self,entry):
        return entry[self.posDictonary["DSSP"]]

    def getHead(self, entry):
        '''
        :param entry: entry with a header on position 1
        :return: the header of an entry
        '''
        return entry[self.posDictonary["HEAD"]]

    def getIDfromFile(self, entry):
        return entry[0]

    def getHeadfromFile(self, entry):
        return entry[1]

    def writeEntryToOutfile(self, entry, outfile):
        '''
        writes an entry to an outfile
        :param entry: an entry containig a single Sequence and DSSP
        :param outfile: file to write the entry to
        :return: none
        '''
        [outfile.write(data + '\t') for data in entry]
        outfile.write('\n')

    def isWanted(self, entry, validList):
        '''
        returns whether the entry with chain is present in the list of relavant sequences
        :param entry: an entry containig the ID and chain at the first two positions
        :return: boolean true or false
        TODO: stop hard coding 0 and 1 for chain and entry
        '''
        return entry[0] + entry[1] in validList

    def isValid(self, seq):
        """    Function to validate length and identity of a sequence. Combines isProtseq and correctSize.
        :param seq: Sequence used for validation
        :return: Returns True if both conditions are met.
        """
        if self.correctSize(seq) & self.isProtseq(seq):
            return True
        return False

    def isProtseq(self, string):
        """    Function to check for uppercases. AAs are always uppercases.
        :param string: String to be checked
        :return: True if the string is an AA-seq.
        """
        if string.isupper():
            return True
        return False

    def correctSize(self, seq):
        """     checks if an input string is greater than or equal to a given size.
        :param seq: Amino acid or DSSSP sequence
        :return: True if of correct size
        """
        param = 50
        if len(seq) >= param:
            return True
        return False

    def hasDSSP(self, inputlist):
        """    Checks whether Sequence is present in the string
        :param inputlist: List with strings
        :return: True if you have a DSSP
        """
        return np.any((['DSSP' in x for x in inputlist]))

    def changeLettersOfDSSp(self, entry):
        '''
        :param entry: entry DSSP
        :param changeWhat: list of characters
        :param changeTo: list of characters to change to
        return: none
        TODO: make it so that order of changes does not matter (now thought about input needed)
        '''

        for w, t in zip(self.changeWhat, self.changeTo):
            entry = entry.replace(w, t)
        return entry

# ----- testcode
if __name__ == '__main__':
    preproc = preProcessClass(outputFileName='simpleDSSPv2_all', debugFilename='debug', DatabaseName="data/PDBfinder2.txt",
                  changeLetters=True, wantedFile="data/wanted.txt")

    preproc.makeTextWithWantedInfo()