from countClass import countClass
from scoreClass import scoreClass
from mkname import createPdbs
from changeIDs import changeID
import pickle
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

import json

def readAllIDs(jsonFile, fastaFile):
    with open(jsonFile) as json_file:
        data = json.load(json_file)

    with open(fastaFile) as fastaFile:
        fasta_dict = {}

        for line in fastaFile:
            if line[0] == ">":
                fasta = line[1:].strip()
            else:
                indx = data[fasta][0]
                fasta_dict[fasta] = [line.strip(), line[indx[0]-1:indx[1]], indx]
    return fasta_dict

def createDSSP(seq, fillChar = 'H'):
    return ''.join([fillChar for _ in seq])

if __name__ == '__main__':
    fasta_dict = readAllIDs("/mnt/D_disk/My_data/Jarek/Documents/Studie/2019-2020/CF-check-master/data/disprot_regions.json", "/mnt/D_disk/My_data/Jarek/Documents/Studie/2019-2020/CF-check-master/data/seqs_disorder_content_all_above_0.6.fasta")

    N=2
    c = countClass(N,"./preprocess_results/simpleDSSPv2_all.txt",
        "countresults/simpleDSSPv2.txt",
        skipSymbols = True)
    # expectedVals, observedVals, scoreMatrix, scoreDict = c.createPrefParams()
    # with open('scoreDict.p', 'wb') as fp:
    #     pickle.dump(scoreDict, fp)

    with open('scoreDict.p', 'rb') as fp:
        scoreDict = pickle.load(fp)


    s = scoreClass(N, scoreDict, "./preprocess_results/simpleDSSPv2_all.txt",
                   "./analysis_results/joanna.txt",
                   removeDAA=True) #NOTE, WE ONLY USE IT TO LOAD DATA
    _, entryScores, _, _ = s.calcScores()
    allprot = s.allProteins
    allprot = [p[3] for p in allprot ]
    randomProts = allprot

    median = np.median(entryScores)
    entryScores -= median

    unfoldedScores_H = []
    unfoldedScores_E = []
    unfoldedScores_T = []
    for ID in fasta_dict.keys():
        chain = ""
        HEAD = "JoannaProt"
        fullSeq = fasta_dict[ID]
        unfoldSeq = fasta_dict[ID][1]

        for scoreArray, dssp in zip([unfoldedScores_H, unfoldedScores_E, unfoldedScores_T], ["H", "E", "T"]):
            DSSP = createDSSP(unfoldSeq, fillChar = dssp)
            entry = [ID, chain, HEAD, unfoldSeq, DSSP]
            _, score, _, _ = s.calcScorePerID(entry, scoreDict)
            scoreArray += [score - median]

    randoms_H = []
    randoms_E = []
    randoms_T = []
    for seq in randomProts:
        chain = ""
        ID = "NUP"
        HEAD = "NOPPES"

        for scoreArray, dssp in zip([randoms_H, randoms_E, randoms_T], ["H", "E", "T"]):
            DSSP = createDSSP(seq, fillChar = dssp)
            entry = [ID, chain, HEAD, seq, DSSP]
            _, score, _, _ = s.calcScorePerID(entry, scoreDict)
            scoreArray += [score - median]

    sns.distplot(np.array(entryScores), hist=False, kde=True,
                 kde_kws={'shade': True, 'linewidth': 2,  "color": 'gray'},
                 hist_kws={"alpha": 0.3, "color": 'gray'},
                 label="allProt")
    sns.distplot(np.array(unfoldedScores_H ), hist=False, kde=True,
                 kde_kws={'shade': True, 'linewidth': 2, "color": 'b'},
                 hist_kws={"alpha": 0.3},
                 label="H_given")
    sns.distplot(np.array(unfoldedScores_E ), hist=False, kde=True,
                 kde_kws={'shade': True, 'linewidth': 2, "color": 'g'},
                 hist_kws={"alpha": 0.3},
                 label="E_given")
    sns.distplot(np.array(unfoldedScores_T ), hist=False, kde=True,
                 kde_kws={'shade': True, 'linewidth': 2, "color": 'm'},
                 hist_kws={"alpha": 0.3},
                 label="T_given")
    sns.distplot(np.array(randoms_H), hist=False, kde=True,
                 kde_kws={'shade': True, 'linewidth': 2, 'linestyle':'--', "color": 'b'},
                 hist_kws={"alpha": 0.3},
                 label="H_allProt")
    sns.distplot(np.array(randoms_E), hist=False, kde=True,
                 kde_kws={'shade': True, 'linewidth': 2, 'linestyle':'--', "color": 'g'},
                 hist_kws={"alpha": 0.3},
                 label="E_allProt")
    sns.distplot(np.array(randoms_T), hist=False, kde=True,
                 kde_kws={'shade': True, 'linewidth': 2, 'linestyle':'--', "color": 'm'},
                 hist_kws={"alpha": 0.3},
                 label="T_allProt")
    plt.savefig("disorder_basic.pdf", dpi=500)
    plt.xlabel("preference score")
    plt.ylabel("normalized frequency")
    plt.title("Disorder content scores")
    plt.savefig("disorder_labels.pdf", dpi=500)
    plt.legend()
    plt.savefig("disorder_all.png", dpi=500)
    plt.show()