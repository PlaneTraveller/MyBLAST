#!/usr/bin/env python3

#==============================================================================#
# This is an implementation of BLASTp algorithm                                #
#
# For testing purposes, the database is the proteome of whipworm (Trichuris
# trichiura) and the query is a protein from E. coli.
#==============================================================================#

#==============================================================================#
# Importing modules
from Bio import SeqIO
from Bio.Seq import Seq
import json
import blosum as bl
import math

subMat = bl.BLOSUM(62)
AA = ['C','S','T','A','G','P','D','E','Q','N','H','R','K','M','I','L','V','W','Y','F']

workDir = "./data/"

dbPath = workDir
dbFile = "whipworm" + "_db.json"

oriDBPath = workDir
oriDBFile = "whipworm" + ".faa"

qrPath = workDir
qrFile = "query.fasta"

resPath = workDir
resFile = "results" + ".csv"

k = 3
T = 13
E = 0.001

#==============================================================================#
# File Reading
hashedDB = {}

# Reading the json file from database
with open(dbPath+dbFile, "r") as f:
    # Final dictionary: {kmer: [{entryID: [positions]}
    hashedDB = json.load(f)




#==============================================================================#
# Functions
def makeSeeds(record):
    pieces = []
    seeds = {}
    querySeq = str(record.seq)
    # querySeq = str(record)

    # What I should really do is to ind. with protein kmers
    # Final dictionary: {kmer: [positions]}
    for i in range(len(querySeq)-k+1):
        pieces.append([str(querySeq[i:i+k]), i])


    for i in range(len(pieces)):
        seq = pieces[i][0]
        loc = pieces[i][1]
        seedList = findSub(seq, T)
        for seed in seedList:
            if seed in seeds:
                seeds[seed].append(loc)
            else:
                seeds[seed] = [loc]

    return seeds


def findSub(target, T):
    success = []

    for aa in AA:
        for bb in AA:
            for cc in AA:
                score = subMat[target[0] + aa]
                score += subMat[target[1] + bb]
                score += subMat[target[2] + cc]
                if score >= T:
                    success.append(aa + bb + cc)
    return success

# tmpseq = Seq("CSTCSTAGPDQNHRRKMLVWYF")
# print(makeSeeds(tmpseq))


def blast(record):
    print("Blasting...")
    # Reading the original database using SeqIO
    # with open(oriDBPath+oriDBFile, "r") as f:
        # oriDB = SeqIO.to_dict(f, "fasta")
    oriDB = SeqIO.index(oriDBPath+oriDBFile, "fasta")

    queryID = record.id
    querySeq = str(record.seq)
    seeds = makeSeeds(record)

    HSP = findHSP(seeds, querySeq, hashedDB, oriDB)
    return(HSP)


def findHSP(seeds, querySeq, db, oriDB):
    HSP = {}
    print("Finding HSPs...")
    for kmer in seeds:
    # Final dictionary: {kmer: [positions]}
        if kmer in db:
            # Final dictionary: {kmer: {entryID: [positions]}
            for entry in db[kmer]:
                for matchPos in db[kmer][entry]:
                    dbEntry = str(oriDB[entry].seq)
                    for seedPos in seeds[kmer]:
                        hspScore = elongation(dbEntry, querySeq, matchPos, seedPos, subMat)
                        if isSignificant(hspScore, len(dbEntry), len(querySeq)):
                            if entry not in HSP:
                                HSP[entry] = [matchPos, seedPos, hspScore]
                                print("HSP found!")
    return HSP


def isSignificant(hspScore, m, n):
    lam = 0.318
    k = 0.13
    # h = 0.4

    bitscore = (lam * hspScore - math.log(k)) / math.log(2)
    expect = m * n * 2**(-bitscore)
    # expect = k*m*n*math.exp(-lam*hspScore)

    # print(expect)
    if expect <= E:
        return True
    else:
        return False


def elongation(dbEntry, querySeq, matchPos, seedPos, subMat):
    score = 0
    maxscore = score
    # elongating to the right
    dbRange = len(dbEntry)
    qrRange = len(querySeq)

    for i in range(0, min(dbRange-matchPos, qrRange-seedPos)):
        a = i + matchPos
        b = i + seedPos
        sc = subMat[dbEntry[a]+querySeq[b]]
        # if sc > 0:
        if score+sc >= maxscore * 0.8:
            score += sc
            # tst.append(dbEntry[a])
        else:
            break
        if sc > 0: maxscore = score

    # tst.append("   ")
    for i in range(0, -min(matchPos, seedPos)-1, -1):
        a = i + matchPos
        b = i + seedPos
        sc = subMat[dbEntry[a]+querySeq[b]]
        # if score+sc >= sc * 0.8:
        if score+sc >= maxscore * 0.8:
        # if sc > 0:
            score += sc
            # tst.append(dbEntry[a])
        else:
            break
        if sc > 0: maxscore = score

    return score

# # testing the elongation function
# tst = []
# # print(subMat["AH"])
# print(elongation("CSTCSTAGPDQNHRRKMLVWYF",
#                "ABCSTCSTAGPDQNHRRKMLVWYF", 3, 5, subMat))
# print(tst)

#==============================================================================#
# Main

def main():
    # Reading the query from file
    with open(qrPath+qrFile, "r") as f:
        print("Reading query...")
        qr = SeqIO.parse(f, "fasta")
        HSP = []

        for record in qr:
            print(record.id)
            HSP.append(blast(record))

        print(HSP)

if __name__ == "__main__":
    main()


#==============================================================================#
# Sample results
# [{'CDW60865.1': [95, 121, 956.0], 'CDW57800.1': [493, 22, 57.0], 'CDW60670.1': [24, 55, 54.0]}, {'CDW60470.1': [115, 157, 648.0]}, {'CDW58513.1': [0, 0, 1992.0], 'CDW56129.1': [14, 218, 57.0], 'CDW60852.1': [39, 218, 56.0], 'CDW60899.1': [47, 218, 54.0], 'CDW57425.1': [18, 221, 54.0]}]
