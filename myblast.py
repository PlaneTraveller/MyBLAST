#!/usr/bin/env python3

#==============================================================================#
# This is My implementation of BLASTp algorithm using BLOSUM62 matrix.
#
# For testing purposes, the database is the proteome of whipworm (Trichuris
# trichiura) and the query is a protein from E. coli.
#
# TODO: Parallelize the code
# TODO: Validity check
# TODO: Argparse
# TODO: Local alignment algorithm
# TODO: Result outputting
#==============================================================================#

#==============================================================================#
# Importing modules
from Bio import SeqIO
from Bio.Seq import Seq
import json
import blosum as bl
import math
from multiprocessing import Pool
from multiprocessing import Process
import argparse
import csv


AA = ['C','S','T','A','G','P','D','E','Q','N','H','R','K','M','I','L','V','W','Y','F']

#==============================================================================#
# Argument parsing
def parseArgs():
    parser = argparse.ArgumentParser(description="My implementation of BLASTp")
    parser.add_argument(
        "-k",
        "--kmer",
        type=int,
        default=3,
        help="Kmer length",
    )
    parser.add_argument(
        "-T",
        "--threshold",
        type=int,
        default=13,
        help="Threshold for kmer",
    )
    parser.add_argument(
        "-E",
        "--expect",
        type=float,
        default=0.001,
        help="Expectation value",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of threads",
    )
    parser.add_argument(
        "-w",
        "--workdir",
        type=str,
        default="./data/",
        help="Working directory",
    )
    parser.add_argument(
        "-d",
        "--database",
        type=str,
        # default="whipworm",
        help="Database file name",
    )
    parser.add_argument(
        "-q",
        "--query",
        type=str,
        # default="query.fasta",
        help="Query file name",
    )
    parser.add_argument(
        "-r",
        "--result",
        type=str,
        default="results.csv",
        help="Result file name",
    )
    parser.add_argument(
        "-s",
        "--submat",
        type=int,
        default=62,
        help="Substitution matrix",
    )


    return parser.parse_args()


#==============================================================================#
# Parameters
options = parseArgs()

k = options.kmer
T = options.threshold
E = options.expect
threads = options.threads

workDir = options.workdir

dbPath = workDir
dbFile = options.database.split(".")[0] + "_db.json"
oriDBPath = workDir
oriDBFile = options.database
qrPath = workDir
qrFile = options.query
resPath = workDir
resFile = options.result


subMat = bl.BLOSUM(options.submat)


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
            # Parellelization Handel

            # procs = []
            # for entry in db[kmer]:
            #     proc = Process(target=mpWorker, args=(entry, HSP, kmer, oriDB, db, seeds, querySeq))
            #     procs.append(proc)
            #     proc.start()

            # for proc in procs:
            #     proc.join()

            # with Pool(threads) as p:
            #     p.map(mpWorker, db[kmer])

            for entry in db[kmer]:
                mpWorker(entry, HSP, kmer, oriDB, db, seeds, querySeq)
    return HSP

def mpWorker(entry, HSP, kmer, oriDB, db, seeds, querySeq):
    if entry in HSP: return
    for matchPos in db[kmer][entry]:
        dbEntry = str(oriDB[entry].seq)
        for seedPos in seeds[kmer]:
            hspScore = elongation(dbEntry, querySeq, matchPos, seedPos, subMat)
            if isSignificant(hspScore, len(dbEntry), len(querySeq)):
                # if entry not in HSP:
                HSP[entry] = [matchPos, seedPos, hspScore]
                print("HSP found!")

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
        HSP = {}

        for record in qr:
            print(record.id)
            HSP[record.id] = blast(record)

        with open(resPath+resFile, "w") as f:
            f.write("Query,Hit,MatchPos,SeedPos,Score\n")
            for query in HSP:
                for hit in HSP[query]:
                    f.write(query + "," + hit + "," + ",".join(str(e) for e in HSP[query][hit])+ "\n")

        print(HSP)

if __name__ == "__main__":
    main()


#==============================================================================#
# Sample results
# HSP = [{'CDW60865.1': [95, 121, 956.0], 'CDW57800.1': [493, 22, 57.0], 'CDW60670.1': [24, 55, 54.0]}, {'CDW60470.1': [115, 157, 648.0]}, {'CDW58513.1': [0, 0, 1992.0], 'CDW56129.1': [14, 218, 57.0], 'CDW60852.1': [39, 218, 56.0], 'CDW60899.1': [47, 218, 54.0], 'CDW57425.1': [18, 221, 54.0]}]
