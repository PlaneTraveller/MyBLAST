#!/usr/bin/env python3

#==============================================================================#
# This is an implementation of protein database processing for BLAST algorithm #
#==============================================================================#

#==============================================================================#
# Importing modules
from Bio import SeqIO
from Bio.Seq import Seq
import json
import argparse


#==============================================================================#
# Argument parser
def parseArgs():
    parser = argparse.ArgumentParser(description="My implementation of BLASTp database constructor")
    parser.add_argument(
        "-k",
        "--kmer",
        type=int,
        default=3,
        help="Kmer length",
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


    return parser.parse_args()


#==============================================================================#
# Parameters
options = parseArgs()

k = options.kmer
dataDir = options.workdir
dbFile = options.database

writeDir = dataDir
writeFile = dbFile.split(".")[0] + "_db.json"



#==============================================================================#
# Hashing
# Final dictionary: {kmer: [{entryID: [positions]}
with open(dataDir+dbFile, "r") as f:
    parsed_db = SeqIO.parse(f, "fasta")
    db = {}
    for record in parsed_db:
        entryID = record.id
        entrySeq = record.seq

        for i in range(len(entrySeq)-k+1):
            kmer = entrySeq[i:i+k]
            kmerStr = str(kmer)
            if kmer in db:
                if entryID in db[kmer]:
                    db[kmerStr][entryID].append(i)
                else:
                    db[kmerStr][entryID] = [i]
            else:
                db[kmerStr] = {entryID: [i]}
        # print(db)

    # Writing the dictionary to a file
    with open(writeDir+writeFile, "w") as f:
        f.write(json.dumps(db, indent = 4))



