#!/usr/bin/env python3

#==============================================================================#
# This is an implementation of protein database processing for BLAST algorithm #
#==============================================================================#

#==============================================================================#
# Importing modules
from Bio import SeqIO
from Bio.Seq import Seq
import json

dataDir = "./data/"
dbFile = "whipworm.faa"

writeDir = "./data/"
writeFile = dbFile.split(".")[0] + "_db.json"

k = 3

# Hashing the database into a dictionary
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





    
























#==============================================================================#
# Importing modules
#==============================================================================#
# Importing modules
#==============================================================================#
# Importing modules
#==============================================================================#
# Importing modules
#==============================================================================#
# Importing modules
#==============================================================================#
# Importing modules
