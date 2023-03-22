#!/usr/bin/env python
import sys
import os
import random
import string
import re

def main():
    #1 : .variant file prefix, 2: .perbase file prefix
    #$1.perbase is output from bam-readcount
    #they all should be sorted
    homedir=os.getenv("HOME")
    curdir=os.getcwd()
    pydir=os.path.dirname(os.path.realpath(__file__))
    
    VARFILE = open(sys.argv[1] + ".variant" , "r")
    RAWFILE = open(sys.argv[2] + ".perbase" , "r")
    line = 0

    posToVar   = dict()
    for var in VARFILE:
        varRow = var.rstrip().split("\t")
        if varRow[0] == "#CHROM":
            continue
        line = line + 1
        varChrom = varRow[0]
        varPosition = varRow[1]
        varRefBase = varRow[2].upper()
        varAltBase = varRow[3].upper()
        posToVar[varChrom + ":" + varPosition] = varRefBase + ":" + varAltBase

    for raw in RAWFILE:
        rawRow = raw.rstrip().split("\t")
        rawChrom = rawRow[0]
        rawPosition = rawRow[1]
        rawRefBase = rawRow[2].upper()
        rawTotalCount = rawRow[3]

        varRefBase = posToVar[rawChrom + ":" + rawPosition].split(":")[0];
        varAltBase = posToVar[rawChrom + ":" + rawPosition].split(":")[1];
        if(varRefBase != rawRefBase):
            raise("fatal error, position not matched")
        
        rawInfoCount = dict()
        for i in range(4,9): #-'base' is not included, column 5 to 10, =,A,C,G,T,N
            thisrawBase = rawRow[i]
            rawBaseInfo = thisrawBase.rstrip().split(":")
            rawInfoCount[rawBaseInfo[0]] = rawBaseInfo[1]
        refCount = rawInfoCount[varRefBase]
        def readRawInfoCount(x):
            return rawInfoCount[x]
        altCount = list(map(readRawInfoCount,  varAltBase.rstrip().split(",")))
        
        del rawInfoCount[varRefBase]
        for alt in varAltBase.rstrip().split(","):
            del rawInfoCount[alt]
        etcInfo = ""
        for base in list(rawInfoCount.keys()):
            etcInfo = etcInfo + base + ":" + rawInfoCount[base] + ","
        print("\t".join((rawChrom, rawPosition, rawRefBase, varAltBase,rawTotalCount, refCount, ",".join(map(str,altCount)), etcInfo)))
        

    VARFILE.close()
    RAWFILE.close()

if __name__ == "__main__":
        main()
