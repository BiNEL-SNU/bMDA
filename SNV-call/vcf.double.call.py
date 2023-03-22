#!/usr/bin/env python
import sys
import os
import random
import string
import re
import numpy as np

homedir=os.getenv("HOME")
curdir=os.getcwd()
pydir=os.path.dirname(os.path.realpath(__file__))


def main():
    ## this script performs inter sample double calling from vcf files obtained by GATK, varscan, MuTect tools (Lab note KJH35). note that this script applies only to SNP
    #1 : bed files of triple called sites, 2: list of target files 3: output file name 4: threshold
    #bed file must be sorted in chromosome coordinates

    ##initialize parameters
    nTH = int(sys.argv[4])
    depTH = 5
    listCaller = ["GATK", "varscan", "MuTect"]
    listGenotype = ["gt", "ref", "alt", "refCount", "altCount"]
    numGT = len(listGenotype)
    dictFileSuffix = dict()
    dictFileSuffix["GATK"] = "GATK.vcf"
    dictFileSuffix["varscan"] = "varscan.snp.Somatic.hc.filtered"
    dictFileSuffix["MuTect"] = "mutect1.final.vcf"

    ##setting extra parameters

    FILEtcsites = open(curdir+"/"+sys.argv[1] , "r")
    FILElist = open(curdir+"/"+sys.argv[2] , "r")
    FILEOUT = open(curdir + "/" + sys.argv[3] + ".variant", "w")

    listTarget = FILElist.read().strip().split("\n")
    numTarget = len(listTarget)

    posToIndex = dict()
    idxToPos = []
    line = 0
    for sites in FILEtcsites:
        thisSite = sites.strip().split("\t")
        posToIndex[thisSite[0] + ":" + thisSite[1] + ":" + thisSite[2]] = line
        idxToPos.append(thisSite[0] + ":" + thisSite[1] + ":" + thisSite[2])
        line = line + 1
    numSites = line
    print(f"number of double called sites: {numSites}")

    #4-dimensional gtArray, dim1 : chr position, dim2: sample, dim3: variant caller (GATK, varscan, MuTect), dim4: genotype info (gt, ref, alt, ref countt, alt count)
    gtArray = np.chararray((numSites, numTarget, len(listCaller) + 1, numGT), unicode=True, itemsize=100)
    for i in range(numSites):
        for j in range(numTarget):
            for k in range(len(listCaller) + 1):
                for l in range(numGT):
                    gtArray[i,j,k,l] = ""


    for target in listTarget:
        targetIndex = listTarget.index(target)
        for caller in range(len(listCaller)):
            nameCaller = listCaller[caller]
            result = []
            if nameCaller == "GATK":
                result = vcfParser(target + "." + dictFileSuffix[nameCaller], target)
            elif nameCaller == "varscan":
                result = varscanParser(target + "." + dictFileSuffix[nameCaller])
            elif nameCaller == "MuTect":
                result = vcfParser(target + "." + dictFileSuffix[nameCaller], target)
            else:
                raise Exception("Wrong variant caller name")
            # print(nameCaller, result)
            for line in result:
                pos = line[0] + ":" + str(int(line[1]) - 1) + ":" + line[1]
                if pos in posToIndex:
                    gtArray[posToIndex[pos], targetIndex, listCaller.index(nameCaller), 0] = line[4]
                    gtArray[posToIndex[pos], targetIndex, listCaller.index(nameCaller), 1] = line[2]
                    gtArray[posToIndex[pos], targetIndex, listCaller.index(nameCaller), 2] = line[3]
                    gtArray[posToIndex[pos], targetIndex, listCaller.index(nameCaller), 3]= line[5][0]
                    gtArray[posToIndex[pos], targetIndex, listCaller.index(nameCaller), 4] = line[5][1]

    countGATK = 0;
    countVarscan = 0;
    countMuTect = 0;
    countCS = 0
    listCS = []
    for i in range(numSites):  ##compute confident sites
        # print(list(gtArray[i, :, 0, 0]).count("0/1"))
        if ((list(gtArray[i, :, 0, 0]).count("0/1") + list(gtArray[i, :, 0, 0]).count("1/1")) >= nTH):
            countGATK += 1
        if ((list(gtArray[i, :, 1, 0]).count("0/1") + list(gtArray[i, :, 1, 0]).count("1/1")) >= nTH):
            countVarscan += 1
        if ((list(gtArray[i, :, 2, 0]).count("0/1") + list(gtArray[i, :, 2, 0]).count("1/1")) >= nTH):
            countMuTect += 1
        if ((list(gtArray[i, :, 0, 0]).count("0/1") + list(gtArray[i, :, 0, 0]).count("1/1")) >= nTH) | ((list(gtArray[i, :, 1, 0]).count("0/1") + list(gtArray[i, :, 1, 0]).count("1/1")) >= nTH) | ((list(gtArray[i, :, 2, 0]).count("0/1") + list(gtArray[i, :, 2, 0]).count("1/1")) >= nTH):
            countCS += 1;
            listCS.append(i)

    print(countGATK, countVarscan, countMuTect, countCS)

    ## among confident sites, samples called by at least one caller is considered called (sites not in confident sites is out of interest)
    for i in range(numSites):
        for j in range(numTarget):
            this = list(gtArray[i,j, list(range(len(listCaller))), 0])
            if (this.count("0/1") + this.count("1/1") > 0):
                gtArray[i,j,len(listCaller), 0] = "1"
            else:
                gtArray[i,j,len(listCaller), 0] = "0"
            if gtArray[i,j, 0, 1] != '':
                for l in range(numGT)[1:]:
                    gtArray[i,j,len(listCaller), l] = gtArray[i,j, 0, l]  ##by default use info of GATK (if the information exist)
            elif gtArray[i,j, 1, 1] != '':
                for l in range(numGT)[1:]:
                    gtArray[i,j,len(listCaller), l] = gtArray[i,j, 1, l]  ##use second caller if not exist
            elif gtArray[i,j, 2, 1] != '':
                for l in range(numGT)[1:]:
                    gtArray[i,j,len(listCaller), l] = gtArray[i,j, 2, l]  ##use third caller if not exist

    gtArrayCS = gtArray[listCS, :, len(listCaller), :]  ##3 dimension. dim1 : chr position, dim2: samples, dim3: genotype
    idxToPosCS = list(np.array(idxToPos)[listCS])
    posToIndexCS = dict()
    for n in range(len(idxToPosCS)):
        posToIndexCS[idxToPosCS[n]] = n;
    FILEOUT.write("CHROM\tPOS\tREF\tALT\t" + "\t".join(listTarget)+ "\n")
    for i in range(countCS):
        thisPos = idxToPosCS[i].split(":")
        thisREF = sorted(list(gtArrayCS[i,:, 1]))[numTarget - 1]
        if sum(1 for i in list(set(list(gtArrayCS[i,:, 2]))) if i in ['A','G','C','T']) != 1: ##this filter removes ambiguous double alt site (different sample A, T + inside sample A,T)
            listCS.remove(posToIndex[idxToPosCS[i]])
            continue
        thisALT = sorted(list(gtArrayCS[i,:, 2]))[numTarget - 1]
        FILEOUT.write("\t".join([thisPos[0], thisPos[2], thisREF, thisALT]))
        for j in range(numTarget):
            FILEOUT.write("\t" + gtArrayCS[i,j, 0])
        FILEOUT.write("\n")
    gtArrayCS = gtArray[listCS, :, len(listCaller), :]  ##3 dimension. dim1 : chr position, dim2: samples, dim3: genotype
    idxToPosCS = list(np.array(idxToPos)[listCS])
    posToIndexCS = dict()
    for n in range(len(idxToPosCS)):
        posToIndexCS[idxToPosCS[n]] = n;

#    print idxToPos[0], listTarget[0], gtArrayCS[9, :, :]




    FILEtcsites.close()
    FILElist.close()
    FILEOUT.close()

def vcfParser(nameFile, sample):
    FILEvcf = open(curdir + "/" + nameFile)

    result = []
    line = 0
    thisGT = -1
    vcfheader = re.compile(r'^#')
    for vcf in FILEvcf:
        vcfRow = vcf.rstrip().split("\t")
        if vcfheader.match(vcf[0]):
            if vcfRow[0] == "#CHROM":
                thisGT = vcfRow.index(sample)
            continue
        line = line + 1
        vcfChrom = vcfRow[0]
        vcfPosition = vcfRow[1]
        vcfRefBase = vcfRow[3].upper()
        vcfAltBase = vcfRow[4].upper()
        vcfVarQ = vcfRow[5]
        vcfGTinfo = vcfRow[8].split(":")
        vcfGenotype = vcfRow[thisGT].split(":")
        vcfGT = vcfGenotype[vcfGTinfo.index("GT")]
        vcfAD = vcfGenotype[vcfGTinfo.index("AD")].split(",")
        vcfDP = vcfGenotype[vcfGTinfo.index("DP")]
        result.append([vcfChrom, vcfPosition, vcfRefBase, vcfAltBase, vcfGT, vcfAD])
    numLine = line
    FILEvcf.close()
    return result

def varscanParser(nameFile):
    FILEvs = open(curdir + "/" + nameFile)
    result = []
    line = 0
    for vs in FILEvs:
        vsRow = vs.rstrip().split("\t")
        if vsRow[0] == "chrom":
            continue
        line = line + 1
        Chrom = vsRow[0]
        Position = vsRow[1]
        RefBase = vsRow[2].upper()
        AltBase = vsRow[3].upper()
        AD = [vsRow[8], vsRow[9]]
        if vsRow[10] == "100%":  ##this needs minor revision
            GT = "1/1"
        else:
            GT = "0/1"
        result.append([Chrom, Position, RefBase, AltBase, GT, AD])
    numLine = line
    FILEvs.close()
    return result

if __name__ == "__main__":
        main()
