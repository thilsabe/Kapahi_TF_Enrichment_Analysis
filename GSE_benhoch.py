from GSEutils import *

import os
import sys
import math
import numpy
import time
import random
import re
import urllib
import bisect
import statsmodels.stats.multitest as multi
import RibotagExtract as RE

from operator import itemgetter, attrgetter, methodcaller

def print_and_log(msg, logf):
    print(msg)
    logf.write("%s\n" % msg)

def read_TF_binding_sites(TFChipSeqf, line, lineNum, nTF, TFon):

    scorePos = int(SETUP['SCOREPOS'])
    genePos = int(SETUP['GENEPOS'])
    TFPos = int(SETUP['TFPOS'])
    startPos = int(SETUP['STARTPOS'])
    endPos = int(SETUP['ENDPOS'])
    annotPos = int(SETUP['ANNOTPOS'])
    maxPos = max(scorePos, genePos, TFPos, startPos, endPos, annotPos)

    TFbindingLocs = {}
    sLine = line.strip().split('\t')
    if len(sLine) < maxPos:
       print_and_log("error in TF ChipSeq file format\n%s" % line, logf)
       sys.exit()
    TF = sLine[TFPos]
    TFon += 1
    prevTF = TF
    while line and prevTF == TF:
        lineNum += 1
        if TFon >= nTF:
            print("starting at "+str(TFon))
            gene = sLine[genePos].strip()
            try:
                score = float(sLine[scorePos])
            except:
                print("error in peak score")
                print(line, sLine[scorePos])
                sys.exit()
            start = sLine[startPos]
            end = sLine[endPos]
            location = sLine[annotPos]
            m = True
            if SETUP['GENOMIC_FEATURE'] == "PROMOTER":
                m = re.match('promoter\-TSS.*', location)
            if m:
                if gene.upper() not in TFbindingLocs:
                    TFbindingLocs[gene.upper()] = [[start, end, score]]
                else:
                    TFbindingLocs[gene.upper()].append([start, end, score])
        line = TFChipSeqf.readline()
        sLine = line.strip().split('\t')
        prevTF = TF
        while line and len(sLine) <= maxPos+1:
            print_and_log("error in TF ChipSeq file format in line %d\n'%s'" % (lineNum, line), logf)
            line = TFChipSeqf.readline()
            sLine = line.strip().split('\t')
           # prevTF = TF
        if line:
            TF = sLine[TFPos]

    #print "number of peaks for %s: " % TF, lineNum
    return prevTF, TFbindingLocs, line, lineNum, TFon

def max_score(geneTFpeaks):

    maxScore = geneTFpeaks[0][2]
    if len(geneTFpeaks) > 1:
        for peak in geneTFpeaks[1:]:
            if peak[2] > maxScore:
               maxScore = peak[2]
    return maxScore

'''
    get_pvalue arguments:
       the score calculated for the gene set in question
       the number of genes from the gene set that have peak values for this transcription factor
       Chip-seq peak values for a specific transcripthion factot
       the number of genes in this gene set
'''
def get_pvalue(score, num, TFpeaks, setSize, bckgGenes):

     def calc_pvalue_hist(rands, observation):

         hist, bin_edges = numpy.histogram(rands, bins=100)
         index = bisect.bisect(bin_edges, observation)
         #print "hist:", hist
         #print "bins:", bins
         #print "index, min-bins, max-bins:",index, min(bin_edges), max(bin_edges)
         tail = 0
         if index == 0:
            return 1.0
         for i in range(index-1, len(hist)):
            tail += hist[i]
         pvalue = pvalue/float(SETUP['PERMUTATIONS'])
         return pvalue

     def calc_pvalue(rands, observation):

         tail = sum(i >= observation for i in rands)
         if SETUP['EXECMODE'] == 'DEB':
            print("Number or simulations with score exceeding gene-set score: = %d" % tail)
         pvalue = tail/float(SETUP['PERMUTATIONS'])
         return pvalue

     #print "score, setSize, number of genes in the gene set that have peak values:", score, setSize, num
     if setSize == 0:
        return 0
     random.seed()
     norm_randScores = []
     randScores = []
     randNums = []
     # this is the bootstrapping loop
     for i in range(0, int(SETUP['PERMUTATIONS'])):
          random.seed(i+random.random())
          randScore = 0
          randNum = 0
          # from the background gene list select random set of genes of the size of the gene list in question
          randSample = random.sample(bckgGenes, setSize)
          for j in range(0, setSize):
              randGene = randSample[j]
              if randGene in TFpeaks:
                 randScore += max_score(TFpeaks[randGene])
                 randNum += 1
          randScores.append(randScore)    # accumulates for each iteration the score for this random gene set
          randNums.append(randNum)        # accumulates for each iteration the number of genes in this random gene set that have peaks
     if SETUP['EXECMODE'] == 'DEB':
         logf.write("random gene-sets scores:\n")
         for rs in randScores:
            logf.write("%.4f " % rs)
         logf.write("\n")
     score_pvalue = calc_pvalue(randScores, score)
     num_pvalue = calc_pvalue(randNums, num)
     return score_pvalue, num_pvalue

def get_set_score(geneSet, TFpeaks):

    numOfGenes = 0
    totalScore = 0
    geneList = []
    for gene in geneSet:
        if gene in TFpeaks:
            score = max_score(TFpeaks[gene])
            totalScore += abs(score)
            numOfGenes +=1
            geneList.append(RE.FBnumtoComName(gene))
            if SETUP['EXECMODE'] == 'DEB':
                logf.write("\t%s" % gene)
    if SETUP['EXECMODE'] == 'DEB':
       logf.write("\n")

    if numOfGenes == 0:
        norm_score = 0
    else:
        norm_score = totalScore/float(numOfGenes)

    return totalScore, numOfGenes, norm_score, geneList

def enrichment_analysis(TF, TFpeaks, allGenes, geneSet, bckgGenes, outf):

     geneSetSize = len(geneSet)
     score, numOfGenes, normScore, geneList = get_set_score(geneSet, TFpeaks)
     #geneList = [RE.FBnumtoComName(g.strip()) for g in geneList]
     pvalue_score, pvalue_num = get_pvalue(score, numOfGenes, TFpeaks, geneSetSize, bckgGenes)
     if SETUP['EXECMODE'] == 'DEB':
         print("%s: genes: %d, score: %f, normalized score: %f" % (TF, numOfGenes, score, normScore))
         print("%s: genes: %d, genes in backgroung: %d, pvalue(num): %f, score: %f, pvalue(score): %f\n" % (TF, len(TFpeaks), numOfGenes, pvalue_num, score, pvalue_score))
     outf.write("%s\t%d\t%d\t%f\t%f\t%f\t%s\n" % (TF, len(TFpeaks), numOfGenes, pvalue_num, score, pvalue_score, ", ".join(geneList)))
     outf.flush()

def write_output_header(outf, SETUP, origGeneSetSize, geneSetSize):

    outf.write("#Time:\t%s\n" % time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime()))   #localtime(time.time()))
    outf.write("#gene-set-file:\t%s\n" % SETUP['GENESET'])
    outf.write("#gene-set-name:\t%s\n" % SETUP['GENESETNAME'])
    outf.write("#original-gene-set-size\t%d\n" % origGeneSetSize)
    outf.write("#gene-set-size (after removing duplications and genes not in background list\t%d\n" % geneSetSize)
    outf.write("#Background-file:\t%s\n" % SETUP['BACKGROUND'])
    outf.write("#Permutations:\t%s\n" % SETUP['PERMUTATIONS'])
    outf.write("TF\tTF-genes\tTF-genes-in-gene-set\tpvalue-num\tscore\tpvalue-score\ttarget-genes-from-input\n")

def TFmulti(TFanalysisOutFile, fdr=0.05):

    TFanalysis_out=[]
    inputfile = str(TFanalysisOutFile)[(str(TFanalysisOutFile).find("'")+1):(str(TFanalysisOutFile).find(".out'")+4)]#str(TFanalysisOutFile).split(' ')[2][1:-2]
    with open(inputfile, 'r') as f:
        for line in f:
            TFanalysis_out.append(line.strip().split('\t'))
        TFanalysis=TFanalysis_out[8:len(TFanalysis_out)] #Removing TFanalysis-specific data and headers for pvalue correction below
        scorepvals = [float(TF[5]) for TF in TFanalysis]

        scorepvals_cor_bon = multi.multipletests(scorepvals, method = "bonferroni", alpha = fdr)
        scorepvals_cor_bh = multi.multipletests(scorepvals, method = "fdr_bh", alpha = fdr)
        TFanalysis_out[7].append('Benjamini-Hochberg Corrected Score Pvalue')
        TFanalysis_out[7].append('Benjamini-Hochberg Significant?')
        TFanalysis_out[7].append('Bonferroni Corrected Score Pvalue')
        TFanalysis_out[7].append('Bonferroni Significant?')
        for tf in range(0, len(TFanalysis)):
            TFanalysis_out[tf+8].append(scorepvals_cor_bh[1][tf])
            TFanalysis_out[tf+8].append(scorepvals_cor_bh[0][tf])
            TFanalysis_out[tf+8].append(scorepvals_cor_bon[1][tf])
            TFanalysis_out[tf+8].append(scorepvals_cor_bon[0][tf])

    with open(inputfile[:-4] + "_MultiCorrect.out", "w") as file:
        for nested_list in TFanalysis_out:
            for e in nested_list:
                 file.write(str(e) + '\t')
            file.write('\n')
    return TFanalysis_out

#To put into a single directory defined in genePipeline, I should be able to define the below code and pass the outputpath to it
def GSEanalysis(FDR, nTF, outputpath):
    geneSetf, TFChipSeqf, bckgf, outf, logf, SETUP = get_setup()

    bckgGenes, lineNum, duplicates = read_gene_set(bckgf, SETUP['BACKGROUND'], SETUP['BACKGROUNDNAME'], logf)
    origGeneSet, lineNum, duplicates = read_gene_set(geneSetf, SETUP['GENESET'], SETUP['GENESETNAME'], logf)
    geneSet = remove_genes_not_in_bckgrnd(origGeneSet, bckgGenes.keys(), logf)
    print_and_log("%d unique genes in gene-set also in background list" % (len(geneSet)), logf)
    if SETUP['EXECMODE'] == 'DEB':
        logf.write("Unique genes in gene-set and in background list: ")
        for gene in geneSet:
            logf.write(" %s")
        logf.write("\n")


    write_output_header(outf, SETUP, len(origGeneSet), len(geneSet))

    line = TFChipSeqf.readline() # header
    line = TFChipSeqf.readline()
    lineNum = 1
    TFon = 0
    while line:
         TF, peaks, line, lineNum, TFon = read_TF_binding_sites(TFChipSeqf, line, lineNum, nTF, TFon)
         print_and_log("%s:\t%d genes" % (TF, len(peaks.keys())), logf)
         if TFon >= nTF:
             enrichment_analysis(TF, peaks, bckgGenes, geneSet, bckgGenes.keys(), outf)

    TFanalysis_out = TFmulti(outf, FDR)
    return TFanalysis_out