from __future__ import print_function
import sys
if sys.version_info.major < 3:
    input = raw_input

import sys
import time
import os
import RibotagExtract as RE

setupFileName = 'GSE.setup'

def open_file(fileName, mode, msg):

    fp = ''
    flag = False
    if mode == 'w':
       if os.path.exists(fileName):
          print("%s file %s already exist" % (msg, fileName))
          overwrite = input("do you want to overwrite it's content? [y|n]")
          if not (overwrite == 'y'):
              sys.exit()
          try:
              fp = open(fileName, "w")
          except:
              print("%s file %s cannot be open" % (msg, fileName))
              flag = True
       else:
          fp = open(fileName, "w")
       return flag, fp

    elif mode == 'r':
       try:
          fp = open(fileName, 'r')
       except:
          print("%s file %s does not exist" % (msg, fileName))
          flag = True
       return flag, fp
    else:
       print("Script error: mode '%s' is not valid" % mode)
       sys.exit()


def help_message():

        print("setup file %s must exist in the local directory" % setupFileName)
        print("The format of each line in the setup file must be: keyword-TAB-value")
        print("Required keywords are: EXECMODE, PERMUTATIONS, GENESET, GENESETNAME, TFCHIPSEQ, BACKGROUND, BACKGROUNDNAME, GENOMIC_FEATURE, TFPOSE, GENEPOS, ANNOTPOS, STARTPOS, ENDPOS, SCOREPOS, OUTFILE, LOGFILE")
        print("Please correct and run again")
        sys.exit()

SETUP = {}
INTKEYWORDS = {'PERMUTATIONS', 'TFPOS', 'GENEPOS', 'STARTPOS', 'ANNOTPOS', 'ENDPOS','SCOREPOS'}
STRKEYWORDS = {'EXECMODE', 'GENESET', 'GENESETNAME', 'TFCHIPSEQ', 'BACKGROUND', 'BACKGROUNDNAME','GENOMIC_FEATURE', 'OUTFILE', 'LOGFILE'}

def get_setup():

    try:
        setUpf = open(setupFileName, "r")
    except:
        help_message()
    lines = setUpf.readlines()
    if len(lines) < 7:
        help_message()
    error = False
    lineNum = 0

    for line in lines:
        lineNum += 1
        sLine = line.strip().split('\t')
        if len(sLine) < 2:
             print("Line number %d in %s is not in the correct format. %s" % (lineNum, setupFileName, line))
             error = True
        keyword = sLine[0].strip()
        if keyword not in INTKEYWORDS.union(STRKEYWORDS):
            print("keyword '%s' in %s is not a valid keyword" % (keyword, setupFileName))
            sys.exit()
        value = sLine[1].strip()
        SETUP[keyword] = value
        print(keyword, value)

        if keyword == 'EXEC_MODE':
            if (value != 'DEB') and (value != 'NODEB'):
                print("Unrecognized value for 'EXEC_MODE' keyword in setup file %s. Value can be 'DEB' or 'NODEB' but it is %s" % (setupFileName, value))
                error = True

        elif keyword == 'GENESET':
            flag, geneSetf = open_file(value, 'r', 'Gene-set')
            if flag: error = True

        elif keyword == 'TFCHIPSEQ':
            flag, TFChipSeqf = open_file(value, 'r', 'Chipseq')
            if flag: error = True

        elif keyword == 'BACKGROUND':
            flag, bckgf = open_file(value, 'r', 'Background')
            if flag: error = True

        elif keyword in INTKEYWORDS:
            try:
                value = int(value)
            except:
                print("%s keyword must follow by an integer" % keyword)
                error = True

        elif keyword == 'GENOMIC_FEATURE':
            if value != 'ALL' and value != 'PROMOTER':
                print("Unrecognized value for 'GENOMIC_FEATURE' keyword in setup file %s. Value can be 'ALL' or 'PROMOTER' but it is" % (setupFileName, value))
                error = True

        elif keyword == 'OUTFILE':
            flag, outf = open_file(value, 'w', 'Output')
            if flag: error = True

        elif keyword == 'LOGFILE':
            flag, logf = open_file(value, 'w', 'Log')
            if flag: error = True

    if error:
        print("There are errors in the setup file.")
        print("Please correct and run again.")
        sys.exit()

    for keyword in INTKEYWORDS.union(STRKEYWORDS):
        if not keyword in SETUP:
           print("missing keyword '%s' in setup file 'geneSetEnrichment.setup'" % keyword)
           sys.exit()
    return geneSetf, TFChipSeqf, bckgf, outf, logf, SETUP



def print_and_log(msg, logf):
    print(msg)
    logf.write("%s\n" % msg)


def read_gene_set(genef, fileName, geneListName, logf):

     print_and_log("... Reading %s genes from file %s" % (geneListName, fileName), logf)
     duplicates = 0
     lineNum = 0
     genes = {}
     lines = genef.readlines()
     for line in lines:
         if line == '':
             continue
         lineNum += 1
         sLine = line.split()
         try:
             gene = sLine[0].strip()
         except:
             print("error reding gene name from file")
             print("line %d: '%s'" % (lineNum, line))
             sys.exit()
         uGene = gene.upper()
         if uGene in genes:
            print_and_log("duplicated gene %s in list %s" % (gene, geneListName), logf)
            duplicates += 1
         else:
            genes[uGene] = gene

     print_and_log("%d rows\n%d duplictes\n%d unique genes\n" % (lineNum, duplicates, len(genes)), logf)

     return genes, lineNum, duplicates

def remove_genes_not_in_bckgrnd(geneSet_all, bckgGenes, logf):

    geneSet = []
    for gene in geneSet_all:
        if gene in bckgGenes:
            geneSet.append(gene)
        else:
            print_and_log("gene '%s' in gene set is not in the background set, therefor it is removed from the set" % geneSet_all[gene], logf)
    return geneSet


if __name__ == "__main__":

    print(".....")
    print("This is not a stand-alone program.")
    print("It is intended to be imported by a main script, like 'GSE.py'.")
    print(".....")
