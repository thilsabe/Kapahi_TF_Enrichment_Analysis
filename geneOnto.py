# -*- coding: utf-8 -*-
"""
Created on Sun May  3 10:38:43 2020

@author: tyler
"""


from goatools.go_enrichment import GOEnrichmentStudy
from goatools import obo_parser
import Bio.UniProt.GOA as GOA
from IPython.display import Image
import gzip
import RibotagExtract as RE
import os

def geneOnto(inputGeneList, genesetname, save=0, methods = ["bonferroni", "sidak", "holm", "fdr"], fdr_thresh = 0.01):
    #inputGeneList can be either a list of gene names or the name of a file containing gene names
    #save = 0 (default), 1 to write output GO to current directory
    #methods = ['fdr'] (default), can be "bonferroni", "sidak", "holm", "fdr"
    #fdr_thresh = 0.01 (default), used to set threshold for inclusion in return variable
    #Converting input genes to FB numbers
    if os.path.isfile(str(inputGeneList)):
        filestem = inputGeneList[:-4]
        inputGeneList_FB = []
        with open(inputGeneList, 'r') as f:
            for line in f:
                if line.startswith('Gene'):
                    continue
                else:
                    inputGeneList_FB.append(RE.geneFBnum(line.strip()).lower())
    else:
        filestem = genesetname
        inputGeneList_FB = [RE.geneFBnum(g) for g in inputGeneList]
    
    #Setting up GO structure
    GOstructure = obo_parser.GODag('go-basic.obo')
    
    #Now setting up fly annotation file
    # filename = <LOCATION OF GAF FILE>
    
    geneAssocFile = 'gene_association.fb.gz' #GAF file with fly annotations and FB numbers from flybase.org
    with gzip.open(geneAssocFile, 'rt') as fly_FB_fp:
        geneAssociations = {}  # Initialise the dictionary of functions
        for entry in fly_FB_fp:
            if entry.startswith('!'):
                continue
            else:
                temp = entry.strip().split('\t')
                if temp[1].lower() not in geneAssociations:
                    geneAssociations[temp[1].lower()] = set([temp[4]])
                else:
                    geneAssociations[temp[1].lower()].add(temp[4])

    GOgene_dict = {x: geneAssociations[x]
                   for x in geneAssociations 
                   if x in inputGeneList_FB}
    
    popFB = geneAssociations.keys()
    assocFB = {}
    for x in geneAssociations:
        assocFB[x] = geneAssociations[x]
    
    studyFB = GOgene_dict.keys()
    #methods = ["bonferroni", "sidak", "holm", "fdr"]
    gFB = GOEnrichmentStudy(popFB, assocFB, GOstructure,
                             propagate_counts=True,
                             alpha=0.05,
                             methods=methods)
    gFB_res = gFB.run_study(studyFB)
    sFB_fdr = []
    for x in gFB_res:
        if x.p_fdr <= fdr_thresh:
            sFB_fdr.append([x.goterm.id, x.name, x.depth, x.p_fdr, x.study_items])
    if save == 1:
        gFB.wr_xlsx(filestem+"flyGOanalysis.xlsx", gFB_res)
    
    go_id2 = 'GO:1903818'
    rec = GOstructure[go_id2]    
        
    return sFB_fdr