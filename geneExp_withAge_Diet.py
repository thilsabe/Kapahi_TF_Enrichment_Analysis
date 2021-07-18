# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 2019
Code adapted from Project Code
@author: Tyler
General Code Steps:
Step 1: Import the file containing normalized expression values for genes on two diets over 4 timepoints.
Step 2: Convert list of input genes to their flybase codes by comparing against a dictionary of common names downloaded/identified on flybase.org (download if not possible to query directly from the site).
Step 3: Correlate values of input list genes, returning the correlation coefficients.
Step 4: If there is a strong correlation between input genes, will run a linear regression and return the slope.
Step 5: Cluster analysis returning a list of genes that group with the input list. Not sure which analysis, yet.
"""
#from sklearn import linear_model
import io
import sys
import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import RibotagExtract as RE
import csv
import scipy as sp

def geneExp_Age(targetGenes, geneListname, save):
    #Diverting print outputs
    text_trap = io.StringIO()
    sys.stdout = text_trap
    #To restore output: sys.stdout = sys.__stdout__
    geneExpData = 'OXR1 RNAi RNA-seq data_FB_genePipeline.csv' #geneExpData should be a matrix with gene FB#s listed with expression data at different timepoints
    #targetGenes =  #targetGenes should be a list of genes, FB#s or gene symbols as listed on flybase.org

    #Converting targetGene list into FB#s, if given as gene symbols
    targetGenes_FB = {}
    for gene in targetGenes:
        if gene.lower().startswith('fbgn'):
            targetGenes_FB[RE.FBnumtoComName(str(gene))] = gene.lower()
            if (targetGenes_FB[RE.FBnumtoComName(str(gene))] == 'NotFound'):
                targetGenes_FB.pop(gene, None)
        else:
            targetGenes_FB[str(gene)] = RE.geneFBnum(gene)
            if (targetGenes_FB[str(gene)] == 'NotFound'):
                targetGenes_FB.pop(gene, None)

    #Opening gene expression data
    expvals = open(geneExpData, 'r')
    expvals_clean = expvals.read().splitlines()
    expvals.close()
    expvals_clean = [gene.split(',') for gene in expvals_clean]
    for i in range(1, len(expvals_clean)):
        expvals_clean[i][0] = RE.geneFBnum(expvals_clean[i][0])
    expvals = {}
    expvals_targets = {}
    for gene in expvals_clean:#[1:len(expvals_clean)]: #Getting expression values and averaging them for each time point
        if gene[0] == 'Gene':
            continue
        if gene[0] == 'NotFound':
            continue
        for timepoint in range(1,len(gene)):
#            if gene[timepoint] == '':
#                gene[timepoint] = None
            try:
                gene[timepoint] = float(gene[timepoint])
            except:
                continue
        if gene[0] == 'FB#':
            expvals[gene[0]] = gene[1:len(gene)]
            #expvals_avg[gene[0]] = ['10', '20', '25', '30', '40']
        else:
            expvals[gene[0]] = [float(value) for value in gene[1:len(gene)]]
            #expvals_avg[gene[0]] = [(sum(gene[1:3])/float(len(gene[1:3]))), (sum(gene[4:6])/float(len(gene[4:6]))), (sum(gene[7:9])/float(len(gene[7:9]))), (sum(gene[10:12])/float(len(gene[10:12]))), (sum(gene[13:15])/float(len(gene[13:15])))]
    for target in targetGenes_FB:
        try:
            expvals_targets[targetGenes_FB[target]] = expvals[targetGenes_FB[target]]
            #expvals_targets_avg[target] = expvals_avg[targetGenes_FB[target]]
        except:
            pass

    #Now fitting expression values over time using a multiple variable linear regression in the statsmodels (sm) package.
    abd_target = {}
    abd_nonTarget = {}
#
#    abd_target_AL = []
#    abd_nonTarget_AL = []
#    abd_target_DR = []
#    abd_nonTarget_DR = []
#    for a in expvals_targets.keys():
#        abd_target_DR.extend(expvals_targets[a][0:4])
#        abd_target_AL.extend(expvals_targets[a][4:8])
#
#    for a in expvals.keys():
#        abd_nonTarget_DR.extend(expvals_targets[a][0:4])
#        abd_nonTarget_AL.extend(expvals_targets[a][4:8])
    for a in expvals_targets.keys():
        abd_target[a] = {}
        abd_target[a]['DR'] = expvals_targets[a][0:4]
        abd_target[a]['AL'] = expvals_targets[a][4:8]

    for a in expvals.keys():
        if a == 'NotFound':
            continue
        abd_nonTarget[a] = {}
        abd_nonTarget[a]['DR'] = expvals[a][0:4]
        abd_nonTarget[a]['AL'] = expvals[a][4:8]
#
#    time_target = []
#    time_target.extend([X]*len(abd_target_DR))
#    time_nonTarget = []
#    time_nonTarget.extend([X]*len(expvals.keys()))
#    gene_target = []
#    gene_nonTarget = []
#    for v in expvals_targets.keys():
#        gene_target.extend([v]*len(X))
#    for v in expvals.keys():
#        gene_nonTarget.extend([v]*len(X))
#    df_targets_DR = pd.DataFrame({'Day': abd_target, 'Abundance': abd_target_DR, 'Gene': gene_target})
#    df_targets_AL = pd.DataFrame({'Day': time_target, 'Abundance': abd_target_AL, 'Gene': gene_target})
#    df_nonTarget = pd.DataFrame({'Day': time_nonTarget, 'Abundance': abd_nonTarget, 'Gene': gene_nonTarget})

    LinVals = {}
    NonLinVals = {}
    #df_targets_range = pd.DataFrame(columns=['Day', 'Abundance', 'Gene'])
    X = [1, 7, 14, 21]
    for g in abd_target.keys():#range(0, len(df_targets_DR['Abundance']), 15):
        #Cycling through abundance data for target genes, range scaling abundance values at all time points for that gene to between 0 and 1 before fitting
        Y_DR = abd_target[g]['DR']#np.interp(abd_target[g]['DR'], (min(abd_target[g]['DR']), max(abd_target[g]['DR'])), (0, 1))
        Y_AL = abd_target[g]['AL']#np.interp(abd_target[g]['AL'], (min(abd_target[g]['AL']), max(abd_target[g]['AL'])), (0, 1))
        tempdf_DR = pd.DataFrame({'Day': X, 'Abundance': Y_DR, 'Gene': g})
        tempdf_AL = pd.DataFrame({'Day': X, 'Abundance': Y_AL, 'Gene': g})
        #df_targets_range = pd.concat([df_targets_range, tempdf], ignore_index=True)
        modelOLS_DR = sm.OLS(Y_DR, X).fit()
        modelOLS_AL = sm.OLS(Y_AL, X).fit()
        modelRLM_DR = sm.RLM(Y_DR, X).fit()
        modelRLM_AL = sm.RLM(Y_AL, X).fit() #RLM is robust=True, which accounts for outliers when fitting
        LinVals[RE.FBnumtoComName(g)] = {}
        LinVals[RE.FBnumtoComName(g)] = [modelOLS_DR.params[0], modelOLS_DR.pvalues[0], modelOLS_AL.params[0], modelOLS_AL.pvalues[0], modelRLM_DR.params[0], modelRLM_DR.pvalues[0], modelRLM_AL.params[0], modelRLM_AL.pvalues[0]]
        if (LinVals[RE.FBnumtoComName(g)][1] <= 0.05):
            fig, ax = plt.subplots()
            fig = sm.graphics.plot_fit(modelOLS_DR, 0, ax=ax)
            #fig = sm.graphics.plot_fit(modelOLS_DR, 0, ax=ax)
            ax.set_ylabel("Abundance")
            ax.set_xlabel("Day")
            ax.set_title("RLM Linear Regression of Gene on DR: "+ RE.FBnumtoComName(g))
            ax.set_title("OLS Linear Regression of Gene on DR: "+ RE.FBnumtoComName(g))
            if save == 1:
                fig.savefig(RE.FBnumtoComName(g)+'_ExpressionWithAge_DR.png', bbox_inches='tight')
        if (LinVals[RE.FBnumtoComName(g)][3] <= 0.05):
            fig, ax = plt.subplots()
            fig = sm.graphics.plot_fit(modelOLS_AL, 0, ax=ax)
            ax.set_ylabel("Abundance")
            ax.set_xlabel("Day")
            ax.set_title("OLS Linear Regression of Gene on AL: "+ RE.FBnumtoComName(g))
            if save == 1:
                fig.savefig(RE.FBnumtoComName(g)+'_ExpressionWithAge_AL.png', bbox_inches='tight')
    csv_columns = ['Gene','DR OLS Coef','DR OLS Pval', 'AL OLS Coef', 'AL OLS Pval', 'DR RLM Coef', 'DR RLM Pval', 'AL RLM Coef', 'AL RLM Pval']
    if save == 1:
        with open(geneListname+'_AgeRelatedExpressionTable.csv', "w", newline='') as f:
            w = csv.DictWriter(f, fieldnames=csv_columns)
            w.writeheader()
    #        for k,d in sorted(geneOutput_export.items()):
    #            w.writerow(mergedict({'Gene': k},d))
            for key, value in sorted(LinVals.items()):
                row = {csv_columns[0]: key, csv_columns[1]: value[0], csv_columns[2]: value[1], csv_columns[3]: value[2], csv_columns[4]: value[3], csv_columns[5]: value[4], csv_columns[6]: value[5], csv_columns[7]: value[6], csv_columns[8]: value[7]}
                w.writerow(row)
    return (LinVals)