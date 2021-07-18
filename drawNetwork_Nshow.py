# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 14:46:09 2019

@author: THilsabeck
"""
import RibotagExtract as RE
import geneOnto as GeneO
import pandas as pd
import numpy as np
from operator import itemgetter
import datetime
import os
import io
import sys
#from py2cytoscape.data.cyrest_client import CyRestClient
#from py2cytoscape import cyrest
import networkx as nx
#import graphviz
#from IPython.display import Image
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Circle
import matplotlib.patches as patches
import numpy as np
from matplotlib.font_manager import FontProperties
import statistics as stats
now = datetime.datetime.now()

def get_coordinates_in_circle(n, r): #Code from https://www.howtobuildsoftware.com/index.php/how-do/6zf/python-networkx-set-initial-positions-of-fixed-nodes-in-networkx-spring-graph
    thetas = [2*np.pi*(float(i)/n) for i in range(n)]
    returnlist = [(r*np.cos(theta),r*np.sin(theta)) for theta in thetas]
    return returnlist

def drawNetwork(inputGeneList, save=0, Nshow = 10, nwk_type='spring', nwk_dist=0):
    #Nshow (default=10) - integer or 'all'; number of nodes per GO term to show
    #nwk_type (default='geneOnto') - layout of nodes to use: geneOnto=run gene ontology and group by deepest sig GO term, spring=spring force layout, circle = concentric circle layout
    #nwk_dist (default=0) - type of nodes to include: 0=input genes and sig TFs, 1=0+primary known interactors, 2=0+1+known interactors of primary group
    TFout_sig = []
    TFNetwork_df = []
    riboGenes = []
    inputGenes = []
    if os.path.isfile(str(inputGeneList)):
        outputpath = inputGeneList[:-4]+"_geneAnalysis_"+now.strftime("%Y-%m-%d")
        with open(inputGeneList, 'r') as f:
            for line in f:
                if line.strip().split(',')[0] in ['Gene', 'TF', 'Target']:
                    continue
                else:
                    TF=line.strip().split(',')[0]
                    TFout_sig.append(RE.geneFBnum(TF))
                    riboGenes.append(RE.geneFBnum(TF))
                    for t in line.strip().split(','):
                        if t == TF:
                            continue
                        else:
                            riboGenes.append(RE.geneFBnum(t.lower().lstrip().replace('"', '')))
                            inputGenes.append(RE.geneFBnum(t.lower().lstrip().replace('"', '')))
                            TFNetwork_df.append((TF, RE.geneFBnum(TF), t.lstrip().replace('"', ''), RE.geneFBnum(t.lower().lstrip().replace('"', '')), '', str(1))) #TF gene name, TF gene FB#, target gene name, target gene FB#, Interaction type, gene source: TF=lightblue, Primary interaction = black, secondary interaction = green
    #os.mkdir(os.getcwd()+"\\"+outputpath)
        GENESET = inputGeneList[:-4]
        GENESETNAME = inputGeneList[:-4]
    elif (type(inputGeneList)==list):
        outputpath = "analysisGeneSet_geneAnalysis_"+now.strftime("%Y-%m-%d")
        #os.mkdir(os.getcwd()+"\\"+outputpath)
        with open("analysisGeneSet_FB"+now.strftime("%Y-%m-%d")+".txt", "w+") as f:
            for item in inputGeneList:
                f.write("%s\n" % RE.geneFBnum(item))
        GENESET = "analysisGeneSet_FB.txt"
        GENESETNAME = "analysisGeneSet_FB"
    else:
        print ("ERROR: Genes not given as filename or list.")
     
    #Now creating the network plot using skeleton code from http://jonathansoma.com/lede/algorithms-2017/classes/networks/networkx-graphs-from-source-target-dataframe/
    FBgeneInterFile = 'gene_genetic_interactions_fb_2019_04.txt' #from flybase for d6
    FBgeneInters=pd.read_csv(FBgeneInterFile,sep="\t")
    PriIntersKwn = []
    SecIntersKwn = []
    for gene in riboGenes:#ribo_out.keys():
        for i in list(np.where(FBgeneInters['Starting_gene(s)_FBgn']==RE.geneFBnum(gene))[0]):
            PriIntersKwn.append((FBgeneInters['Starting_gene(s)_symbol'][i], FBgeneInters['Starting_gene(s)_FBgn'][i], FBgeneInters['Interacting_gene(s)_symbol'][i], FBgeneInters['Interacting_gene(s)_FBgn'][i], FBgeneInters['Interaction_type'][i], str(2)))
            for s in list(np.where(FBgeneInters['Starting_gene(s)_FBgn']==FBgeneInters['Interacting_gene(s)_FBgn'][i])[0]):
                SecIntersKwn.append((FBgeneInters['Starting_gene(s)_symbol'][s], FBgeneInters['Starting_gene(s)_FBgn'][s], FBgeneInters['Interacting_gene(s)_symbol'][s], FBgeneInters['Interacting_gene(s)_FBgn'][s], FBgeneInters['Interaction_type'][s], str(3)))
    if nwk_dist == 0:
        next
    elif nwk_dist == 1:
        for i in PriIntersKwn:
            riboGenes.append(i[2])
    elif nwk_dist == 2:
        for i in PriIntersKwn:
            riboGenes.append(i[2])
        for i in SecIntersKwn:
            riboGenes.append(i[2])
    riboGenes = list(dict.fromkeys(riboGenes)) #Removing duplicates
    ribo_out = RE.riboext(riboGenes, GENESETNAME+'_network', save)#, outputpath)
    ribo_out.pop('Comparison Stats', None)
    #LinVals = GA.geneExp_Age(riboGenes, GENESETNAME, save)
    GO_fdr = GeneO.geneOnto(riboGenes, GENESETNAME, save)
    #Creating colormap for TF analysis genes and their target genes based on tissue with highest FoldChange
    largestFoldChange = {}
    for node in ribo_out:
        try:
            values = []
            for tissue in ['Fatbody', 'Germ', 'Gut', 'Heart', 'Muscle', 'Neuron', 'Tubule']:
                if ribo_out[node][tissue+ ' FoldChange'] == 0:
                    values.append(0.999)
                else:
                    try:
                        values.append(ribo_out[node][tissue+ ' FoldChange'])
                    except:
                        values.append(0.999)
            values = [abs(x-1) for x in values] #Determing distance from 1 (no change) so max() will pick value in either direction
            if values.index(max(values)) == 0:
                largestFoldChange.update({RE.geneFBnum(node):'gold'}) #For Fatbody
            elif values.index(max(values)) == 1:
                largestFoldChange.update({RE.geneFBnum(node):'green'}) #For Germ
            elif values.index(max(values)) == 2:
                largestFoldChange.update({RE.geneFBnum(node):'brown'}) #For Gut
            elif values.index(max(values)) == 3:
                largestFoldChange.update({RE.geneFBnum(node):'red'}) #For Heart
            elif values.index(max(values)) == 4:
                largestFoldChange.update({RE.geneFBnum(node):'silver'}) #For Muscle
            elif values.index(max(values)) == 5:
                largestFoldChange.update({RE.geneFBnum(node):'blue'}) #For Neuron
            elif values.index(max(values)) == 6:
                largestFoldChange.update({RE.geneFBnum(node):'magenta'}) #For Tubule
        except:
            continue
    #Adding primary and secondary interactions to TFnetwork list
    if nwk_dist == 0:
        TFNetwork_df.extend([i for i in PriIntersKwn if i[3] in riboGenes])
    elif nwk_dist == 1:
        TFNetwork_df.extend(PriIntersKwn)
    elif nwk_dist == 2:
        TFNetwork_df.extend(PriIntersKwn)
        TFNetwork_df.extend(SecIntersKwn)
    TFNetwork_df = [x for x in TFNetwork_df if x[1] != 'NotFound']
    TFNetwork_df = [x for x in TFNetwork_df if x[3] != 'NotFound']
    #Adding largest tissue foldchange for each gene: neuron=0, fatbody=1, germ=2,gut=3,heart=4,muscle=5,tubule=6
    interdf = pd.DataFrame(TFNetwork_df, columns=['source', 'sourceFB', 'target', 'targetFB', 'interactionType', 'geneSource'])#, 'largest tissue foldchange'])#pd.read_csv('testnetwork.txt', sep='\t')
    Pridf = pd.DataFrame(PriIntersKwn, columns=['source', 'sourceFB', 'target', 'targetFB', 'interactionType', 'geneSource'])#, 'largest tissue foldchange'])#pd.read_csv('testnetwork.txt', sep='\t')
    Secdf = pd.DataFrame(SecIntersKwn, columns=['source', 'sourceFB', 'target', 'targetFB', 'interactionType', 'geneSource'])#, 'largest tissue foldchange'])#pd.read_csv('testnetwork.txt', sep='\t')
    inputdf = pd.DataFrame(inputGenes, columns=['Input_gene'])
    inputdf['Input_geneFB'] = [RE.geneFBnum(gene) for gene in inputGenes]
    inputdf = pd.DataFrame([x for x in inputdf['Input_geneFB'] if x != 'NotFound'], columns = ['Input_geneFB'])
    mapping_labels = dict()
    for s in TFNetwork_df:
        mapping_labels.update({s[1]:s[0]})
        mapping_labels.update({s[3]:s[2]})
    #sources = list(interdf.source)#.unique())
    #targets = list(interdf.target)#.unique())
    if len(TFNetwork_df) < 100:
        figsize = len(TFNetwork_df)*2
    elif len(TFNetwork_df) > 1000:
        figsize = 300
    else:
        figsize = 250
    fig =plt.figure(figsize=(figsize,figsize))
    ax = plt.gca()
    # 1. Create the graph
    g = nx.from_pandas_edgelist(interdf, source='sourceFB', target='targetFB')

    # 2. Create a layout for our nodes
    #Code from swatchai on Stackoverflow.com: https://stackoverflow.com/questions/55750436/group-nodes-together-in-networkx
    # All this code should replace original `pos=nx.spring_layout(graph)`
    #colors = ['gold', 'green', 'brown', 'red', 'silver', 'blue', 'magenta']
#    nodes_by_color = {}
#    for val in colors:
#        nodes_by_color[val] = {}
#        fnode = 1
#        for node in g:
#            try: #Need to catch errors because gene names are still not consistent between ribotag dataset and FB name dataset, and ribotag is not complete, either.
#                if largestFoldChange[node] == val:
#                    if fnode == 1:
#                        nodes_by_color[val] = [node]
#                        fnode = 0
#                    else:
#                        nodes_by_color[val].append(node)
#            except:
#                continue

    nodes_by_connectivity = {}
    for val in g._node:
        nodes_by_connectivity[val] = g.degree(val)
    sortedConnectivity=sorted(nodes_by_connectivity, key=nodes_by_connectivity.__getitem__)
    GOterms = [(GO[1], GO[2]) for GO in GO_fdr]
    nodes_by_GO = {GO[1]: GO[4] for GO in GO_fdr}
    if nwk_type.lower() == 'circle':
        layout = nx.circular_layout(g)
        radii = [30,45,60,75,90,105,120]#[:num_radii]  # for concentric circles
    #Setting nodes based on connectivity
        NinRadii = {r: 0 for r in radii}
        for ea in layout.keys():
            new_r = 1
            if g.degree(ea) <= g.degree(sortedConnectivity[round(len(sortedConnectivity)*(1/7))]):
                new_r = radii[6]
                NinRadii[radii[6]] += 1
            elif g.degree(sortedConnectivity[round(len(sortedConnectivity)*(1/7))]) < g.degree(ea) <= g.degree(sortedConnectivity[round(len(sortedConnectivity)*(2/7))]):
                new_r = radii[5]
                NinRadii[radii[5]] += 1
            elif g.degree(sortedConnectivity[round(len(sortedConnectivity)*(2/7))]) < g.degree(ea) <= g.degree(sortedConnectivity[round(len(sortedConnectivity)*(3/7))]):
                new_r = radii[4]
                NinRadii[radii[4]] += 1
            elif g.degree(sortedConnectivity[round(len(sortedConnectivity)*(3/7))]) < g.degree(ea) <= g.degree(sortedConnectivity[round(len(sortedConnectivity)*(4/7))]):
                new_r = radii[3]
                NinRadii[radii[3]] += 1
            elif g.degree(sortedConnectivity[round(len(sortedConnectivity)*(4/7))]) < g.degree(ea) <= g.degree(sortedConnectivity[round(len(sortedConnectivity)*(5/7))]):
                new_r = radii[2]
                NinRadii[radii[2]] += 1
            elif g.degree(sortedConnectivity[round(len(sortedConnectivity)*(5/7))]) < g.degree(ea) <= g.degree(sortedConnectivity[round(len(sortedConnectivity)*(6/7))]):
                new_r = radii[1]
                NinRadii[radii[1]] += 1
            elif g.degree(sortedConnectivity[round(len(sortedConnectivity)*(6/7))]) < g.degree(ea):
                new_r = radii[0]
                NinRadii[radii[0]] += 1
            layout[ea][0] = new_r# *= new_r   # relayoutition nodes as concentric circles

        #Setting outter ring to last empty ring
        empty = []
        for r in NinRadii.keys():
            if NinRadii[r] == 0:
                empty.append(r)
        for r in empty:
            NinRadii.pop(r)
        for ea in layout.keys():
            if layout[ea][0] == list(NinRadii.keys())[-1]:
                layout[ea][0] = list(NinRadii.keys())[-2]+15
        NinRadii[list(NinRadii.keys())[-2]+15] = NinRadii[list(NinRadii.keys())[-1]]
        #NinRadii.pop(120)
        #NinRadii[list(NinRadii.keys())[-2]+15] = NinRadii.pop(list(NinRadii.keys())[-1])


    #Evenly spacing nodes in each circle
        for r in NinRadii.keys():
            if r == 30 and NinRadii[r] <= 6:
                circular_positions = get_coordinates_in_circle(NinRadii[r], 6)
                for i,n in enumerate([n for n in layout.keys() if layout[n][0] == r]):
                    layout[n] = circular_positions[i]
            else:
                circular_positions = get_coordinates_in_circle(NinRadii[r], r)
                for i,n in enumerate([n for n in layout.keys() if layout[n][0] == r]):
                    layout[n] = circular_positions[i]

    #Using spring_layout and attempted to set fixed nodes using top 1/7 connected genes
    elif nwk_type.lower() == 'spring':
        fixed_nodes = [n for n in g.nodes() if g.degree(n) > g.degree(sortedConnectivity[round(len(sortedConnectivity)*(6/7))])]
        layout = nx.spring_layout(g)
        circular_positions = get_coordinates_in_circle(len(fixed_nodes), 30)
        i = 0
        for p in layout.keys():
            if p in fixed_nodes:
                layout[p] = circular_positions[i]
                i += 1
        layout = nx.drawing.layout.kamada_kawai_layout(g)#nx.nx_agraph.graphviz_layout(g, prog='twopi')#nx_agraph.pygraphviz_layout(g, prog='neato')#spectral_layout(g)#spring_layout(g,k=2,iterations=100)
    
    #Placing nodes based on deepest gene ontology group
    elif nwk_type.lower() == 'geneonto':
        #Below code modified from user swatchai at https://stackoverflow.com/questions/55750436/group-nodes-together-in-networkx
        # --- Begin_myhack ---
        # All this code should replace original `pos=nx.spring_layout(graph)`
        layout = nx.circular_layout(g)   # replaces your original pos=...
        # prep center points (along circle perimeter) for the clusters
        GOgrpPos = {}
        for ea in layout.keys():
            GOgrpPos[ea] = {}
            GOgrpPos[ea]['deepestGO'] = 'Not in Significant GO'
            GOgrpPos[ea]['depth'] = 0
            GOgrpPos[ea]['posx'] = 0
            for GO in GOterms:
                if ea in nodes_by_GO[GO[0]]:
                    if GO[1] >= GOgrpPos[ea]['depth']:
                        GOgrpPos[ea]['deepestGO'] = GO[0]
                        GOgrpPos[ea]['depth'] = GO[1]
                else:
                    GOgrpPos[ea]['posx']+=1
        GOdeepestTF = list(set([GOgrpPos[ea]['deepestGO'] for ea in layout.keys() if ea in TFout_sig]))
        GOdeepest = list(set([GOgrpPos[ea]['deepestGO'] for ea in layout.keys()]))
        GOdeepest = [GO for GO in GOdeepest if GO not in GOdeepestTF]
        GOdeepestN = {}
        GOdeepestList = {}
        for goterm in GOgrpPos:
            if GOgrpPos[goterm]['deepestGO'] in GOdeepestN:
                GOdeepestN[GOgrpPos[goterm]['deepestGO']] += 1
                GOdeepestList[GOgrpPos[goterm]['deepestGO']].add(goterm)
            else:
                GOdeepestN[GOgrpPos[goterm]['deepestGO']] = 1
                GOdeepestList[GOgrpPos[goterm]['deepestGO']] = set([goterm])
        GO_angs = np.linspace(0, 2*np.pi, 1+len(GOdeepest))
        GO_angsTF = np.linspace(0, 2*np.pi, 1+len(GOdeepestTF))
        Node_angs = {}
        NodeRepos = {}
        GOdeepestNshow = {}
        GOrad = 0 #2.5*figsize/(len(GOdeepest)*np.pi)#3.5     # radius of non TF GO circle
        for GO in GOdeepest:
            if type(Nshow) == str:
                Nshow = GOdeepestN[GO]
            Node_angs[GO] = {}
            # if GOdeepestN[GO] < 15:
            Node_angs[GO]['positions'] = [np.linspace(0, 2*np.pi, 1+Nshow)]#GOdeepestN[GO])]
            if GOdeepestN[GO] < Nshow:
                Node_angs[GO]['radii'] = [GOdeepestN[GO]/(3.5*np.pi)] #Setting radii for each GO circle
            else:
                Node_angs[GO]['radii'] = [Nshow/(3.5*np.pi)] #Setting radii for each GO circle
            # else:
            #     count = 0
            #     for c in [10, 11, 12, 13, 14, 15, 16]:#range(0, GOdeepestN[GO], 15):
            #         if count == 0:
            #             Node_angs[GO]['positions'] = [np.linspace(0, 2*np.pi, c)]
            #             Node_angs[GO]['radii'] = [c/(3.5*np.pi)]
            #         else:
            #             Node_angs[GO]['positions'].append(np.linspace(0, 2*np.pi, c))
            #             Node_angs[GO]['radii'].append((c)/(3.5*np.pi)) #Setting radii for each GO circle
            #         count += 15
            NodeRepos[GO] = []
            GOconnectivity = { key:value for key,value in nodes_by_connectivity.items() if key in GOdeepestList[GO]}
            GOdeepestNshow[GO] = sorted(GOconnectivity, key=GOconnectivity.__getitem__, reverse=True)[:Nshow]
            Node_angs[GO]['radii'] = [sum([nodes_by_connectivity[fb] for fb in GOdeepestNshow[GO]])/(2*np.pi)] #Setting radii for each GO circle based on size/connectivity of included nodes
            GOrad += 2*Node_angs[GO]['radii'][0]
        GOradTF = 0 # radius of whole TF GO circle
        for GO in GOdeepestTF:
            if type(Nshow) == str:
                Nshow = GOdeepestN[GO]
            Node_angs[GO] = {}
            # if GOdeepestN[GO] < 15:
            Node_angs[GO]['positions'] = [np.linspace(0, 2*np.pi, 1+Nshow)]#GOdeepestN[GO])]
            NodeRepos[GO] = []
            GOconnectivity = { key:value for key,value in nodes_by_connectivity.items() if key in GOdeepestList[GO]}
            GOdeepestNshow[GO] = sorted(GOconnectivity, key=GOconnectivity.__getitem__, reverse=True)[:Nshow]
            Node_angs[GO]['radii'] = [sum([nodes_by_connectivity[fb] for fb in GOdeepestNshow[GO]])/(10*2*np.pi)] #Setting radii for each GO circle based on size/connectivity of included nodes
            GOradTF += 10*Node_angs[GO]['radii'][0]
        GOdeepestNshowList = [item for sublist in GOdeepestNshow for item in GOdeepestNshow[sublist]]
        GOrad = GOrad/(2*np.pi)
        GOradTF = GOradTF/(len(GOdeepestTF)*2*np.pi)
        if GOrad <= GOradTF:
            GOrad = GOradTF+5
        GOrepos = {}
        GOrepos['pos'] = []
        GOreposTF = {}
        GOreposTF['pos'] = []
        #***Need to make sure each ring is centered on the right GO cluster position
        for ea in GO_angs: #Setting GO cluster positions
            if ea > 0:
                GOrepos['pos'].append(np.array([GOrad*np.cos(ea), GOrad*np.sin(ea)])) # location of each cluster
        for ea in GO_angsTF: #Setting GO cluster positions
            if ea > 0:
                GOreposTF['pos'].append(np.array([GOradTF*np.cos(ea), GOradTF*np.sin(ea)])) # location of each cluster
        for got in GOdeepest: #Setting node positions in GO clusters
            for r, vals in enumerate(Node_angs[got]['positions']):
                    for ea in vals:
                        if ea > 0:
                            NodeRepos[got].append(np.array([GOrepos['pos'][GOdeepest.index(got)][0] + Node_angs[got]['radii'][r]*np.cos(ea), GOrepos['pos'][GOdeepest.index(got)][1] + Node_angs[got]['radii'][r]*np.sin(ea)]))
        for got in GOdeepestTF: #Setting node positions in GO clusters
            for r, vals in enumerate(Node_angs[got]['positions']):
                    for ea in vals:
                        if ea > 0:
                            NodeRepos[got].append(np.array([GOreposTF['pos'][GOdeepestTF.index(got)][0] + Node_angs[got]['radii'][r]*np.cos(ea), GOreposTF['pos'][GOdeepestTF.index(got)][1] + Node_angs[got]['radii'][r]*np.sin(ea)]))
        #Code to check layout of nodes:
            # for i in GOrepos['pos']:
            #     plt.scatter(i[0], i[1])
            # for j in NodeRepos:
            #     for g in NodeRepos[j]:
            #         plt.scatter(g[0], g[1])
        for i in GOdeepestList:
            for ea in GOdeepestNshow[i]:
                layout[ea] = NodeRepos[i][GOdeepestNshow[i].index(ea)]#GOgrpPos[ea]['deepestGO']][list(GOdeepestNshow[GOgrpPos[ea]['deepestGO']]).index(ea)]
        # --- End_myhack ---
        #Removing genes not in Nshow list from layout:
        remove = [n for n in layout.keys() if n not in GOdeepestNshowList]
        for key in remove:
            layout.pop(key, None)
            g.remove_node(key)
            mapping_labels.pop(key, None)
            # interdf
            # inputdf
            # Pridf
            # Secdf
        
    # 3. Draw the parts we want
    nx.draw_networkx_edges(g, layout, width=4, edge_color='#AAAAAA')
    #Highlighting input genes
    Inputgenes = [node for node in g.nodes() if node in inputdf.Input_geneFB.unique() and node in GOdeepestNshowList and node not in TFout_sig]
    Inputsize = [g.degree(node) * 1000 for node in g.nodes() if node in Inputgenes and node in GOdeepestNshowList]
    Inputnodes = nx.draw_networkx_nodes(g, layout, nodelist=Inputgenes, node_size=Inputsize, node_color='yellow', linewidths=5)#AAAAAA'
    Inputnodes.set_edgecolor(pd.Series([largestFoldChange[key] for key in largestFoldChange.keys() if key in Inputgenes]))
    #Highlighting TFs
    TFsources = [node for node in g.nodes() if node in interdf.sourceFB[interdf.geneSource == "1"].unique() and node in GOdeepestNshowList]
    TFsize = [g.degree(node) * 1000 for node in g.nodes() if node in TFsources and node in GOdeepestNshowList]
    #Source_color_map = [interdf.geneSource[interdf.sourceFB == node] for node in TFsources]
    TFnodes = nx.draw_networkx_nodes(g, layout, nodelist=TFsources, node_size=TFsize, node_color='darkturquoise', linewidths=5)
    TFnodes.set_edgecolor(pd.Series([largestFoldChange[key] for key in largestFoldChange.keys() if key in TFsources and key in GOdeepestNshowList]))
    #nx.draw_networkx_labels(g, pos=nx.kamada_kawai_layout(g), labels=interdf.source.unique())
    TFtargets = [node for node in g.nodes() if node in interdf.targetFB[interdf.geneSource == "1"].unique() and node in GOdeepestNshowList]
    Targetsize = [g.degree(node) * 1000 for node in g.nodes() if node in TFtargets and node in GOdeepestNshowList]
    #Target_color_map = [interdf.geneSource[interdf.targetFB == node] for node in TFtargets]
    Targetnodes = nx.draw_networkx_nodes(g, layout, nodelist=TFtargets, node_size=Targetsize, node_color='yellow', linewidths=5)#AAAAAA'
    Targetnodes.set_edgecolor(pd.Series([largestFoldChange[key] for key in largestFoldChange.keys() if key in TFtargets and key in GOdeepestNshowList]))
    #nx.draw_networkx_labels(g, pos=nx.kamada_kawai_layout(g), labels=interdf.target.unique())
    if nwk_dist == 1:
        Pritargets = [node for node in g.nodes() if node in interdf.targetFB[interdf.geneSource == "2"].unique() and node in GOdeepestNshowList]#Pridf.sourceFB]
        Prisize = [g.degree(node) * 500 for node in g.nodes() if node in Pritargets and node in GOdeepestNshowList]
        Prinodes = nx.draw_networkx_nodes(g, layout, nodelist=Pritargets, node_size=Prisize, node_color='springgreen', linewidths=5)#, linewidths=5)
        Prinodes.set_edgecolor(pd.Series([largestFoldChange[key] for key in largestFoldChange.keys() if key in Pritargets and key in GOdeepestNshowList]))
        #Prinodes.set_edgecolor('blue')
    elif nwk_dist == 2:
        Pritargets = [node for node in g.nodes() if node in interdf.targetFB[interdf.geneSource == "2"].unique() and node in GOdeepestNshowList]#Pridf.sourceFB]
        Prisize = [g.degree(node) * 500 for node in g.nodes() if node in Pritargets and node in GOdeepestNshowList]
        Prinodes = nx.draw_networkx_nodes(g, layout, nodelist=Pritargets, node_size=Prisize, node_color='springgreen', linewidths=5)#, linewidths=5)
        Prinodes.set_edgecolor(pd.Series([largestFoldChange[key] for key in largestFoldChange.keys() if key in Pritargets]))
        #Prinodes.set_edgecolor('blue')
        Sectargets = [node for node in g.nodes() if node in interdf.targetFB[interdf.geneSource == "3"].unique() and node in GOdeepestNshowList]#Secdf.sourceFB]
        Secsize = [g.degree(node) * 250 for node in g.nodes() if node in Sectargets and node in GOdeepestNshowList]
        Secnodes = nx.draw_networkx_nodes(g, layout, nodelist=Sectargets, node_size=Secsize, node_color='pink', linewidths=5)
        Secnodes.set_edgecolor(pd.Series([largestFoldChange[key] for key in largestFoldChange.keys() if key in Sectargets and key in GOdeepestNshowList]))
        #Secnodes.set_edgecolor('purple')
    enhancers = [(source, target) for source in interdf.sourceFB[interdf.interactionType=='enhanceable'] for target in interdf.targetFB[interdf.interactionType=='enhanceable']]
    enhancers = [e for e in enhancers if e in g.edges()]
    nx.draw_networkx_edges(g, layout, edgelist=enhancers, arrows=True, width=2, edge_color='blue')
    suppressors = [(source, target) for source in interdf.sourceFB[interdf.interactionType=='suppressible'] for target in interdf.targetFB[interdf.interactionType=='suppressible']]
    suppressors = [s for s in suppressors if s in g.edges()]
    nx.draw_networkx_edges(g, layout, edgelist=suppressors, arrows=True, width=2, edge_color='red')
    nx.draw_networkx_labels(g, layout, labels=mapping_labels, font_size=0.6*figsize)
#    high_degree_targets = [node for node in g.nodes() if node in interdf.targetFB.unique() and g.degree(node) > 50]
#    nx.draw_networkx_nodes(g, layout, nodelist=high_degree_targets, node_size=size, node_color='#fc8d62')
    #source_dict = dict(zip(TFsources, TFsources))
    #target_dict = dict(zip(TFtargets, TFtargets))
    #nx.draw_networkx_labels(g, layout, labels=source_dict)
    #nx.draw_networkx_labels(g, layout, labels=target_dict)

    # 4. Turn off the axis because I know you don't want it
    #Below circling and annotation code modified from farenorth on StackOverflow (https://stackoverflow.com/questions/37489874/how-do-i-put-a-circle-with-annotation-in-matplotlib)
    if nwk_type.lower() == 'geneonto':
        font = FontProperties()
        font.set_weight('bold')
        font.set_size(0.6*figsize)
        for i, GO in enumerate(GOdeepest):
            if GO == 'Not in Significant GO':
                circle_rad = Node_angs[GO]['radii'][-1]*1.35  # This is the radius, in points
                #ax.plot(repos[GOdeepest.index(GO)], 'o',
                #        ms=circle_rad * 2, mec='b', mfc='none', mew=2)
        ####    #Change radius (outter bound) for circle to just beyond last node position for circle
                circle = plt.Circle(GOrepos['pos'][GOdeepest.index(GO)], radius=circle_rad,lw=10.,ec='r',fill=False)
                ax.add_artist(circle)
                circle.set_clip_box(ax.bbox)
                # xyt=(60,60)
                # if GO_angs[i] >= np.pi:
                #     circle_rad = -circle_rad
                #     xyt = (0,0)
                plt.annotate(GO, xy=(GOrepos['pos'][GOdeepest.index(GO)][0]+circle_rad, GOrepos['pos'][GOdeepest.index(GO)][1]),
                             xytext = (60,60),
                            textcoords='offset points',
                            color='r', size=0.6*figsize,#'large',
                            fontproperties=font,
                            #horizontalalignment='left',
                            arrowprops=dict(
                                arrowstyle='simple,tail_width=0.3,head_width=0.8,head_length=0.8',
                                facecolor='r', shrinkB=circle_rad * 2)
                )
            else:
                circle_rad = Node_angs[GO]['radii'][-1]*1.35  # This is the radius, in points
                #ax.plot(repos[GOdeepest.index(GO)], 'o',
                #        ms=circle_rad * 2, mec='b', mfc='none', mew=2)
        ####    #Change radius (outter bound) for circle to just beyond last node position for circle
                circle = plt.Circle(GOrepos['pos'][GOdeepest.index(GO)], radius=circle_rad,lw=10.,ec='b',fill=False)
                ax.add_artist(circle)
                circle.set_clip_box(ax.bbox)
                # xyt=(60,60)
                # if GO_angs[i] >= np.pi:
                #     circle_rad = -circle_rad
                #     xyt=(0,0)
                plt.annotate(GO, xy=(GOrepos['pos'][GOdeepest.index(GO)][0]+circle_rad, GOrepos['pos'][GOdeepest.index(GO)][1]),
                             xytext = (60,60),
                            textcoords='offset points',
                            color='b', size=0.6*figsize,#'large',
                            fontproperties=font,
                            #horizontalalignment='left',
                            arrowprops=dict(
                                arrowstyle='simple,tail_width=0.3,head_width=0.8,head_length=0.8',
                                facecolor='b', shrinkB=circle_rad * 2)
                )
        for i, GO in enumerate(GOdeepestTF):
            circle_rad = Node_angs[GO]['radii'][-1]*1.35  # This is the radius, in points
            #ax.plot(repos[GOdeepest.index(GO)], 'o',
            #        ms=circle_rad * 2, mec='b', mfc='none', mew=2)
    ####    #Change radius (outter bound) for circle to just beyond last node position for circle
            circle = plt.Circle(GOreposTF['pos'][GOdeepestTF.index(GO)], radius=circle_rad,lw=10.,ec='r',fill=False)
            ax.add_artist(circle)
            circle.set_clip_box(ax.bbox)
            # xyt=(60,60)
            # if GO_angs[i] >= np.pi:
            #     circle_rad = -circle_rad
            #     xyt=(0,0)
            plt.annotate(GO, xy=(GOreposTF['pos'][GOdeepestTF.index(GO)][0]+circle_rad, GOreposTF['pos'][GOdeepestTF.index(GO)][1]),
                         xytext = (60,60),
                        textcoords='offset points',
                        color='b', size=0.6*figsize,#'large',
                        fontproperties=font,
                        #horizontalalignment='left',
                        arrowprops=dict(
                            arrowstyle='simple,tail_width=0.3,head_width=0.8,head_length=0.8',
                            facecolor='b', shrinkB=circle_rad * 2)
            )
                
        # labels = GOdeepest

        # # image = mpimg.imread("Figure 2020-05-07 150300.png") # just a image of your graph
        # # plt.imshow(image)
        # ax = plt.gca()
        
        # # set your own radius and centers of circles in loop, like here
        # r = 11; c = (157,177)
        # circ1 = patches.Circle(c,2*r,lw=3.,ec='b',fill=False)
        # ax.add_artist(circ1)
        # circ1.set_clip_box(ax.bbox)
        
        # # annotate circles
        # # I have one circle but for your array pos_annotation_node
        # # you need 'i' to extract proper position
        # for i,label in enumerate(labels):
        #     annot_array_end = (c[0], c[1]+2*(-1)**i*r)
        #     annot_text_pos = (c[0]+3,c[1]+5*(-1)**i*r)
        #     ax.annotate(label, 
        #       xy= annot_array_end, 
        #       xytext=annot_text_pos,
        #       color='b',
        #       fontproperties=font,
        #       arrowprops=dict(fc='b', shrink=.005)
        #     )
        
        # plt.show()
    plt.axis('off')

    plt.title("Regulatory Network of "+GENESETNAME+" Input Genes", fontsize = figsize)

    # 5. Tell matplotlib to show it
    plt.show()
    if save == 1:
        if nwk_dist == 0:
            fig.savefig(GENESETNAME+'_RegulatoryNetwork_Basic'+now.strftime("%Y-%m-%d")+'.png')
        if nwk_dist == 1:
            fig.savefig(GENESETNAME+'_RegulatoryNetwork_Primary'+now.strftime("%Y-%m-%d")+'.png')
        if nwk_dist == 2:
            fig.savefig(GENESETNAME+'_RegulatoryNetwork_PriAndSec'+now.strftime("%Y-%m-%d")+'.png')
#