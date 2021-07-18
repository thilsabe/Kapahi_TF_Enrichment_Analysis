# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 16:13:29 2019

@author: THilsabeck
"""
#import sys
import GSE_benhoch as GSE #_path as GSE#_combscript as GSE#
import RibotagExtract as RE
import geneExp_withAge_Diet as GA
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
import statistics as stats
now = datetime.datetime.now()

#For when I need to run from linux terminal
#first_arg = sys.argv[1]
#second_arg = sys.argv[2]

def get_coordinates_in_circle(n, r): #Code from https://www.howtobuildsoftware.com/index.php/how-do/6zf/python-networkx-set-initial-positions-of-fixed-nodes-in-networkx-spring-graph
    thetas = [2*np.pi*(float(i)/n) for i in range(n)]
    returnlist = [(r*np.cos(theta),r*np.sin(theta)) for theta in thetas]
    return returnlist

def genePipeline(inputGeneList, save=0, TF_FDR=0.05, MultiTest='Bon', EXECMODE='NODEB', PERMUTATIONS=1000, TFCHIPSEQ='lbl.AllTFchIPseq_FB.txt', STARTPOS=2, ENDPOS=3, GENEPOS=15, TFPOS=19, ANNOTPOS=7, SCOREPOS=26, BACKGROUND='All', GENOMIC_FEATURE='ALL', nTF=1, nwk_type='circle', nwk_dist=0):
    #Diverting print outputs
    text_trap = io.StringIO()
    sys.stdout = text_trap
    #To restore output: sys.stdout = sys.__stdout__
    #MultiTest= 'Bon' Bonferroni (default), 'BH' Benjamini-Hochberg, 'none' selects any with p-val-score <= 0.05, and 'PK' which takes top 20; for selecting TFs from TF analysis
    #save = 0 (no, default) or 1 (yes). Specifies whether or not files should be saved.
    #nTF = 1 (start at first TF, default) or other integer number specifying which TF to start analysis on. Make sure to rename original TFanalysis_out file, otherwise will be overwritten!
    #nwk_dist = 0 (only include TFs and their targets in network diagram, default) or 1 (include primary known interactors) or 2 (include primary and secondary interactors, which are primary interactors of the initial primary interactors)
    #nwk_type = 'circle' (default) (nodes setup in concentric circles based on connectivity, default) or 'spring' (use a spring-weighted layout for the network nodes)
    if (MultiTest == 'BH'):
        MTC=-3
    elif (MultiTest == 'Bon'):
        MTC=-1
    elif (MultiTest == 'none'):
        MTC=-5
    elif (MultiTest == 'PK'):
        MTC=0
    else:
        print ("Incorrect multiple testing correction provided")
    #EXECMODE=NODEB or DEB
    #BACKGROUND=All (default), Fatbody, Germline, Gut, Heart, Muscle, Neuron, or Tubule
    BACKGROUNDNAME=BACKGROUND
    BACKGROUND = BACKGROUND+"_Background_FB.txt"
    #nwk_dist = 0, 1 (default), or 2. If 0, will not extend network, 1 extends to first known interactors, 2 extends to interactors of first interactors
    #GENOMIC_FEATURE=ALL OR PROMOTER
    if os.path.isfile(str(inputGeneList)):
        outputpath = inputGeneList[:-4]+"_geneAnalysis_"+now.strftime("%Y-%m-%d")
        with open(inputGeneList, 'r') as f:
            with open(inputGeneList[:-4]+'_FB.txt', 'w+') as w:
                for line in f:
                    print (line)
                    if line.startswith('Gene'):
                        continue
                    else:
                        w.write("%s\n" % RE.geneFBnum(line.strip()))
        #os.mkdir(os.getcwd()+"\\"+outputpath)
        OUTFILE = inputGeneList[:-4]+'_FB.out'
        LOGFILE = inputGeneList[:-4]+'_FB.log'
        GENESET = inputGeneList[:-4]+'_FB.txt'
        GENESETNAME = inputGeneList[:-4]+'_FB'
    elif (type(inputGeneList)==list):
        outputpath = "analysisGeneSet_geneAnalysis_"+now.strftime("%Y-%m-%d")
        #os.mkdir(os.getcwd()+"\\"+outputpath)
        with open("analysisGeneSet_FB.txt", "w+") as f:
            for item in inputGeneList:
                f.write("%s\n" % RE.geneFBnum(item))
        GENESET = "analysisGeneSet_FB.txt"
        GENESETNAME = "analysisGeneSet_FB"
        OUTFILE = GENESET[:-4]+'.out'
        LOGFILE = GENESET[:-4]+'.log'
    else:
        print ("ERROR: Genes not given as filename or list.")
    #Creating and writing GSE.setup file for TF analysis
    GSEsetup = [EXECMODE, PERMUTATIONS, GENESET, GENESETNAME, TFCHIPSEQ, STARTPOS, ENDPOS, GENEPOS, TFPOS, ANNOTPOS,SCOREPOS, BACKGROUND, BACKGROUNDNAME, GENOMIC_FEATURE, OUTFILE, LOGFILE]
    with open("GSE.setup", "w+") as f:
        fieldN = 0
        for field in ['EXECMODE', 'PERMUTATIONS', 'GENESET', 'GENESETNAME', 'TFCHIPSEQ', 'STARTPOS', 'ENDPOS', 'GENEPOS', 'TFPOS', 'ANNOTPOS','SCOREPOS', 'BACKGROUND', 'BACKGROUNDNAME', 'GENOMIC_FEATURE', 'OUTFILE', 'LOGFILE']:
            if (type(GSEsetup[fieldN])==str):
                f.write(field+"\t"+GSEsetup[fieldN]+"\n")
            else:
                f.write(field+"\t"+str(GSEsetup[fieldN])+"\n")
            fieldN+=1

    TFanalysis_out = GSE.GSEanalysis(TF_FDR, nTF, outputpath)
    TFout_sig = []
    TFNetwork_df = []
    if (MultiTest == 'BH' or 'Bon'): #Default
        for nested_list in TFanalysis_out[8:]:
            if (str(nested_list[MTC]).lower()=='true'):
                TFout_sig.append(nested_list[0])
                for t in nested_list[-5].split(','):
                    TFNetwork_df.append((nested_list[0], RE.geneFBnum(nested_list[0]), t.lstrip().replace('"', ''), RE.geneFBnum(t.lower().lstrip().replace('"', '')), '', str(1))) #TF gene name, TF gene FB#, target gene name, target gene FB#, Interaction type, gene source: TF=lightblue, Primary interaction = black, secondary interaction = green
                if (nested_list[MTC-1]==0):
                    print ("Warning: "+nested_list[0]+"p-value is 0, suggesting something is off. Likely too few permutations.")
    elif (MultiTest == 'none'):
        for nested_list in TFanalysis_out[8:]:
            if (nested_list[MTC]<=0.05):
                TFout_sig.append(nested_list[0])
                for t in nested_list[-5].split(','):
                    TFNetwork_df.append((nested_list[0], RE.geneFBnum(nested_list[0]), t.lstrip().replace('"', ''), RE.geneFBnum(t.lower().lstrip().replace('"', '')), '', str(1))) #TF gene name, TF gene FB#, target gene name, target gene FB#, Interaction type, gene source: TF=lightblue, Primary interaction = black, secondary interaction = green
                if (nested_list[MTC]==0):
                    print ("Warning: "+nested_list[0]+"p-value is 0, suggesting something is off. Likely too few permutations.")
    elif (MultiTest == 'PK'):
        TFs = TFanalysis_out[8:]
        TFs.sort(key=itemgetter(5))
        for nested_list in TFs[:20]:
            TFout_sig.append(nested_list[0])
            for t in nested_list[-5].split(','):
                    TFNetwork_df.append((nested_list[0], RE.geneFBnum(nested_list[0]), t.lstrip().replace('"', ''), RE.geneFBnum(t.lower().lstrip().replace('"', '')), '', str(1))) #TF gene name, TF gene FB#, target gene name, target gene FB#, Interaction type, gene source: TF=lightblue, Primary interaction = black, secondary interaction = green
    #Now adding significant TFs from above to the input list for the ribotag analysis
    inputGenes = []
    riboGenes = []
    with open(GENESET, 'r') as f:
        for line in f:
            if line.startswith('Gene'):
                continue
            else:
                if line.lower().startswith('fbgn'):
                    riboGenes.append(RE.FBnumtoComName(line.strip()))#TFanalysis_out.append(line.strip().split('\t'))
                    inputGenes.append(RE.FBnumtoComName(line.strip()))
                else:
                    riboGenes.append(RE.FBnumtoComName(RE.geneFBnum(line.strip()))) #Probably a more efficient way to get the same name for any input gene
                    inputGenes.append(RE.FBnumtoComName(RE.geneFBnum(line.strip()))) #Probably a more efficient way to get the same name for any input gene
    riboGenes.extend(TFout_sig)
    #ribo_out = RE.riboext(riboGenes, GENESETNAME)#, outputpath)
    #ribo_out.pop('Comparison Stats', None)
    #ribo_out = ribo_out.drop(['Comparison Stats'])
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
    LinVals = GA.geneExp_Age(riboGenes, GENESETNAME, save)
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
        next
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
        figsize = len(TFNetwork_df)
    elif len(TFNetwork_df) > 1000:
        figsize = 250
    else:
        figsize = 100
    fig=plt.figure(figsize=(figsize,figsize))
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
    if nwk_type.lower() == 'circle':
        layout = nx.circular_layout(g)
        radii = [30,45,60,75,90,105,120]#[:num_radii]  # for concentric circles

    #Setting nodes based on tissue
    #    for ea in layout.keys():
    #        new_r = 1
    #        if ea in nodes_by_color['gold']:
    #            new_r = radii[0]
    #        elif ea in nodes_by_color['green']:
    #            new_r = radii[1]
    #        elif ea in nodes_by_color['brown']:
    #            new_r = radii[2]
    #        elif ea in nodes_by_color['red']:
    #            new_r = radii[3]
    #        elif ea in nodes_by_color['silver']:
    #            new_r = radii[4]
    #        elif ea in nodes_by_color['blue']:
    #            new_r = radii[5]
    #        elif ea in nodes_by_color['magenta']:
    #            new_r = radii[6]
    #        else:
    #            new_r = radii[7]#pass
    #        layout[ea] *= new_r   # relayoutition nodes as concentric circles
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

    # 3. Draw the parts we want
    nx.draw_networkx_edges(g, layout, width=4, edge_color='#AAAAAA')
    TFsources = [node for node in g.nodes() if node in interdf.sourceFB[interdf.geneSource == "1"].unique()]
    TFsize = [g.degree(node) * 1000 for node in g.nodes() if node in TFsources]
    #Source_color_map = [interdf.geneSource[interdf.sourceFB == node] for node in TFsources]
    TFnodes = nx.draw_networkx_nodes(g, layout, nodelist=TFsources, node_size=TFsize, node_color='darkturquoise', linewidths=5)
    TFnodes.set_edgecolor(pd.Series([largestFoldChange[key] for key in largestFoldChange.keys() if key in TFsources]))
    #nx.draw_networkx_labels(g, pos=nx.kamada_kawai_layout(g), labels=interdf.source.unique())
    TFtargets = [node for node in g.nodes() if node in interdf.targetFB[interdf.geneSource == "1"].unique()]
    Targetsize = [g.degree(node) * 1000 for node in g.nodes() if node in TFtargets]
    #Target_color_map = [interdf.geneSource[interdf.targetFB == node] for node in TFtargets]
    Targetnodes = nx.draw_networkx_nodes(g, layout, nodelist=TFtargets, node_size=Targetsize, node_color='yellow', linewidths=5)#AAAAAA'
    Targetnodes.set_edgecolor(pd.Series([largestFoldChange[key] for key in largestFoldChange.keys() if key in TFtargets]))
    #Highlighting input genes
    Inputgenes = [node for node in g.nodes() if node in inputdf.Input_geneFB.unique()]
    Inputsize = [g.degree(node) * 1000 for node in g.nodes() if node in Inputgenes]
    Inputnodes = nx.draw_networkx_nodes(g, layout, nodelist=Inputgenes, node_size=Inputsize, node_color='yellow', linewidths=5)#AAAAAA'
    Inputnodes.set_edgecolor(pd.Series([largestFoldChange[key] for key in largestFoldChange.keys() if key in Inputgenes]))
    #nx.draw_networkx_labels(g, pos=nx.kamada_kawai_layout(g), labels=interdf.target.unique())
    if nwk_dist == 1:
        Pritargets = [node for node in g.nodes() if node in interdf.targetFB[interdf.geneSource == "2"].unique()]#Pridf.sourceFB]
        Prisize = [g.degree(node) * 500 for node in g.nodes() if node in Pritargets]
        Prinodes = nx.draw_networkx_nodes(g, layout, nodelist=Pritargets, node_size=Prisize, node_color='springgreen', linewidths=5)#, linewidths=5)
        Prinodes.set_edgecolor(pd.Series([largestFoldChange[key] for key in largestFoldChange.keys() if key in Pritargets]))
        #Prinodes.set_edgecolor('blue')
    elif nwk_dist == 2:
        Pritargets = [node for node in g.nodes() if node in interdf.targetFB[interdf.geneSource == "2"].unique()]#Pridf.sourceFB]
        Prisize = [g.degree(node) * 500 for node in g.nodes() if node in Pritargets]
        Prinodes = nx.draw_networkx_nodes(g, layout, nodelist=Pritargets, node_size=Prisize, node_color='springgreen', linewidths=5)#, linewidths=5)
        Prinodes.set_edgecolor(pd.Series([largestFoldChange[key] for key in largestFoldChange.keys() if key in Pritargets]))
        #Prinodes.set_edgecolor('blue')
        Sectargets = [node for node in g.nodes() if node in interdf.targetFB[interdf.geneSource == "3"].unique()]#Secdf.sourceFB]
        Secsize = [g.degree(node) * 250 for node in g.nodes() if node in Sectargets]
        Secnodes = nx.draw_networkx_nodes(g, layout, nodelist=Sectargets, node_size=Secsize, node_color='pink', linewidths=5)
        Secnodes.set_edgecolor(pd.Series([largestFoldChange[key] for key in largestFoldChange.keys() if key in Sectargets]))
        #Secnodes.set_edgecolor('purple')
    enhancers = [(source, target) for source in interdf.sourceFB[interdf.interactionType=='enhanceable'] for target in interdf.targetFB[interdf.interactionType=='enhanceable']]
    enhancers = [e for e in enhancers if e in g.edges()]
    nx.draw_networkx_edges(g, layout, edgelist=enhancers, arrows=True, width=2, edge_color='blue')
    suppressors = [(source, target) for source in interdf.sourceFB[interdf.interactionType=='suppressible'] for target in interdf.targetFB[interdf.interactionType=='suppressible']]
    suppressors = [s for s in suppressors if s in g.edges()]
    nx.draw_networkx_edges(g, layout, edgelist=suppressors, arrows=True, width=2, edge_color='red')
    nx.draw_networkx_labels(g, layout, labels=mapping_labels, font_size=0.4*figsize)
#    high_degree_targets = [node for node in g.nodes() if node in interdf.targetFB.unique() and g.degree(node) > 50]
#    nx.draw_networkx_nodes(g, layout, nodelist=high_degree_targets, node_size=size, node_color='#fc8d62')
    #source_dict = dict(zip(TFsources, TFsources))
    #target_dict = dict(zip(TFtargets, TFtargets))
    #nx.draw_networkx_labels(g, layout, labels=source_dict)
    #nx.draw_networkx_labels(g, layout, labels=target_dict)

    # 4. Turn off the axis because I know you don't want it
    plt.axis('off')

    plt.title("Regulatory Network of Input Genes", fontsize = figsize)

    # 5. Tell matplotlib to show it
    plt.show()
    if save == 1:
        if nwk_dist == 0:
            fig.savefig(GENESETNAME+'_RegulatoryNetwork_Basic.png')
        if nwk_dist == 1:
            fig.savefig(GENESETNAME+'_RegulatoryNetwork_Primary.png')
        if nwk_dist == 2:
            fig.savefig(GENESETNAME+'_RegulatoryNetwork_PriAndSec.png')
#
#    cytoscape=cyrest.cyclient()
#    cy = CyRestClient()
#    cy.session.delete()
#    cytoscape.session.new()
#    nx_network=nx.from_pandas_edgelist(df, 'source', 'target', edge_attr=True)
#    net_module=cy.network.create_from_networkx(nx_network)
#    cy.layout.apply(name='circular', network=net_module)
#    net_module.update_node_table(df=df, network_key_col='source', data_key_col='source')
#    cy.layout.apply(name='circular', network=net_module) #'hierarchical'
#    network = cy.network.create(name='Test', collection='TestNetwork')
#    for source in testnetwork['source']:
#        network.add_node(source)
#
#    cyrest.network.create_attribute(self, column=['Test Attribute'], listType=['string'], network=['testnetwork']) #Attempting to add an attribute, get error on self
#    df_network=network.get_node_table()
#    cytoscape.vizmap.apply(styles="default")
    return (TFanalysis_out, ribo_out)
#    geneSetf, TFChipSeqf, bckgf, outf, logf, SETUP = GSE.get_setup(outputpath)
#
#    bckgGenes, lineNum, duplicates = GSE.read_gene_set(bckgf, SETUP['BACKGROUND'], SETUP['BACKGROUNDNAME'], logf)
#    origGeneSet, lineNum, duplicates = GSE.read_gene_set(geneSetf, SETUP['GENESET'], SETUP['GENESETNAME'], logf)
#    geneSet = GSE.remove_genes_not_in_bckgrnd(origGeneSet, bckgGenes.keys(), logf)
#    GSE.print_and_log("%d unique genes in gene-set also in background list" % (len(geneSet)), logf)
#    if SETUP['EXECMODE'] == 'DEB':
#        logf.write("Unique genes in gene-set and in background list: ")
#        for gene in geneSet:
#            logf.write(" %s")
#        logf.write("\n")
#
#    setupFileName = 'GSE.setup'
#    GSE.write_output_header(outf, SETUP, len(origGeneSet), len(geneSet))
#
#    line = TFChipSeqf.readline() # header
#    line = TFChipSeqf.readline()
#    lineNum = 1
#    while line:
#         TF, peaks, line, lineNum = GSE.read_TF_binding_sites(TFChipSeqf, line, lineNum)
#         GSE.print_and_log("%s:\t%d genes" % (TF, len(peaks.keys())), logf)
#         GSE.enrichment_analysis(TF, peaks, bckgGenes, geneSet, bckgGenes.keys(), outf)
#
#    GSE.TFmulti(outf, outputpath, TF_FDR)