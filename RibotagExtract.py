import pandas as pd
from scipy import stats
#from scipy.stats import shapiro
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import csv
import io
import sys
#%matplotlib inline

text_trap = io.StringIO()
sys.stdout = text_trap
#To restore output: sys.stdout = sys.__stdout__
FBFile = 'fb_annotation_synonym_fb_2019_03_dmel.txt' #Combined using fb_annotation and fb_synonyms data for d6 on flybase
FBcodeList=pd.read_csv(FBFile,sep="\t")
for i in range(0, len(FBcodeList['current_symbol'])):
    try:
        FBcodeList['current_symbol'][i] = FBcodeList['current_symbol'][i].lower()
    except:
        continue

    #mpl.use('agg')
#Extraction code from online used to extract tissue information from each input gene
def gen_dict_extract(var, key):
    if isinstance(var, dict):
        for k, v in var.items():
            if k == key:
                yield v
            if isinstance(v, (dict, list)):
                yield from gen_dict_extract(v, key)
    elif isinstance(var, list):
        for d in var:
            yield from gen_dict_extract(d, key)

def mergedict(a,b):
    a.update(b)
    return a
    ###Extract information from ribotag dataset

def geneFBnum(gene):
    FBnum=[] #defined early to catch errors for longer lists
    if gene.lower().startswith('fbgn'):
        FBnum = gene.strip()
    else:
        try:
            FBnum = FBcodeList['##primary_FBid'][list(FBcodeList['current_symbol']).index(gene.lower().strip())]
        except:
            try:
                FBnum = FBcodeList['##primary_FBid'][list(FBcodeList['annotation_ID']).index(gene.strip())]
            except:
                try:
                    found = 0
                    place = 0
                    for s in list(FBcodeList['symbol_synonyms']):
                        if found == 1:
                            break
                        try:
                            if np.isnan(s):
                                place += 1
                                continue
                        except Exception:
                            pass
                        for g in s.split(sep=','):
                            if gene.lower() == g.lower():
                                FBnum = FBcodeList['##primary_FBid'][place]
                                print (gene+' FB# was found at position '+str(place)+'.')
                                found = 1
                        else:
                            place += 1
                except:
                    FBnum = 'NotFound'
                    print (str(gene)+" FB# not found. Need to update FBid in script OR name likely not on flybase or is wrong.")
    if FBnum == []:
        FBnum = 'NotFound'
    return (FBnum)

def FBnumtoComName(FB):
    try:
        comName = FBcodeList['current_symbol'][list(FBcodeList['##primary_FBid']).index(FB.lower())]
    except:
        comName = 'NotFound'
        print (str(FB)+" not found. Need to update FBid in database OR is wrong.")
    return (comName)

def riboext(geneList, geneListname, save=1):#, outputpath):
    #Provide gene list file name for a tab separated text file with column heading Gene
    #GeneFileName='cnstMtdGenes.csv'#BetaQuantCase_TFanalysis.csv'#D21to95_betaCaseQuantGenes.csv'#D21to95_ALDRdeltabetaCaseQuantGenes.csv'#D21to95_betaCaseGenes.csv'#PurineGenes.csv' #Variance_Genes.csv' #MLE_betaGenes.csv'#'D21to95_AL_alphaGenes.csv'#'
    #save = 0 (no) or 1 (yes)

    #Loading requiresite files and initiating output dictionary
    geneList = pd.DataFrame(geneList, columns=['Gene'])#geneList=pd.read_csv(GeneFileName)
    neuron=pd.read_csv('Neuron_Results_all.csv')
    fatbody=pd.read_csv('Fat body_Results_all.csv')
    germ=pd.read_csv('Germ line_Results_all.csv')
    gut=pd.read_csv('Gut Results_all.csv')
    heart=pd.read_csv('Heart_Results_all.csv')
    muscle=pd.read_csv('Muscle_Results_all.csv')
    tubule=pd.read_csv('Tubule_Results_all.csv')
    geneOutput = dict()
    tissues = ['Neuron', 'Gut', 'Germ', 'Fatbody', 'Heart', 'Muscle', 'Tubule']

    #Running through each gene from input list, provided as gene symbol or FBgn number, and getting
    for gene in geneList['Gene']:
        geneOutput[str(gene.strip())] = {}
        geneOutput[str(gene.strip())]['FB#'] = geneFBnum(gene)
#        if gene.startswith('FBgn'):
#            geneOutput[str(gene.strip())]['FB#'] = gene.strip()
#        else:
#            try:
#                geneOutput[str(gene.strip())]['FB#'] = FBcodeList['##primary_FBid'][list(FBcodeList['current_symbol']).index(gene.strip())]
#                #print ('tried 1')
#            except:
#                try:
#                    geneOutput[str(gene.strip())]['FB#'] = FBcodeList['##primary_FBid'][list(FBcodeList['annotation_ID']).index(gene.strip())]
#                    #print ('tried 2')
#                except:
#                    try:
#                        #print ('tried 3')
#                        #geneOutput[str(gene)]['FB#'] = FBcodeList['##primary_FBid'][list(FBcodeList['symbol_synonyms']).index(gene)]
#                        found = 0
#                        place = 0
#                        for s in list(FBcodeList['symbol_synonyms']):
#                            if found == 1:
#                                break
#                            #print ("Position: "+str(place))
#                            if gene in s.split(sep=','):
#                                geneOutput[str(gene)]['FB#'] = FBcodeList['##primary_FBid'][place]
#                                print (gene+' FB# was found at position '+str(place)+'.')
#                                found = 1
#                            else:
#                                place += 1
#                    except:
#                        print (str(gene)+" FB# not found. Need to update FBid in script OR name likely not on flybase or wrong.")
        try:
            geneOutput[str(gene)]['Neuron Pval'] = neuron['P value'][list(neuron['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Neuron FoldChange'] = neuron['Fold change'][list(neuron['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Neuron AL Expression'] = neuron['AL'][list(neuron['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Neuron DR Expression'] = neuron['DR'][list(neuron['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
        except:
            print (str(gene)+" not found in neuron list.")
        try:
            geneOutput[str(gene)]['Gut Pval'] = gut['P value'][list(gut['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Gut FoldChange'] = gut['Fold change'][list(gut['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Gut AL Expression'] = gut['Myo1A AL'][list(gut['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Gut DR Expression'] = gut['Myo1A DR'][list(gut['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
        except:
            print (str(gene)+" not found in gut list.")
        try:
            geneOutput[str(gene)]['Germ Pval'] = germ['P value'][list(germ['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Germ FoldChange'] = germ['Fold change'][list(germ['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Germ AL Expression'] = germ['NGT AL'][list(germ['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Germ DR Expression'] = germ['NGT DR'][list(germ['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
        except:
            print (str(gene)+" not found in germ list.")
        try:
            geneOutput[str(gene)]['Fatbody Pval'] = fatbody['P value'][list(fatbody['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Fatbody FoldChange'] = fatbody['Fold change'][list(fatbody['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Fatbody AL Expression'] = fatbody['S106 AL'][list(fatbody['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Fatbody DR Expression'] = fatbody['S106 DR'][list(fatbody['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
        except:
            print (str(gene)+" not found in fatbody list.")
        try:
            geneOutput[str(gene)]['Heart Pval'] = heart['P value'][list(heart['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Heart FoldChange'] = heart['Fold change'][list(heart['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Heart AL Expression'] = heart['Hand4.2 AL'][list(heart['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Heart DR Expression'] = heart['Hand4.2 DR'][list(heart['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
        except:
            print (str(gene)+" not found in heart list.")
        try:
            geneOutput[str(gene)]['Muscle Pval'] = muscle['P value'][list(muscle['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Muscle FoldChange'] = muscle['Fold change'][list(muscle['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Muscle AL Expression'] = muscle['MHC AL'][list(muscle['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Muscle DR Expression'] = muscle['MHC DR'][list(muscle['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
        except:
            print (str(gene)+" not found in muscle list.")
        try:
            geneOutput[str(gene)]['Tubule Pval'] = tubule['P value'][list(tubule['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Tubule FoldChange'] = tubule['Fold change'][list(tubule['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Tubule AL Expression'] = tubule['AL C42'][list(tubule['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
            geneOutput[str(gene)]['Tubule DR Expression'] = tubule['DR C42'][list(tubule['GENE_ID']).index(geneOutput[str(gene)]['FB#'])]#FBcodeList['primary_FBgn#'][list(FBcodeList['##gene_symbol']).index(gene)])]
        except:
            print (str(gene)+" not found in tubule list.")

    if save == 1:
        #Running tissue comparison stats
        geneOutput['Comparison Stats'] = {}
        for  tissue in tissues:#['Neuron', 'Gut', 'Fatbody', 'Heart', 'Muscle', 'Tubule']:
            for compt in tissues:#['Neuron', 'Gut', 'Fatbody', 'Heart', 'Muscle', 'Tubule']:
                if (tissue == compt):
                    next
                elif (tissues.index(compt)<=tissues.index(tissue)):
                    next
                else:
                    try:
                        if (stats.shapiro(list(gen_dict_extract(geneOutput, tissue+' FoldChange')))[1] > 0.05) & (stats.shapiro(list(gen_dict_extract(geneOutput, compt+' FoldChange')))[1] > 0.05): #If both are normal, run parametric t-test, else run non-parametric mannwhitneyu
                            #First value is test statistic, second is p-value. Stored as a tuple under Parametric Comparisons dictionary
                            geneOutput['Comparison Stats'][tissue+' and '+compt] = tuple(stats.ttest_ind(list(gen_dict_extract(geneOutput, tissue+' FoldChange')),list(gen_dict_extract(geneOutput, compt+' FoldChange'))))
                        else:
                            #First value is test statistic, second is p-value. Stored as a tuple under Nonparametric Comparisons dictionary
                            geneOutput['Comparison Stats'][tissue+' and '+compt] = tuple(stats.mannwhitneyu(list(gen_dict_extract(geneOutput, tissue+' FoldChange')),list(gen_dict_extract(geneOutput, compt+' FoldChange'))))
                    except:
                        print ("Error during stats, possibly due to lack of significant TFs.")
        #Paired nonparametric ANOVA: h, p2 = stats.kruskal(list(gen_dict_extract(geneOutput, 'Neuron FoldChange')),list(gen_dict_extract(geneOutput, 'Gut FoldChange')),list(gen_dict_extract(geneOutput, 'Germ FoldChange')),list(gen_dict_extract(geneOutput, 'Fatbody FoldChange')),list(gen_dict_extract(geneOutput, 'Heart FoldChange')),list(gen_dict_extract(geneOutput, 'Muscle FoldChange')),list(gen_dict_extract(geneOutput, 'Tubule FoldChange')))
        #Unpaired nonparametric ANOVA: t, p3 = stats.ranksums(list(gen_dict_extract(geneOutput, 'Neuron FoldChange')),list(gen_dict_extract(geneOutput, 'Gut FoldChange')))
        #Using the log-transformed values makes for a stronge interpretation of the data. Sticking with nontransformed, since a fold change of 1 would be no change, and anything else would be the proportional change. Logged code follows: h3, p3 = stats.f_oneway(np.log(list(gen_dict_extract(geneOutput, 'Neuron FoldChange'))),np.log(list(gen_dict_extract(geneOutput, 'Gut FoldChange'))),np.log(list(gen_dict_extract(geneOutput, 'Germ FoldChange'))),np.log(list(gen_dict_extract(geneOutput, 'Fatbody FoldChange'))),np.log(list(gen_dict_extract(geneOutput, 'Heart FoldChange'))),np.log(list(gen_dict_extract(geneOutput, 'Muscle FoldChange'))),np.log(list(gen_dict_extract(geneOutput, 'Tubule FoldChange'))))
        #sns.distplot(np.log(list(gen_dict_extract(geneOutput, 'Neuron FoldChange'))))#, kde = False)
        #g.set(xlabel="log10(P-Value)", ylabel="Density", title = "Histogram of Ranged Non-Target Gene P-values - RLM Linear Regression")

        #Creating boxplot and plotting significance bars  based on above stats for tissues with significantly different expression patters in gene set:
        fig = plt.figure(figsize=(9, 6))
        ax = fig.add_subplot(111)
        ax.boxplot([list(gen_dict_extract(geneOutput, 'Neuron FoldChange')),list(gen_dict_extract(geneOutput, 'Gut FoldChange')),list(gen_dict_extract(geneOutput, 'Germ FoldChange')),list(gen_dict_extract(geneOutput, 'Fatbody FoldChange')),list(gen_dict_extract(geneOutput, 'Heart FoldChange')),list(gen_dict_extract(geneOutput, 'Muscle FoldChange')),list(gen_dict_extract(geneOutput, 'Tubule FoldChange'))])

        #Extracting p-values for all comparisons and putting stars on boxplot comparisons with pvalues <=0.05
        comppvals = []
        #Flattening list of fold change lists then determining max in order to set height of significance bars (y); h is the bar spacer distance, 'k' is color black
        sigN = -1 #offsetting location on boxplot. Was originally placing line one box ahead
        try:
            y, h, col = max([item for sublist in [list(gen_dict_extract(geneOutput, 'Neuron FoldChange')),list(gen_dict_extract(geneOutput, 'Gut FoldChange')),list(gen_dict_extract(geneOutput, 'Germ FoldChange')),list(gen_dict_extract(geneOutput, 'Fatbody FoldChange')),list(gen_dict_extract(geneOutput, 'Heart FoldChange')),list(gen_dict_extract(geneOutput, 'Muscle FoldChange')),list(gen_dict_extract(geneOutput, 'Tubule FoldChange'))] for item in sublist]) + 1, 0.5, 'k'
        except:
            print ('Error in figure setup')
        for c in geneOutput['Comparison Stats']:
            comppvals.append(geneOutput['Comparison Stats'][c][1])
            if geneOutput['Comparison Stats'][c][1] <= 0.05:
                #print (comppvals[-1])
                sigN += 1
                x1 = tissues.index(c.split(' ', 2)[0])+1 #Setting x1 equal to first tissue in comparison, which corresponds to position on boxplot
                x2 = tissues.index(c.split(' ', 2)[2])+1 #Setting x2 equal to second tissue in comparison
                plt.plot([x1, x1, x2, x2], [y+sigN, y+h+sigN, y+h+sigN, y+sigN], lw=1.5, c=col)
                if comppvals[-1] <= 0.0005:
                    plt.text((x1+x2)*.5, y+h+sigN, "***", ha='center', va='bottom', color=col)
                elif comppvals[-1] <= 0.005:
                    plt.text((x1+x2)*.5, y+h+sigN, "**", ha='center', va='bottom', color=col)
                else:
                    plt.text((x1+x2)*.5, y+h+sigN, "*", ha='center', va='bottom', color=col)

        ax.set_xticklabels(['Neuron', 'Gut', 'Germ', 'Fatbody', 'Heart', 'Muscle', 'Tubule'])
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.set(ylabel="DR/AL Fold Change", title = "Ribotag Tissue-specific fold changes of "+geneListname)
        fig.savefig(geneListname+'_TissueFoldChanges_Boxplot.png', bbox_inches='tight')

    #Removing any genes that were not found or only in a single tissue. This is because of errors when trying to create column indexes for heatmaps.
    geneOutput_export = geneOutput
    if save == 1:
        for key in list(geneOutput.keys()):
            if len(geneOutput[key]) <= 1:
                del geneOutput[key]
        #Specifying Index and Column titles for heatmaps
        #FCIndex = ['Neuron Fold Change', 'Gut Fold Change', 'Germ Fold Change', 'Fatbody Fold Change', 'Heart Fold Change', 'Muscle Fold Change', 'Tubule Fold Change']
        #ALDRIndex = ['Neuron AL Expression', 'Neuron DR Expression', 'Gut AL Expression', 'Gut DR Expression', 'Germ AL Expression', 'Germ DR Expression', 'Fatbody AL Expression', 'Fatbody DR Expression', 'Heart AL Expression', 'Heart DR Expression', 'Muscle AL Expression', 'Muscle DR Expression', 'Tubule AL Expression', 'Tubule DR Expression']

        #Creating dataframe for heatmap creation, setting Y and X axis names to indexes and column names defined above. Dropping comparison stats in figures.
        #Need to figure out how to sort by row with maximum sum so genes are grouped
        geneOutput_df = pd.DataFrame.from_dict(geneOutput, orient = 'index').drop(['Comparison Stats']).fillna(0)
        geneOutput_df.sum(0).drop(['FB#']).max()
        Cols = [k for k, v in geneOutput.items()][:-1] #Used to calculate figure dimensions below

        #Now creating clustermaps and boxplot figures
        fig = sns.clustermap(geneOutput_df.filter(items=['Neuron FoldChange', 'Gut FoldChange', 'Germ FoldChange', 'Fatbody FoldChange', 'Heart FoldChange', 'Muscle FoldChange', 'Tubule FoldChange']).transpose(), figsize=(len(Cols)/2,5), annot = True)
        fig.savefig(geneListname+'_TissueFoldChanges_Clustermap.png')#, bbox_inches='tight') #Saving figure with name which includes input file name

        fig = sns.clustermap(geneOutput_df.filter(items=['Neuron AL Expression', 'Neuron DR Expression', 'Gut AL Expression', 'Gut DR Expression', 'Germ AL Expression', 'Germ DR Expression', 'Fatbody AL Expression', 'Fatbody DR Expression', 'Heart AL Expression', 'Heart DR Expression', 'Muscle AL Expression', 'Muscle DR Expression', 'Tubule AL Expression', 'Tubule DR Expression']).transpose(), figsize=(len(Cols)/2,5), annot = True)
        fig.savefig(geneListname+'_TissueExpression_ALDR_Clustermap.png')#, bbox_inches='tight') #Saving figure with name which includes input file name

    if save == 1:
        #Nested dictionary writing script from stackoverflow: https://stackoverflow.com/questions/29400631/python-writing-nested-dictionary-to-csv
        fields = ['Gene', 'FB#', 'Fatbody AL Expression', 'Fatbody DR Expression', 'Fatbody FoldChange', 'Fatbody Pval', 'Germ AL Expression', 'Germ DR Expression', 'Germ FoldChange', 'Germ Pval', 'Gut AL Expression', 'Gut DR Expression', 'Gut FoldChange', 'Gut Pval', 'Heart AL Expression', 'Heart DR Expression', 'Heart FoldChange', 'Heart Pval', 'Muscle AL Expression', 'Muscle DR Expression', 'Muscle FoldChange', 'Muscle Pval', 'Neuron AL Expression', 'Neuron DR Expression', 'Neuron FoldChange', 'Neuron Pval', 'Tubule AL Expression', 'Tubule DR Expression', 'Tubule FoldChange', 'Tubule Pval', 'Fatbody and Heart', 'Neuron and Heart', 'Germ and Heart', 'Heart and Muscle', 'Neuron and Tubule', 'Germ and Tubule', 'Gut and Muscle', 'Gut and Heart', 'Gut and Fatbody', 'Muscle and Tubule', 'Heart and Tubule', 'Fatbody and Tubule', 'Neuron and Fatbody', 'Neuron and Muscle', 'Fatbody and Muscle', 'Germ and Muscle', 'Gut and Germ', 'Germ and Fatbody', 'Gut and Tubule', 'Neuron and Germ', 'Neuron and Gut']

        with open(geneListname+'_RibotagTable.csv', "w", newline='') as f:
            w = csv.DictWriter(f, fields)
            w.writeheader()
            for k,d in sorted(geneOutput_export.items()):
                w.writerow(mergedict({'Gene': k},d))
    return geneOutput_export