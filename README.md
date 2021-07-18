# Kapahi_TF_Enrichment_Analysis
Code determines the enrichment of input genes among targets from 445 publicly available ChIP-Seq data sets of Drosophila melanogaster stage 3 larvae. Original enrichment code written by Tal Oron at the Buck Institute, modifications and additions made by Tyler Hilsabeck in the Kapahi and Brem Labs at the Buck Institute.
Code Descriptions:
GSE_benhoch.py - Modified enrichment code that includes multiple testing correction
GSEutils.py - Support methods called in GSE_benhoch.py
genePipeline.py - Pipeline that integrates reading of an input gene list and its analysis
RibotagExtract.py - Support methods called in genePipeline.py
drawNetwork_Nshow.py - Displays significant TFs and their targets in 3 different orientations
geneOnto.py - Gene ontology analysis
All_Background_FB.txt - Background set of Drosophila genes. Can be paired down to tissue-specific gene sets and used for tissue-specific backgrounds
fb_annotation_synonym_fb_2019_03_dmel.txt - File containing Drosophila gene names and their respective FB numbers, from flybase.org
gene_genetic_interactions_fb_2019_04.txt - File containing Drosophila gene known interactors, from flybase.org
gene_group_data_fb_2020_02.tsv - File containing gene ontology terms
