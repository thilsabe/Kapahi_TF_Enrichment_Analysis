# Transcription Factor Enrichment Analysis Pipeline

A Python pipeline for identifying enriched transcription factors (TFs) 
upstream of a gene set of interest, using 445 publicly available ChIP-Seq 
datasets from *Drosophila melanogaster* stage 3 larvae.

Developed in the **Kapahi and Brem Labs at the Buck Institute for Research 
on Aging**. Original enrichment code written by Tal Oron (Buck Institute); 
modifications, extensions, and pipeline integration by Tyler A.U. Hilsabeck.

---

## Overview

Given an input gene list, this pipeline:

1. Maps gene names/symbols to FlyBase IDs
2. Tests for enrichment of those genes among the known targets of 445 TFs 
   from publicly available *Drosophila* ChIP-Seq data
3. Applies multiple testing correction to identify statistically significant 
   TF-target enrichments
4. Performs gene ontology (GO) analysis on the input gene set
5. Visualizes significant TFs and their targets as an interactive network

This approach is useful for identifying upstream transcriptional regulators 
of a gene set of interest — for example, differentially expressed genes from 
an RNA-seq experiment or hits from a genetic screen — and placing them in 
the context of known TF-target relationships and genetic interactions in 
*Drosophila*.

---

## Repository Contents

### Core Pipeline

| File | Description |
|---|---|
| `genePipeline.py` | Main entry point. Reads an input gene list and coordinates downstream analysis |
| `GSE_benhoch.py` | TF enrichment analysis with multiple testing correction (modified from original Tal Oron code) |
| `GSEutils.py` | Utility/support functions called by `GSE_benhoch.py` |
| `RibotagExtract.py` | Support methods called by `genePipeline.py` |
| `geneOnto.py` | Gene ontology (GO) analysis of input gene sets |
| `drawNetwork_Nshow.py` | Visualizes significant TFs and their targets in 3 network orientations |
| `geneExp_withAge_Diet.py` | Analysis of gene expression patterns across age and dietary conditions |

### Reference Data

| File | Description | Source |
|---|---|---|
| `All_Background_FB.txt` | Background gene set for enrichment testing. Can be subsetted to tissue-specific gene sets for context-specific backgrounds | FlyBase |
| `fb_annotation_synonym_fb_2019_03_dmel.txt` | *Drosophila* gene names and FlyBase IDs | FlyBase (2019) |
| `gene_genetic_interactions_fb_2019_04.txt` | Known genetic interactions in *Drosophila* | FlyBase (2019) |
| `gene_group_data_fb_2020_02.tsv` | Gene ontology terms | FlyBase (2020) |

---

## Usage
```python
# Run the full pipeline with an input gene list
python genePipeline.py --input your_gene_list.txt

# Draw network of significant TFs and targets
python drawNetwork_Nshow.py
```

**Input format:** A plain text file with one gene name or FlyBase ID per line.

**Background:** By default, `All_Background_FB.txt` is used as the background 
gene set. For tissue-specific analyses, this can be subsetted to genes 
expressed in the tissue of interest.

---

## Methods

Enrichment of input genes among TF targets is tested using a hypergeometric 
test across 445 *Drosophila melanogaster* ChIP-Seq datasets. Multiple testing 
correction is applied via the Benjamini-Hochberg procedure. Significant TF 
enrichments are visualized as force-directed network graphs showing TFs, 
their targets within the input gene list, and known genetic interactions.

---

## Requirements

- Python 3.x
- Standard scientific Python stack: `numpy`, `scipy`, `pandas`, `matplotlib`
- NetworkX (for network visualization)

---

## Biological Context

This pipeline was developed to support upstream regulator discovery in the 
context of aging and dietary restriction research in *Drosophila*. By 
identifying which transcription factors are significantly enriched among 
gene sets of interest, it provides mechanistic hypotheses about the 
transcriptional programs underlying diet- and age-dependent phenotypes.

---

## Acknowledgments

Original enrichment code written by **Tal Oron** at the Buck Institute for 
Research on Aging. Extended and integrated into the current pipeline by 
**Tyler A.U. Hilsabeck** in the Kapahi and Brem Labs at the Buck Institute.

ChIP-Seq reference data sourced from publicly available *Drosophila* 
datasets via FlyBase (flybase.org).
