# scTCR

This repository contains the code used in our single cell sequencing paper: "Loss of B cell Tolerance is TCR Dependent"

## Installation guide
Install dependencies listed below

## Demo
Model data is included in the `data` directory.

## Instructions
1. Download `data` and `code` directories
2. Set `data` as the working directory
3. Download `gex` and `vdj` data from GSE157649, available here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157649
4. Download `mcr` data from GSE157649, available here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157649
5. Run each script in numerical order.  These scripts will generate the figures presented in our manuscript.\
\
NB: Expected total run time is 3-5 days

## System requirements and software


[GLIPH](http://50.255.35.37:8080/) v2\
[cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi) v7.0.0

<ins>Python (v3.7.4) packages</ins>\
[AlphaFold](https://github.com/deepmind/alphafold) v2.1.1\
[GIANA](https://github.com/s175573/GIANA) v4.1\
[tcrdist3](https://github.com/kmayerb/tcrdist3) v0.2.2\
[DeepTCR](https://github.com/sidhomj/DeepTCR) v2.1.0

<ins>R (v4.1.1) packages</ins>\
harmony_0.1.0\
Seurat_4.2.0        
clusterProfiler_4.6.0\
SPIA_2.50.0\
qgraph_1.9.2 \
ggfortify_0.4.15 \
ggpubr_0.4.0.999 \
ggplot2_3.3.6 \
pheatmap_1.0.12 


