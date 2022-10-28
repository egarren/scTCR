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
4. Run each script in numerical order in either Rstudio or Python.  These scripts will generate the figures presented in our manuscript.\
\
NB: Expected total run time is 3-5 days

## System requirements and software

<ins>Python (v3.7.4) packages</ins>\
[AlphaFold](https://github.com/deepmind/alphafold) v2.1.1\
[GIANA](https://github.com/s175573/GIANA) v4.1\
[tcrdist3](https://github.com/kmayerb/tcrdist3) v0.2.2\
[DeepTCR](https://github.com/sidhomj/DeepTCR) v2.1.0


[GLIPH](http://50.255.35.37:8080/) v2\

<ins>R (v4.1.1) packages</ins>\
ggalluvial_0.12.3\
harmony_0.1.0\
Rcpp_1.0.9\
sp_1.5-0               
SeuratObject_4.1.2   
Seurat_4.2.0          
ggseqlogo_0.1    
ggforce_0.4.1     
phylotools_0.2.4   
ape_5.6-2          
viridis_0.6.2      
viridisLite_0.4.1   
VennDiagram_1.7.3    
futile.logger_1.4.3  
limma_3.46.0        
EnhancedVolcano_1.13.2\
plotly_4.10.0       
biomaRt_2.46.3        
plyr_1.8.7         
forcats_0.5.2       
stringr_1.4.1         
purrr_0.3.4       
readr_2.1.3        
tidyverse_1.3.2       
qgraph_1.9.2       
cowplot_1.1.1      
ggrepel_0.9.1         
ggpubr_0.4.0.999    
immunarch_0.8.0     
patchwork_1.1.2       
dtplyr_1.2.2       
dplyr_1.0.10       
data.table_1.14.2     
vegan_2.6-4        
lattice_0.20-45     
permute_0.9-7         
xlsx_0.6.5          
RColorBrewer_1.1-3   
scales_1.2.1          
gridExtra_2.3       
ggplot2_3.3.6       
dendextend_1.16.0     
corrplot_0.92       
pheatmap_1.0.12      
tidyr_1.2.1           
tibble_3.1.8     



