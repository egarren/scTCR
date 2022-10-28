rm(list=ls())
library(Seurat)
library(cowplot)
library(ggplot2)
library(plyr)
library(dplyr)
library(biomaRt)
library(plotly)
library(scales)
library(EnhancedVolcano)
library(data.table)
library(ggpubr)
library(limma)
library(VennDiagram)
library(viridis)
library(pheatmap)
library(phylotools) 
library(ggforce)
library(harmony)



#load data
metadata<-read.csv("Tcell.metadata.csv", header=T) #define path
metadata[] <- lapply(metadata, as.character)
mice<-metadata$T.sampleID
for (i in mice){
  path<-paste0("/data/cellranger/",i,"/outs/per_sample_outs/",i,"/count/sample_filtered_feature_bc_matrix")
  seurat.object<-CreateSeuratObject(counts = Read10X(data.dir = path,gene.column=1), min.cells = 3, min.features  = 200, project = i, assay = "RNA")
  assign(i,seurat.object)
}


##Add Count Matrices (eg. all 564 samples, all AID samples)
seurat.list<-lapply(mice,get)
T_all<-merge(seurat.list[[1]], y=seurat.list[2:length(seurat.list)],add.cell.ids=mice,project="564vsAID")
# add sample metadata
T_all<-AddMetaData(T_all, metadata=rownames(T_all@meta.data),col.name = "cellID")
T_all@meta.data<-left_join(x = T_all@meta.data, y = metadata, by = c("orig.ident"="T.sampleID"),keep=F)
rownames(T_all@meta.data)<-T_all@meta.data$cellID
rm(list=setdiff(ls(), c("T_all", "metadata")))


##QC and Filter
#adding feature (gene) metadata
T_all@assays[["RNA"]]@meta.features$original_ensembl<-rownames(T_all@assays[["RNA"]]@meta.features)
genes.meta<-getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "start_position", "end_position", "chromosome_name", 
                      "percentage_gene_gc_content", "external_gene_name", "gene_biotype","go_id","name_1006"),filters=
                    "ensembl_gene_id",values=list(rownames(T_all@assays[["RNA"]]@meta.features)),
                  mart=useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="ensembl.org"),useCache=F) #useast.
m <- match(T_all@assays[["RNA"]]@meta.features$original_ensembl, genes.meta$ensembl_gene_id)
T_all@assays[["RNA"]]@meta.features<-cbind(T_all@assays[["RNA"]]@meta.features,genes.meta[m,])
#renaming features to gene symbols
rownames(T_all@assays[["RNA"]]@meta.features) <- make.names(T_all@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
rownames(T_all@assays[["RNA"]]@data) <- make.names(T_all@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
rownames(T_all@assays[["RNA"]]@counts) <- make.names(T_all@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
# #TCR and Ig identification and filtering
ig_list <- c("IG_C_gene", "IG_C_pseudogene", "IG_D_gene", "IG_D_pseudogene", "IG_J_gene", "IG_LV_gene", 
             "IG_pseudogene", "IG_V_gene", "IG_V_pseudogene")
tr_list <-c("TR_V_gene", "TR_V_pseudogene", "TR_D_gene", "TR_J_gene", "TR_J_pseudogene", "TR_C_gene")
T_all <- subset(T_all, features=rownames(T_all[!(T_all@assays[["RNA"]]@meta.features$gene_biotype %in% ig_list),]))
T_all <- subset(T_all, features=rownames(T_all[!(T_all@assays[["RNA"]]@meta.features$gene_biotype %in% tr_list),]))
T_all <- subset(T_all, features=rownames(T_all[!(is.na(T_all@assays[["RNA"]]@meta.features$mgi_symbol)),]))
#eliminate 18S rRNA contaminating genes (not necessary with new ensembl)
rRNA.genes<-c("AY036118","Gm42418")
T_all <- subset(T_all, features=rownames(T_all[!(T_all@assays[["RNA"]]@meta.features$mgi_symbol %in% rRNA.genes),]))
T_all <- subset(T_all, features=rownames(T_all[!(rownames(T_all) %in% rRNA.genes),]))
# mitochondrial DNA percentage
T_all <- PercentageFeatureSet(T_all, pattern = "^mt", col.name = "percent.mt")
T_all <- PercentageFeatureSet(T_all, pattern = "^Rpl", col.name = "percent.Rpl")
T_all <- PercentageFeatureSet(T_all, pattern = "^Rps", col.name = "percent.Rps")
#cell cycle regression
m_cc<-readRDS("m_cc.rds") #define path
T_all <- CellCycleScoring(T_all, s.features = m_cc$s.genes, g2m.features = m_cc$g2m.genes, set.ident = TRUE)
T_all$CC.Difference <- T_all$S.Score - T_all$G2M.Score
#HSP regression
HSP_genes<-genes.meta[genes.meta$go_id=="GO:0034605",]$mgi_symbol
T_all <- AddModuleScore(object = T_all,features = list(HSP_genes),name = 'HSP.score')
Idents(T_all)<-"condition"
VlnPlot(object = T_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                                     "percent.Rpl","percent.Rps","HSP.score1"),pt.size=0, ncol = 3)
ggsave2("vln.QC.png",device="png")
dim(T_all)
T_all<- subset(x = T_all, subset = nFeature_RNA > 800 & nFeature_RNA < 3500 & nCount_RNA>2000 & 
                  percent.mt >  -Inf & percent.mt < 7 & percent.Rpl < 20 & percent.Rps < 15 &
                 S.Score <0.15 & G2M.Score<0.15) 
dim(T_all)
rm(list=setdiff(ls(), c("T_all", "metadata","genes.meta")))


#Harmony integration
T_all <- T_all %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress=c("nCount_RNA", "percent.mt","HSP.score1","percent.Rpl","percent.Rps"),verbose = FALSE) %>% 
  RunPCA(pc.genes = T_all@var.genes, npcs = 20, verbose = FALSE)
T_all.combined<-RunHarmony(T_all,c("condition","batch","mouse_ID"),plot_convergence=T)
T_all.combined <- RunUMAP(T_all.combined, reduction = "harmony", dims = 1:20) #change based on elbow plot
T_all.combined <- FindNeighbors(T_all.combined, reduction = "harmony", dims = 1:20) #change based on elbow plot
T_all.combined <- FindClusters(T_all.combined, resolution = 0.3) #adjust resolution (bigger=more clusters), initially used 0.2
DimPlot(T_all.combined, reduction = "umap")
for (i in c("condition","Phase","mouse_ID","gender","batch","seurat_clusters")){
  Idents(T_all.combined)<-i
  DimPlot(T_all.combined, reduction = "umap",pt.size=0.1)
  ggsave2(paste0(i,".umap.png"),width=6, height=5,device="png")
}
for (i in c("nCount_RNA","nFeature_RNA","percent.mt","percent.Rps","percent.Rpl","CC.Difference","HSP.score1","S.Score","G2M.Score")){
  FeaturePlot(T_all.combined, features= i,split.by = "condition",pt.size=0.1, order=T)
  ggsave2(paste0(i,".umap.png"),width=10, height=5,device="png")
}
rm(list=setdiff(ls(), c("T_all.combined", "metadata","genes.meta")))
save.image("merged.RData")

