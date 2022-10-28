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
my.ttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

#
load("merged.RData")

## Cluster Analysis
T_all.markers <- FindAllMarkers(object = T_all.combined, test.use = "MAST")
write.csv(T_all.markers,"cluster.markers.csv")
top6 <- T_all.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
top4 <- T_all.markers %>% group_by(cluster) %>% top_n(4, avg_log2FC)
top2 <- T_all.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
top1 <- T_all.markers %>% group_by(cluster) %>% top_n(1, avg_log2FC)
DoHeatmap(object = T_all.combined, features = top6$gene, label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2("cluster.heatmap.png",width=4, height=8,device="png")
T_all.combined@meta.data$my.clusters <- Idents(T_all.combined)  # Store cluster identities in object@meta.data$my.clusters
cluster.counts<-table(T_all.combined@meta.data$my.clusters, T_all.combined@meta.data$condition)
clusters<-names(table(T_all.combined$my.clusters))

#set metadata and subsets
mice<-names(table(T_all.combined$mouse_ID))
Idents(T_all.combined) <- "condition"
m564.combined<-subset(T_all.combined, idents='m564')
AID.combined<-subset(T_all.combined, idents='AID')
AID.mouseID<-metadata[metadata$condition=="AID",]$mouse_ID
m564.mouseID<-metadata[metadata$condition=="m564",]$mouse_ID

##Condition DE
Idents(T_all.combined) <- "condition" #setting idents to condition metadata
T_all.combined.autoimmune.response <- FindMarkers(T_all.combined, ident.1 = "m564", ident.2 = "AID", min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
avg.T_all.combined <- as.data.frame(log1p(AverageExpression(T_all.combined)$RNA))
avg.T_all.combined$gene <- rownames(avg.T_all.combined)
top30 <- head(T_all.combined.autoimmune.response, n = 30)
write.csv(T_all.combined.autoimmune.response ,"564vsAID.csv")
T_all.combined$celltype.condition <- paste(T_all.combined$my.clusters, T_all.combined$condition, sep = "_") #adding metadata identifier


#condition DE (by mouse)
Idents(T_all.combined) <- "mouse_ID"
avg.T_all.ms <- as.data.frame(log1p(AverageExpression(T_all.combined)$RNA))
avg.T_all.ms$gene <- rownames(avg.T_all.ms)
avg.T_all.ms$AID.avg<-rowMeans(avg.T_all.ms[,AID.mouseID])
avg.T_all.ms$m564.avg<-rowMeans(avg.T_all.ms[,m564.mouseID])
avg.T_all.ms$log2FC<-log2(exp(avg.T_all.ms$m564.avg-avg.T_all.ms$AID.avg))
for(i in 1:nrow(avg.T_all.ms)) {
  avg.T_all.ms$ttest[i] <- my.ttest(avg.T_all.ms[i,AID.mouseID],avg.T_all.ms[i,m564.mouseID])
}


###Within Cluster DE
for (i in clusters){
  Idents(T_all.combined) <- "celltype.condition" #setting idents to new metadata column
  df <- FindMarkers(T_all.combined, ident.1 = paste0(i,"_m564"), ident.2 = paste0(i,"_AID"), 
                    min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
  df2 <- head(df, n = 30)
  assign(paste0("cluster",i,"top30"),df2)
  write.csv(df,paste0("cluster",i,".autoimmune.markers.csv"))
  assign(paste0("cluster",i,".autoimmune.response"),df)
  Idents(T_all.combined) <- "my.clusters"
  temp<- subset(T_all.combined, idents = i)
  Idents(temp) <- "condition"
  df <- as.data.frame(log1p(AverageExpression(temp)$RNA))
  df$gene <- rownames(df)
  assign(paste0("avg.cluster",i),df)
  assign(paste0("cluster",i),temp)
}
save.image("analyzed.RData")
