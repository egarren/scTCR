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
library(ggalluvial)
library(tidyverse)

load("analyzed.RData")

Idents(T_all.combined) <- "my.clusters"
DoHeatmap(object = T_all.combined, features = top4$gene, label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells

#rename clusters
new.cluster.ids <- c("TFR","Sostdc1","TFH-Tcf1","TFH-Exhausted","TFH-Activated","TFH-CM","TFH-Effector","TFH-ISG","NA")
names(new.cluster.ids) <- levels(T_all.combined)
T_all.combined <- RenameIdents(T_all.combined, new.cluster.ids)
T_all.combined@meta.data$my.clusters2 <- Idents(T_all.combined)  
T_all.combined<-subset(T_all.combined,idents=c("TFR","Sostdc1","TFH-Tcf1","TFH-Exhausted","TFH-Activated","TFH-CM","TFH-Effector","TFH-ISG"))
T_all.combined@meta.data$my.clusters2  <- factor(T_all.combined@meta.data$my.clusters2, levels = c("TFR","Sostdc1","TFH-Tcf1","TFH-Exhausted","TFH-Activated","TFH-CM","TFH-Effector","TFH-ISG"))


#QC graphs
p<-VlnPlot(object = T_all.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                                                 "percent.Rpl","HSP.score1","CC.Difference"), 
           combine=F,pt.size=0.001)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]+ NoLegend()  +theme(axis.title=element_blank())
}
cowplot::plot_grid(plotlist = p,ncol=2)
ggsave2("QC.vln.png",width=5, height=10,device="png")
FeatureScatter(object = T_all.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave2("RNA.scatter.png",width=5, height=5,device="png")
FeatureScatter(object = T_all.combined, feature1 = "Cxcr5", feature2 = "Pdcd1")
ggsave2("gene.scatter.png",width=5, height=5,device="png")
# CellScatter(object = T_all.combined, cell1 = "AGTCTACTAGGGTG", cell2 = "CACAGATGGTTTCT")
ElbowPlot(object = T_all.combined)
ggsave2("elbow.png",width=5, height=5,device="png")

# ## tSNE/UMAP
# # Visualization
DimPlot(T_all.combined, reduction = "umap")
DimPlot(T_all.combined, reduction = "umap",pt.size=0.001)+ #NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                   axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
ggsave2("umap.png",width=7, height=5,device="png")
ggsave2("umap.highres.tiff",width=8, height=6,dpi=300,device="tiff")
DimPlot(T_all.combined, reduction = "umap", split.by = "condition")+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
ggsave2("umap2.png",width=12, height=5,device="png")
#by mouse and condition
Idents(T_all.combined) <- "condition"
m564.combined<-subset(T_all.combined, idents='m564')
AID.combined<-subset(T_all.combined, idents='AID')
Idents(m564.combined) <- "my.clusters"
DimPlot(m564.combined, reduction = "umap", split.by = "mouse_ID")+ NoLegend()+NoAxes()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ggsave2("umap.m564.mouse.png",width=10, height=2.25,device="png")
Idents(AID.combined) <- "my.clusters"
DimPlot(AID.combined, reduction = "umap", split.by = "mouse_ID")+ NoLegend()+NoAxes()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ggsave2("umap.AID.mouse.png",width=10, height=2.25,device="png")

## Cluster Analysis
Idents(T_all.combined) <- "my.clusters2"
# DefaultAssay(T_all.combined) <- "integrated" 
T_all.combined<-BuildClusterTree(T_all.combined)
png("cluster.tree.png",width=4,height=6,units="in",res=300)
PlotClusterTree(T_all.combined,font=1)
dev.off()
# # Cluster heatmap
top8 <- T_all.markers %>% group_by(cluster) %>% top_n(8, avg_log2FC)
# DefaultAssay(T_all.combined) <- "integrated"
DoHeatmap(object = T_all.combined, features = top8$gene[1:64], label = F)  #slim.col.label to TRUE prints cluster IDs instead of cells, ,size=5,angle=45,hjust=0.5
ggsave2("cluster.heatmap.png",width=6, height=8,device="png")
DoHeatmap(object = T_all.combined, features = top4$gene[1:32], label = F)  #slim.col.label to TRUE prints cluster IDs instead of cells, ,size=5,angle=45,hjust=0.5
ggsave2("cluster.heatmap2.png",width=6, height=4,device="png")
# DefaultAssay(T_all.combined) <- "RNA"
#UMAP by cluster marker
FeaturePlot(object = T_all.combined, features = "Tbx21", cols = c("grey", "blue"), reduction = "umap")
FeaturePlot(object = T_all.combined, features = "Tbx21", cols = c("grey", "blue"), split.by="condition",reduction = "umap")
gene.list<-c("Foxp3","Tnfsf8","Tcf7","Ext1","Klf2","Sox4","Ccl4","Ifit1")
p<-FeaturePlot(object = T_all.combined, features = gene.list, 
               cols = c("grey", "blue"), reduction = "umap",combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
cowplot::plot_grid(plotlist=p,ncol=4)
ggsave2("umap.clustermarkers.png",width=12,height=6,device="png")
#vln cluster marker
p<-VlnPlot(T_all.combined, features = gene.list,group.by = "my.clusters2",pt.size = 0, combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]+ NoLegend()+theme(axis.title=element_blank())
}
cowplot::plot_grid(plotlist = p,ncol=4)
ggsave2("vlnplot.clustermarkers.png",width=18, height=9,device="png")
ggsave2("vlnplot.clustermarkers.png",width=9, height=4.5,device="png")
#vln cluster marker
p<-RidgePlot(T_all.combined, features = gene.list,group.by = "my.clusters2", combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]+ NoLegend()+theme(axis.title=element_blank())
}
cowplot::plot_grid(plotlist = p,ncol=4)
ggsave2("ridgeplot.clustermarkers.png",width=15, height=6,device="png")
#dot plot
# markers.to.plot <- c("Foxp3","Bcl6","Cxcr5","Cd69","Icos","Ctla4","Cd74","Klf2","Sox4","Ccl5","Nkg7")
# DotPlot(T_all.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
#         split.by = "condition") + RotatedAxis()
# top2 <- T_all.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top2$gene<-make.unique(top2$gene,sep="--")
DotPlot(T_all.combined, features = rev(top2$gene), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "condition")+ RotatedAxis()
top6$gene<-make.unique(top6$gene,sep="--")
DotPlot(T_all.combined, features = top6$gene)+ 
  coord_flip()+theme(legend.title = element_blank(),axis.title=element_blank())+ RotatedAxis()
ggsave2("dotplot.png",width=6, height=12,device="png")
raw.value<-ggplot_build(DotPlot(T_all.combined, features = top6$gene)+ 
                          coord_flip()+theme(legend.title = element_blank(),axis.title=element_blank())+ RotatedAxis())$data
write.csv(raw.value, file="rawvalues.csv")
# top1 <- T_all.markers %>% group_by(cluster) %>% top_n(1, avg_logFC)

# Cluster frequency comparison
ggplot(T_all.combined@meta.data, aes(x = my.clusters,fill = condition)) +  
  geom_bar(aes(y = (..count..)/sum(..count..)),position='dodge')+ 
  scale_y_continuous(labels = percent)+ 
  labs(y="Frequency", x = "Cluster")
ggsave2("clusterfreq.png",width=4, height=2.5,device="png")
cluster.Counts <- table(T_all.combined@meta.data$my.clusters2,T_all.combined@meta.data$condition)
cluster.prop <- as.data.frame(scale(cluster.Counts,scale=colSums(cluster.Counts),center=FALSE)*100) 
ggplot(cluster.prop, aes(fill=Var1,y=Freq, x=Var2,alluvium=Var1,stratum=Var1)) + 
  geom_lode()+geom_flow()+geom_stratum(alpha=0) +theme_classic()+
  theme(legend.title = element_blank(),axis.title=element_blank())+ 
  scale_x_discrete(labels= c("WT","564Igi"))#+ theme(legend.position = "none")
ggsave2("cluster.prop.flow.png",width=4, height=4,device="png")
write.csv(cluster.prop,file="rawvalue.csv")
#cluster freq by mouse
mice<-dimnames(table(T_all.combined@meta.data$mouse_ID))[[1]]
freqlist=list()
for (i in mice) {
  freq<-as.data.frame(table(T_all.combined@meta.data[T_all.combined@meta.data$mouse_ID==i,]$my.clusters2, 
                            T_all.combined@meta.data[T_all.combined@meta.data$mouse_ID==i,]$condition))
  freq$pct<-100*(freq$Freq/sum(freq$Freq))
  freq$mouse<-i
  freqlist[[i]]<-freq
}
cluster.freq.mice<-do.call(rbind,freqlist)
setnames(cluster.freq.mice, old=c("Var1","Var2"), new=c("cluster","condition"))
cluster.freq.mice$condition<-as.character(cluster.freq.mice$condition)
cluster.freq.mice$condition[cluster.freq.mice$condition%in% c("AID")]<-"WT"
cluster.freq.mice$condition[cluster.freq.mice$condition%in% c("m564")]<-"564Igi"
cluster.freq.mice$condition<-as.factor(cluster.freq.mice$condition)
cluster.freq.mice$condition <- factor(cluster.freq.mice$condition, levels = c("WT","564Igi"))
ggbarplot(cluster.freq.mice, x = "cluster", y = "pct",add = c("mean_se", "jitter"),palette=c("black","red"),
          color = "condition",xlab=F,position = position_dodge(0.8), legend="right")+
  stat_compare_means(aes(group = condition),label = "p.format", method="t.test",size=2,label.y=60)+#stat_compare_means(method="anova",label.y=0)+
  labs(y = "Frequency")+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))
ggsave2("cluster.freq.mice.png",width=4.5, height=3,device="png")


#Volcano plot
res<-lapply(ls(pattern = ".autoimmune.response"), get)
names(res)<-ls(pattern = ".autoimmune.response")
thresh_p_val_adj <- 0.05
thresh_lfc <- 0.2
plt_ls <- lapply(res, function(x) {
  x %>% rownames_to_column(var = "gene") %>% 
    filter(!(grepl("Rps", gene) | grepl("Rpl", gene)| grepl("mt.", gene)| grepl("H2.", gene))) %>% 
    mutate(up_in = ifelse(p_val_adj >= thresh_p_val_adj, "NS", ifelse(avg_log2FC > thresh_lfc, "564Igi", ifelse(avg_log2FC < -thresh_lfc, "B6", "NS"))),
           rank_pval = rank(p_val_adj), rank_lfc_inc = rank(avg_log2FC), rank_lfc_dec = rank(-abs(avg_log2FC)),
           lp = -log10(p_val_adj)) %>% 
    arrange(-abs(avg_log2FC))
})
graphs<-c(ls(pattern = ".autoimmune.response")[1:4])#,"T_all.combined.autoimmune.response"
plist <- lapply(graphs, function(x){
  plt_df <- plt_ls[[x]]
  if (sum(plt_df$up_in != "NS") >= 30) {
    plt_df <- plt_df %>% mutate(gene_label = ifelse(up_in != "NS" & (rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), gene, NA))
  } else {
    plt_df <- plt_df %>% mutate(gene_label = ifelse((up_in != "NS" | rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), gene, NA))
  }
  if (x == "cluster0.autoimmune.response") {
    plt_df$lp <- pmin(plt_df$lp, 50)
    plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 0.5)
  } else if (x == "cluster1.autoimmune.response") {
    plt_df$lp <- pmin(plt_df$lp, 100)
    plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 0.5)
  } else if (x == "cluster2.autoimmune.response") {
    plt_df$lp <- pmin(plt_df$lp, 50)
    plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 1)
  } else if (x == "cluster3.autoimmune.response") {
    plt_df$lp <- pmin(plt_df$lp, 50)
    plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 1)
  } else if (x == "T_all.combined.autoimmune.response") {
    plt_df$lp <- pmin(plt_df$lp, 100)
    plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 1)
  } else {
    stop("Invalid tissue!!!")
  }
  
  p <- ggplot(plt_df, aes(x = avg_log2FC, y = lp, color = up_in, label = gene_label)) +
    geom_point(size = 0.5) +
    geom_text_repel(size = 2, segment.size = 0.1, seed = 1, show.legend=F) +
    scale_color_manual(values = c("564Igi" = hue_pal()(2)[2], "B6" = hue_pal()(2)[1], "NS" = "grey80")) +
    scale_x_continuous(limits = c(-max(plt_df$avg_log2FC), max(plt_df$avg_log2FC)), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, max(plt_df$lp)), expand = expansion(mult = c(0.05, 0.01))) +
    # labs(title = x, caption = paste0("p_val_adj < ", thresh_p_val_adj, ", fold change > ", sprintf("%.2f", 2^thresh_lfc))) +
    labs(title = x, color = "Up-regulated in:", y = NULL, x =NULL) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),legend.position="bottom",legend.title=element_text(size=8),
          plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 10, color = "black"))
  return(p)
})
ggarrange(plotlist=lapply(plist, function(x){x + theme(legend.position = "none")}),nrow=1,common.legend=TRUE,legend="bottom")
ggsave2("DEG.png",width=8, height=2.33,device="png")

#MA plot
T_all.combined@meta.data$all<-"all"
Idents(T_all.combined) <- "all" #setting idents to condition metadata
df <- as.data.frame(log1p(AverageExpression(T_all.combined)$RNA))
colnames(df)<-c("logCPM")
df$gene <- rownames(df)
df2<-get("T_all.combined.autoimmune.response")
df2$gene<-rownames(df2)
df3<-left_join(df,df2,by="gene")
rownames(df3)<-df3$gene
glist <- c("Def8","Id3")
plt_df3 <- df3 %>%
  # rownames_to_column(var = "gene") %>%
  filter(!(grepl("Rps", gene) | grepl("Rpl", gene))) %>% 
  mutate(up_in = ifelse(avg_log2FC > 0.5 & p_val_adj <0.05, "564Igi",
                        ifelse(avg_log2FC < -0.5 & p_val_adj < 0.05, "B6", "NS")) %>%
           factor(levels = c("B6", "564Igi", "NS"))) %>%
  mutate(label = ifelse((gene %in% glist & up_in != "NS") | abs(avg_log2FC) > 2, gene, NA))
## set colors
col_de <- c(hue_pal()(2), "black")
names(col_de) <- c("B6", "564Igi", "NS")
ggplot(plt_df3, aes(x = logCPM, y = avg_log2FC, color = up_in)) +
  geom_point() +
  geom_label_repel(aes(label = label)) +
  scale_color_manual(values = col_de) +
  labs(x = "Log Average Expression", y = "Log2 Fold Change (564Igi vs B6)") +
  theme_bw() +
  theme(legend.position = "none")


save.image("graphed.RData")

