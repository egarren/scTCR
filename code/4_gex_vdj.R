rm(list=ls())
library(cowplot)
library(ggplot2)
library(plyr)
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
library(ggseqlogo)
library(immunarch)
library(cowplot)
library(ggrepel)
library(Seurat)
library(dplyr)
library(tidyverse)
my.ttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
data_concater <- function(x){
  x<- levels(factor(x))
  paste(x, collapse = "+")
}
data_concater2 <- function(x){
  x<- levels(factor(x))
  paste(x[1])
}
data_concater3 <- function(x){
  x<- levels(factor(x))
  if(length(x)>1){paste(x[1:2], collapse = "+")}else{paste(x[1])}
}

load("graphed.RData")
load("clone.data.RData")

#paired data
df<-clone.data
df.TRA<-df[df$chain=="TRA",]
df.TRB<-df[df$chain=="TRB",]
cells<-inner_join(df.TRB[,c("barcode","cdr3","v_gene","j_gene")],
                  df.TRA[,c("barcode","cdr3","v_gene","j_gene")],by="barcode")
colnames(cells)<-gsub(".x",".b",colnames(cells))
colnames(cells)<-gsub(".y",".a",colnames(cells))
# cells<-left_join(cells,unique(clone.data[,c("barcode","mouse_ID","condition","my.clusters2")]),by="barcode")
cells$cdr3_ba<-paste0(cells$cdr3.b,"+",cells$cdr3.a)
clone_sizes<-as.data.frame(table(cells$cdr3_ba))
colnames(clone_sizes)<-c("cdr3_ba","clone.size")
cells<-left_join(cells,clone_sizes,by="cdr3_ba",keep=F) #add in clone size info
cells<-cells[order(cells$clone.size,decreasing=T),] #keep largest clone for duplicate barcodes
cells<-cells[!duplicated(cells$barcode),] #get rid of duplicate barcodes
rownames(cells)<-cells$barcode
T_all.combined@meta.data$orig.barcode<-rownames(T_all.combined@meta.data)
T_all.combined@meta.data$orig.barcode<-substr(T_all.combined@meta.data$orig.barcode,1,nchar(T_all.combined@meta.data$orig.barcode)-2)
rownames(T_all.combined@meta.data)<-T_all.combined@meta.data$orig.barcode
T_all.combined<-AddMetaData(T_all.combined, metadata=cells)
T_all.combined@meta.data$condition2<-"B6"
T_all.combined@meta.data$condition2[T_all.combined@meta.data$condition=="m564"]<-"564Igi"
T_all.combined@meta.data$condition2 <- factor(x = T_all.combined@meta.data$condition2, levels = c("B6", "564Igi"))
head(T_all.combined@meta.data)
sum(is.na(T_all.combined@meta.data$cdr3_ba))

#key objects
T_all.combined <- RenameCells(T_all.combined, new.names = T_all.combined@meta.data$cellID)
Idents(T_all.combined) <- "condition"
m564.combined<-subset(T_all.combined, idents='m564')
AID.combined<-subset(T_all.combined, idents='AID')
# identifying top clones
m564.top.clones<-names(head(sort(table(m564.combined@meta.data$cdr3_ba), decreasing = T), n=9))
AID.top.clones<-names(head(sort(table(AID.combined@meta.data$cdr3_ba), decreasing = T), n=9))
save.image("vdj.RData")

#stats for paper
median(m564.combined@meta.data$nFeature_RNA)
table(m564.combined@meta.data$productive)
median(AID.combined@meta.data$nFeature_RNA)
table(AID.combined@meta.data$productive)


#mapping expanded clones by cdr3aa
Tall.vdj.clones<-T_all.combined

Tall.vdj.clones<-subset(x = T_all.combined, subset = !is.na(cdr3_ba))
FeaturePlot(Tall.vdj.clones, features= "clone.size",pt.size=0.1, order=T,cols=c("grey","darkgreen"))+ #NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),
        axis.title=element_blank(),plot.title=element_blank(),legend.position="none")
ggsave2("umap.combined.clone.size.png",width=3, height=3,device="png")
p<-FeaturePlot(Tall.vdj.clones, features= "clone.size",split.by = "condition2",pt.size=0.1, 
               order=T,cols=c("grey","darkgreen"),combine=F) #NoLegend()+
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] +NoAxes()+#NoLegend()+
    theme(panel.border = element_rect(colour = "black", size=1),plot.title=element_blank(),
          axis.title.y.right=element_blank(),
          axis.line=element_blank())
}
cowplot::plot_grid(p[[1]],p[[2]],ncol=2)
ggsave2("umap.clone.size.png",width=9, height=3.5,device="png")

##Mapping individual clones within condition
plot.list <- list()
#AID
#mapping top clones
Idents(AID.combined) <- "cdr3_ba"
# plot.list <- list()
for (i in unique(x = names(head(sort(table(AID.combined@meta.data$cdr3_ba), decreasing = T), n=4)))) {#
  plot.list[[i]] <- DimPlot(
    object = AID.combined, cols.highlight="forestgreen",
    cells.highlight = Cells(subset(AID.combined, idents=i)),sizes.highlight=0.3) + 
    NoLegend() + NoAxes()+ggtitle(i)+
    theme(plot.title = element_text(size = 5,hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
#564
#mapping top clones
Idents(m564.combined) <- "cdr3_ba"
for (i in unique(x = names(head(sort(table(m564.combined@meta.data$cdr3_ba), decreasing = T), n=4)))) {#
  plot.list[[i]] <- DimPlot(
    object = m564.combined, cols.highlight="forestgreen",
    cells.highlight = Cells(subset(m564.combined, idents=i)),
    sizes.highlight=0.3) + NoLegend() + NoAxes()+ggtitle(i)+
    theme(plot.title = element_text(size = 5,hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
# CombinePlots(plots = plot.list, ncol = 3)
# ggsave2("564.topclones.umap.png",width=6, height=6,device="png")
CombinePlots(plots = plot.list, ncol = 4)
ggsave2("topclones.umap.png",width=10, height=5,device="png")


#mapping selected clones
TCR.select<-c("CASSRDWGGLEYF+CAAGSPNYSNNRLTL","CASSAGLGSSYEQYF+CAAGASSGSWQLIF","CASSEGLGSSYEQYF+CAAGASSGSWQLIF",
              "CASSGLGGYAEQFF+CAPPATNAYKVIF","CASSIGNNNQAPLF+CAMRDTNAYKVIF","CASSLGGYEQYF+CAAEFANKMIF",
              "CASSQGGADQDTQYF+CALSDEDYANKMIF")
Idents(AID.combined) <- "cdr3_ba"
Idents(m564.combined) <- "cdr3_ba"
plot.list <- list()
for (i in TCR.select) {#
   temp.plot<-try(DimPlot(
    object = AID.combined, cols.highlight="forestgreen",
    cells.highlight = Cells(subset(AID.combined, idents=i)),sizes.highlight=1) + 
    NoLegend() + NoAxes()+ggtitle(i)+
    theme(plot.title = element_text(size = 7,hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)))
  if(!is(temp.plot,"try-error")){
    plot.list[[paste0("AID",i)]] <-temp.plot
  }else{
    plot.list[[paste0("AID",i)]] <-DimPlot(
      object = AID.combined,cols="grey",group.by="condition") + 
      NoLegend() + NoAxes()+ggtitle(i)+
      theme(plot.title = element_text(size = 7,hjust = 0.5),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  }
}
for (i in TCR.select) {#
   temp.plot <- try(DimPlot(
     object = m564.combined, cols.highlight="forestgreen",
     cells.highlight = Cells(subset(m564.combined, idents=i)),sizes.highlight=1) + 
       NoLegend() + NoAxes()+ggtitle(i)+
       theme(plot.title = element_text(size = 7,hjust = 0.5),
             panel.border = element_rect(colour = "black", fill=NA, size=0.5)))
   if(!is(temp.plot,"try-error")){
     plot.list[[paste0("m564",i)]] <-temp.plot
   }else{
     plot.list[[paste0("m564",i)]] <-DimPlot(
       object = m564.combined,cols="grey",group.by="condition") + 
       NoLegend() + NoAxes()+ggtitle(i)+
       theme(plot.title = element_text(size = 7,hjust = 0.5),
             panel.border = element_rect(colour = "black", fill=NA, size=0.5))
   } 
}
# CombinePlots(plots = plot.list, ncol = 3)
# ggsave2("564.topclones.umap.png",width=6, height=6,device="png")
CombinePlots(plots = plot.list, ncol = length(TCR.select))
ggsave2("selectclones.umap2.png",width=22, height=6,device="png")


# DE b/w conditions in individual public clones
m564clonefreq<-as.data.frame(table(m564.combined@meta.data$cdr3_ba))
AIDclonefreq<-as.data.frame(table(AID.combined@meta.data$cdr3_ba))
pubclonefreq<-inner_join(m564clonefreq,AIDclonefreq,by="Var1")
public.clones<-unique(pubclonefreq[(pubclonefreq$Freq.x>5|pubclonefreq$Freq.y>5)&pubclonefreq$Freq.x>2&pubclonefreq$Freq.y>2,])
DefaultAssay(T_all.combined) <- "RNA"
Idents(T_all.combined) <- "cellID"
for(i in public.clones$Var1){
  barcodes<-T_all.combined@meta.data$cellID[T_all.combined@meta.data$cdr3_ba %in% i]
  df<-subset(T_all.combined, idents=barcodes)
  Idents(df) <- "condition" #setting idents to condition metadata
  df.autoimmune.response <- try(FindMarkers(df, ident.1 = "m564", ident.2 = "AID",min.pct=0,logfc.threshold = -Inf,test.use = "MAST"))
  if(!is(df.autoimmune.response,"try-error")){
    # df.autoimmune.response$log2FC<-log2(exp(df.autoimmune.response$avg_logFC))
    assign(paste0(i,".cdr3.autoimmune.response"),df.autoimmune.response)
    avg.df <- as.data.frame(log1p(AverageExpression(df)$RNA))
    avg.df$gene <- rownames(df)
    assign(paste0(i,".cdr3.avg.exp"),avg.df)
  }
}
save.image("temp.vdj.analyzed2.RData")




# identifying expanded clones
m564.expanded <- names(table(m564.combined@meta.data$cdr3_ba))[table(m564.combined@meta.data$cdr3_ba) > 20]
m564.unexpanded <- names(table(m564.combined@meta.data$cdr3_ba))[table(m564.combined@meta.data$cdr3_ba) ==1]
AID.expanded <- names(table(AID.combined@meta.data$cdr3_ba))[table(AID.combined@meta.data$cdr3_ba) > 20]
AID.unexpanded <- names(table(AID.combined@meta.data$cdr3_ba))[table(AID.combined@meta.data$cdr3_ba) ==1]

#expanded DE (by condition)
Idents(m564.combined) <- "cdr3_ba"
expandedDE.m564<-FindMarkers(m564.combined, ident.1 = m564.expanded, ident.2 = m564.unexpanded, 
                             min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
expandedDE.m564$gene<-rownames(expandedDE.m564)
# expandedDE.m564$log2FC<-log2(exp(expandedDE.m564$avg_logFC))
Idents(AID.combined) <- "cdr3_ba"
expandedDE.AID<-FindMarkers(AID.combined, ident.1 = AID.expanded, ident.2 = AID.unexpanded, 
                            min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
expandedDE.AID$gene<-rownames(expandedDE.AID)
# expandedDE.AID$log2FC<-log2(exp(expandedDE.AID$avg_logFC))
exp.vs.expandedDE <-merge(expandedDE.m564[,c("gene","avg_log2FC")], expandedDE.AID[,c("gene","avg_log2FC")], by="gene")
rownames(exp.vs.expandedDE)<-exp.vs.expandedDE$gene
save.image("temp.vdj.analyzed3.RData")


#by clusters
for (i in names(table(T_all.combined$my.clusters))){
  #cluster subset
  Idents(m564.combined) <- "my.clusters"
  m564.clust.combined<-subset(m564.combined, idents=i)
  Idents(AID.combined) <- "my.clusters"
  AID.clust.combined<-subset(AID.combined, idents=i)
  # identifying expanded clones
  if(i %in% c("0","1","2")){
    m564.clust.expanded <- names(table(m564.clust.combined@meta.data$cdr3_ba))[table(m564.clust.combined@meta.data$cdr3_ba) > 5]
    AID.clust.expanded <- names(table(AID.clust.combined@meta.data$cdr3_ba))[table(AID.clust.combined@meta.data$cdr3_ba) > 5]
  }else{
    m564.clust.expanded <- names(table(m564.clust.combined@meta.data$cdr3_ba))[table(m564.clust.combined@meta.data$cdr3_ba) > 1]
    AID.clust.expanded <- names(table(AID.clust.combined@meta.data$cdr3_ba))[table(AID.clust.combined@meta.data$cdr3_ba) > 1]
  }
  m564.clust.unexpanded <- names(table(m564.clust.combined@meta.data$cdr3_ba))[table(m564.clust.combined@meta.data$cdr3_ba) ==1]
  AID.clust.unexpanded <- names(table(AID.clust.combined@meta.data$cdr3_ba))[table(AID.clust.combined@meta.data$cdr3_ba) ==1]
  #expanded DE (by condition)
  Idents(m564.clust.combined) <- "cdr3_ba"
  Idents(AID.clust.combined) <- "cdr3_ba"
  expandedDE.m564.clust<-try(FindMarkers(m564.clust.combined, ident.1 = m564.clust.expanded, ident.2 = m564.clust.unexpanded, 
                               min.pct=0,logfc.threshold = -Inf,test.use = "MAST"), silent=TRUE)
  expandedDE.AID.clust<-try(FindMarkers(AID.clust.combined, ident.1 = AID.clust.expanded, ident.2 = AID.clust.unexpanded, 
                                    min.pct=0,logfc.threshold = -Inf,test.use = "MAST"), silent=TRUE)
  if(!is(expandedDE.m564.clust,"try-error")&!is(expandedDE.AID.clust,"try-error")){
    expandedDE.m564.clust$gene<-rownames(expandedDE.m564.clust)
    # expandedDE.m564.clust$log2FC<-log2(exp(expandedDE.m564.clust$avg_logFC))
    expandedDE.AID.clust$gene<-rownames(expandedDE.AID.clust)
    # expandedDE.AID.clust$log2FC<-log2(exp(expandedDE.AID.clust$avg_logFC))
    clust.exp.vs.exp <-merge(expandedDE.m564.clust[,c("gene","avg_log2FC")], expandedDE.AID.clust[,c("gene","avg_log2FC")], by="gene")
    rownames(clust.exp.vs.exp)<-clust.exp.vs.exp$gene
    assign(paste0("clust",i,".expandedDE.m564"),expandedDE.m564.clust)
    assign(paste0("clust",i,".expandedDE.AID"),expandedDE.AID.clust)
    assign(paste0("clust",i,".exp.vs.expandedDE"),clust.exp.vs.exp)
  }
  save.image(paste0("temp.clust.",i,".RData"))
}

load("temp.clust.7.RData")
# save(T_all.combined,file="T_all.RData")
save(list=c(ls(pattern="DE"),ls(pattern=".markers"),ls(pattern=".response"),"clusters"),file="DE.dfs.RData")
save.image("vdj.analyzed.RData")

#Volcano Plots

#publication
res<-lapply(ls(pattern = ".autoimmune.response")[1], get)
names(res)<-ls(pattern = ".autoimmune.response")[1]
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
graphs<-c(ls(pattern = ".autoimmune.response")[1])#,"T_all.combined.autoimmune.response"
for(i in graphs){
  plt_df <- plt_ls[[i]]
  if (sum(plt_df$up_in != "NS") >= 30) {
    plt_df <- plt_df %>% mutate(gene_label = ifelse(up_in != "NS" & (rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), gene, NA))
  } else {
    plt_df <- plt_df %>% mutate(gene_label = ifelse((up_in != "NS" | rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), gene, NA))
  }
  if (i == "cluster0.autoimmune.response") {
    plt_df$lp <- pmin(plt_df$lp, 50)
    plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 0.5)
  } else {
    plt_df$lp <- pmin(plt_df$lp, 50)
    plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 5)
  }
  
 ggplot(plt_df, aes(x = avg_log2FC, y = lp, color = up_in, label = gene_label)) +
    geom_point(size = 0.5) +
    geom_text_repel(size = 2, segment.size = 0.1, seed = 1, show.legend=F) +
    scale_color_manual(values = c("564Igi" = hue_pal()(2)[2], "B6" = hue_pal()(2)[1], "NS" = "grey80")) +
    scale_x_continuous(limits = c(-max(plt_df$avg_log2FC), max(plt_df$avg_log2FC)), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, max(plt_df$lp)), expand = expansion(mult = c(0.05, 0.01))) +
    # labs(title = x, caption = paste0("p_val_adj < ", thresh_p_val_adj, ", fold change > ", sprintf("%.2f", 2^thresh_lfc))) +
    labs(title = i, color = "Up-regulated in:", y = NULL, x =NULL) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),legend.position="bottom",legend.title=element_text(size=8),
          plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 10, color = "black"))
  ggsave2(paste0(i,".volcano.png"),width=3, height=3,device="png")
}

#clone heatmap
public.clones
clone<-public.clones$Var1[1]
Idents(T_all.combined) <- "cellID"
barcodes<-T_all.combined@meta.data$cellID[T_all.combined@meta.data$cdr3_ba %in% clone]
df<-subset(T_all.combined, idents=barcodes)
df@meta.data$condition<-as.character(df@meta.data$condition)
df@meta.data$condition[df@meta.data$condition%in% c("AID")]<-"B6"
df@meta.data$condition[df@meta.data$condition%in% c("m564")]<-"564Igi"
df@meta.data$condition<-as.factor(df@meta.data$condition)
df@meta.data$condition <- factor(df@meta.data$condition, levels = c("B6","564Igi"))
Idents(df) <- "condition" #setting idents to condition metadata
deg.df<-get(paste0(clone,".cdr3.autoimmune.response"))
genes<-rownames(deg.df)[deg.df$p_val_adj<0.1]
genes<-head(rownames(deg.df),n=50)
DoHeatmap(object = df, features = genes, label = F)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2(paste0(clone,".heatmap.png"),width=5, height=4,device="png")


#Scatter Plots
for(i in c(ls(pattern=".cdr3.avg.exp"))){
  df<-get(i)
  if(!is(df,"try-error")){
    colnames(df)[apply(sapply(c("AID","avg_log2FC.y"), function (y) sapply(colnames(df), 
                                                                       function (x) grepl(y, x))), 1, any)]<-"AID"
    colnames(df)[apply(sapply(c("m564","avg_log2FC.x"), function (y) sapply(colnames(df), 
                                                                        function (x) grepl(y, x))), 1, any)]<-"m564"
    df$sig<-"unsig"
    df$sig[abs(df$m564-df$AID)>0.2]<-"DE"
    # up<-rownames(df[abs(df$m564)>0.29&abs(df$AID)>0.29&!is.na(df$gene),])
    # diff<-rownames(df[abs(df$m564-df$AID)>0.3&!is.na(df$gene),])
    gene.list<-rownames(df[abs(df$m564-df$AID)>0.2&!is.na(df$gene),])
    if(length(gene.list)>30){gene.list<-rownames(df[abs(df$m564-df$AID)>0.5&!is.na(df$gene),])}
    if(length(gene.list)>30){gene.list<-rownames(df[abs(df$m564-df$AID)>1&!is.na(df$gene),])}
    # gene.list<-c(up,diff)
    df$sig <- factor(df$sig, levels = c("unsig","DE"))
    p2 <- ggplot(df, aes(m564,AID)) + geom_point(aes(colour=sig),size=0.5)+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) + 
      ggtitle(i)+theme_classic()+
      labs(x="564Igi",y="B6")+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
            axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
      # geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))+#geom_abline(intercept = 0, slope = 1)+
      scale_color_manual(values = c("black","red"))
    if(length(gene.list)>0){p2 <- LabelPoints(plot = p2, points = gene.list, repel = TRUE,xnudge=0,ynudge=0,size=2.5,segment.size=0.1)}
    p2
    ggsave2(paste0(i,".scatter.png"),width=3, height=3,device="png")
  }
}


for(i in c(ls(pattern="exp.vs.expandedDE"),ls(pattern="expandedDE.ms"))){
  df<-get(i)
  if(!is(df,"try-error")){
    colnames(df)[apply(sapply(c("AID","avg_log2FC.y"), function (y) sapply(colnames(df),
                                                                       function (x) grepl(y, x))), 1, any)]<-"AID"
    colnames(df)[apply(sapply(c("m564","avg_log2FC.x"), function (y) sapply(colnames(df),
                                                                        function (x) grepl(y, x))), 1, any)]<-"m564"
    df$sig<-"unsig"
    df$sig[abs(df$m564)>0.1&abs(df$AID)>0.1]<-"co.DE"
    df$sig[abs(df$m564)>0.1&abs(df$AID)<=0.1]<-"564.DE"
    df$sig[abs(df$m564)<=0.1&abs(df$AID)>0.1]<-"AID.DE"
    up<-rownames(df[abs(df$m564)>0.1&abs(df$AID)>0.1&!is.na(df$gene),])
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.2&abs(df$AID)>0.2&!is.na(df$gene),])}
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.3&abs(df$AID)>0.3&!is.na(df$gene),])}
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.4&abs(df$AID)>0.4&!is.na(df$gene),])}
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.5&abs(df$AID)>0.5&!is.na(df$gene),])}
    if(length(up)>10){up<-rownames(df[abs(df$m564)>1&abs(df$AID)>1&!is.na(df$gene),])}
    diff<-rownames(df[abs(df$m564-df$AID)>0.3&!is.na(df$gene),])
    # gene.list<-rownames(df[abs(df$m564-df$AID)>0.1&!is.na(df$gene),])
    if(length(diff)>10){diff<-rownames(df[abs(df$m564-df$AID)>0.5&!is.na(df$gene),])}
    if(length(diff)>10){diff<-rownames(df[abs(df$m564-df$AID)>1&!is.na(df$gene),])}
    if(length(diff)>10){diff<-rownames(df[abs(df$m564-df$AID)>2&!is.na(df$gene),])}
    if(length(diff)>10){diff<-rownames(df[abs(df$m564-df$AID)>3&!is.na(df$gene),])}
    if(length(diff)>10){diff<-rownames(df[abs(df$m564-df$AID)>4&!is.na(df$gene),])}
    if(length(diff)>10){diff<-rownames(df[abs(df$m564-df$AID)>5&!is.na(df$gene),])}
    gene.list<-c(up,diff,"Pdcd1")
    df$sig <- factor(df$sig, levels = c("unsig","co.DE","564.DE","AID.DE"))
    p2 <- ggplot(df, aes(m564,AID)) + geom_point(aes(colour=sig),size=0.5)+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) +
      ggtitle(i)+theme_classic()+
      labs(x=bquote(~Log[2]~ (frac("564Igi expanded","564Igi unexpanded"))),y=bquote(~Log[2]~ (frac("WT expanded","WT unexpanded"))))+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
            axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
      geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))+#geom_abline(intercept = 0, slope = 1)+
      scale_color_manual(values = c("grey","black","red","orange"))
    if(length(gene.list)>0){p2 <- LabelPoints(plot = p2, points = gene.list, repel = TRUE,xnudge=0,ynudge=0,size=2.5,segment.size=0.1)}
    p2
    ggsave2(paste0(i,".scatter.png"),width=3, height=3,device="png")
  }
}
  