rm(list=ls())
library(pheatmap)
library(reticulate)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)
library(data.table)
library(RColorBrewer)
library(Seurat)


meta.list<-list()
for(i in c("AID","m564")){
  file.names<-list.files(path=paste0("../scTfh.data/",i))
  if(i=="AID"){cond<-"B6"}else{cond<-"564Igi"}
  meta.list[[i]]<-data.frame(condition2=cond,file.name=file.names)
}
meta<-do.call("rbind",meta.list)

#ML feature heatmap
df<-read.csv("unsupervised.rep.features.csv",header=T,row.names=1)
deep.feats<-t(df)
pheatmap(deep.feats)
annot.col<-data.frame(file.name=colnames(deep.feats))
annot.col<-merge(annot.col,meta,by="file.name")
names(annot.col)[names(annot.col) == "condition2"] <- "BMChimera"
row.names(annot.col) <- annot.col$file.name
annot.col$file.name <- NULL
p<-pheatmap(deep.feats,annotation_col=annot.col,fontsize_col=3,treeheight_row=2,treeheight_col=2,
            annotation_names_col=F,show_colnames=F,show_rownames=F,
         annotation_colors=list(BMChimera=c("B6"=hue_pal()(2)[1],"564Igi"=hue_pal()(2)[2])),silent=T)
png("DeepTCR.unsup.rep.heatmap.png", width = 4, height = 4, res = 200,units="in")
grid.arrange(p$gtable)
dev.off()

##PCA on samples
mtcars.pca <- prcomp(df, center = TRUE,scale. = TRUE)
summary(mtcars.pca)
d1<-meta
rownames(meta)<-meta$file.name
#PCA plots
autoplot(mtcars.pca,data=d1,colour='condition2')+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_color_manual(values=c("B6"=hue_pal()(2)[1],"564Igi"=hue_pal()(2)[2]))
ggsave2("DeepTCR.PCA.unsup.rep.condition.samples.png",width=4.7, height=3.2,device="png")

pd <- import("pandas")
df<- pd$read_pickle("./unsupervised.rep/KNN_sample.pkl")

##UMAP plots
#create metadata info
#add metadata
load("../../clone_ba.data.RData")
df2<-cells
df2$clone_ID<-paste0(df2$cdr3_ba,"+",df2$mouse_ID)
df<-unique(df2[,c("cdr3_ba","condition","mouse_ID","clone_ID")])
df$top.clust<-NA
for(i in 1:nrow(df)){
  ID<-df$clone_ID[i]
  temp.cells<-df2[df2$clone_ID==ID,]
  df$top.clust[i]<-names(sort(table(temp.cells$my.clusters2), decreasing = T))[1]
}
df$top.clust[df$my.clusters==""]<-NA
temp_count<-as.data.frame(table(df2$clone_ID))
colnames(temp_count)<-c("clone_ID","clone_ID_count")
df<-left_join(df,temp_count,by="clone_ID",keep=F)
d1<-df


#load unsupervised umap from DeepTCR
cluster.df<-read.csv("unsupervised.rep.clusters.csv",header=T)
cluster.df$mouse_ID<-gsub("\\..*","",cluster.df$Sample)
cluster.df$clone_ID<-paste0(cluster.df$Beta_Sequences,"+",cluster.df$Alpha_Sequences,"+",cluster.df$mouse_ID)
cluster.df$cluster<-paste0("DeepTCR_",cluster.df$cluster)

df<-read.csv("unsup.rep.umap.csv",header=F)
colnames(df)<-c("x","y","condition","mouse","cdr3b","cdr3a","freq","count")
df$mouse_ID<-gsub("\\..*","",df$mouse)
df$clone_ID<-paste0(df$cdr3b,"+",df$cdr3a,"+",df$mouse_ID)
df$condition2<-"B6"
df$condition2[df$condition=="m564"]<-"564Igi"
df$condition2 <- factor(x = df$condition2, levels = c("B6", "564Igi"))
df<-left_join(df,unique(d1[,c("clone_ID","top.clust","clone_ID_count")]),by="clone_ID")
df<-left_join(df,cluster.df[,c("clone_ID","cluster")],by="clone_ID")

df2<-df[df$clone_ID_count>10,]

ggplot(df2,aes(x, y,size=clone_ID_count,color=condition2))+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="right")+
  geom_point()+
  labs(x="UMAP1",y="UMAP2")+
  scale_color_manual(values=c("B6"=hue_pal()(2)[1],"564Igi"=hue_pal()(2)[2]))
ggsave2("DeepTCR.umap.condition.png",width=4.7, height=3.2,device="png")
ggplot(df2,aes(x, y,size=clone_ID_count,color=condition2))+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="none")+
  geom_point()+
  labs(x="UMAP1",y="UMAP2")+
  scale_color_manual(values=c("B6"=hue_pal()(2)[1],"564Igi"=hue_pal()(2)[2]))
ggsave2("DeepTCR.umap.condition.nolegend.png",width=3.5, height=3.2,device="png")

ggplot(df2,aes(x, y,size=clone_ID_count,color=mouse_ID))+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="right")+
  geom_point()+
  labs(x="UMAP1",y="UMAP2")
ggsave2("DeepTCR.umap.mouse.png",width=4.7, height=3.2,device="png")
ggplot(df2,aes(x, y,size=clone_ID_count,color=mouse_ID))+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="none")+
  geom_point()+
  labs(x="UMAP1",y="UMAP2")
ggsave2("DeepTCR.umap.mouse.nolegend.png",width=3.5, height=3.2,device="png")

ggplot(df2,aes(x, y,size=clone_ID_count,color=top.clust))+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="right")+
  geom_point()+
  scale_color_manual(values=c("TFR"=hue_pal()(8)[1],"Sostdc1"=hue_pal()(8)[2],"TFH-Tcf1"=hue_pal()(8)[3],"TFH-Exhausted"=hue_pal()(8)[4],
                              "TFH-Activated"=hue_pal()(8)[5],"TFH-CM"=hue_pal()(8)[6],"TFH-Effector"=hue_pal()(8)[7],"TFH-ISG"=hue_pal()(8)[8]))+
  labs(x="UMAP1",y="UMAP2")
ggsave2("DeepTCR.umap.clust.png",width=4.7, height=3.2,device="png")
ggplot(df2,aes(x, y,size=clone_ID_count,color=top.clust))+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="none")+
  geom_point()+
  scale_color_manual(values=c("TFR"=hue_pal()(8)[1],"Sostdc1"=hue_pal()(8)[2],"TFH-Tcf1"=hue_pal()(8)[3],"TFH-Exhausted"=hue_pal()(8)[4],
                              "TFH-Activated"=hue_pal()(8)[5],"TFH-CM"=hue_pal()(8)[6],"TFH-Effector"=hue_pal()(8)[7],"TFH-ISG"=hue_pal()(8)[8]))+
  labs(x="UMAP1",y="UMAP2")
ggsave2("DeepTCR.umap.clust.nolegend.png",width=3.5, height=3.2,device="png")

ggplot(df2,aes(x, y,size=clone_ID_count,color=cluster))+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="none")+
  geom_point()+
  labs(x="UMAP1",y="UMAP2")
ggsave2("DeepTCR.umap.knnclust.png",width=3.5, height=3.2,device="png")




