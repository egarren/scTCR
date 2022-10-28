library(ape) 
library(RColorBrewer)
library(scales)
library(ggbiplot)
library(ggfortify)
library(cluster)
library(cowplot)
library(ggplot2)
library(viridis)
library(gridExtra)
library(data.table)
library(dplyr)
library(pheatmap)

load("clone_ba.data.RData")

#load TCRdist groups
tr.clone_df<-read.csv("nn_vs_k.csv",header=T)
tr.clone_df[c("count","clone_id","K_neighbors","nsubject")]<-sapply(tr.clone_df[c("count","clone_id","K_neighbors","nsubject")],as.numeric)
tr.clone_df$epitope[tr.clone_df$epitope=="AID"]<-"B6"
tr.clone_df$epitope[tr.clone_df$epitope=="m564"]<-"564Igi"
tr.clone_df$epitope<-factor(tr.clone_df$epitope,levels=c("B6","564Igi"))

#plot probability of generation
p<-ggplot(tr.clone_df, aes(x= pgen_cdr3_b_aa_nlog10, y =K_neighbors,color = nsubject)) +
  geom_point() + xlim(0,20) + scale_color_viridis()+
  facet_wrap('epitope', scales = "free_y")+ ylab('Neighbors') + xlab('-Log10 Pgen')+theme_bw()+
  theme(strip.text = element_text(face = "bold"))
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- hue_pal()(2)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
png("TCRdist.pgen.png", width = 6, height = 3, res = 200,units="in")
grid.arrange(g, newpage = T)
dev.off()

##comparing repertoires
rw_beta<-fread("rw_beta.csv",header=F)
#count number of clones of AID each 564 clone is similar to
near<-rowSums(rw_beta<50)
#heatmap
m<-as.matrix(rw_beta)
p <- pheatmap(m, 
              cluster_rows = TRUE,cluster_cols = TRUE,clustering_method = "ward.D2",
              show_rownames = F, show_colnames=F,treeheight_row=0,treeheight_col=0,
              # color = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
              color = colorRampPalette(c("red","white", "navy"))(100),
              silent = T)
save.image("temp.TCRdist1.3.RData")
png("TCRdist.rep.compare3.png", width = 5, height = 7, res = 200,units="in")
grid.arrange(p$gtable)
dev.off()
#figure out axes
table(cells$condition)
dim(m)


## adding repertoire annotation
d1<-tr.clone_df
df<-unique(cells[,c("cdr3.b","condition")])
clone_condition<-data.frame(cdr3.b=unique(df$cdr3.b),condition="B6")
clone_condition$condition[clone_condition$cdr3.b %in% df$cdr3.b[df$condition=="m564"]]<-"564Igi"
inter<-intersect(df$cdr3.b[df$condition=="AID"],df$cdr3.b[df$condition=="m564"])
clone_condition$condition[clone_condition$cdr3.b %in% inter]<-"Both"
d1<-left_join(unique(d1[,c("cdr3_b_aa","clone_id","nsubject")]),clone_condition,by=c("cdr3_b_aa"="cdr3.b"),keep=F)
d1$cdr3b_TCRdist_id<-paste0(d1$cdr3_b_aa,d1$clone_id)
rownames(d1)<-d1$cdr3b_TCRdist_id
#add meta info
d1$TCRdist<-paste0("TCRdist_",d1$clone_id)
temp_count<-as.data.frame(table(cells$cdr3.b))
colnames(temp_count)<-c("cdr3_b_aa","cdr3b_count")
d1<-left_join(d1,temp_count,by="cdr3_b_aa",keep=F)

##clustering
pw_beta<-fread("pw_beta.csv",header=F)


#distance matrix to PCA object
Mat<-pw_beta
rownames(Mat)<-d1$cdr3b_TCRdist_id

##PCA
mtcars.pca<-cmdscale(Mat)
summary(mtcars.pca)



#add counts
res<-as.data.frame(mtcars.pca)
colnames(res)<-c("MDS1","MDS2")
res$cdr3b_TCRdist_id<-rownames(res)
res<-left_join(res,d1,by="cdr3b_TCRdist_id",keep=F)
res$cdr3b_count<-as.numeric(res$cdr3b_count)

#add most common cluster
res$top.clust<-NA
for(i in 1:nrow(res)){
  ID<-res$cdr3_b_aa[i]
  temp.cells<-cells[cells$cdr3.b==ID,]
  res$top.clust[i]<-names(sort(table(temp.cells$my.clusters2), decreasing = T))[1]
}

res<-res[res$nsubject>1,]
#MDS plots
ggplot(res, aes(x=MDS1,y=MDS2,color=condition,size=cdr3b_count))+geom_point()+theme_bw()+#labs(x="MDS1",y="MDS2")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="right")+
  scale_color_manual(values=c("B6"=hue_pal()(2)[1],"564Igi"=hue_pal()(2)[2], "Both"="grey70" ))
ggsave2("TCRdist.MDS.condition.png",width=4.7, height=3.2,device="png")
ggplot(res, aes(x=MDS1,y=MDS2,color=condition,size=cdr3b_count))+geom_point()+theme_bw()+#labs(x="MDS1",y="MDS2")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="none")+
  scale_color_manual(values=c("B6"=hue_pal()(2)[1],"564Igi"=hue_pal()(2)[2], "Both"="grey70" ))
ggsave2("TCRdist.MDS.condition.nolegend.png",width=3.5, height=3.2,device="png")

ggplot(res, aes(x=MDS1,y=MDS2,color=top.clust,size=cdr3b_count))+geom_point()+theme_bw()+#labs(x="MDS1",y="MDS2")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="right")+
  scale_color_manual(values=c("TFR"=hue_pal()(8)[1],"Sostdc1"=hue_pal()(8)[2],"TFH-Tcf1"=hue_pal()(8)[3],"TFH-Exhausted"=hue_pal()(8)[4],
                              "TFH-Activated"=hue_pal()(8)[5],"TFH-CM"=hue_pal()(8)[6],"TFH-Effector"=hue_pal()(8)[7],"TFH-ISG"=hue_pal()(8)[8]))
ggsave2("TCRdist.MDS.cluster.png",width=4.7, height=3.2,device="png")
ggplot(res, aes(x=MDS1,y=MDS2,color=top.clust,size=cdr3b_count))+geom_point()+theme_bw()+#labs(x="MDS1",y="MDS2")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="none")+
  scale_color_manual(values=c("TFR"=hue_pal()(8)[1],"Sostdc1"=hue_pal()(8)[2],"TFH-Tcf1"=hue_pal()(8)[3],"TFH-Exhausted"=hue_pal()(8)[4],
                              "TFH-Activated"=hue_pal()(8)[5],"TFH-CM"=hue_pal()(8)[6],"TFH-Effector"=hue_pal()(8)[7],"TFH-ISG"=hue_pal()(8)[8]))
ggsave2("TCRdist.MDS.cluster.nolegend.png",width=3.5, height=3.2,device="png")


ggplot(res, aes(x=MDS1,y=MDS2,color=TCRdist,size=cdr3b_count))+geom_point()+theme_bw()+#labs(x="MDS1",y="MDS2")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="none")
ggsave2("TCRdist.MDS.group.png",width=3.5, height=3.2,device="png")


#TREE
#pairwise distance matrix
# Mat.cor=cor(t(ddc.s[, 4:99]))
# Mat.dist=sqrt(1-Mat.cor)
tr = ape::nj(Mat)
save.image("temp.TCRdist3.RData")
# tr$tip.label = gsub('\\.[1-9]{1}$','', tr$tip.label)
## A list of high-contrast colors
# col.gr<-colorRampPalette(c("red","green"))(length(unique(ddc.s[,2])))
# names(col.gr)<-unique(ddc.s[,2])
col.gr<-c(hue_pal()(2),"grey")
names(col.gr)<-c("B6","564Igi","Both")
png("TCRdist.tree.png", width = 3, height = 3, res = 200,units="in")
par(mar=c(0,0,0,0))
plot.phylo(tr, font=0.1, cex=0.1, 
           # tip.color=col.gr[as.character(ddc.s[,2])], ## color the CDR3s according to GIANA group
           tip.color=col.gr[as.character(ddc.s[,3])], ## color the CDR3s according to antigen-specificity
           align.tip.label=TRUE, x.lim=c(0,1.1))
dev.off()


