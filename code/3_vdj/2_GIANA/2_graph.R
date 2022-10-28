library(ape) 
library(dplyr)
library(RColorBrewer)
library(scales)
library(ggbiplot)
library(ggfortify)
library(cluster)
library(cowplot)

load("clone_ba.data.RData")

d1=read.table('tutorial--RotationEncodingBL62.txt', header=F,sep='\t',stringsAsFactors = F)
Mat=read.table('tutorial--RotationEncodingBL62.txt_EncodingMatrix.txt', header=F,sep='\t',stringsAsFactors = F)

## adding antigen annotation
df<-unique(cells[,c("cdr3.b","condition")])
clone_condition<-data.frame(cdr3.b=unique(df$cdr3.b),condition="B6")
clone_condition$condition[clone_condition$cdr3.b %in% df$cdr3.b[df$condition=="m564"]]<-"564Igi"
inter<-intersect(df$cdr3.b[df$condition=="AID"],df$cdr3.b[df$condition=="m564"])
clone_condition$condition[clone_condition$cdr3.b %in% inter]<-"Both"
d1<-left_join(unique(d1[,1:2]),clone_condition,by=c("V1"="cdr3.b"),keep=F)



## adding isometric coordinates
Mat=Mat[,c(1,4:99)]
Mat=unique(Mat)
rownames(Mat)=Mat[,1]
Mat=Mat[,2:97]
d1=cbind(d1, Mat[d1[,1],])
rownames(d1)<-d1[,1]


##PCA
mtcars.pca <- prcomp(d1[,4:99], center = TRUE,scale. = TRUE)
summary(mtcars.pca)
# ggbiplot::ggbiplot(mtcars.pca)
# ggbiplot::ggbiplot(mtcars.pca,ellipse=TRUE,groups=d1$condition)#labels=rownames(d1), 

#add counts
save<-d1
# d1<-save
d1$GIANA<-paste0("GIANA_",d1$V2)
temp_count<-as.data.frame(table(cells$cdr3.b))
colnames(temp_count)<-c("V1","cdr3b_count")
d1<-left_join(d1,temp_count,by="V1",keep=F)


#add most common cluster
d1$top.clust<-NA
for(i in 1:nrow(d1)){
  ID<-d1$V1[i]
  temp.cells<-cells[cells$cdr3.b==ID,]
  d1$top.clust[i]<-names(sort(table(temp.cells$my.clusters2), decreasing = T))[1]
}


#k means plots
autoplot(kmeans(d1[,4:99], 20), data=d1,colour='condition',size="cdr3b_count")+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="right")+
  scale_color_manual(values=c("B6"=hue_pal()(2)[1],"564Igi"=hue_pal()(2)[2], "Both"="grey70" ))
ggsave2("GIANA.PCA.condition.png",width=4.7, height=3.2,device="png")
autoplot(kmeans(d1[,4:99], 20), data=d1,colour='condition',size="cdr3b_count")+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="none")+
  scale_color_manual(values=c("B6"=hue_pal()(2)[1],"564Igi"=hue_pal()(2)[2], "Both"="grey70" ))
ggsave2("GIANA.PCA.condition.nolegend.png",width=3.5, height=3.2,device="png")

autoplot(kmeans(d1[,4:99], 20), data=d1,colour='top.clust',size="cdr3b_count")+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="right")+
  scale_color_manual(values=c("TFR"=hue_pal()(8)[1],"Sostdc1"=hue_pal()(8)[2],"TFH-Tcf1"=hue_pal()(8)[3],"TFH-Exhausted"=hue_pal()(8)[4],
                              "TFH-Activated"=hue_pal()(8)[5],"TFH-CM"=hue_pal()(8)[6],"TFH-Effector"=hue_pal()(8)[7],"TFH-ISG"=hue_pal()(8)[8]))
ggsave2("GIANA.PCA.cluster.png",width=4.7, height=3.2,device="png")
autoplot(kmeans(d1[,4:99], 20), data=d1,colour='top.clust',size="cdr3b_count")+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="none")+
  scale_color_manual(values=c("TFR"=hue_pal()(8)[1],"Sostdc1"=hue_pal()(8)[2],"TFH-Tcf1"=hue_pal()(8)[3],"TFH-Exhausted"=hue_pal()(8)[4],
                              "TFH-Activated"=hue_pal()(8)[5],"TFH-CM"=hue_pal()(8)[6],"TFH-Effector"=hue_pal()(8)[7],"TFH-ISG"=hue_pal()(8)[8]))
ggsave2("GIANA.PCA.cluster.png",width=3.5, height=3.2,device="png")

autoplot(kmeans(d1[,4:99], 20),data=d1,colour='GIANA',size="cdr3b_count")+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="none")
ggsave2("GIANA.PCA.group.png",width=3.5, height=3.2,device="png")


#filter
# ## Choose a subset of TCRs 
select<-head(names(table(cells$cdr3.b)),n=200)
ddc.s<-d1[which(d1[,1] %in% select),]

#TREE
#pairwise distance matrix
Mat.cor=cor(t(ddc.s[, 4:99]))
Mat.dist=sqrt(1-Mat.cor)
tr = ape::nj(Mat.dist)
# tr$tip.label = gsub('\\.[1-9]{1}$','', tr$tip.label)
## A list of high-contrast colors
col.gr<-colorRampPalette(c("red","green"))(length(unique(ddc.s[,2])))
names(col.gr)<-unique(ddc.s[,2])
col.gr<-c(hue_pal()(2),"grey")
names(col.gr)<-c("B6","564Igi","Both")
png("GIANA.tree.png", width = 3, height = 3, res = 200,units="in")
par(mar=c(0,0,0,0))
plot.phylo(tr, font=0.5, cex=0.5, 
           # tip.color=col.gr[as.character(ddc.s[,2])], ## color the CDR3s according to GIANA group
           tip.color=col.gr[as.character(ddc.s[,3])], ## color the CDR3s according to antigen-specificity
           align.tip.label=TRUE, x.lim=c(0,1.1))
dev.off()






