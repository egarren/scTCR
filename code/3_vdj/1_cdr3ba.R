library(dplyr)
library(tibble)
library(tidyr)
library(pheatmap)
library(corrplot)
library(dendextend)
library(ggplot2)
library(gridExtra)
library(scales)
library(RColorBrewer)
library(xlsx)
library(vegan)
library(data.table)
quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n),na.rm=T)
  breaks[!duplicated(breaks)]
}


load("clone_ba.data.RData")

metadata<-read.csv("/home/eha8/metadata/Tcell.metadata.csv", header=T) #define path
metadata[] <- lapply(metadata, as.character)

#load cdr3 data
mice<-metadata$T.sampleID
for (i in mice){
  csv.path<-paste0("/n/scratch3/users/e/eha8/Rsessions/20221022_scTfh/cellranger/",i,"/outs/per_sample_outs/",i,"/vdj_t/filtered_contig_annotations.csv")
  fasta.path<-paste0("/n/scratch3/users/e/eha8/Rsessions/20221022_scTfh/cellranger/",i,"/outs/per_sample_outs/",i,"/vdj_t/filtered_contig.fasta")
  df<-merge(read.csv(csv.path, header=T),read.fasta(fasta.path),by.x="contig_id",by.y="seq.name")
  colnames(df)[colnames(df)=="length"]<-"vdj_length"
  # df$T.sampleID<-i
  assign(paste0(i,".vdj"),df)
}

#rename barcodes
for (i in mice){
  assign(paste0(i,".vdj"),barcoder(get(paste0(i,".vdj")),prefix=paste0(i,"_")))
}
Tall.vdj<-rbind(T1.vdj,T2.vdj,T3.vdj,T4.vdj,T7.vdj,T8.vdj,T11.vdj,T12.vdj,T34.vdj,T36.vdj)
save.image("temp.vdj.only.RData")

load("graphed.RData")
#extract and add Seurat metadata
clone.data<-Tall.vdj[Tall.vdj$productive=="true"&!is.na(Tall.vdj$productive)&Tall.vdj$chain %in% c("TRA","TRB"),] 
clone.data<-transform(clone.data, vab.freq = ave(seq(nrow(clone.data)), v_gene, FUN=length))
T_all.combined@meta.data$orig.barcode<-rownames(T_all.combined@meta.data)
T_all.combined@meta.data$orig.barcode<-substr(T_all.combined@meta.data$orig.barcode,1,nchar(T_all.combined@meta.data$orig.barcode)-2)
clone.data<-left_join(x = clone.data, y = subset(T_all.combined@meta.data, select=-c(ms.barcoder)),by = c("barcode"="orig.barcode"),keep=F)
clone.data<-left_join_replace(clone.data,metadata,by="ms.barcoder")
clone.data$my.clusters2<-factor(clone.data$my.clusters2, levels = c("TFR","Sostdc1","TFH-Tcf1","TFH-Exhausted","TFH-Activated","TFH-CM","TFH-Effector","TFH-ISG"))
clone.data.seurat<-clone.data[!is.na(clone.data$my.clusters),] 
save(list=c(ls(pattern="clone."),"metadata","Tall.vdj","clusters",ls(pattern="mice")),file="clone.data.RData")


df<-clone.data
df.TRA<-df[df$chain=="TRA",]
df.TRB<-df[df$chain=="TRB",]
cells<-inner_join(df.TRB[,c("barcode","cdr3","v_gene","j_gene")],
                  df.TRA[,c("barcode","cdr3","v_gene","j_gene")],by="barcode")
colnames(cells)<-gsub(".x",".b",colnames(cells))
colnames(cells)<-gsub(".y",".a",colnames(cells))
cells<-left_join(cells,unique(clone.data[,c("barcode","mouse_ID","condition","my.clusters2")]),by="barcode")
cells$cdr3_ba<-paste0(cells$cdr3.b,"+",cells$cdr3.a)
save(list=c(ls(pattern="cells"),"metadata","Tall.vdj","clusters",ls(pattern="mice")),file="clone_ba.data.RData")

####Between Mice comparisons
clone_sizes<-as.data.frame(table(cells$cdr3_ba))
colnames(clone_sizes)<-c("cdr3_ba","Clone Size")
all_clonotypes<-unique(cells$cdr3_ba)

#clonotype sharing between mice heatmap (all clones)
clone_size_tab <- table(cells$cdr3_ba,cells$mouse_ID)
total_size <- rowSums(clone_size_tab)
sample_count<-rowSums(clone_size_tab>0)
# min_total_size <- 2
plt_mtx <- clone_size_tab
plt_mtx_scale<-plt_mtx/rowSums(plt_mtx)*100
plt_mtx_scale_cap <- pmin(plt_mtx_scale, quantile(plt_mtx_scale, 0.94))
annot.col<-data.frame(BMChimera=metadata$condition)
annot.col$BMChimera[annot.col$BMChimera %in% "AID"]<-"B6"
annot.col$BMChimera[annot.col$BMChimera %in% "m564"]<-"564Igi"
rownames(annot.col)<-metadata$mouse_ID
annot.row<-clone_sizes[clone_sizes$cdr3_ba %in% rownames(plt_mtx_scale_cap), ] %>% remove_rownames %>% column_to_rownames(var="cdr3_ba")
p <- pheatmap(plt_mtx_scale_cap, 
              cluster_rows = TRUE,cluster_cols = TRUE,clustering_method = "ward.D2",show_rownames = F, treeheight_row=0,treeheight_col=0,
              color = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
              # color = colorRampPalette(c("white", "red"))(100),
              fontsize = 8,fontsize_col = 8,fontsize_row=8,
              annotation_col=annot.col,annotation_names_col=F,annotation_names_row=F,
              annotation_colors=list(BMChimera=c("B6"=hue_pal()(2)[1],"564Igi"=hue_pal()(2)[2])),
              # annotation_row=annot.row,
              silent = T)
png("clonotype.mice.heatmap.all.png", width = 5, height = 7, res = 200,units="in")
grid.arrange(p$gtable)
dev.off()

#clonotype sharing between mice heatmap
clone_size_tab <- table(cells$cdr3_ba,cells$mouse_ID)
total_size <- rowSums(clone_size_tab)
sample_count<-rowSums(clone_size_tab>0)
min_total_size <- 2
plt_mtx <- clone_size_tab[total_size > min_total_size & sample_count>1,]
plt_mtx_scale<-plt_mtx/rowSums(plt_mtx)*100
plt_mtx_scale_cap <- pmin(plt_mtx_scale, quantile(plt_mtx_scale, 0.94))
annot.col<-data.frame(BMChimera=metadata$condition)
annot.col$BMChimera[annot.col$BMChimera %in% "AID"]<-"B6"
annot.col$BMChimera[annot.col$BMChimera %in% "m564"]<-"564Igi"
rownames(annot.col)<-metadata$mouse_ID
annot.row<-clone_sizes[clone_sizes$cdr3_ba %in% rownames(plt_mtx_scale_cap), ] %>% remove_rownames %>% column_to_rownames(var="cdr3_ba")
p <- pheatmap(plt_mtx_scale_cap, 
              cluster_rows = TRUE,cluster_cols = TRUE,clustering_method = "ward.D2",show_rownames = T, 
              color = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
              # color = colorRampPalette(c("white", "red"))(100),
              fontsize = 8,fontsize_col = 8,fontsize_row=8,treeheight_col=3,treeheight_row=3,
              annotation_col=annot.col,annotation_names_col=F,annotation_names_row=F,
              annotation_colors=list(BMChimera=c("B6"=hue_pal()(2)[1],"564Igi"=hue_pal()(2)[2])),
              annotation_row=annot.row,
              silent = T)
png("clonotype.mice.heatmap.png", width = 8, height = 6, res = 200,units="in")
grid.arrange(p$gtable)
dev.off()


#TRSS cor plot
#repertoire similarity between mice
mice_vec<-unique(cells$mouse_ID)
p_mtx_mice <- sapply(mice_vec, function (x) {
  cluster_clonotypes <- cells$cdr3_ba[cells$mouse_ID == x]
  p <- table(cluster_clonotypes)[all_clonotypes]
  p[is.na(p)] <- 0
  p <- p / sum(p)
})
sim_scores <- sapply(mice_vec, function(x) {
  sapply(mice_vec, function (y){
    sum(sqrt(p_mtx_mice[,x] * p_mtx_mice[,y]))
  })
})
color_vec <-unique(metadata[,c("condition","mouse_ID")])
color_vec$color<-hue_pal()(2)[1]
color_vec$color[color_vec$condition=="m564"]<-hue_pal()(2)[2]
bh_dist <- -log10(sim_scores)
bh_dist <- pmin(bh_dist, 10*max(bh_dist[!is.infinite(bh_dist)]))
hc <- hclust(as.dist(bh_dist), method = "ward.D")
sim_scores <- sim_scores[(hc$order), (hc$order)]
diag(sim_scores) <- NA
tl_col <- color_vec$color[order(match(color_vec$mouse_ID,rownames(sim_scores)))]
sim_scores_cap <- pmin(sim_scores, quantile(sim_scores, 0.95, na.rm = T))
#plot
png("clonotype.mice.TRSS.png", width = 3, height = 3, res = 200,units="in")
layout(matrix(c(1,1, rep(2,8)), nrow = 10, ncol = 1, byrow = TRUE))
par(mar=c(0, 2.7, 3.25, 3)) #must adjust to align dendrogram
dend <- as.dendrogram(hc)
dend %>% plot()
corrplot(sim_scores_cap*100,
         is.corr = F,
         type="upper", order="original",method = "circle",
         outline = T,
         col = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
         # color = colorRampPalette(c("white", "red"))(100),
         mar = c(0, 0, 0, 0),
         tl.col = tl_col,
         tl.pos = 'td', #tl.cex = 2,
         na.label = " ")
dev.off()


#size pie charts (condition vs mouseID)
for(i in c("AID","m564")){
  clone_size_df <- cells %>% 
    filter(condition==i) %>%
    group_by(mouse_ID, cdr3_ba, condition) %>% 
    tally(name = "clone_size") %>% 
    mutate(size_cat = cut(clone_size, breaks = c(0, 1, 4, 9, 19, 49, Inf), 
                          labels = c("1", "2", "5", "10","20","50+"))) %>% 
    mutate(size_cat = factor(size_cat, levels = rev(c("1", "2", "5", "10","20","50+"))))
  cluster_size_df <- cells %>% 
    filter(condition==i) %>%
    group_by(mouse_ID) %>% tally(name = "cluster_size")
  cluster_size_df$cluster_size2<-paste0("n=",cluster_size_df$cluster_size)
  p<- ggplot(clone_size_df) +
    geom_bar(aes(x = 1, y = clone_size, fill = size_cat),
             stat = "identity", size = 0, color = NA, position = position_fill()) +
    geom_text(aes(x = 1, y = 0.4, label = cluster_size2), data = cluster_size_df,size=2.5, hjust = -0.7, vjust =6) +#
    coord_polar(theta = "y") +
    facet_grid( ~ mouse_ID, switch = "y") +
    scale_fill_manual(values = rev(brewer.pal(7, "Reds"))[1:6]) +
    ggtitle("") +
    labs(x = "", y = "", fill = "") +
    theme(panel.background = element_rect(fill = "white", color = "black"),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          strip.background=element_rect(fill="white"),
          # strip.text.y.left = element_text(angle = 0),
          # strip.text.x = element_text(size = 15),
          strip.text.x = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12))
  png(paste0(i,".clonotype.mice.pie.png"), width = 10, height = 4, res = 200,units="in")
  grid.arrange(p, newpage = T)
  dev.off()
}

## expansion score
scores <- list()
for (i in unique(cells$mouse_ID)) {
  count_mtx <- with(cells[cells$mouse_ID == i , ],{table(cdr3_ba)}) %>% #, condition
    as.numeric() %>% matrix(ncol = 1, byrow = F)
  scores[[i]] <- apply(count_mtx, 2, function(x) {
    xx <- x[x > 0]
    p <- xx / sum(xx)
    1 + sum(p * log2(p)) / log2(length(p))
  })
}
scores
df<-do.call(rbind.data.frame,scores)
df$mouse_ID<-names(scores)
colnames(df)<-c("exp.score","mouse_ID")
df<-left_join(df,metadata[,c("mouse_ID","condition")],keep=F)
df<-df[order(df$condition),]
write.xlsx(df,"exp.scores.mice.xlsx",showNA=F,row.names=F)

#diversity
df.list <- list()
df.cluster<-cells.seurat
for(k in c("m564","AID")){
  meta2<-metadata[metadata$condition==k,]
  meta<-data.frame(mouse_ID=names(table(meta2$mouse_ID))) #blank dataframe of mice (to fill in upcoming loop)
  for(l in c("shannon","simpson","invsimpson")){
    meta[[l]]<-0
    for (j in 1:nrow(meta)) {
      ms.clust<-df.cluster[df.cluster$mouse_ID==meta$mouse_ID[j],]
      meta[[l]][j]<-diversity(table(ms.clust$cdr3_ba),index=l)
    }
  }
  df.list[[k]]<-meta
}
df<-rbindlist(df.list)
df3<-left_join(as.data.frame(df),metadata[,c("mouse_ID","condition")],keep=F, by="mouse_ID")
df3<-df3[order(df3$condition),]
write.xlsx(df3,"div.mice.xlsx",showNA=F,row.names=F)

#################################

###Between cluster comparisons
cells.seurat<-cells[!is.na(cells$my.clusters2),]
all_clonotypes.seurat<-unique(cells.seurat$cdr3_ba)
clone_sizes.seurat<-as.data.frame(table(cells.seurat$cdr3_ba))
colnames(clone_sizes.seurat)<-c("cdr3_ba","Clone Size")
df<-unique(cells.seurat[,c("cdr3_ba","condition")])
clone_condition<-data.frame(cdr3_ba=unique(df$cdr3_ba),condition="B6")
clone_condition$condition[clone_condition$cdr3_ba %in% df$cdr3_ba[df$condition=="m564"]]<-"564Igi"
inter<-intersect(df$cdr3_ba[df$condition=="AID"],df$cdr3_ba[df$condition=="m564"])
clone_condition$condition[clone_condition$cdr3_ba %in% inter]<-"Both"

#clonotype sharing between clusters heatmap (all clones)
clone_size_tab <- table(cells.seurat$cdr3_ba,cells.seurat$my.clusters2)
total_size <- rowSums(clone_size_tab)
sample_count<-rowSums(clone_size_tab>0)
plt_mtx <- clone_size_tab
plt_mtx_scale<-plt_mtx/rowSums(plt_mtx)*100
plt_mtx_scale_cap <- pmin(plt_mtx_scale, quantile(plt_mtx_scale, 0.94))
annot.col<-data.frame(Cluster=unique(cells.seurat$my.clusters2))
rownames(annot.col)<-annot.col$Cluster
annot.row<-clone_sizes.seurat[clone_sizes.seurat$cdr3_ba %in% rownames(plt_mtx_scale_cap), ] %>% remove_rownames %>% column_to_rownames(var="cdr3_ba")
p <- pheatmap(plt_mtx_scale_cap, 
              cluster_rows = TRUE,cluster_cols = TRUE,clustering_method = "ward.D2",show_rownames = F, treeheight_row=0,treeheight_col=0,
              color = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
              # color = colorRampPalette(c("white", "red"))(100),
              fontsize = 8,fontsize_col = 8,fontsize_row=8,
              annotation_col=annot.col,annotation_names_col=F,annotation_names_row=F,
              annotation_colors=list(Cluster=c("TFR"=hue_pal()(8)[1],"Sostdc1"=hue_pal()(8)[2],"TFH-Tcf1"=hue_pal()(8)[3],"TFH-Exhausted"=hue_pal()(8)[4],
                                               "TFH-Activated"=hue_pal()(8)[5],"TFH-CM"=hue_pal()(8)[6],"TFH-Effector"=hue_pal()(8)[7],"TFH-ISG"=hue_pal()(8)[8])),
              # annotation_row=annot.row,
              silent = T)
png("clonotype.cluster.heatmap.all.png", width = 5, height = 7, res = 200,units="in")
grid.arrange(p$gtable)
dev.off()

#clonotype sharing between cluster heatmap (public clones)
clone_size_tab <- table(cells.seurat$cdr3_ba,cells.seurat$my.clusters2)
total_size <- rowSums(clone_size_tab)
sample_count<-rowSums(clone_size_tab>0)
min_total_size <- 2
plt_mtx <- clone_size_tab[total_size > min_total_size & sample_count>1,]
plt_mtx_scale<-plt_mtx/rowSums(plt_mtx)*100
plt_mtx_scale_cap <- pmin(plt_mtx_scale, quantile(plt_mtx_scale, 0.94))
annot.col<-data.frame(Cluster=unique(cells.seurat$my.clusters2))
rownames(annot.col)<-annot.col$Cluster
annot.row<-clone_sizes.seurat[clone_sizes.seurat$cdr3_ba %in% rownames(plt_mtx_scale_cap), ] 
annot.row<-merge(annot.row,clone_condition,by="cdr3_ba",keep=F)
names(annot.row)[names(annot.row) == "condition"] <- "BMChimera"
annot.row <-annot.row %>% remove_rownames %>% column_to_rownames(var="cdr3_ba")
p <- pheatmap(plt_mtx_scale_cap, 
              cluster_rows = TRUE,cluster_cols = TRUE,clustering_method = "ward.D2",show_rownames = F, 
              color = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
              fontsize = 8,fontsize_col = 8,fontsize_row=8,treeheight_col=5,treeheight_row=5,
              annotation_col=annot.col,annotation_names_col=F,annotation_names_row=F,
              annotation_colors=list(Cluster=c("TFR"=hue_pal()(8)[1],"Sostdc1"=hue_pal()(8)[2],"TFH-Tcf1"=hue_pal()(8)[3],"TFH-Exhausted"=hue_pal()(8)[4],
                                                    "TFH-Activated"=hue_pal()(8)[5],"TFH-CM"=hue_pal()(8)[6],"TFH-Effector"=hue_pal()(8)[7],"TFH-ISG"=hue_pal()(8)[8]),
                                     BMChimera=c("B6"=hue_pal()(2)[1],"564Igi"=hue_pal()(2)[2],"Both"=alpha("grey",0.5))),
              annotation_row=annot.row,
              silent = T)
png("clonotype.cluster.heatmap.png", width = 6, height = 6, res = 200,units="in")
grid.arrange(p$gtable)
dev.off()

#clonotype sharing between cluster heatmap (select clones)
TCR.select<-c("CASSRDWGGLEYF+CAAGSPNYSNNRLTL","CASSAGLGSSYEQYF+CAAGASSGSWQLIF","CASSEGLGSSYEQYF+CAAGASSGSWQLIF",
              "CASSGLGGYAEQFF+CAPPATNAYKVIF","CASSIGNNNQAPLF+CAMRDTNAYKVIF","CASSLGGYEQYF+CAAEFANKMIF",
              "CASSQGGADQDTQYF+CALSDEDYANKMIF")
sum("CASSQGGADQDTQYF" %in% cells$cdr3.b)
cells.NAlabel<-cells
cells.NAlabel$my.clusters2<-as.character(cells.NAlabel$my.clusters2)
cells.NAlabel$my.clusters2[is.na(cells.NAlabel$my.clusters2)]<-"NA"
clone_size_tab <- table(cells.NAlabel$cdr3_ba,cells.NAlabel$my.clusters2)
plt_mtx <- clone_size_tab[rownames(clone_size_tab) %in% TCR.select,]
empty<-colSums(plt_mtx)==0
plt_mtx<-plt_mtx[ , !empty]
plt_mtx<-plt_mtx[,-which(colnames(plt_mtx) %in% c("NA"))]
plt_mtx_scale<-plt_mtx/rowSums(plt_mtx)*100
plt_mtx_scale_cap <- pmin(plt_mtx_scale, quantile(plt_mtx_scale, 0.94))
annot.col<-data.frame(Cluster=unique(cells.NAlabel$my.clusters2))
rownames(annot.col)<-annot.col$Cluster
annot.row<-clone_sizes[clone_sizes$cdr3_ba %in% rownames(plt_mtx_scale_cap), ] 
annot.row<-merge(annot.row,clone_condition,by="cdr3_ba",keep=F)
names(annot.row)[names(annot.row) == "condition"] <- "BMChimera"
annot.row <-annot.row %>% remove_rownames %>% column_to_rownames(var="cdr3_ba")
p <- pheatmap(plt_mtx_scale_cap, 
              cluster_rows = TRUE,cluster_cols = TRUE,clustering_method = "ward.D2",show_rownames = T, 
              color = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
              fontsize = 8,fontsize_col = 8,fontsize_row=8,treeheight_col=5,treeheight_row=5,
              annotation_col=annot.col,annotation_names_col=F,annotation_names_row=F,
              annotation_colors=list(Cluster=c("TFR"=hue_pal()(8)[1],"Sostdc1"=hue_pal()(8)[2],"TFH-Tcf1"=hue_pal()(8)[3],"TFH-Exhausted"=hue_pal()(8)[4],
                                               "TFH-Activated"=hue_pal()(8)[5],"TFH-Effector"=hue_pal()(8)[7],"TFH-ISG"=hue_pal()(8)[8],"NA"="grey"),
                                     BMChimera=c("B6"=hue_pal()(2)[1],"564Igi"=hue_pal()(2)[2],"Both"=alpha("grey",0.5))),
              annotation_row=annot.row,
              silent = T)
png("select.clonotype.cluster.heatmap.png", width = 6, height = 3.5, res = 200,units="in")
grid.arrange(p$gtable)
dev.off()


#TRSS cor plot
cells.seurat$my.clusters2<-as.character(cells.seurat$my.clusters2)
mice_vec<-unique(cells.seurat$my.clusters2)
p_mtx_mice <- sapply(mice_vec, function (x) {
  cluster_clonotypes <- cells.seurat$cdr3_ba[cells.seurat$my.clusters2 == x]
  p <- table(cluster_clonotypes)[all_clonotypes.seurat]
  p[is.na(p)] <- 0
  p <- p / sum(p)
})
sim_scores <- sapply(mice_vec, function(x) {
  sapply(mice_vec, function (y){
    sum(sqrt(p_mtx_mice[,x] * p_mtx_mice[,y]))
  })
})
color_vec<-as.data.frame(sort(table(cells.seurat$my.clusters2),decreasing=T))
color_vec$color<-hue_pal()(8)
bh_dist <- -log10(sim_scores)
bh_dist <- pmin(bh_dist, 10*max(bh_dist[!is.infinite(bh_dist)]))
hc <- hclust(as.dist(bh_dist), method = "ward.D")
sim_scores <- sim_scores[(hc$order), (hc$order)]
diag(sim_scores) <- NA
tl_col <- color_vec$color[order(match(color_vec$Var1,rownames(sim_scores)))]
sim_scores_cap <- pmin(sim_scores, quantile(sim_scores, 0.95, na.rm = T))
#plot
png("clonotype.cluster.TRSS.png", width = 8, height = 6, res = 200,units="in")
layout(matrix(c(1,1, rep(2,8)), nrow = 10, ncol = 1, byrow = TRUE))
par(mar=c(0,15.5, 7, 13.2)) #must adjust to align dendrogram
dend <- as.dendrogram(hc)
dend %>% plot()
## corrplot
corrplot(sim_scores_cap,
         is.corr = F,
         type="upper", order="original",method = "circle",
         outline = T,
         col = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
         mar = c(0, 0, 0, 0),
         tl.col = tl_col,
         tl.pos = 'td', #tl.cex = 2,
         na.label = " ")
dev.off()



#size pie charts (condition vs mouseID)
clone_size_df <- cells.seurat %>% 
  group_by(my.clusters2, cdr3_ba, condition) %>% 
  tally(name = "clone_size") %>% 
  mutate(size_cat = cut(clone_size, breaks = c(0, 1, 4, 9, 19, 49, Inf), 
                        labels = c("1", "2", "5", "10","20","50+"))) %>% 
  mutate(size_cat = factor(size_cat, levels = rev(c("1", "2", "5", "10","20","50+"))))
cluster_size_df <- cells.seurat %>% 
  group_by(my.clusters2, condition) %>% tally(name = "cluster_size")
cluster_size_df$cluster_size2<-paste0("n=",cluster_size_df$cluster_size)
condition.labs<-c("B6","564Igi")
names(condition.labs)<-c("AID","m564")
p <- ggplot(clone_size_df) +
  geom_bar(aes(x = 1, y = clone_size, fill = size_cat),
           stat = "identity", size = 0, color = NA, position = position_fill()) +
  geom_text(aes(x = 1, y = 0.4, label = cluster_size2), data = cluster_size_df,size=2.5, hjust = -0.7, vjust =6) +#
  coord_polar(theta = "y") +
  facet_grid(condition ~ my.clusters2, switch = "y",labeller=labeller(condition=condition.labs)) +#
  scale_fill_manual(values = rev(brewer.pal(7, "Reds"))[1:6]) +
  ggtitle("") +
  labs(x = "", y = "", fill = "") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        strip.background=element_rect(fill="white"),
        strip.text.y.left = element_text(angle = 0),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))
# p
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- hue_pal()(8)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
png("clonotype.cluster.pie.png", width = 15, height = 5, res = 200,units="in")
grid.arrange(g, newpage = T)
dev.off()
## expansion score
scores <- list()
for (i in unique(cells.seurat$my.clusters2)) {
  count_mtx <- with(cells.seurat[cells.seurat$my.clusters2 == i , ],{table(cdr3_ba, condition)}) %>% 
    as.numeric() %>% matrix(ncol = 2, byrow = F)
  scores[[i]] <- apply(count_mtx, 2, function(x) {
    xx <- x[x > 0]
    p <- xx / sum(xx)
    1 + sum(p * log2(p)) / log2(length(p))
  })
}
df<-do.call(rbind.data.frame,scores)
df$cluster<-names(scores)
colnames(df)<-c("exp.score.b6","exp.score.564","cluster")
df2<-as.data.frame(table(cells.seurat$my.clusters2))
df3<-left_join(df2,df,by=c("Var1"="cluster"))
write.xlsx(df3,"exp.scores.cluster.xlsx",showNA=F,row.names=F)
## expansion score by mouse
scores <- list()
for(j in unique(cells.seurat$mouse_ID)){
  for (i in unique(cells.seurat$my.clusters2)) {
    count_mtx <- with(cells.seurat[cells.seurat$my.clusters2 == i &cells.seurat$mouse_ID == j, ],{table(cdr3_ba)}) %>% 
      as.numeric() %>% matrix(ncol = 1, byrow = F)
    scores[[paste0(j,"_",i)]] <- apply(count_mtx, 2, function(x) {
      xx <- x[x > 0]
      p <- xx / sum(xx)
      1 + sum(p * log2(p)) / log2(length(p))
    })
  }
}
df<-do.call(rbind.data.frame,scores)
df$name<-names(scores)
colnames(df)<-c("exp.score","name")
df$mouse_ID<-sub("_.*", "", df$name)
df$cluster<-sub(".*_", "", df$name)
df2<-spread(subset(df,select=-c(name)), key=cluster, value=exp.score)
df2<-df2[,c("mouse_ID",names(table(cells.seurat$my.clusters2)))]
df3<-left_join(df2,metadata[,c("mouse_ID","condition")],keep=F, by="mouse_ID")
df3<-df3[order(df3$condition),]
write.xlsx(df3,"exp.scores.cluster.mice.xlsx",showNA=F,row.names=F)
#diversity
df.list <- list()
for (i in unique(cells.seurat$my.clusters2)){
  df.cluster<-cells.seurat[cells.seurat$my.clusters2==i,] #gating on given cluster
  for(k in c("m564","AID")){
    meta2<-metadata[metadata$condition==k,]
    meta<-data.frame(mouse_ID=names(table(meta2$mouse_ID))) #blank dataframe of mice (to fill in upcoming loop)
    # if(k=="AID"){meta$condition<-"WT"}else{meta$condition<-"564Igi"}
    meta$my.clusters2<-i
    for(l in c("shannon","simpson","invsimpson")){
      meta[[l]]<-0
      for (j in 1:nrow(meta)) {
        ms.clust<-df.cluster[df.cluster$mouse_ID==meta$mouse_ID[j],]
        meta[[l]][j]<-diversity(table(ms.clust$cdr3_ba),index=l)
      }
    }
    df.list[[paste0(i,k)]]<-meta
  }
}
df<-rbindlist(df.list)
df2<-as.data.frame(spread(subset(df,select=c("mouse_ID","my.clusters2","shannon")), key=my.clusters2, value=shannon))
df2<-df2[,c("mouse_ID",names(table(cells.seurat$my.clusters2)))]
df3<-left_join(df2,metadata[,c("mouse_ID","condition")],keep=F, by="mouse_ID")
df3<-df3[order(df3$condition),]
write.xlsx(df3,"div.cluster.mice.xlsx",showNA=F,row.names=F)


####export data
head(cells)
clone_sizes<-as.data.frame(table(cells$cdr3_ba))
colnames(clone_sizes)<-c("cdr3_ba","clone.size")
cells2<-left_join(cells,clone_sizes,by="cdr3_ba",keep=F) #add in clone size info

## TCRdist
# TRgene_list<-list()
# for(i in c("v_gene.a","j_gene.a","v_gene.b","j_gene.b")){
#   vect<-sort(unique(cells[[i]]))
#   length(vect)<-10000
#   TRgene_list[[i]]<-vect
# }
# df<-do.call(cbind.data.frame, TRgene_list)
# write.csv(df,"TRgenes.csv",row.names=F,na="")
# ### must manually match to TR genes in 'alphabeta_gammadelta_db.tsv'
df<-read.csv(file="TRgene.match.csv")
dict<-df[df$cellranger_id!="",c("cellranger_id","id")]
dict$cellranger_id<-paste0(dict$cellranger_id,"*01")
df<-cells[,c("mouse_ID","v_gene.a","j_gene.a","cdr3.a","v_gene.b","j_gene.b","cdr3.b","condition")]
colnames(df)<-c("subject","v_a_gene","j_a_gene","cdr3_a_aa","v_b_gene","j_b_gene","cdr3_b_aa","epitope")
df$count<-1
for(i in c("v_a_gene","j_a_gene","v_b_gene","j_b_gene")){
  df[[i]]<-paste0(df[[i]],"*01")
}
df$v_a_gene<-dplyr::recode(df$v_a_gene,!!!setNames(dict$id, dict$cellranger_id))
imgt.cells<-df
write.csv(df,"dash.csv",row.names=F)

##GIANA
#export all clones
df<-unique(imgt.cells[,c("cdr3_b_aa","v_b_gene","j_b_gene")])
write.table(df, file = "tutorial.txt", sep = "\t",row.names = FALSE,col.names=F,quote=F)
#export clones by condition
for(i in unique(imgt.cells$epitope)){
  df<-imgt.cells[imgt.cells$epitope==i,]
  df<-unique(df[,c("cdr3_b_aa","v_b_gene","j_b_gene")])
  write.table(df, file = paste0(i,"_tutorial.txt"), sep = "\t",row.names = FALSE,col.names=F,quote=F)
}
##use mouse genes (replace GIANA inbuilt IMGT human file)
library(ampir)
df<-read.csv(file="TRgene.match.csv")
df<-df[df$organism=="mouse" & df$chain=="B"& df$region=="V",c("id","aligned_protseq")]
df_to_faa(df, file = "Imgt_Human_TRBV.fasta")


##GLIPH2
df<-unique(cells2[,c("cdr3.b","v_gene.b","j_gene.b","cdr3.a","mouse_ID","clone.size")])
colnames(df)<-c("CDR3b","TRBV","TRBJ","CDR3a","subject.condition","count")
write.table(df,file="GLIPH.txt",sep="\t",row.names=F,col.names=F,quote=F)
#submit file to: http://50.255.35.37:8080/

##DeepTCR
df<-clone.data
df.TRA<-df[df$chain=="TRA",]
df.TRB<-df[df$chain=="TRB",]
cells2<-inner_join(df.TRB[,c("barcode","cdr3","v_gene","j_gene","d_gene")],
                  df.TRA[,c("barcode","cdr3","v_gene","j_gene")],by="barcode")
colnames(cells2)<-gsub(".x",".b",colnames(cells2))
colnames(cells2)<-gsub(".y",".a",colnames(cells2))
cells2<-left_join(cells2,unique(clone.data[,c("barcode","mouse_ID","condition","my.clusters2")]),by="barcode")
cells2$cdr3_ba<-paste0(cells2$cdr3.b,"+",cells2$cdr3.a)

setwd("/n/scratch3/users/e/eha8/Rsessions/20221022_scTfh")
dir.create("./deepTCR")
setwd("./deepTCR")
dir.create("./scTfh.data")
setwd("./scTfh.data")
for(i in names(table(cells2$condition))){
  df<-cells2[cells2$condition==i,]
  dir.create(paste0("./",i))
  setwd(paste0("./",i))
  for(j in names(table(df$mouse_ID))){
    df2<-df[df$mouse_ID==j,]
    msclone_sizes<-as.data.frame(table(df2$cdr3_ba))
    colnames(msclone_sizes)<-c("cdr3_ba","clone.size")
    cells3<-left_join(df2,msclone_sizes,by="cdr3_ba",keep=F) #add in clone size info
    cells4<-unique(cells3[,c("cdr3.b","v_gene.b","j_gene.b","d_gene","cdr3.a","v_gene.a","j_gene.a","mouse_ID","clone.size")])
    cells4$clone.size<-as.numeric(as.character(cells4$clone.size))
    write.table(cells4,file=paste0(j,".deepTCR.tsv"),sep="\t",row.names=F,col.names=F,quote=F)
  }
  setwd("../")
}
setwd("../")

#making count matrix (for ML regression)
df<-as.data.frame.matrix(table(cells$cdr3_ba,cells$condition))
df$cdr3_ba<-rownames(df)
df2<-left_join(df,unique(cells[,c("cdr3_ba","cdr3.b","cdr3.a")]),by="cdr3_ba")
df3<-data.frame(alpha=df2$cdr3.a,beta=df2$cdr3.b,AID=df2$AID,m564=df2$m564)
df3$alpha<-gsub("[^GALMFWKQESPVICYHRNDT]+", "", df3$alpha)
df3$beta<-gsub("[^GALMFWKQESPVICYHRNDT]+", "", df3$beta)
write.csv(df3, "BMchim.counts_regression.csv",row.names=F,col.names=T,quote=F)

