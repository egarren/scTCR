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
library(immunarch)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(qgraph)
library(dplyr)
library(tidyverse)
gm_mean = function(x, na.rm=TRUE){  #geometric mean functions for graphing
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
} 
quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n),na.rm=T)
  breaks[!duplicated(breaks)]
}


load("clone_ba.data.RData")


for(i in c("GLIPH2","TCRdist",'GIANA',"CDR3b","deepTCR","CDR3ba")){
  if(i=="GLIPH2"){
    res<-read.csv("P7436_LL1J9NOFBB.csv",header=T)
    res$index<-paste0("GLIPH_",res$index)
    res<-res[,1:19]
    res <- res %>% 
      filter((grepl("global", type) & (Fisher_score < 0.05) & (pattern != "single")))  # only keep significant global clusters
    length(unique(res$index[grepl("global", res$type)]))
    res$cdr3_ba<-paste0(res$TcRb,"+",res$TcRa)
    res$TCRgroup<-res$index
    cells.res<-left_join(cells,res,by="cdr3_ba")
  }
  if(i=="TCRdist"){
    res<-read.csv("nn_vs_k.csv",header=T)
    res<-unique(res[,c("cdr3_b_aa","clone_id")])
    colnames(res)<-c("cdr3.b","TCRdist")
    res$TCRdist<-paste0("TCRdist_",res$TCRdist)
    res$TCRgroup<-res$TCRdist
    cells.res<-left_join(cells,res,by="cdr3.b")
  }
  if(i=="GIANA"){
    res<-read.table('tutorial--RotationEncodingBL62.txt', header=F,sep='\t',stringsAsFactors = F)
    res<-res[,1:2]
    colnames(res)<-c("cdr3.b","GIANA")
    res$GIANA<-paste0("GIANA_",res$GIANA)
    res$TCRgroup<-res$GIANA
    cells.res<-left_join(cells,res,by="cdr3.b")
  }
  if(i=="CDR3b"){
    cells.res<-cells
    cells.res$TCRgroup<-cells.res$cdr3.b
  }
  if(i=="CDR3ba"){
    cells.res<-cells
    cells.res$TCRgroup<-cells.res$cdr3_ba
  }
  if(i=="deepTCR"){
    cluster.df<-read.csv("./deepTCR/results/unsupervised.rep.clusters.csv",header=T)
    cluster.df$mouse_ID<-gsub("\\..*","",cluster.df$Sample)
    cluster.df$clone_ID<-paste0(cluster.df$Beta_Sequences,"+",cluster.df$Alpha_Sequences,"+",cluster.df$mouse_ID)
    cluster.df$TCRgroup<-paste0("DeepTCR_",cluster.df$cluster)
    cells2<-cells
    cells2$clone_ID<-paste0(cells$cdr3_ba,"+",cells2$mouse_ID)
    cells.res<-left_join(cells2,cluster.df[,c("clone_ID","TCRgroup")],by="clone_ID",keep=F)

  }
  save(list=c("cells.res"),file=paste0(i,".clone_ba.data.RData"))
  
  ####Network plot
  ###CDR3 Network Analysis by condition
  df<-cells.res
  temp_count<-as.data.frame(table(df$TCRgroup))
  colnames(temp_count)<-c("TCRgroup","TCRgroup_count")
  df<-left_join(df,temp_count,by="TCRgroup",keep=F)
  df<-df[df$TCRgroup_count>5,]
  # extract link between diag and patients ID for network
  tmp1=data.frame(Diagnosis=df$condition,ID=df$mouse_ID,stringsAsFactors = F)
  tmp1=tmp1[which(duplicated(paste(tmp1[,1],tmp1[,2],sep="_"))==F),]
  colnames(tmp1)=c("x","y")
  # extract link between clonotype_id.long and patients ID for network
  tmp2=data.frame(df$mouse_ID,df$TCRgroup,stringsAsFactors = F)
  tmp2=tmp2[which(duplicated(paste(tmp2[,1],tmp2[,2],sep="_"))==F),]
  colnames(tmp2)=c("x","y")
  # create table for qgraph
  toQgraph=rbind(tmp1,tmp2)
  # color of Dx
  l=unique(c(toQgraph[,1],toQgraph[,2]))
  col.tmp=rep(adjustcolor("grey85",.3),length(l))
  col.tmp[which(l %in% c("AID"))]<-adjustcolor(hue_pal()(2)[1],.8)
  col.tmp[which(l %in% c("m564"))]<-adjustcolor(hue_pal()(2)[2],.8)
  # color of subjects
  col.tmp[which(l %in% unique(df$mouse_ID[which(df$condition=="AID")]))]<-adjustcolor(hue_pal()(2)[1],.3)
  col.tmp[which(l %in% unique(df$mouse_ID[which(df$condition=="m564")]))]<-adjustcolor(hue_pal()(2)[2],.3)
  # size of nodes
  size.tmp=rep(1,length(l))
  size.tmp[which(l %in% c("AID","m564"))]<-15
  size.tmp[which(l %in% c(names(table(df$mouse_ID))))]<-10
  for(j in which(size.tmp==1)){
    size.tmp[j]=log10(gm_mean(df$TCRgroup_count[which(df$TCRgroup==l[j])]))
  }
  # color of edges
  line.tmp=rep("grey",nrow(toQgraph))
  line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("AID"))]<-hue_pal()(2)[1]
  line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("m564"))]<-hue_pal()(2)[2]
  labels.cex.tmp=l
  labels.cex.tmp=rep(0.00000000001,length(l))
  labels.cex.tmp[which(l %in% c("AID","m564"))]<-1.5
  labels.cex.tmp[which(l %in% c(names(table(df$mouse_ID))))]<-1
  #graph
  dim(toQgraph)
  toQgraph$x[toQgraph$x %in% "AID"]<-"WT"
  toQgraph$x[toQgraph$x %in% "m564"]<-"564Igi"
  for(k in seq(0.85,1.1,0.02)){
    png(paste0(i,".rep",k,".network.by.mouse.png"),width=5,height=5,units="in",res=400)
    qgraph(toQgraph,color=col.tmp,vsize=size.tmp/2,edge.color=line.tmp,edge.width=0.25,
           labels=F,label.cex=labels.cex.tmp,directed=F,repulsion=k)
    dev.off()
  }
  
  
  ####Between Mice comparisons
  clone_sizes<-as.data.frame(table(cells.res$TCRgroup))
  clone_sizes$Freq<-log10(clone_sizes$Freq)
  colnames(clone_sizes)<-c("TCRgroup","log(Size)")
  all_clonotypes<-unique(cells.res$TCRgroup)
  
  #adding nsubject
  cells.res$number_subject<-NA
  for(j in 1:nrow(cells.res)){
    ID<-cells.res$TCRgroup[j]
    cells.res$number_subject[j]<-length(unique(cells.res$mouse_ID[cells.res$TCRgroup %in% ID]))
  }
  
  #adding n cdr3
  cells.res$number_unique_cdr3<-NA
  for(j in 1:nrow(cells.res)){
    ID<-cells.res$TCRgroup[j]
    cells.res$number_unique_cdr3[j]<-length(unique(cells.res$cdr3_ba[cells.res$TCRgroup %in% ID]))
  }
  
  ##EXP vs EXP plot
  df<-cells.res
  imm.list<-list()
  for(j in names(table(df$mouse_ID))){ #creating immunarch list from collapsed data
    df2<-df[df$mouse_ID==j&!is.na(df$TCRgroup),]
    df2<-transform(df2,Clones=ave(seq(nrow(df2)),TCRgroup,FUN=length))
    df2<-unique(df2[,c("TCRgroup","Clones")])
    df2$Proportion<-df2$Clones/sum(df2$Clones)
    imm.list[[j]]<-tibble(Clones=df2$Clones, Proportion=df2$Proportion, CDR3.aa=df2$TCRgroup )
  }
  immdata = repLoad("./rep")
  immdata$meta$Sample<-immdata$meta$mouse_ID
  pr = pubRep(imm.list, "aa", .coding = T, .verbose = F)
  pr.AID = pubRepFilter(pr, immdata$meta, c(condition = "AID"))
  pr.564 = pubRepFilter(pr, immdata$meta, c(condition = "m564"))
  pr.AID[is.na(pr.AID)]<-0
  pr.564[is.na(pr.564)]<-0
  pr.AID[["avgfreq.AID"]] = rowMeans(public_matrix(pr.AID), na.rm = T)
  pr.564[["avgfreq.564"]] = rowMeans(public_matrix(pr.564), na.rm = T)
  pr.AID[["sum.AID"]] = rowSums(public_matrix(pr.AID)[,1:5], na.rm = T)+1
  pr.564[["sum.564"]] = rowSums(public_matrix(pr.564)[,1:5], na.rm = T)+1
  pr.res.GLIPH = dplyr::full_join(as.data.frame(pr.AID), as.data.frame(pr.564), by = c("CDR3.aa")) #,"J.name"
  pr.res.GLIPH[is.na(pr.res.GLIPH)]<-0
  pr.res.GLIPH$sum.AID[pr.res.GLIPH$sum.AID==0]<-1
  pr.res.GLIPH$sum.564[pr.res.GLIPH$sum.564==0]<-1
  pr.res.GLIPH[["Samples.sum"]] = pr.res.GLIPH[["Samples.x"]] + pr.res.GLIPH[["Samples.y"]]
  pr.res.GLIPH[["freq.ratio"]] = apply(pr.res.GLIPH[, c("avgfreq.AID", "avgfreq.564")],1, function(x) log10(x[1])/log10(x[2]))
  pr.res.GLIPH[["log2FC"]]= apply(pr.res.GLIPH[, c("sum.AID", "sum.564")],1, function(x) log2((x[2])/(x[1])))
  pr.res.GLIPH<-left_join(pr.res.GLIPH,unique(cells.res[,c("TCRgroup","number_subject","number_unique_cdr3")]), by=c("CDR3.aa"="TCRgroup"))
  pr.res.GLIPH<-unique(pr.res.GLIPH[,c("sum.564","sum.AID","CDR3.aa","number_subject","number_unique_cdr3","log2FC")])
  labels<-pr.res.GLIPH[abs(pr.res.GLIPH$log2FC)>5|pr.res.GLIPH$sum.564>100,]#,]
  labels<-unique(labels[,c("CDR3.aa","sum.564","sum.AID")])
  p<-ggplot(pr.res.GLIPH,aes(x = sum.564, y =  sum.AID,size=number_subject))+theme_classic()+
    scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(1, 10))+
    geom_point(alpha=0.25,aes(colour=number_unique_cdr3))+ #tfr.score
    guides(fill = guide_legend(override.aes = list(size = 7)))+#,color=F)+ #,size=F
    labs(x = "564Igi", y = "B6", size="Samples",color="Clones")+
    theme(legend.direction = "vertical", legend.box = "horizontal")
  if(i!="CDR3ba"){
    p<-p+scale_color_gradient(low="grey",high="darkgreen",trans="log",breaks=c(1,5,10,15,25),labels=format(c(1,5,10,15,25)))#
  }else{
    p<-p+scale_color_gradient(low="grey",high="darkgreen",trans="log")
  }
  p
  ggsave2(paste0(i,".overlap.scatter.png"),width=6, height=4,device="png")
  
  #clonotype sharing between mice heatmap (all clones)
  clone_size_tab <- table(cells.res$TCRgroup,cells.res$mouse_ID)
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
  annot.row<-clone_sizes[clone_sizes$TCRgroup %in% rownames(plt_mtx_scale_cap), ] %>% remove_rownames %>% column_to_rownames(var="TCRgroup")
  p <- pheatmap(plt_mtx_scale_cap, 
                cluster_rows = TRUE,cluster_cols = TRUE,clustering_method = "ward.D2",show_rownames = F, treeheight_row=3,treeheight_col=3,
                color = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
                fontsize = 8,fontsize_col = 8,fontsize_row=8,
                annotation_col=annot.col,annotation_names_col=F,annotation_names_row=F,
                annotation_colors=list(BMChimera=c("B6"=hue_pal()(2)[1],"564Igi"=hue_pal()(2)[2])),
                annotation_row=annot.row,
                silent = T)
  png(paste0(i,".mice.heatmap.all.png"), width = 5, height = 7, res = 200,units="in")
  grid.arrange(p$gtable)
  dev.off()
  p <- pheatmap(plt_mtx_scale_cap, 
                cluster_rows = TRUE,cluster_cols = TRUE,clustering_method = "ward.D2",show_rownames = F,show_colnames=F, treeheight_row=3,treeheight_col=3,
                color = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
                fontsize = 8,fontsize_col = 8,fontsize_row=8,
                annotation_col=annot.col,annotation_names_col=F,annotation_names_row=F,
                annotation_colors=list(BMChimera=c("B6"=hue_pal()(2)[1],"564Igi"=hue_pal()(2)[2])),
                silent = T)
  png(paste0(i,".mice.heatmap.all2.png"), width = 3.5, height = 3, res = 200,units="in")
  grid.arrange(p$gtable)
  dev.off()
  
  #TRSS cor plot
  #repertoire similarity between mice
  mice_vec<-unique(cells.res$mouse_ID)
  p_mtx_mice <- sapply(mice_vec, function (x) {
    cluster_clonotypes <- cells.res$TCRgroup[cells.res$mouse_ID == x]
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
  png(paste0(i,".mice.TRSS.png"), width = 3, height = 3, res = 200,units="in")
  layout(matrix(c(1,1, rep(2,8)), nrow = 10, ncol = 1, byrow = TRUE))
  # dendrogram
  par(mar=c(0, 2.7, 3.25, 3)) #must adjust to align dendrogram
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
  for(j in c("AID","m564")){
    clone_size_df <- cells.res %>% 
      filter(condition==j) %>%
      group_by(mouse_ID, TCRgroup, condition) %>% 
      tally(name = "clone_size") %>% 
      mutate(size_cat = cut(clone_size, breaks = c(0, 4, 9, 49, 99, 199, Inf), 
                            labels = c("1", "5", "10","50","100","200+"))) %>% 
      mutate(size_cat = factor(size_cat, levels = rev(c("1", "5", "10","50","100","200+"))))
    cluster_size_df <- cells.res %>% 
      filter(condition==j) %>%
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
            strip.text.x = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 12))
    png(paste0(i,".",j,".mice.pie.png"), width = 10, height = 4, res = 200,units="in")
    grid.arrange(p, newpage = T)
    dev.off()
  }
  
  ## expansion score
  scores <- list()
  for (j in unique(cells.res$mouse_ID)) {
    count_mtx <- with(cells.res[cells.res$mouse_ID == j , ],{table(TCRgroup)}) %>% #, condition
      as.numeric() %>% matrix(ncol = 1, byrow = F)
    scores[[j]] <- apply(count_mtx, 2, function(x) {
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
  write.xlsx(df,paste0(i,".exp.scores.mice.xlsx"),showNA=F,row.names=F)
  
  #diversity
  df.list <- list()
  df.cluster<-cells.res
  for(k in c("m564","AID")){
    meta2<-metadata[metadata$condition==k,]
    meta<-data.frame(mouse_ID=names(table(meta2$mouse_ID))) #blank dataframe of mice (to fill in upcoming loop)
    for(l in c("shannon","simpson","invsimpson")){
      meta[[l]]<-0
      for (j in 1:nrow(meta)) {
        ms.clust<-df.cluster[df.cluster$mouse_ID==meta$mouse_ID[j],]
        meta[[l]][j]<-diversity(table(ms.clust$TCRgroup),index=l)
      }
    }
    df.list[[k]]<-meta
  }
  df<-rbindlist(df.list)
  df3<-left_join(as.data.frame(df),metadata[,c("mouse_ID","condition")],keep=F, by="mouse_ID")
  df3<-df3[order(df3$condition),]
  write.xlsx(df3,paste0(i,".div.mice.xlsx"),showNA=F,row.names=F)
  
  
  
  #################################
  
  ###Between cluster comparisons
  cells.res.seurat<-cells.res[!is.na(cells.res$my.clusters2),]
  all_clonotypes.seurat<-unique(cells.res.seurat$TCRgroup)
  clone_sizes.seurat<-as.data.frame(table(cells.res.seurat$TCRgroup))
  clone_sizes.seurat$Freq<-log10(clone_sizes.seurat$Freq)
  colnames(clone_sizes.seurat)<-c("TCRgroup","log(Size)")
  df<-unique(cells.res.seurat[,c("TCRgroup","condition")])
  clone_condition<-data.frame(TCRgroup=unique(df$TCRgroup),condition="B6")
  clone_condition$condition[clone_condition$TCRgroup %in% df$TCRgroup[df$condition=="m564"]]<-"564Igi"
  inter<-intersect(df$TCRgroup[df$condition=="AID"],df$TCRgroup[df$condition=="m564"])
  clone_condition$condition[clone_condition$TCRgroup %in% inter]<-"Both"
  
  #clonotype sharing between clusters heatmap (all clones)
  clone_size_tab <- table(cells.res.seurat$TCRgroup,cells.res.seurat$my.clusters2)
  total_size <- rowSums(clone_size_tab)
  sample_count<-rowSums(clone_size_tab>0)
  # min_total_size <- 2
  plt_mtx <- clone_size_tab
  plt_mtx_scale<-plt_mtx/rowSums(plt_mtx)*100
  plt_mtx_scale_cap <- pmin(plt_mtx_scale, quantile(plt_mtx_scale, 0.94))
  annot.col<-data.frame(Cluster=unique(cells.res.seurat$my.clusters2))
  rownames(annot.col)<-annot.col$Cluster
  # color_key<-as.list(hue_pal()(8))
  # names(color_key)<-names(table(cells.res.seurat$my.clusters2))
  annot.row<-clone_sizes.seurat[clone_sizes.seurat$TCRgroup %in% rownames(plt_mtx_scale_cap), ] %>% remove_rownames %>% column_to_rownames(var="TCRgroup")
  p <- pheatmap(plt_mtx_scale_cap, 
                cluster_rows = TRUE,cluster_cols = TRUE,clustering_method = "ward.D2",show_rownames = F, treeheight_row=5,treeheight_col=2,
                color = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
                # color = colorRampPalette(c("white", "red"))(100),
                fontsize = 8,fontsize_col = 8,fontsize_row=8,
                annotation_col=annot.col,annotation_names_col=F,annotation_names_row=F,
                annotation_colors=list(Cluster=c("TFR"=hue_pal()(8)[1],"Sostdc1"=hue_pal()(8)[2],"TFH-Tcf1"=hue_pal()(8)[3],"TFH-Exhausted"=hue_pal()(8)[4],
                                                 "TFH-Activated"=hue_pal()(8)[5],"TFH-CM"=hue_pal()(8)[6],"TFH-Effector"=hue_pal()(8)[7],"TFH-ISG"=hue_pal()(8)[8])),
                # annotation_row=annot.row,
                silent = T)
  png(paste0(i,".cluster.heatmap.all.png"), width = 5, height = 7, res = 200,units="in")
  grid.arrange(p$gtable)
  dev.off()
  
  #TRSS cor plot
  cells.res.seurat$my.clusters2<-as.character(cells.res.seurat$my.clusters2)
  mice_vec<-unique(cells.res.seurat$my.clusters2)
  p_mtx_mice <- sapply(mice_vec, function (x) {
    cluster_clonotypes <- cells.res.seurat$TCRgroup[cells.res.seurat$my.clusters2 == x]
    p <- table(cluster_clonotypes)[all_clonotypes.seurat]
    p[is.na(p)] <- 0
    p <- p / sum(p)
  })
  sim_scores <- sapply(mice_vec, function(x) {
    sapply(mice_vec, function (y){
      sum(sqrt(p_mtx_mice[,x] * p_mtx_mice[,y]))
    })
  })
  color_vec<-as.data.frame(sort(table(cells.res.seurat$my.clusters2),decreasing=T))
  color_vec$color<-hue_pal()(8)
  bh_dist <- -log10(sim_scores)
  bh_dist <- pmin(bh_dist, 10*max(bh_dist[!is.infinite(bh_dist)]))
  hc <- hclust(as.dist(bh_dist), method = "ward.D")
  sim_scores <- sim_scores[(hc$order), (hc$order)]
  diag(sim_scores) <- NA
  tl_col <- color_vec$color[order(match(color_vec$Var1,rownames(sim_scores)))]
  sim_scores_cap <- pmin(sim_scores, quantile(sim_scores, 0.95, na.rm = T))
  #plot
  png(paste0(i,".clonotype.cluster.TRSS.png"), width = 4, height = 4, res = 200,units="in")
  layout(matrix(c(1,1, rep(2,8)), nrow = 10, ncol = 1, byrow = TRUE))
  # dendrogram
  par(mar=c(0,7.6, 5, 3.8)) #must adjust to align dendrogram
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
  clone_size_df <- cells.res.seurat %>% 
    group_by(my.clusters2, TCRgroup, condition) %>% 
    tally(name = "clone_size") %>% 
    mutate(size_cat = cut(clone_size, breaks = c(0, 4, 9, 49, 99, 199, Inf), 
                          labels = c("1", "5", "10","50","100","200+"))) %>% 
    mutate(size_cat = factor(size_cat, levels = rev(c("1", "5", "10","50","100","200+"))))
  cluster_size_df <- cells.res.seurat %>% 
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
  for (l in stripr) {
    j <- which(grepl('rect', g$grobs[[l]]$grobs[[1]]$childrenOrder))
    g$grobs[[l]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  png(paste0(i,".clonotype.cluster.pie.png"), width = 15, height = 5, res = 200,units="in")
  grid.arrange(g, newpage = T)
  dev.off()
  ## expansion score
  scores <- list()
  for (j in unique(cells.res.seurat$my.clusters2)) {
    count_mtx <- with(cells.res.seurat[cells.res.seurat$my.clusters2 == j , ],{table(TCRgroup, condition)}) %>% 
      as.numeric() %>% matrix(ncol = 2, byrow = F)
    scores[[j]] <- apply(count_mtx, 2, function(x) {
      xx <- x[x > 0]
      p <- xx / sum(xx)
      1 + sum(p * log2(p)) / log2(length(p))
    })
  }
  df<-do.call(rbind.data.frame,scores)
  df$cluster<-names(scores)
  colnames(df)<-c("exp.score.b6","exp.score.564","cluster")
  df2<-as.data.frame(table(cells.res.seurat$my.clusters2))
  df3<-left_join(df2,df,by=c("Var1"="cluster"))
  write.xlsx(df3,paste0(i,".exp.scores.cluster.xlsx"),showNA=F,row.names=F)
  ## expansion score by mouse
  scores <- list()
  for(j in unique(cells.res.seurat$mouse_ID)){
    for (k in unique(cells.res.seurat$my.clusters2)) {
      count_mtx <- with(cells.res.seurat[cells.res.seurat$my.clusters2 == k &cells.res.seurat$mouse_ID == j, ],{table(TCRgroup)}) %>% 
        as.numeric() %>% matrix(ncol = 1, byrow = F)
      scores[[paste0(j,"_",k)]] <- apply(count_mtx, 2, function(x) {
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
  df2<-df2[,c("mouse_ID",names(table(cells.res.seurat$my.clusters2)))]
  df3<-left_join(df2,metadata[,c("mouse_ID","condition")],keep=F, by="mouse_ID")
  df3<-df3[order(df3$condition),]
  write.xlsx(df3,paste0(i,".exp.scores.cluster.mice.xlsx"),showNA=F,row.names=F)
  #diversity
  df.list <- list()
  for (j in unique(cells.res.seurat$my.clusters2)){
    df.cluster<-cells.res.seurat[cells.res.seurat$my.clusters2==j,] #gating on given cluster
    for(k in c("m564","AID")){
      meta2<-metadata[metadata$condition==k,]
      meta<-data.frame(mouse_ID=names(table(meta2$mouse_ID))) #blank dataframe of mice (to fill in upcoming loop)
      # if(k=="AID"){meta$condition<-"WT"}else{meta$condition<-"564Igi"}
      meta$my.clusters2<-j
      for(l in c("shannon","simpson","invsimpson")){
        meta[[l]]<-0
        for (m in 1:nrow(meta)) {
          ms.clust<-df.cluster[df.cluster$mouse_ID==meta$mouse_ID[m],]
          meta[[l]][m]<-diversity(table(ms.clust$TCRgroup),index=l)
        }
      }
      df.list[[paste0(j,k)]]<-meta
    }
  }
  df<-rbindlist(df.list)
  df2<-as.data.frame(spread(subset(df,select=c("mouse_ID","my.clusters2","shannon")), key=my.clusters2, value=shannon))
  df2<-df2[,c("mouse_ID",names(table(cells.res.seurat$my.clusters2)))]
  df3<-left_join(df2,metadata[,c("mouse_ID","condition")],keep=F, by="mouse_ID")
  df3<-df3[order(df3$condition),]
  write.xlsx(df3,paste0(i,".div.cluster.mice.xlsx"),showNA=F,row.names=F)
}


