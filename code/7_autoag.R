rm(list=ls())
library(ape) 
library(dplyr)
library(RColorBrewer)
library(scales)
library(ggfortify)
library(cluster)
library(cowplot)
library(readxl)
library(pheatmap)
library(gridExtra)
library(corrr)
library(corrplot)
library(factoextra)
library(purrr)
library(tibble)
library(Seurat)
library(M3C)
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}


## load antigen MFIs
for(i in c("IgG","IgG2c","IgG2a")){ 
  df<-read.csv(paste0("autoag_data_",i,".csv"), check.names=F)
  colnames(df)<-c("mouse",paste0(colnames(df[,-1])," (",i,")"))
  df$mouse<-as.numeric(df$mouse)
  autoag<-left_join(autoag,df,by="mouse")
}

antigens<-colnames(select(autoag,-c(mouse)))

# load metadata
meta<-read.csv("mouse_meta.csv")

#load antigen metadata
autoag.meta<-read.csv("autoag_meta2.csv")
autoag.meta$Ag<-paste0(autoag.meta$Ag," (",autoag.meta$Class,")")

#Differential analysis (wilcoxan)
df<-left_join(autoag,meta,by="mouse")
rownames(df)<-df$mouse
thresh_p.adj <- 0.05
thresh_lfc <- 1
table(df$tx2)
for(j in c("BMchim.564-Icos.RAG-TCRA","BMchim.564-Icos.RAG-TCRB","BMchim.564-Icos.RAG-TCRE","BMchim.564")){
  for(k in c("564homo","564homo.Icos","BMchim.564")){
    if(j!=k){
      #ttest loop
      df2<-df[df$tx2 %in% c(k,j),]
      ttest.pval.df<-data.frame(test=antigens)
      ttest.pval.df$tx2<-NA
      ttest.list<-list()
      for(i in antigens){
        test<-try(wilcox.test(df2[[i]] ~ df2$tx2, data = df2))
        test2<-try(t.test(df2[[i]] ~ df2$tx2, data = df2))
        if(!is(test,"try-error")&!is(test2,"try-error")&!is.na(test$p.value)){
          ttest.list[[i]]<-c(p.val=test$p.value,test2$estimate)
          ttest.pval.df$tx2[ttest.pval.df$test==i]<-test$p.value
        }
      }
      ttest.df<-as.data.frame(purrr::map_dfr(lapply(ttest.list,unlist), ~as_tibble(t(.))))
      rownames(ttest.df)<-names(ttest.list)
      ttest.df$p.adj<-ttest.df$p.val
      ttest.df<-ttest.df[order(ttest.df$p.adj),]
      ttest.df$log2FC<-log2(ttest.df[[3]]/ttest.df[[2]])
      assign(paste0(k,"vs",j,".ttest.df"),ttest.df)
      write.csv(ttest.df,paste0(k,"vs",j,".ttest.csv"))
      #volcano plot
      plt_df<-ttest.df %>% tibble::rownames_to_column(var = "antigen") %>% 
        mutate(up_in = ifelse(p.adj >= thresh_p.adj, "NS", ifelse(log2FC > thresh_lfc, "group2", ifelse(log2FC < -thresh_lfc, "group1", "NS"))),
               rank_pval = rank(p.adj), rank_lfc_inc = rank(log2FC), rank_lfc_dec = rank(-abs(log2FC)),
               lp = -log10(p.adj)) %>% 
        arrange(-abs(log2FC))
      if (sum(plt_df$up_in != "NS") >= 30) {
        plt_df <- plt_df %>% mutate(antigen_label = ifelse(up_in != "NS" & (rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), antigen, NA))
      } else {
        plt_df <- plt_df %>% mutate(antigen_label = ifelse((up_in != "NS" | rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), antigen, NA))
      }
      plt_df$lp <- pmin(plt_df$lp, 10) # pval bound
      plt_df$log2FC <- pmin(plt_df$log2FC, 4) #lfc bound
      plt_df$log2FC <- pmax(plt_df$log2FC, -4) #lfc bound
      group2.color="red"
      if(j=="BMchim.564-Icos.RAG-TCRA"){group2.color="#70AD47"}
      if(j=="BMchim.564-Icos.RAG-TCRB"){group2.color="#7030A0"}
      if(j=="BMchim.564-Icos.RAG-TCRE"){group2.color="#FFC000"}
      if(j=="BMchim.564"){group2.color="#00BFC4"}
      group1.color="black"
      if(k=="564homo"){group1.color="black"}
      if(k=="564homo.Icos"){group1.color="#4472C4"}
      if(k=="BMchim.564"){group1.color="#00BFC4"}
      max_x<-max(abs(plt_df$log2FC))
      ggplot(plt_df, aes(x = log2FC, y = lp, color = up_in, label = antigen_label)) +
        geom_point(size = 0.5) +
        ggrepel::geom_text_repel(size = 2, segment.size = 0.1, seed = 1, show.legend=F) +
        scale_color_manual(values = c("group2" = group2.color, "group1" = group1.color, "NS" = "grey80")) +
        scale_x_continuous(limits = c(-max_x, max_x), expand = expansion(mult = c(0.01, 0.01))) +
        scale_y_continuous(limits = c(0, max(plt_df$lp)), expand = expansion(mult = c(0.05, 0.01))) +
        labs(title = paste0(k," vs ", j), color = "Up-regulated in:", y = NULL, x =NULL) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),legend.position="bottom",legend.title=element_text(size=8),
              plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
              axis.title = element_text(size = 10, color = "black"))
      ggsave2(paste0(k,"vs",j,".volcano.png"),width=3, height=3,device="png")
    }
  }
}

#export for prism
df<-df[order(df$tx2),]
select.ags<-c("Ag1","Ag2","Ag3")
select.ags<-antigens
write.csv(df[,c("tx2",select.ags)],"autoag.prism.csv")


##heatmap
mtx<-select(autoag,-c(mouse))
rownames(mtx)<-autoag$mouse
#annotations
annot.row<-autoag.meta
rownames(annot.row)<-annot.row$Ag
annot.row<-select(annot.row,c(Type,Class))
annot.col<-meta
rownames(annot.col)<-meta$mouse
annot.col<-select(annot.col,c(tx2))#Gender,age.round,
annot_colors<-list(tx2=c("564homo"="black","BMchim.564-Icos.RAG-TCRA"="#70AD47","BMchim.564-Icos.RAG-TCRB"="#7030A0",
                         "BMchim.564-Icos.RAG-TCRE"="#FFC000","BMchim.564"="#00BFC4","564homo.Icos"="#4472C4"),
                   Class=c("IgG"="#E7872B","IgG2a"="#5AAA46","IgG2c"="#317EC2"))
#plot
p <- pheatmap(t(mtx), 
              cluster_rows = TRUE,cluster_cols = TRUE,clustering_method = "ward.D2",
              scale="row",
              show_rownames = T, treeheight_row=5,treeheight_col=2,show_colnames=F,
              color = colorRampPalette(c("navy","white","firebrick3"))(100),
              fontsize = 8,fontsize_col = 8,fontsize_row=5,
              annotation_names_row=F, annotation_col=annot.col,annotation_names_col=F,
              annotation_colors=annot_colors,
              annotation_row=annot.row,
              silent = T)
png("autoag.heatmap.png", width = 15, height = 7, res = 200,units="in")
grid.arrange(p$gtable)
dev.off()


##PCA
df<-left_join(autoag,meta,by="mouse")
rownames(df)<-df$mouse
d1<-select(df,-c(mouse))
num_data<-select(d1,-c(Gender,age.round,tx2))
num_data<-num_data[ , which(apply(num_data, 2, var) != 0)] #getting rid of zero variance
mtcars.pca <- prcomp(num_data, center = TRUE,scale. = TRUE)
summary(mtcars.pca)

#PCA plots
autoplot(mtcars.pca,data=d1,colour='tx2')+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.title=element_blank())+
  scale_color_manual(values=c("564homo"="black","BMchim.564-Icos.RAG-TCRA"="#70AD47","BMchim.564-Icos.RAG-TCRB"="#7030A0",
                            "BMchim.564-Icos.RAG-TCRE"="#FFC000","BMchim.564"="#00BFC4","564homo.Icos"="#4472C4"))
ggsave2("autoag.PCA.condition.png",width=5.5, height=3.2,device="png")

#compute loadings and contributions
loadings <- mtcars.pca$rotation
sdev <- mtcars.pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
var.cos2 <- var.coord^2
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
var.contrib<-as.data.frame(var.contrib)
var.contrib$ag<-rownames(var.contrib)
for(i in c("PC1","PC2")){
  temp<-var.contrib[order(-var.contrib$PC1),c("ag",i)]
  temp<-head(temp,n=10)
  ggplot(temp,aes(x=reorder(rownames(temp),!!sym(i)), y=!!sym(i))) +
    geom_point( color="blue", size=4, alpha=0.6)+
    geom_segment( aes(x=rownames(temp), xend=rownames(temp), y=0, yend=!!sym(i))) +#,color='skyblue'
    labs(x="",y=i)+ theme(panel.border = element_rect(colour = "black", size=1,fill=NA), panel.background = element_blank())+
    coord_flip() 
  ggsave2(paste0(i,".png"),width=4, height=3,device="png")
}



