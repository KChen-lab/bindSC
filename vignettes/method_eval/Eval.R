


anchoring_dis <-function(dt1,dt2){
  dt1 <- as.matrix(dt1)
  dt2 <- as.matrix(dt2)
  d <- dt1[,1]
  X2 <- dt2
  for(i in seq(1,nrow(dt1),1)){
    X1 <- t(as.matrix(dt1[i,]))
    
    X3 <- t(as.matrix(dt2[i,]))
    d_11 <- cdist(X1,X3)
    d_tol <- cdist(X1, X2)
    d[i] <- d_11/max(as.vector(d_tol))*2
  }
  return(d)
}


calc_silhouette_coef<-function(x = NULL, clst = NULL){
  #dt$cluster<-as.character(dt$cluster)
  clst <- as.numeric(as.factor(clst))
  out<-silhouette(clst, dist(x))
  return(out)
}

calc_alignment_score<-function(x = NULL,  clst = NULL, k = NULL){
  
  #colnames(x)[n] <- "cluster"
  out<-get.knn(x, k= k)
  tp<-out$nn.index
  alignment<-rep(0,dim(x)[1])
  for(i in seq(1,dim(tp)[1],1)){
    flag<- clst[tp[i,]]
    alignment[i]<-Entropy(table(flag), base=exp(1))
  }
  return(alignment)
}

get_UMAP <- function(dt = NULL){
  a <- umap(dt)
  return(a$layout)
}

score_plot <- function(plt_dt = NULL){
  method_color <- c("red","lightblue3","slateblue1","lightsalmon","palegreen4")
  plt_dt$score<-factor(plt_dt$score, levels =c("Silhouette score",
    "Alignment mixing score","Anchoring distance"))
  p<-ggboxplot(plt_dt, x = "method", y = "value", fill = "method",
               outlier = FALSE, order =c("bindSC","Seurat","LIGER","Harmony"), outlier.shape=NA,
               palette =method_color[-2])+ facet_wrap(~score, scales = "free") + 
    theme(axis.text.x=element_text(angle = -60, hjust = 0))+ 
    theme(legend.position = "none") + xlab("")
  return(p)
}

alignmentScore <-function(coembed = NULL, merged_dt = NULL){
  
  
  method_lst <- c("bindSC","Seurat","LIGER","Harmony")
  align_score <- c()
  for(method in method_lst){
    pos <- merged_dt$method==method
    plt_dt <- coembed[pos,]
    clst <- merged_dt$celltype[pos]
    tech <- merged_dt$tech[pos]
    silhouette <- calc_silhouette_coef(x=plt_dt, clst = clst)
    align_mix <- calc_alignment_score(x=plt_dt, clst = tech, k = 10)
    anchor_dis <- anchoring_dis(dt1 = plt_dt[tech=="A",], dt2 = plt_dt[tech=="B",])
    tp <- data.frame("score"=c(rep("Silhouette score", length(silhouette[,3])),
                               rep("Alignment mixing score", length(align_mix)),
                               rep("Anchoring distance", length(anchor_dis))),
                     "value" = c(silhouette[,3], align_mix, anchor_dis))
    tp$method <- method
    align_score <- rbind(align_score, tp)
  }
  return(align_score)
  
  
}


umapPlot <- function(merged_dt = NULL, tech_label=NULL){
  
  merged_dt$method <- factor(merged_dt$method, levels=c("bindSC","Seurat","LIGER","Harmony"))
  
  i<-0 
  j<-0
  fig <- list()
  #tech_label<-c("scRNA-seq","scATAC-seq")
  method_lst <- c("bindSC","Seurat","LIGER","Harmony")
  for(method in method_lst){
    plt_dt <- merged_dt[merged_dt$method==i,]
    j <- j + 1
    i<-0
    for(tech in c("A","B")){
      i <- i + 1
      tp <- merged_dt[merged_dt$method==method,]
      plt_dt_method <- merged_dt[merged_dt$method==method & merged_dt$tech==tech,]
      label_pos<-aggregate(. ~ celltype, plt_dt_method, median)
      #label_pos_method <- label_pos[label_pos$method==method & label_pos$tech==tech,]
      
      p<-ggplot(plt_dt_method, aes(x = UMAP1, y = UMAP2, color=celltype)) + 
        geom_point(alpha = 0.5, size =0.5)  + 
        scale_colour_manual(values = paletteDiscrete(plt_dt_method$celltype))  +
        geom_text_repel(data = label_pos, aes(label=celltype),color="black",point.padding = 0.2, box.padding =0.2)+
        theme_classic() + xlim(c(min(tp$UMAP1), max(tp$UMAP1))) + 
        ylim(c(min(tp$UMAP2), max(tp$UMAP2))) 
      
      
      if(i==1){p <- p + ggtitle(method_lst[j])}
      if(j==1){ 
        p <- p +  ylab(tech_label[i])+
          theme(axis.text=element_blank(),
                axis.title.x=element_blank(),
                axis.ticks=element_blank()) 
      }
      else{
        p <- p + ylab("")+
          theme(axis.title.x=element_blank(),
                axis.text=element_blank(),
                axis.ticks=element_blank())
        
      }
      p <- p + theme(legend.position = "none")  
      fig[[j+i*4-4]] <- p
    }
  }
  
  library(gridExtra)
  library(ggpubr)
  
  p <- ggarrange(plotlist = fig, ncol=4, nrow=2)
  return(p)
}

















