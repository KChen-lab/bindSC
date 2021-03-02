
library(umap)
library(stringr)
library(cluster)
library(FNN)
library(ggpubr)
library(DescTools)


calc_silhouette_coef<-function(dt){
  if(dim(dt)[1]<=600){


  }else{
    dt$Time<-as.character(dt$Time)
    dt$Time[dt$Time=="A549_0h"]<-0
    dt$Time[dt$Time=="A549_1h"]<-1
    dt$Time[dt$Time=="A549_3h"]<-1
    dt$Time[dt$Time==3]<-1
   
  }
  dt$Time<-as.numeric(dt$Time)
  out<-silhouette(dt$Time, dist(cbind(dt$UMAP1,dt$UMAP2)))
  return(out)
}

calc_alignment_score<-function(dt){
  
  out<-get.knn(cbind(dt$UMAP1, dt$UMAP2), k=20)
  tp<-out$nn.index
  alignment<-rep(0,dim(dt)[1])
  for(i in seq(1,dim(tp)[1],1)){
    flag<-dt$Technology[tp[i,]]
    alignment[i]<-Entropy(table(flag), base=exp(1))
  }
  return(alignment)
}

find_cellmatch<-function(dt,lst,size){
  out<-get.knn(cbind(dt$UMAP1, dt$UMAP2), k=size)
  tp<-out$nn.index
  find<-0
  tol<-0
  for(i in seq(1,dim(tp)[1],1)){
    if(lst[i]>0){
      tol<-tol+1
      if(lst[i]%in%tp[i,]){
        find<-find+1
      }
    }
  }
  out<-find/tol
  return(out)
}

  
# Get anchors from goldstandard


umap_plot <- function(out, cell_type, data_type){

  #tmp <- umap(rbind(out$u, out$r))
  tmp <- out
  plt_dt <-data.frame("UMAP1"=tmp$layout[,1],"UMAP2"=tmp$layout[,2], 
                   "cell_type"=cell_type,  
                   "datatype"=data_type)

  colnames(plt_dt)<-c("UMAP1","UMAP2","Time","Technology")
  plt_dt$Time<-factor(plt_dt$Time)

  p1<-ggscatter(plt_dt, x = "UMAP1", y = "UMAP2",
     color = "Time", palette = c("darkseagreen4","lightpink","darkorchid1"), 
     repel = FALSE,size=0.5, alpha=0.5,legend.title = "", title="bindSC", font.title=16) 

  p2<-ggscatter(plt_dt, x = "UMAP1", y = "UMAP2",
     color = "Technology", palette = c("indianred3", "#E7B800"), repel = FALSE,size=0.5, alpha=0.5,legend.title = "",title="", font.title=16)

  p <- ggarrange(p1,p2,nrow = 1, common.legend = FALSE)
  out <- c()
  out$plt_dt <- plt_dt
  out$plt <- p 
  return(out)
}



method_compare <- function(plt_dt, dataset){

    bindSC <- plt_dt

    seurat<-read.csv(paste("../../inst/",dataset,"/seurat.csv",sep=""))
    liger<-read.csv(paste("../../inst/",dataset,"/liger.csv",sep=""))

    #seurat <- read.csv("/home/jdou1/project/integraModel/dev_6/A549/ins/A549/seurat.csv")


    #liger <- read.csv("/home/jdou1/project/integraModel/dev_6/A549/ins/A549/liger.csv")


    bindSC<-cbind(seurat$id,bindSC)
    #colnames(bindSC)<-c("id","UMAP1","UMAP2",  "Time","Technology")
    seurat$Technology<-bindSC$Technology
    seurat$Time<-bindSC$Time
    #colnames(seurat)<-c("id","UMAP1","UMAP2","Technology","Time")
    liger$Technology<-bindSC$Technology
    liger$Time<-bindSC$Time
    #colnames(liger)<-c("id","UMAP1","UMAP2","Technology","Time")

    p1<-ggscatter(bindSC, x = "UMAP1", y = "UMAP2",
       color = "Time", palette = c("darkseagreen4","lightpink","darkorchid1"), 
       repel = FALSE,size=0.5, alpha=0.5,legend.title = "", title="bindSC", font.title=16) 

    p2<-ggscatter(bindSC, x = "UMAP1", y = "UMAP2",
       color = "Technology", palette = c("indianred3", "#E7B800"), repel = FALSE,size=0.5, alpha=0.5,legend.title = "",title="", font.title=16)

    p3<-ggscatter(seurat, x = "UMAP1", y = "UMAP2",
       color = "Time", palette = c("darkseagreen4","lightpink","darkorchid1"), 
       repel = FALSE,size=0.5, alpha=0.5,legend.title = "", title="Seurat", font.title=16) 
    p4<-ggscatter(seurat, x = "UMAP1", y = "UMAP2",
       color = "Technology", palette = c("indianred3", "#E7B800"), repel = FALSE,size=0.5, alpha=0.5,legend.title = "")

    p5<-ggscatter(liger, x = "UMAP1", y = "UMAP2",
       color = "Time", palette = c("darkseagreen4","lightpink","darkorchid1"), 
       repel = FALSE,size=0.5, alpha=0.5,legend.title = "", title="Liger") 
    p6<-ggscatter(liger, x = "UMAP1", y = "UMAP2",
       color = "Technology", palette = c("indianred3", "#E7B800"), repel = FALSE,size=0.5, alpha=0.5,legend.title = "")  

 
    p1<-p1+ rremove("ticks") + rremove("xylab")+ rremove("axis.text") 
    p2<-p2+ rremove("ticks") + rremove("xylab")+ rremove("axis.text")
    p3<-p3+ rremove("ticks") + rremove("xylab")+ rremove("axis.text")
    p4<-p4+ rremove("ticks") + rremove("xylab")+ rremove("axis.text")
    p5<-p5+ rremove("ticks") + rremove("xylab")+ rremove("axis.text")
    p6<-p6+ rremove("ticks") + rremove("xylab")+ rremove("axis.text")

    a<-ggarrange(p1, p3, p5, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
    b<-ggarrange(p2, p4, p6, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
    f1 <- ggarrange(a,b,nrow = 2, common.legend = FALSE) 
    #tiff("Coembedding.tiff", width=8, height =6, units = 'in', res=300, compression="lzw")
    #print(figure)
    #dev.off()


    ID<-str_replace(seurat$id,"_1","")
    ID<-str_replace(ID,"_2","")
    rd<-rep(0,length(ID))
    for(i in seq(1,length(ID),1)){
      lst<-which(ID==ID[i],arr.ind = TRUE)
      if(length(lst)>1){
        rd[i]<-lst[1]
        if(lst[1]==i){rd[i]<-lst[2]}
      }
    }

    bindSC_silhoutte<-calc_silhouette_coef(bindSC)
    seurat_silhoutte<-calc_silhouette_coef(seurat)
    liger_silhoutte<-calc_silhouette_coef(liger)
    dt_plt<-data.frame("value"=c(bindSC_silhoutte[,3],seurat_silhoutte[,3], liger_silhoutte[,3]), 
        "method"=c(rep("bindSC",length(bindSC_silhoutte[,3])),
                   rep("Seurat",length(seurat_silhoutte[,3])),
                   rep("Liger",length(liger_silhoutte[,3]))))

    library(ggpubr)
    my_comparisons<-list(c("bindSC","Liger"),c("Liger","Seurat"))
    f2 <- ggboxplot(dt_plt, x="method",y="value", fill="method", color="gray",
                  palette = c("firebrick3", "#FFA07A", "palegreen"),
                  xlab="", ylab="Silhoutte coefficient", width=0.5,outlier.shape=NA) + 
      stat_compare_means(comparisons = my_comparisons, label.y=c(0.75,0.65)) +  ylim(-0.3,0.8)

    #tiff("SilhoutteCor.tiff", width=4, height =5, units = 'in', res=300, compression="lzw")
    #print(p4)
    #dev.off()

    # alignment score 

    bindSC_alignment<-calc_alignment_score(bindSC)
    seurat_alignment<-calc_alignment_score(seurat)
    liger_alignment<-calc_alignment_score(liger)
    dt_plt<-data.frame("value"=c(bindSC_alignment,seurat_alignment, liger_alignment), 
        "method"=c(rep("bindSC",length(bindSC_alignment)),
                   rep("Seurat",length(seurat_alignment)),
                   rep("Liger",length(liger_alignment))))
    library(ggpubr)
    my_comparisons<-list(c("bindSC","Liger"),c("Liger","Seurat"))
    f3 <- ggboxplot(dt_plt, x="method",y="value", fill="method", color="gray",
                  palette = c("firebrick3", "#FFA07A", "palegreen"),
                  xlab="", ylab="Alignment score", width=0.5,outlier.shape=NA) + 
      stat_compare_means(comparisons = my_comparisons) 

    #tiff("AlignmentScore.tiff", width=4, height =5, units = 'in', res=300, compression="lzw")
    #print(p5)
    #dev.off()

    # cell_anchor accuracy 
    dt_plt<-c()
    if(dataset=="A549"){step <- c(5,25,100,200,500)}
    if(dataset=="Sim"){step <- c(1,5,10,15,20)}

    for(i in step){
      out1<-find_cellmatch(bindSC,rd,i)
      out2<-find_cellmatch(liger,rd,i)
      out3<-find_cellmatch(seurat,rd,i)
      dt_plt<-rbind(dt_plt, c(out1,i,"bindSC"))
      dt_plt<-rbind(dt_plt, c(out2,i,"Seurat"))
      dt_plt<-rbind(dt_plt, c(out3,i,"Liger"))
      dt_plt<-rbind(dt_plt, c(i/length(rd),i,"Random"))
    }
    colnames(dt_plt)<-c("anchor_accuracy","neighbor_size","Method")
    dt_plt<-as.data.frame(dt_plt)
    dt_plt$anchor_accuracy<-as.numeric(as.character(dt_plt$anchor_accuracy))
    dt_plt$neighbor_size<-factor(dt_plt$neighbor_size,levels=step)
    dt_plt$Method<-factor(dt_plt$Method,levels=c("bindSC","Liger","Seurat","Random"))
    f4 <- ggline(dt_plt, "neighbor_size", "anchor_accuracy", color = "Method",
     palette = c("firebrick3", "#FFA07A", "palegreen","grey")) + xlab("Neighbor size") + 
      ylab("Anchor accuracy")
    #tiff("Anchor_accuracy.tiff", width=5, height =5, units = 'in', res=300, compression="lzw")
    #print(p6)
    #dev.off()  
    out<-c()
    out$coembedding <- f1
    f2 <- ggarrange(f2,f3, f4, nrow = 1, common.legend = TRUE)
    out$eval <- f2
    return(out)
}