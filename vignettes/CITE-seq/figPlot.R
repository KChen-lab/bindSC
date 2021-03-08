

myplot<-function(dt=NULL,title=NULL,xmin=0,xmax=0, ymin=0, ymax=0){
  

    tp <- dt[,c("celltype.l1","UMAP1","UMAP2")]
    colnames(tp)[1]<-"celltype"
    tp$cluster <- tp$celltype
    label_pos<-aggregate(. ~ cluster, tp[,c("cluster","UMAP1","UMAP2")], median) 
    
    tx <- title
    if(xmin==0){
      xmin<-min(tp$UMAP1)
      xmax<-max(tp$UMAP1)
      ymin<-min(tp$UMAP2)
      ymax<-max(tp$UMAP2)
    }
    p <- ggplot(tp, aes(x = UMAP1, y = UMAP2,  color = cluster)) + 
      xlim(c(xmin,xmax)) + ylim(c(ymin,ymax))+
      geom_point(alpha = 0.5, size =0.25)   + ggtitle(tx) + 
      theme(plot.title = element_text(size = 40, face = "bold")) + 
      scale_colour_manual(values = paletteDiscrete(tp$cluster)) + 
      geom_label_repel(data = label_pos, 
                       aes(label = cluster),
                       color = "black",
                       point.padding =0.2,
                       box.padding = 0.2)  + theme_void() + 
      theme(legend.position = "none") + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", size=1))
    return(p)

 
}
