#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

cisSim<-function(dt,rate){
  result<-list()
  gene<-rownames(dt)
  seq<-seq(1,length(gene),1)
  rd<-matrix(0, length(gene),1)
  for(i in seq(1,length(gene),1)){
    if(runif(1)<rate){rd[i,1]<-1}
  }
  sel<-seq[rd[,1]>0]
  sel_permu<-sample(sel)
  interaction<-cbind(sel,sel_permu)
  colnames(interaction)<-c("gene1","gene2")
  interaction<-data.frame(gene1=sel, gene2=sel_permu)
  gene_order<-gene
  for(i in seq(1,length(gene),1)){
    if(rd[i,1]==0){gene_order[i]<-seq[i]}
    else{ gene_order[i]<-interaction$gene2[interaction$gene1==seq[i]]}
  }
  gene_order<-as.numeric(gene_order)
  tp<-dt[["RNA"]][gene_order]
  rownames(tp)<-rownames(dt)
  out<-CreateSeuratObject(counts=tp)
  out<-AddMetaData(out, dt@meta.data)
  result[[1]]<-out
  result[[2]]<-gene_order
  return(result)
}


library("liger")
library("Seurat")
library("gridExtra")
library("splatter")
library("scater")
library("ggplot2")
library("umap")

rate<-args[1]
method<-args[2]

print(rate)
print(method)

params<-newSplatParams()
params<-setParam(params, "nGenes",400)
params.groups <- newSplatParams(batchCells = 300, nGenes = 400)
sim3 <- splatSimulateGroups(params.groups,
                            group.prob = c(0.33, 0.33, 0.34),
                            verbose = FALSE)
sim3 <- logNormCounts(sim3)
sim <- runPCA(sim3)



sim <- as.Seurat(sim, counts = "counts", data = "logcounts")
rate<-0.9
sim_new<-cisSim(sim,rate)
sim_new<-sim_new[[1]]
A<-as.matrix(sim[["RNA"]][])
B<-as.matrix(sim_new[["RNA"]][])


sim<-c()
#sim$X <- A[1:800, seq(1,1000,2)]
#sim$Z0 <- B[1:800,1:600]
#sim$Y <- A[seq(1,1000,2),1:600]

sim$X <- A
sim$Z0 <- B
sim$Y <- A


sim$X <- sim$X[rowSums(sim$X)>0,]
sim$Y <- sim$Y[rowSums(sim$Y)>0,]
sim$Z0 <- sim$Z0[rowSums(sim$Z0)>0,]

f <- intersect(rownames(sim$X), rownames(sim$Z0))
sim$X <- sim$X[f,]
sim$Z0 <- sim$Z0[f,]

sim$X_meta <- sim_new@meta.data[colnames(sim$X),]
sim$Y_meta <- sim_new@meta.data[colnames(sim$Y),]
saveRDS(sim,file="sim.rds")




cell_type<-c(sim$X_meta$Group,sim$Y_meta$Group)
type<-c(rep("A",ncol(sim$X)), rep("B",ncol(sim$Y)))






if(method=="liger"){
  colnames(sim$X) <- paste(colnames(sim$X),"_ref",sep="")
  RNA <- CreateSeuratObject(counts = sim$X-min(sim$X) , project = "test", min.features =0, min.cells = 0)
  ATAC <- CreateSeuratObject(counts = sim$Z0 - min(sim$Z0), project = "test", min.features =0, min.cells = 0)

  ligerex = seuratToLiger(list(RNA, ATAC), names = c('rna', 'atac'))
  #ligerex = liger::normalize(ligerex)
  ligerex@scale.data$rna <- t(as.matrix(ligerex@raw.data$rna))
  ligerex@scale.data$atac <- t(as.matrix(ligerex@raw.data$atac))

  ligerex@var.genes<-intersect(rownames(ligerex@raw.data$rna), rownames(ligerex@raw.data$atac))
  #ligerex = scaleNotCenter(ligerex)
  ligerex = optimizeALS(ligerex, k = 6) 
  ligerex = quantileAlignSNF(ligerex,  dist.use="kd_tree")
  tmp <- umap(ligerex@H.norm)
  #ligerex <- liger::runUMAP(ligerex, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)
  #tmp<-ligerex@tsne.coords

  tmp<-data.frame(UMAP1=tmp$layout[,1],UMAP2=tmp$layout[,2],Datatype=type,cell_type=cell_type)
  write.csv(tmp,file=paste("liger",rate,".csv",sep=""))

   
}

if(method=="cca"){
  if(1){
      genes.use <- f
      atac<-CreateSeuratObject(counts=sim$Z0)
      atac<-AddMetaData(atac, sim$Y_meta)
      rna<-CreateSeuratObject(counts=sim$X)
      rna<-AddMetaData(rna, sim$X_meta)
      atac <- ScaleData(atac, features = genes.use, do.scale = FALSE)
      atac <- RunPCA(atac, features = genes.use, verbose = FALSE)

      transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = f, 
          reference.assay = "RNA", query.assay = "RNA", reduction = "cca")
      celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$Group, weight.reduction = atac[["pca"]], dims=1:5)
      atac <- AddMetaData(atac, metadata = celltype.predictions)

      refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

      imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["pca"]], dims=1:5)
      tmp <- atac
      atac[["RNA"]] <- imputation
      coembed <- merge(x =rna, y = atac)


      # Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
      # datasets
      coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
      coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
      coembed <- RunUMAP(coembed, dims = 1:30)
      tmp<-coembed@reductions$umap@cell.embeddings
  }


  if(0){
      cca.results <- RunCCA(
        object1 = X1,
        object2 = Y, 
        standardize = TRUE,
        num.cc = 50,
        verbose = verbose,
      )
      tmp <- umap(cca.results$ccv)
      tmp<-tmp$layout
  }

  tmp<-data.frame(UMAP1=tmp[,1], UMAP2=tmp[,2], cell_type=cell_type, Datatype=type)
  tmp$cell_type<-as.factor(tmp$cell_type)
  x_limit<-c(quantile(tmp$UMAP1,0.005), quantile(tmp$UMAP1,1-0.005))
  y_limit<-c(quantile(tmp$UMAP2,0.005), quantile(tmp$UMAP2,1-0.005))
  write.csv(tmp,file=paste("seurat",rate,".csv",sep=""))
   
}
