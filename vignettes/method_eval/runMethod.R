#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("harmony")
library("bindSC")
library("liger")
library("irlba")
library("Seurat")
library("gridExtra")
library("splatter")
library("scater")
library("ggplot2")
library("umap")



runHarmony <- function(X=NULLL,Z0 = NULL, K =NULL){
  
  memory_start <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14]) 
  time_start <- Sys.time()
  
  gene.use <- rownames(X)
  dt2<-CreateSeuratObject(counts=Z0)
  dt1<-CreateSeuratObject(counts=X)
  merged <- merge(x = dt1, y = dt2, add.cell.ids= c("dt1", "dt2"), project = "eval")
  
  merged@meta.data$type <- c(rep("X",dim(X)[2]), rep("Y",dim(Z0)[2]))
  
  merged <- NormalizeData(merged) %>% FindVariableFeatures() %>% ScaleData() 
  
  merged <-RunPCA(merged, verbose = FALSE, features = VariableFeatures(merged))
  
  merged <- RunHarmony(merged, group.by.vars = "type")

  res<- merged@reductions$harmony@cell.embeddings[,seq(1,K,1)]
  
  
  #memory_usage <- sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt
  time_end <- Sys.time()
  time_usage <-(time_end - time_start)
  #print(c(time_usage, memory_usage))
  
  return(res)
  
}


runLIGER <- function(X=NULLL,Z0 = NULL, K =NULL){
  
  memory_start <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14]) 
  time_start <- Sys.time()
  
  
  colnames(X) <- paste(colnames(X),"_ref",sep="")
  RNA <- CreateSeuratObject(counts = X-min(X) , project = "eval", min.features =0, min.cells = 0)
  ATAC <- CreateSeuratObject(counts = Z0- min(Z0), project = "eval", min.features =0, min.cells = 0)
  
  ligerex = seuratToLiger(list(RNA, ATAC), names = c('dt1', 'dt2'))
  #ligerex = liger::normalize(ligerex)
  ligerex@scale.data$dt1 <- t(as.matrix(ligerex@raw.data$dt1))
  ligerex@scale.data$dt2 <- t(as.matrix(ligerex@raw.data$dt2))
  
  ligerex@var.genes<-intersect(rownames(ligerex@raw.data$dt1), rownames(ligerex@raw.data$dt2))
  #ligerex = scaleNotCenter(ligerex)
  ligerex = optimizeALS(ligerex, k = K) 
  ligerex = quantileAlignSNF(ligerex,  dist.use="kd_tree")
  res <- ligerex@H.norm
  
  
  
 # memory_usage <- sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt
  time_end <- Sys.time()
  time_usage <-(time_end - time_start)
  #print(c(time_usage, memory_usage))
  
  return(res)

}

runBindSC <- function(X=NULL, Y=NULL, Z0=NULL, 
            X.clst=NULL, Y.clst=NULL, K=NULL, alpha=NULL, lambda = NULL){
           
            memory_start <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14]) 
            time_start <- Sys.time()
            
             res <- BiCCA( X = X ,
             Y = Y,
             Z0 =Z0, 
             X.clst = X.clst,
             Y.clst = Y.clst,
             alpha = alpha, 
             lambda = lambda,
             K = K,
             temp.path  = "out",
             num.iteration = 50,
             tolerance = 0.01,
             save = TRUE,
             parameter.optimize = FALSE,
             block.size = 1000)
      
             #memory_usage <- sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt
             time_end <- Sys.time()
             time_usage <-(time_end - time_start)
             #print(c(time_usage, memory_usage))
             return(rbind(res$u, res$r))

}



runSeurat <- function(X=NULL, Z0=NULL, K = NULL){
  
  memory_start <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14]) 
  time_start <- Sys.time()
  
  dt1<-CreateSeuratObject(counts=X)
  dt2<-CreateSeuratObject(counts=Z0)
  gene.use <- rownames(X)
  dt2 <- ScaleData(dt2, features = gene.use, do.scale = TRUE)
  dt2 <- RunPCA(dt2, features = gene.use, verbose = FALSE)
  
  transfer.anchors <- FindTransferAnchors(reference = dt1, query = dt2, features = gene.use,  dims = seq(1,K,1), 
                                          reference.assay = "RNA", query.assay = "RNA", reduction = "cca")
  
  
  refdata <- GetAssayData(dt1, assay = "RNA", slot = "data")[gene.use, ]
  
  imputation <- TransferData(anchorset = transfer.anchors, 
                             refdata = refdata, weight.reduction = dt2[["pca"]], dims=seq(1,K,1))
  tmp <- dt2
  dt2[["RNA"]] <- imputation
  coembed <- merge(x =dt1, y = dt2)
  
  
  # Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
  # datasets
  coembed <- ScaleData(coembed, features = gene.use, do.scale = FALSE)
  coembed <- RunPCA(coembed, features = gene.use, verbose = FALSE)
 
  #memory_usage <- sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt
  time_end <- Sys.time()
  time_usage <-(time_end - time_start)
  #print(c(time_usage, memory_usage))
  res <- coembed@reductions$pca@cell.embeddings[,seq(1,K,1)]
  return(res)
  
}


