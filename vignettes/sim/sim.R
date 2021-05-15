#!/usr/bin/env Rscript
library("liger")
library("irlba")
library("Seurat")
library("gridExtra")
library("splatter")
library("scater")
library("ggplot2")
library("umap")
library("optparse")


addnoise <- function(dt=NULL, snr=NULL){
  for(i in seq(1, nrow(dt),1)){
    rowSd <- sd(dt[i,])
    tp <- matrix(rnorm(length(dt[i,])), nrow=1,ncol = length(dt[i,]))*rowSd*snr
    dt[i,] <- dt[i,] + tp
  }
  return(dt)
}

getPCA <- function(in1=NULL){
  in1 <- CreateSeuratObject(counts = in1)
  in1 <- ScaleData(in1,features = rownames(in1))
  in1 <- RunPCA(in1, features = rownames(in1))
  myumap <- umap(in1@reductions$pca@cell.embeddings)
  return(myumap$layout)
}


getGeneScore <- function(dt1 = NULL, dt2=NULL, c = NULL, gene_loci=NULL,win_size=3){
  
  
  gene_mat <- dt1
  for(i in seq(1,nrow(gene_mat),1)){
    # identify cis-window along the genome based on gene_loci info, using win_size=3
    start<-max(gene_loci[i]-win_size, 0)
    end <- min(gene_loci[i]+win_size,nrow(dt2))
    gene_activity <-colSums(dt2[seq(start,end,1),])
    norm1 <- sum(gene_mat[i,]^2)
    norm2 <- sum(gene_activity^2)
    if(abs(norm1*norm2)>0){
      gene_mat[i,] <- (1-c)*gene_mat[i,]/norm1 + c* gene_activity/norm2
    }
  }
  return(gene_mat)
}

sim_XYZ <- function(cell_perct=NULL, noise=NULL, nGene=500, nPeak=1500, nCell=2000){
  
  params<-newSplatParams()
  params<-setParam(params, "nGenes", nGene)
  params.groups <- newSplatParams(batchCells = nCell, nGenes = nGene)
  sim_tp <- splatSimulateGroups(params.groups,
                                group.prob = cell_perct,
                                verbose = FALSE)
  sim_tp <- logNormCounts(sim_tp)
  sim <- runPCA(sim_tp)
  sim <- as.Seurat(sim, counts = "counts", data = "logcounts")
  
  X0 <- sim[["RNA"]][]
  x.clst <- sim@meta.data$Group
  y.clst <- sim@meta.data$Group
  
  # generate trans-peak matrix 
  nrow <- dim(X0)[1]
  ncol <- dim(X0)[2]
  
  # randomly sampling unit to generate the peak profile 

  pos <- sample(nrow, nPeak-nrow, replace = TRUE)
  peak_mat  <- addnoise(X0[pos,], snr = 0)
  # add original gene expression matrix 
  Y <- rbind(peak_mat, X0)

  
  # insert the gene location in the peak matrix 
  gene_loci <- sort(sample(nPeak-nrow, nrow, replace=TRUE))
  

  # build gene expression matrix X by summing cis and trans peak signal
  X <- X0
  dt_x <-getPCA(in1=X)
  tp <- data.frame("UMAP1"=dt_x[,1],"UMAP2"=dt_x[,2],"celltype"=x.clst)
  
  # check umap from gene expression 
  p1<-ggplot(tp, aes(x = UMAP1, y = UMAP2, color=celltype)) + 
    geom_point(alpha = 0.5, size =0.5)  
  
  
  # peak matrix Y
  dt_x <-getPCA(in1=Y)
  tp <- data.frame("UMAP1"=dt_x[,1],"UMAP2"=dt_x[,2],"celltype"=x.clst)
  # check umap from peak profile 
  p3<-ggplot(tp, aes(x = UMAP1, y = UMAP2, color=celltype)) + 
    geom_point(alpha = 0.5, size =0.5)  
  
  # model gene activity Z0 by adding noise 
  Z0 <- getGeneScore(dt1 = X0,dt2 = Y , c = noise, gene_loci=gene_loci)
  print(cor(as.vector(X), as.vector(Z0)))
  rownames(X) <- paste0("gene",seq(1,nrow(X),1))
  colnames(X) <- paste0("a",seq(1,ncol(X),1))
  rownames(Y) <- paste0("peak",seq(1,nrow(Y),1))
  colnames(Y) <- paste0("b",seq(1,ncol(Y),1))
  rownames(Z0) <- paste0("gene",seq(1,nrow(Z0),1))
  colnames(Z0) <- paste0("b",seq(1,ncol(Z0),1))
  # check umap from gene activity Z0
  dt_x <-getPCA(in1=Z0)
  tp <- data.frame("UMAP1"=dt_x[,1],"UMAP2"=dt_x[,2],"celltype"=x.clst)
  p2<-ggplot(tp, aes(x = UMAP1, y = UMAP2, color=celltype)) + 
    geom_point(alpha = 0.5, size =0.5)  
  
  simData <- list()
  simData$X <- X
  simData$Y <- Y
  simData$Z0 <-  Z0
  simData$x.clst <- x.clst
  simData$y.clst <- y.clst
  simData$p1 <- p1
  simData$p2 <- p2
  simData$p3 <- p3
  return(simData)
}


option_list = list(
  make_option(c("--SNR"), type="double", default=0.1, 
              help="noise level assigned on gene activity matrix", metavar="character"),
  make_option(c("--nGene"), type="integer", default=500, 
              help="number of genes", metavar="character"),
  make_option(c("--nPeak"), type="integer", default=1500, 
              help="number of peaks", metavar="character"),
  make_option(c("--nCell"), type="integer", default=500, 
              help="number of cells", metavar="character"),
  make_option(c("--out"), type="character", default="out", 
              help="output directory name [default= %default]", metavar="character")
) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


#prob_vec <- c(0.3,0.3,0.3,0.2,0.2,0.2,0.1,0.1,0.1,0.05,0.05,0.05) # 12 populaiton with different sample size
prob_vec <- c(0.3,0.3,0.3)  # three populaitons with balanced sample size
prob_vec <- prob_vec/sum(prob_vec)

# The noise level is assigned on gene activity matrix and cis/trans effect matrix

#noise <- 0.1
#nGene <- 500
#nPeak <- 1500
#nCell <- 2000

print(opt)
noise <-opt$SNR
popSize <- length(prob_vec)
out_dir <- paste0("Ncell",opt$nCell,"_PopSize",opt$popSize,"_nGene",opt$nGene,"_nPeak",opt$nPeak,"_noise",noise)
simData <- sim_XYZ(cell_perct= prob_vec, noise=noise, nGene=opt$nGene, nPeak=opt$nPeak, nCell=opt$nCell)
#simData <- sim_XYZ(cell_perct = prob_vec, noise=noise)
saveRDS(simData,file=paste0(opt$out,"/" , out_dir,".RDS"))




