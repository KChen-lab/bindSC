
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' L2 normalize the columns (or rows) of a given matrix
#' @param mat Matrix to cosine normalize
#' @param MARGIN Perform normalization over rows (1) or columns (2)
#'
#'
#' @return returns l2-normalized matrix
#'
#'



L2Norm <- function(mat, MARGIN = 1){
  normalized <- sweep(
    x = mat,
    MARGIN = MARGIN,
    STATS = apply(
      X = mat,
      MARGIN = MARGIN,
      FUN = function(x){
        sqrt(x = sum(x ^ 2))
      }
    ),
    FUN = "/"
  )
  normalized[!is.finite(x = normalized)] <- 0
  return(normalized)
}


#' Normalize the columns of a given matrix
#' @param mat Matrix to normalize
#'
#' @return returns column-normalized matrix
#'
#'

colScale <- function(X){

    X <- scale(X, TRUE, apply(X,2,sd))
    return(X)
}


#' Normalize the rows of a given matrix
#' @param mat Matrix to normalize
#'
#' @return returns row-normalized matrix
#'
#'

rowScale <- function(X){
    tp <- t(X) 
    tp <- scale(tp, TRUE, apply(tp,2,sd))
    X  <- t(tp)
    return(X)
}


FeatureScale <- function(X){
  for(i in seq(1, dim(X)[1],1)){
    X[i,] <- X[i,] - mean(X[i,])
    X[i,] <- X[i,]/sd(X[i,])
  }
  return(X)
}
#' Matrix Crossproduct
#' Given a list of matrices x1, x2, ... , xn as arguments, return a matrix cross-product. 
#' x1%*%x2%*%xn
#' @param lst list of matrix
#'
#' @return returns crossproduct
#'
#'

is.sparseMatrix <- function(x) is(x, 'sparseMatrix') 


FindVariableFeature<-function(x=null, label=null,  top_n = NULL, T = NULL ){

  # We use seurat function to DE analysis 
  dt <- CreateSeuratObject(counts = x)
  meta <- data.frame("cluster"=label)
  rownames(meta) <- colnames(x)
  dt <- Seurat::AddMetaData(dt, metadata = meta)
  #dt <- NormalizeData(dt, normalization.method = "LogNormalize", scale.factor = 10000)
  markers <- FindAllMarkers(dt, only.pos = TRUE, min.pct = 0.5 , logfc.threshold =T, test.use="t")
}



PseudoHeatmapPlot <- function(mat1 = NULL, mat2 = NULL, mat3 =NULL,add_anno = NULL, levels=NULL, order=NULL){


    overlap <- intersect(colnames(mat1), colnames(mat2))
    overlap <- sort(overlap)
    mat1 <- mat1[,overlap]
    mat2 <- mat2[,overlap]
    mat3 <- mat3[,overlap]

    col_group <- colsplit(overlap,"_",c("id1","id2"))[,1]
    mat1 <- FeatureScale(mat1+1-min(mat1))
    mat2 <- FeatureScale(mat2+1-min(mat2))
    mat3 <- FeatureScale(mat3)
  
    mat <- rbind(mat1, mat2, mat3)
   
    colData <- data.frame("cluster"=factor(col_group, levels=order))
    new_order <- order(colData$cluster)
    colData$cluster <- colData$cluster[order(colData$cluster)]
    #print(order(colData$cluster))
 
    row_anno <- factor(c(rep("Gene expression", nrow(mat1)),
                                        rep("Chromatin accessibility", nrow(mat2)),
                                        rep("TF activity", nrow(mat3))),
                       levels=c("Gene expression", "Peak", "TF activity"))
 

    label_feature <- intersect(add_anno, rownames(mat))
    pos <- match(label_feature, rownames(mat))
    print(mat[1:5,1:5])
    print(dim(mat))
    print(pos)
    h <- ArchRHeatmap(mat = as.matrix(mat[,new_order]) ,colData = colData, draw=FALSE, name="Level",
                    clusterCols = FALSE, clusterRows = TRUE, useRaster = TRUE,
                    customRowLabel=pos,
                    limits = c(-1,1),
                    split = row_anno, 
                    colAnnoPerRow = 30, scale=TRUE)
    return(h)
   
}


PseudoCellProfile <- function(x = NULL, resolution1 = 0.25, resolution2 = 1){
  x<- PseudoCell(x, resolution1 = resolution1, resolution2 = resolution2)
  n_X <- dim(x$X)[2]
  n_Y <- dim(x$Y)[2]
  pseudoCell_X  <- PseudoProfiling(x = x$X, id = x$pseudo_cell[seq(1,n_X,1)])
  pseudoCell_Y  <- PseudoProfiling(x = x$Y, id = x$pseudo_cell[seq(n_X+1,n_X+n_Y,1)])
  x$pseudoCell_X  <- pseudoCell_X
  x$pseudoCell_Y  <- pseudoCell_Y
  return(x)
}

PseudoProfiling <- function(x = NULL, id = NULL, bin=1000){
  res<-c()
  n_row <- nrow(x)
  for(i in seq(1,n_row,bin)){

    end <- min(i+bin-1, n_row)
    tp <- as.matrix(x[seq(i,end,1), ])
    mat<-as.data.frame(t(tp))
    mat$id <- id
    mat <- as.data.frame(mat)
    mat_median<-aggregate(. ~ id,mat, mean)
    if(i==1){ res<-mat_median}
    else{
      res<-cbind(res, mat_median[,-1])
    }
  }
  rownames(res) <- res$id
  res<-t(res[,-1])

  return(res)
}



PseudoCell <- function( x = NULL, resolution1 = NULL, resolution2 =NULL){
  coembedding <- L2Norm(as.matrix(rbind(x$u, x$r)))
  joint_cluster  <- knn_cluster(coembedding, k=20, resolution  = resolution1)
  
  # Get pseudo_cell for each joint cluster 
  pseudo_id <- rep("0",dim(coembedding)[1])

  impuClst <- unique(joint_cluster)
  for(i in impuClst){
    pos <- which(joint_cluster==i,  arr.ind=TRUE)
    if(length(pos)>20){
      tp <- knn_cluster(coembedding[pos,], k=5, resolution = resolution2)
      pseudo_id[pos] <- paste(i,"_",tp,sep="")
    }
    else{
      pseudo_id[pos] <- paste(i,"_",0,sep="")
    }    
  }
  x$impuClst <- joint_cluster
  x$pseudo_cell <- pseudo_id
  return(x)
}


RunModularityClustering <- function(
  SNN = matrix(),
  modularity = 1,
  resolution = 0.8,
  algorithm = 1,
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  print.output = TRUE,
  temp.file.location = NULL,
  edge.file.name = NULL
) {
  #edge_file <- edge.file.name %||% ''
  edge_file <- ''
  clusters <- RunModularityClusteringCpp(
    SNN,
    modularity,
    resolution,
    algorithm,
    n.start,
    n.iter,
    random.seed,
    print.output,
    edge_file
  )
  return(clusters)
}


RunModularityClusteringCpp <- function(SNN, modularityFunction, resolution, algorithm, nRandomStarts, nIterations, randomSeed, printOutput, edgefilename) {
    .Call('_Seurat_RunModularityClusteringCpp', PACKAGE = 'Seurat', SNN, modularityFunction, resolution, algorithm, 
      nRandomStarts, nIterations, randomSeed, printOutput, edgefilename)
}


knn_cluster <- function(x, k=20, prune.SNN =1/15, resolution = 0.25){

  l2norm <- L2Norm(x)
  snn <- swne::CalcSNN(t(l2norm), k = k, prune.SNN = prune.SNN)
  #snn<- as.matrix(snn)
  #snn <- as(snn, "sparseMatrix")
  #ids <- runLeiden(SNN=snn, resolution = resolution)

  ids <- RunModularityClustering(SNN = snn, resolution = resolution) 
  return(ids)
}



paraSel_plot <- function(dt = NULL){
  dt <- as.data.frame(dt)
  p1 <- ggplot( dt, aes(x=alpha, y = lambda, fill=alignment )) + scale_fill_gradient(low="white", high="blue") + geom_tile() + theme_classic()
  p2 <- ggplot( dt, aes(x=alpha, y = lambda, fill=silhouette1)) + scale_fill_gradient(low="white", high="blue") + geom_tile() + theme_classic()
  p3 <- ggplot( dt, aes(x=alpha, y = lambda, fill=silhouette2)) + scale_fill_gradient(low="white", high="blue") + geom_tile() + theme_classic()
  #p4 <- ggplot( dt, aes(x=alpha, y = lambda, fill=cost_all)) + scale_fill_gradient(low="white", high="blue") + geom_tile() + theme_classic()
  p <- ggarrange(p2, p3, p1, ncol=3, nrow=1)
  return(p)
}

score_plot <- function(dt = NULL){
  dt<-data.frame(dt)
  #dt$silhouette1 <- (dt$silhouette1 - min(dt$silhouette1))/(max(dt$silhouette1)-min(dt$silhouette1))
  #dt$silhouette2 <- (dt$silhouette2 - min(dt$silhouette2))/(max(dt$silhouette2)-min(dt$silhouette2))
  #dt$alignment <- (dt$alignment - min(dt$alignment))/(max(dt$alignment)-min(dt$alignment))
  p <- ggplot() + xlab("alpha") + ylab("score") + 
  geom_point(data= dt, aes(x=alpha,y = silhouette1, colour = "darkred")) + 
  geom_point(data=dt, aes(x=alpha,y = silhouette2, colour="steelblue")) + 
  geom_point(data=dt, aes(x=alpha,y = alignment, colour="green")) + 
  geom_line(data= dt, aes(x=alpha,y = silhouette1, colour = "darkred")) + 
  geom_line(data=dt, aes(x=alpha,y = silhouette2, colour="steelblue")) + 
  geom_line(data=dt, aes(x=alpha,y = alignment, colour="green")) + 
  scale_color_discrete(name = "score", 
    labels = c("silhouette (A)","alignment", "silhouette (B)")) + theme_bw() 
   
  return(p)
}



runUMAP <- function(x=NULL){
  tmp <- umap(x)
  return(tmp$layout)
}


PC_optimize <- function(x=NULL){

  x  <- x/dim(x)[2]
  y  <- crossprod_e(x,x)
  dim_max <- min(c(dim(x)[1], dim(x)[2]))
  dim_max <- min(c(500, dim_max-1))
  SVD <- irlba(y, nv = dim_max)
  p_value <- rep(1, dim(SVD$u)[2])
  FLAG <- TRUE
  k <- 0
  for(i in seq(1, dim(SVD$u)[2],1)){
    p_value[i] <- ad.test(SVD$u[,i])$p.value
    if(p_value[i]>0.1 && FLAG){
      k <- i
      FLAG = FALSE
    }
  }

  # identify dimension used
  dim <- (as.integer(k/5) + 1)*5
  out <- c()
  out$p_value  <- p_value[1:dim]
  out$k <- dim
  out$u <- SVD$u[,1:dim]
  plt <- data.frame("PC"=1:5, "pvalue"=0-log10(out$p_value))
  out$plot <- ggscatter(plt, x="PC", y="pvalue", ylab = "-log10(pvalue)",
        size=2, alpha=1, font.title=16) + 
     geom_hline(yintercept = 1, linetype = 2,color="red")
  return(out)
}




UMAP_plot <- function(meta = NULL, color=NULL, xlim=NULL, alpha = 0.1,
                      ylim=NULL, showLabel=TRUE, title=NULL, mylabel=NULL){

  library(ggrepel)
  #meta <- x@meta
  dt <- data.frame("UMAP1" = meta[,"UMAP1"], "UMAP2" = meta[,"UMAP2"], label=meta[,color])
  dt$label <- as.factor(dt$label)
  
  x_min <- xlim[1]
  x_max <- xlim[2]
  y_min <- ylim[1]
  y_max <- ylim[2]

  label_pos<-aggregate(. ~ label, dt, median)
  p1 <- ggplot(dt, aes(x = UMAP1, y = UMAP2,  color = label)) + 
    geom_point(alpha = alpha, size =0.5)   + 
    scale_colour_manual(values = mylabel) + 
    xlim(x_min, x_max) + 
    ylim(y_min, y_max) + 
    theme_classic()  +  theme(legend.position = "none") + 
    geom_text_repel(data = label_pos, repel = TRUE,
                    aes(label = label), color="black", fontface="bold",
                    alpha = 0.75,box.padding = 0.5, point.padding = 0.1) + 
    NoLegend() + theme(axis.text=element_blank(), axis.title=element_blank(),
                       axis.ticks=element_blank()) 
  
  
    
  return(p1)

}






plotIteration <- function(x = NULL){

    dt<-data.frame("iter"=seq(1, length(x$cost_l),1),
                "cost_left" = x$cost_l, 
                "cost_right" = x$cost_r,
                "cost_z0" = x$cost_z0, 
                "cost_all" = x$cost_all,
                "delta" = x$delta
                )

    p1 <- ggscatter(dt, x = "iter", y = "cost_left", conf.int = FALSE, color = "#00AFBB") + 
        xlab("Iteration time") + ylab("obj. on left CCA")

    p2 <- ggscatter(dt, x = "iter", y = "cost_right",, conf.int = FALSE, color = "#00AFBB") + 
        xlab("Iteration time") + ylab("obj. on right CCA")

    p3 <- ggscatter(dt, x = "iter", y = "cost_z0", conf.int = FALSE, color = "#00AFBB") + 
        xlab("Iteration time") + ylab("obj. on ||Z-Z0||")

    p4 <- ggscatter(dt, x = "iter", y = "cost_right",, conf.int = FALSE, color = "#00AFBB") + 
        xlab("Iteration time") + ylab("obj. on all three terms")

    p5 <- ggscatter(dt, x = "iter", y = "delta", conf.int = FALSE, color = "#00AFBB") + 
        xlab("Iteration time") + ylab("Relative change of Z")


    p <- ggarrange(p1, p2, p3, p4, p5 ,nrow = 2, ncol=3) 

    return(p)
} 




celltype_assign <- function(train_x=NULL, train_y = NULL, test = NULL, test_clst = NULL){

  res <- list()
  clst_num <- length(unique(train_y))
  n  <- Entropy(rep(1, clst_num)/clst_num) 
  train <- svm(train_x, as.factor(train_y), probability = TRUE)
  b5 <- predict(train, test ,probability = TRUE)
  b5 <- attr(b5, "probabilities")
  res$prob <- b5
  res$clst <- test_clst
  
  res$score <- aggregate(. ~clst, res, mean)
  res$alignment <- calcEntropy(res$score)
  #print(res$score)
  return(res)
}


calcEntropy<-function(x=NULL){
  y<-matrix(0, nrow(x),1)
  for(i in seq(1,nrow(x),1)){
    y[i] <- Entropy(x[i, seq(2, ncol(x),1)])
  }
  return(y)
}


knn_cluster <-function(x, k=20, prune.SNN =1/15, resolution = 0.5){

  #res <- irlba(crossprod(x,x), nu =15, nv =15)
  snn <- swne::CalcSNN(t(x), k = k, prune.SNN = prune.SNN)
  ids <- RunModularityClustering(SNN = snn, resolution = resolution) 
  return(ids)
}



RunModularityClustering <- function(
  SNN = matrix(),
  modularity = 1,
  resolution = 0.8,
  algorithm = 1,
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  print.output = TRUE,
  temp.file.location = NULL,
  edge.file.name = NULL
) {
  #edge_file <- edge.file.name %||% ''
  edge_file <- ''
  clusters <- RunModularityClusteringCpp(
    SNN,
    modularity,
    resolution,
    algorithm,
    n.start,
    n.iter,
    random.seed,
    print.output,
    edge_file
  )
  return(clusters)
}


RunModularityClusteringCpp <- function(SNN, modularityFunction, resolution, algorithm, nRandomStarts, nIterations, randomSeed, printOutput, edgefilename) {
    .Call('_Seurat_RunModularityClusteringCpp', PACKAGE = 'Seurat', SNN, modularityFunction, resolution, algorithm, 
      nRandomStarts, nIterations, randomSeed, printOutput, edgefilename)
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



#' Run Leiden clustering algorithm
#' This code is modified from Tom Kelly (https://github.com/TomKellyGenetics/leiden), where we added more parameters (seed.use and n.iter) to run the Python version. In addition, we also take care of the singleton issue after running leiden algorithm.
#' @description Implements the Leiden clustering algorithm in R using reticulate to run the Python version. Requires the python "leidenalg" and "igraph" modules to be installed. Returns a vector of partition indices.
#' @param SNN An adjacency matrix compatible with \code{\link[igraph]{igraph}} object.
#' @param seed.use set seed
#' @param n.iter number of iteration
#' @param initial.membership arameters to pass to the Python leidenalg function defaults initial_membership=None
#' @param node.sizes Parameters to pass to the Python leidenalg function
#' @param partition_type Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition (see the Leiden python module documentation for more details)
#' @param resolution A parameter controlling the coarseness of the clusters
#' @param weights defaults weights=None
#' @return A parition of clusters as a vector of integers
##'
##'
#' @keywords graph network igraph mvtnorm simulation
#' @importFrom reticulate import r_to_py
##' @export

runLeiden <- function(SNN = matrix(), resolution = 1, partition_type = c(
  'RBConfigurationVertexPartition',
  'ModularityVertexPartition',
  'RBERVertexPartition',
  'CPMVertexPartition',
  'MutableVertexPartition',
  'SignificanceVertexPartition',
  'SurpriseVertexPartition'),
seed.use = 42L,
n.iter = 10L,
initial.membership = NULL, weights = NULL, node.sizes = NULL) {
  if (!reticulate::py_module_available(module = 'leidenalg')) {
    stop("Cannot find Leiden algorithm, please install through pip (e.g. pip install leidenalg).")
  }

  #import python modules with reticulate
  leidenalg <- import("leidenalg", delay_load = TRUE)
  ig <- import("igraph", delay_load = TRUE)

  resolution_parameter <- resolution
  initial_membership <- initial.membership
  node_sizes <- node.sizes
  #convert matrix input (corrects for sparse matrix input)
  adj_mat <- as.matrix(SNN)

  #compute weights if non-binary adjacency matrix given
  is_pure_adj <- all(as.logical(adj_mat) == adj_mat)
  if (is.null(weights) && !is_pure_adj) {
    #assign weights to edges (without dependancy on igraph)
    weights <- t(adj_mat)[t(adj_mat)!=0]
    #remove zeroes from rows of matrix and return vector of length edges
  }

  ##convert to python numpy.ndarray, then a list
  adj_mat_py <- r_to_py(adj_mat)
  adj_mat_py <- adj_mat_py$tolist()

  #convert graph structure to a Python compatible object
  GraphClass <- if (!is.null(weights) && !is_pure_adj){
    ig$Graph$Weighted_Adjacency
  } else {
    ig$Graph$Adjacency
  }
  snn_graph <- GraphClass(adj_mat_py)

  #compute partitions
  partition_type <- match.arg(partition_type)
  part <- switch(
    EXPR = partition_type,
    'RBConfigurationVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$RBConfigurationVertexPartition,
      initial_membership = initial.membership, weights = weights,
      resolution_parameter = resolution,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'ModularityVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$ModularityVertexPartition,
      initial_membership = initial.membership, weights = weights,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'RBERVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$RBERVertexPartition,
      initial_membership = initial.membership, weights = weights, node_sizes = node.sizes,
      resolution_parameter = resolution_parameter,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'CPMVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$CPMVertexPartition,
      initial_membership = initial.membership, weights = weights, node_sizes = node.sizes,
      resolution_parameter = resolution,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'MutableVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$MutableVertexPartition,
      initial_membership = initial.membership,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'SignificanceVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$SignificanceVertexPartition,
      initial_membership = initial.membership, node_sizes = node.sizes,
      resolution_parameter = resolution,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'SurpriseVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$SurpriseVertexPartition,
      initial_membership = initial.membership, weights = weights, node_sizes = node.sizes,
      n_iterations = n.iter,
      seed = seed.use
    ),
    stop("please specify a partition type as a string out of those documented")
  )
  partition <- part$membership+1
  idents <- partition

  if (min(table(idents)) == 1) {
    idents <- assignSingletons(idents, SNN)
  }
  idents <- factor(idents)
  names(idents) <- row.names(SNN)
  return(idents)
}

# Group single cells that make up their own cluster in with the cluster they are most connected to.
#
# @param idents  clustering result
# @param SNN     SNN graph
# @return        Returns scAI object with all singletons merged with most connected cluster

assignSingletons <- function(idents, SNN) {
  # identify singletons
  singletons <- c()
  for (cluster in unique(idents)) {
    if (length(which(idents %in% cluster)) == 1) {
      singletons <- append(x = singletons, values = cluster)
    }
  }
  #singletons = names(table(idents))[which(table(idents)==1)]
  # calculate connectivity of singletons to other clusters, add singleton to cluster it is most connected to
  cluster_names <- unique(x = idents)
  cluster_names <- setdiff(x = cluster_names, y = singletons)
  connectivity <- vector(mode="numeric", length = length(x = cluster_names))
  names(x = connectivity) <- cluster_names
  for (i in singletons) {
    print(i)
    for (j in cluster_names) {
      subSNN = SNN[
        which(idents %in% i),
        which(idents %in% j)
        ]
      if (is.object(x = subSNN)) {
        connectivity[j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
      } else {
        connectivity[j] <- mean(x = subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(x = connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
    which(idents %in% i)[which(idents %in% i)] <- closest_cluster
  }
  if (length(x = singletons) > 0) {
    message(paste(
      length(x = singletons),
      "singletons identified.",
      length(x = unique(idents)),
      "final clusters."
    ))
  }
  return(idents)
}



##########################################################################################
# Plot Aesthetics Objects and Methods
##########################################################################################

#' List of color palettes that can be used in plots
#' 
#' A collection of some original and some borrowed color palettes to provide appealing color aesthetics for plots in ArchR
#' 
#' @export
bindSCPalettes <- list(

  #DISCLOSURE: This is a collection of palettes that includes some original palettes and some palettes originally
  #implemented by others in other packages.
  #They are included here for convenience because they help improve plot aesthetics.

  #NOTE: all palettes included in the "Primarily Continuous Palettes" section should also work for discrete usage but not vice versa.
  #Each continuous palette has been ordered by color to generate a visually appealing discrete palette.
  
  #---------------------------------------------------------------
  # Primarily Discrete Palettes
  #---------------------------------------------------------------
  
  #20-colors
  stallion = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D"),

  stallion2 = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767"),

  calm = c("1"="#7DD06F", "2"="#844081", "3"="#688EC1", "4"="#C17E73", "5"="#484125", "6"="#6CD3A7", "7"="#597873","8"="#7B6FD0", "9"="#CF4A31", "10"="#D0CD47",
          "11"="#722A2D", "12"="#CBC594", "13"="#D19EC4", "14"="#5A7E36", "15"="#D4477D", "16"="#403552", "17"="#76D73C", "18"="#96CED5", "19"="#CE54D1", "20"="#C48736"),

  kelly = c("1"="#FFB300", "2"="#803E75", "3"="#FF6800", "4"="#A6BDD7", "5"="#C10020", "6"="#CEA262", "7"="#817066", "8"="#007D34", "9"="#F6768E", "10"="#00538A",
          "11"="#FF7A5C", "12"="#53377A", "13"="#FF8E00", "14"="#B32851", "15"="#F4C800", "16"="#7F180D", "17"="#93AA00", "18"="#593315", "19"="#F13A13", "20"="#232C16"),

  #16-colors
  bear = c("1"="#faa818", "2"="#41a30d","3"="#fbdf72", "4"="#367d7d",  "5"="#d33502", "6"="#6ebcbc", "7"="#37526d",
           "8"="#916848", "9"="#f5b390", "10"="#342739", "11"="#bed678","12"="#a6d9ee", "13"="#0d74b6",
           "14"="#60824f","15"="#725ca5", "16"="#e0598b"),
  
  #15-colors
  ironMan = c("9"='#371377',"3"='#7700FF',"2"='#9E0142',"10"='#FF0080', "14"='#DC494C',"12"="#F88D51","1"="#FAD510","8"="#FFFF5F","4"='#88CFA4',
           "13"='#238B45',"5"="#02401B", "7"="#0AD7D3","11"="#046C9A", "6"="#A2A475", "15"='grey35'),
  
  circus = c("1"="#D52126", "2"="#88CCEE", "3"="#FEE52C", "4"="#117733", "5"="#CC61B0", "6"="#99C945", "7"="#2F8AC4", "8"="#332288",
             "9"="#E68316", "10"="#661101", "11"="#F97B72", "12"="#DDCC77", "13"="#11A579", "14"="#89288F", "15"="#E73F74"),

  #12-colors
  paired = c("9"="#A6CDE2","1"="#1E78B4","3"="#74C476","12"="#34A047","11"="#F59899","2"="#E11E26",
               "10"="#FCBF6E","4"="#F47E1F","5"="#CAB2D6","8"="#6A3E98","6"="#FAF39B","7"="#B15928"),
  
  #11-colors
  grove = c("11"="#1a1334","9"="#01545a","1"="#017351","6"="#03c383","8"="#aad962","2"="#fbbf45","10"="#ef6a32","3"="#ed0345","7"="#a12a5e","5"="#710162","4"="#3B9AB2"),
  
  #7-colors
  summerNight = c("1"="#2a7185", "2"="#a64027", "3"="#fbdf72","4"="#60824f","5"="#9cdff0","6"="#022336","7"="#725ca5"),
  
  #5-colors
  zissou = c("1"="#3B9AB2", "4"="#78B7C5", "3"="#EBCC2A", "5"="#E1AF00", "2"="#F21A00"), #wesanderson
  darjeeling = c("1"="#FF0000", "2"="#00A08A", "3"="#F2AD00", "4"="#F98400", "5"="#5BBCD6"), #wesanderson
  rushmore = c("1"="#E1BD6D", "5"="#EABE94", "2"="#0B775E", "4"="#35274A" , "3"="#F2300F"), #wesanderson
  captain = c("1"="grey","2"="#A1CDE1","3"="#12477C","4"="#EC9274","5"="#67001E"),

  #---------------------------------------------------------------
  # Primarily Continuous Palettes
  #---------------------------------------------------------------
  
  #10-colors
  horizon = c("1"='#000075',"4"='#2E00FF', "6"='#9408F7', "10"='#C729D6', "8"='#FA4AB5', "3"='#FF6A95', "7"='#FF8B74', "5"='#FFAC53', "9"='#FFCD32', "2"='#FFFF60'),
  
  #9-colors
  horizonExtra =c("1"="#000436","4"="#021EA9","6"="#1632FB","8"="#6E34FC","3"="#C732D5","9"="#FD619D","7"="#FF9965","5"="#FFD32B","2"="#FFFC5A"),
  blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D"),
  sambaNight = c("6"='#1873CC',"2"='#1798E5',"8"='#00BFFF',"5"='#4AC596',"1"='#00CC00',"4"='#A2E700',"9"='#FFFF00',"7"='#FFD200',"3"='#FFA500'), #buencolors
  solarExtra = c("5"='#3361A5', "7"='#248AF3', "1"='#14B3FF', "8"='#88CEEF', "9"='#C1D5DC', "4"='#EAD397', "3"='#FDB31A',"2"= '#E42A2A', "6"='#A31D1D'),  #buencolors
  whitePurple = c("9"='#f7fcfd',"6"='#e0ecf4',"8"='#bfd3e6',"5"='#9ebcda',"2"='#8c96c6',"4"='#8c6bb1',"7"='#88419d',"3"='#810f7c',"1"='#4d004b'),
  whiteBlue = c("9"='#fff7fb',"6"='#ece7f2',"8"='#d0d1e6',"5"='#a6bddb',"2"='#74a9cf',"4"='#3690c0',"7"='#0570b0',"3"='#045a8d',"1"='#023858'),
  whiteRed = c("1"="white", "2"="red"),
  comet = c("1"="#E6E7E8","2"="#3A97FF","3"="#8816A7","4"="black"),

  #7-colors
  greenBlue = c("4"='#e0f3db',"7"='#ccebc5',"2"='#a8ddb5',"5"='#4eb3d3',"3"='#2b8cbe',"6"='#0868ac',"1"='#084081'),
  
  #6-colors
  beach = c("4"="#87D2DB","1"="#5BB1CB","6"="#4F66AF","3"="#F15F30","5"="#F7962E","2"="#FCEE2B"),
  
  #5-colors
  coolwarm = c("1"="#4858A7", "4"="#788FC8", "5"="#D6DAE1", "3"="#F49B7C", "2"="#B51F29"),
  fireworks = c("5"="white","2"="#2488F0","4"="#7F3F98","3"="#E22929","1"="#FCB31A"),
  greyMagma = c("2"="grey", "4"="#FB8861FF", "5"="#B63679FF", "3"="#51127CFF", "1"="#000004FF"),
  fireworks2 = c("5"="black", "2"="#2488F0","4"="#7F3F98","3"="#E22929","1"="#FCB31A"),
  purpleOrange = c("5"="#581845", "2"="#900C3F", "4"="#C70039", "3"="#FF5744", "1"="#FFC30F")
)

#' Optimized discrete color palette generation
#'
#' This function assesses the number of inputs and returns a discrete color palette that is tailored to provide the most
#' possible color contrast from the designated color set.
#'
#' @param values A character vector containing the sample names that will be used. Each entry in this character vector will be
#' given a unique color from the designated palette set.
#' @param set The name of a color palette provided in the `bindSCPalettes` list object.
#' @param reverse A boolean variable that indicates whether to return the palette colors in reverse order.
#' @export
paletteDiscrete <- function(
  values = NULL,
  set = "stallion",  
  reverse = FALSE
  ){

  
  values <- gtools::mixedsort(values)
  n <- length(unique(values))
  pal <- bindSCPalettes[[set]]
  palOrdered <- pal[gtools::mixedsort(names(pal))] #mixed sort gets 1,2,3,4..10,11,12

  if(n > length(palOrdered)){
    message("Length of unique values greater than palette, interpolating..")
    palOut <- colorRampPalette(pal)(n)
  }else{
    palOut <- palOrdered[seq_len(n)]
  }
  
  if(reverse){
    palOut <- rev(palOut)
  }

  names(palOut) <- unique(values)

  return(palOut)
  
}

#' Continuous Color Palette
#'
#' @param set The name of a color palette provided in the `bindSCPalettes` list object.
#' @param n The number of unique colors to generate as part of this continuous color palette.
#' @param reverse A boolean variable that indicates whether to return the palette colors in reverse order.
#' @export
paletteContinuous <- function(
  set = "solarExtra", 
  n = 256, 
  reverse = FALSE
  ){

  
  pal <- bindSCPalettes[[set]]
  palOut <- colorRampPalette(pal)(n)
  
  if(reverse){  
    palOut <- rev(palOut)
  }

  return(palOut)
  
}



plot_geneScoreChange <- function(X=NULL, Z0=NULL,Z_impu=NULL){

  X<- (X-min(X))/(max(X)-min(X))
  Z0<- (Z0-min(Z0))/(max(Z0)-min(Z0))
  Z_impu<- (Z_impu-min(Z_impu))/(max(Z_impu)-min(Z_impu))
  before_bindSC <- data.frame("True"=as.vector(X), "Est"=as.vector(Z0))
 
  val <-cor(before_bindSC$True, before_bindSC$Est)
  p21 <- ggplot(before_bindSC, aes(x=True, y=Est) ) + geom_bin2d() + theme_classic() + 
    xlab("True") + ylab("Initilized gene score matrix") + 
    stat_cor(label.y=1) + 
    scale_fill_gradient(low="lightblue1",high="darkblue",trans="log10")

  after_bindSC <- data.frame("True"=as.vector(X), "Est"=as.vector(Z_impu))
  val <-cor(after_bindSC$True, after_bindSC$Est)
  p22 <- ggplot(after_bindSC, aes(x=True, y=Est) ) + geom_bin2d() + theme_classic() + 
    xlab("True") + ylab("Imputed gene score matrix") + 
    stat_cor(label.y=1) + 
    scale_fill_gradient(low="lightblue1",high="darkblue",trans="log10")

  p4 <- ggarrange(p21, p22, ncol=2)
  return(p4)

}


