#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Bi-order Canonical Correlation Analysis (BiCCA)
#'
#' Assuming there are two matrices X (M features by K samples) and Y (N features by L samples), we want to 
#' find the correlation in feature and sample levels between X and Y. 
#' Standard CCA could not handle this case because it requires that there should be at least one 
#' dimension shared by two datasets. The BiCCA function introudces one gene activity matrix Z (M features by L
#' samples) to bridge X with Y. The gene score matrix Z is solved by maximalizing correlation between (X, Z) 
#' in cell level and correlation between (Z, Y) in feature level simultaneously. 

#'
#' @param X First matrix (M features by K samples)
#' @param Y Second matrix (N features by L samples)
#' @param Z0 Transition matrix (M features by L samples) 
#' @param alpha Couple coefficient to assign weight on initilized gene activity matrix Z0 [default 0.5]. alpha is set to [0 1). 
#' If alpha is set to 0, the initilized Z0 will not be involved after the 2nd iteration. 
#' @param lambda Couple coefficient to assign a weight factor to multi-objective function [default 0.5]. lambda is set to (0 1).
#' If lambda is close to 0, the first modality will not be involved during the iteration process. 
#' If lambda is close to 1, the second modality will not be involved during the iteration process. In this case, biCCA will reduce to the traditional CCA.
#' @param X.clst Cluster ID vector for the first matrix 
#' @param Y.clst Cluster ID vector for the second matrix 
#' @param K Number of canonical vectors to calculate in low-dimension space [default 10]
#' @param num.iteration Maximal number of iteration [default 100]
#' @param tolerance  Relative change ratio for frobenius norm of matrix Z during iteration [default 0.01]
#' @param block.size Block size for SVD decomposition. This will reduce time and memory usage when cell number is large
#' @param save  Save the temporay files [default FALSE]
#' @param temp.path Folder that is used to store temporary files. Only works when option save is set to be TRUE [default NULL]
#'
#' @return Returns the object with list:
#'
#' @return * `Z_est` - contains the updated gene score matrix (M features by L samples)
#' @return * `u` - contains the canonical correlation vectors for X (K samples by K factor) 
#' @return * `r` - contains the canonical correlation vectors for Z (sample level)(L samples by K factor) 
#' @return * `s` - contains the canonical correlation vectors for Z (feature level)(M features by K factor) 
#' @return * `v` - contains the canonical correlation vectors for Y (N features by K factor)
#' @return * `cost_l` - the cost function on the first modality penalty term 
#' @return * `cost_r` - the cost function on the second modality penalty term 
#' @return * `cost_z0` - the cost function on the the initialized gene score penalty
#' @return * `delta` - relative change ratio for frobenius norm of matrix Z during iteration
#'
#'
#' @export
#' @import umap irlba progress Matrix rdist ggplot2 

#'
#' @author Jinzhuang Dou, \email{Dou1@mdanderson.org} \email{jinzhuangdou198706@gmail.com}

#'
#' @examples

#' X <- sim$X
#' Y <- sim$Y
#' Z0 <- sim$Z0
#'
#' out <- BiCCA(X=X, Z0=Z0, Y=Y, 
#'       K = 5,
#'       alpha = 0.5,
#'       lambda = 0,
#'       X.clst = sim$X_meta$Group,
#'       Y.clst = sim$Y_meta$Group,  
#'       num.iteration = 100, 
#'       temp.path = "./tp",
#'       tolerance = 0.01, 
#'       save = TRUE, 
#'       block.size = 500)



BiCCA <- function(
    X=NULL, 
    Z0=NULL, 
    Y=NULL, 
    X.clst = NULL, 
    Y.clst = NULL,
    parameter.optimize = FALSE, 
    alpha = NULL,  
    lambda = NULL, 
    K = 5, 
    temp.path=NULL,
    num.iteration = 100, 
    tolerance = 0.01,  
    save = FALSE, 
    calc.score=FALSE,
    block.size = 0){

    start_time <- Sys.time()
    if(parameter.optimize == FALSE){
        message(paste(start_time, " Started!"))
    }


    preCheck(X=X, 
             Z0=Z0, 
             Y=Y, 
             K = K, 
             alpha = alpha, 
             lambda = lambda,
             num.iteration = num.iteration, 
             X.clst = X.clst, 
             Y.clst = Y.clst, 
             temp.path = temp.path, 
             tolerance = tolerance, 
             save = save, 
             block.size = block.size)

    current_time <- Sys.time()
    
    INF <- 10^(10)
    X_Ncell <- dim(X)[2]
    X_Ngene <- dim(X)[1]
    Y_Ncell <- dim(Y)[2]
    Y_Nloci <- dim(Y)[1]
    Z_Ncell <- dim(Z0)[2]
    Z_Ngene <- dim(Z0)[1]


    if(parameter.optimize == FALSE){
        message(paste(current_time, "  Dimension Check: X[",X_Ngene,"x",X_Ncell,"]", 
                                    " Y[",Y_Nloci,"x",Y_Ncell,"]",
                                    " Z0[",Z_Ngene,"x",Z_Ncell,"]", sep=""))
    }
    if(parameter.optimize==TRUE){
	   calc.score <- TRUE
    }


    gc()

    if(save==TRUE && !is.na(temp.path)){
        if (!dir.exists(temp.path)){
                dir.create(temp.path)
                if(parameter.optimize==FALSE){
                    message(paste(current_time, "  The output will be saved in ", temp.path, sep=""))
                }
        }
     } 


    # initialization 
    Z0 <- as.matrix(Z0)
    Z0 <- Z0/norm(Z0, type="F")
    cnt <- 0 
    delta <- INF
    Z_old <- Z0
    Z_in  <- Z0
    rd_delta <- c()
    rd_cost_l <- c()
    rd_cost_r <- c()
    rd_cost_z0 <- c()
    rd_cost_all <- c()
    rd_cell_alignment <- list()
    rd_cell_alignment_mean <- c()
    
    if(parameter.optimize==FALSE){
        pb <- progress_bar$new(
              format = "  Iterating [:bar] :percent eta: :eta",
              total = num.iteration, clear = FALSE, width= 60)
    }
    

    if(parameter.optimize == FALSE){
        message(paste(current_time, "  Decomposing started!",sep=""))
    }

    scale_3 <- norm(Z_in, type = "F")

    for(cnt in seq(1,num.iteration,1)){
        gc()
        if(parameter.optimize == FALSE) {pb$tick()}
    
        if(block.size==0){
          in1 <- crossprod(X, Z0)
          cca_l <-  irlba(in1, nv = K , nu = K)
        }
        else{
          cca_l <- irlba_block(x = X, y=Z0, k= K, blocksize = block.size)
        }

        xu <- crossprod(t(X), cca_l$u)
        zv <- crossprod(t(Z0), cca_l$v)
        z_l <- crossprod(t(xu), t(cca_l$v))
       
        if(block.size==0){
           in2 <- crossprod(t(Z0), t(Y))
           cca_r <- irlba(in2, nv = K, nu = K)
        }
        else{
           cca_r <- irlba_block(x = t(Z0), y=t(Y), k = K, blocksize = block.size)
        }

        zs <- crossprod(Z0, cca_r$u)
        yt <- crossprod(Y, cca_r$v) 
        z_r <- crossprod(t(cca_r$u), t(yt))
        scale_r <- norm(z_r, type = "F")

        scale_x  <- norm(z_l, type = "F")
        scale_y  <- norm(z_r, type = "F")
        cost_l  <- sum(diag(crossprod(xu, zv)))/scale_x
        cost_r  <- sum(diag(crossprod(zs, yt)))/scale_y
  


        Z0 <- (z_l*lambda/scale_x + z_r*(1-lambda)/scale_y)*(1-alpha) + Z_in*alpha/scale_3
        delta <- norm(Z0 - Z_old, type = "F")/norm(Z_old, type = "F")
        cost_z0 <-  norm(Z0 - Z_in, type = "F")/scale_3
        Z_old  <- Z0
        cost_all <- (cost_l*lambda + cost_r*(1-lambda))*(1-alpha) + cost_z0 * alpha
        
        rd_cost_l <- c(rd_cost_l, cost_l)
        rd_cost_r <- c(rd_cost_r, cost_r)
        rd_delta <- c(rd_delta, delta)
        rd_cost_z0 <- c(rd_cost_z0, cost_z0)
        rd_cost_all <- c(rd_cost_all, cost_all)
       
        cca_l$u <- L2Norm(cca_l$u)
        cca_l$v <- L2Norm(cca_l$v)

        if(save==TRUE && !is.na(temp.path)){
            if (!dir.exists(temp.path)){
                dir.create(temp.path)
            } 

            set.seed(123)
  
            out<-list(
              "u" = cca_l$u, 
              "r" = cca_l$v,
              "s" = cca_r$u,
              "v" = cca_r$v,  
  	      "d" = cca_l$d,
              "Z_est" = Z0,   
              "cost_l" =  cost_l,
              "cost_r" =  cost_r,
              "cost_z0" = cost_z0,
              "delta" = delta)

            saveRDS(out, file=paste(temp.path,"/iteration",cnt,".rds",sep=""))
        }

        if(delta < tolerance){break}

    }

    if(parameter.optimize == FALSE){
        message("\n")
    }
    current_time <- Sys.time()
    if(delta >= tolerance){
        if(parameter.optimize == FALSE){
            message(paste(current_time, " The decomposition is not converged. The output may not be optimal.
                     You can increase num.iteration again."))
        }
    }
    else{
        if(parameter.optimize == FALSE){
             message(paste(current_time, " Done! The decomposition is converged."))
        }
    }
    if(calc.score){
    
     rd_cell_alignment <- calc_score(dt1 = cca_l$u, dt2 = cca_l$v, label1 = X.clst,  label2 = Y.clst)
    }
    else{
      rd_cell_alignment <- NULL
    }

   
    if(!is.na(colnames(X)[1])) {rownames(cca_l$u) <- colnames(X)}
    if(!is.na(colnames(Y)[1])) {rownames(cca_l$v) <- colnames(Y)}

    out<-list(
              "u" = cca_l$u, 
              "r" = cca_l$v,
	      "d" = cca_l$d,
              "s" = cca_r$u,
              "v" = cca_r$v,
              "Z_est" = Z0,   
              "cost_l" =  rd_cost_l,
              "cost_r" =  rd_cost_r,
              "cost_z0" = rd_cost_z0,
              "cost_all" = rd_cost_all,
              "delta" = rd_delta,
              "score" = rd_cell_alignment )

    saveRDS(out, file=paste(temp.path,"/out_final.rds",sep=""))
    return(out)
}





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


impuZ <- function(X=NULL, bicca=NULL){
    out <- X%*%bicca$u%*%t(bicca$r)
    return(out)
}


dimReduce <- function(dt1=NULL, dt2=NULL, K=50){
  if( is.null(dt1)){
      stop("Please input dt1!")
  }
  else{
    res <- list()
    if(is.null(dt2)){  
        out <- irlba(crossprod(dt1,dt1), nu =K, nv =K)
        rownames(out$u) <- colnames(dt1)
        colnames(out$u) <- paste("C",seq(1,K,1), sep="")
        res$dt1 <- out$u
        return(res$dt1)
    }
    else{
        out <- irlba(crossprod(dt1,dt2), nu =K, nv =K)
        rownames(out$u) <- colnames(dt1)
        rownames(out$v) <- colnames(dt2)
        colnames(out$u) <- paste("C",seq(1,K,1), sep="")
        colnames(out$v) <- paste("C",seq(1,K,1), sep="")
        res$dt1 <- out$u
        res$dt2 <- out$v
        return(res)
    }
  }
}


BiCCA_para_opt <-function(X=NULL, 
    Z0=NULL, 
    Y=NULL, 
    X.clst = NULL, 
    Y.clst = NULL,
    temp.path=NULL,
    num.iteration = 100, 
    tolerance = 0.01,  
    save = TRUE, 
    tp = "out",
    alpha.lst = NULL,
    lambda.lst = NULL,
    ncore = NCORES,
    block.size = NULL,
    calc.score = FALSE,
    K.lst = NULL){

    start_time <- Sys.time()
    message(paste(start_time, " Started!"))

    current_time <- Sys.time()


    preCheck(X=X, 
         Z0=Z0, 
         Y=Y, 
         K = K.lst, 
         alpha = alpha.lst, 
         lambda = lambda.lst,
         num.iteration = num.iteration, 
         X.clst = X.clst, 
         Y.clst = Y.clst, 
         temp.path = tp, 
         tolerance = tolerance, 
         save = save, 
         block.size = block.size)
    


    INF <- 10^(10)
    M <- 200
    X_Ncell <- dim(X)[2]
    X_Ngene <- dim(X)[1]
    Y_Ncell <- dim(Y)[2]
    Y_Nloci <- dim(Y)[1]
    Z_Ncell <- dim(Z0)[2]
    Z_Ngene <- dim(Z0)[1]

    
    message(paste(current_time, "  Dimension Check: X[",X_Ngene,"x",X_Ncell,"]", 
                                " Y[",Y_Nloci,"x",Y_Ncell,"]",
                                " Z0[",Z_Ngene,"x",Z_Ncell,"]", sep=""))

    
 
    #message(paste(current_time, "  KNN clustering on X ",sep=""))
    #sink("temp.txt")
    #invisible(capture.output (X_clst <- knn_cluster(t(X), resolution = 0.5)))
    #sink()

    tol <- length(K.lst) * length(alpha.lst) * length(lambda.lst)
    pb <- progress_bar$new(
          format = "  Parameter optimization [:bar] :percent eta: :eta",
          total = tol, clear = FALSE, width= 60)

    rd <- c()
    for(K in K.lst){
        for(alpha in alpha.lst){
            for(lambda in lambda.lst){
          
                #message(paste(current_time, "  ","Decomposing with [K=",K.lst,",alpha=",alpha, ",lambda=", lambda,  "] ...", sep=""))
                skip_to_next <- FALSE
                pb$tick()
                tp.dir <- paste(tp,K,alpha,lambda, sep="_")
                tryCatch( 
                    res <- BiCCA(X=X, 
                        Z0=Z0, Y=Y,
                        X.clst = X.clst,
                        Y.clst = Y.clst,
                        K = K, 
                        num.iteration = num.iteration, 
                        temp.path = tp.dir,
                        tolerance = tolerance,
                        alpha = alpha, 
                        lambda  = lambda, 
                        parameter.optimize = TRUE, 
                        save = TRUE), 
                        error = function(e) { 
                            print(e)
                            skip_to_next <<- TRUE}
                        )
               
                rd <- rbind(rd, 
                    c(K, alpha, lambda,  
                        mean(res$score$silhouette1[,3]), 
                        mean(res$score$silhouette2[,3]),  
                        mean(res$score$alignment),
                        tail(res$cost_all, n=1)
                    )
                )
                if(skip_to_next) { next }  
            }
        }
    }

   
    colnames(rd)<-c("K","alpha", "lambda" , "silhouette1", "silhouette2",  "alignment","cost_all")

    message(paste(current_time, " Parameter optimization is done!"))
    return(rd)

}





# R code for split-and-merge approach for Single Value Decomposition(SVD)
# X is our data variable
# svd(t(X)%*%Y)

 irlba_block <- function(x=NULL, y=NULL, blocksize=NULL, k = NULL){

    #x <- as.matrix(t(x))
    blocksize <- min(c(blocksize, ncol(x)))
    x <- t(x)

    m <- nrow(x)
    Y <- NULL               #Initializing the end Y matrix
    U1 <- NULL              #Initializing the U bar matrix form X
    xi <- NULL
    X1 <- NULL

    i <- 1    
    part_len <- blocksize

    #loop through the data
    block_cnt <- 0
    T <- k
    #T <- nrow(x)
    while ( i <= m ){
      block_cnt <- block_cnt + 1
      
      if( i+part_len-1 < m ){
        xi <- x[i:(i+part_len-1), ]%*%y       
      }
      else if(i != m) {
        xi <- x[i:m, ]%*%y
      }
      else{
        xi <-matrix((x[i,]), nrow = 1) %*%y     
      }
      fin<-NULL

      xi.svd <- irlba(xi, nv=T)
      #xi.svd <- svd(xi)                     #Performing svd on individual submatrices
      if(is.null(U1)){                      #When first submatrix is being used
        #fin<-xi.svd$u
        U1 <- xi.svd$u
      }
      else                                 #We need to add U matrices diagonally
      {
          U1<-rbind(U1,xi.svd$u)
      }
      #U1 <- fin                  #Updating U bar matrix
      d<-diag(xi.svd$d, nrow = length(xi.svd$d))
      yi <- d %*% t(xi.svd$v)    #Creating y=v*d
      Y <- rbind(Y,yi)           #Updating the y matrix
      i <- i+part_len
    }

    y.svd <- svd(Y)
    #y.svd <- svd(as.matrix(Y))
    #Assigning final matrices to a random variable to output.
    X1$u <- NULL
    for(j in seq(1,block_cnt,1)){
      row_start <- (j-1)*blocksize + 1
      row_end <- j*blocksize
      if(row_end >= nrow(U1)){
          row_end <- nrow(U1)
      }
      col_start <- (j-1)*T + 1 
      col_end <- j*T
      if(col_end >= nrow(y.svd$u)){
          col_end <- nrow(y.svd$u)
      }
      tp <- U1[seq(row_start,row_end,1),]%*%y.svd$u[seq(col_start,col_end,1),]
      X1$u <- rbind(X1$u, tp)
    }

    #X1$u <- U1 %*% y.svd$
    X1$u <- X1$u[,seq(1,k,1)]
    X1$v <- y.svd$v[,seq(1,k,1)]
    #X1$d <- diag(y.svd$d, nrow = length(y.svd$d))
    X1$d <- y.svd$d[seq(1,k,1)]
     #Returning the three matrices through X1 
    return(X1)
  
}








label_transfer <- function(dt1 = NULL, X.clst = NULL, dt2 = NULL, k=3){
    
    dis<-cdist(as.matrix(dt2), as.matrix(dt1))
    # using k-neighbour to update coembeddings 
    #k <- 3
    for(i in seq(1, nrow(dis),1)){
          b<-sort(dis[i,])
          d<-rep(0, length(b))
          d[which(dis[i,]<=b[k])]<-1
          dis[i,]<-d/k
    }
    U <- dis%*%dt1
        

    model <- svm(dt1, as.factor(X.clst), probability=TRUE)
    pred <- predict(model, U, probability=TRUE)
    out <- attr(pred, "probabilities")
   
    celltype <- rep("unknown", nrow(out))
    prob_max <- rep(0, nrow(out))
    names <- colnames(out)
    for(i in seq(1,nrow(out),1)){
          tp <- out[i,]
          pos <- which(tp==max(tp),arr.ind = TRUE)
          celltype[i] <- names[pos]
          prob_max[i] <- tp[pos]
    }
    prob <- data.frame("celltype"=celltype, "Prob_max"=prob_max)
    rownames(prob) <- rownames(dt2)
    return(prob)
}

calc_score <- function(dt1 = NULL, dt2 = NULL, label1 = NULL, label2 = NULL ){

    dt3 <- rbind(dt1,dt2)
    dt1 <- as.matrix(dt1)
    dt2 <- as.matrix(dt2)
    
    out <-c()
    out$silhouette1 <- calc_silhouette_coef(x = dt1, clst = label1)
    out$silhouette2 <- calc_silhouette_coef(x = dt2, clst = label2)
    label3 <- c(rep("A", nrow(dt1)), rep("B", nrow(dt2)))
  
    out$alignment <- calc_alignment_score(x = dt3, k = 20, clst = label3)
    return(out)
}


# dt1 bindSC co-embedding: on RNA cells 
# dt2 bindSC co-embedding: on Merfish cells 
# X : gene expression matrix on RNA cells 
# k:  neighbor size
GeneImp <- function(dt1 = NULL, X = NULL, dt2 = NULL, k=NULL){
    
    dis<-cdist(as.matrix(dt2), as.matrix(dt1))
   # using k-neighbour to update coembeddings 
  
    for(i in seq(1, nrow(dis),1)){
          b<-sort(dis[i,])
          d<-rep(0, length(b))
          d[which(dis[i,]<=b[k])]<-1
          dis[i,]<-d/k
    }
    U <- dis%*%X
    return(U)
        
}



preCheck <- function(X=NULL, 
                    Z0=NULL, 
                    Y=NULL,  
                    K = NULL,
                    alpha = NULL,
                    lambda = NULL,
                    num.iteration = NULL,
                    X.clst = NULL,
                    Y.clst = NULL,
                    temp.path = NULL, 
                    tolerance = NULL, 
                    save = NULL,
                    block.size = NULL,
                    ncore = NCORES) {

    if (!is.matrix(X) && !is.sparseMatrix(X)){ff
        stop("Please input X as a matrix!")  
    }
    if (!is.matrix(Z0) && !is.sparseMatrix(Z0)){
        stop("Please input Z0 as a matrix!")
    }
    if (!is.matrix(Y) && !is.sparseMatrix(Y)){
        stop("Please input Y as a matrix!")
    }
    if(block.size <=100 & block.size>0){
        stop("Please set bloc.size more than 100 or set it to 0 (the block SVD mode is disable).")
    }
    if(length(alpha)==1){
        if(alpha < 0 || alpha >1){
            stop("Please set alpha in [0 1]!")
        }
        if(lambda < 0 || lambda >1){
            stop("Please set lambda in [0 1]!")
        }
    }
    if(dim(X)[1]!=dim(Z0)[1]){
        stop("X and Z0 should have the same row number!")
    }
    if(dim(Z0)[2]!=dim(Y)[2]){
        stop("Z0 and Y should have the same column number!")
    }
    if(dim(X)[2]!=length(X.clst)){
        stop("X and X.clst should have matched sample size!")
    }
    if(dim(Y)[2]!=length(Y.clst)){
        stop("Y and Y.clst should have matched sample size!")
    }

    if(sum(is.na(X.clst))>0){
        stop("X.clst has NA values!")
    }

    if(sum(is.na(Y.clst))>0){
        stop("Y.clst has NA values!")
    }


    tp <- min(dim(X)[2], dim(Z0)[2])
    if(!all.equal(K, as.integer(K)) || K<0 || K >tp){
        stop(paste("The K should be integer ranging from 2 to", tp))
    }

    tp <- min(dim(Z0)[1], dim(Y)[1])

    if(length(K)==1){
        if(!all.equal(K, as.integer(K)) || K<0 || K >tp){
            stop(paste("The K should be integer ranging from 2 to", tp))
        }
    }

    if(save && is.null(temp.path)){
        stop("Please specify the folder to store temporary output")
    }
    return(TRUE)
}
