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
#' dimension shared by two datasets. The BiCCA function introudces one transition matrix Z (M features by L
#' samples) to bridge X with Y. The transition matrix Z is solved by maximalizing correlation between (X, Z) 
#' in sample level and correlation between (Z, Y) in feature level simultaneously. Then sample/feature level 
#' correlation can be obtained by applying standard CCA on (X, Z) and (Y, Z), respectively.


#'
#' @param X First matrix (M features by K samples)
#' @param Y Second matrix (N features by L samples)
#' @param Z0 Transition matrix (M features by L samples) 
#' @param alpha Couple coefficient to assign weight on initilized gene activity matrix Z0 [default 0.5]. alpha is set to [0 1). 
#' If alpha is set to 0, the initilized Z0 will not be involved after the 2nd iteration. 
#' @param lambda Couple coefficient to assign a weight factor to multi-objective function [default 0.5]. lambda is set to (0 1).
#' If lambda is set to 0, the first modality will not be involved during the iteration process. 
#' If lambda is set to 1, the second modality will not be involved during the iteration process. In this case, biCCA will reduces to traditional CCA.
#' @param X.clst Cluster ID vector for the first matrix 
#' @param Y.clst Cluster ID vector for the second matrix 
#' @param K Number of canonical vectors to calculate in low-dimension space [default 10]
#' @param num.iteration Maximal number of iteration [default 100]
#' @param tolerance  Relative change ratio for frobenius norm of matrix Z during iteration [default 0.05]
#' @param block.size Sample/feature size for each block, only works when bigMemory is set to TRUE
#' @param save  Save the temporay files [default FALSE]
#' @param temp.path Folder that is used to store temporary files. Only works when option save is set to be TRUE [default NULL]
#'
#' @return Returns the object with list:
#'
#' @return * `Z` - contains the estimated transition matrix (M features by L samples)
#' @return * `u` - contains the canonical correlation vectors for X (K samples by K factor) 
#' @return * `r` - contains the canonical correlation vectors for Z (sample level)(L samples by K factor) 
#' @return * `s` - contains the canonical correlation vectors for Z (feature level)(M features by K factor) 
#' @return * `v` - contains the canonical correlation vectors for Y (N features by K factor)
#' @return * `delta` - relative change ratio for frobenius norm of matrix Z during iteration
#'
#'
#' @export
#' @import umap irlba progress Matrix rdist ggplot2 bigstatsr

#'
#' @author Jinzhuang Dou, \email{Dou1@mdanderson.org} \email{jinzhuangdou198706@gmail.com}

#'
#' @examples

#' X <- sim$X
#' Y <- sim$Y
#' Z0 <- sim$Z0
#'
#' out <- BiCCA(X=X, Z0=Z0, Y=Y, 
#'       K = 5, X.clst =  sim$X_meta$Group,
#'       Y.clst =  sim$Y_meta$Group,
#'       alpha = 0.5,
#'       lambda = 0,
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
    block.size = 5000){

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
    
    #if(parameter.optimize == FALSE){
    #    message(paste(current_time, "  KNN clustering on X ",sep=""))
    #    sink("temp.txt")
    #    invisible(capture.output (X_clst <- knn_cluster(X, resolution = 0.5)))
    #    sink()
    #}
   

    if(parameter.optimize == FALSE){
        message(paste(current_time, "  Decomposing started!",sep=""))
    }

    scale_3 <- norm(Z_in, type = "F")

    for(cnt in seq(1,num.iteration,1)){
        gc()
        if(parameter.optimize == FALSE) {pb$tick()}
    
        in1 <- crossprod(X, Z0)
        cca_l <-  irlba(in1, nv = K , nu = K)
        xu <- crossprod(t(X), cca_l$u)
        zv <- crossprod(t(Z0), cca_l$v)
        z_l <- crossprod(t(xu), t(cca_l$v))
       
        #cell_align <- celltype_assign(train_x  = cca_l$u, train_y = X.clst, test = cca_l$v, test_clst = Y.clst)
        #rd_cell_alignment[[cnt]] <- cell_align$score
        #rd_cell_alignment_mean <- c(rd_cell_alignment_mean, mean(cell_align$alignment))

        

        in2 <- crossprod(t(Z0), t(Y))
        cca_r <- irlba(in2, nv = K, nu = K)
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
    rd_cell_alignment <- calc_score(dt1 = cca_l$u, dt2 = cca_l$v, label1 = X.clst,  label2 = Y.clst)
    out<-list(
              "u" = cca_l$u, 
              "r" = cca_l$v,
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



RunModularityClusteringCpp <- function(SNN, modularityFunction, resolution, algorithm, nRandomStarts, nIterations, randomSeed, printOutput, edgefilename) {
    .Call('_Seurat_RunModularityClusteringCpp', PACKAGE = 'Seurat',
     SNN, modularityFunction, resolution, algorithm, nRandomStarts, 
     nIterations, randomSeed, printOutput, edgefilename)
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
    if(block.size <=100){
        stop("Please set bloc.size more than 100!")
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
