library(tidyverse)
library(stringr)
library(cowplot)
library(monocle)
library(Matrix)
library(data.table)
library(Seurat)
library(ggpubr)
library(gridExtra)
mybg<-theme_bw() +theme(panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line()) +    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())

get_I<-function(peak.matrix,
  annotation.file,
  seq.levels = c(1:22, "X", "Y"),
  include.body = TRUE,
  upstream = 2000,
  downstream = 0,
  verbose = TRUE){
  if (!PackageCheck('GenomicRanges', error = FALSE)) {
    stop("Please install GenomicRanges from Bioconductor.")
  }
  if (!PackageCheck('rtracklayer', error = FALSE)) {
    stop("Please install rtracklayer from Bioconductor.")
  }
  # convert peak matrix to GRanges object
  peak.df <- rownames(x = peak.matrix)
  peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, pattern = ":", replacement = "-"), split = "-"))
  peak.df <- as.data.frame(x = peak.df)
  colnames(x = peak.df) <- c("chromosome", 'start', 'end')
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)

  # if any peaks start at 0, change to 1
  # otherwise GenomicRanges::distanceToNearest will not work
  BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1

  # get annotation file, select genes
  gtf <- rtracklayer::import(con = annotation.file)
  gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = 'coarse')

  # change seqlevelsStyle if not the same
  if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
    GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
  }
  gtf.genes <- gtf[gtf$type == 'gene']

  # Extend definition up/downstream
  if (include.body) {
    gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
  } else {
    gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, upstream = upstream, downstream = downstream)
  }
  gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
  keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]
  peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
  gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]

  # Some GTF rows will not have gene_name attribute
  # Replace it by gene_id attribute
  gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]

  peak.ids$gene.name <- gene.ids$gene_name
  peak.ids <- as.data.frame(x = peak.ids)
  peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
  annotations <- peak.ids[, c('peak', 'gene.name')]
  colnames(x = annotations) <- c('feature', 'new_feature')
  return(annotations)
}






construct_cds_ATAC <- function(UMI, df_cell, df_gene)
    {
    # accept a UMI count matrix, a df_cell data frame, a df_gene data frame, and then construct
    # a cds object for monocle downstream analysis
    # and then estimate the dispersion and size
    row.names(df_gene) = df_gene$peak
    row.names(df_cell) = df_cell$sample
    
    pd = new("AnnotatedDataFrame", data = df_cell)
    fd = new("AnnotatedDataFrame", data = df_gene)
    
    colnames(UMI) = df_cell$sample
    row.names(UMI) = df_gene$peak
    
    row.names(pd) = colnames(UMI)
    row.names(fd) = row.names(UMI) 
    
    cds = newCellDataSet(UMI, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
    return(cds)
}

construct_cds_RNA <- function(UMI, df_cell, df_gene)
    {
    # accept a UMI count matrix, a df_cell data frame, a df_gene data frame, and then construct
    # a cds object for monocle downstream analysis
    # and then estimate the dispersion and size
    row.names(df_gene) = df_gene$gene_short_name
    row.names(df_cell) = df_cell$sample
    
    pd = new("AnnotatedDataFrame", data = df_cell)
    fd = new("AnnotatedDataFrame", data = df_gene)
    
    colnames(UMI) = df_cell$sample
    row.names(UMI) = df_gene$gene_short_name
    
    row.names(pd) = colnames(UMI)
    row.names(fd) = row.names(UMI) 
    
    cds = newCellDataSet(UMI, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
    return(cds)
}

# here I define a function that accept a vector of peak name in the regression output, then it remove the first "X"
# and then generate a data frame with the chr, start, end location
df_peak_coor2 <- function(peak_names)
    {
    df_peak_name = data.frame(id = peak_names)
    df_peak_name = ((df_peak_name
                     %>% select(id) 
                     %>% separate(id, into = c("chr", "start", "stop"), sep = "-")))
    return(df_peak_name)
}


# accept a cds, and aggregate nearby peak within distance

make_bin_col <- function(cds, distance) {
    coords_string_df <- df_peak_coor2((fData(cds))$peak)
    coords_string_df$start = as.numeric(coords_string_df$start)
    coords_string_df$stop = as.numeric(coords_string_df$stop)
    
    coords_ranges <- GenomicRanges::makeGRangesFromDataFrame(coords_string_df)
    coords_range_merge <- GenomicRanges::reduce(coords_ranges, min.gapwidth = distance)

    merge_df <- data.frame(seqnames=GenomicRanges::seqnames(coords_range_merge),
                           starts=GenomicRanges::start(coords_range_merge),
                           ends=GenomicRanges::end(coords_range_merge))
    
    merge_df$name <- paste(merge_df$seqnames, merge_df$starts, merge_df$ends, sep="-")

    overlaps <- GenomicRanges::findOverlaps(coords_ranges,
                                            coords_range_merge,
                                            select="first")
    overlaps <- as.data.frame(overlaps)

    merge_df <- merge_df[overlaps$overlaps,]
    merge_df$name
}


#' Converts sparse matrix to data.table
#' @import data.table
#' @export

sparse_to_datatable <- function(sparse) {
  dgt_mat <- as(Matrix::t(sparse), "dgTMatrix")
  dt <- data.table::data.table(cell = dgt_mat@Dimnames[[1]][dgt_mat@i+1], site=dgt_mat@Dimnames[[2]][dgt_mat@j+1], val = dgt_mat@x)
  data.table::setkey(dt, site, cell)
  dt
}





aggregate_nearby_peaks <- function(cds, distance = 1000) {
    
    
  fData(cds)$bin <- make_bin_col(cds, distance)
  cds <- cds[!is.na(fData(cds)$bin),]

  exprs_dt <- sparse_to_datatable(Matrix(exprs(cds), sparse = TRUE))
  bin_info <- data.table::data.table(site = (fData(cds))$peak,
                                     bin = fData(cds)$bin)
  data.table::setkey(bin_info, site)
  data.table::setkey(exprs_dt, site)
  exprs_dt <- merge(exprs_dt, bin_info)

  data.table::setkey(exprs_dt, cell, bin)
  genomic_bins <- exprs_dt[,sum(val), by="cell,bin"]
  out <- Matrix::sparseMatrix(j=as.numeric(factor(genomic_bins$cell)),
                              i=as.numeric(factor(genomic_bins$bin)),
                              x=genomic_bins$V1)

  fdf <- data.frame(dhs = levels(factor(genomic_bins$bin)),
                    row.names = levels(factor(genomic_bins$bin)))
  pdf <- data.frame(cells = levels(factor(genomic_bins$cell)),
                    row.names = levels(factor(genomic_bins$cell)))
  fdf$bin <- NULL
  #pdf <- pdf[row.names(pData(cds)),]
  pdf <- cbind(pdf, pData(cds)[rownames(pdf), ])
  pdf$pdf <- NULL

  fd <- new("AnnotatedDataFrame", data = fdf)
  pd <- new("AnnotatedDataFrame", data = pdf)

  if (class(exprs(cds)) == "dgCMatrix") {
    compart_cds <-  newCellDataSet(as(out, "sparseMatrix"),
                                    phenoData = pd,
                                    featureData = fd,
                                    expressionFamily=negbinomial.size(),
                                    lowerDetectionLimit=1)
  } else {
    compart_cds <-  newCellDataSet(as.matrix(out),
                                   phenoData = pd,
                                   featureData = fd,
                                   expressionFamily=negbinomial.size(),
                                   lowerDetectionLimit=1)
  }
  return(compart_cds)
}


normalize_expr_data <- function(cds,
                                norm_method = c("log", "vstExprs", "none"),
                                pseudo_expr = 1,
                                relative_expr = TRUE){
  FM <- exprs(cds)
  use_for_ordering <- NULL
  # If the user has selected a subset of genes for use in ordering the cells
  # via setOrderingFilter(), subset the expression matrix.
  if (is.null(fData(cds)$use_for_ordering) == FALSE &&
      nrow(subset(fData(cds), use_for_ordering == TRUE)) > 0) {
    FM <- FM[fData(cds)$use_for_ordering, ]
  }

  norm_method <- match.arg(norm_method)
  if (cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {

    # If we're going to be using log, and the user hasn't given us a pseudocount
    # set it to 1 by default.
    if (is.null(pseudo_expr)){
      if(norm_method == "log")
        pseudo_expr = 1
      else
        pseudo_expr = 0
    }

    checkSizeFactors(cds)

    if (norm_method == "vstExprs") {
      if (relative_expr == FALSE)
        message("Warning: relative_expr is ignored when using norm_method == 'vstExprs'")

      if (is.null(fData(cds)$use_for_ordering) == FALSE &&
          nrow(subset(fData(cds), use_for_ordering == TRUE)) > 0) {
        VST_FM <- vstExprs(cds[fData(cds)$use_for_ordering,], round_vals = FALSE)
      }else{
        VST_FM <- vstExprs(cds, round_vals = FALSE)
      }

      if (is.null(VST_FM) == FALSE) {
        FM <- VST_FM
      }
      else {
        stop("Error: set the variance-stabilized value matrix with vstExprs(cds) <- computeVarianceStabilizedValues() before calling this function with use_vst=TRUE")
      }
    }else if (norm_method == "log") {
      # If we are using log, normalize by size factor before log-transforming
      
      if (relative_expr)
        FM <- Matrix::t(Matrix::t(FM)/sizeFactors(cds))

      if(is.null(pseudo_expr))
        pseudo_expr <- 1 
      FM <- FM + pseudo_expr
      FM <- log2(FM)
    }else if (norm_method == "none"){
      # If we are using log, normalize by size factor before log-transforming
      FM <- Matrix::t(Matrix::t(FM)/sizeFactors(cds))
      FM <- FM + pseudo_expr
    }
  }else if (cds@expressionFamily@vfamily == "binomialff") {
    if (norm_method == "none"){
      #If this is binomial data, transform expression values into TF-IDF scores.
      ncounts <- FM > 0
      ncounts[ncounts != 0] <- 1
      FM <- Matrix::t(Matrix::t(ncounts) * log(1 + ncol(ncounts)/rowSums(ncounts)))
    }else{
      stop("Error: the only normalization method supported with binomial data is 'none'")
    }
  }else if (cds@expressionFamily@vfamily == "Tobit") {
    FM <- FM + pseudo_expr
    if (norm_method == "none"){

    }else if (norm_method == "log"){
      FM <- log2(FM)
    }else{
      stop("Error: the only normalization methods supported with Tobit-distributed (e.g. FPKM/TPM) data are 'log' (recommended) or 'none'")
    }
  }else if (cds@expressionFamily@vfamily == "uninormal") {
    if (norm_method == "none"){
      FM <- FM + pseudo_expr
    }else{
      stop("Error: the only normalization method supported with gaussian data is 'none'")
    }
  }
  # if(norm_method != "none")
    #normalize_expr_data
  return (FM)
}

checkSizeFactors <- function(cds)
{
  if (cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size"))
  {
    if (is.null(sizeFactors(cds))){
      stop("Error: you must call estimateSizeFactors() before calling this function.")
    }
    if (sum(is.na(sizeFactors(cds))) > 0){
      stop("Error: one or more cells has a size factor of NA.")
    }
  }
}


FindIntegrationMatrix <- function(
  object,
  assay = NULL,
  integration.name = 'integrated',
  features.integrate = NULL,
  verbose = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = 'neighbors')
  nn.cells1 <- neighbors$cells1
  nn.cells2 <- neighbors$cells2
  anchors <- GetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'anchors'
  )
  if (verbose) {
    message("Finding integration vectors")
  }
  features.integrate <- features.integrate %||% rownames(
    x = GetAssayData(object = object, assay = assay, slot = "data")
  )
  data.use1 <- t(x = GetAssayData(
    object = object,
    assay = assay,
    slot = "data")[features.integrate, nn.cells1]
  )
  data.use2 <- t(x = GetAssayData(
    object = object,
    assay = assay,
    slot = "data")[features.integrate, nn.cells2]
  )
  anchors1 <- nn.cells1[anchors[, "cell1"]]
  anchors2 <- nn.cells2[anchors[, "cell2"]]
  data.use1 <- data.use1[anchors1, ]
  data.use2 <- data.use2[anchors2, ]
  integration.matrix <- data.use2 - data.use1
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'integration.matrix',
    new.data = integration.matrix
  )
  return(object)
}





