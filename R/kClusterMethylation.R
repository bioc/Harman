#' @title Cluster beta values with a set value for k
#' @description This function is part of a set of three functions to be run in
#' series. \code{\link{discoverClusteredMethylation}} takes a matrix of
#' methylation beta values (typically from the Illumina Infinium Methylation
#' Assay) and clusters the data across a range of ks specified buy the user.
#' 
#' Then the data is reclustered again across the the best two candidate values
#' for k (determined by the rate of change in Bayesian information criterion),
#' and minimum cluster size and distance filters are employed. If both clusters
#' meet these filters, then the higher value of k is returned.  
#' This function should be run on uncorrected data that ideally has slides
#' removed which are prone to batch effect. This will bias towards finding
#' clusters that are driven by biological factors such as X-chromosome
#' inactivation and allele-specific methylation.
#' 
#' The output of this function is input for the
#' \code{\link{kClusterMethylation}} function which extracts cluster membership
#' and statistics on variance for a given matrix of beta values. It might be
#' useful to discover clusters on samples less prone to clustering due to batch
#' effect or cellular heterogeneity and then recluster all the data for set
#' values of k via the \code{\link{kClusterMethylation}} function.
#' 
#' Finally, a comparison of differences of uncorrected to
#' batch-corrected beta values can be made using \code{\link{clusterStats}}.
#' This function generates a data.frame containing log variance ratio and mean
#' beta differences to clusters after correction.
#' @param betas a matrix of methylation beta values with samples in columns and
#' betas in rows.
#' @param row_ks vector, the value of k for each row in betas. The names
#' of \code{row_ks} must match the rownames of \code{betas}.
#' @param cores integer, specifies the number of cores to use for computation.
#' @param printInfo logical, whether to print information during computation or
#' not.
#' @return  A \code{kClusters} S3 object.
#' @details Betas values should be of type \code{double} with samples in
#' columns and betas in rows. The betas need to be bounded between 0 and 1.
#' The matrix is typically exported from a \code{\link[minfi]{GenomicRatioSet}},
#' \code{\link[minfi]{GenomicMethylSet}} or \code{\link[minfi]{MethylSet}}
#' object via the \code{getBeta} S4 accessor method.
#' @seealso \code{\link{discoverClusteredMethylation}},
#' \code{\link{clusterStats}}
#' @examples
#' library(HarmanData)
#' data(episcope)
#' bad_batches <- c(1, 5, 9, 17, 25)
#' is_bad_sample <- episcope$pd$array_num %in% bad_batches
#' myK <- discoverClusteredMethylation(episcope$original[, !is_bad_sample])
#' mykClust = kClusterMethylation(episcope$original, row_ks=myK)
#' res = clusterStats(pre_betas=episcope$original,
#'                    post_betas=episcope$harman,
#'                    kClusters = mykClust)
#' all.equal(episcope$ref_md$meandiffs_harman, res$meandiffs)
#' all.equal(episcope$ref_lvr$var_ratio_harman, res$log2_var_ratio)
#' @importFrom methods is
#' @importFrom Ckmeans.1d.dp Ckmeans.1d.dp
#' @importFrom parallel mclapply
#' @export
kClusterMethylation <- function(betas, row_ks, cores=1, printInfo=FALSE) {
  
  ######  Sanity checks  #####
  if(!methods::is(betas, "matrix")) {
    stop(paste("Require 'betas' as class matrix, not class \'",
               class(betas), "\'.", sep=""))
  }
  if(typeof(betas) != "double") {
    stop(paste("'betas' is type \'", typeof(betas),
               "\', needs to be of type \'double\'.", sep=""))
  }
  if(anyNA(betas)) {
    stop("The matrix cannot have NA values.")
  }
  range_bs <- range(betas)
  if(range_bs[1] < 0 | range_bs[2] > 1) {
    the_range <- paste(range_bs, collapse = ", ")
    stop(paste("The range of value in the 'betas' matrix is (", the_range,
               "), this is outside the range of possible legitimate beta values (0, 1).",
               sep=""))
  }
  rm(range_bs)
  
  # Check that the length of row_ks is equal to betas and that it's named
  if(!identical(rownames(betas), names(row_ks))) {
    stop('Betas rownames and names of row_ks are not identical')
  }
  
  # Check that not all clusters are 1
  if(sum(row_ks) == length(row_ks)) {
    warning('All rows have a cluster size of 1')
  }
    
  ######  Clustering  #####
  # For computational efficiency, only redo clustering on probes with >1 cluster
  ps <- row_ks[row_ks > 1]
  
  # Pre-fill ones matrix
  clusters <- matrix(1, nrow = nrow(betas), ncol = ncol(betas),
                     dimnames = dimnames(betas))
  size <- NA
  if(length(ps) > 0) {
    # Need drop=FALSE as it's a special case for 1 row
    bs <- betas[names(ps), ,drop=FALSE]
    if(printInfo) {cat("Clustering ", length(ps), " probes on ",
                       cores, " core(s)", sep="") }
    x <- parallel::mclapply(seq_along(ps), function(i) {
      Ckmeans.1d.dp::Ckmeans.1d.dp(bs[i, ], k=ps[i])[c("cluster", "centers",
                                                       "size", "withinss",
                                                       "tot.withinss")]
    }, mc.cores=cores)
    names(x) <- names(ps)
    if(printInfo) { cat(", done.\n") }
    
    myclusters <- t(sapply(x, "[[", "cluster"))
    colnames(myclusters) <- colnames(bs)
    # Stuff in values
    clusters[rownames(myclusters), ] <- myclusters
    size <- lapply(x, function(a) a[c("centers","size", "withinss", "tot.withinss")])
  }
  kClusters <- list(clusters = clusters,
                    size = size)
  class(kClusters) <- "kClusters"
  kClusters
}