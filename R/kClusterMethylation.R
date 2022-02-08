#' @title Cluster beta values with a set value for k
#' @description TODO
#' 
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
#' library(Harman)
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
  
  ######  Clustering  #####
  # For computational efficiency, only redo clustering on probes with >1 cluster
  ps <- row_ks[row_ks > 1]
  bs <- betas[names(ps), ]
  if(printInfo) {cat("Clustering ", length(ps), " probes on ",
                   cores, " core(s)", sep="") }
  x <- parallel::mclapply(1:length(ps), function(i) {
    Ckmeans.1d.dp::Ckmeans.1d.dp(bs[i, ], k=ps[i])[c("cluster", "centers",
                                                     "size", "withinss",
                                                     "tot.withinss")]
  }, mc.cores=cores)
  names(x) <- names(ps)
  if(printInfo) { cat(", done.\n") }
  
  myclusters <- t(sapply(x, "[[", "cluster"))
  colnames(myclusters) <- colnames(bs)
  #stopifnot(identical(rownames(myclusters), rownames(bs)))
  # Pre-fill ones matrix
  clusters <- matrix(1, nrow = nrow(betas), ncol = ncol(betas),
                     dimnames = dimnames(betas))
  # Stuff in values
  clusters[rownames(myclusters), ] <- myclusters
  size <- lapply(x, function(a) a[c("centers","size", "withinss", "tot.withinss")])
  kClusters <- list(clusters = clusters,
                    size = size)
  class(kClusters) <- "kClusters"
  kClusters
}