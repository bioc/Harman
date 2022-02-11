#' @title Discover clustered beta values
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
#' @param ks integer, the range of k's to consider for clustering (defaults to
#' 1:10). 
#' @param min_cluster_size integer, the minimum number of samples in a cluster
#' (defaults to 5).
#' @param min_cluster_dist numeric, the minimum beta difference in cluster
#' centroids (defaults to 0.1).
#' @param max_clusters integer, the maximum value of k returned (defaults to 4).
#' @param cores integer, specifies the number of cores to use for computation.
#' @param printInfo logical, whether to print information during computation or
#' not.
#' @return  A named vector containing the optimal value for k.
#' @details Betas values should be of type \code{double} with samples in
#' columns and betas in rows. The betas need to be bounded between 0 and 1.
#' The matrix is typically exported from a \code{\link[minfi]{GenomicRatioSet}},
#' \code{\link[minfi]{GenomicMethylSet}} or \code{\link[minfi]{MethylSet}}
#' object via the \code{getBeta} S4 accessor method.
#' @seealso \code{\link{kClusterMethylation}}, \code{\link{clusterStats}}
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
#' @importFrom Ckmeans.1d.dp Ckmeans.1d.dp
#' @importFrom parallel mclapply
#' @importFrom methods is
#' @export
discoverClusteredMethylation <- function(betas, ks=1:10, min_cluster_size=5,
                                         min_cluster_dist=0.1, max_clusters=4,
                                         cores=1, printInfo=FALSE) {
  
  ######  Internal function  #####
  testCluster <- function(b, ks, min_cluster_size, min_cluster_dist) {
    ckset <- list(Ckmeans.1d.dp::Ckmeans.1d.dp(b, ks[1]),
                  Ckmeans.1d.dp::Ckmeans.1d.dp(b, ks[2]))
    names(ckset) <- ks
    # Are there at least min_cluster_size number of samples in a cluster
    has_min_n <- sapply(ckset, function(x) min(x$size) >= min_cluster_size)
    # If more than one cluster, the min gap must be at least min_cluster_dist
    has_min_gap <- sapply(ckset, function(x) {
      if(length(x$centers) > 1) {
        min(diff(x$centers)) >= min_cluster_dist
      } else {
        TRUE
      }
    })
    num_clust <- 1
    # Only run if something qualifies
    if(sum(has_min_n & has_min_gap) > 0) {
      # Very rarely returns -Inf, so build edge case for this
      num_clust <- max(1, max(ks[has_min_n & has_min_gap]))
    }
    num_clust
  }
  
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
  
  if(ncol(betas) <= min_cluster_size) {
    stop("Minimum cluster size must be less than the number of beta value columns")
  }
  
  if(min_cluster_dist >= 1 | min_cluster_dist <= 0) {
    stop("Minimum cluster distance must be between 0 and 1.")
  }
  
  ######  Clustering  #####
  # Ckmeans.1d.dp will throw a warning to increase k if max k is used.
  # Suppress that warning
  if(printInfo) { cat("First clustering") }
  suppressWarnings(
    unbiased_clusters <- parallel::mclapply(seq_len(nrow(betas)), function(p) {
      myclust <- Ckmeans.1d.dp::Ckmeans.1d.dp(betas[p, ], k=ks)
      list(num_clusters=length(myclust$centers),
           BIC=myclust$BIC)
    }, mc.cores = cores)
  )
  names(unbiased_clusters) <- rownames(betas)
  
  # Cycle over clusters and find the best two options by BIC
  ktext_to_int <- ks
  names(ktext_to_int) <- names(unbiased_clusters[[1]]$BIC)
  if(printInfo) { cat(", now find top two k") }
  # get the two top ks by BIC
  num_clusters <- sapply(unbiased_clusters, function(x) x$num_clusters)
  top_ks <- lapply(unbiased_clusters[num_clusters > 1], function(x) {
    bic <- x$BIC
    bic <- bic[!is.nan(bic)]
    bic_diffs <- diff(c(0, bic))
    highest_bic_diffs <- sort(bic_diffs, decreasing = TRUE)[1:2]
    as.integer(ktext_to_int[names(highest_bic_diffs)])
  })
  rm(ktext_to_int)
  
  # Recluster for multiple cluster probes and specify distance and number cutoff
  if(printInfo) {cat(", done. Reclustering") }
  bs_small <- betas[names(top_ks), ]
  best_k <- parallel::mclapply(seq_along(top_ks), function(i) {
    testCluster(bs_small[i, ], top_ks[[i]],
                min_cluster_size, min_cluster_dist)
  }, mc.cores=cores)
  names(best_k) <- names(top_ks)
  best_k <- do.call("c", best_k)
  if(printInfo) { cat(", done.\n") }
  
  num_clusters[names(best_k)] <- best_k
  # Reduce maximum k to max_clusters
  num_clusters[num_clusters > max_clusters] <- max_clusters
  num_clusters
}