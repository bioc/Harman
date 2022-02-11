#' @title Compute LVR and meandiff statistics for beta values after batch
#' correction
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
#' @param pre_betas a matrix of methylation beta values \bold{prior to}
#' correction.
#' @param post_betas a matrix of methylation beta values \bold{after}
#' correction.
#' @param kClusters a kClusters S3 object
#' @return  A \code{data.frame} containing clustering stats.
#' @details Betas values should be of type \code{double} with samples in
#' columns and betas in rows. The betas need to be bounded between 0 and 1.
#' The matrix is typically exported from a \code{\link[minfi]{GenomicRatioSet}},
#' \code{\link[minfi]{GenomicMethylSet}} or \code{\link[minfi]{MethylSet}}
#' object via the \code{getBeta} S4 accessor method.
#' @seealso \code{\link{kClusterMethylation}},
#' \code{\link{discoverClusteredMethylation}}
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
#' @importFrom matrixStats rowVars
#' @export
clusterStats <- function(pre_betas, post_betas, kClusters) {
  
  ######  Sanity checks  #####
  if(!methods::is(pre_betas, "matrix")) {
    stop(paste("Require 'pre_betas' as class matrix, not class \'",
               class(pre_betas), "\'.", sep=""))
  }
  if(!methods::is(post_betas, "matrix")) {
    stop(paste("Require 'post_betas' as class matrix, not class \'",
               class(post_betas), "\'.", sep=""))
  }
  if(!methods::is(kClusters, "kClusters")) {
    stop(paste("Require 'kClusters' as class kClusters, not class \'",
               class(kClusters), "\'.", sep=""))
  }
  
  ######  Internal functions  #####
  sumSquares <- function(x) {
    x <- x[!is.na(x)]
    # the 'classic' approach
    sum((x - mean(x))^2)
  }
  
  computeStats <- function(bs) {
    sample_n = ncol(bs) - 1
    x <- lapply(rownames(bs), function(p) {
      df <- as.data.frame(table(pre_clusters[p, ]))
      df <- cbind(p, df)
      names(df) <- c("probe_id", "cluster_id", "size")
      df$centers <- tapply(bs[p, ], INDEX = pre_clusters[p, ], mean)
      df$withinss <- tapply(bs[p, ], INDEX = pre_clusters[p, ], sumSquares)
      df
    })
    names(x) <- rownames(bs)
    list("clust" = do.call("rbind", x),
         "tot_ssq" = sapply(x, function(p) sum(p$withinss)))
  }
  
  ######  Compute stats  #####
  # Cluster-wise stats
  pre_clusters <- kClusters$clusters
  pre_ssq <- computeStats(pre_betas)
  post_ssq <- computeStats(post_betas)
  
  # Aggregate into probe-wise stats
  df <- data.frame("probe_id"=rownames(pre_clusters),
                   "num_clusters"=apply(pre_clusters, 1, function(x) length(unique(x))),
                   "pre_total_withinss"=pre_ssq$tot_ssq,
                   "post_total_withinss"=post_ssq$tot_ssq,
                   "pre_withinvar"=pre_ssq$tot_ssq / (ncol(pre_clusters) - 1),
                   "post_withinvar"=post_ssq$tot_ssq / (ncol(pre_clusters) - 1))
  df$pre_var <- matrixStats::rowVars(pre_betas)
  df$post_var <- matrixStats::rowVars(post_betas)
  df$log2_var_ratio <- log2(df$post_withinvar / df$pre_withinvar)
  df$meandiffs=rowMeans(abs(post_betas - pre_betas))
  df
}