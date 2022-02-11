library(Harman)
library(HarmanData)
library(RUnit)


#####  Global Data  #####

data(episcope)


#####  Tests  #####

test.discoverClusteredMethylation <- function() {
  betas <- matrix(runif(2600), nrow = 26, ncol=100)
  rownames(betas) <- letters
  colnames(betas) <- paste("cg", 1:ncol(betas), sep="")
  betas_na <- betas
  betas_range_low <- betas
  betas_range_high <- betas
  betas_na[1,1] <- NA
  betas_range_low[1,1] <- -0.1
  betas_range_high[26,5] <- 1.1
  
  res <- discoverClusteredMethylation(betas, max_clusters = 2)
  checkTrue(max(res) >= 2)
  checkException(discoverClusteredMethylation(betas, min_cluster_size = 100))
  checkException(discoverClusteredMethylation(betas, min_cluster_dist = 1))
  checkException(discoverClusteredMethylation(betas_na))
  checkException(discoverClusteredMethylation(betas_range_low))
  checkException(discoverClusteredMethylation(betas_range_high))
}


test.clusterStats <- function() {
  m <- matrix()
  not_m <- data.frame()
  kclust = list()
  class(kclust) <- "kClusters"
  not_kclust = list()
  checkException(clusterStats(not_m, m, kclust))
  checkException(clusterStats(m, not_m, kclust))
  checkException(clusterStats(m, m, not_kclust))
}


test.kClusterMethylation <- function() {
  betas <- matrix(runif(130), nrow = 26, ncol=5)
  rownames(betas) <- letters
  colnames(betas) <- paste("cg", 1:ncol(betas), sep="")
  betas_na <- betas
  betas_range_low <- betas
  betas_range_high <- betas
  betas_na[1,1] <- NA
  betas_range_low[1,1] <- -0.1
  betas_range_high[26,5] <- 1.1
  ks <- rep(1, nrow(betas))
  names(ks) <- rownames(betas)
  ks[5] <- 2
  
  checkEquals(kClusterMethylation(betas[1, , drop=FALSE],
                                  row_ks = ks[1])$size, NA)
  checkException(kClusterMethylation(betas_na, row_ks = ks))
  checkException(kClusterMethylation(betas_range_low, row_ks = ks))
  checkException(kClusterMethylation(betas_range_high, row_ks = ks))

  betas_dummy <- matrix(c(0.1, 0.2, 0.3, 0.7, 0.8, 0.9), nrow=1)
  rownames(betas_dummy) <- letters[1]
  res <- kClusterMethylation(betas_dummy, row_ks = c("a" = 2))
  checkIdentical(res$size$a$size, c(3,3))
  checkEquals(res$size$a$centers, c(0.2,0.8))
}


test.cluster.math <- function() {
  bad_batches <- c(1, 5, 9, 17, 25)
  is_bad_sample <- episcope$pd$array_num %in% bad_batches
  
  mydiscoverClust <- discoverClusteredMethylation(episcope$original[, !is_bad_sample], printInfo = TRUE)
  mykClust <- kClusterMethylation(betas=episcope$original, row_ks=mydiscoverClust)
  
  res <- clusterStats(pre_betas=episcope$original,
                      post_betas=episcope$harman,
                      kClusters = mykClust)
  checkEquals(episcope$ref_md$meandiffs_harman, res$meandiffs)
  checkEquals(episcope$ref_lvr$var_ratio_harman, res$log2_var_ratio)
}

