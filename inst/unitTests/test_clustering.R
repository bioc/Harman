library(Harman)
library(HarmanData)
library(RUnit)


#####  Global Data  #####

data(EpiSCOPE_sample)


#####  Tests  #####

test.clusters <- function() {
  bad_batches <- c(1, 5, 9, 17, 25)
  is_bad_sample <- episcope$pd$array_num %in% bad_batches
  
  mydiscoverClust <- discoverClusteredMethylation(episcope$original[, !is_bad_sample], verbose = TRUE)
  mykClust = kClusterBetas(episcope$original, k=mydiscoverClust)
  
  res = ClusterStatistics(pre_betas=episcope$original,
                          post_betas=episcope$harman,
                          clusters = mykClust)
  checkEquals(episcope$ref_md$meandiffs_harman, res$meandiffs)
  checkEquals(episcope$ref_lvr$var_ratio_harman, res$log2_var_ratio)
}
