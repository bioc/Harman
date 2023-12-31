#' @title Harman batch correction method
#' @description Harman is a PCA and constrained optimisation based technique 
#' that maximises the removal of batch effects from datasets, with the
#' constraint that the probability of overcorrection (i.e. removing genuine
#' biological signal along with batch noise) is kept to a fraction which is set
#' by the end-user (Oytam et al, 2016;
#' \url{http://dx.doi.org/10.1186/s12859-016-1212-5}).
#' 
#' Harman expects unbounded data, so for example, with HumanMethylation450
#' arrays do not use the Beta statistic (with values constrained between 0 and
#' 1), instead use the logit transformed M-values.
#' @param datamatrix matrix or data.frame, the data values to correct with
#' samples in columns and data values in rows. Internally, a data.frame will be
#' coerced to a matrix. Matrices need to be of type \code{integer} or
#' \code{double}.
#' @param expt vector or factor with the experimental variable of interest
#' (variance to be kept).
#' @param batch vector or factor with the batch variable (variance to be
#' removed).
#' @param limit numeric, confidence limit. Indicates the limit of confidence in
#' which to stop removing a batch effect. Must be between \code{0} and \code{1}.
#' @param numrepeats integer, the number of repeats in which to run the
#' simulated batch mean distribution estimator using the random selection
#' algorithm. (N.B. 32 bit Windows versions may have an upper limit of 300000
#' before catastrophic failure)
#' @param randseed integer, the seed for random number generation.
#' @param forceRand logical, to enforce Harman to use a random selection
#' algorithm to compute corrections. Force the simulated mean code to
#' use random selection of scores to create the simulated batch mean (rather
#' than full explicit calculation from all permutations).
#' @param printInfo logical, whether to print information during computation or
#' not.
#' @return  A \code{harmanresults} S3 object:
#' \describe{
#'   \item{factors}{\code{A data.frame} of the \code{expt} and \code{batch}
#'   vectors}
#'   \item{parameters}{The harman runtime parameters. See \code{\link{harman}}
#'   for details}
#'   \item{stats}{Confidence intervals and the degree of correction for each
#'   principal component}
#'   \item{center}{The centering vector returned by \code{\link{prcomp}} with
#'   \code{center=TRUE}}
#'   \item{rotation}{The matrix of eigenvectors (by column) returned from
#'   \code{\link{prcomp}}}
#'   \item{original}{The original PC scores returned by \code{\link{prcomp}}}
#'   \item{corrected}{The harman corrected PC scores}
#' }
#' @details The \code{datamatrix} needs to be of type \code{integer} or
#' \code{numeric}, or alternatively a data.frame that can be coerced into one
#' using \code{\link{as.matrix}}. The matrix is to be constructed with data
#' values (typically microarray probes or sequencing counts) in rows and samples
#' in columns, much like the `assayData` slot in the canonical Bioconductor
#' \code{eSet} object, or any object which inherits from it. The data should
#' have normalisation and any other global adjustment for noise reduction
#' (such as background correction) applied prior to using Harman.
#' 
#' For converge, the number of simulations, \code{numrepeats} parameter should
#' probably should be at least 100,000. The underlying principle of Harman rests
#' upon PCA, which is a parametric technique. This implies Harman should be
#' optimal when the data is normally distributed. However, PCA is known to be
#' rather robust to very non-normal data.
#' 
#' The output \code{harmanresults} object may be presented to summary and data
#' exploration functions such as \code{\link{plot.harmanresults}} and
#' \code{\link{summary.harmanresults}} as well as the
#' \code{\link{reconstructData}} function which creates a corrected matrix of
#' data with the batch effect removed.
#' @seealso \code{\link{reconstructData}},
#' \code{\link{pcaPlot}}, \code{\link{arrowPlot}}
#' @references Oytam et al (2016) BMC Bioinformatics 17:1.
#' DOI: 10.1186/s12859-016-1212-5
#' @examples
#' library(HarmanData)
#' data(OLF)
#' expt <- olf.info$Treatment
#' batch <- olf.info$Batch
#' olf.harman <- harman(olf.data, expt, batch)
#' plot(olf.harman)
#' olf.data.corrected <- reconstructData(olf.harman)
#'
#' ## Reading from a csv file
#' datafile <- system.file("extdata", "NPM_data_first_1000_rows.csv.gz",
#' package="Harman")
#' infofile <- system.file("extdata", "NPM_info.csv.gz", package="Harman")
#' datamatrix <- read.table(datafile, header=TRUE, sep=",", row.names="probeID")
#' batches <- read.table(infofile, header=TRUE, sep=",", row.names="Sample")
#' res <- harman(datamatrix, expt=batches$Treatment, batch=batches$Batch)
#' arrowPlot(res, 1, 3)
#' @importFrom methods is
#' @importFrom stats runif
#' @export
harman <- function(datamatrix, expt, batch, limit=0.95, numrepeats=100000L,
                   randseed, forceRand=FALSE, printInfo=FALSE) {
  
  ######  Coerce a data.frame to a matrix  ##### 
  if(is.data.frame(datamatrix)) {
    datamatrix <- as.matrix(datamatrix)
  }
  
  ######  Sanity checks  #####
  if(!methods::is(datamatrix, "matrix")) {
    stop(paste("Require 'datamatrix' as a matrix or data.frame, not class \'",
               class(datamatrix), "\'.", sep=""))
  }
  
  if(!typeof(datamatrix) %in% c("integer", "double")) {
    stop(paste("'datamatrix' is type \'", typeof(datamatrix), "\',
               needs to be of type \'integer\' or \'double\'.", sep=""))
  }
  
  if(!is.vector(expt) && !is.factor(expt)) {
    stop(paste("Require 'expt' to be a vector or factor, not class \'",
               class(expt), "\'.", sep=""))
  }
  if(!is.vector(batch) && !is.factor(batch)) {
    stop(paste("Require 'batch' to be a vector or factor, not class \'",
               class(batch), "\'.", sep=""))
  }
  if(is.vector(expt)) {
    if(sum(is.na(expt)) > 0 ||
       sum(is.nan(expt)) > 0 ||
       sum(is.null(expt)) > 0) {
      stop("Cannot have NA, NaN or NULL as 'expt' levels.")
    }
  }
  if(is.vector(batch)) {
    if(sum(is.na(batch)) > 0 ||
       sum(is.nan(batch)) > 0 ||
       sum(is.null(batch)) > 0) {
      stop("Cannot have NA, NaN or NULL as 'batch' levels.")
    }
  }
  if(length(expt) != length(batch)) stop("'expt' and 'batch' vectors not the
                                         same length.")
  if(!is.numeric(limit) || limit < 0 || limit > 1) {
    stop(paste("'limit' needs to be a number between 0 and 1, not \"", limit,
               "\".", sep=""))
  }
  if(!is.numeric(numrepeats)) stop("'numrepeats' needs to be numeric.")
  
  if (!missing(randseed)) {
    if(!is.numeric(randseed)) stop("'randseed' needs to be numeric.")
  } else {
    randseed <- stats::runif(1, 0, 1e9)
  }
  
  strict <- FALSE
  #  Sanity checking to see if the expt vector length is equal to the
  #  number of matrix columns
  if(length(expt) != ncol(datamatrix)) {
    msg <- "'expt' vector not equal to the number of datamatrix columns."
    if(strict == FALSE) {
      warning(msg)
    } else {
      stop(msg)
    }
  }
  
  # Coerce expt and batch to factors
  expt <- factor(expt)
  batch <- factor(batch)
  
  if(length(levels(expt)) < 2 || length(levels(batch)) < 2) {
    stop("Require more than one experimental factor and/or batch for experiment
         structure")
  }
  
  #####  PCA  #####  
  
  # Don't shift the original data into the .RunHarman function as we just need
  # the PCs to kick it off.
  if(printInfo == TRUE) cat('Performing PCA... ')
  pca <- harmanScores(datamatrix)
  pc_data_scores <- pca$scores
  # Try and free up RAM, but keep the sample names first.
  sample_names <- dimnames(datamatrix)[[2]]
  rm(datamatrix)
  gc()
  if(printInfo == TRUE) cat('done.\n')
  
  #####  Call Rcpp wrapper function  #####
  
  # Form group construct by converting all expt and batch names to an integer
  group <- as.matrix(data.frame(expt=as.integer(expt), batch=as.integer(batch)))
  rownames(group) <- sample_names

  
  if(printInfo == TRUE) cat('Now calling the Rcpp layer.\n')
  res <- .callHarman(pc_data_scores, group, limit, numrepeats, randseed,
                     forceRand, printInfo)
  
    
  #####  Form S3 object  #####
  
  parameters <- list(limit=limit, numrepeats=numrepeats, randseed=randseed,
                     forceRand=forceRand)
  factors <- data.frame(expt=expt, batch=batch)
  rownames(factors) <- sample_names
  factors$expt.numeric <- group[, 'expt']
  factors$batch.numeric <- group[, 'batch']
  dim_names <- paste('PC',seq_len(length(res[["confidence_vector"]])), sep='')
  stats <- data.frame(dimension=dim_names,
                      confidence=res[["confidence_vector"]],
                      correction=res[["correction_vector"]])
  
  # Corrected scores are returned as a numeric array
  # Need to coerce them into a matrix
  corrected <- matrix(res[["corrected_scores"]],
                      nrow=nrow(pc_data_scores),
                      ncol=ncol(pc_data_scores),
                      dimnames=dimnames(pc_data_scores))

  results <- list(factors=factors,
                  parameters=parameters,
                  stats=stats,
                  center=pca$center,
                  rotation=pca$rotation,
                  original=pc_data_scores,
                  corrected=corrected)
  
  # Define the S3 results class
  results <- structure(results, class = "harmanresults")
  results
}