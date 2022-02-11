#' @name harmanresults
#' @title Harman results object
#' @description The S3 object returned after running \code{\link{harman}}.
#' @slot factors \code{A data.frame} of the \code{expt} and \code{batch}
#' vectors.
#' @slot parameters The harman runtime parameters. See \code{\link{harman}}
#' for details.
#' @slot stats Confidence intervals and the degree of correction for each
#' principal component.
#' @slot center The centering vector returned by \code{\link{prcomp}} with
#' \code{center=TRUE}.
#' @slot rotation The matrix of eigenvectors (by column) returned from
#' \code{\link{prcomp}}.
#' @slot original The original PC scores returned by \code{\link{prcomp}}.
#' @slot corrected The harman corrected PC scores.
#' @details \code{harmanresults} is the S3 object used to store the results from
#' \code{\link{harman}}.
#' This object may be presented to summary and data exploration functions such
#' as \code{\link{plot.harmanresults}}
#' and \code{\link{summary.harmanresults}} as well as the
#' \code{\link{reconstructData}} function which creates a corrected matrix of
#' data with the batch effect removed.
#' @seealso \code{\link{harman}}, \code{\link{reconstructData}},
#' \code{\link{pcaPlot}}, \code{\link{arrowPlot}}
#' @examples
#' ## HarmanResults
#' library(HarmanData)
#' data(OLF)
#' expt <- olf.info$Treatment
#' batch <- olf.info$Batch
#' olf.harman <- harman(as.matrix(olf.data), expt, batch)
#' plot(olf.harman)
#' summary(olf.harman)
#' pcaPlot(olf.harman, pc_x=2, pc_y=3)
#' pcaPlot(olf.harman, pc_x=2, pc_y=3, colBy='expt', pch=1)
#' olf.data.corrected <- reconstructData(olf.harman)
NULL


#' @title     Summarizing harman results.
#' @description Summary method for class \code{\link{harmanresults}}.
#' @param object An object of class \code{harmanresults}.
#' @param ... further parameters.
#' @return    Returns an object of class \code{summary.harmanresults}.
#' @seealso  \code{\link{harmanresults}}
#' @examples
#' library(HarmanData)
#' data(OLF)
#' expt <- olf.info$Treatment
#' batch <- olf.info$Batch
#' olf.harman <- harman(olf.data, expt, batch)
#' summary(olf.harman)
#' @export
summary.harmanresults <- function(object, ...) {
  # 1) % of the variance removed.
  # 2) Sequence of corrections from the 1st to last PC
  # 3) Confidence threshold
  # 4) Whether is was balanced or not
  # 5) Parameters
  # 6) Batch structure
  
  ans <- list()
  ans$totals <- c(length(levels(object$factors$expt)),
                  length(levels(object$factors$batch))
  )
  names(ans$totals) <- c('expt', 'batch')
  ans$factor_table <- table(expt=object$factors$expt,
                            batch=object$factors$batch)
  ans$parameters <- object$parameters
  ans$correction <- object$stats$correction
  names(ans$correction) <- object$stats$dimension
  class(ans) <- c("summary.harmanresults")
  ans
}
#' @title Printing Harmanresults summaries.
#' @param x an object of class \code{summary.harmanresults}, usually, a result
#' of a call to \code{summary.harmanresults}.
#' @param ... further parameters.
#' @description Print method for \code{summary.harmanresults}.
#' @return    Prints summary information from an object of class
#' \code{summary.harmanresults}.
#' @export
print.summary.harmanresults <- function(x, ...) {
  
  cat('Total factor levels:\n\n')
  print(x$totals)
  cat('\nExperiment x Batch Design:\n\n')
  print(x$factor_table)
  cat('\nSimulation parameters:\n')
  cat(x$parameters$numrepeats, ' simulations (with seed of ',
      x$parameters$randseed, '). ForceRand is ', x$parameters$forceRand, '.\n',
      sep='')
  cat('\nHarman results with confidence limit of ', x$parameters$limit, ':\n',
      sep='')
  print(x$correction)
  cat('\nTop batch-effected PCs:\n')
  top <- sort(x$correction)
  top <- top[top < 1]
  ntop <- min(5, length(top))
  if(ntop > 0) {
    print(top[seq_len(ntop)])
  } else {
    cat('None, all have no correction with this confidence limit\n')
  }
}


#' @title     Plot method for harman
#' @description Plot method for instances of \code{\link{harmanresults}}.
#' @param     x An instance of \code{harmanresults}.
#' @param     ... further plotting parameters.
#' @return    None
#' @seealso  \code{\link{harmanresults}} \code{\link{pcaPlot}}
#' @examples 
#' library(HarmanData)
#' data(OLF)
#' expt <- olf.info$Treatment
#' batch <- olf.info$Batch
#' olf.harman <- harman(olf.data, expt, batch)
#' plot(olf.harman)
#' @importFrom graphics par
#' @export
plot.harmanresults <- function(x, ...) {
  # set xlim and ylim  to be the same
  this_pc_x <- 1
  this_pc_y <- 2
  
  params <- list(...)
  
  if('pc_x' %in% names(params)) {
    this_pc_x <- params[['pc_x']]
  } else {
    if(length(params) >= 1 && is.numeric(params[[1]])) {
      this_pc_x <- params[[1]]
    }
  }
  if('pc_y' %in% names(params)) {
    this_pc_y <- params[['pc_y']]
  } else {
    if(length(params) >= 2 && is.numeric(params[[2]])) {
      this_pc_y <- params[[2]]
    }
  }
  
  xrange <- range(c(x$original[, this_pc_x], x$corrected[, this_pc_x]))
  yrange <- range(c(x$original[, this_pc_y], x$corrected[, this_pc_y]))
  old_mfrow <- graphics::par()$mfrow
  graphics::par(mfrow=c(1, 2))
  pcaPlot(x, this='original', main='Original', xlim=xrange, ylim=yrange,
          legend=TRUE, ...)
  pcaPlot(x, this='corrected', main='Corrected', xlim=xrange, ylim=yrange,
          legend=FALSE, ...)
  graphics::par(mfrow=old_mfrow)
}
