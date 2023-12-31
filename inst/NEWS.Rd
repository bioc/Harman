\name{NEWS}
\title{News for Package \pkg{Harman}}
\encoding{UTF-8}

\section{Changes in version 1.23.1}{
  \itemize{
		\item Code syntax tweaks to remove warnings raised in later version of R.
		\item Increase the forceRand unit test from 10,000 to 30,000 iterations.
		This solves an occasional unit test failure due to sampling.
		}
}

\section{Changes in version 1.15.0}{
  \itemize{
		\item Previously, if no random seed is specified to the \code{harman}
		function, a seed of 0 was used. Now a random seed between 0 and 1e9 is
		picked at random.
		}
}

\section{Changes in version 1.2.0}{
  \itemize{
		\item Dynamically resizes legends depending upon the number of batches in
		\code{pcaPlot}
		\item A custom prcomp function to get the appropriate scores. The standard R
		\code{prcomp} function did not work in instances where the number of samples
		(matrix columns) was greater than the number of assays (matrix rows), so a
		special case is needed for less assays than samples. We need to use u'
		instead of v from the SVD. Presently, this is under developement, so an
		error is thrown if rows < cols.
		\item \code{prcompPlot} now has an argument for scaling, which defaults to
		FALSE. Previously, scaling was always TRUE. This new default makes
		prcompPlots plots agree with the 'original' plots of plotting
		\code{harmanresults} objects.
		\item Extensive updating of the vignette with a new comparison to ComBat
		from \pkg{sva}.
  	 }
}

\section{Version 1.0.2}{
  \itemize{
		\item First public version on Bioconductor.
	  }
}
