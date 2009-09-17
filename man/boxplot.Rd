\name{boxplot}
\alias{boxplot,evaloutput-method}
\title{Make a boxplot of the classifier evaluation}
\description{
  This method displays the slot \code{scores} of performance
  scores of an object of class \code{\link{evaloutput}}.
}
\arguments{
  \item{x}{An object of class \code{\link{evaloutput}}.}
  \item{\dots}{Further graphical parameters passed to the
               classical \code{boxplot} function.}
}
\value{The only return is a boxplot.}

\references{
Slawski, M. Daumer, M.  Boulesteix, A.-L. (2008)
CMA - A comprehensive Bioconductor package for supervised classification with high dimensional data.
\emph{BMC Bioinformatics 9: 439} }

\author{Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}}

\seealso{\code{\link{evaluation}}}
\keyword{multivariate}
