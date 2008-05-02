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
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\code{\link{evaluation}}}
\keyword{multivariate}
