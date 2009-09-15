\name{summary}
\alias{summary,evaloutput-method}
\title{Summarize classifier evaluation}
\description{This method principally does nothing more than
applying the pre-implemented \code{summary()} function to the slot \code{score}
of an object of class \code{\link{evaloutput}}. One then obtains the usual
five-point-summary, consisting of minimum and maximum, lower and upper quartile
and the median. Additionally, the mean is also shown.}

\arguments{
  \item{object}{An object of class \code{\link{evaloutput}}.}
  \item{\dots}{Further arguments passed to the pre-implemented
          \code{summary} function.}
}

\value{No return.}
\note{That the results normally differ for different evaluation
      schemes ("iterationwise" or "observationwise").}
\author{Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}}
\seealso{\code{\link{evaluation}}, \code{\link{compare}}, \code{\link{obsinfo}}.}
\keyword{multivariate}
