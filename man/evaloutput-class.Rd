\name{evaloutput-class}
\docType{class}
\alias{evaloutput-class}
\alias{evaloutput}
\alias{obsinfo,evaloutput-method}
\alias{show,evaloutput-method}

\title{"evaloutput"}
\description{Object returned by the method \code{\link{evaluation}}.}
\section{Slots}{
	 \describe{
    \item{\code{score}:}{A numeric vector of performance scores whose length depends
                         on \code{"scheme"}, s.below. It equals the number of
                         iterations (number of different datasets) if
                         \code{"scheme = iterationwise"} and the number
                         of all observations in the complete dataset otherwise.
                         As not necessarily all observation must be predicted
                         at least one time, \code{score} can also contain \code{NA}s
                         for those observations not classified at all.}
    \item{\code{measure}:}{performance measure used, s. \code{\link{evaluation}}.}
    \item{\code{scheme}:}{scheme used, s. \code{\link{evaluation}}}
    \item{\code{method}:}{name of the classifier that has been evaluated.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(evaloutput-object)} for brief information.}
    \item{summary}{Use \code{summary(evaloutput-object)} to apply the
                   classic \code{summary()} function to the slot \code{score},
                   s. \code{\link{summary,evaloutput-method}}}
    \item{boxplot}{Use \code{boxplot(evaloutput-object)} to display
                   a boxplot of the slot \code{score}, s. \code{\link{boxplot,evaloutput-method}}.}
    \item{obsinfo}{Use \code{obsinfo(evaloutput-object, threshold)} to display all
                   observations consistenly correctly or incorrectly classified
                   (depending on the value of the argument \code{threshold}), s.
                   \code{\link{obsinfo}}.}
	 }
}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\code{\link{evaluation}}}
\keyword{multivariate}
