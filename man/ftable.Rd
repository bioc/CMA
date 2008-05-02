\name{ftable}
\alias{ftable,cloutput-method}
\title{Cross-tabulation of predicted and true class labels}
\description{
 An object of class \code{\link{cloutput}} contains (among others)
 the slot \code{y} and \code{yhat}. The former contains the true,
 the last the predicted class labels. Both are cross-tabulated in
 order to obtain a so-called confusion matrix. Counts out of the
 diagonal are misclassifications.}
\arguments{
  \item{x}{An object of class \code{\link{cloutput}}}
  \item{\dots}{Currently unused argument.}
}

\value{No return.}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{For more advanced evaluation: \code{\link{evaluation}}}
\keyword{multivariate}

