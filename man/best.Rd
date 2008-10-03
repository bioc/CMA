\name{best}
\alias{best,tuningresult-method}
\alias{best}
\title{Show best hyperparameter settings}
\description{
  In this package hyperparameter tuning is performed by
  an inner cross-validation step for each \code{learningset}.
  A grid of values is tried and evaluated in terms of the
  misclassification rate, the results are saved in an object
  of class \code{\link{tuningresult}}. This method displays
  (separately for each \code{learningset}) the hyperparameter/
  hyperparameter combination that showed the best results. Note
  that this must not be unique; in this case, only one combination
  is displayed.
}
\usage{
best(object, ...)
}
\arguments{
  \item{object}{An object of class \code{\link{tuningresult}}.}
  \item{\dots}{Currently unused argument.}
}

\value{A list with elements equal to the number of different learningsets.
       Each element contains the best hyperparameter combination and
       the corresponding misclassification rate.}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de}

        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\code{tune}}
\keyword{multivariate}
