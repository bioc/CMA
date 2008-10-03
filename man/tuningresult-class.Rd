\name{tuningresult-class}
\docType{class}
\alias{tuningresult-class}
\alias{tuningresult}
\alias{show,tuningresult-method}
\title{"tuningresult"}
\description{Object returned by the function \code{\link{tune}}}
\section{Slots}{
	 \describe{
    \item{\code{hypergrid}:}{A \code{data.frame} representing the
                             grid of values that were tried and
                             evaluated. The number of columns equals
                             the number of tuned hyperparameters
                             and the number rows equals the number
                             of all possible combinations of the
                             discrete grids.}
    \item{\code{tuneres}:}{A list whose lengths equals the number
                           of different \code{learningsets} for which
                           tuning has been performed and whose
                           elements are numeric vectors with length
                           equal to the number of rows of \code{hypergrid} (s.above),
                           containing the misclassifcation
                           rate belonging to the respective
                           hyperparameter/hyperparameter combination.
                           In order to to get an overview
                           about the best hyperparmeter/hyperparameter
                           combination, use the convenience method
                           \code{\link{best}}}
    \item{\code{method}:}{Name of the classifier that has been tuned.}
    \item{\code{fold}:}{Number of cross-validation fold used for
                        tuning, s. argument of the same name in
                        \code{\link{tune}}}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(tuninresult-object)} for brief information.}
    \item{best}{Use \code{best(tuningresult-object)} to see which
                hyperparameter/hyperparameter combination has performed
                best in terms of the misclassification rate,
                s. \code{\link{best,tuningresult-method}}}
    \item{plot}{Use \code{plot(tuningresult-object, iter, which)}
                to display the performance of hyperparameter/hyperparameter
                combinations graphically, either as one-dimensional or
                as two-dimensional (contour) plot,
                s. \code{\link{plot,tuningresult-method}}}
	 }
}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de}

        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\code{\link{tune}}}
\keyword{multivariate}
