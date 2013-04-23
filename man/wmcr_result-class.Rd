\name{wmcr.result-class}
\docType{class}
\alias{wmcr.result-class}
\alias{wmcr.result}
\alias{show,wmcr.result-method}

\title{"wmcr.result"}
\description{Object returned by function \code{weighted.mcr}.}
\section{Slots}{
	 \describe{
    \item{\code{corrected.mcr}:}{The corrected misclassification rate for
      the best method.}
    \item{\code{best.method}:}{The method which performed best in the
      subsampling approach.}
    \item{\code{mcrs}:}{Misclassification rates of all classifiers used
      in the subsampling approach.}
   \item{\code{weights}:}{The weights used for the different classifiers
   in the correction method.}
    \item{\code{cov}:}{Estimated covariance matrix for the
    misclassification rates of the different classifiers.}
    \item{\code{uncorrected.mcr}}{The uncorrected misclassification rate of the
      best method.}
    \item{\code{ranges}}{Minimum and maximal mean misclassification rates as
    well as the theoretical bound for nested cross validation (averaging
    over foldwise minima or maxima respectively).}
  \item{\code{mcr.m}}{matrix of resampling fold errors, columns
    correspond to the fold errors of a single classifier}
   \item{\code{shrinkage}}{a logical value indicating whether shrinkage
     (WMCS) has been aplied.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(wmcr.result-object)} for brief information}
   
}}

\author{Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}}

\seealso{\code{weighted.mcr}}


\keyword{tuning bias, selection bias, corrected misclassification rate}