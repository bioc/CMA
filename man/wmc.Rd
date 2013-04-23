\name{wmc}
\alias{wmc}
\title{Tuning / Selection bias correction based on matrix of subsampling
fold errors}
\description{Perform tuning / selection bias correction for a matrix of subsampling fold errors. 
 }
\usage{
wmc(mcr.m,n.tr,n.ts,shrinkage=F)
}

\arguments{
  \item{mcr.m}{A matrix of resampling fold errors. Columns correspond
  the the fold errors of a single classifier.}
  \item{n.tr}{Number of observations in the resampling training sets.}
  \item{n.ts}{Number of observations in the resampling test sets.}
  \item{shrinkage}{A logical value indicating whether shrinkage (WMCS)
  shall be applied.}
  }

\details{The algorithm tries to avoid the additional computational costs
           of a nested cross validation by estimating the corrected
           misclassification rate of the best classifier by a weighted mean
           of all classifiers included in the subsampling approach.}
\value{A list containing the corrected misclassification rate, the index
           of the best method and a logical value indicating whether
           shrinkage has been applied.}

\references{Bernau Ch., Augustin, Th. and Boulesteix, A.-L. (2011): Correcting the optimally selected resampling-based error rate: A smooth analytical alternative to nested cross-validation. Department of Statistics: Technical Reports, Nr. 105.}

\author{Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}
        
       }

\seealso{\code{\link{weighted.mcr}},\code{\link{classification}},\code{\link{GeneSelection}}, \code{\link{tune}}, \code{\link{evaluation}},
}

\keyword{tuning bias, selection bias, corrected misclassification rate}