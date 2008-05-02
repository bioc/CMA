\name{learningsets-class}
\docType{class}
\alias{learningsets-class}
\alias{learningsets}
\alias{show,learningsets-method}
\title{"learningsets"}
\description{An object returned from \code{\link{GenerateLearningsets}}
             which is usually passed as arguments to \code{\link{GeneSelection}},
             \code{\link{tune}} and \code{\link{classification}}.}

\section{Slots}{
	 \describe{
    \item{\code{learnmatrix}:}{ A matrix of dimension \code{niter x
	ntrain}. Each row contains the indices of those observations
      representing the learningset for one iteration. If \code{method =
	CV}, zeros appear due to rounding issues.}
    \item{\code{method}:}{The method used to generate the \code{learnmatrix},
                          s.\code{\link{GenerateLearningsets}}}
    \item{\code{ntrain}:}{Number of observations in one learning set.If
      \code{method = CV}, this number is not attained for all
      iterations, due to rounding issues.}
    \item{\code{iter}:}{Number of iterations (different learningsets)
      that are stored in \code{learnmatrix}.}
  }
}
\section{Methods}{
 \item{show}{Use \code{show(learningsets-object)} for brief information.}
}



\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}




\seealso{\code{\link{GenerateLearningsets}}, \code{\link{GeneSelection}},
         \code{\link{tune}}, \code{\link{classification}}}

\keyword{multivariate}
