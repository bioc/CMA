\name{genesel-class}
\docType{class}
\alias{genesel-class}
\alias{genesel}
\alias{show,genesel-method}

\title{"genesel"}
\description{Object returned from a call to \code{\link{GeneSelection}}}
\section{Slots}{
	 \describe{
    \item{\code{rankings}:}{A list of matrices.
                            For the two-class case and the multi-class case
                            where a genuine multi-class method has been used
                            for variable selection, the length of the list is one.
                            Otherwise, it is named according to the different
                            binary scenarios (e.g. \code{1 vs 3}). Each list
                            element is a \code{matrix} with rows corresponding  to
                            iterations (different \code{learningsets}) and columns
                            to variables.
                            Each row thus contains an index vector representing the order of the
                            variables with respect to their variable importance
                            (s. slot \code{importance})}
    \item{\code{importance}:}{A list of matrices, with the same structure as
                              described for the slot \code{rankings}.
                              Each row of these matrices are ordered according to
                              \code{rankings} and contain the variable importance
                              measure (absolute value of test statistic or regression
                              coefficient).}
    \item{\code{method}:}{Name of the method used for variable selection, s. \code{\link{GeneSelection}}.}
    \item{\code{scheme}:}{The scheme used in the case of a non-binary
      response, one of \code{"pairwise"}, \code{"one-vs-all"} or \code{"multiclass"}}.
                
                
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(genesel-object)} for brief information}
    \item{toplist}{Use \code{toplist(genesel-object, k=10, iter = 1)} to display
                   the top first 10 variables and their variable importance
                   for the first iteration (first learningset),
                    s.\code{\link{toplist}}.}
    \item{plot}{Use \code{plot(genesel-object, k=10, iter=1)} to display
                a barplot of the variable importance of the top first 10
                variables, s. \code{\link{plot,genesel-method}}}
	 }
}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\code{\link{GeneSelection}}}
\keyword{multivariate}
