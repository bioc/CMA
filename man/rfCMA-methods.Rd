\name{rfCMA-methods}
\docType{methods}
\alias{rfCMA-methods}
\alias{rfCMA,matrix,numeric,missing-method}
\alias{rfCMA,matrix,factor,missing-method}
\alias{rfCMA,data.frame,missing,formula-method}
\alias{rfCMA,ExpressionSet,character,missing-method}
\title{Classification based on Random Forests}
\description{Random Forests were proposed by Breiman (2001)
             and are implemented in the package \code{randomForest}.

             In this package, they can as well be used to rank variables
             according to their importance, s. \code{GeneSelection}.}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For references, further argument and output information, consult
\code{\link{rfCMA}}}
\keyword{multivariate}