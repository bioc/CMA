\name{knnCMA-methods}
\docType{methods}
\alias{knnCMA-methods}
\alias{knnCMA,matrix,numeric,missing-method}
\alias{knnCMA,matrix,factor,missing-method}
\alias{knnCMA,data.frame,missing,formula-method}
\alias{knnCMA,ExpressionSet,character,missing-method}
\title{Nearest Neighbours}
\description{
 Ordinary \code{k} nearest neighbours algorithm from the
very fast implementation in the package \code{class}
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{knnCMA}}.
}
\keyword{multivariate}