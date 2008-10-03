\name{scdaCMA-methods}
\docType{methods}
\alias{scdaCMA-methods}
\alias{scdaCMA,matrix,numeric,missing-method}
\alias{scdaCMA,matrix,factor,missing-method}
\alias{scdaCMA,data.frame,missing,formula-method}
\alias{scdaCMA,ExpressionSet,character,missing-method}
\title{Shrunken Centroids Discriminant Analysis}
\description{
  The nearest shrunken centroid classification algorithm is
  detailly described in Tibshirani et al. (2002).

  It is widely known under the name PAM (prediction analysis for microarrays),
  which can also be found in the package \code{pamr}.
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For references, further argument and output information, consult
\code{\link{scdaCMA}}.
}
\keyword{multivariate}