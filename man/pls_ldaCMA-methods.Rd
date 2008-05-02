\name{pls_ldaCMA-methods}
\docType{methods}
\alias{pls_ldaCMA-methods}
\alias{pls_ldaCMA,matrix,numeric,missing-method}
\alias{pls_ldaCMA,matrix,factor,missing-method}
\alias{pls_ldaCMA,data.frame,missing,formula-method}
\alias{pls_ldaCMA,ExpressionSet,character,missing-method}
\title{Partial Least Squares combined with Linear Discriminant Analysis}
\description{
 -This method constructs a classifier that extracts
  Partial Least Squares components that are plugged into
  Linear Discriminant Analysis.
  The Partial Least Squares components are computed by the package
  \code{plsgenomics}.
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{pls_ldaCMA}.}}
\keyword{multivariate}
