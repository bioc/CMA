\name{pls_lrCMA-methods}
\docType{methods}
\alias{pls_lrCMA-methods}
\alias{pls_lrCMA,matrix,numeric,missing-method}
\alias{pls_lrCMA,matrix,factor,missing-method}
\alias{pls_lrCMA,data.frame,missing,formula-method}
\alias{pls_lrCMA,ExpressionSet,character,missing-method}
\title{Partial Least Squares followed by logistic regression}
\description{
 This method constructs a classifier that extracts
  Partial Least Squares components that form the the covariates
  in a binary logistic regression model.
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
\code{\link{pls_lrCMA}}
}
\keyword{multivariate}