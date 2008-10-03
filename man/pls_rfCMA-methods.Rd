\name{pls_rfCMA-methods}
\docType{methods}
\alias{pls_rfCMA-methods}
\alias{pls_rfCMA,matrix,numeric,missing-method}
\alias{pls_rfCMA,matrix,factor,missing-method}
\alias{pls_rfCMA,data.frame,missing,formula-method}
\alias{pls_rfCMA,ExpressionSet,character,missing-method}
\title{Partial Least Squares followed by random forests}
\description{
  This method constructs a classifier that extracts
  Partial Least Squares components used to generate Random Forests, s. \code{\link{rfCMA}}.
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
\code{\link{pls_rfCMA}}.
}
\keyword{multivariate}