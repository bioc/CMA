\name{pnnCMA-methods}
\docType{methods}
\alias{pnnCMA-methods}
\alias{pnnCMA,matrix,numeric,missing-method}
\alias{pnnCMA,matrix,factor,missing-method}
\alias{pnnCMA,data.frame,missing,formula-method}
\alias{pnnCMA,ExpressionSet,character,missing-method}
\title{Probabilistic Neural Networks}
\description{
             Probabilistic Neural Networks is the term Specht (1990) used
             for a Gaussian kernel estimator for the conditional class
             densities.
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For references, further argument and output information, consult
\code{\link{pnnCMA}}.
}
\keyword{multivariate}