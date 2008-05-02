\name{shrinkldaCMA-methods}
\docType{methods}
\alias{shrinkldaCMA-methods}
\alias{shrinkldaCMA,matrix,numeric,missing-method}
\alias{shrinkldaCMA,matrix,factor,missing-method}
\alias{shrinkldaCMA,data.frame,missing,formula-method}
\alias{shrinkldaCMA,ExpressionSet,character,missing-method}
\title{Shrinkage linear discriminant analysis}
\description{Linear Discriminant Analysis combined with the
             James-Stein-Shrinkage approach of Schaefer and Strimmer (2005)
             for the covariance matrix.\cr
             Currently still an experimental version.
             For \code{S4} method information, see \link{shrinkldaCMA-methods}}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{shrinkldaCMA}}.
}
\keyword{multivariate}

