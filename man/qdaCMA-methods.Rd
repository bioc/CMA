\name{qdaCMA-methods}
\docType{methods}
\alias{qdaCMA-methods}
\alias{qdaCMA,matrix,numeric,missing-method}
\alias{qdaCMA,matrix,factor,missing-method}
\alias{qdaCMA,data.frame,missing,formula-method}
\alias{qdaCMA,ExpressionSet,character,missing-method}
\title{Quadratic Discriminant Analysis}
\description{
 Performs a quadratic discriminant analysis under the assumption
of a multivariate normal distribution in each classes without restriction
concerning the covariance matrices. The function \code{qda} from
the package \code{MASS} is called for computation.
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{qdaCMA}}.}
\keyword{multivariate}