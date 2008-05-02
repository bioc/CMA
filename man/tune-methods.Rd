\name{tune-methods}
\docType{methods}
\alias{tune-methods}
\alias{tune,matrix,numeric,missing-method}
\alias{tune,matrix,factor,missing-method}
\alias{tune,data.frame,missing,formula-method}
\alias{tune,ExpressionSet,character,missing-method}
\title{Hyperparameter tuning for classifiers}
\description{
  Performs hyperparameter tuning for the following signatures:
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{tune}}.}
\keyword{multivariate}
