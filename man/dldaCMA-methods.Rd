\name{dldaCMA-methods}
\docType{methods}
\alias{dldaCMA-methods}
\alias{dldaCMA,matrix,numeric,missing-method}
\alias{dldaCMA,matrix,factor,missing-method}
\alias{dldaCMA,data.frame,missing,formula-method}
\alias{dldaCMA,ExpressionSet,character,missing-method}
\title{Diagonal Discriminant Analysis}
\description{
 Performs a diagonal discriminant analysis under the assumption
of a multivariate normal distribution in each classes (with equal,
diagonally structured) covariance matrices.
The method is also known under the name 'naive Bayes' classifier.
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{dldaCMA}}.
}
\keyword{multivariate}