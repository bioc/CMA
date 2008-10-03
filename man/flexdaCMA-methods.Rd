\name{flexdaCMA-methods}
\docType{methods}
\alias{flexdaCMA-methods}
\alias{flexdaCMA,matrix,numeric,missing-method}
\alias{flexdaCMA,matrix,factor,missing-method}
\alias{flexdaCMA,data.frame,missing,formula-method}
\alias{flexdaCMA,ExpressionSet,character,missing-method}
\title{Flexible Discriminant Analysis}
\description{
  This method is experimental.

 It is easy to show that, after appropriate scaling of the predictor matrix \code{X},
 Fisher's Linear Discriminant Analysis is equivalent to Discriminant Analysis
 in the space of the  fitted values from the linear regression of the
 \code{nlearn x K} indicator matrix of the class labels on \code{X}.
 This gives rise to 'nonlinear discrimant analysis' methods that expand
 \code{X} in a suitable, more flexible basis. In order to avoid overfitting,
 penalization is used. In the implemented version, the linear model is replaced
 by a generalized additive one, using the package \code{mgcv}.
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{flexdaCMA}}.
}
\keyword{multivariate}