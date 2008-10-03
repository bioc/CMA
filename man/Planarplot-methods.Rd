\name{Planarplot-methods}
\docType{methods}
\alias{Planarplot-methods}
\alias{Planarplot,matrix,numeric,missing-method}
\alias{Planarplot,matrix,factor,missing-method}
\alias{Planarplot,data.frame,missing,formula-method}
\alias{Planarplot,ExpressionSet,character,missing-method}
\title{Visualize Separability of different classes}
\description{Given two variables, the methods trains a classifier
  (argument \code{classifier}) based on these two variables
  and plots the resulting class regions, learning- and test
  observations in the plane.

  Appropriate variables are usually found by \code{\link{GeneSelection}}.}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{Planarplot}}.}
\keyword{multivariate}