\name{LassoCMA-methods}
\docType{methods}
\alias{LassoCMA-methods}
\alias{LassoCMA,matrix,numeric,missing-method}
\alias{LassoCMA,matrix,factor,missing-method}
\alias{LassoCMA,data.frame,missing,formula-method}
\alias{LassoCMA,ExpressionSet,character,missing-method}
\title{L1 penalized logistic regression}
\description{
 The Lasso (Tibshirani, 1996) is one of the most popular
 tools for simultaneous shrinkage and variable selection. Recently,
 Friedman, Hastie and Tibshirani (2008) have developped and algorithm to
 compute the entire solution path of the Lasso for an arbitrary
 generalized linear model, implemented in the package \code{glmnet}.
 The method can be used for variable selection alone, s. \code{\link{GeneSelection}}
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For references, further argument and output information, consult
\code{\link{LassoCMA}}.
}
\keyword{multivariate}
