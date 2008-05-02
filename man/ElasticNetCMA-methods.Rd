\name{ElasticNetCMA-methods}
\docType{methods}
\alias{ElasticNetCMA-methods}
\alias{ElasticNetCMA,matrix,numeric,missing-method}
\alias{ElasticNetCMA,matrix,factor,missing-method}
\alias{ElasticNetCMA,data.frame,missing,formula-method}
\alias{ElasticNetCMA,ExpressionSet,character,missing-method}
\title{Classfication and variable selection by the ElasticNet}
\description{
 Zou and Hastie (2004) proposed a combined L1/L2 penalty
 for regularization and variable selection.
 The Elastic Net penalty encourages a grouping
 effect, where strongly correlated predictors tend to be in or out of the model together.
 The computation is done with the function \code{glmpath} from the package
 of the same name.
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For references, further argument and output information, consult
\code{\link{ElasticNetCMA}}
}
\keyword{multivariate}