\name{gbmCMA-methods}
\docType{methods}
\alias{gbmCMA-methods}
\alias{gbmCMA,matrix,numeric,missing-method}
\alias{gbmCMA,matrix,factor,missing-method}
\alias{gbmCMA,data.frame,missing,formula-method}
\alias{gbmCMA,ExpressionSet,character,missing-method}
\title{Tree-based Gradient Boosting}
\description{
 Roughly speaking, Boosting combines 'weak learners'
 in a weighted manner in a stronger ensemble. This
 method calls the function \code{gbm.fit} from the
 package \code{gbm}. The 'weak learners' are
 simple trees that need only very few splits (default: 1).
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{gbmCMA}.}
}
\keyword{multivariate}