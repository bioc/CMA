\name{plrCMA-methods}
\docType{methods}
\alias{plrCMA-methods}
\alias{plrCMA,matrix,numeric,missing-method}
\alias{plrCMA,matrix,factor,missing-method}
\alias{plrCMA,data.frame,missing,formula-method}
\alias{plrCMA,ExpressionSet,character,missing-method}
\title{L2 penalized logistic regression}
\description{
 High dimensional logistic regression combined with an L2-type (Ridge-)penalty.
 Multiclass case is also possible.
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{plrCMA}}.
}
\keyword{multivariate}