\name{nnetCMA-methods}
\docType{methods}
\alias{nnetCMA-methods}
\alias{nnetCMA,matrix,numeric,missing-method}
\alias{nnetCMA,matrix,factor,missing-method}
\alias{nnetCMA,data.frame,missing,formula-method}
\alias{nnetCMA,ExpressionSet,character,missing-method}
\title{Feed-Forward Neural Networks}
\description{
This method provides access to the function
\code{nnet} in the package of the same name that trains
Feed-forward Neural Networks with one hidden layer.
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{nnetCMA}}.
}
\keyword{multivariate}
