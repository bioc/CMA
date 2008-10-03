\name{pknnCMA-methods}
\docType{methods}
\alias{pknnCMA-methods}
\alias{pknnCMA,matrix,numeric,missing-method}
\alias{pknnCMA,matrix,factor,missing-method}
\alias{pknnCMA,data.frame,missing,formula-method}
\alias{pknnCMA,ExpressionSet,character,missing-method}
\title{Probabilistic nearest neighbours}
\description{
 Nearest neighbour variant that replaces the simple
  voting scheme by a weighted one (based on euclidean distances).
  This is also used to compute class probabilities.
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{pknnCMA}.}
}
\keyword{multivariate}
