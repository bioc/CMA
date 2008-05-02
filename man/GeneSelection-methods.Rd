\name{GeneSelection-methods}
\docType{methods}
\alias{GeneSelection-methods}
\alias{GeneSelection,matrix,numeric,missing-method}
\alias{GeneSelection,matrix,factor,missing-method}
\alias{GeneSelection,data.frame,missing,formula-method}
\alias{GeneSelection,ExpressionSet,character,missing-method}
\title{General method for variable selection with various methods}
\description{
  Performs gene selection for the following signatures:
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{GeneSelection}}.
}
\keyword{multivariate}