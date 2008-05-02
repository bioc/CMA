\name{fdaCMA-methods}
\docType{methods}
\alias{fdaCMA-methods}
\alias{fdaCMA,matrix,numeric,missing-method}
\alias{fdaCMA,matrix,factor,missing-method}
\alias{fdaCMA,data.frame,missing,formula-method}
\alias{fdaCMA,ExpressionSet,character,missing-method}
\title{Fisher's Linear Discriminant Analysis}
\description{
 Fisher's Linear Discriminant Analysis constructs a subspace of
 'optimal projections' in which classification is performed.
 The directions of optimal projections are computed by the
 function \code{cancor} from the package \code{stats}. For
 an exhaustive treatment, see e.g. Ripley (1996).\cr
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For references, further argument and output information, consult
\code{\link{fdaCMA}}.
}
\keyword{multivariate}