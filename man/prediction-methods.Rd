\name{prediction-methods}
\docType{methods}
\alias{prediction-methods}
\alias{prediction,matrix,ANY,matrix,missing-method}
\alias{prediction,data.frame,missing,data.frame,formula-method}
\alias{prediction,ExpressionSet,character,ExpressionSet,missing-method}
\title{General method for predicting class lables of new observations}
\description{
  Perform prediction signatures:
}
\section{Methods}{
\describe{

\item{X.tr = "matrix", X.new="matrix", y.tr='any',f = "missing"}{signature 1}

\item{X.tr = "data.frame", X.new="data.frame", y.tr = "missing", f = "formula"}{signature 2}

\item{X.tr = "ExpressionSet",X.new = "ExpressionSet", y.tr = "character", f = "missing"}{signature 3}



}
For further argument and output information, consult
\code{\link{classification}}.}
\keyword{multivariate}
