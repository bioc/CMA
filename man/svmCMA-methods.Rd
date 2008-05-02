\name{svmCMA-methods}
\docType{methods}
\alias{svmCMA-methods}
\alias{svmCMA,matrix,numeric,missing-method}
\alias{svmCMA,matrix,factor,missing-method}
\alias{svmCMA,data.frame,missing,formula-method}
\alias{svmCMA,ExpressionSet,character,missing-method}
\title{Support Vector Machine}
\description{Calls the function \code{svm} from the package \code{e1071}
that provides an interface to the award-winning LIBSVM routines.}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{svmCMA}}.}
\keyword{multivariate}