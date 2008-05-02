\name{compBoostCMA-methods}
\docType{methods}
\alias{compBoostCMA-methods}
\alias{compBoostCMA,matrix,numeric,missing-method}
\alias{compBoostCMA,matrix,factor,missing-method}
\alias{compBoostCMA,data.frame,missing,formula-method}
\alias{compBoostCMA,ExpressionSet,character,missing-method}
\title{Componentwise Boosting}
\description{Roughly speaking, Boosting combines 'weak learners'
            in a weighted manner in a stronger ensemble.\cr
            'Weak learners' here consist of linear functions
            in one component (variable), as proposed by
            Buehlmann and Yu (2003).\cr
            It also generates sparsity and can as well be as
            used for variable selection alone. (s. \code{\link{GeneSelection}}.)
}
\section{Methods}{
\describe{

\item{X = "matrix", y = "numeric", f = "missing"}{signature 1}

\item{X = "matrix", y = "factor", f = "missing"}{signature 2}

\item{X = "data.frame", y = "missing", f = "formula"}{signature 3}

\item{X = "ExpressionSet", y = "character", f = "missing"}{signature 4}
}
For further argument and output information, consult
\code{\link{compBoostCMA}.}
}
\keyword{multivariate}
