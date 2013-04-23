\name{weighted.mcr-methods}
\docType{methods}
\alias{weighted.mcr-methods}
\alias{weighted.mcr,character,character,numeric,character,matrix,numeric-method}
\alias{weighted.mcr,character,character,numeric,character,matrix,factor-method}
\alias{weighted.mcr,character,character,missing,character,matrix,factor-method}

\title{General method for tuning / selection bias correction}
\description{
  Perform tuning / selection bias correction in subsampling for the
  following signatures:
}
\section{Methods}{
\describe{

\item{classifiers="character",parameters="character",nbgenes="numeric",sel.method="character",X = "matrix", y = "numeric"}{signature 1}

\item{classifiers="character",parameters="character",nbgenes="numeric",sel.method="character",X = "matrix", y = "factor"}{signature 2}

\item{classifiers="character",parameters="character",nbgenes="missing",sel.method="character",X = "matrix", y = "factor"}{signature 3}

}
For further argument and output information, consult
\code{\link{weighted.mcr}}.}
\keyword{tuning bias, selection bias, corrected misclassification rate}
