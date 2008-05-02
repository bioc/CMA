\name{filter}
\alias{ttest}
\alias{welchtest}
\alias{wilcoxtest}
\alias{ftest}
\alias{kruskaltest}
\alias{limmatest}
\alias{golubcrit}
\alias{rfe}
\title{Filter functions for Gene Selection}
\description{The functions listed above are usually not called by the
             user but via \code{\link{GeneSelection}}.}
\usage{
ttest(X, y, learnind, ...)
welchtest(X, y, learnind, ...)
ftest(X, y, learnind,...)
kruskaltest(X, y, learnind,...)
limmatest(X, y, learnind,...)
golubcrit(X, y, learnind,...)
rfe(X, y, learnind,...)
}
\arguments{
 \item{X}{A \code{numeric} matrix of gene expression values.}
 \item{y}{A \code{numeric} vector of class labels.}
 \item{learnind}{An index vector specifying the observations that
                  belong to the learning set.}
 \item{...}{Currently unused argument.}}
                  
\value{An object of class \code{\link{varseloutput}}.}
\keyword{multivariate}

