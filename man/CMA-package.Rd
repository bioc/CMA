\name{CMA-package}
\alias{CMA-package}
\alias{CMA}
\docType{package}
\title{Synthesis of microarray-based classification}
\description{The aim of the package is to provide a user-friendly
environment for the evaluation of classification methods using
gene expression data. A strong focus is on combined variable selection,
hyperparameter tuning, evaluation, visualization and comparison of (up to now) 21
classification methods from three main fields: Discriminant Analysis,
Neural Networks and Machine Learning. Although the package has been
created with the intention to be used for Microarray data, it can as well
be used in various \code{(p > n)}-scenarios.}
\details{
\tabular{ll}{
Package: & CMA \\
Type: & Package \\
Version: & 0.8.5 \\
Date: & 2008-2-13 \\
License: & GPL (version 2 or later) \\
}
Most Important Steps for the workflow are:
\describe{
\item{1.}{Generate evaluation datasets using \code{\link{GenerateLearningsets}}}
\item{2.}{(Optionally): Perform variable selection using \code{\link{GeneSelection}}}
\item{3.}{(Optionally): Peform hyperparameter tuning using \code{\link{tune}}}
\item{4.}{Perform \code{\link{classification}} using 1.-3.}
\item{5.}{Repeat 2.-4. based on 1. for several methods:
          \code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}},
         \code{\link{pnnCMA}}, \code{\link{qdaCMA}}, \code{\link{rfCMA}},
         \code{\link{scdaCMA}}, \code{\link{shrinkldaCMA}}, \code{\link{svmCMA}}}
\item{6.}{Evaluate the results from 5. using \code{\link{evaluation}} and make a comparison
          by calling \code{\link{compare}}}
      }
      }
\author{
Martin Slawski \email{martin.slawski@campus.lmu.de}

Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}

Maintainer: Martin Slawski \email{martin.slawski@campus.lmu.de}.
}
\keyword{multivariate}
