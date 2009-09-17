\name{roc}
\alias{roc,cloutput-method}
\alias{roc}
\title{Receiver Operator Characteristic}
\description{The empirical Receiver Operator Characteristic (ROC)
             is widely used for the evaluation of diagnostic tests,
             but also for the evaluation of classfiers.
             In this implementation, it can only be used for the binary
             classification case.
             The input are a numeric vector of class probabilities
             (which play the role of a test result) and the true class labels.
             Note that misclassifcation performance can (partly widely) differ
             from the Area under the ROC (AUC). This is due to the fact that
             misclassifcation rates are always computed for the threshold
             'probability = 0.5'.}

\arguments{
  \item{object}{An object of \code{\link{cloutput}}.}
  \item{plot}{Should the ROC curve be plotted ? Default is \code{TRUE}.}
  \item{\dots}{Argument to specifiy further graphical options.}
}
\value{The empirical area under the curve (AUC).}
\references{

Slawski, M. Daumer, M.  Boulesteix, A.-L. (2008)
CMA - A comprehensive Bioconductor package for supervised classification with high dimensional data.
\emph{BMC Bioinformatics 9: 439}

}
\author{Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}}
\seealso{\code{\link{evaluation}}}
\keyword{multivariate}
