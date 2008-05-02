\name{cloutput-class}
\docType{class}
\alias{cloutput-class}
\alias{cloutput}
\alias{show,cloutput-method}

\title{"cloutput"}
\description{Object returned by one of the classifiers (functions ending with \code{CMA})}
\section{Slots}{
	 \describe{
    \item{\code{learnind}:}{Vector of indices that indicates which observations
                            where used in the learning set.}
    \item{\code{y}:}{Actual (true) class labels of predicted observations.}
    \item{\code{yhat}:}{Predicted class labels by the classifier.}
    \item{\code{prob}:}{A \code{numeric} \code{matrix} whose rows
                        equals the number of predicted observations (length of \code{y}/\code{yhat})
                        and whose columns equal the number of different classes in the learning set.
                        Rows add up to one.
                        Entry \code{j,k} of this matrix contains the probability for the \code{j}-th
                        predicted observation to belong to class \code{k}.
                        Can be a matrix of \code{NA}s, if the classifier used does not
                        provide any probabilities}
    \item{\code{method}:}{Name of the classifer used.}
    \item{\code{mode}:}{\code{character}, one of \code{"binary"} (if the number of classes in the learning set is two)
                          or \code{multiclass} (if it is more than two).}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(cloutput-object)} for brief information}
    \item{ftable}{Use \code{ftable(cloutput-object)} to obtain a confusion matrix/cross-tabulation
                  of \code{y} vs. \code{yhat}, s. \code{\link{ftable,cloutput-method}}.}
    \item{plot}{Use \code{plot(cloutput-object)} to generate a probability plot of the matrix
                \code{prob} described above, s. \code{\link{plot,cloutput-method}}}
    \item{roc}{Use \code{roc(cloutput-object)} to compute the empirical ROC curve and the
               Area Under the Curve (AUC) based on the predicted probabilities, s.\code{\link{roc,cloutput-method}}}
	 }
}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\code{clvarseloutput}
        \code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}},
         \code{\link{pnnCMA}}, \code{\link{qdaCMA}}, \code{\link{rfCMA}},
         \code{\link{scdaCMA}}, \code{\link{shrinkldaCMA}}, \code{\link{svmCMA}}}
\keyword{multivariate}
