\name{clvarseloutput-class}
\docType{class}
\alias{clvarseloutput-class}
\alias{clvarseloutput}
\title{"clvarseloutput"}
\description{Object returned by all classifiers that can peform variable selection
             or compute variable importance. These are:
             \itemize{
             \item Random Forest, s. \code{\link{rfCMA}},
             \item Componentwise Boosting, s. \code{\link{compBoostCMA}},
             \item LASSO-logistic regression, s. \code{\link{LassoCMA}},
             \item ElasticNet-logistic regression,
	     s. \code{\link{ElasticNetCMA}}
	     }.
             Objects of class \code{clvarseloutput} extend both the class
             \code{cloutuput} and \code{varsel}, s. below.}
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
    \item{\code{varsel}:}{\code{numeric} vector of variable importance measures (for Random Forest) or
                          absolute values of regression coefficients (for the other three methods mentionned above)
                          (from which the majority will be zero).}
  }
}
\section{Extends}{
Class \code{"\linkS4class{cloutput}"}, directly.
Class \code{"\linkS4class{varseloutput}"}, directly.
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
\seealso{\code{\link{rfCMA}}, \code{\link{compBoostCMA}}, \code{\link{LassoCMA}}, \code{\link{ElasticNetCMA}}}
\keyword{multivariate}
