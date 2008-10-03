\name{compare}
\alias{compare}
\title{Compare different classifiers}
\description{
  Classifiers can be evaluated separately using the method
  \code{\link{evaluation}}. Normally, several classifiers
  are used for the same dataset and their performance is
  compared. This comparison procedure is essentially facilitated by
  this method.
  For S4 method information, s. \code{\link{compare-methods}}}

\usage{
compare(clresultlist, measure = c("misclassification", "sensitivity", "specifity", "average probability", "brier score", "auc"), aggfun = mean, plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clresultlist}{A list of lists (!) of objects of class \code{\link{cloutput}} or \code{\link{clvarseloutput}}.
                      Each inner list is usually returned by \code{\link{classification}}.
                      Additionally, the different list elements of the outer list
                      should have been created by different classifiers, s.
                      also example below.}
  \item{measure}{A character vector containing one or more of the elements listed below.
                 By default, all measures are computed, using \code{\link{evaluation}}
                 with \code{scheme = "iterationwise"}.
                 Note that \code{"sensitivity", "specifity", "auc"} cannot be computed
                 for the multiclass case.
                 \describe{
                 \item{\code{"misclassification"}}{The missclassifcation rate.}
                 \item{\code{"sensitivity"}}{The sensitivity or 1-false negative rate.
                                             Can only be computed for binary classifcation.}
                 \item{\code{"specifity"}}{The specifity or 1-false positive rate.
                                           Can only be computed for binary classification.}
                 \item{\code{"average probability"}}{The average probability assigned to the
                                                     correct class. Requirement is that
                                                     the used classifier provides probability
                                                     estimations. The optimum performance is 1.}
                 \item{\code{"brier score"}}{The Brier Score is generally defined as
                                             \code{<sum over all observation i> <sum over all classes k> (I(y_i=k)-P(k))^2},
                                             with \code{I()} denoting the indicator function and \code{P(k)} the estimated
                                             probability for class \code{k}. The optimum performance is 0.}
                 \item{\code{"auc"}}{The Area under the Curve (AUC) belonging to the empirical ROC curve
                                     computed from the estimated probabilities and the true class labels.
                                     Can only be computed for binary classification and if \code{"scheme = iterationwise"}, s. below.
                                     S. also \code{roc,cloutput-method}.}
                 }}
  \item{aggfun}{Function that determines how performance among different iterations are aggregared.
                Default is \code{mean}, other possible choices are quantiles.}
  \item{plot}{Should the performance of different classifiers be visualized by a joint boxplot ?
               Default is \code{FALSE}.}
  \item{\dots}{Further arguments passed to \code{boxplot} in the case that \code{plot = TRUE}.}
}
\note{If more than one measure is computed and \code{plot = TRUE}, one separate
     plot is created for each of them.}
\value{A \code{data.frame} with rows corresponding to the compared classifiers
       and columns to the performance measures, aggregated by \code{aggfun}, s. above.}
\references{ Dudoit, S., Fridlyand, J., Speed, T. P. (2002)

            Comparison of discrimination methods for the classification of tumors
            using gene expression data.

            \emph{Journal of the American Statistical Association 97, 77-87}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de}

        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
        
\seealso{\code{\link{classification}}, \code{\link{evaluation}}}
\examples{
\dontrun{
### compare the performance of several discriminant analysis methods
### for the Khan dataset:
data(khan)
khanX <- as.matrix(khan[,-1])
khanY <- khan[,1]
set.seed(27611)
fiveCV10iter <- GenerateLearningsets(y=khanY, method = "CV", fold = 5, niter = 2, strat = TRUE)
### candidate methods:  DLDA, LDA, QDA, pls_LDA, sclda
class_dlda <- classification(X = khanX, y=khanY, learningsets = fiveCV10iter, classifier = dldaCMA)
### peform GeneSlection for LDA, FDA, QDA (using F-Tests):
genesel_da <- GeneSelection(X=khanX, y=khanY, learningsets = fiveCV10iter, method = "f.test")
###
class_lda <- classification(X = khanX, y=khanY, learningsets = fiveCV10iter, classifier = ldaCMA, genesel= genesel_da, nbgene = 10)

class_qda <- classification(X = khanX, y=khanY, learningsets = fiveCV10iter, classifier = qdaCMA, genesel = genesel_da, nbgene = 2)

### We now make a comparison concerning the performance (sev. measures):
### first, collect in a list:
dalike <- list(class_dlda, class_lda, class_qda)
### use pre-defined compare function:
comparison <- compare(dalike, plot = TRUE, measure = c("misclassification", "brier score", "average probability"))
print(comparison)
}
}
\keyword{multivariate}
