\name{evaluation}
\alias{evaluation}
\title{Evaluation of classifiers}
\description{
  The performance of classifiers can be evaluted by six different
  measures and two different schemes that are described more precisely
  below.\cr
  For S4 method information, s. \code{\link{evaluation-methods}}.
}
\usage{
evaluation(clresult, cltrain = NULL, cost = NULL, y = NULL, measure = c("misclassification", "sensitivity", "specifity", "average probability", "brier score", "auc", "0.632", "0.632+"),
                     scheme = c("iterationwise", "observationwise", "classwise"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clresult}{A list of objects of class \code{\link{cloutput}} or \code{\link{clvarseloutput}}
                  }
  \item{cltrain}{An object of class \code{\link{cloutput}} in which the
                  \emph{whole} dataset was used as learning set. Only
                   used if \code{method = "0.632"} or \code{method = "0.632+"}
                   in order to obtain an estimation for the resubsitution
                   error rate.}
  \item{cost}{An optional cost matrix used if \code{measure = "misclassification"}.
              If it is not specified (default), the cost is the usual indicator
              loss. Otherwise, entry \code{i,j} of \code{cost} quantifies the
              loss when the true class is class \code{i-1} and the predicted
              class is \code{j-1}, provided the conventional coding \code{0,...,K-1}
              in the case of \code{K} classes is used. Usually, the matrix contains
              only non-negative entries with zeros on the diagonal, but this is
              not obligatory. Make sure that the dimension of the matrix matches
              the number of classes.}
  \item{y}{A vector containing the true class labels. Only needed if \code{scheme = "classwise"}.}
  \item{measure}{Peformance measure to be used:
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
                 \item{\code{"0.632"}}{The 0.632 estimator (s. reference) for the misclassification rate (applied iteration- or)
                                     observationwise, if bootstrap learning sets have been used. Note that
                                     \code{cltrain} must be provided.}
                 \item{\code{"0.632+"}}{The 0.632+ estimator (s. reference) for the misclassification rate (applied iteration- or)
                                     observationwise, if bootstrap learning sets have been used. Note that
                                     \code{cltrain} must be provided.}
  
                 }
                 }
                \item{scheme}{\describe{
                              \item{\code{"iterationwise"}}{The performance measures listed above are computed for
                                                            each different iteration, i.e. each different \code{learningset}}

                              \item{\code{"observationwise"}}{The performance measures listed above (except for \code{"auc"})
                                                              are computed separately for each observation classified one
                                                              or several times, depending on the \code{learningset} scheme.}
                              \item{\code{"classwise"}}{The performance measures (exceptions: \code{"auc", "0.632", "0.632+"})
                                                        are computed separately for each class, averaged over both iterations
                                                        and observations.}
                              }
                            }
      }
\value{An object of class \code{\link{evaloutput}}.}
\references{Efron, B. and Tibshirani, R. (1997). 
            Improvements on cross-validation: The .632+ bootstrap method.\cr 
            \emph{Journal of the American Statistical Association, 92, 548-560}.}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\code{\link{evaloutput}}, \code{\link{classification}}, \code{\link{compare}}}
\examples{
### simple linear discriminant analysis example using bootstrap datasets:
### datasets:
data(golub)
golubY <- golub[,1]
### extract gene expression from first 10 genes
golubX <- as.matrix(golub[,2:11])
### generate 25 bootstrap datasets
set.seed(333)
bootds <- GenerateLearningsets(y = golubY, method = "bootstrap", ntrain = 30, niter = 10, strat = TRUE)
### run classification()
ldalist <- classification(X=golubX, y=golubY, learningsets = bootds, classifier=ldaCMA)
### Evaluation:
eval_iter <- evaluation(ldalist, scheme = "iter")
eval_obs <- evaluation(ldalist, scheme = "obs")
show(eval_iter)
show(eval_obs)
summary(eval_iter)
summary(eval_obs)
### auc with boxplot
eval_auc <- evaluation(ldalist, scheme = "iter", measure = "auc")
boxplot(eval_auc)
### which observations have often been misclassified ?
obsinfo(eval_obs, threshold = 0.75)
}

\keyword{multivariate}

