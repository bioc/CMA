\name{GenerateLearningsets}
\alias{GenerateLearningsets}
\title{Repeated Divisions into learn- and tets sets}
\description{
  Due to very small sample sizes, the classical division learnset/testset
  does not give accurate information about the classification performance.
  Therefore, several different divisions should be used and aggregated.
  The implemented methods are discussed in Braga-Neto and Dougherty (2003)
  and Molinaro et al. (2005) whose terminology is adopted.

  This function is usually the basis for all deeper analyses.
}
\usage{
GenerateLearningsets(n, y, method = c("LOOCV", "CV", "MCCV", "bootstrap"),
                     fold = NULL, niter = NULL, ntrain = NULL, strat = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{The total number of observations in the available data set. May be \code{missing}
           if \code{y} is provided instead.}
  \item{y}{A vector of class labels, either \code{numeric} or a \code{factor}.
            \emph{Must} be given if \code{strat=TRUE} or \code{n} is not specified.}
  \item{method}{Which kind of scheme should be used to generate divisions
                into learning sets and test sets ? Can be one of the following:
                 \describe{
                 \item{"LOOCV"}{Leaving-One-Out Cross Validation.}
                 \item{"CV"}{(Ordinary) Cross-Validation. Note that \code{fold}
                             must as well be specified.}
                 \item{"MCCV"}{Monte-Carlo Cross Validation, i.e.
                               random divisions into learning sets with
                               \code{ntrain}(s.below) observations
                               and tests sets with \code{ntrain} observations.}
                 \item{"bootstrap"}{Learning sets are generated
                                    by drawing \code{ntrain} times with replacement from
                                    all observations. Those not drawn not all
                                    form the test set.}
                }
                }
  \item{fold}{Gives the number of CV-groups. Used only when \code{method="CV"}}
  \item{niter}{Number of iterations (s.\code{details).}}
  \item{ntrain}{Number of observations in the learning sets. Used
                only when \code{method="MCCV"} or \code{method="bootstrap"}.}
  \item{strat}{Logical. Should stratified sampling be performed,
                  i.e. the proportion of observations from each class in the learning
                  sets be the same as in the whole data set ?

                  Does not apply for \code{method = "LOOCV"}.}		
}

\details{
\itemize{
  \item When \code{method="CV"}, \code{niter} gives the number of times
  the whole CV-procedure is repeated. The output matrix has then \code{fold}x\code{niter} rows.
  When \code{method="MCCV"} or \code{method="bootstrap"}, \code{niter} is simply the number of considered
  learning sets.
\item Note that \code{method="CV",fold=n} is equivalent to \code{method="LOOCV"}.
}
}

\value{An object of class \code{\link{learningsets}}}
\references{Braga-Neto, U.M., Dougherty, E.R. (2003).

            Is cross-validation valid for small-sample microarray classification ?

            \emph{Bioinformatics, 20(3), 374-380}

            Molinaro, A.M., Simon, R., Pfeiffer, R.M. (2005).

            Prediction error estimation: a comparison of resampling
            methods.

            \emph{Bioinformatics, 21(15), 3301-3307}
            
            Slawski, M. Daumer, M.  Boulesteix, A.-L. (2008)
CMA - A comprehensive Bioconductor package for supervised classification with high dimensional data.
\emph{BMC Bioinformatics 9: 439}
            
            }


\author{Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}
        
        Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}}

\seealso{\code{\link{learningsets}}, \code{\link{GeneSelection}}, \code{\link{tune}},
         \code{\link{classification}}}
\examples{
# LOOCV
loo <- GenerateLearningsets(n=40, method="LOOCV")
show(loo)
# five-fold-CV
CV5 <- GenerateLearningsets(n=40, method="CV", fold=5)
show(loo)
# MCCV
mccv <- GenerateLearningsets(n=40, method = "MCCV", niter=3, ntrain=30)
show(mccv)
# Bootstrap
boot <- GenerateLearningsets(n=40, method="bootstrap", niter=3)
# stratified five-fold-CV
set.seed(113)
classlabels <- sample(1:3, size = 50, replace = TRUE, prob = c(0.3, 0.5, 0.2))
CV5strat <- GenerateLearningsets(y = classlabels, method="CV", fold=5, strat = TRUE)
show(CV5strat)
}

\keyword{multivariate}
