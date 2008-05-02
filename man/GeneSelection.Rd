\name{GeneSelection}
\alias{GeneSelection}
\title{General method for variable selection with various methods}
\description{
 For different learning data sets as defined by the argument \code{learningsets},
 this method ranks the genes from the most relevant to the less relevant using
 one of various 'filter' criteria or provides a sparse collection of variables
 (Lasso, ElasticNet, Boosting). The results are typically used for variable selection for
 the classification procedure that follows.\cr
 For S4 class information, s. \code{GeneSelection-methods}.
}
\usage{
GeneSelection(X, y, f, learningsets, method = c("t.test", "welch.test", "wilcox.test", "f.test", "kruskal.test", "limma", "rfe", "rf", "lasso", "elasticnet", "boosting", "golub"), scheme, trace = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Gene expression data. Can be one of the following:
           \itemize{
           \item A \code{matrix}. Rows correspond to observations, columns to variables.
           \item A \code{data.frame}, when \code{f} is \emph{not} missing (s. below).
           \item An object of class \code{ExpressionSet}.
       }
     }
  \item{y}{Class labels. Can be one of the following:
           \itemize{
           \item A \code{numeric} vector.
           \item A \code{factor}.
           \item A \code{character} if \code{X} is an \code{ExpressionSet}.
           \item \code{missing}, if \code{X} is a \code{data.frame} and a
           proper formula \code{f} is provided.
       }
     }
       
  \item{f}{A two-sided formula, if \code{X} is a \code{data.frame}. The
           left part correspond to class labels, the right to variables.}
  \item{learningsets}{An object of class \code{\link{learningsets}}. May
                      be missing, then the complete datasets is used as
                      learning set.}
  \item{method}{A character specifying the method to be used:
                \describe{
                \item{\code{t.test}}{two-sample t.test (equal variances for both classes assumed).}
                \item{\code{welch.test}}{Welch modification of the t.test (unequal variances for both classes).}
                \item{\code{wilcox.test}}{Wilcoxon rank sum test.}
                \item{\code{f.test}}{F test belonging to the linear hypothesis that the mean is the same for all classes. Usually used for the
                                    multiclass scheme, is equivalent to \code{method = t.test} in the two-class case.}
                \item{\code{kruskal.test}}{Multi-class generalization of the Wilcoxon rank sum test and the nonparametric
                                           pendant to the F test, respectively.}
                \item{\code{limma}}{'Moderated t' statistic for the two-class case and 'moderated F' statistic for
                                     the multiclass case, described in Smyth (2003). Requires the package
                                     \code{limma}.}
                \item{\code{rfe}}{One-step Recursive Feature Elimination, based on the Support Vector Machine.
                                    The method is decribed in Guyon et al. (2002). Requires the package \code{e1071}.
                                    Take care that appropriate hyperparameters are passed by the \code{...} argument.}
                \item{\code{rf}}{Random Forest Variable Importance Measure. Requires the package \code{randomForest}}
                \item{\code{lasso}}{\code{L1} penalized logistic regression leads to sparsity with respect
                                      to the variables used. Calls the function \code{\link{LassoCMA}}, which requires
                                      the package \code{glmpath}.
                                      \bold{warning}: Take care that appropriate hyperparameters are passed by the \code{...} argument.}
                \item{\code{elasticnet}}{Penalized logistic regression with both \code{L1} and \code{L2} penalty, claimed
                                         by Zhou and Hastie (2004) to select 'variable groups'. Calls
                                         the function \code{\link{ElasticNetCMA}}, which requires the package \code{glmpath}}.
                                         \bold{warning}: Take care that appropriate hyperparameters are passed by the \code{...} argument.}
                \item{\code{boosting}}{Componentwise boosting (Buehlmann and Yu, 2003) has been shown to mimic
                                       the LASSO (Efron et al., 2004; Buehlmann and Yu, 2006). Calls the function \code{\link{compBoostCMA}}
                                       Take care that appropriate hyperparameters are passed by the \code{...} argument.}
                \item{\code{golub}}{The (theoretically unfounded) variable selection criterion used by Golub et al. (1999), s. \code{\link{golub}}.}

    }            
  \item{scheme}{The scheme to be used in the case of a non-binary response. Must be one
                of \code{"pairwise"},\code{"one-vs-all"} or \code{"multiclass"}. The
                last case only makes sense if \code{method} is one of \code{f.test, limma, rf, boosting},
                which can directly be applied to the multi class case.}
  \item{trace}{Should the progress be traced ? Default is \code{TRUE}.}
  \item{\dots}{Further arguments passed to the function performing variable selection, s. \code{method}.}
}
\note{most of the methods described above are only apt for the binary classification case. The only
     ones that can be used without restriction in the multiclass case are
     \itemize{
     \item \code{f.test}
     \item \code{kruskal.test}
     \item \code{rf}
     \item \code{boosting}
   }
     For the rest, pairwise or one-vs-all schemes are used.}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{An object of class \code{\link{genesel}}.}
\references{Smyth, G. K., Yang, Y.-H., Speed, T. P. (2003).\cr
            Statistical issues in microarray data analysis.\cr
            \emph{Methods in Molecular Biology 224, 111-136}.

            Guyon, I., Weston, J., Barnhill, S., Vapnik, V. (2002).\cr
            Gene Selection for Cancer Classification using support
            vector machines.
            \emph{Journal of Machine Learning Research, 46, 389-422}

            Zhou, H., Hastie, T. (2004).\cr
            Regularization and variable selection via the elastic net.
            \emph{Journal of the Royal Statistical Society B, 67(2),301-320}

            Buelmann, P., Yu, B. (2003).\cr
            Boosting with the L2 loss: Regression and Classification.\cr
            \emph{Journal of the American Statistical Association, 98, 324-339}
            
            Efron, B., Hastie, T., Johnstone, I., Tibshirani, R. (2004).\cr
            Least Angle Regression.\cr
            \emph{Annals of Statistics, 32:407-499}

            Buehlmann, P., Yu, B. (2006).\cr
            Sparse Boosting.\cr
            \emph{Journal of Machine Learning Research, 7- 1001:1024}
            }
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{\code{\link{filter}}, \code{\link{GenerateLearningsets}}, \code{\link{tune}},
         \code{\link{classification}}}
\examples{
# load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression from first 10 genes
golubX <- as.matrix(golub[,-1])
### Generate five different learningsets
set.seed(111)
five <- GenerateLearningsets(y=golubY, method = "CV", fold = 5, strat = TRUE)
### simple t-test:
selttest <- GeneSelection(golubX, golubY, learningsets = five, method = "t.test")
### show result:
show(selttest)
toplist(selttest, k = 10, iter = 1)
plot(selttest, iter = 1)
}
\keyword{multivariate}
