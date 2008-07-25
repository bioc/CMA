\name{ElasticNetCMA}
\alias{ElasticNetCMA}
\title{Classfication and variable selection by the ElasticNet}
\description{
 Zou and Hastie (2004) proposed a combined L1/L2 penalty
 for regularization and variable selection.
 The Elastic Net penalty encourages a grouping
 effect, where strongly correlated predictors tend to be in or out of the model together.
 The computation is done with the function \code{glmpath} from the package
 of the same name.\cr
 The method can be used for variable selection alone, s. \code{\link{GeneSelection}}.\cr
 For \code{S4} method information, see \code{ElasticNetCMA-methods}.}
\usage{
ElasticNetCMA(X, y, f, learnind, norm.fraction = 0.1, lambda2=1e-3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Gene expression data. Can be one of the following:
           \itemize{
           \item A \code{matrix}. Rows correspond to observations, columns to variables.
           \item A \code{data.frame}, when \code{f} is \emph{not} missing (s. below).
           \item An object of class \code{ExpressionSet}.
            \bold{note}: by default, the predictors are scaled to have unit variance
           and zero mean. Can be changed by passing \code{standardize = FALSE} via
            the \code{...} argument. 
	   }
           }
  \item{y}{Class labels. Can be one of the following:
           \itemize{
           \item A \code{numeric} vector.
           \item A \code{factor}.
           \item A \code{character} if \code{X} is an \code{ExpressionSet} that
                 specifies the phenotype variable.
           \item \code{missing}, if \code{X} is a \code{data.frame} and a
           proper formula \code{f} is provided.
	   }
           \bold{WARNING}: The class labels will be re-coded to
           range from \code{0} to \code{K-1}, where \code{K} is the
           total number of different classes in the learning set.
           }
  \item{f}{A two-sided formula, if \code{X} is a \code{data.frame}. The
           left part correspond to class labels, the right to variables.}
  \item{learnind}{An index vector specifying the observations that
                  belong to the learning set. May be \code{missing};
                  in that case, the learning set consists of all
                  observations and predictions are made on the
                  learning set.}
  \item{norm.fraction}{L1 Shrinkage intensity, expressed as the fraction
                       of the coefficient L1 norm compared to the
                       maximum possible L1 norm (corresponds to \code{fraction = 1}).
                       Lower values correspond to higher shrinkage.
                       Note that the default (0.1) need not produce good
                       results, i.e. tuning of this parameter is recommended.}
  \item{lambda2}{L2 shrinkage intensity, a postive real number. The
                 default (0.001) need not produce good results.}
  \item{\dots}{Further arguments passed to the function \code{glmpath}
               from the package of the same name.}
}
\note{For a strongly related method, s. \code{\link{LassoCMA}}.\cr
      Up to now, this method can only be applied to binary classification.}
\value{An object of class \code{\link{clvarseloutput}}.}
\references{Zhou, H., Hastie, T. (2004).\cr
            Regularization and variable selection via the elastic net.
            \emph{Journal of the Royal Statistical Society B, 67(2),301-320}

             Young-Park, M., Hastie, T. (2007)\cr
             L1-regularization path algorithm for generalized linear models.\cr
             \emph{Journal of the Royal Statistical Society B, 69(4), 659-677}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\code{\link{compBoostCMA}}, \code{\link{dldaCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}},
         \code{\link{pnnCMA}}, \code{\link{qdaCMA}}, \code{\link{rfCMA}},
         \code{\link{scdaCMA}}, \code{\link{shrinkldaCMA}}, \code{\link{svmCMA}}}
\examples{
### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression
golubX <- as.matrix(golub[,-1])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run ElasticNet - penalized logistic regression (no tuning)
result <- ElasticNetCMA(X=golubX, y=golubY, learnind=learnind, norm.fraction = 0.2, lambda2=0.01)
show(result)
ftable(result)
plot(result)
}
\keyword{multivariate}