\name{flexdaCMA}
\alias{flexdaCMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Flexible Discriminant Analysis}
\description{
 This method is experimental.
 
 It is easy to show that, after appropriate scaling of the predictor matrix \code{X},
 Fisher's Linear Discriminant Analysis is equivalent to Discriminant Analysis
 in the space of the  fitted values from the linear regression of the
 \code{nlearn x K} indicator matrix of the class labels on \code{X}.
 This gives rise to 'nonlinear discrimant analysis' methods that expand
 \code{X} in a suitable, more flexible basis. In order to avoid overfitting,
 penalization is used. In the implemented version, the linear model is replaced
 by a generalized additive one, using the package \code{mgcv}.

 For \code{S4} method information, s. \code{\link{flexdaCMA-methods}}.
}
\usage{
flexdaCMA(X, y, f, learnind, comp = 1, plot = FALSE, ...)
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
  \item{comp}{Number of discriminant coordinates (projections) to compute.
              Default is one, must be smaller than or equal to \code{K-1}, where
              \code{K} is the number of classes.}
  \item{plot}{Should the projections onto the space spanned by the optimal
              projection directions be plotted ? Default is \code{FALSE}.}
  \item{\dots}{Further arguments passed to the function \code{gam} from
               the package \code{mgcv}.}
}

\note{Excessive variable selection has usually to performed before
      \code{flexdaCMA} can be applied in the \code{p > n} setting.
      Recall that the original predictor dimension is even enlarged,
      therefore, it should be applied only with very few variables.}

\value{An object of class \code{\link{cloutput}}.}
\references{Ripley, B.D. (1996)

            Pattern Recognition and Neural Networks.

            \emph{Cambridge University Press}}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de}

        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{\code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{gbmCMA}},
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
### extract gene expression from first 5 genes
golubX <- as.matrix(golub[,2:6])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run flexible Discriminant Analysis
result <- flexdaCMA(X=golubX, y=golubY, learnind=learnind, comp = 1)
### show results
show(result)
ftable(result)
plot(result)
}
\keyword{multivariate}