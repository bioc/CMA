\name{gbmCMA}
\alias{gbmCMA}
\title{Tree-based Gradient Boosting}
\description{Roughly speaking, Boosting combines 'weak learners'
            in a weighted manner in a stronger ensemble. This
            method calls the function \code{gbm.fit} from the
            package \code{gbm}. The 'weak learners' are
            simple trees that need only very few splits (default: 1).

            For \code{S4} method information, see \code{\link{gbmCMA-methods}}.}
\usage{
gbmCMA(X, y, f, learnind, models=FALSE,...)
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
                  \item{models}{a logical value indicating whether the model object shall be returned }
  \item{\dots}{Further arguments passed to the function \code{gbm.fit}
               from the package of the same name. Worth mentionning are
               \describe{
               \item{\code{ntrees}}{Number of trees to fit (size of the ensemble),
                                    defaults to 100. This parameter should
                                    be optimized using \code{\link{tune}}.}
               \item{\code{shrinkage}}{The learning rate (default is
		 0.001). Usually fixed to a very low value.}
               \item{\code{distribution}}{Loss function to be used. Default is \code{"bernoulli"},
                           i.e. \code{LogitBoost}, a (less robust) alternative
                           is \code{"adaboost"}.}
               \item{\code{interaction.depth}}{Number of splits used by the 'weak learner'
                 (single decision tree). Default is \code{1}.}
                 
	     }
                }
}

\note{Up to now, this method can only be applied to binary classification.}
\value{An onject of class \code{\link{cloutput}}.}
\references{ Ridgeway, G. (1999).

            The state of boosting.

            \emph{Computing Science and Statistics, 31:172-181}

            Friedman, J. (2001).

            Greedy Function Approximation: A Gradient Boosting Machine.

            \emph{Annals of Statistics 29(5):1189-1232}.
            }


\author{Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}}

\seealso{\code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}},
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
### run tree-based gradient boosting (no tuning)
gbmresult <- gbmCMA(X=golubX, y=golubY, learnind=learnind, n.trees = 500)
show(gbmresult)
ftable(gbmresult)
plot(gbmresult)
}

\keyword{multivariate}
