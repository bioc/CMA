\name{compBoostCMA}
\alias{compBoostCMA}
\title{Componentwise Boosting}
\description{Roughly speaking, Boosting combines 'weak learners'
            in a weighted manner in a stronger ensemble.

            'Weak learners' here consist of linear functions
            in one component (variable), as proposed by
            Buehlmann and Yu (2003).

            It also generates sparsity and can as well be as
            used for variable selection alone. (s. \code{\link{GeneSelection}}).

            For \code{S4} method information, see \code{\link{compBoostCMA-methods}.}
}
\usage{
compBoostCMA(X, y, f, learnind, loss = c("binomial", "exp", "quadratic"), mstop = 100, nu = 0.1, ...)
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
  \item{loss}{Character specifying the loss function - one of \code{"binomial"} (LogitBoost),
              \code{"exp"} (AdaBoost), "quadratic"(L2Boost).}
  \item{mstop}{Number of boosting iterations, i.e. number of updates
               to perform. The default (100) does not necessarily produce
               good results, therefore usage of \code{\link{tune}} for
               this argument is highly recommended.}
  \item{nu}{Shrinkage factor applied to the update steps, defaults to 0.1.
            In most cases, it suffices to set \code{nu} to a very low value
            and to concentrate on the optimization of \code{mstop}.}
  \item{\dots}{Currently unused arguments.}
}
\details{The method is partly based on code from the package \code{mboost}
         from T. Hothorn and P. Buehlmann.

         The algorithm for the multiclass case is described
         in Lutz and Buehlmann (2006) as 'rowwise updating'.}
\value{An object of class \code{\link{clvarseloutput}}.}
\references{  Buelmann, P., Yu, B. (2003).

              Boosting with the L2 loss: Regression and Classification.

              \emph{Journal of the American Statistical Association, 98, 324-339}

              Buehlmann, P., Hothorn, T.

              Boosting: A statistical perspective.

              \emph{Statistical Science (to appear)}
              
              Lutz, R., Buehlmann, P. (2006).

              Boosting for high-multivariate responses in high-dimensional linear regression.

              \emph{Statistica Sinica 16, 471-494}.
}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de}

        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{\code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
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
### run componentwise (logit)-boosting (not tuned)
result <- compBoostCMA(X=golubX, y=golubY, learnind=learnind, mstop = 500)
### show results
show(result)
ftable(result)
plot(result)
### multiclass example:
### load Khan data
data(khan)
### extract class labels
khanY <- khan[,1]
### extract gene expression
khanX <- as.matrix(khan[,-1])
### select learningset
set.seed(111)
learnind <- sample(length(khanY), size=floor(ratio*length(khanY)))
### run componentwise multivariate (logit)-boosting (not tuned)
result <- compBoostCMA(X=khanX, y=khanY, learnind=learnind, mstop = 1000)
### show results
show(result)
ftable(result)
plot(result)
}
\keyword{multivariate}