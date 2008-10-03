\name{nnetCMA}
\alias{nnetCMA}
\title{Feed-forward Neural Networks}
\description{This method provides access to the function
\code{nnet} in the package of the same name that trains
Feed-forward Neural Networks with one hidden layer.

For \code{S4} method information, see \link{nnetCMA-methods}
}
\usage{
nnetCMA(X, y, f, learnind, eigengenes = FALSE, ...)
}
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
  \item{eigengenes}{Should the training be performed be in the space of
                    eigengenes obtained from a singular value decomposition
                    of the Gene expression data matrix ? Default is \code{FALSE};
                    in this case, variable selection is necessary to reduce
                    the number of weights that have to be optimized.}
  \item{\dots}{Further arguments passed to the function \code{nnet}
               from the package of the same name.

               Important parameters are:
               \itemize{
               \item \code{"size"}, i.e. the number of units in the hidden layer
               \item \code{"decay"} for weight decay.
	     }
	   }
	   }
\note{\itemize{
\item Excessive variable selection is usually necessary if \code{eigengenes = FALSE}
\item Different runs of this method on the same dataset not necessarily produce
      the same results due to the fact that optimization for Feed-Forward
      Neural Networks is rather difficult and depends on the choice of
      (normally randomly chosen) starting values for the network
      weights.
      }
}
\value{An object of class \code{\link{cloutput}}.}
\references{Ripley, B.D. (1996)

            Pattern Recognition and Neural Networks.

            \emph{Cambridge University Press}}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de}

        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{\code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
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
### extract gene expression from first 10 genes
golubX <- as.matrix(golub[,2:11])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run nnet (not tuned)
nnetresult <- nnetCMA(X=golubX, y=golubY, learnind=learnind, size = 3, decay = 0.01)
### show results
show(nnetresult)
ftable(nnetresult)
plot(nnetresult)
### in the space of eigengenes (not tuned)
golubXfull <-  as.matrix(golubX[,-1])
nnetresult <- nnetCMA(X=golubXfull, y=golubY, learnind = learnind, eigengenes = TRUE,
                      size = 3, decay = 0.01)
### show results
show(nnetresult)
ftable(nnetresult)
plot(nnetresult)
}
\keyword{multivariate}
