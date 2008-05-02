\name{pknnCMA}
\alias{pknnCMA}
\title{Probabilistic Nearest Neighbours}
\description{
  Nearest neighbour variant that replaces the simple
  voting scheme by a weighted one (based on euclidean distances).
  This is also used to compute class probabilities.\cr
  For \code{S4} class information, see \link{pknnCMA-methods}.
}
\usage{
pknnCMA(X, y, f, learnind, beta = 1, k = 1, ...)
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
                  belong to the learning set. Must not be missing
                  for this method.}
  \item{beta}{Slope parameter for the logistic function which
              is used for the computation of class probabilities.
              The default value (1) need not produce reasonable
              results and can produce warnings.}
  \item{k}{Number of nearest neighbours to use.}
  \item{\dots}{Currently unused argument.}
}
\details{The algorithm is as follows:
         \itemize{
         \item Determine the \code{k} nearest neighbours
         \item For each class represented among these, compute the average
               euclidean distance.
         \item The negative distances are plugged into the logistic function
               with parameter \code{beta}.
               \item Classify into the class with highest probability.
	     }
}
\value{An object of class \code{\link{cloutput}}.}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}},
         \code{\link{pnnCMA}}, \code{\link{qdaCMA}}, \code{\link{rfCMA}},
         \code{\link{scdaCMA}}, \code{\link{shrinkldaCMA}}, \code{\link{svmCMA}}}
\examples{
### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression from first 10 genes
golubX <- as.matrix(golub[,-1])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run probabilistic k-nearest neighbours
result <- pknnCMA(X=golubX, y=golubY, learnind=learnind, k = 3)
### show results
show(result)
ftable(result)
plot(result)
}
\keyword{multivariate}
