\name{pnnCMA}
\alias{pnnCMA}
\title{Probabilistic Neural Networks}
\description{Probabilistic Neural Networks is the term Specht (1990) used
             for a Gaussian kernel estimator for the conditional class
             densities.\cr
             For \code{S4} method information, see \link{pnnCMA-methods}.}
\usage{
pnnCMA(X, y, f, learnind, sigma = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Gene expression data. Can be one of the following:
           \itemize{
           \item A \code{matrix}. Rows correspond to observations, columns to variables.
           \item A \code{data.frame}, when \code{f} is \emph{not} missing (s. below).
           \item An object of class \code{ExpressionSet}.\cr
           Each variable (gene) will be scaled for unit variance and zero mean.
           }}
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
                  belong to the learning set. For this method, this
                  must \emph{not} be \code{missing}.}
  \item{sigma}{Standard deviation of the Gaussian Kernel used.\cr
               This hyperparameter should be tuned, s. \code{\link{tune}}.
               The default is \code{1}, but this generally  does not
               lead to good results. Actually, this method reacts
               very sensitively to the value of sigma. Take care
               if warnings appear related to the particular choice.}
}


\value{An object of class \code{\link{cloutput}}.}
\note{There is actually no strong relation of this method to Feed-Forward
      Neural Networks, s. \code{\link{nnetCMA}}.}
\references{Specht, D.F. (1990).\cr
            Probabilistic Neural Networks.
            \emph{Neural Networks, 3, 109-118}.}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
        

\seealso{\code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}},
         \code{\link{qdaCMA}}, \code{\link{rfCMA}},
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
### run PNN
pnnresult <- pnnCMA(X=golubX, y=golubY, learnind=learnind, sigma = 3)
### show results
show(pnnresult)
ftable(pnnresult)
plot(pnnresult)
}
         
\keyword{multivariate}