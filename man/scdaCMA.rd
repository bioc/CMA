\name{scdaCMA}
\alias{scdaCMA}
\title{Shrunken Centroids Discriminant Analysis}
\description{
The nearest shrunken centroid classification algorithm is
detailly described in Tibshirani et al. (2002).

It is widely known under the name PAM (prediction analysis for microarrays),
which can also be found in the package \code{pamr}.

For \code{S4} method information, see \link{scdaCMA-methods}.
}
\usage{
scdaCMA(X, y, f, learnind, delta = 0.5, ...)
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
  \item{delta}{The shrinkage intensity for the class centroids -
               a hyperparameter that must be tuned. The default
               \code{0.5} not necessarily produces good results.}
  \item{\dots}{Currently unused argument.}
}



\value{An object of class \code{\link{cloutput}}.}

\note{The results can differ from those obtained by
      using the package \code{pamr}.}

\references{ Tibshirani, R., Hastie, T., Narasimhan, B., and Chu, G., (2003).

             Class prediction by nearest shrunken centroids with applications to DNA microarrays.

             \emph{Statistical Science, 18, 104-117}}
             
\author{Martin Slawski \email{martin.slawski@campus.lmu.de}

        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
        

\seealso{\code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}},
         \code{\link{pnnCMA}}, \code{\link{qdaCMA}}, \code{\link{rfCMA}},
         \code{\link{shrinkldaCMA}}, \code{\link{svmCMA}}}

\examples{
### load Khan data
data(khan)
### extract class labels
khanY <- khan[,1]
### extract gene expression
khanX <- as.matrix(khan[,-1])
### select learningset
set.seed(111)
learnind <- sample(length(khanY), size=floor(2/3*length(khanY)))
### run Shrunken Centroids classfier, without tuning
scdaresult <- scdaCMA(X=khanX, y=khanY, learnind=learnind)
### show results
show(scdaresult)
ftable(scdaresult)
plot(scdaresult)}
\keyword{multivariate}