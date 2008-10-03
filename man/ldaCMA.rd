\name{ldaCMA}
\alias{ldaCMA}
\title{Linear Discriminant Analysis}
\description{Performs a linear discriminant analysis under the assumption
of a multivariate normal distribution in each classes (with equal, but
generally structured) covariance matrices. The function \code{lda} from
the package \code{MASS} is called for computation.

For \code{S4} method information, see \link{ldaCMA-methods}.
}
\usage{
ldaCMA(X, y, f, learnind, ...)
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
  \item{\dots}{Further arguments to be passed to \code{lda} from the
               package \code{MASS}}
}
\note{Excessive variable selection has usually to performed before
      \code{ldaCMA} can be applied in the \code{p > n} setting.
      Not reducing the number of variables can result in an error
      message.}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{An object of class \code{\link{cloutput}}.}
\references{McLachlan, G.J. (1992).

           Discriminant Analysis and Statistical Pattern Recognition.

           \emph{Wiley, New York}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de}

        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{\code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{LassoCMA}}, \code{\link{nnetCMA}},
         \code{\link{pknnCMA}}, \code{\link{plrCMA}}, \code{\link{pls_ldaCMA}},
         \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}}, \code{\link{pnnCMA}},
         \code{\link{qdaCMA}}, \code{\link{rfCMA}}, \code{\link{scdaCMA}},
         \code{\link{shrinkldaCMA}}, \code{\link{svmCMA}}}
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
### run LDA
ldaresult <- ldaCMA(X=golubX, y=golubY, learnind=learnind)
### show results
show(ldaresult)
ftable(ldaresult)
plot(ldaresult)
### multiclass example:
### load Khan data
data(khan)
### extract class labels
khanY <- khan[,1]
### extract gene expression from first 10 genes
khanX <- as.matrix(khan[,2:11])
### select learningset
set.seed(111)
learnind <- sample(length(khanY), size=floor(ratio*length(khanY)))
### run LDA
ldaresult <- ldaCMA(X=khanX, y=khanY, learnind=learnind)
### show results
show(ldaresult)
ftable(ldaresult)
plot(ldaresult)
}
\keyword{multivariate}
