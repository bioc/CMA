\name{plrCMA}
\alias{plrCMA}
\title{L2 penalized logistic regression}
\description{High dimensional logistic regression combined with an
             L2-type (Ridge-)penalty.
             Multiclass case is also possible.
             For \code{S4} method information, see \link{plrCMA-methods}}
\usage{
plrCMA(X, y, f, learnind, lambda = 0.01, scale = TRUE, ...)
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
  \item{lambda}{Parameter governing the amount of penalization.
                This hyperparameter should be \code{\link{tune}}d.}
  \item{ scale }{Scale the predictors as specified by \code{X} to have unit variance
                 and zero mean.}
  \item{\dots}{Currently unused argument.}
}

\value{An object of class \code{\link{cloutput}}.}
\references{Zhu, J., Hastie, T. (2004).
            Classification of gene microarrays by penalized logistic regression.\cr
            \emph{Biostatistics 5:427-443}.}
\author{Special thanks go to \cr
        Ji Zhu (University of Ann Arbor, Michigan)
 \cr
        Trevor Hastie (Stanford University)
 \cr
        who provided the basic code that was then adapted by

        Martin Slawski \email{martin.slawski@campus.lmu.de}, \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}.}

\seealso{code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}},
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
### run penalized logistic regression (no tuning)
plrresult <- plrCMA(X=golubX, y=golubY, learnind=learnind)
### show results
show(plrresult)
ftable(plrresult)
plot(plrresult)
### multiclass example:
### load Khan data
data(khan)
### extract class labels
khanY <- khan[,1]
### extract gene expression from first 10 genes
khanX <- as.matrix(khan[,-1])
### select learningset
set.seed(111)
learnind <- sample(length(khanY), size=floor(ratio*length(khanY)))
### run penalized logistic regression (no tuning)
plrresult <- plrCMA(X=khanX, y=khanY, learnind=learnind)
### show results
show(plrresult)
ftable(plrresult)
plot(plrresult)
}
\keyword{multivariate}
