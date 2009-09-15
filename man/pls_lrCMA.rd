\name{pls_lrCMA}
\alias{pls_lrCMA}
\title{Partial Least Squares followed by logistic regression}
\description{
  This method constructs a classifier that extracts
  Partial Least Squares components that form the the covariates
  in a binary logistic regression model.
  The Partial Least Squares components are computed by the package
  \code{plsgenomics}.

  For \code{S4} method information, see \code{\link{pls_lrCMA-methods}}.
}
\usage{
pls_lrCMA(X, y, f, learnind, comp = 2, lambda = 1e-4, plot = FALSE)
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
  \item{comp}{Number of Partial Least Squares components to extract.
              Default is 2 which can be suboptimal, depending on the
              particular dataset. Can be optimized using \code{\link{tune}}.}
  \item{lambda}{Parameter controlling the amount of L2 penalization for logistic
                regression, usually taken to be a small value in order to
                stabilize estimation in the case of separable data.}
  \item{plot}{If \code{comp <= 2}, should the classification space of the
              Partial Least Squares components be plotted ? Default is \code{FALSE}.}
}

\value{An object of class \code{\link{cloutput}}.}

\note{Up to now, only the two-class case is supported.}

\references{Boulesteix, A.L.,  Strimmer, K. (2007).

            Partial least squares: a versatile tool for the analysis of high-dimensional genomic data.

            \emph{Briefings in Bioinformatics 7:32-44.}}

\author{Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}}


\seealso{\code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_rfCMA}},
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
### run PLS, combined with logistic regression
result <- pls_lrCMA(X=golubX, y=golubY, learnind=learnind)
### show results
show(result)
ftable(result)
plot(result)
}
\keyword{multivariate}