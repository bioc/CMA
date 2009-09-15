\name{pls_ldaCMA}
\alias{pls_ldaCMA}
\title{Partial Least Squares combined with Linear Discriminant Analysis}
\description{
  This method constructs a classifier that extracts
  Partial Least Squares components that are plugged into
  Linear Discriminant Analysis.
  The Partial Least Squares components are computed by the package
  \code{plsgenomics}.

  For \code{S4} method information, see \code{\link{pls_ldaCMA-methods}}.
}
\usage{
pls_ldaCMA(X, y, f, learnind, comp = 2, plot = FALSE)
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
  \item{plot}{If \code{comp <= 2}, should the classification space of the
              Partial Least Squares components be plotted ? Default is \code{FALSE}.}
}

\value{An object of class \code{\link{cloutput}}.}
\references{Nguyen, D., Rocke, D. M., (2002).

            Tumor classifcation by partial least squares using microarray gene expression data.

            \emph{Bioinformatics 18, 39-50}
            
            Boulesteix, A.L.,  Strimmer, K. (2007).

            Partial least squares: a versatile tool for the analysis of high-dimensional genomic data.

            \emph{Briefings in Bioinformatics 7:32-44.}
            }

\author{Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}}

\seealso{\code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}},
         \code{\link{pnnCMA}}, \code{\link{qdaCMA}}, \code{\link{rfCMA}},
         \code{\link{scdaCMA}}, \code{\link{shrinkldaCMA}}, \code{\link{svmCMA}}}

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
plsresult <- pls_ldaCMA(X=khanX, y=khanY, learnind=learnind, comp = 4)
### show results
show(plsresult)
ftable(plsresult)
plot(plsresult)}

\keyword{multivariate}
