\name{shrinkldaCMA}
\alias{shrinkldaCMA}
\title{Shrinkage linear discriminant analysis}
\description{Linear Discriminant Analysis combined with the
             James-Stein-Shrinkage approach of Schaefer and Strimmer (2005)
             for the covariance matrix.

             Currently still an experimental version.

             For \code{S4} method information, see \link{shrinkldaCMA-methods}}
\usage{
shrinkldaCMA(X, y, f, learnind, models=FALSE, ...)
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
  \item{\dots}{Further arguments to be passed to \code{cov.shrink} from the
               package \code{corpcor}}
 }

\value{An object of class \code{\link{cloutput}}.}

\note{This is still an experimental version.

      Covariance shrinkage is performed by calling functions
      from the package \code{corpcor}.

      Variable selection is \emph{not} necessary.}

\references{ Schaefer, J., Strimmer, K. (2005).

             A shrinkage approach to large-scale covariance estimation and implications for functional genomics.

             \emph{Statististical Applications in Genetics and Molecular Biology, 4:32. }
             }
             
\author{Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}}


\seealso{\code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}},
         \code{\link{pnnCMA}}, \code{\link{qdaCMA}}, \code{\link{rfCMA}},
         \code{\link{scdaCMA}}, \code{\link{svmCMA}}.}
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
### run  shrinkage-LDA
result <- shrinkldaCMA(X=golubX, y=golubY, learnind=learnind)
### show results
show(result)
ftable(result)
plot(result)}
\keyword{multivariate}
