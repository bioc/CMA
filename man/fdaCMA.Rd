\name{fdaCMA}
\alias{fdaCMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fisher's Linear Discriminant Analysis}
\description{
 Fisher's Linear Discriminant Analysis constructs a subspace of
 'optimal projections' in which classification is performed.
 The directions of optimal projections are computed by the
 function \code{cancor} from the package \code{stats}. For
 an exhaustive treatment, see e.g. Ripley (1996).

 For \code{S4} method information, see \code{\link{fdaCMA-methods}.}
}
\usage{
fdaCMA(X, y, f, learnind, comp = 1, plot = FALSE,models=FALSE)
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
           \bold{WARNING}: The class labels will be re-coded to
           range from \code{0} to \code{K-1}, where \code{K} is the
           total number of different classes in the learning set.
	   }
           }
  \item{f}{A two-sided formula, if \code{X} is a \code{data.frame}. The
           left part correspond to class labels, the right to variables.}
  \item{learnind}{An index vector specifying the observations that
                  belong to the learning set. May be \code{missing};
                  in that case, the learning set consists of all
                  observations and predictions are made on the
                  learning set.}
  \item{comp}{Number of discriminant coordinates (projections) to compute.
              Default is one, must be smaller than or equal to \code{K-1}, where
              \code{K} is the number of classes.}
  \item{plot}{Should the projections onto the space spanned by the optimal
              projection directions be plotted ? Default is \code{FALSE}.}
              \item{models}{a logical value indicating whether the model object shall be returned }
}

\note{Excessive variable selection has usually to performed before
      \code{fdaCMA} can be applied in the \code{p > n} setting.
      Not reducing the number of variables can result in an error
      message.}

\value{An object of class \code{\link{cloutput}}.}

\references{Ripley, B.D. (1996)

            Pattern Recognition and Neural Networks.

            \emph{Cambridge University Press}}

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
### run FDA
fdaresult <- fdaCMA(X=golubX, y=golubY, learnind=learnind, comp = 1, plot = TRUE)
### show results
show(fdaresult)
ftable(fdaresult)
plot(fdaresult)
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
### run FDA
fdaresult <- fdaCMA(X=khanX, y=khanY, learnind=learnind, comp = 2, plot = TRUE)
### show results
show(fdaresult)
ftable(fdaresult)
plot(fdaresult)
}
\keyword{multivariate}