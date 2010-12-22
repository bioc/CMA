\name{pls_rfCMA}
\alias{pls_rfCMA}
\title{Partial Least Squares followed by random forests}
\description{
  This method constructs a classifier that extracts
  Partial Least Squares components used to generate Random Forests, s. \code{\link{rfCMA}}.

  For \code{S4} method information, see \code{\link{pls_rfCMA-methods}}.
}
\usage{
pls_rfCMA(X, y, f, learnind, comp = 2 * nlevels(as.factor(y)), seed = 111,models=FALSE, ...)
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
  \item{comp}{Number of Partial Least Squares components to extract.
              Default ist two times the number of different classes.}
  \item{seed}{Fix Random number generator seed to \code{seed}. This is
              useful to guarantee reproducibility of the results,
              due to the random component in the random Forest.}
              \item{models}{a logical value indicating whether the model object shall be returned }
  \item{\dots}{Further arguments to be passed to \code{randomForests} from the
               package of the same name.}
}

\value{An object of class \code{\link{cloutput}}. }

\references{Boulesteix, A.L.,  Strimmer, K. (2007).

            Partial least squares: a versatile tool for the analysis of high-dimensional genomic data.

            \emph{Briefings in Bioinformatics 7:32-44.}}

\author{Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}}

\seealso{\code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}},
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
### run PLS, combined with Random Forest
#result <- pls_rfCMA(X=golubX, y=golubY, learnind=learnind)
### show results
#show(result)
#ftable(result)
#plot(result)
}

\keyword{multivariate}
