\name{Planarplot}
\alias{Planarplot}
\title{Visualize Separability of different classes}
\description{
  Given two variables, the methods trains a classifier
  (argument \code{classfier}) based on these two variables
  and plots the resulting class regions, learning- and test
  observations in the plane.

  Appropriate variables are usually found by \code{\link{GeneSelection}}.

  For S4 method information, s. \code{\link{Planarplot-methods}}.
}
\usage{
Planarplot(X, y, f, learnind, predind, classifier, gridsize = 100, ...)
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
           }
  \item{f}{A two-sided formula, if \code{X} is a \code{data.frame}. The
           left part correspond to class labels, the right to variables.}
  \item{learnind}{An index vector specifying the observations that
                  belong to the learning set. May be \code{missing};
                  in that case, the learning set consists of all
                  observations and predictions are made on the
                  learning set.}
  \item{predind}{A vector containing \emph{exactly} two indices that
                 denote the two variables used for classification.}
  \item{classifier}{Name of function ending with \code{CMA} indicating
                    the classifier to be used.}
  \item{gridsize}{The gridsize used for two-dimensional plotting.

                  For both variables specified in \code{predind},
                  an equidistant grid of size \code{gridsize} is created.
                  The resulting two grids are then combined to obtain
                  \code{gridsize^2} points in the real plane which are
                  used to draw the class regions. Defaults to 100
                  which is usually a reasonable choice, but takes
                  some time.}
  \item{\dots}{Further argument passed to \code{classifier}.}
}

\value{No return.}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de}

        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}.
        Idea is from the \code{MLInterfaces} package, contributed
        by Jess Mar, Robert Gentleman and Vince Carey.}
\seealso{\code{GeneSelection},
         \code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}},
         \code{\link{pnnCMA}}, \code{\link{qdaCMA}}, \code{\link{rfCMA}},
         \code{\link{scdaCMA}}, \code{\link{shrinkldaCMA}}, \code{\link{svmCMA}}}
\examples{
### simple linear discrimination for the golub data:
data(golub)
golubY <- golub[,1]
golubX <- as.matrix(golub[,-1])
golubn <- nrow(golubX)
set.seed(111)
learnind <- sample(golubn, size=floor(2/3*golubn))
Planarplot(X=golubX, y=golubY, learnind=learnind, predind=c(2,4),
           classifier=ldaCMA)

}
\keyword{multivariate}
