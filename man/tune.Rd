\name{tune}
\alias{tune}
\title{Hyperparameter tuning for classifiers}
\description{Most classifiers implemented in this package depend on one
or even several hyperparameters (s. details) that should be optimized
to obtain good (and comparable !) results. As tuning scheme, we propose
three fold Cross-Validation on each \code{learningset} (for fixed selected
variables). Note that \code{learningsets} usually do not contain the
complete dataset, so tuning involves a second level of splitting the dataset.
Increasing the number of folds leads to larger datasets (and possibly to higher accuracy),
but also to higher computing times.\cr
For S4 method information, s. \code{link{tune-methods}}}
\usage{
tune(X, y, f, learningsets, genesel, genesellist = list(), nbgene, classifier, fold = 3, strat = FALSE, grids = list(), trace = TRUE, ...)
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
     }
  \item{f}{A two-sided formula, if \code{X} is a \code{data.frame}. The
           left part correspond to class labels, the right to variables.}
  \item{learningsets}{An object of class \code{\link{learningsets}}. May
                      be missing, then the complete datasets is used as
                      learning set.}
  \item{genesel}{Optional (but usually recommended) object of class
                 \code{\link{genesel}} containing variable importance
                 information for the argument \code{learningsets}}
  \item{genesellist}{In the case that the argument \code{genesel} is missing,
                     this is an argument list passed to \code{\link{GeneSelection}}.
                     If both \code{genesel} and \code{genesellist} are missing,
                     no variable selection is performed.}
  \item{nbgene}{Number of best genes to be kept for classification, based
                on either \code{genesel} or the call to \code{\link{GeneSelection}}
                using \code{genesellist}. In the case that both are missing,
                this argument is not necessary.
                \bold{note}:
		\itemize{
  
                \item If the gene selection method has been one of
                      \code{"lasso", "elasticnet", "boosting"},
                      \code{nbgene} will be reset to \code{min(s, nbgene)}
                      where \code{s} is the number of nonzero coefficients.
                \item if the gene selection scheme has been \code{"one-vs-all", "pairwise"}
                      for the multiclass case, there exist several rankings.
                      The top \code{nbgene} will be kept of \emph{each} of them,
                      so the number of effective used genes will sometimes be much
                      larger.
		    }
		}
  \item{classifier}{Name of function ending with \code{CMA} indicating
                    the classifier to be used.}
  \item{fold}{The number of cross-validation folds used within each \code{learningset}.
              Default is 3. Increasing \code{fold} will lead to higher computing times.}
  \item{strat}{Should stratified cross-validation according to the class proportions
               in the complete dataset be  used ? Default is \code{FALSE}.}
  \item{grids}{A named list. The names correspond to the arguments to be tuned,
               e.g. \code{k} (the number of nearest neighbours) for \code{\link{knnCMA}},
               or \code{cost} for \code{\link{svmCMA}}. Each element is a \code{numeric}
               vector defining the grid of candidate values. Of course, several hyperparameters
               can be tuned simultaneously (though requiring much time). By
               default, \code{grids} is an empty list. In that case, a pre-defined
               list will be used, s. details.}
  \item{trace}{Should progress be traced ? Default is \code{TRUE}.}
  \item{\dots}{Further arguments to be passed to \code{classifier}, of course
               not one of the arguments to be tuned (!).}
}
\note{The computation time can be enormously high. Note that for each different
      \code{learningset}, the classifier must be trained \code{fold} times
      \code{number of possible different hyperparameter combinations} times.
      E.g. if the number of the learningsets is fifty, \code{fold = 3} and
      two hyperparameters (each with 5 candidate values) are tuned, 50x3x25=3750
      training iterations are necessary !}
\details{The following default settings are used, if the arguments \code{grids}
         is an empty list:
         \describe{
         \item{\code{\link{gbmCMA}}}{\code{n.trees = c(50, 100, 200, 500, 1000)}}
         \item{\code{\link{compBoostCMA}}}{\code{mstop = c(50, 100, 200, 500, 1000)}}
         \item{\code{\link{LassoCMA}}}{\code{norm.fraction = seq(from=0.1, to=0.9, length=9)}}
         \item{\code{\link{ElasticNetCMA}}}{\code{norm.fraction = seq(from=0.1, to=0.9, length=5), lambda2 = 2^{-(5:1)}}}
         \item{\code{\link{plrCMA}}}{\code{lambda = 2^{-4:4}}}
         \item{\code{\link{pls_ldaCMA}}}{\code{comp = 1:10}}
         \item{\code{\link{pls_lrCMA}}}{\code{comp = 1:10}}
         \item{\code{\link{pls_rfCMA}}}{\code{comp = 1:10}}
         \item{\code{\link{rfCMA}}}{\code{mtry = ceiling(c(0.1, 0.25, 0.5, 1, 2)*sqrt(ncol(X))), nodesize = c(1,2,3)}}
         \item{\code{\link{knnCMA}}}{\code{k=1:10}}
         \item{\code{\link{pknnCMA}}}{\code{k = 1:10}}
         \item{\code{\link{scdaCMA}}}{\code{delta = c(0.1, 0.25, 0.5, 1, 2, 5)}}
         \item{\code{\link{pnnCMA}}}{\code{sigma = c(2^{-2:2})}},
         \item{\code{\link{nnetCMA}}}{\code{size = 1:5, decay = c(0, 2^{-(4:1)})}}
         \item{\code{\link{svmCMA}}, \code{kernel = "linear"}}{\code{cost = c(0.1, 1, 5, 10, 50, 100, 500)}}
         \item{\code{\link{svmCMA}}, \code{kernel = "radial"}}{\code{cost = c(0.1, 1, 5, 10, 50, 100, 500), gamma = 2^{-2:2}}}

         \item{\code{\link{svmCMA}}, \code{kernel = "polynomial"}}{\code{cost = c(0.1, 1, 5, 10, 50, 100, 500), degree = 2:4}}
         }
}
\value{An object of class \code{\link{tuningresult}}}

\references{
Slawski, M. Daumer, M.  Boulesteix, A.-L. (2008)
CMA - A comprehensive Bioconductor package for supervised classification with high dimensional data.
\emph{BMC Bioinformatics 9: 439}

}

\author{Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}
        
        Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}}

\seealso{\code{\link{tuningresult}}, \code{\link{GeneSelection}}, \code{\link{classification}}}
\examples{
\dontrun{
### simple example for a one-dimensional grid, using compBoostCMA.
### dataset
data(golub)
golubY <- golub[,1]
golubX <- as.matrix(golub[,-1])
### learningsets
set.seed(111)
lset <- GenerateLearningsets(y=golubY, method = "CV", fold=5, strat =TRUE)
### tuning after gene selection with the t.test
tuneres <- tune(X = golubX, y = golubY, learningsets = lset,
              genesellist = list(method = "t.test"),
              classifier=compBoostCMA, nbgene = 100,
              grids = list(mstop = c(50, 100, 250, 500, 1000)))
### inspect results
show(tuneres)
best(tuneres)
plot(tuneres, iter = 3)
}
}
\keyword{multivariate}
