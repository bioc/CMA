\name{classification}
\alias{classification}
\title{General method for classification with various methods}
\description{
  Most general function in the package, providing an interface
  to perform variable selection, hyperparameter tuning and
  classification in one step. Alternatively, the first two
  steps can be performed separately and can then be plugged into
  this function.\cr
  For S4 method information, s. \code{\link{classification-methods}}.
}
\usage{
classification(X, y, f, learningsets, genesel, genesellist = list(), nbgene, classifier, tuneres, tuninglist = list(), trace = TRUE, ...)
}
\arguments{
  \item{X}{Gene expression data. Can be one of the following:
           \itemize{
           \item A \code{matrix}. Rows correspond to observations, columns to variables.
           \item A \code{data.frame}, when \code{f} is \emph{not} missing (s. below).
           \item An object of class \code{ExpressionSet}.}
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
           total number of different classes in the learning set.}
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
  \item{tuneres}{Analogous to the argument \code{genesel} - object of
                 class \code{\link{tuningresult}} containing information
                 about the best hyperparameter choice for the argument
                 \code{learningsets}.}
  \item{tuninglist}{Analogous to the argument \code{genesellist}. In the
                    case that the argument \code{tuneres} is missing, this
                    in argument list passed to \code{\link{tune}}.
                    If both \code{tuneres} and \code{tuninglist} are missing,
                    no variable selection is performed.
                    \bold{warning}: Note that if a user-defined hyperparameter
                    grid is passed, this will result in a list within a list:
                    \code{tuninglist = list(grids=list(argname = c())}, s. example.
                    \bold{warning}: Contrary to \code{\link{tune}}, if
                    \code{tuninglist} is an empty list (default),
                    \emph{no} hyperparameter tuning
                    will be performed at all. To use pre-defined hyperparameter
                    grids, the argument is \code{tuninglist = list(grids = list())}.}
  \item{trace}{Should progress be traced ? Default is \code{TRUE}.}
  \item{\dots}{Further arguments passed to the function \code{classifier}.}
}
\details{For details about hyperparameter tuning, consult \code{\link{tune}}.}
\value{A list of objects of class \code{\link{cloutput}} and \code{\link{clvarseloutput}},
       respectively; its length equals the number of different \code{learningsets}.
       The single elements of the list can convenienly be combined using
       the \code{\link{join}} function. The results can be analyzed and
       evaluated by various measures using the method \code{\link{evaluation}}.}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\code{\link{GeneSelection}}, \code{\link{tune}}, \code{\link{evaluation}},
         \code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}},
         \code{\link{pnnCMA}}, \code{\link{qdaCMA}}, \code{\link{rfCMA}},
         \code{\link{scdaCMA}}, \code{\link{shrinkldaCMA}}, \code{\link{svmCMA}}}
\examples{
### a simple k-nearest neighbour example
### datasets
\dontrun{plot(x)
data(golub)
golubY <- golub[,1]
golubX <- as.matrix(golub[,-1])
### learningsets
set.seed(111)
lset <- GenerateLearningsets(y=golubY, method = "CV", fold=5, strat =TRUE)
### 1. GeneSelection
selttest <- GeneSelection(golubX, golubY, learningsets = lset, method = "t.test")
### 2. tuning
tunek <- tune(golubX, golubY, learningsets = lset, genesel = selttest, nbgene = 20, classifier = knnCMA)
### 3. classification
knn1 <- classification(golubX, golubY, learningsets = lset, genesel = selttest,
                       tuneres = tunek, nbgene = 20, classifier = knnCMA)
### steps 1.-3. combined into one step:
knn2 <- classification(golubX, golubY, learningsets = lset,
                       genesellist = list(method  = "t.test"), classifier = knnCMA,
                       tuninglist = list(grids = list(k = c(1:8))), nbgene = 20)
### show and analyze results:
knnjoin <- join(knn2)
show(knn2)
eval <- evaluation(knn2, measure = "misclassification")
show(eval)
summary(eval)
boxplot(eval)
}
}
\keyword{multivariate}
