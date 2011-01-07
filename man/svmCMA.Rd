%changes !!?
\name{svmCMA}
\alias{svmCMA}
\title{Support Vector Machine}
\description{Calls the function \code{svm} from the package \code{e1071}
that provides an interface to the award-winning LIBSVM routines.
For \code{S4} method information, see \link{svmCMA-methods}
}
\usage{
svmCMA(X, y, f, learnind, probability, models=FALSE,seed=341,...) %!!?
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

  \item{probability}{logical indicating whether the model should allow for
          probability predictions.}
          \item{seed}{Fix random number generator for reproducibility.}
          \item{models}{a logical value indicating whether the model object shall be returned }
  \item{\dots}{Further arguments to be passed to \code{svm} from the
               package \code{e1071}}
               
  
}
\note{Contrary to the default settings in \code{e1071:::svm}, the used
      kernel is a linear kernel which has turned to be out a better
      default setting in the small sample, large number of predictors - situation,
      because additional nonlinearity is mostly not necessary there. It
      additionally avoids the tuning of a further kernel parameter \code{gamma},
      s. help of the package \code{e1071} for details.\cr
      Nevertheless, hyperparameter tuning concerning the parameter \code{cost} must
      usually be performed to obtain reasonale results, s. \code{\link{tune}}.}
\value{An object of class \code{\link{cloutput}}.}
\references{
            Boser, B., Guyon, I., Vapnik, V.  (1992) \cr
            A training algorithm for optimal margin classifiers.\cr
            \emph{Proceedings of the fifth annual workshop on Computational learning theory, pages 144-152, ACM Press}.

            Chang, Chih-Chung and Lin, Chih-Jen :
            LIBSVM: a library for Support Vector Machines
            \url{http://www.csie.ntu.edu.tw/~cjlin/libsvm}
            
            Schoelkopf, B., Smola, A.J. (2002) \cr
            Learning with kernels.
            \emph{MIT Press, Cambridge, MA.}
            }
                    
\author{Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}
        
        Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}}


\seealso{\code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}},
         \code{\link{pnnCMA}}, \code{\link{qdaCMA}}, \code{\link{rfCMA}},
         \code{\link{scdaCMA}}, \code{\link{shrinkldaCMA}}}
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
### run _untuned_linear SVM
svmresult <- svmCMA(X=golubX, y=golubY, learnind=learnind,probability=TRUE)
### show results
show(svmresult)
ftable(svmresult)
plot(svmresult)}
\keyword{multivariate}