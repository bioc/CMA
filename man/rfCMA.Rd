\name{rfCMA}
\alias{rfCMA}
\title{Classification based on Random Forests}
\description{Random Forests were proposed by Breiman (2001)
             and are implemented in the package \code{randomForest}.

             In this package, they can as well be used to rank variables
             according to their importance, s. \code{GeneSelection}.

             For \code{S4} method information, see \link{rfCMA-methods}}
\usage{
rfCMA(X, y, f, learnind, varimp = TRUE, seed = 111, models=FALSE,type=1,scale=FALSE,importance=TRUE, ...)}
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
  \item{varimp}{Should additional information for variable selection be provided ? Defaults to \code{TRUE}.}
  \item{seed}{Fix Random number generator seed to \code{seed}. This is
              useful to guarantee reproducibility of the results.}
              \item{models}{a logical value indicating whether the model object shall be returned }
 \item{type}{Parameter passed to function \code{importance}. Either 1 or 2, specifying the type of importance measure
          (1=mean decrease in accuracy, 2=mean decrease in node
          impurity).}
 \item{scale}{Parameter passed to function \code{importance}. For permutation based measures, should the measures be
          divided by their standard errors?}
 \item{importance}{Parameter passed to function \code{randomForest}.Should importance of predictors be assessed by permutation?}
  \item{\dots}{Further arguments to be passed to \code{randomForest} from the
               package of the same name.}
}


\value{If \code{varimp}, then an object of class \code{\link{clvarseloutput}} is returned,
       otherwise an object of class \code{\link{cloutput}}}
\references{Breiman, L. (2001)

            Random Forest.

            \emph{Machine Learning, 45:5-32.}}

\author{Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}}

\seealso{\code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}},
         \code{\link{pnnCMA}}, \code{\link{qdaCMA}},
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
### run random Forest
#rfresult <- rfCMA(X=khanX, y=khanY, learnind=learnind, varimp = FALSE)
### show results
#show(rfresult)
#ftable(rfresult)
#plot(rfresult)}


\keyword{multivariate}