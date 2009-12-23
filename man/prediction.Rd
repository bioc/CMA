\name{prediction}
\alias{prediction}
\title{General method for predicting classes of new observations}
\description{
  This method constructs the given classifier using the specified training data, gene selection and tuning results.. Subsequently, class labels are predicted for new observations.\cr
  For S4 method information, s. \code{\link{classification-methods}}.
}
\usage{
prediction(X.tr,y.tr,X.new,f,classifier,genesel,models=F,nbgene,tuneres,...)
}
\arguments{
  \item{X.tr}{Training gene expression data. Can be one of the following:
           \itemize{
           \item A \code{matrix}. Rows correspond to observations, columns to variables.
           \item A \code{data.frame}, when \code{f} is \emph{not} missing (s. below).
           \item An object of class \code{ExpressionSet}.}
       }
       
  \item{X.new}{ gene expression data. Can be one of the following:
           \itemize{
           \item A \code{matrix}. Rows correspond to observations, columns to variables.
           \item A \code{data.frame}, when \code{f} is \emph{not} missing (s. below).
           \item An object of class \code{ExpressionSet}.}}
  \item{y.tr}{Class labels of training observation. Can be one of the following:
           \itemize{
           \item A \code{numeric} vector.
           \item A \code{factor}.
           \item A \code{character} if \code{X} is an \code{ExpressionSet} that
                 specifies the phenotype variable.
           \item \code{missing}, if \code{X} is a \code{data.frame} and a
           proper formula \code{f} is provided.
	 }
           \bold{WARNING}: The class labels will be re-coded for classifier construction to
           range from \code{0} to \code{K-1}, where \code{K} is the
           total number of different classes in the learning set. 
           }
  \item{f}{A two-sided formula, if \code{X} is a \code{data.frame}. The
           left part correspond to class labels, the right to variables.}
  \item{genesel}{Optional (but usually recommended) object of class
                 \code{\link{genesel}} containing variable importance
                 information for the argument \code{learningsets}.
                 In this case the object contains a single variable selection. Appropriate 
                 \code{\link{genesel}}-objects can be obtained using the function \code{\link{genesel}}
                 without learningset and setting \code{X=X.tr} and \code{y=y.tr} (i.e. corresponding to the training data of this function).
                 }
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
                 \code{learningsets}. Appropriate tuning-objects can be obtained using the function 
                 \code{\link{tune}} without learningsets and setting parameters \code{X=X.tr}, \code{y=y.tr} and \code{genesel=genesel}
                 (i.e. using the same training data and gene selection as in this function)}
  

  \item{models}{a logical value indicating whether the model object shall be returned }
  \item{\dots}{Further arguments passed to the function \code{classifier}.}
}
\details{This function builds the specified classifier and predicts the class labels of new observations. Hence,
its usage differs from those of most other prediction functions in R.}
\value{A object of class \code{\link{predoutput-class}}; Predicted classes can be seen by \code{show(predoutput)}}
       
\references{
Slawski, M. Daumer, M.  Boulesteix, A.-L. (2008)
CMA - A comprehensive Bioconductor package for supervised classification with high dimensional data.
\emph{BMC Bioinformatics 9: 439} }

\author{Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}
        
        }

\seealso{\code{\link{GeneSelection}}, \code{\link{tune}}, \code{\link{evaluation}},
         \code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}},
         \code{\link{pnnCMA}}, \code{\link{qdaCMA}}, \code{\link{rfCMA}},
         \code{\link{scdaCMA}}, \code{\link{shrinkldaCMA}}, \code{\link{svmCMA}}\code{\link{classification}}}
\examples{
### a simple k-nearest neighbour example
### datasets
\dontrun{plot(x)
data(golub)
golubY <- golub[,1]
golubX <- as.matrix(golub[,-1])
###Splitting data into training and test set
X.tr<-golubX[1:30]
X.new<-golubX[31:39]
y.tr<-golubY[1:30]
### 1. GeneSelection
selttest <- GeneSelection(X=X.tr, y=y.tr, method = "t.test")
### 2. tuning
tunek <- tune(X.tr, y.tr, genesel = selttest, nbgene = 20, classifier = knnCMA)
### 3. classification
pred <- prediction(X.tr=X.tr,y.tr=y.tr,X.new=X.new, genesel = selttest,
                       tuneres = tunek, nbgene = 20, classifier = knnCMA)
### show and analyze results:
show(pred)

}
}
\keyword{multivariate}
