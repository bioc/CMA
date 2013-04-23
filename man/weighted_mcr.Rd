\name{weighted.mcr}
\alias{weighted.mcr}
\title{Tuning / Selection bias correction}
\description{Performs subsampling for several classifiers or a single
  classifiers with different tuning parameter values or numbers of
  selected genes. Eventually, a specific procedure for correcting for
  the tuning or selection bias, which is caused by optimal selection of
  classifiers or tuning parameters, is applied. 
 }
\usage{
weighted.mcr(classifiers,parameters,nbgenes,sel.method,X,y,portion,niter=100,shrinkage=F)
}
\arguments{
  \item{classifiers}{A character vector of the several CMA classifiers
    that shall be used. If the same classifier shall be used with
    different tuning parameters it must appear several times in this
    vector.}

  \item{parameters}{A character containing the tuning parameter values
    corresponding to the classification methods in \code{classifiers}.
   Must have the same length as \code{classifiers}.}
  \item{nbgenes}{A numeric vector indicating how many variables
    shall be selected by \code{sel.method} for the corresponding
    classifier. 
    Must have the same length as \code{classifiers}.}
  \item{sel.method}{The CMA-method (represented as a string) that shall be applied for variable
    selection. If this parameter is set to \code{'none'} no variable
    selection is performed.}
  \item{X}{The  matrix of gene expression data. Can be one of the
    following. Rows correspond to observations, columns to variables.}
  
  \item{y}{Class labels. Can be one of the following:
           \itemize{
           \item A \code{numeric} vector.
           \item A \code{factor}.
          }
           \bold{WARNING}: The class labels will be re-coded to
           range from \code{0} to \code{K-1}, where \code{K} is the
           total number of different classes in the learning set.}
  \item{portion}{A numeric value which indicates the portion of
    observations that will be used for training the classifiers.}
  \item{niter}{The number of subsampling iterations.}
  \item{shrinkage}{A logical value indicating whether shrinkage (WMCS) shall be applied.}
}
\details{The algorithm tries to avoid the additional computational costs
           of a nested cross validation by estimating the corrected
           misclassification rate of the best classifier by a weighted mean
           of all classifiers included in the subsampling approach.}
\value{An object of class \code{wmcr.result} which provides the
    corrected and uncorrected misclassification rate of the best
    classifier as well as weights and misclassifcation rates for all
    classifiers used in the subsampling approach.}
       
\references{Bernau Ch., Augustin, Th. and Boulesteix, A.-L. (2011): Correcting the optimally selected resampling-based error rate: A smooth analytical alternative to nested cross-validation. Department of Statistics: Technical Reports, Nr. 105.}

\author{Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}
        
       }

\seealso{\code{\link{wmc}},\code{\link{classification}},\code{\link{GeneSelection}}, \code{\link{tune}}, \code{\link{evaluation}},
        }
\examples{
#inputs
classifiers<-rep('knnCMA',7)
nbgenes<-rep(50,7)
parameters<-c('k=1','k=3','k=5','k=7','k=9','k=11','k=13')
portion<-0.8
niter<-100
data(golub)
X<-as.matrix(golub[,-1])         
y<-golub[,1]
sel.method<-'t.test'
#function call
wmcr<-weighted.mcr(classifiers=classifiers,parameters=parameters,nbgenes=nbgenes,sel.method=sel.method,X=X,y=y,portion=portion,niter=niter)
}
\keyword{tuning bias, selection bias, corrected misclassification rate}
