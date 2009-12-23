\name{predoutput-class}
\docType{class}
\alias{predoutput-class}
\alias{predoutput}
\alias{show,predoutput-method}

\title{"predoutput"}
\description{Object returned by the function \code{prediction}}
\section{Slots}{
	 \describe{
    \item{\code{Xnew}:}{Gene Expression matrix of new observations}
    \item{\code{yhat}:}{Predicted class labels for the new data.}
    \item{\code{model}:}{List containing the constructed classifier.}
 
  }
}
\section{Methods}{
  \describe{
    \item{show}{Returns predicted class labels for the new data.}
 
	 }
}

\author{Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}}

\seealso{
        \code{\link{compBoostCMA}}, \code{\link{dldaCMA}}, \code{\link{ElasticNetCMA}},
         \code{\link{fdaCMA}}, \code{\link{flexdaCMA}}, \code{\link{gbmCMA}},
         \code{\link{knnCMA}}, \code{\link{ldaCMA}}, \code{\link{LassoCMA}},
         \code{\link{nnetCMA}}, \code{\link{pknnCMA}}, \code{\link{plrCMA}},
         \code{\link{pls_ldaCMA}}, \code{\link{pls_lrCMA}}, \code{\link{pls_rfCMA}},
         \code{\link{pnnCMA}}, \code{\link{qdaCMA}}, \code{\link{rfCMA}},
         \code{\link{scdaCMA}}, \code{\link{shrinkldaCMA}}, \code{\link{svmCMA}}}
\keyword{multivariate}