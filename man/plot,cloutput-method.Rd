\name{plot}
\alias{plot,cloutput-method}
\alias{plot,cloutput-method}
\alias{plot,cloutput,missing-method}
\title{Probability plot}
\description{A popular way of visualizing the output of classifier
             is to plot, separately for each class, the predicted
             probability of each predicted observations for the respective class.
             For this purpose, the plot area is divided into \code{K}
             parts, where \code{K} is the number of classes.
             Predicted observations are assigned, according to their
             true class, to one of those parts. Then, for each part
             and each predicted observation, the predicted probabilities
             are plotted, displayed by coloured dots, where each
             colour corresponds to one class.}
             
\arguments{
  \item{x}{An object of class \code{\link{cloutput}} whose slot
           \code{probmatrix} does not contain any missing value,
           i.e. probability estimations are provided by the classifier.}
  \item{main}{A title for the plot (character).}
}

\note{The plot usually only makes sense if a sufficiently large numbers
      of observations has been classified. This is usually achieved
      by running the classifier on several \code{\link{learningsets}}
      with the method \code{\link{classification}}. The output can
      then be processed via \code{\link{join}} to obtain an object
      of class \code{\link{cloutput}} to which this method can be applied.}
\value{No return.}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\code{\link{cloutput}}}
\keyword{multivariate}