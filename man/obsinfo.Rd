\name{obsinfo}
\alias{obsinfo}
\title{Classifiability of observations}
\description{
  Some observations are harder to classify than others.
  It is frequently of interest to know which observations
  are consistenly misclassified; these are candiates for
  outliers or wrong class labels.
}

\arguments{
  \item{object}{An object of class \code{\link{evaluation}}, generated
                with \code{scheme = "observationwise"}}
  \item{threshold}{threshold value of (observation-wise) performance measure,
                   s. \code{\link{evaluation}} that has to be exceeded
                   in order to speak of consistent misclassification.
                   If \code{measure = "average probability"}, then
                   values \emph{below} \code{threshold} are regarded
                   as consistent misclassification. Note that the
                   default values 1 is not sensible in that case}
  \item{show}{Should the information be printed ? Default is \code{TRUE}.}
}
\details{As not all observation must have been classified at least once,
         observations not classified at all are also shown.}
\value{A list with two components
  \item{misclassification}{A \code{data.frame} containing the indices
                           of consistenly misclassfied observations and
                           the corresponding performance measure.}
  \item{notclassified}{The indices of those observations not classfied at all,
                       s. details.}
}

\author{Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}}
\seealso{\code{\link{evaluation}}}
\keyword{multivariate}
