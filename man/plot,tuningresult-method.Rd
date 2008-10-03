\name{plot}
\alias{plot,tuningresult-method}
\alias{plot,tuningresult,missing-method}
\title{Visualize results of tuning}
\description{
After hyperparameter tuning using \code{\link{tune}} it is useful
to see which choice of hyperparameters is suitable and how good the
performance is.}
\arguments{
  \item{x}{An object of class \code{\link{tuningresult}}.}
  \item{iter}{Iteration number (\code{learningset}) for which tuning
              results should be displayed.}
  \item{which}{Character vector (maximum length is two) naming
               the arguments for which tuning results should
               be display. Default is \code{NULL}; if the
               number of tuned hyperparameter is less or equal
               than two, then the results for these hyperparameters
               will be plotted. If this number is two, then a \code{contour}
               plot will be made, otherwise a simple line segment plot.
               If the number of tuned hyperparameters exceeds two, then
               which may not be \code{NULL}.}
  \item{\dots}{Further graphical options passed either to \code{plot}
               or \code{contour}.}
}
\note{Frequently, several hyperparameter (combinations) perform "best",
      s. also the remark in \code{\link{best}}.}
\value{no return.}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de}

        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\code{\link{tune}}, \code{\link{tuningresult}}}
\keyword{multivariate}
