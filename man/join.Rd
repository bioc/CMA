\name{join}
\alias{join}
\title{Combine list elements returned by the method classification}
\description{The method \code{classification}  returns a list of
             class \code{\link{cloutput}} or \code{\link{clvarseloutput}}.
             It is often more convenient to work with an object of class
             \code{\link{cloutput}} instead with a whole list, e.g.
             because the convenience method defined for that class can
             be used.\cr
             For S4 method information, s. \code{\link{join-methods}}}
\usage{
join(cloutputlist)
}
\arguments{
  \item{cloutputlist}{A list of objects of classes \code{\link{cloutput}}
                      or \code{\link{clvarseloutput}}, usually that
                      returned by a call to the method \code{\link{classification}}.
                      The only requirement for a succesful join is that the
                      used dataset and classfier are the same for each
                      list element.}
}
\value{An object of class \code{\link{cloutput}}.
        \bold{warning:}If the elements of
       \code{cloutputlist} have originally been of class \code{\link{clvarseloutput}},
       the slot \code{varsel} will be dropped !}


\note{The result of the \code{join} method is incompatible with the methods
      \code{\link{evaluation}}, \code{\link{compare}}. These require the lists returned by
      \code{\link{classification}}.}
\seealso{\code{\link{classification}}, \code{\link{evaluation}}}
\keyword{multivariate}
