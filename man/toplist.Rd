\name{toplist}
\alias{toplist}
\alias{toplist,genesel-method}
\title{Display  'top' variables}
\description{This is a convenient method to get quick
             access to the most important variables, based
             on the result of call to \code{\link{GeneSelection}}.}
\usage{
toplist(object, k = 10, iter = 1, show = TRUE, ...)
}
\arguments{
  \item{object}{An object of \code{\link{genesel}}.}
  \item{k}{Number of top genes for which information should be displayed.
             Defaults to 10.}
  \item{iter}{teration number (\code{learningset}) for which tuning
              results should be displayed.}
  \item{show}{Should the results be printed ? Default is \code{TRUE}.}
  \item{\dots}{Currently unused argument.}
}

\value{The type of output depends on the gene selection scheme. For
       the multiclass case, if gene selection has been run with
       the \code{"pairwise"} or \code{"one-vs-all"} scheme, then the
       output will be a list of \code{data.frames}, each containing
       the gene indices plus variable importance for the top \code{k}
       genes. The list elements are named according to the binary
       scenarios (e.g., \code{1 vs. 3}).
       Otherwise, a single \code{data.frame} is returned.
}
\author{Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}}
\seealso{\code{\link{genesel}}, \code{\link{GeneSelection}}, \code{\link{plot,genesel-method}}}
\keyword{multivariate}
