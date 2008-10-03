\name{golub}
\alias{golub}
\docType{data}
\title{ALL/AML dataset of Golub et al. (1999)}
\description{s. below}
\usage{data(golub)}
\format{
  A data frame with 38 observations and 3052 variables. The first
  column (named \code{golub.cl}) contains the tumor classes
  (ALL = acute lymphatic leukaemia, AML = acute myeloid leukaemia).\\
    \code{golub.cl}: a factor with levels \code{ALL} \code{AML}.\\
    \code{X2-X3051}: Gene expression values.
  %}
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
  Adopted from the dataset in the package \code{multtest}.
}
\references{Golub, T., Slonim, D. K., Tamayo, P., Huard, C., Gaasenbeek, M., Mesirov, J. P.,

           Coller, H., Loh, M. L., Downing, J., Caligiuri, M. A., Bloomfeld, C. D., Lander, E. S. (1999).

           Molecular classification of cancer: class discovery and class prediction by gene expression monitoring.

           \emph{Science 286, 531-537}.
}
\examples{
data(golub)
}
\keyword{datasets}
