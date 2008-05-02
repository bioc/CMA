\name{khan}
\alias{khan}
\docType{data}
\title{Small blue round cell tumor dataset of Khan et al. (2001)}
\description{s. below}
\usage{data(khan)}
\format{
  A data frame with 63 observations on the following 2309 variables.
  The first column (named \code{khanY}) contains the tumor classes
  (BL = Burkitt Lymphoma, EWS = Ewing Sarcoma, NB = Neuro Blastoma,
   RMS = Rhabdomyosarcoma).\\
  %\describe{
    \code{khanY}: a factor with levels \code{BL} \code{EWS} \code{NB} \code{RMS} \\
    \code{X2-X2309}: Gene expression values.
  }

\source{Adopted from the dataset in the package \code{pamr}.}
\references{Khan, J., Wei, J. S., Ringner, M., Saal, L. H., Ladanyi, M., Westermann, F., Berthold,
F., Schwab, M., Antonescu, C. R., Peterson, C., Meltzer, P. S., (2001). \cr
Classification and diagnostic prediction of cancers using gene expression profiling and artificial neural
networks. \cr
\emph{Nature Medicine 7, 673-679}.}


\examples{
data(khan)}
\keyword{datasets}
