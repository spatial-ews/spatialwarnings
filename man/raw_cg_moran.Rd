% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/task_generic_indic.R
\name{raw_cg_moran}
\alias{raw_cg_moran}
\title{Moran's Index at lag of 1}
\usage{
raw_cg_moran(mat, subsize = 1)
}
\arguments{
\item{mat}{A matrix}

\item{subsize}{logical. Dimension of the submatrix used to coarse-grain the 
original matrix (set to 1 for no coarse-graining).}
}
\value{
The Moran's I index measuring autocorrelation at lag 1 as a 
  named vector
}
\description{
This functions computes the Moran's spatial correlation index 
  (with lag one) on a matrix.
}
\examples{

\donttest{
data(serengeti)
raw_cg_moran(serengeti[[1]], subsize = 1)
}

}
\references{
Dakos, V., van Nes, E. H., Donangelo, R., Fort, H., & 
Scheffer, M. (2010). Spatial correlation as leading indicator of 
catastrophic shifts. Theoretical Ecology, 3(3), 163-174.

Legendre, P., & Legendre, L. F. J. (2012). Numerical Ecology.
Elsevier Science.
}
\seealso{
\code{\link{generic_sews}}
}
