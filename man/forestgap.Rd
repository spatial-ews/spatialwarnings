% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{forestgap}
\alias{forestgap}
\alias{forestgap.pars}
\title{A list of binary matrices and their associated parameters}
\format{
A list of logical matrices which are the end results of simulations 
  from Kubo's Forest Gap model along a gradient of increasing values of 
  stress (see references).

The parameters used for the simulations, as a data frame.
}
\source{
Generated using the implementation of Kubo's model in caspr 0.2.0 
  \url{https://github.com/fdschneider/caspr}.
}
\usage{
forestgap

forestgap.pars
}
\description{
A list of binary matrices and their associated parameters
}
\details{
Kubo's forest gap model has three parameters, \eqn{\alpha}{alpha} 
  that controls the reproductive rate of trees, \eqn{d}{d} controls the 
  non-spatialized mortality and \eqn{\delta}{delta} the increased mortality 
  due to the presence of a neighboring gap.
}
\references{
Kubo, T., Iwasa, Y., & Furumoto, N. (1996). Forest spatial dynamics with gap
  expansion: Total gap area and gap size distribution. Journal of Theoretical
  Biology, 180(3), 229-246. \doi{10.1006/jtbi.1996.0099}
}
\keyword{datasets}
