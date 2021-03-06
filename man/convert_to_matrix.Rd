% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_matrix.R
\name{convert_to_matrix}
\alias{convert_to_matrix}
\title{Convert an object to a matrix}
\usage{
convert_to_matrix(object, ...)
}
\arguments{
\item{object}{An object (typically, a matrix or a list of matrices)}

\item{...}{Additional arguments (currently ignored)}
}
\description{
This function is mainly for internal use by the 
  \code{\link{spatialwarnings}} package to convert objects before they are 
  processed by \code{*_sews} functions.
}
\details{
This generic function is here so that other packages can extend it. 
For example, \href{https://github.com/spatial-ews/spatialwarningsGis}{spatialwarningsGis} will provide methods so that GIS objects can be handled 
(e.g. \code{RasterLayer} from package \code{raster}).
}
\examples{

# this does nothing
convert_to_matrix(serengeti[2:3]) 

}
