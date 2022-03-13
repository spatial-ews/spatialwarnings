# 
# This file defines a generic function to convert matrices. spatialwarningsGis
# will add methods to this generic function so that foreign gis formats can be 
# directly used in spatialwarnings. 
# 
#'
#'
#' @title Convert an object to a matrix
#' 
#' @description This function is mainly for internal use by the 
#'   \code{\link{spatialwarnings}} package to convert objects before they are 
#'   processed by \code{*_sews} functions. 
#' 
#' @param object An object (typically, a matrix or a list of matrices)
#' 
#' @param ... Additional arguments (currently ignored)
#' 
#' @details 
#'   This generic function is here so that other packages can extend it. 
#' For example, \href{https://github.com/spatial-ews/spatialwarningsGis}{spatialwarningsGis} will provide methods so that GIS objects can be handled 
#' (e.g. \code{RasterLayer} from package \code{raster}).
#' 
#' @examples 
#' 
#' # this does nothing
#' convert_to_matrix(serengeti[2:3]) 
#' 
#'@export
convert_to_matrix <- function(object, ...) { 
  UseMethod("convert_to_matrix")
}

#'@export
convert_to_matrix.matrix <- function(object, ...) { 
  object
}

#'@export 
convert_to_matrix.factor <- function(object, ...) { 
  if ( is.matrix(object) ) { 
    convert_to_matrix.matrix(object) 
  } else { 
    stop("This factor object has no dimensions. I cannot convert it to a matrix.")
  }
}

#'@export
convert_to_matrix.list <- function(object, ...) { 
  lapply(object, convert_to_matrix) 
}

