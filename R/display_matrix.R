# 
# This file contains a function to vizualise matrices as those we use in
#   spatialwarnings 
# 

#' @title Plot a matrix 
#' 
#' @description Display a matrix or a list of matrices in a plot 
#' 
#' @param object A matrix, a list of matrices, an object produced by 
#'   \code{*_sews} functions or \code{indictest()}
#' 
#' @param palette A color palette to use in the plot. It can be any color 
#'   palette understood by \link[ggplot2]{scale_fill_brewer}. 
#' 
#' @param along A vector of values used in facet headers. If \code{object} is 
#'   a matrix, this is ignored. 
#' 
#' @param ... Other arguments are ignored. 
#' 
#' @details This function will plot a matrix using ggplot2, using the provided 
#'   palette. Its use-case is very much like \code{image()}, but its produces 
#'   nicer plots by default (\code{image()} is much faster, however). 
#' 
#' @return A ggplot2 object, which is printed when this function is used 
#'   interactively. 
#' 
#' @examples 
#' 
#' # display_matrix works with single matrices or lists of matrices 
#' display_matrix(serengeti[2:3])
#' 
#' # display_matrix is compatible with "*_sews" objects
#' indics <- compute_indicator(serengeti[2:3], raw_moran) 
#' display_matrix(indics)
#' 
#'@export
display_matrix <- function(object, palette = "RdYlBu", along = NULL, ...) { 
  UseMethod("display_matrix")
}

#'@export
display_matrix.matrix <- function(object, palette = "RdYlBu", 
                                  along = NULL, ...) { 
  check_mat(object, NAerror = FALSE) 
  
  if ( is.numeric(object) ) { 
    fillscale <- scale_fill_distiller(palette = palette)
  } else { 
    fillscale <- scale_fill_brewer(palette = palette)
  }
  
  ggplot(tabularize(object), aes_string(x = "col", y = "row")) + 
    geom_raster(aes_string(fill = "value")) + 
    fillscale + 
    scale_y_reverse() + 
    theme_spwarnings() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank()) + 
    coord_fixed()
}

#'@export 
display_matrix.list <- function(object, palette = "RdYlBu", 
                                along = NULL, ...) { 
                                  
  # Convert and check objects
  object <- lapply(object, convert_to_matrix)
  lapply(object, check_mat, NAerror = FALSE)

  # Convert all matrices to data frames
  all_tabs <- Map(function(n, o) { data.frame(matrixn = n, tabularize(o)) }, 
                  seq_along(object), object)
  all_tabs <- do.call(rbind, all_tabs)
  
  # Add the along vector 
  check_suitable_for_plots(all_tabs, along)
  if ( ! is.null(along) ) { 
    all_tabs[ ,"along"] <- along[all_tabs[ ,"matrixn"]]
  } else { 
    all_tabs[ ,"along"] <- all_tabs[ ,"matrixn"]
  }
  
  # Choose proper color scale 
  if ( is.numeric(all_tabs[ ,"value"]) ) { 
    fillscale <- scale_fill_distiller(palette = palette)
  } else { 
    fillscale <- scale_fill_brewer(palette = palette)
  }
  
  ggplot(all_tabs, aes_string(x = "col", y = "row")) + 
    geom_raster(aes_string(fill = "value")) + 
    facet_wrap( ~ along ) + 
    fillscale + 
    scale_y_reverse() + 
    theme_spwarnings() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank()) + 
    coord_fixed()
}

#'@export 
display_matrix.sews_result_list <- function(object, palette = "RdYlBu", 
                                            along = NULL, ...) { 
  
  # Extract original data from object 
  orig_dat <- lapply(object, function(o) o[["orig_data"]])
  display_matrix.list(orig_dat, palette = palette, along = along, ...)
}

#'@export 
display_matrix.sews_result_single <- function(object, palette = "RdYlBu", 
                                              along = NULL, ...) { 
  
  # Extract original data from object 
  display_matrix.matrix(object[["orig_data"]], palette = palette, 
                        along = NULL, ...)
}

#'@export
display_matrix.default <- function(object, palette = "RdYlBu", 
                                   along = NULL, ...) { 
  display_matrix( convert_to_matrix(object) )
}


tabularize <- function(mat) { 
  data.frame(expand.grid(row = seq.int(nrow(mat)), 
                         col = seq.int(ncol(mat))), 
             value = as.vector(mat))
}

