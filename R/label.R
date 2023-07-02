#' @title Labelling of unique patches and detection of percolation. 
#' 
#' @description Label each patch with a number in a binary matrix 
#' 
#' @param mat A binary matrix
#' 
#' @param nbmask Either "moore" for 8-way neighborhood, "von_neumann" for four-way 
#'   neighborhood, or a square matrix with an odd number of lines and columns that 
#'   describes which neighbors are to be considered around a cell. Default 
#'   is 4-way or Von Neumann neighborhood (the neighborhood of a cell comprises the cell 
#'   above, below, on the right and on the left of the target cell). 
#' 
#' @param wrap Whether to wrap around lattice boundaries (`TRUE`/`FALSE`), 
#'   effectively using periodic boundaries.
#'
#' @return A matrix containing ID numbers for each connected patch. Default 
#'   parameters assume 4-cell neighborhood and periodic boundaries. The 
#'   distribution of patch sizes is returned as the attribute "psd" and the 
#'   percolation status as "percolation" (whether a TRUE patch has a width 
#'   or height equal to the size of the matrix). 
#' 
#' @details The \code{label} function "labels" the patches of a binary 
#'   (TRUE/FALSE) matrix. It returns a matrix of similar height and width, 
#'   with integer values representing the ID of each unique patch (contiguous
#'    cells). Empty cells are labelled as \code{NA}.
#' 
#' @seealso \code{\link{patchsizes}}, \code{\link{patchdistr_sews}}
#' 
#' @examples 
#' 
#' data(forestgap)
#' 
#' rmat <- matrix(rnorm(100) > .1, ncol = 10)
#' display_matrix(label(rmat))
#' 
#' # With 8-way neighborhood mask and no wrapping around borders
#' display_matrix(label(rmat, "moore", wrap = FALSE))
#' 
#' # On real data: 
#' display_matrix(label(forestgap[[5]], nbmask8, wrap = FALSE))
#' 
#' @export
label <- function(mat, 
                  nbmask = "von_neumann", 
                  wrap = FALSE) {
  
  nbmask <- parse_neighbor_spec(nbmask)
  
  if ( ! is.logical(mat) ) { 
    stop('Labelling of patches requirese a logical matrix',
         '(TRUE/FALSE values): please convert your data first.')
  }
  
  if ( ! is.matrix(mat) ) { 
    stop('The input object must be a matrix')
  }
  
  # The matrix is full
  if ( all(mat) ) { 
    result <- matrix(1, nrow = nrow(mat), ncol = ncol(mat)) 
    attr(result, "psd") <- prod(dim(mat))
    attr(result, "percolation") <- TRUE
    return(result)
  } 
  
  # The matrix is empty
  if ( !any(mat) ) { 
    result <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat)) 
    attr(result, "psd") <- integer(0)
    attr(result, "percolation") <- FALSE
    return(result)
  }
  
  # The matrix is a row/column vector 
  if ( ncol(mat) == 1 || nrow(mat) == 1 ) { 
    vec <- as.vector(mat)
    result <- cumsum( c(vec[1] > 0, diff(vec)) == 1 ) * vec
    result <- ifelse(result > 0, result, NA)
    # If we wrap, then we need to merge the two patches at the end of the vector
    if ( wrap && !is.na(head(result, 1)) && !is.na(tail(result, 1)) ) { 
      result[ result == tail(result, 1) ] <- head(result, 1)
    }
    
    # PSD is the just the number of times each unique values appears in the 
    # result vector. 
    attr(result, "psd") <- tabulate(result)
    # Adjust dimensions
    dim(result) <- dim(mat)
    
    # If there is a patch, then it necessarily has the width or height 
    # of the matrix, so percolation is present. 
    attr(result, "percolation") <- any(mat)
    return(result)
  }
  
  # Otherwise we scan for patches
  label_cpp(mat, nbmask, wrap)
}

#' @rdname label
#'
#' @description \code{percolation()} detects whether percolation occurs in the
#'   matrix (i.e. a patch has a width or a height equal to the size of the 
#'   matrix)
#' 
#' @export
percolation <- function(mat, nbmask = "von_neumann") { 
  
  # We never wrap for percolation, by definition. 
  patches <- label(mat, nbmask, wrap = FALSE)
  return(attr(patches, "percolation"))
}


#' @title Get patch sizes.
#' 
#' @description Get the distribution of patch sizes from a logical matrix
#' 
#' @param mat A logical matrix or a list of such matrices.
#' 
#' @param merge Controls whether the obtained patch size distributions are to 
#'   be pooled together if \code{mat} is a list of matrices. 
#' 
#' @param nbmask Either "moore" for 8-way neighborhood, "von_neumann" for four-way 
#'   neighborhood, or a square matrix with an odd number of lines and columns that 
#'   describes which neighbors are to be considered around a cell. Default 
#'   is 4-way or Von Neumann neighborhood (the neighborhood of a cell comprises the cell 
#'   above, below, on the right and on the left of the target cell). 
#' 
#' @param wrap Whether to wrap around lattice boundaries (`TRUE`/`FALSE`), 
#'   effectively using periodic boundaries.
#' 
#' @return If mat is a logical matrix, then the function returns a vector of 
#'   patch sizes. If mat is a list of logical matrices, then it returns 
#'   a list of vectors of patch sizes: this list is flattened if merge is TRUE.
#' 
#' @seealso \code{\link{label}}
#' 
#' @examples
#' 
#' data(forestgap)
#' patchsizes(forestgap[[5]]) # Use a single matrix
#' 
#' # Compute the average patch size of each matrix
#' list_patches <- patchsizes(forestgap) # get the patch size for each matrix
#' print( sapply(list_patches, mean)) # print the average patch size 
#' 
#' # Example with 8-way neighborhood
#' patchsizes(forestgap[[5]], nbmask = "moore")
#' 
#' # Same neighborhood as above, but specified in matrix form 
#' moore_nb <- matrix(c(1, 1, 1, 
#'                      1, 0, 1, 
#'                      1, 1, 1), 
#'                    nrow = 3, ncol = 3, byrow = TRUE)
#' patchsizes(forestgap[[5]], nbmask = moore_nb) 
#' 
#' @export
patchsizes <- function(mat, 
                       merge = FALSE,
                       nbmask = "von_neumann", 
                       wrap = FALSE) { 
  
  if ( is.list(mat)) { 
    result <- lapply(mat, patchsizes, merge, nbmask, wrap) 
    if (merge) { 
      # This always works even if there is only one element
      result <- do.call(c, result)
      names(result) <- NULL
    }
    return(result)
  }
  
  if ( ! is.logical(mat) ) { 
    stop('Computing patch-size distributions requires a logical matrix',
         '(TRUE/FALSE values): please convert your data first.')
  }
  
  # We use the label function -> it returns patch sizes as attributes
  map <- label(mat, nbmask, wrap)
  
  return(attr(map, "psd"))
}


parse_neighbor_spec <- function(spec) { 
  
  moore_nb <- matrix(c(1, 1, 1, 
                       1, 0, 1, 
                       1, 1, 1), 
                     nrow = 3, ncol = 3, byrow = TRUE)
  von_neumann <- matrix(c(0, 1, 0, 
                          1, 0, 1, 
                          0, 1, 0), 
                        nrow = 3, ncol = 3, byrow = TRUE)
  
  if ( is.matrix(spec) ) { 
    
    if( nrow(spec) == 3 && ncol(spec) == 3 ) { 
      return(spec)
    } 
  }
  
  if ( is.numeric(spec) && spec == 4 ) { 
    spec <- "von_neumann"
  }
  
  if ( is.numeric(spec) && spec == 8 ) { 
    spec <- "moore"
  }
  
  if ( is.character(spec) && spec == "moore" ) { 
    return(moore_nb)
  }
  
  if ( is.character(spec) && spec == "von_neumann" ) { 
    return(von_neumann)
  }
  
  stop("Unknown neighborhood specification")
}
