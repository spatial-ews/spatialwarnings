# 
# This contains a function to compute clustering on matrices
# 

#' @title Clustering of pairs
#'
#' @description Compute clustering of pairs in a landscape
#'
#' @param mat A matrix, usually with discrete values (logical, integers, etc.) 
#' 
#' @param wrap Whether space should be considered to wrap at the edges
#' 
#' @param use_8_nb Set to \code{TRUE} to use an 8-way neighborhood. The default is set 
#'   to FALSE, which uses a 4-way neighborhood
#' 
#' @param divide_by_null By default, this function computes the observed number of pairs, 
#'   and divides it by the expected number of pairs under random spatial structure.
#'   Set this parameter to \code{FALSE} to return the raw proportion of pairs. 
#' 
#' @details 
#' 
#'   The clustering of pairs is defined as the density of pairs, i.e. 
#'     the proportion of all pairs of cells that share the same state, divided by the 
#'     null expectation given a random, homogeneous spatial structure. 
#'   
#'   For example, let's consider a matrix with two states, 'a' and 'b'. This function 
#'     will count all pairs of cells 'a-a' or 'b-b', and divide this by the total
#'     number of pairs. This proportion is then again divided by the probability of
#'     obtaining these proportion of pairs under the assumption of no spatial structure
#'     (random mixing of states in the matrix). 
#'     
#'   Clustering is equal to one when there is no spatial structure. It is above one when 
#'     two states tend to be next to each other (i.e. cluster) more than expected by 
#'     chance. Values below one means that those two states tend to be neighbors less 
#'     than expected by chance. 
#'   
#'   If you are only interested in the proportion of pairs of a given state, and you 
#'     do not want to divide by the null expectation, set divide_by_null = FALSE. This 
#'     may happen if you have a more complicated null model than a random, homogeneous
#'     landscape. 
#'   
#' @return 
#' 
#'   A vector with the requested clutering values, whose names are equal to 
#'     each state (unique value) found in the original matrix, preceded by 'clust_' 
#'     (to make sure names are compatible with other functions).
#' 
#' @examples 
#' 
#' # The clustering of a random matrix is close to one
#' ls <- 100 # lattice size
#' mm <- matrix(sample(c("sp1", "sp2", "sp3", "sp4"), size = ls^2, replace = TRUE), 
#'              nrow = ls, ncol = ls)
#' clust <- raw_clustering(mm, wrap = TRUE, use_8_nb = TRUE) 
#' print(clust)
#' 
#' # Compute clustering along the gradient for the serengeti dataset 
#' data(forestgap)
#' clust_indic <- compute_indicator(serengeti, raw_clustering, 
#'                                  wrap = TRUE, use_8_nb = FALSE)
#' # The interesting one is the clustering of state 1 (TRUE in the original matrix), 
#' # which shows similar trends to the variance (as a generic indicator, see 
#' # ?generic_sews)
#' plot(clust_indic, along = serengeti.rain) 
#' # Add null trend
#' clust_test <- indictest(clust_indic, nrep = 19)
#' plot(clust_test, along = serengeti.rain) 
#' 
#'@export
raw_clustering <- function(mat, 
                           wrap = TRUE, 
                           use_8_nb = FALSE,
                           divide_by_null = TRUE) { 
  if ( ! is.matrix(mat) ) { 
    stop("raw_clustering only accepts a single matrix as input.\nUse compute_indicator(<list>, raw_clustering) to make it work on several matrices at once")
  }

  
  if ( is.logical(mat) ) { 
    ilevels <- c("FALSE", "TRUE")
    mat <- matrix(as.integer(mat), nrow = nrow(mat), ncol = ncol(mat))
  }
  
  # If the matrix has floats, convert it to character internally 
  if ( is.numeric(mat) || is.character(mat) ) { 
    dims <- dim(mat)
    mat <- as.factor(mat)
    dim(mat) <- dims
  }
  
  if ( is.factor(mat) ) { 
    ilevels <- levels(mat)
    mat <- matrix(as.integer(mat) - 1, nrow = nrow(mat), ncol = ncol(mat))
  }
  
  counts <- clustering_core(mat, length(ilevels), 
                            wrap = wrap, use_8_nb = use_8_nb)
  
  ncells <- prod(dim(mat))
  nstates <- length(ilevels)
  
  # The number of pairs depends on whether we consider eight neighbors or only four. See 
  # rationale in clustering.cpp on how to count pairs 
  total_pairs <- ifelse(use_8_nb, 4 * ncells, 2 * ncells)
  if ( ! wrap ) { 
    if ( use_8_nb ) { 
      total_pairs <- total_pairs - ( 4 + nrow(mat) + ncol(mat) + 
                                     (nrow(mat)-1)*2 + (ncol(mat)-1)*2 ) 
    } else { 
      total_pairs <- total_pairs - ( 2*nrow(mat) + 2*ncol(mat) ) 
    }
  }
  
  rho_pairs <- counts[ ,1] / total_pairs
  
  rho_null <- ( counts[ ,2] / ncells ) * ( counts[ ,2] - 1 ) / ( ncells - 1 )
  
  if ( divide_by_null ) { 
    clusts <- rho_pairs / rho_null
  }
  
  names(clusts) <- paste0("clust_", ilevels)
  
  return(clusts)
}
