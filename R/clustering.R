#
# This contains a function to compute clustering on matrices
#

#' @title Clustering of pairs
#'
#' @description Compute the number of pairs of neighbor cells in a landscape and derive 
#'   from it clustering indices
#'
#' @param mat A matrix, usually with discrete values (logical, integers, etc.)
#'
#' @param wrap Whether space should be considered to wrap at the edges
#'
#' @param use_8_nb Set to \code{TRUE} to use an 8-way neighborhood. The default is set
#'   to FALSE, which uses a 4-way neighborhood
#'
#' @param prop Return the counts for all pairs in the matrix (FALSE), or
#'   the proportions for each pair (out of all possible pairs)
#' 
#' @details
#'
#'   The clustering of pairs is defined as the density of pairs, i.e.
#'     the proportion of all neighboring pairs of cells that share the same state, 
#'     divided by the null expectation given a random, homogeneous spatial structure.
#'
#'   For example, let's consider a matrix with two states, 'a' and 'b'.
#'     \code{raw_clustering} will count all pairs of cells 'a-a' or 'b-b', and
#'     divide this by the total number of pairs. This proportion is then again
#'     divided by the probability of obtaining these proportion of pairs under
#'     the assumption of no spatial structure (random mixing of cells in the
#'     matrix). 
#'
#'   Clustering is equal to one when there is no spatial structure. It is
#'     above one when two states are found next to each other (i.e. cluster)
#'     more than expected by chance. Values below one means that those two
#'     states tend to be neighbors less frequently than expected by chance. 
#'
#'   If you are only interested in the proportion of pairs for each combination
#'     of states, you can use the function \code{pair_counts}, which returns
#'     a matrix with as many rows and columns as there are states in the matrix, and 
#'     contains the counts for all possible pairs of cells found in the matrix 
#'     (or their relative proportions).
#'
#' @return
#'
#'   A vector with the requested clutering values for \code{raw_clustering},
#'     whose names are equal to each state (unique value) found in the
#'     original matrix, preceded by 'clust_' (to make sure names are
#'     compatible with other functions).
#'
#'  \code{pair_counts} returns the counts of pairs of states in the matrix,
#'     or the proportion of each pair, depending on the value of \code{prop}
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
#' \donttest{
#' 
#' data(forestgap)
#' clust_indic <- compute_indicator(serengeti, raw_clustering,
#'                                  wrap = TRUE, use_8_nb = FALSE)
#' # The interesting one is the clustering of state 0 (FALSE in the original matrix),
#' # which corresponds to grassland pixels, which get more and more clustered with 
#' # increasing rainfall (see also ?generic_sews for how that compares with generic 
#' # indicators)
#' plot(clust_indic, along = serengeti.rain)
#' 
#' # Add null trend
#' clust_test <- indictest(clust_indic, nulln = 19)
#' plot(clust_test, along = serengeti.rain)
#' 
#' # Show the proportion of each pairs of states in the matrix... 
#' pair_counts(serengeti[[5]])
#' 
#' # ... or the total count
#' pair_counts(serengeti[[5]], prop = FALSE)
#' 
#' }
#' 
#'@export
raw_clustering <- function(mat,
                           wrap = TRUE,
                           use_8_nb = FALSE) {

  mat_fmt <- format_mat_for_pair_count(mat)
  mat <- mat_fmt[[1]]
  ilevels <- mat_fmt[[2]]

  ncells <- prod(dim(mat))

  # Count the number of neighboring pairs of each state
  total_pairs <- pair_counts(mat, wrap = wrap, use_8_nb = use_8_nb,
                             prop = TRUE)
  pair_props <- diag(total_pairs)

  state_props <- sapply(seq_along(ilevels)-1, function(lvl) {
    mean(mat == lvl)
  })

  # N1 = number of cells in state 1
  # n = total number of cells
  # P1 = proportion of cells in state 1 = N1/n
  # P(two time N1) = ( N1 / n * (N1-1) / (n-1) ) =
  #   = P1 * ( P1 * n - 1 ) / ( n - 1)
  pair_nulls <- state_props * ( state_props * ncells - 1 ) / ( ncells - 1 )
  clust_pairs <- pair_props / pair_nulls

  names(clust_pairs) <- paste0("clust_", ilevels)

  return(clust_pairs)
}


#'@export
#'@rdname raw_clustering
pair_counts <- function(mat,
                        wrap = TRUE,
                        use_8_nb = FALSE,
                        prop = TRUE) {
  mat_fmt <- format_mat_for_pair_count(mat)
  mat <- mat_fmt[[1]]
  ilevels <- mat_fmt[[2]]

  pairs <- pair_counts_internal(mat, length(ilevels), wrap, use_8_nb)
  rownames(pairs) <- colnames(pairs) <- ilevels

  if ( prop ) {
    total_pairs <- total_pair_count(mat, wrap, use_8_nb)
    pairs <- pairs / total_pairs
  }

  # Make a semi-triangular matrix
  pairs[lower.tri(pairs)] <- pairs[upper.tri(pairs)] + pairs[lower.tri(pairs)]
  pairs[upper.tri(pairs)] <- NA

  return(pairs)
}


# Format a matrix before feeding it to pair-counting cpp functions
format_mat_for_pair_count <- function(mat) {

   if ( ! is.matrix(mat) ) {
     msg <- paste0("clustering functions only accept a single matrix as input. ",
                   "Use compute_indicator(<list>, fun) to use several matrices ",
                   "at once")
     stop(msg)
  }

  if ( is.logical(mat) ) {
    dims <- dim(mat)
    mat <- factor(c("FALSE", "TRUE")[mat+1], levels = c("TRUE", "FALSE"))
    dim(mat) <- dims
  }

  # If the matrix has floats, convert it to factor. We don't know the full
  # range of levels in mat, so just use as.factor()
  if ( is.numeric(mat) || is.character(mat) ) {
    dims <- dim(mat)
    mat <- as.factor(mat)
    dim(mat) <- dims
  }

  # mat is always a factor at this point
  lvls <- levels(mat)

  if ( is.factor(mat) ) {
    mat <- matrix(as.integer(mat) - 1, nrow = nrow(mat), ncol = ncol(mat))
  }

  return( list(mat, lvls) )
}

# Get the number of pairs for a given matrix size and options. There is
# more information about pair counting in clustering.cpp
total_pair_count <- function(mat, wrap, use_8_nb) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  ncells <- nr * nc
  total_pairs <- ifelse(use_8_nb,
                        4 * ncells,
                        2 * ncells)
  if ( ! wrap ) {
    if ( use_8_nb ) {
      total_pairs <- total_pairs - ( 2 + nr + nc + (nr-1)*2 + (nc-1)*2 )
    } else {
      total_pairs <- total_pairs - (nr + nc) # ok
    }
  }

  return(total_pairs)
}
