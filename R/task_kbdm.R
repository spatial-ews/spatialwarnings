
#'
#' @title Indicator based on Kolmogorov Complexity 
#' 
#' @description 
#' 
#'   Computes the Kolmogorov Complexity on a set of matrices, 
#'     using the Block Decomposition Method. 
#'
#' @details 
#' 
#    The Kolmogorov complexity of a given matrix has been suggested to 
#'     be a useful indicator to anticipate transitions in model ecological 
#'     systems (Dakos and Soler-Toscano, 2017). When close to the transition 
#'     critical point, the complexity is expected to decrease. 
#'   
#'   The Kolmogorov complexity cannot be computed directly for large strings 
#'     (i.e. matrices). However, the complexity of smaller submatrices can be 
#'     estimated, then combined to obtain an approximation of the complexity 
#'     of the whole matrix. This method, the Block Decomposition Method is 
#'     implemented in this indicator following Dakos and Soler-Toscano (2017). 
#' 
#' @return 
#' 
#'   \code{kbdm_sews} returns an object of class \code{simple_sews_single} 
#'     (a list) if mat is a single matrix, and an object of class 
#'     \code{simple_sews_list} if mat is a list of matrices. These objects can 
#'     be used with generic methods indictest (to test significance) or plot 
#'     (to display trends), see also the examples below. 
#'   
#' @references 
#'   
#'   Dakos, V., and F. Soler-Toscano. 2017. Measuring complexity to infer 
#'   changes in the dynamics of ecological systems under stress. Ecological 
#'   Complexity 32:144-155. 
#'   
#' @param mat A logical matrix (\code{TRUE}/\code{FALSE} values) or a list of logical 
#'   matrices
#' 
#' @param subsize A submatrix size to carry out the Block Decomposition Method
#'   (must be between 1 and 3)
#' 
#' @seealso \code{\link{raw_kbdm}}, \code{\link[acss]{acss}}, 
#'   \code{\link{indictest}}, to test the significance of indicator values. 
#' 
#' @examples 
#' 
#' \donttest{ 
#' 
#' kbdm_result <- kbdm_sews(serengeti, subsize = 3)
#' plot(kbdm_result, along = serengeti.rain)
#' 
#' kbdm_test <- indictest(kbdm_result, nulln = 49)
#' plot(kbdm_test, along = serengeti.rain)
#' 
#' # Plot deviation to null expectation
#' plot(kbdm_test, along = serengeti.rain, what = "z_score") 
#' 
#' }
#' 
#' 
#'@export 
kbdm_sews <- function(mat, subsize = 3) { 
    compute_indicator(mat, raw_kbdm, subsize = subsize, 
                      taskname = "Kbdm Complexity")
}



# This function takes a matrix and a returns a single value.
#' 
#' @title Kolmogorov complexity of a matrix 
#' 
#' @description Compute the Kolmogorov complexity of a matrix using the 
#'   Block Decomposition Method (requires the \code{acss} package).
#' 
#' @details 
#' 
#'   The Kolmogorov complexity cannot be computed directly for large strings 
#'     (i.e. matrices). However, the complexity of smaller submatrices can be 
#'     estimated, then combined to obtain an approximation of the complexity 
#'     of the whole matrix. This method, the Block Decomposition Method is 
#'     implemented in this function. See also \code{\link{kbdm_sews}} 
#'     for more details. 
#' 
#' @return The KBDM numeric value as a named vector
#' 
#' @param mat A logical matrix (with TRUE/FALSE values)
#' 
#' @param subsize A submatrix size to carry out the Block Decomposition Method
#'   (must be between 1 and 3)
#' 
#' @seealso \code{\link{kbdm_sews}}, \code{\link[acss]{acss}}
#' 
#' @examples 
#' 
#' \donttest{ 
#' raw_kbdm(forestgap[[1]], subsize = 3)
#' }
#' 
#'@export
raw_kbdm <- function(mat, subsize) {
  
  if ( ! requireNamespace("acss") ) { 
    stop(paste0('Computation of kbdm requires the package acss. Install it ', 
                'using install.packages("acss")'))
  }
  
  if ( subsize > 3 || subsize < 1 ) { 
    stop("subsize must be between 1 and 3.")
  }
  
  # Split matrix
  xs <- seq(1, nrow(mat), by = subsize)
  ys <- seq(1, ncol(mat), by = subsize)
  allblockns <- expand.grid(seq.int(length(xs)-1),
                            seq.int(length(ys)-1))
  all_substr <- Map(function(xblockn, yblockn) {
                      dict <- as.vector(mat[xs[xblockn]:(xs[xblockn+1]-1),
                                            ys[yblockn]:(ys[yblockn+1]-1)]) 
                      dict <- as.integer(dict) 
                      dict <- paste(dict, collapse = "")
                    },
                    allblockns[ ,1], allblockns[ ,2]) 
  all_substr <- unlist(all_substr)
  
  # Summarize the substrings
  counts <- table(all_substr)
  counts <- data.frame(string = names(counts),
                       multip = as.vector(counts),
                       kctm = acss::acss(names(counts), alphabet = 2)[ ,1])
  
  # Compute Kbdm
  return( c(kbdm = with(counts, sum(log2(multip) + kctm))) )
}
  
