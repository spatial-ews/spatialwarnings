# A function that computes the "simple" flowlength as described in Rodriguez 
#   et al. 
# 
#' @title Flowlength connectivity indicator (uniform topography)
#' 
#' @description Measures the connectivity of runoff-source areas as determined 
#'   by vegetation patterns and (uniform) topography
#' 
#' @details 
#'   
#' This function computes Flowlength, a simple metric that measures the 
#'   potential hydrological connectivity of runoff-source areas 
#'   (e.g., bare soil) considering vegetation cover, vegetation patterns and 
#'   topography. Flowlength is defined as the average length of all the 
#'   potential runoff pathways in the target area. Thus, a higher value of 
#'   the index indicates a higher hydrologic connectivity of runoff source 
#'   areas. This function is designed for an idealized uniform hillslope 
#'   (e.g., with constant slope angle, the direction of maximum slope being 
#'   from the top to the bottom of the input matrices). 
#' 
#' The deviations of Flowlength from its expected values under random or 
#'   aggregated-pattern null models can be used as an indicator of imminent 
#'   transition to a degraded state (Rodriguez et al. 2017) in the context 
#'   of arid drylands. An increased deviation of flowlength compared to its 
#'   null values is expected as a possible transition gets closer. This deviation
#'   can be computed using any null model provided by \code{spatialwarnings} (see 
#'   \code{\link{indictest}} for more details), but a specific null model is 
#'   provided for Flowlength based on a much-faster analytical approximation, 
#'   using the argument \code{null_method = "approx_rand"} when calling 
#'   \code{indictest} (see examples below). 
#' 
#' In general, Flowlength can be used as indicator of dryland functional status 
#'   by assessing potential water and soil losses in patchy landscapes 
#'   (Mayor et al. 2008, Moreno-de las Heras et al. 2012, Mayor et al. 2013, 
#'   Okin et al. 2015). Finally, the combination of observed and expected 
#'   Flowlength under null models for random or aggregated vegetation cover 
#'   can be used for assessing the cover-independent role of bare-
#'   soil connectivity (Rodriguez et al. 2018).
#'   
#' @references
#' 
#' Rodriguez, F., A. G. Mayor, M. Rietkerk, and S. Bautista. 2017. A null model 
#'   for assessing the cover-independent role of bare soil connectivity as 
#'   indicator of dryland functioning and dynamics. Ecological Indicators.
#' 
#' Mayor, A.G., Bautista, S., Small, E.E., Dixon, M., Bellot, J., 2008. 
#'   Measurement of the connectivity of runoff source areas as determined by 
#'   vegetation pattern and topography: a tool for assessing potential water 
#'   and soil losses in drylands. Water Resour. Res. 44, W10423.
#' 
#' Mayor, A.G., Kefi, S., Bautista, S., Rodriguez, F., Carteni, F., Rietkerk, 
#'   M., 2013. Feedbacks between vegetation pattern and resource loss 
#'   dramatically decrease ecosystem resilience and restoration potential in 
#'   a simple dryland model. Landsc. Ecol. 28, 931-942.
#' 
#' Moreno-de las Heras, M., Saco, P.M., Willgoose, G.R., Tongway, D.J., 2012. 
#'   Variations in hydrological connectivity of Australian semiarid landscapes 
#'   indicate abrupt changes in rainfall-use efficiency of vegetation. 
#'   J. Geophys. Res. 117, G03009.
#' 
#' Okin, G.S., Moreno-de las Heras, M., Saco, P.M., Throop, H.L., Vivoni, E.R., 
#'   Parsons, A.J., Wainwright, J., Peters, D.P.C., 2015. Connectivity in 
#'   dryland landscapes: shifting concepts of spatial interactions. 
#'   Front. Ecol. Environ. 13 (1), 20-27.
#' 
#' @param mat The input matrix (must be a logical matrix)
#' 
#' @param slope The slope of the area documented by the matrix (in degree). 
#' 
#' @param cell_size The horizontal size of a cell in the matrix (as viewed 
#'   from above). 
#' 
#' @return A `simple_sews` object containing the flow length value, among 
#'   other things, see \code{\link{simple_sews_object}} for more information.
#' 
#' @seealso 
#'   \code{\link{raw_flowlength_uniform}}, 
#'   \code{\link{indictest}} to test the significance of indicator values. 
#' 
#' @examples 
#' 
#' \donttest{ 
#' fl_result <- flowlength_sews(arizona, slope = 20, cell_size = 1)
#' 
#' # Compute the Z-score (standardized deviation to null distribution) and plot 
#' #   its variations along the gradient. This Z-score is suggested by 
#' #   Rodriguez et al. (2017) as an indicator of degradation. 
#' fl_test <- indictest(fl_result, nulln = 19)
#' plot(fl_test, what = "z_score")
#' 
#' # Use the analytical approximation suggested in Rodriguez et al. (2017), 
#' # instead of permuting the original values in the matrix (much faster)
#' fl_test <- indictest(fl_result, null_method = "approx_rand")
#' plot(fl_test, what = "z_score")
#' 
#' }
#' 
#'@export
flowlength_sews <- function(mat,        # Input matrix
                            slope = 20, # Slope (in degrees)
                            cell_size = 1) { # Cell size
  
  fl_result <- compute_indicator(mat, raw_flowlength_uniform, 
                                 slope = slope, 
                                 cell_size = cell_size, 
                                 taskname = "Flow length (uniform topography)")
  
  # Make sure we adjust the class to add flowlength-related information. This is 
  # used later on so that indictest() will use specific methods written for the 
  # flowlength. 
  if ( inherits(fl_result, "sews_result_list") ) { 
    class(fl_result) <- c("flowlength_sews_list", class(fl_result))
    for ( i in seq_along(fl_result) ) { 
      class(fl_result[[i]]) <- c("flowlength_sews_single", class(fl_result[[i]]))
    }
  } else { 
    class(fl_result) <- c("flowlength_sews_single", class(fl_result))
  }
  
  return(fl_result)
}

#'@export
indictest.flowlength_sews_list <- function(x, 
                                           nulln = 999, 
                                           null_method = "perm", 
                                           null_control = NULL, 
                                           ...) { 
  
  results <- future_lapply_seed(x, indictest.flowlength_sews_single, 
                                nulln, null_method, null_control, 
                                ...)
  
  # Add matrixn value with correct number
  for ( nb in seq_along(results) ) { 
    results[[nb]][['matrixn']] <- nb
  }
  
  class(results) <- c('flowlength_sews_test_list', 
                      'simple_sews_test_list', 'sews_test', 
                      'sews_result_list')
  
  return(results)
}

#'@export
indictest.flowlength_sews_single <- function(x, 
                                             nulln = 999, 
                                             null_method = "perm", 
                                             null_control = NULL, 
                                             ...) { 
  
  # Read null_control values we may use to compute null expectation for FL
  qinf <- ifelse(is.null(null_control[["qinf"]]), .05, null_control[["qinf"]]) 
  qsup <- ifelse(is.null(null_control[["qsup"]]), .95, null_control[["qsup"]])
  
  # These methods use an analytical approximation to compute the expected values
  # of flowlength given an homogeneous matrix, or (TODO) one that is somewhat
  # aggregated (see Rodriguez et al., Ecological Indicators).
  if ( ( ! is.function(null_method) ) && ( null_method == "approx_rand" ) ) { 
    
    rho <- mean(x[["orig_data"]])
    L <- nrow(x[["orig_data"]])
    slope <- x[["fun.args"]][["slope"]]
    cell_size <- x[["fun.args"]][["cell_size"]]
    ds <- cell_size / cos(slope * pi / 180)
    
    # Eq. 1 in the above paper (with ds = 1)
    E_fl <- ds * (1 - rho)*(rho * L - (1 - rho)*(1 - (1 - rho)^L)) / (rho^2*L)
    E_fl2 <- ds^2*( (1 - rho)*(rho^2*(1 - rho)*(L+1)^2 + rho^2*L - 6 * (1 - rho) + 
                        (1 - rho)^(L+1)*(rho^2*(2*L^2 - 1) + 6 * rho * L + 6)) ) / 
                        (rho^4 * L^2)
    V_fl <- E_fl2 - E_fl^2
    
    # NOTE: the values above are for a column vector, we need to compute 
    # the variance for the average of FL for all columns, which is different
    V_fl <- V_fl / ( ncol(x[["orig_data"]]) )
    
    # Now add the values in the resulting test object 
    # "nulldistr"   "null_mean"   "null_sd"     "null_qsup"   "null_qinf"  
    # "z_score"     "pval"        "null_method" "nulln"       "get_nullmat"
    ans <- x
    
    # These are undefined because we never computed a null distribution 
    # ans[["nulldistr"]] <- NULL
    # ans[["nulln"]] <- NULL
    # ans[["get_nullmat"]] <- NULL
    
    ans[["null_mean"]] <- E_fl
    ans[["null_sd"]] <- sqrt(V_fl) 
    ans[["null_method"]] <- null_method
    ans[["null_qsup"]] <- qnorm(qsup, mean = E_fl, sd = sqrt(V_fl))
    ans[["null_qinf"]] <- qnorm(qinf, mean = E_fl, sd = sqrt(V_fl))
    # P-value close to zero when indicator is above the null distribution
    ans[["pval"]] <- pnorm(x[["value"]], mean = E_fl, sd = sqrt(V_fl), 
                           lower.tail = FALSE)
    ans[["z_score"]] <- ( x[["value"]] - E_fl ) / sqrt(V_fl)
    
    class(ans) <- c("flowlength_sews_test_single", 
                    'simple_sews_test_single', 'sews_test', 
                    'sews_result_single')
    
    return(ans)
  }
  
  # Catch-all: we revert to the default, usual way of doing it 
  ans <- indictest.simple_sews_single(x, nulln, null_method, null_control)
  class(ans) <- c("flowlength_sews_test_single", 
                  'simple_sews_test_single', 'sews_test', 
                  'sews_result_single')
  
  return(ans)
}

#
# A function that computes the "simple" flowlength as described in Rodriguez 
#   et al. 
# 
#' @title Flow length (uniform slope)
#' 
#' @description Compute a simple approximation of the flow length assuming a 
#'   constant slope 
#' 
#' @details 
#'   
#'   This function computes the Flowlength of a given matrix, using a 
#'     uniform approximation (the slope is constant across the whole matrix, 
#'     with maximum slope being from the top of the matrix to its bottom), 
#'     as per Rodriguez et al. (2017). See \code{\link{flowlength_sews}} for 
#'     more details. 
#'   
#' @references
#' 
#' Rodriguez, F., A. G. Mayor, M. Rietkerk, and S. Bautista. 2017. A null model 
#'   for assessing the cover-independent role of bare soil connectivity as 
#'   indicator of dryland functioning and dynamics. Ecological Indicators.
#' 
#' @param mat The input matrix (must be a logical matrix)
#' 
#' @param slope The slope of the area documented by the matrix (in degrees). 
#' 
#' @param cell_size The horizontal size of a cell in the matrix (as viewed 
#'   from above). 
#' 
#' @return A named vector of length 1 containing the flow length numerical 
#'   value 
#' 
#' @seealso \code{\link{flowlength_sews}}
#' @seealso 
#'   \code{\link{indictest}}, to test the significance of indicator values. 
#' 
#' @examples 
#' 
#' \donttest{ 
#' raw_flowlength_uniform(arizona[[1]], slope = 20, cell_size = 1)
#' }
#' 
#' 
#'@export
raw_flowlength_uniform <- function(mat,        # Input matrix
                                   slope, # Slope (in degrees)
                                   cell_size) { # Cell size
  
  if ( is.vector(mat) ) { 
    mat <- matrix(mat, ncol = 1, nrow = length(mat))
  }
  
  # Mat must be a logical matrix 
  if ( ! is.matrix(mat) ) { 
    stop('The supplied object is not a matrix')
  }
  
  # Mat must be a logical matrix 
  if ( ! is.logical(mat) ) { 
    stop('The supplied matrix must be logical to compute the flow length')
  }
  
  # Take the negative of mat, so that ones represent empty cells
  nmat <- ! mat 
  
  # Cell size along the slope 
  p_slope <- cell_size / cos(slope * pi / 180)
  
  # Compute flow length without taking slope into account (C++ function)
  flmean <- fl_internal(nmat)
  
  # Adjust value for slope
  fl <- flmean * p_slope
  
  return(c(flowlength = fl))
}
