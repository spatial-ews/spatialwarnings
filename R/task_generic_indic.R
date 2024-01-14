# 
# 
# This file contains the workflow function for generic EWS
# 
# 
#' @title Generic Spatial Early-Warning signals
#' 
#' @description Computation, significance assessment and display of spatial 
#'   generic early warning signals (Moran's I, variance and skewness)
#' 
#' @param mat A matrix (quantitative data), a binary matrix (TRUE/FALSE data), 
#'   or a list of those
#' 
#' @param subsize The subsize used for the coarse-graining phase (see Details)
#'   
#' @param abs_skewness Should the absolute skewness be used instead of its 
#'   raw values ? 
#' 
#' @param moranI_coarse_grain Should the input matrix be coarse-grained before
#'   computing the Moran's I indicator value ?
#' 
#' 
#' @return 
#' 
#' \code{generic_sews} returns an object of class \code{simple_sews_single}
#'   (actually a list) if \code{mat} is a single matrix or an object of class 
#'   \code{simple_sews_list} if \code{mat} is a list. You probably want to use some  
#'   of the methods written for these complicated objects instead of extracting 
#'   values directly (they are displayed using \code{print(<object>)}).
#' 
#' \code{indictest} returns an object of class \code{generic_test} (actually 
#'   a data.frame). 
#' 
#' \code{plot} methods return ggplot objects, usually immediately displayed 
#'   when R is being used interactively.
#' 
#' @details 
#' 
#' The Generic Early warning signal are based on the property of a 
#'   dynamical system to "slow down" when approaching a critical point, 
#'   that is take more time to return to equilibrium after a perturbation. This 
#'   is expected to be reflected in several spatial characteristics: the 
#'   variance, the spatial autocorrelation (at lag-1) and the skewness. This 
#'   function provides a convenient workflow to compute these indicators, 
#'   assess their significance and display the results. 
#' 
#' Before computing the actual indicators, the matrix can be "coarse-grained". 
#'   This process reduces the matrix by averaging the nearby cells using 
#'   a square window defined by the \code{subsize} parameter. This makes spatial  
#'   variance and skewness reflect actual spatial patterns when working with 
#'   binary (\code{TRUE}/\code{FALSE} data), but is optional when using 
#'   continuous data. Keep in mind that it effectively reduces the size of 
#'   the matrix by approximately \code{subsize} on each dimension. 
#'   
#' The significance of generic early-warning signals can be estimated by 
#'   reshuffling the original matrix (function \code{indictest}). Indicators 
#'   are then recomputed on the shuffled matrices and the values obtained are 
#'   used as a null distribution. P-values are obtained based on the rank of 
#'   the observed value in the null distribution. A small P-value means 
#'   that the indicator is significantly above the null values, as expected 
#'   before a critical point. 
#'
#' The \code{plot} method can displays the results graphically. A text summary 
#'   can be obtained using the \code{summary} method. 
#' 
#' @references 
#' 
#'   Kefi, S., Guttal, V., Brock, W.A., Carpenter, S.R., Ellison, A.M., 
#'   Livina, V.N., et al. (2014). Early Warning Signals of Ecological 
#'   Transitions: Methods for Spatial Patterns. PLoS ONE, 9, e92097.
#'   
#'   Dakos, V., van Nes, E. H., Donangelo, R., Fort, H., & 
#'   Scheffer, M. (2010). Spatial correlation as leading indicator of 
#'   catastrophic shifts. Theoretical Ecology, 3(3), 163-174.
#'   
#'   Guttal, V., & Jayaprakash, C. (2008). Spatial variance and spatial 
#'   skewness: leading indicators of regime shifts in spatial ecological 
#'   systems. Theoretical Ecology, 2(1), 3-12. 
#'   
#' @seealso 
#'   \code{\link{indictest}}, to test the significance of indicator values. 
#'   Individual indicators: \code{\link{raw_cg_moran}}
#'     \code{\link{raw_cg_variance}}, \code{\link{raw_cg_skewness}}, 
#'     \code{\link{simple_sews}}
#' 
#' @examples
#' 
#' data(serengeti)
#' gen_indic <- generic_sews(serengeti, subsize = 5, 
#'                            moranI_coarse_grain = TRUE)
#' 
#' # Display results
#' summary(gen_indic)
#' 
#' # Display trends along the varying model parameter
#' plot(gen_indic, along = serengeti.rain)
#' 
#' # Compute significance (long)
#' \donttest{
#' gen_test <- indictest(gen_indic, nulln = 199)
#' 
#' print(gen_test)
#' 
#' # Display the trend, now with a grey ribbon indicating the 5%-95% quantile
#' # range of the null distribution
#' plot(gen_test, along = serengeti.rain)
#' 
#' # Display the effect size compared to null distribution 
#' plot(gen_test, along = serengeti.rain, what = "z_score")
#' 
#' # Note that plot() method returns a ggplot object that can be modified
#' # for convenience
#' if ( require(ggplot2) ) { 
#'   plot(gen_test, along = serengeti.rain) + 
#'     geom_vline(xintercept = 733, color = "red", linetype = "dashed") +
#'     xlab('Annual rainfall') + 
#'     theme_minimal()
#' }
#' 
#' }
#' 
#' @export
generic_sews <- function(mat, 
                         subsize = 4, 
                         abs_skewness = FALSE, 
                         moranI_coarse_grain = FALSE) { 
  compute_indicator(mat, raw_generic_indic, 
                    subsize = subsize, 
                    abs_skewness = abs_skewness, 
                    moranI_coarse_grain = moranI_coarse_grain, 
                    taskname = "Generic indicators")
}


# Build the right indicator function (closure) that take into accounts the 
#   above options. 
raw_generic_indic <- function(mat, 
                              subsize, 
                              abs_skewness, 
                              moranI_coarse_grain) { 
  
  # Warn if the matrix is continuous but we will coarse grain anyway
  if ( is.numeric(mat) && subsize > 1 ) { 
    warning(paste("Input matrix has continous values but will be", 
                  "coarse-grained anyway. Set subsize=1 to disable coarse", 
                  "graining."))
  }
  
  if ( is.logical(mat) && subsize == 1 ) { 
    warning(paste("Input matrix is binary but no coarse-graining will be",
                  "performed."))
  }
  
  # We do coarse-graining only once for the whole matrix
  mat_cg <- coarse_grain(mat, subsize)
  
  # Handle the case where there is a single unique value in the matrix
  if ( length(unique(as.vector(mat_cg))) == 1 ) { 
    return( c(variance = var(as.vector(mat_cg)), 
              skewness = NA_real_, 
              moran    = NA_real_, 
              mean     = mean(mat)) )
  }
  
  skewness_value <- cpp_skewness(mat_cg)
  if (abs_skewness) { 
    skewness_value <- abs(skewness_value)
  }
  
  if (moranI_coarse_grain) { 
    moran_value <- raw_moran(mat_cg) 
  } else { 
    moran_value <- raw_moran(mat)
  }
  
  c(variance = var(as.vector(mat_cg)),
    skewness = skewness_value,
    moran    = moran_value,
    mean     = mean(mat))
}


#' @title Spatial variance indicator
#'
#' @description This functions computes the spatial variance of a matrix. 
#' 
#' @param mat A matrix. Its values can be logical, with \code{FALSE} (empty) 
#'   or \code{TRUE} (occupied) values. The entries can also be continuous 
#'   (like NDVI or EVI data). 
#' 
#' @param subsize  Dimension of the submatrix used to coarse-grain the 
#'   original matrix. This must be an integer less than size of the full 
#'   matrix. Coarse-graining reduces the size of the matrix by a factor 
#'   \code{subsize} in each dimension of the matrix. Variance is calculated 
#'   on the coarse-grained matrix. 
#' 
#' @return The variance of the coarse-grained matrix as a named vector 
#' 
#' @details
#' 
#' Spatial variance is a measure of fluctuations in space. Based on the theory 
#' of critical slowing down, when systems approach critical points
#' they are expected to show increased fluctuations in space. Thus, increasing 
#' spatial variance is proposed as an early warning signal of impending 
#' critical transitions. 
#' 
#' Many high resolution spatial data are classified as FALSE (empty) 
#' or TRUE (occupied). In such cases, spatial variance captures just 
#' the variance in data, but not that of spatial structure. 
#' To resolve the issue, this function employs a method called coarse-graining, 
#' proposed in Kefi et al (2014), and described in detail in 
#' Sankaran et al. (2017). One must specify a subsize above one for 
#' binary valued data sets to obtain meaningful values. 
#' 
#' \code{subsize} has to be an integer. It has to be less than or equal to 
#' half of matrix size (N). \code{subsize} must also be preferably a 
#' divisor of N. If it is not a divisor of N, the remainder rows and columns 
#' are discarded when computing coarse-graining matrices. 
#' 
#' Null model evaluations are also done on coarse-grained matrices. 
#' 
#' @references 
#' 
#' Guttal, V., and Jayaprakash, C. (2009). Spatial variance and 
#' spatial skewness: leading indicators of regime shifts in spatial 
#' ecological systems. Theoretical Ecology, 2(1), 3-12.
#' 
#' Kefi, S., Guttal, V., Brock, W.A., Carpenter, S.R., Ellison, A.M., 
#' Livina, V.N., et al. (2014). Early Warning Signals of Ecological 
#' Transitions: Methods for Spatial Patterns. PLoS ONE, 9, e92097.
#' 
#' Sankaran, S., Majumder, S., Kefi, S., and Guttal, V. (2017). Implication 
#' of being discrete and spatial in detecting early warning signals
#' of regime shifts. Ecological Indicators. 
#' 
#' @seealso \code{\link{generic_sews}}
#' 
#' @examples
#' 
#' \donttest{
#' data(serengeti)
#' raw_cg_variance(serengeti[[1]])
#' compute_indicator(serengeti, fun = raw_cg_variance, subsize = 5)
#' }
#' 
#'@export
raw_cg_variance <- function(mat, subsize = 5) { 
  if ( ! is.matrix(mat) ) { 
    stop("raw_cg_variance only accepts a single matrix as input.")
  }
  
  c(variance = var( as.vector( coarse_grain(mat, subsize) ) ))
}


#' @title Moran's Index at lag of 1
#'
#' @description This functions computes the Moran's spatial correlation index 
#'   (with lag one) on a matrix. 
#'
#' @references 
#'
#' Dakos, V., van Nes, E. H., Donangelo, R., Fort, H., & 
#' Scheffer, M. (2010). Spatial correlation as leading indicator of 
#' catastrophic shifts. Theoretical Ecology, 3(3), 163-174.
#'
#' Legendre, P., & Legendre, L. F. J. (2012). Numerical Ecology.
#' Elsevier Science.
#'
#' @param mat A matrix
#' 
#' @param subsize logical. Dimension of the submatrix used to coarse-grain the 
#'   original matrix (set to 1 for no coarse-graining).
#' 
#' @return The Moran's I index measuring autocorrelation at lag 1 as a 
#'   named vector
#'
#' @seealso \code{\link{generic_sews}}
#' 
#' @examples 
#' 
#' \donttest{
#' data(serengeti)
#' raw_cg_moran(serengeti[[1]], subsize = 1)
#' }
#' 
#'@export
raw_cg_moran <- function(mat, subsize = 1) { 
  if ( ! is.matrix(mat) ) { 
    stop("raw_cg_moran only accepts a single matrix as input.")
  }
  c(moran = raw_moran( coarse_grain(mat, subsize) ))
}



#' @title Skewness indicator
#'
#' @description Compute the spatial skewness of spatial data (a matrix).
#'   
#' @param mat A matrix. The matrix 
#'   values can be logical, with \code{FALSE} (empty) or \code{TRUE} (occupied) 
#'   values. The entries can also be continuous (like NDVI or EVI data). 
#' 
#' @param subsize Dimension of the submatrix used to coarse-grain the 
#'   original matrix. This must be an integer less than size of the full 
#'   matrix. Coarse-graining reduces the size of the matrix by a factor 
#'   \code{subsize} in each dimension of the matrix. Skewness is calculated 
#'   on the coarse-grained matrix. 
#' 
#' @param absolute Should the function return the absolute value or raw value 
#'   of skewness ?
#' 
#' @return The spatial skewness of the matrix as a named vector
#' 
#' @details
#' 
#' Spatial skewness is a measure of fluctuations in space; specifically, it 
#' measures if fluctuations are getting biased (skewed) in one direction. Based 
#' on the theory of critical slowing down, when systems approach critical 
#' points they are expected to show increased fluctuations in space. Thus, 
#' increasing spatial skewness is proposed as an early warning signal of 
#' impending critical transitions. 
#' 
#' Computing spatial skewness is straightforward. However, detecting trends of 
#' skewness that correspond to critical slowing down can be tricky, especially 
#' if data come from discrete classification of state variable.
#' 
#' Many high resolution spatial data are classified as FALSE (empty) 
#' or TRUE (occupied by plant). In such cases, spatial skewness captures just 
#' the skewness in data, but not that of spatial structure. 
#' To resolve the issue, this function employs a method called coarse-graining, 
#' proposed in Kefi et al (2014), and described in detail in 
#' Sankaran et al. (2017). One must specify a subsize above one for 
#' binary valued data sets to obtain meaningful values. 
#' 
#' \code{subsize} has to be an integer. It has to be less than or equal to 
#' half of matrix size (N). \code{subsize} must also be preferably a 
#' divisor of N. If it is not a divisor of N, the remainder rows and columns 
#' are discarded when computing coarse-graining matrices. 
#' 
#' Null model evaluations are also done on coarse-grained matrices. 
#' 
#' @references 
#' 
#' Guttal, V., and Jayaprakash, C. (2009). Spatial variance and 
#' spatial skewness: leading indicators of regime shifts in spatial 
#' ecological systems. Theoretical Ecology, 2(1), 3-12.
#' 
#' Kefi, S., Guttal, V., Brock, W.A., Carpenter, S.R., Ellison, A.M., 
#' Livina, V.N., et al. (2014). Early Warning Signals of Ecological 
#' Transitions: Methods for Spatial Patterns. PLoS ONE, 9, e92097.
#' 
#' Sankaran, S., Majumder, S., Kefi, S., and Guttal, V. (2017). Implication of 
#' being discrete and spatial in detecting early warning signals of regime 
#' shifts. Ecological indicators. 
#' 
#' @seealso \code{\link{generic_sews}}
#' 
#' @examples 
#' 
#' data(serengeti)
#' \donttest{
#' raw_cg_skewness(serengeti[[1]])
#' compute_indicator(serengeti, fun = raw_cg_skewness, subsize = 5)
#' }
#' 
#'@export
raw_cg_skewness <- function(mat, subsize = 5, absolute = TRUE) { 
  if ( ! is.matrix(mat) ) { 
    stop("raw_cg_skewness only accepts a single matrix as input.")
  }
  
  a <- cpp_skewness( coarse_grain(mat, subsize) )
  if (absolute) { 
    a <- abs(a) 
  } 
  c(skewness = a)
}
