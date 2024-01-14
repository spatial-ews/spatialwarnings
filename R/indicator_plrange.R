#
# This file contains the indicator function based on the power-law range
#   (how much of the data behaves like a power law?)
#
#' @title Power-law range indicator
#' 
#' @description Compute the power-law range of a matrix 
#' 
#' @param mat A logical matrix, or a list of logical matrices 
#'
#' @param merge Controls whether the patch-size distributions of the input 
#'   matrices are merged together before computing the power-law range. Setting 
#'   this value to \code{TRUE} makes the function return a single value even 
#'   if multiple matrices are given as input. 
#'
#' @param xmin_bounds A vector of two integer values, defining a range in which 
#'   to search for the best xmin (see Details). 
#' 
#' @details 
#' 
#' Some ecosystems show typical changes in their patch-size 
#' distribution as they become more and more degraded. In particular, an 
#' increase in the truncation of the patch-size distribution (PSD) is expected 
#' to occur. The power-law range (PLR) measures the truncation of the PSD 
#' in a single value (see also \code{\link{patchdistr_sews}} for more details). 
#' 
#' To compute the PLR, power-laws are fitted with a variable 
#' minimum patch size (xmin) and the one with the lowest Kolmogorov-Smirnov
#' distance to the empirical distribution is retained. PLR is then computed 
#' using this best-fitting xmin: 
#' 
#' \deqn{\frac{log(x_{max}) - log(x_{min})}{log(x_{max}) - log(x_{smallest})}}{ (log(xmax) - log(xmin))/(log(xmax) - log(xsmallest))}
#' 
#' where \eqn{x_{max}}{x_max} is the maximum observed patch size, and 
#' \eqn{x_{smallest}}{x_smallest} is the minimum observed patch size. 
#' 
#' @references 
#' 
#' Clauset, A., Shalizi, C. R., & Newman, M. E. (2009). 
#'   Power-law distributions in empirical data. SIAM review, 51(4), 661-703.
#' 
#' Berdugo, M., Kefi, S., Soliveres, S. & Maestre, F.T. (2017). Plant spatial 
#' patterns identify alternative ecosystem multifunctionality states in 
#' global drylands. Nature in Ecology and Evolution.
#' 
#' @return A data.frame with columns minsize, maxsize which are the observed 
#'   minimum and maximum patch sizes. The estimated \eqn{x_{min}}{x_min} and the 
#'   value of the power-law range. If multiple matrices were provided, then 
#'   a list of data.frames is returned
#' 
#' @seealso \code{\link{patchdistr_sews}}
#' 
#' @examples
#' \donttest{
#' forestgap.plr <- indicator_plrange(forestgap) 
#' do.call(rbind, forestgap.plr) # convert results to data.frame
#' 
#' # Restrict to small xmins 
#' forestgap.plr2 <- indicator_plrange(forestgap, xmin_bounds = c(1, 10)) 
#' do.call(rbind, forestgap.plr2) 
#' }
#'@export
indicator_plrange <- function(mat,
                              merge = FALSE,
                              xmin_bounds = NULL) {
  
  if ( !merge && is.list(mat) ) {
      # Returns a list of lists
    return( lapply(mat, indicator_plrange, merge, xmin_bounds) )
  }
  
  psd <- patchsizes(mat, merge = merge)
  
  if ( is.null(xmin_bounds) ) {
    xmin_bounds <- range(psd)
  }
  
  # If there are not enough patches to work with -> return NA
  if ( length(unique(psd)) <= 2 ) { 
    warning('Not enough different patch sizes to estimate xmin: returning NA')
    result <- data.frame(min(psd), max(psd), NA, NA)
  } else {
    # Compute xmin and range
    plrange_result <- plrange(psd, xmin_bounds) # returns xmin also
    result <- data.frame(min(psd), max(psd), plrange_result)
  }
  
  names(result) <- c("minsize", "maxsize", "xmin_est", "plrange")
  return(result)
}
#
# This file contains the indicator function based on the power-law range
#   (how much of the data behaves like a power law?)
#
#' @title Power-law range indicator
#' 
#' @description Compute the power-law range of a matrix 
#' 
#' @param mat A logical matrix, or a list of logical matrices 
#'
#' @param xmin_bounds A vector of two integer values, defining a range in which 
#'   to search for the best xmin (see Details). 
#' 
#' 
#' @details 
#' 
#' Some ecosystems show typical changes in their patch-size 
#' distribution as they become more and more degraded. In particular, an 
#' increase in the truncation of the patch-size distribution (PSD) is expected 
#' to occur. The power-law range (PLR) measures the truncation of the PSD 
#' in a single value (see also \code{\link{patchdistr_sews}} for more details). 
#' 
#' To compute the PLR, power-laws are fitted with a variable 
#' minimum patch size (xmin) and the one with the lowest Kolmogorov-Smirnov
#' distance to the empirical distribution is retained. PLR is then computed 
#' using this best-fitting xmin: 
#' 
#' \deqn{\frac{log(x_{max}) - log(x_{min})}{log(x_{max}) - log(x_{smallest})}}{ (log(xmax) - log(xmin))/(log(xmax) - log(xsmallest))}
#' 
#' where \eqn{x_{max}}{x_max} is the maximum observed patch size, and 
#' \eqn{x_{smallest}}{x_smallest} is the minimum observed patch size. 
#' 
#' @references 
#' 
#' Clauset, A., Shalizi, C. R., & Newman, M. E. (2009). 
#'   Power-law distributions in empirical data. SIAM review, 51(4), 661-703.
#' 
#' Berdugo, M., Kefi, S., Soliveres, S. & Maestre, F.T. (2017). Plant spatial 
#' patterns identify alternative ecosystem multifunctionality states in 
#' global drylands. Nature in Ecology and Evolution.
#' 
#' @return A named vector containing the power-law range value
#' 
#' @seealso \code{\link{patchdistr_sews}}
#' 
#' @examples
#' \donttest{
#' forestgap.plr <- raw_plrange(forestgap[[2]]) 
#' 
#' # Restrict to small xmins 
#' forestgap.plr2 <- indicator_plrange(forestgap[[2]], xmin_bounds = c(1, 10)) 
#' }
#'@export
raw_plrange <- function(mat, xmin_bounds = NULL) { 
  
  psd <- patchsizes(mat)
  
  if ( is.null(xmin_bounds) ) {
    xmin_bounds <- range(psd)
  }
  plrange_result <- plrange(psd, xmin_bounds)
  
  return( c(plrange = plrange_result[ ,"plrange"]) )
}

plrange <- function(psd, xmin_bounds) {
  
  if ( any(psd == 0) ) { 
    warning(paste0("You have zeros in your patch size distribution, plrange ", 
                   "will be NaN"))
  }
  
  # If psd is empty, then return NA
  if ( length(unique(psd)) <= 1) {
    return( data.frame(xmin_est = NA_real_, plrange = NA_real_) )
  }
  
  xsmallest <- min(psd)
  xmax <- max(psd)
  xmin <- xmin_estim(psd, bounds = xmin_bounds)
  
  if ( is.na(xmin) ) { # finding xmin failed
    result <- data.frame(NA_real_, NA_real_)
  } else {
    result <- data.frame(xmin, 1 - (log10(xmin) - log10(xsmallest)) /
                                (log10(xmax) - log10(xsmallest)))
  }
  
  names(result) <- c('xmin_est', 'plrange')
  
  return(result)
}
