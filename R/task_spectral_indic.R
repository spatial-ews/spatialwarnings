# 
# 
# This file contains a function that compute, assess and display results from 
#   spectrum-based indicators: reddening of powerspectrum, spectrum density
#   ratio.
# 



#' @title Spectrum-based spatial early-warning signals. 
#' 
#' @description Computation of spatial early warning signals based on spectral
#'   properties.
#' 
#' @param mat The input matrix or a list of matrices. 
#' 
#' @param sdr_low_range The range of values (in proportion) to 
#'   use for the computation of the spectral density ratio.
#'   For example, for the lowest 20% (default value), set sdr_low_range to 
#'   c(0, .2). See also the Details section. 
#'  
#' @param sdr_high_range The range of values (in proportion) to 
#'   use for the computation of the spectral density ratio. For example, for 
#'   the higher 20% (default value), set sdr_high_range to 
#'   c(.8, 1). See also the Details section. 
#' 
#' @param quiet Do not display some warnings
#' 
#' @return 
#' 
#' Function \code{spectral_sews} object of class \code{spectral_sews_list} or 
#'   \code{spectral_sews_single} depending on whether the input was a list of 
#'   matrices or a single matrix. 
#' 
#' Function \code{indictest} 
#' 
#' The \code{plot} methods returns a ggplot object (usually displayed 
#' immediately when called interactively). 
#' 
#' @details Spectral early warning signals are based on the fact that some 
#'   dynamical systems can exhibit an change in some characteristics of their 
#'   spatial structure when approaching a transition. In particular, long-range 
#'   correlations are expected to have an increased importance. 
#' 
#' This is expected to be reflected in the spectrum of the spatial structure
#'   by an increase of the relative importance of lower frequencies over higher 
#'   frequencies ("reddening" of the spectrum).
#'   
#' This task allows computing the radial-spectrum which gives the relative 
#'   importance of each space scale as a function of distance, from 1 to 
#'   \code{N/2} (N being the minimum between the number of rows and columns). 
#'   If the matrix is not square, then it is cropped to biggest square that 
#'   fits within the left side of the matrix. 
#' 
#' Additionally, it summarizes this spectrum into a Spectral 
#'   Density Ratio (SDR), which is the ratio of low frequencies over 
#'   high frequencies of the r-spectrum. The SDR value is expected to increase
#'   before a transition point.
#' 
#' The significance of spectral early-warning signals can be estimated by 
#'   reshuffling the original matrix (function \code{indictest}). Indicators 
#'   are then recomputed on the shuffled matrices and the values obtained are 
#'   used as a null distribution. P-values are obtained based on the rank of 
#'   the observed value in the null distribution. 
#' 
#' The trend of SDR values can be plotted using the \code{plot()} method. 
#'   Alternatively, the spectrum itself can be plotted (with facets 
#'   if multiple input matrices were used) using the \code{plot_spectrum} 
#'   method.
#' 
#' 
#' @references 
#' 
#'   Kefi, S., Guttal, V., Brock, W.A., Carpenter, S.R., Ellison, A.M., 
#'   Livina, V.N., et al. (2014). Early Warning Signals of Ecological 
#'   Transitions: Methods for Spatial Patterns. PLoS ONE, 9, e92097.
#' 
#' @seealso \code{\link{rspectrum}}, \code{\link{plot_spectrum}}, 
#'   \code{\link{raw_sdr}}, \code{\link{extract_spectrum}}
#' 
#' @seealso 
#'   \code{\link{indictest}}, to test the significance of indicator values. 
#' 
#' @examples
#' \dontrun{
#' 
#' data(serengeti) 
#' data(serengeti.rain) 
#' 
#' 
#' spec_indic <- spectral_sews(serengeti, 
#'                             sdr_low_range  = c(0, .2), 
#'                             sdr_high_range = c(.8, 1))
#' 
#' summary(spec_indic)
#' 
#' # Display trends along the varying model parameter
#' plot(spec_indic, along = serengeti.rain)
#' 
#' # Computing spectra many times is expensive, consider setting parallel 
#' # computing using: options(mc.cores = n)
#' 
#' # Assess significance
#' spec_test <- indictest(spec_indic, nulln = 199)
#' summary(spec_test)
#' 
#' # Display the SDR trend, now with a grey ribbon representing 5%-95% 
#' # quantiles of the null distribution
#' plot(spec_test, along = serengeti.rain)
#' 
#' # Add a line highlighting the shift 
#' if (require(ggplot2)) {
#'   plot(spec_test, along = serengeti.rain) + 
#'     geom_vline(xintercept = 590, color = "red", linetype = "dashed")
#' }
#' 
#' 
#' # Display radial-spectra
#' plot_spectrum(spec_indic, along = serengeti.rain)
#' 
#' # Graphics can be modified using ggplot2 functions
#' if (require(ggplot2)) { 
#'   plot_spectrum(spec_indic, along = serengeti.rain) + 
#'     scale_y_log10()
#' }
#' 
#' }
#' @export
spectral_sews <- function(mat, 
                          sdr_low_range  = NULL, 
                          sdr_high_range = NULL, 
                          quiet = FALSE) { 
  
  if ( is.null(sdr_low_range) ) { 
    if ( !quiet ) { 
      warning("Choosing the 20% lowest frequencies for spectral density ratio ",
              "as no range was specified. Use parameter sdr_low_range to ", 
              "choose a different value.")
    }
    sdr_low_range <- c(0, .2)
  }
  
  if ( is.null(sdr_high_range) ) { 
    if ( !quiet ) { 
      warning("Choosing the 20% highest frequencies for spectral density ", 
              "ratio as no range was specified. Use parameter sdr_high_range ", 
              "to choose a different value.")
    }
    sdr_high_range <- c(.8, 1)
  }
  
  # Handle list case
  if ( is.list(mat) ) { 
    results <- lapply(mat, spectral_sews, sdr_low_range, 
                      sdr_high_range, quiet)
    names(results) <- names(mat)
    class(results) <- c('spectral_sews_list',  
                        'spectral_sews', 
                        'simple_sews_list', 
                        'simple_sews', 
                        'sews_result_list')
    
    return(results)
  }
  
  # Convert object to matrix and check if it is suitable for spatialwarnings
  mat <- convert_to_matrix(mat)
  check_mat(mat)
  
  # Now the mat is always a matrix -> process it
  # Check and warn if not square
  if ( nrow(mat) != ncol(mat) ) { 
    if ( !quiet ) { 
      warning('The matrix is not square: only a squared subset of the left ", 
              "part will be used.')
    }
    
    mindim <- min(dim(mat))
    mat <- mat[seq.int(mindim), seq.int(mindim)]
  }
  
  # Compute powerspectrum
  spectrum <- rspectrum(mat)
  
  # Compute SDR
  ranges_absolute <- convert_ranges_to_absolute(mat, sdr_low_range, 
                                                sdr_high_range)
  
  sdr_value <- indicator_sdr_do_ratio(spectrum, 
                                      ranges_absolute[["low"]], 
                                      ranges_absolute[["high"]])
  
  # Return list containing both
  output <- list(value = c(sdr = sdr_value), 
                 spectrum = spectrum, 
                 orig_data = mat, 
                 low_range = ranges_absolute[['low']], 
                 high_range = ranges_absolute[['high']], 
                 taskname = "Spectrum-based indicators")
  class(output) <- c('spectral_sews_single', 
                     'spectral_sews', 
                     'simple_sews_single', 
                     'simple_sews', 
                     'sews_result_single')
  
  return(output)
}

# 
#' @title Spectral Density Ratio (SDR) indicator
#' 
#' @description Compute the ratio of low frequencies over high frequencies
#'   of the r-spectrum. 
#' 
#' @param mat A matrix with continuous values, or a logical matrix 
#'   (TRUE/FALSE). 
#' 
#' @param sdr_low_range The range of values (in proportion) to 
#'   use for the computation of the spectral density ratio.
#'   For example, for the lowest 20\% (default value), set \code{sdr_low_range}
#'   to \code{c(0, .2)}.
#'  
#' @param sdr_high_range The range of values (in proportion) to 
#'   use for the computation of the spectral density ratio. For example, for 
#'   the highest 20\% (default value), set \code{sdr_high_range} to 
#'   \code{c(.8, 1)}. 
#' 
#' @return The SDR values computed on the matrix as a named vector
#' 
#' @details 
#' 
#' SDR measures the increase in long-range correlations before a critical point. 
#'   It is the ratio of the average low frequency value over high frequency 
#'   values. In this implementation, an increase in SDR implies a "reddening" 
#'   of the \link[=rspectrum]{r-spectrum}. See also \code{\link{spectral_sews}} for 
#'   a more complete description. 
#' 
#' Low and high frequencies are averaged in order to compute the SDR. The 
#'   parameters \code{sdr_low_range} and \code{sdr_high_range} control which 
#'   frequencies are selected for averaging. For example 
#'   \code{sdr_low_range = c(0, .2)} (default) uses the lower 20% to compute 
#'   the average of low frequencies. \code{sdr_high_range = c(.8, 1)} uses the 
#'   higher 20% for the average of high frequencies. 
#' 
#' @seealso \code{\link{indictest}}, 
#'   \code{\link{rspectrum}}, \code{\link{plot_spectrum}}, 
#'   \code{\link{spectral_sews}}, \code{\link{extract_spectrum}}
#' 
#' @references 
#' 
#' Carpenter, S.R. & Brock, W.A. (2010). Early warnings of regime shifts in 
#'   spatial dynamics using the discrete Fourier transform. Ecosphere
#' 
#' @examples 
#' 
#' \dontrun{ 
#' data(serengeti)
#' serengeti.sdr <- raw_sdr(serengeti[[1]], 
#'                          sdr_low_range = c(0, 0.2), 
#'                          sdr_high_range = c(0.8, 1))
#' compute_indicator(serengeti, raw_sdr)
#' }
#' 
#' @export
raw_sdr <- function(mat, 
                    sdr_low_range  = NULL, 
                    sdr_high_range = NULL) { 
  
  check_mat(mat) # checks if sensible
  
  if ( is.null(sdr_low_range) ) { 
    warning("Choosing the 20% lowest frequencies for spectral density ratio ",
            "as none was specified. Use parameter sdr_low_range to choose ", 
            "a different value.")
    sdr_low_range <- c(0, .2)
  }
  
  if ( is.null(sdr_high_range) ) { 
    warning("Choosing the 20% highest frequencies for spectral density ratio ",
            "as none was specified. Use parameter sdr_high_range to choose ", 
            "a different value.")
    sdr_high_range <- c(.8, 1)
  }
  
  ranges_absolute <- convert_ranges_to_absolute(mat, sdr_low_range, 
                                                sdr_high_range)
  
  c(sdr = indicator_sdr_core(mat, 
                             ranges_absolute[["low"]],
                             ranges_absolute[["high"]]))
}

indicator_sdr_core <- function(mat, low_range, high_range) { 
  
  # Compute r-spectrum
  spectrum <- rspectrum(mat)
  
  # Compute ratio
  return( indicator_sdr_do_ratio(spectrum, low_range, high_range) )
}

indicator_sdr_do_ratio <- function(spectrum, low_range, high_range) { 
  
  # Compute subsets
  low_subset  <- with(spectrum, dist <= max(low_range)  & 
                                  dist >= min(low_range))
  high_subset <- with(spectrum, dist <= max(high_range) & 
                                  dist >= min(high_range))
  
  # If the number of values to estimate means is very low, then warn. 
  if ( sum(low_subset) < 3 || sum(high_subset) < 3 ) { 
    warning('The number of values used to compute the SDR ratio is very low ', 
            'it may be unreliable')
  }
  
  # Return ratio of means
  with(spectrum, mean(rspec[low_subset]) / mean(rspec[high_subset])) 
}

# Convert ranges from proportional values to absolute distance values
convert_ranges_to_absolute <- function(mat, 
                                       sdr_low_range  = NULL, 
                                       sdr_high_range = NULL) { 
    
  maxdist <- 1 + floor(min(dim(mat)) / 2)
  low_range_absolute  <- sdr_low_range * maxdist
  high_range_absolute <- sdr_high_range * maxdist
  
  return( list(low = low_range_absolute, 
               high = high_range_absolute) )
}
