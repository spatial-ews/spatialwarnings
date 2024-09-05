# 
# 
# This file contains the indictest functions for spectral sews
# 

# 
# Indictest functions for spectral_sews objects.
#'@export
indictest.spectral_sews_list <- function(x, 
                                         nulln = 999, 
                                         null_method = 'perm', 
                                         null_control = NULL, 
                                         ...) { 
  
  # Compute a distribution of null values for SDR
  results <- future_lapply_seed(x, indictest.spectral_sews_single, 
                                nulln, null_method, null_control, ...)
  
  # Format and return output
  class(results) <- c('spectral_sews_test_list', 
                      'spectral_sews_list', 
                      'simple_sews_test_list', 
                      'simple_sews_list', 
                      'sews_result_list')
  
  return(results)
}

#' @export
indictest.spectral_sews_single <- function(x, 
                                           nulln = 999, 
                                           null_method = 'perm', 
                                           null_control = NULL, 
                                           ...) { 
  
  # Build closure passed to compute_indicator_with_null that uses the correct
  #   high and low ranges, and is compatible with the use of matrixn(). 
  sdr_indicf <- function(mat) { 
    spectrum <- rspectrum(mat)
    
    c(sdr = indicator_sdr_do_ratio(spectrum, x[['low_range']], 
                                   x[['high_range']]), 
      spectrum = spectrum[ ,'rspec'])
  }
  
  # Compute null values for SDR
  test_values <- compute_indicator_with_null(x[['orig_data']], 
                                             nulln = nulln, 
                                             indicf = sdr_indicf, 
                                             null_method = null_method, 
                                             null_control = null_control)
  
  # Format results. We import summary stats for SDR results 
  x <- append(x, lapply(test_values[["summary_values"]], function(o) o[1]))
  
  # We extract spectrum values 
  spec <- llply(test_values[["summary_values"]], function(o) o[-1])
  spec <- data.frame(x[["spectrum"]], as.data.frame(spec))
  row.names(spec) <- as.character(seq.int(nrow(spec)))
  x[["spectrum"]] <- spec
  
  # We import information about the null computation in x 
  x <- append(x, test_values[["info"]])
  
  class(x) <- c('spectral_sews_test_single', 
                'spectral_sews_single', 
                'simple_sews_test_single', 
                'simple_sews_single', 
                'sews_result_single')
  
  return(x)
}


# Methods to extract the spectrum from spectral_sews objects
# ---------------------------------------------------------

#' @title Extract the r-spectrum from objects 
#'
#' @description Extract the r-spectrum from objects produced by
#'   \code{spectral_sews}. 
#'
#' @param x An object produced by \code{spectral_sews} or the result of the 
#'   \code{indictest} function called on such object
#' 
#' @param ... Other arguments are ignored 
#' 
#' @return The empirical r-spectrum as a \code{data.frame}
#' 
#' @seealso \code{\link{spectral_sews}}, \code{\link{rspectrum}}
#' 
#' @examples 
#' 
#' # Extract the r-spectrum after computing indicators
#' indics <- spectral_sews(serengeti[2:3])
#' extract_spectrum(indics) 
#' 
#'@export
extract_spectrum <- function(x, ...) { 
  UseMethod("extract_spectrum")
}
#'@export
extract_spectrum.spectral_sews_list <- function(x, ...) { 
  values <- Map(function(n, o) data.frame(matrixn = n, o[["spectrum"]]), 
                seq_along(x), x)
  do.call(rbind, values)
}
#'@export
extract_spectrum.spectral_sews_single <- function(x, ...) { 
  x[["spectrum"]]
}

#' @title Display the r-spectrum of a \code{spectral_sews} object
#' 
#' @description Display the r-spectrum (or multiple spectra) that are contained
#'   in an object returned by \code{\link{spectral_sews}} object (or the result 
#'   of \code{\link{indictest}} applied on such object). 
#' 
#' @param x An object produced by \code{\link{spectral_sews}} or the result 
#'   returned by \code{\link{indictest}} applied on such object
#' 
#' @param along A vector providing numerical or categorical values along 
#'   which the indicator trends will be plotted. If \code{NULL}, then the
#'   indicator values are plotted sequentially in their original order. 
#' 
#' @param log Whether to use a log scale or a linear scale on the y axis
#' 
#' @param display_null Whether to display null information. This argument is 
#'   ignored if \code{x} has not been produced through \code{\link{indictest}} 
#'   (and thus does not contain data regarding the null model)
#' 
#' @param ... Other arguments are ignored 
#' 
#' @seealso \code{\link{rspectrum}}, \code{\link{spectral_sews}}, 
#'   \code{\link{extract_spectrum}}
#' 
#'@export
plot_spectrum <- function(x, 
                          along = NULL, 
                          log = TRUE, 
                          display_null = TRUE, 
                          ...) { 
  UseMethod("plot_spectrum")
}

# Method for indictest output
#'@export
plot_spectrum.spectral_sews_test_list <- function(x, 
                                                  along = NULL, 
                                                  log = TRUE, 
                                                  display_null = TRUE, 
                                                  ...) { 
  
  ggobj <- plot_spectrum.spectral_sews_list(x, along, log)
  
  # Add layers with null model information. We use this method so that the 
  # ribbon appears below the line corresponding to observed values. 
  if ( display_null ) { 
    ggobj$layers <- c(geom_line(aes_string(x = "dist", y = "null_mean"), 
                                color = 'black', alpha = .1), 
                      geom_ribbon(aes_string(x = "dist", ymin = "null_qinf", 
                                            ymax = "null_qsup"), 
                                  fill = 'grey',
                                  group = 1, 
                                  alpha = .8), 
                      ggobj$layers)
  }
    
  return(ggobj) 
}

#'@export
plot_spectrum.spectral_sews_test_single <- function(x, 
                                                  along = NULL, 
                                                  log = TRUE, 
                                                  display_null = TRUE, 
                                                  ...) { 
  
  # Get base plot 
  ggobj <- plot_spectrum.spectral_sews_single(x, along, log)
  
  # Add layers with null model information. We use this method so that the 
  # ribbon appears below the line corresponding to observed values. 
  if ( display_null ) { 
    ggobj$layers <- c(geom_line(aes_string(x = "dist", y = "null_mean"), 
                                color = 'black', alpha = .1), 
                      geom_ribbon(aes_string(x = "dist", ymin = "null_qinf", 
                                            ymax = "null_qsup"), 
                                  fill = 'grey',
                                  group = 1, 
                                  alpha = .8), 
                      ggobj$layers)
  }
  
  return(ggobj)
}

# Method for spectral_sews output (list object)
#'@export
plot_spectrum.spectral_sews_list <- function(x, 
                                             along = NULL, 
                                             log = TRUE, 
                                             display_null = TRUE, 
                                             ...) { 
  
  tab <- extract_spectrum(x)
  check_suitable_for_plots(tab, along)
  
  if ( ! is.null(along) ) { 
    tab[ ,"along"] <- along[tab[ ,"matrixn"]]
  } else { 
    tab[ ,"along"] <- tab[ ,"matrixn"]
  }
  
  p <- ggplot(tab, aes_string(x = "dist", y = "rspec")) + 
    geom_line() + 
    theme_spwarnings() + 
    facet_wrap( ~ along ) + 
    labs(x = "Distance", y = "r-spectrum")
  
  if (log) { 
    p <- p + scale_y_continuous(trans = "log10")
  }
  return(p)
}

# Method for spectral_sews output (single object)
#'@export
plot_spectrum.spectral_sews_single <- function(x, 
                                               along = NULL, 
                                               log = TRUE, 
                                               display_null = TRUE, 
                                               ...) { 
  p <- ggplot(extract_spectrum(x), 
         aes_string(x = "dist", y = "rspec")) + 
    geom_line() + 
    theme_spwarnings() + 
    labs(x = "Distance", y = "r-spectrum")
  
  if (log) { 
    p <- p + scale_y_continuous(trans = "log10")
  }
  
  return(p)
}
