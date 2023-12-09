# 
# 

# 
# 
# 
# Testing method 
#'@export
indictest.variogram_sews_list <- function(x, 
                                          nulln = 999, 
                                          null_method = 'perm', 
                                          null_control = NULL, 
                                          ...) { 
  
  results <- future_lapply_seed(x, indictest.variogram_sews_single, 
                                nulln, null_method, null_control, ...)
  
  # Add matrixn column with correct matrixn number
  for ( nb in seq_along(results) ) { 
    results[[nb]][['matrixn']] <- nb
  }
  
  class(results) <- c('variogram_sews_test_list', 
                      'variogram_sews_list', 
                      'variogram_sews_test', 
                      'simple_sews_test_list', 
                      'sews_result_list')
  return(results)
}
#'@export
indictest.variogram_sews_single <- function(x, 
                                            nulln = 999, 
                                            null_method = 'perm', 
                                            null_control = NULL, 
                                            ...) { 
  
  # This function will produce a vector, with the 4 first values holding 
  # the summarized metrics, and the last ones the values of the variogram. 
  metric_compute <- function(mat) { 
    a <- with(x[["pars"]], fit_variogram(mat, model, nmax, nbins, cutoff))
    unlist(with(a, c(compute_vario_metrics(pars), vario[ ,"gamma"])))
  }
  
  test_values <- compute_indicator_with_null(x[['orig_data']], 
                                             nulln = nulln, 
                                             indicf = metric_compute, 
                                             null_method = null_method, 
                                             null_control = null_control)
  
  # Format results. The first four values are parameters, the rest is the 
  # variogram.
  x <- append(x, llply(test_values[["summary_values"]], function(o) o[1:4]))
  
  # Import variogram data in original object
  vario <- lapply(test_values[["summary_values"]], function(o) o[-(1:4)])
  vario <- data.frame(x[["variogram"]], as.data.frame(vario))
  row.names(vario) <- as.character(seq.int(nrow(vario)))
  x[["variogram"]] <- vario
  
  # Import information
  x <- append(x, test_values[["info"]])
  
  class(x) <- c('variogram_sews_test_single', 
                'variogram_sews_single', 
                'simple_sews_test_single', 
                'sews_result_single')
  return(x)
}

#' @rdname variogram_sews_plot
#' 
#' @param what The trendline to be displayed. Defaults to the indicator's 
#'   values ("value") but other metrics can be displayed. Accepted values are 
#'   "value", "pval", "difference" (obs - null mean), or "z_score" 
#'   ( (obs - null mean) / (null sd) ).
#' 
#' @param display_null Chooses whether a grey ribbon should be added to reflect
#'   the null distribution. Note that it can not be displayed when the trend 
#'   line reflects something else than the indicator values (when \code{what} 
#'   is not set to "value").
#' 
#' @param ... Other arguments are ignored.
#'@export
plot.variogram_sews_test <- function(x, 
                                along = NULL, 
                                what = "value", 
                                display_null = TRUE, 
                                ...) { 
  NextMethod("plot")
}
#'@export
plot.variogram_sews_test_list <- function(x, along = NULL, 
                                          what = "value", 
                                          display_null = TRUE, 
                                          ...) { 
  plot.simple_sews_test_list(x, along, what, display_null, ...)
}

# We need to put that into the plot help page for variogram-based indics
#'@rdname variogram_sews_plot
#'@method plot_variogram variogram_sews_test
#'@export
plot_variogram.variogram_sews_test <- function(x, 
                                along = NULL, 
                                what = "value", 
                                display_null = TRUE, 
                                ...) { 
  NextMethod("plot")
}
#'@export
plot_variogram.variogram_sews_test_list <- function(x, along = NULL, ...) { 
  ggobj <- plot_variogram.variogram_sews_list(x, along = along)
  
  # Extract null values and display them 
  variodf <- extract_variogram(x)
  if ( ! is.null(along) ) { 
    variodf[ ,"along"] <- along[variodf[ ,"matrixn"]]
  } else { 
    variodf[ ,"along"] <- variodf[ ,"matrixn"]
  }
  
  ggobj$layers <- c(geom_ribbon(aes_string(x = "dist", ymin = "null_qinf", 
                                           ymax = "null_qsup"), 
                                data = variodf, 
                                fill = 'grey',
                                group = 1, 
                                alpha = .8), 
                    geom_line(aes_string(x = "dist", y = "null_mean"), 
                              data = variodf, 
                              color = 'black', alpha = .1), 
                    ggobj$layers)
  return(ggobj)
}

#'@export
plot_variogram.variogram_sews_test_single <- function(x, ...) { 
  ggobj <- plot_variogram.variogram_sews_single(x)
  # Extract null values and display them 
  variodf <- extract_variogram(x)
  
  ggobj$layers <- c(geom_ribbon(aes_string(x = "dist", ymin = "null_qinf", 
                                           ymax = "null_qsup"), 
                                data = variodf, 
                                fill = 'grey',
                                group = 1, 
                                alpha = .8), 
                    geom_line(aes_string(x = "dist", y = "null_mean"), 
                              data = variodf, 
                              color = 'black', alpha = .1), 
                    ggobj$layers)
  return(ggobj)
}
