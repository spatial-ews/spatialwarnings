# 
# This file contains functions to test significance and plot simple indicators 
# 

#' @export
indictest.simple_sews <- function(x, 
                                  nulln = 999, 
                                  null_method = 'perm', 
                                  null_control = NULL, 
                                  ...) { 
  NextMethod('indictest')
}
#'@export
indictest.simple_sews_single <- function(x, 
                                         nulln = 999, 
                                         null_method = 'perm', 
                                         null_control = NULL, 
                                         ...) { 
  
  # We do not support low numbers of null simulations
  if ( nulln < 3 ) { 
    stop('The number of null simulations should be above 3 to ', 
         'assess significance')
  }
  
  new_indicf <- function(mat) { 
    do.call(x[['indicf']], c(list(mat), x[['fun.args']]))
  }
  
  # Compute a distribution of null values
  null_values <- compute_indicator_with_null(x[["orig_data"]],
                                             nulln = nulln, 
                                             indicf = new_indicf, 
                                             null_method = null_method, 
                                             null_control = null_control)
  
  # Add data into original object
  x[["nulldistr"]] <- null_values[["nulldistr"]]
  x <- append(x, null_values[["summary_values"]])
  x <- append(x, null_values[["info"]])
  
  class(x) <- c('simple_sews_test_single', 'sews_test', 
                'sews_result_single')
  
  return(x)
}
#'@export
indictest.simple_sews_list <- function(x, 
                                       nulln = 999, 
                                       null_method = 'perm', 
                                       null_control = NULL, 
                                       ...) { 
  
  results <- future_lapply_seed(x, indictest.simple_sews_single, 
                                nulln, null_method, null_control, 
                                ...)
  
  # Add matrixn value with correct number
  for ( nb in seq_along(results) ) { 
    results[[nb]][['matrixn']] <- nb
  }
  
  class(results) <- c('simple_sews_test_list', 'sews_test', 
                      'sews_result_list')
  
  return(results)
}

#'@method as.data.frame simple_sews_test_single
#'@export
as.data.frame.simple_sews_test_single <- function(x, ...) { 
  # Find or create the indicator names
  indicnames <- ifNULLthen(names(x[['value']]), 
                           paste0("indic_", seq_along(x[['value']])))
  
  # Names of components to extract from x 
  test_names <- c("value", "null_mean", "null_sd", "null_qsup",
                  "null_qinf", "z_score", "pval")
  
  # We need to explicitely add a `matrixn` column because it will 
  # be used by funs down the stream. 
  ans <- data.frame(matrixn = 1, 
                    indic = indicnames, 
                    as.data.frame.list(x[test_names]))
  
  return(ans)
}
#'@method as.data.frame simple_sews_test_list
#'@export
as.data.frame.simple_sews_test_list <- function(x, ...) { 
  # Find or create the indicator names
  indicnames <- ifNULLthen(names(x[[1]][['value']]), 
                           paste0("indic_", seq_along(x[[1]][['value']])))

  # Extract the relevant components from the object, and put them in a data 
  # frame. 
  test_names <- c("value", "null_mean", "null_sd", "null_qsup",
                  "null_qinf", "z_score", "pval")
  tab <- Map(function(n, x) { 
    data.frame(matrixn = n, 
               indic = indicnames, 
               as.data.frame.list(x[test_names])) 
  }, seq_along(x), x)
  
  tab <- do.call(rbind, tab)
  
  # Reorder cols so that 'matrixn' is first'
  tab <- data.frame(matrixn = tab[ ,'matrixn'], 
                    tab[ ! names(tab) %in% 'matrixn'])
  tab
}


#'@export
summary.simple_sews_test_single <- function(object, 
                                            taskname = object[["taskname"]], 
                                            ...) { 
  object.list <- list(object)
  summary.simple_sews_test_list(object.list, taskname, ...)
}
#'@export
summary.simple_sews_test_list <- function(object, 
                                          taskname = object[[1]][["taskname"]], 
                                          ...) { 
  
  tab <- as.data.frame(object)
  
  # Format pvalues 
  tab[ ,'stars'] <- pval_stars(tab[ ,'pval'])
  
  tab <- tab[ ,c('matrixn', 'indic', 'value', 'pval', 'stars')]
  
  # Format P-values 
  nulln <- object[[1]][["nulln"]]
  tab[ ,"pval"] <- format.pval(tab[ ,"pval"], digits = DIGITS, eps = 1/nulln)
  
  # Reshape the table 
  tab_pretty <- llply(unique(tab[['indic']]), function(current_indic) { 
    a <- tab[tab[ ,"indic"] == current_indic, c('value', 'pval', 'stars')]
    names(a) <- c(as.character(current_indic), 
                  paste0('pval_', current_indic), 
                  paste0('stars_', current_indic))
    a
  })
  tab_pretty <- do.call(data.frame, tab_pretty)
  tab_pretty <- data.frame(matrixn = seq.int(nrow(tab_pretty)), tab_pretty)
  
  # Build the header to print. Note that we set it to the names of the 
  # data.frame. This is handy because then print.data.frame handles all the 
  # spacing for us, but we effectively create a data.frame with duplicated 
  # column names. I don't know how brittle this is, but probably very. 
  header <- names(tab_pretty) 
  header[grepl('pval_', header)] <- 'P<null'
  header[grepl('stars_', header)] <- ''
  header[grepl('matrixn', header)] <- 'Mat #'
  names(tab_pretty) <- header
  
  cat('Spatial Early-Warning:', taskname, '\n') 
  cat("\n")
  display_size_info(object)

  cat('\n')
  print.data.frame(tab_pretty, row.names = FALSE, digits = DIGITS)
  cat('\n')
  cat(" Null model: ", object[[1]][["null_method"]], " (", 
      object[[1]][["nulln"]], " null matrices used)", "\n", sep = "")
  cat(" Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", '\n')
  cat('\n')
  
  cat("The following methods are available:\n")
  cat(list_methods("simple_sews_test_list"), "\n")
  
  invisible(tab)
}



#'@export
print.simple_sews_test_single <- function(x, ...) { 
  summary.simple_sews_test_single(x, ...)
}
#'@export
print.simple_sews_test_list <- function(x, ...) { 
  summary.simple_sews_test_list(x, ...)
}

#' @title Spatial early-warning signals: display of trends
#' 
#' @param x A \code{simple_sews_test} object (as provided by **_sews functions, such 
#'   as \code{generic_sews()}
#' 
#' @param along A vector providing values over which the indicator trend 
#'   will be plotted. If \code{NULL} then the values are plotted sequentially 
#'   in their original order. 
#' 
#' @details Note that the produced plot is adjusted depending on whether 
#'   \code{along} is numeric or not. 
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
#' 
#'@rdname simple_sews_methods
#'@method plot simple_sews_test
#'@export
plot.simple_sews_test <- function(x, 
                                  along = NULL,
                                  what = "value", 
                                  display_null = TRUE, 
                                  ...) { 
  NextMethod('plot')
}

#'@method plot simple_sews_test_list
#'@export
plot.simple_sews_test_list <- function(x, 
                                       along = NULL, 
                                       what = "value", 
                                       display_null = TRUE, 
                                       ...) { 
  
  # Give the user some slack 
  what <- tolower(what)
  what <- gsub("-", "_", what) # common mistake with the underscore in z-score
  
  # Transform into data.frame
  x.df <- as.data.frame(x)
  
  # Add the difference to the table 
  if ( what == "difference" ) { 
    x.df[ ,"difference"] <- x.df[ ,"value"] - x.df[ ,"null_mean"]
  }
  
  # If along is not provided, then use the matrixn number
  set_default_xlab <- FALSE 
  if ( is.null(along) ) { 
    along <- unique(x.df[ ,"matrixn"])
    set_default_xlab <- TRUE 
  }
  
  check_suitable_for_plots(x.df, along)
  
  plot_data <- data.frame(as.data.frame(x.df), 
                          gradient = along[x.df[ ,'matrixn']])
  
  # Create base ggplot2 object
  # Create base plot object 
  plot <- ggplot2::ggplot(plot_data) + 
            linescale_spwarnings() + 
            theme_spwarnings() + 
            facet_wrap(~ indic, scales = "free_y")
  
  # Check if we really want to add a null ribbon
  add_null <- display_null
  if ( display_null && ! "null_mean" %in% colnames(x.df) ) { 
    warning('Null data was specified to be displayed but could not be found ', 
            'in the provided object')
    add_null <- FALSE
  }
  if ( display_null && what == "value" ) { 
    warning(
      'Display quantiles for null values cannot be shown when the indicator value is ',
      'not displayed'
    )
    add_null <- FALSE
  }
  
  if ( add_null ) { 
    null_data <- data.frame(plot_data,
                            null_ymin = x.df[ ,'null_qinf'],
                            null_ymax = x.df[ ,'null_qsup'])
    
    plot <- plot + 
      ggplot2::geom_ribbon(ggplot2::aes_string(x = 'gradient',
                                               ymin = 'null_ymin',
                                               ymax = 'null_ymax'),
                           data = null_data, 
                           fill = 'grey',
                           group = 1, 
                           alpha = .8) + 
      ggplot2::geom_line(ggplot2::aes_string(x = "gradient", 
                                             y = "null_mean"), 
                         group = 1, 
                         color = 'black', alpha = .1) 
  }
  
  # Add the trend on the graph (Note that we add it over the null ribbon)
  plot <- plot + 
    ggplot2::geom_point(ggplot2::aes_string(x = 'gradient', 
                                            y = what)) + 
    ggplot2::geom_line(ggplot2::aes_string(x = 'gradient', 
                                           y = what))
  
  # Set ylab  
  ylabel <- "Value"
  plot <- plot + ylab(ylabel)
  
  # Set xlab
  if ( set_default_xlab ) { 
    plot <- plot + 
      ggplot2::xlab('Matrix number') 
  } else { 
    plot <- plot + ggplot2::xlab(as.character(substitute(along)))
  }
  
  return(plot)
}
