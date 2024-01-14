# 
# 
# This file contains methods related to objects of the patchdistr task
# 



# Plot methods
# --------------------------------------------------
# 
#' @rdname patchdistr_sews_plot
#' @name patchdistr_sews_plot
#' 
#' @title Early-warning signals based on patch size distributions
#' 
#' @description Plot early-warning signals based on patch size distributions 
#' 
#' @param x An object as produced by \code{\link{patchdistr_sews}}, or 
#'   the result of \code{\link{indictest}} called on such object
#' 
#' @param along A vector providing values along which the indicator trends 
#'   will be plotted. If \code{NULL} then the values are plotted sequentially 
#'   in their original order. 
#' 
#' @param ... Ignored 
#' 
#' @details 
#'   
#'   The \code{plot} function will produce a figure summarizing the changes
#'   in patch size distributions along a set of values. The figure has two 
#'   panels: 
#'   \itemize{ 
#'      \item the upper panel shows the percolation status of empty 
#'        (\code{FALSE}) and occupied cells (\code{TRUE}), and shows the mean 
#'        value (proportion of \code{TRUE} values). The background shows 
#'        the proportion of each type of distribution for each unique values 
#'        of the \code{along} vector. 
#'        
#'      \item the bottom panel displays the power-law range
#'   }
#'  
#'  The \code{plot_distr} function displays each distribution in an 
#'    individual facet, with an overlay of the best distribution fit and a blue 
#'    bar showing the power-law range. If appropriate, a grey ribbon is used to
#'    display the expected distribution given the null expectation (i.e. when 
#'    \code{plot_distr} is called on the results of \code{indictest()}. This 
#'    function can produce quite crowded graphs, but it displays in full the
#'    shape of the distributions, and can be useful e.g. to assess the quality 
#'    of the fits. 
#' 
#' @seealso \code{\link{patchdistr_sews}}
#' 
#' @examples
#' 
#' \donttest{ 
#' data(forestgap)
#' psd_indic <- patchdistr_sews(forestgap)
#' 
#' plot(psd_indic, along = forestgap.pars[ ,"d"]) 
#' 
#' # When along is non-numeric, bars are used for display
#' plot(psd_indic, along = as.factor(forestgap.pars[ ,"d"]))
#' 
#' # Display individual distributions
#' plot_distr(psd_indic, along = forestgap.pars[ ,"d"])
#' 
#' # We can display the distributions along with the null expectation after 
#' # indictest() is run
#' psd_test <- indictest(psd_indic, nulln = 19)
#' plot_distr(psd_test, along = forestgap.pars[ ,"d"])
#' }
#'
#'@method plot patchdistr_sews
#'@export
plot.patchdistr_sews <- function(x, along = NULL, ...) { 
  if ( 'patchdistr_sews_single' %in% class(x) ) { 
    stop('I cannot plot a trend with only one value')
  }
  
  plot.patchdistr_sews_list(x, along)
}

# Note: plot will display how characteristics change along a gradient. To 
#   have a plot of distributions, use plot_distr
# 
#'@method plot patchdistr_sews_list
#'@export
plot.patchdistr_sews_list <- function(x, along = NULL, ...) { 
  
  if ( !is.null(along) && (length(along) != length(x)) ) { 
    stop('The along values are unfit for plotting (size mismatch)')
  }
  
  obj_table <- as.data.frame(x)
  
  # Subset table for best fits
  obj_table <- obj_table[is.na(obj_table[ ,'best']) | obj_table[ ,'best'], ]
  
  if ( all(is.na(obj_table[ ,"best"])) && 
       all(is.na(obj_table[ ,"plrange"])) ) { 
    warning("There is no non-missing data to plot in the provided object")
    return( invisible(NULL) )
  }
  
  # If along is provided, then add it to the table
  xtitle <- deparse(substitute(along))
  if ( is.null(along) ) { 
    along <- as.factor(obj_table[ ,'matrixn'])
    xtitle <- "Matrix number"
  }
  
  obj_table[ ,"along"] <- along
  
  # Now we summarise the obj_table
  alltypes <- na.omit(unique(obj_table[ ,"type"]))
  summ <- ddply(obj_table, 'along',
                function(df) { 
                  type_freqs <- sapply(alltypes, 
                         function(current_type) { 
                           sum(!is.na(df[ ,'type']) & df[ 'type'] == current_type)
                         })
                  type_freqs <- type_freqs / ifelse(sum(type_freqs) > 0, 
                                                    sum(type_freqs), 1)

                  data.frame(type = alltypes, type_freq = type_freqs, 
                             percolation = mean(df[ ,'percolation']), 
                             percolation_empty = mean(df[ ,'percolation_empty']), 
                             cover = mean(df[ ,"cover"]), 
                             # We remove NAs as sometime the power-law range cannot 
                             # be computed (not enough points)
                             plrange = mean(df[ ,'plrange'], na.rm = TRUE))
                  })
  
  # Make a data.frame with dummy facets 
  classif_df <- data.frame(plot_type = "PSD Classif.", summ)
  classif_df[ ,'plrange'] <- NA
  
  plrange_df <- ddply(summ, "along", function(df) { 
                       data.frame(plot_type = "PL-range", 
                                  plrange = df[1, "plrange"])
                      })
  
  summ <- rbind.fill(classif_df, plrange_df)
  
  # Construct base plot
  plot <- ggplot(summ) + 
    theme_spwarnings() + 
    fillscale_spwarnings(name = "PSD type") + 
    scale_linetype_discrete(name = "Cover") + 
    facet_grid( plot_type ~ ., switch = "y") + 
    ylab('') + 
    xlab(xtitle)
  
  if ( ! is.numeric(along) ) { 
    # Note that it is quite tricky to make ggplot produce a line over factors: 
    # `group=1` seems to do the trick. 
    plot <- plot + 
      geom_bar(aes_q(x = ~along, y = ~type_freq, fill = ~type), 
               position = "stack", stat = "identity") + 
      stat_summary(aes_q(x = ~along, y = ~percolation, 
                         linetype = "Perc. (full)", group = 1), 
                   fun = mean, geom = "line") + 
      stat_summary(aes_q(x = ~along, y = ~percolation_empty, 
                         linetype = "Perc. (empty)", group = 1), 
                   fun = mean, geom = "line") + 
      stat_summary(aes_q(x = ~along, y = ~cover, 
                         linetype = "Mean cover", group = 1), 
                   fun = mean, geom = "line") 
  
  } else { 
    plot <- plot + 
      geom_area(aes_q(x = ~along, y = ~type_freq, fill = ~type), 
                position = "stack") + 
      geom_line(aes_q(x = ~along, y = ~percolation, linetype = "Perc. (full)"))  + 
      geom_line(aes_q(x = ~along, y = ~percolation_empty, linetype = "Perc. (empty)"), 
                color = 'black') +
      geom_line(aes_q(x = ~along, y = ~cover, linetype = "Mean cover"), 
                color = 'black') 
      
  }
  
  # Add plrange to the graph
  plot <- plot + 
      geom_point(aes_q(x = ~along, y = ~plrange, group = 1), 
                stat = "identity") + 
      geom_line(aes_q(x = ~along, y = ~plrange, group = 1), 
                stat = "identity") 
  
  return(plot)
}


#' @rdname patchdistr_sews_plot
#' 
# // along arg is already documented in plot() method
#' 
#' @param best_only Plot only the best fit the empirical (inverse cumulative) 
#'   patch-size distribution with an overlay of the estimated fits. 
#' 
#' @param plrange Plot the power-law range 
#'
#'@export
plot_distr <- function(x, along = NULL, best_only = TRUE, plrange = TRUE) { 
  UseMethod('plot_distr')
}

#'@export
plot_distr.patchdistr_sews <- function(x, along = NULL, best_only = TRUE, 
                                        plrange = TRUE) { 
  NextMethod('plot_distr')
}

#'@export
plot_distr.patchdistr_sews_single <- function(x, 
                                              along = NULL, 
                                              best_only = TRUE, 
                                              plrange = TRUE) { 
  
  if ( !is.null(along) ) { 
    warning('Ignoring along parameter when plotting only one patch-size', 
            'distribution')
  }
  
  # Get plottable data.frames
  values <- predict(x, best_only = best_only, xmin_rescale = TRUE)  
  
  # Create base plot 
  plot <- ggplot() + 
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") + 
    theme_spwarnings() +
    labs(x = 'Patch size',
         y = 'Observed frequency (P>=x)',
         color = 'Fitted type')
  
  # If plrange exists, then add it to the plot 
  plrange_dat <- x[['plrange']]
  if ( plrange && !is.na(plrange_dat[ ,"plrange"]) ) { 
    plrange_dat[ ,"xmax"] <- max(values[["obs"]][ ,"patchsize"])
    plot <- plot + 
      geom_segment(aes_q(x = ~xmin_est, y = 1, 
                         xend = ~xmax,  yend = 1), 
                   data = plrange_dat, 
                   arrow = arrow(ends = "both", type = "open", 
                                 length = unit(0.05, "inches")), 
                   color = "blue")
  }
  
  # Add observed values
  plot <- plot + 
    geom_point(aes_string(x = "patchsize", y = "y"), data = values[['obs']]) 
  
  # It can happen that no distribution have been fitted. Check for that 
  # before plotting otherwise it produces an error. 
  if ( any(!is.na(values[['pred']][ ,"y"])) ) { 
    
    # Get data and remove NAs (where no fit has been carried out
    pred_values <- values[["pred"]]
    pred_values <- pred_values[!is.na(pred_values[ ,"y"]), ]
    
    # Use only the y values of the fit that falls within the range of the observed 
    # values. Otherwise fit curves can go very low and we don't see a thing on the plot 
    pred_values <- pred_values[pred_values[ ,"y"] >= min(values[["obs"]]) & 
                               pred_values[ ,"y"] <= max(values[["obs"]]), ]
    
    if ( nrow(pred_values) > 0 ) { 
      plot <- plot + 
        geom_line(aes_string(x = "patchsize", y = "y", color = "type"), 
                  data = pred_values)
    }
    
  } else { 
    warning('No distribution has been fitted to the observed patch size distribution')
  }
  
  return(plot)
}

#'@export
plot_distr.patchdistr_sews_list <- function(x, 
                                            along = NULL, 
                                            best_only = TRUE, 
                                            plrange = TRUE) { 
  
  if ( !is.null(along) && (length(along) != length(x)) ) { 
    stop('The along values are unfit for plotting (size mismatch)')
  }
  
  # Get plottable data.frames
  values <- predict(x, best_only = best_only, xmin_rescale = TRUE)
  
  values[['obs']][ ,"along"]  <- values[['obs']][ ,'matrixn']
  values[['pred']][ ,"along"] <- values[['pred']][ ,'matrixn']
  # Modify matrixn column if necessary and reorder values
  if ( ! is.null(along) ) { 
    values[['obs']][ ,'along']  <- along[values[["obs"]][ ,'matrixn']]
    values[['pred']][ ,'along'] <- along[values[["pred"]][ ,'matrixn']]
  }
  
  plot <- ggplot() + 
    scale_y_log10() +
    scale_x_log10() + 
    facet_wrap( ~ along ) + 
    xlab('Patch size') + 
    ylab('Frequency (P>=x)') + 
    theme_spwarnings()
  
  if ( plrange ) { 
    # Add plrange to the plot. We need to extract info 
    # from the observed psd so that we can place the segment on the plot. 
    plrange_dat <- unique( as.data.frame(x)[ ,c("matrixn", 'xmin_est')] )
    plrange_dat[ ,"along"] <- plrange_dat[ ,"matrixn"]
    if ( ! is.null(along) ) { 
      plrange_dat[ ,"along"] <- along[plrange_dat[ ,"matrixn"]]
    }
    
    patches_minmax <- ddply(values[['obs']], "matrixn", 
                          function(df) { 
                            data.frame(xmax = max(df[ ,"patchsize"]))
                          })
    plrange_dat <- join(plrange_dat, patches_minmax, type = "left", 
                        match = "first", by = "matrixn")
    
    plot <- plot + 
      geom_segment(aes_q(x = ~xmin_est, y = 1, 
                         xend = ~xmax,  yend = 1), 
                    data = plrange_dat, 
                    arrow = arrow(ends = "both", type = "open", 
                                  length = unit(0.05, "inches")), 
                    color = "blue") 
  
  }
  
  # Add observed values
  plot <- plot + 
    geom_point(aes_string(x = "patchsize", y = "y"), 
               data = values[['obs']]) 
  
  # It can happen that no distribution have been fitted. Check for that 
  # before plotting otherwise it produces an error. 
  if ( any(!is.na(values[['pred']][ ,"y"])) ) { 
    # Get data and remove NAs (where no fit has been carried out
    pred_values <- values[["pred"]]
    pred_values <- pred_values[! is.na(pred_values[ ,"type"]) & 
                               ! is.na(pred_values[ ,"y"]), ]
    
    # Use only the y values of the fit that falls within the range of the observed 
    # values. Otherwise fit curves can go very low and we don't see a thing on the plot 
    pred_values <- pred_values[pred_values[ ,"y"] >= min(values[["obs"]]) & 
                               pred_values[ ,"y"] <= max(values[["obs"]]), ]
    if ( any(is.na(pred_values[ ,"type"])) ) browser()
    
    if ( nrow(pred_values) > 0 ) { 
      plot <- plot + 
        geom_line(aes_string(x = "patchsize", y = "y", color = "type"), 
                  data = pred_values)
    }
    
    plot <- plot + 
              geom_line(aes_string(x = "patchsize", 
                                   y = "y", 
                                   color = "type", 
                                   group = "paste(type, matrixn)"), 
                        data = pred_values)

  } else { 
    warning('No distribution has been fitted to the observed patch size distributions')
  }
  
  return(plot)
}





# Predict methods 
# --------------------------------------------------
#' @rdname patchdistr_sews_predict
#' @name patchdistr_sews_predict
#' 
#' @title predict method for patchdistr_sews objects
#' 
#' @description Export the observed and fitted patch size distributions 
#' 
#' @param object An \code{\link{patchdistr_sews}} object 
#' 
#' @param newdata A vector of patch sizes at which the fit is returned (default 
#'   to 200 regularly-spaced values). 
#' 
#' @param ... Additional arguments (ignored)
#' 
#' @param best_only Return values for only the best fit of each element (matrix)
#'   in \code{object}, or return the values for all fitted distribution. 
#' 
#' @param xmin_rescale If the xmin value used for fits is above one, then setting this
#'   to \code{TRUE} will rescale the predicted probabilities so that they align on 
#'   the cumulative distribution of the observed patch sizes 
#' 
#' @return A list with component obs, a data.frame containing the observed 
#'   distribution values and pred, a data.frame containing the fitted 
#'   values. 
#' 
#' @details The function \code{\link{patchdistr_sews}} fits competing 
#'   distribution models to the observed patch size distributions. This 
#'   functions is able to export the observed values and the fitted values 
#'   altogether. 
#' 
#' @examples 
#' 
#' \donttest{ 
#' patch_indics <- patchdistr_sews(forestgap)
#' 
#' predict(patch_indics)
#' 
#' }
#' 
#' @seealso \code{\link{patchdistr_sews}}
#' 
#'@export
predict.patchdistr_sews_single <- function(object, ..., 
                                           newdata = NULL,
                                           best_only = FALSE, 
                                           xmin_rescale = FALSE) { 
  
  # Get observed values
  vals_obs <- cumpsd(object[["psd_obs"]])
  
  # Shapes table
  shptbl <- object[['psd_type']]
  xmin_used <- shptbl[1, "xmin_fit"] # should be the same for all distribs
  stopifnot(length(unique(shptbl[ , "xmin_fit"])) == 1) # enforce it
  
  # Bail if no fit carried out. Note that we need to set classes in the 
  # output pred df otherwise coercion when merging with existing results 
  # goes wrong (e.g. factor converted to integer).
  if ( all(is.na(shptbl[ ,'type'])) ) { 
    result <- list(obs = vals_obs, 
                   pred = data.frame(type = NA_character_, 
                                     patchsize = NA_integer_, 
                                     y = NA_real_))
    return(result)
  }
  
  # Create x vector of values. We predict only values above xmin, as the distribution 
  # is undefined below that. 
  if ( is.null(newdata) ) { 
    newdata <- unique( round(10^(seq(log10(xmin_used), 
                                     log10(max(object[["psd_obs"]])), 
                                     length.out = 256))) )
  }
  
  if ( best_only ) { 
    shptbl <- shptbl[shptbl[ ,'best'], ]
  }
  
  # Get values
  vals_pred <- data.frame()
  for ( type in shptbl[ ,"type"] ) { 
    type_yvals <- switch(type, 
                         pl  = ippl(newdata, 
                                    shptbl[type, "plexpo"], 
                                    shptbl[type, "xmin_fit"]),
                         tpl = iptpl(newdata, 
                                     shptbl[type, "plexpo"], 
                                     shptbl[type, "cutoff"], 
                                     shptbl[type, 'xmin_fit']),
                         exp = ipdisexp(newdata,  
                                        shptbl[type, "cutoff"], 
                                        shptbl[type, "xmin_fit"]),
                         lnorm = ipdislnorm(newdata, 
                                            shptbl[type, "meanlog"], 
                                            shptbl[type, "sdlog"],  
                                            shptbl[type, "xmin_fit"]))
    
    vals_pred <- rbind(vals_pred, 
                       data.frame(type = type, 
                                  patchsize = newdata, 
                                  y = type_yvals))
  }
  
  # If xmin > 1, we need to rescale the values by a factor so that the probability 
  # matches the observed data in the cumulative distribution representation. This is 
  # mainly used when plotting the distributions. 
  if ( xmin_rescale ) { 
    # We use approx because there is no guarantee that the xmin used falls on an 
    # observed patch size
    scale_factor <- approx(vals_obs[ ,"patchsize"], 
                           vals_obs[ ,"y"], 
                           xmin_used)[["y"]]
    
    vals_pred[ ,"y"] <- vals_pred[ ,"y"] * scale_factor
  }
  
  # Return data
  return( list(obs = vals_obs, 
               pred = vals_pred) )
}

#'@export
predict.patchdistr_sews_list <- function(object, ..., 
                                         newdata = NULL, 
                                         best_only = FALSE, 
                                         xmin_rescale = FALSE) { 
  
  dat <- lapply(object, predict.patchdistr_sews_single, 
                newdata = newdata, best_only = best_only, 
                xmin_rescale = xmin_rescale)
  
  # Add id but handle when psd is empty
  add_id <- function(n, x) { 
    if (nrow(x) > 0) { 
      x <- data.frame(matrixn = n, x) 
    } else {
      x <- data.frame(matrixn = n, patchsize = NA_integer_)
    } 
  }
  
  dat <- Map(function(n, x) { 
               x[['obs']]  <- add_id(n, x[["obs"]])
               x[['pred']] <- add_id(n, x[["pred"]])
               x
             }, seq.int(length(dat)), dat)
  
  # Convert to data.frame
  dat <- list(obs  = plyr::ldply(dat, function(x) x[['obs']]), 
              pred = plyr::ldply(dat, function(x) x[['pred']]))
  
  return(dat)
}




# As data.frame methods
# --------------------------------------------------

#'@export
as.data.frame.patchdistr_sews_single <- function(x, ...) { 
  ans <- data.frame(x[['psd_type']], x[['plrange']],
                    # Add the name explicitely for these as they are single values
                    percolation = x[['percolation']], 
                    percolation_empty = x[['percolation_empty']], 
                    cover = x[['cover']], 
                    npatches = x[['npatches']],
                    unique_patches = x[['unique_patches']], 
                    stringsAsFactors = FALSE)
  # We convert the psd type to a character vector for simplicity
  ans[ ,'type'] <- as.character(ans[ ,'type'])
  ans
}

#'@export
as.data.frame.patchdistr_sews_list <- function(x, ...) { 
  
  # Format data
  results <- lapply(x, as.data.frame.patchdistr_sews_single)
  results <- Map(function(n, df) { 
                   data.frame(matrixn = n, df) 
                 }, seq.int(length(results)), results)
  
  # Bind it altogether and return df 
  do.call(rbind, results)
}




# Summary methods
# --------------------------------------------------

prepare_summary_table <- function(x, ...) { 
  
  dat <- as.data.frame(x)
  # Select lines that either have no fit or select the best fit
  dat <- dat[is.na(dat[ ,"best"]) | dat[ ,"best"], ]
  
  # Add formated unique patch column 
  dat[ ,'pretty_patches'] <- with(dat, 
                                  paste0(npatches, "(", unique_patches, ")"))
  # Add formated plrange column
  dat[ ,'plrange'] <- with(dat, 
                           ifelse(is.na(plrange), "", 
                                  paste0(round(plrange * 100), "%")))
  
  # Add formated distribution type column
  dat[ ,'type'] <- ifelse(is.na(dat[ ,'type']), " ", dat[ ,'type'])
  
  # Set the column names to extract
  cols <- c('pretty_patches', 'cover', 'percolation', 'percolation_empty',
            'type', 'plrange')
  pretty_names <- c('N(uniq.)', 'Cover', 'Percl.Full', 'Percl.Empt', 'Type', 'PLR')
  
  # If there is a matrixn column (for when several matrices were 
  #   used at once), then add it to these columns
  if ( "matrixn" %in% names(dat) ) { 
    cols <- c("matrixn", cols)
    pretty_names <- c('Mat. #', pretty_names)
  }
  
  # If there is a pvalue columne, then add it to the output table at the end, 
  # along with the stars
  if ( "plrpval" %in% names(dat) ) { 
    dat[ ,"pvalstars"] <- pval_stars(dat[ ,"plrpval"])
    dat[ ,'pvalstars'] <- ifelse(is.na(dat[ ,'plrpval']), " ", 
                                 as.character(dat[ ,'pvalstars']))
    dat[ ,'plrpval'] <- ifelse(is.na(dat[ ,'plrpval']), " ", 
                               as.character(dat[ ,'plrpval']))
    
    cols <- c(cols, "plrpval", "pvalstars")
    pretty_names <- c(pretty_names, "P<null (PLR)", "")
  }
  
  # Extract data and rename cols
  dat <- dat[, cols]
  names(dat) <- pretty_names
  
  return(dat)
}

#'@export
summary.patchdistr_sews_single <- function(object, ...) { 
  summary.patchdistr_sews_list(object, ...)
}
#'@export
summary.patchdistr_sews_list <- function(object, ...) { 
  dat <- prepare_summary_table(object)
  
  cat('Spatial Early-Warning: Patch-based indicators\n') 
  cat('\n')
  print.data.frame(dat, row.names = FALSE, digits = DIGITS)
  cat('\n')
  cat("The following methods are available: \n")
  cat(list_methods(class(object)), "\n")
  
  invisible(dat)
}



# Print methods
# --------------------------------------------------
# Note this method works for both single and lists objects so we only 
#   define it for patchdistr_sews (and not their *_list *_single 
#   equivalents)

#'@export
print.patchdistr_sews_single <- function(x, ...) { 
  summary.patchdistr_sews_single(x, ...)
}
#'@export
print.patchdistr_sews_list <- function(x, ...) { 
  summary.patchdistr_sews_list(x, ...)
}



# Helper function 
# ---------------
# Get the inverse cumulative distribution of a psd (P(x >= k)). x here are the 
# values at which to evaluate the inverse cumulative psd. 
cumpsd <- function(dat, x = sort(unique(dat))) { 
  N <- length(dat)
  y <- sapply(x, function(k) { sum(dat >= k) / N })
  return( data.frame(patchsize = x, y = y) )
}


