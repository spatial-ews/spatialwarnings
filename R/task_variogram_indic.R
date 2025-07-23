# 
# Variogram-based indicators. 
# 
#' @title Early-Warning signals based on variograms (EXPERIMENTAL)
#' 
#' @description Compute Early-warning signals based on metrics derived 
#'   form semi-variograms. 
#' 
#' @param mat A matrix (TRUE/FALSE values) or a list of matrices 
#' 
#' @param model The variogram model to use, either "sph" (for a spherical
#'   model) or "exp" (for an exponential model)
#' 
#' @param nmax The maximum number of pairs of cells to use when computing the 
#'   variogram
#' 
#' @param nbins Number of distance bins to use to compute the variogram
#' 
#' @param cutoff Maximum distance to consider in the variogram. If NULL, then 
#'   a distance equal to one third of the diagonal of the matrix is used
#' 
#' @return A list object of class "variogram_sews", that can be displayed 
#'   using \code{summary()}, \code{plot()}, etc. Significance of values can 
#'   be tested using \code{\link{indictest}}. 
#' 
#' @details 
#' 
#' During ecosystem degradation and especially before a regime shift occurs in 
#'   some ecosystems, spatial autocorrelation is expected to increase in a
#'   landscape. This increase can be measured based on variograms, which
#'   represent how the difference (variance) between two points in a landscape
#'   varies as a function of distance. 
#' 
#' The approach used to derive variogram-based EWS is to compute the 
#'   empirical variogram of a landscape (represented passed as a matrix of 
#'   values), then fit a variogram model to it. Three 
#'   parameters are then extracted from the variogram model (see Nijp et al.
#'   2019 for a visual description of these parameters): 
#'   \enumerate{ 
#'     \item The \emph{nugget} (intercept)
#'     \item The \emph{partial sill}, i.e. the reduction in semivariance at 
#'       distance zero
#'     \item The \emph{correlation range}, i.e. the distance at which the 
#'       relationship between semivariance and distance flattens
#'   }
#' 
#' Additionally, the \emph{structural variance} is computed as 
#'   (partial sill)/(nugget + partial sill), wich quantifies whether the 
#'   data are spatially structured (structural variance of one), or completely 
#'   unstructured (value of zero). Theoretical work suggests that partial sill,
#'   correlation range and structural variance should increase before a regime
#'   shift occurs in an ecosystem (Nijp et al. 2019). 
#' 
#' This function offers to fit a spherical model or 
#'   an exponential model. The best-fitting model depends on your data, you 
#'   should try different options and review the fits using
#'   \code{\link{plot_variogram}}.
#' 
#' Please note that this part of the package is still experimental and deserves 
#'   more testing. 
#' 
#' @seealso \code{\link{raw_structvar}},
#'   \code{\link{plot_variogram}}, \code{\link{extract_variogram}},
#'   \code{\link[=predict.variogram_sews_list]{predict}} 
#' 
#' @seealso 
#'   \code{\link{indictest}}, to test the significance of indicator values. 
#' 
#' @references 
#' 
#' Nijp, Jelmer J., Arnaud J.A.M. Temme, George A.K. Voorn, Lammert Kooistra,
#'   Geerten M. Hengeveld, Merel B. Soons, Adriaan J. Teuling, and Jakob
#'   Wallinga. (2019) Spatial Early Warning Signals for Impending Regime Shifts:
#'   A Practical Framework for Application in Real-world Landscapes. Global
#'   Change Biology 25 (6): 1905-21. \doi{10.1111/gcb.14591}
#'
#' @examples 
#' 
#' \donttest{
#' serengeti_ews <- variogram_sews(serengeti, 
#'                                 model ="exp")
#' plot(serengeti_ews, along = serengeti.rain)
#' summary(serengeti_ews)
#' 
#' plot_variogram(serengeti_ews)
#' 
#' # nulln should be set to a higher values for meaningful results
#' serengeti_test <- indictest(serengeti_ews, nulln = 9)
#' plot(serengeti_test) # gray ribbons indicate the null indicator values 
#' summary(serengeti_test)
#' }
#' 
#'@export
variogram_sews <- function(mat, 
                           model = "sph", 
                           nmax = 1e5, 
                           nbins = 32, 
                           cutoff = NULL) { 
  
  if ( ! model %in% c("sph", "exp") ) { 
    stop(paste0("The variogram model must be one of 'sph' (spherical model) ", 
                "or 'exp' (exponential model)"))
  }
  
  # Handle list case
  if ( is.list(mat) ) { 
    results <- lapply(mat, variogram_sews, model, nmax, nbins, cutoff)
    names(results) <- names(mat)
    class(results) <- c('variogram_sews_list',  
                        'simple_sews_list', 
                        'sews_result_list')
    attr(results, 'indicname') <- "Variogram-based indicators"
    return(results)
  }
  
  # Returns a list with components 'vario' (the empirical variogram) and 
  # 'fit' (its fit)
  vario <- fit_variogram(mat, model = model, nmax = nmax, nbins = nbins, 
                         cutoff = cutoff)
  metrics <- compute_vario_metrics(vario[['pars']])
  
  # Return list containing both
  output <- list(variogram = vario[['vario']], 
                 fit = vario[['pars']], 
                 value = metrics, 
                 orig_data = mat, 
                 pars = list(model = model, 
                             nmax = nmax, 
                             nbins = nbins, 
                             cutoff = cutoff), 
                 taskname = "Variogram-based indicators")
  
  class(output) <- c('variogram_sews_single', 
                     'simple_sews_single', 
                     'sews_result_single')
  attr(output, 'indicname') <- "Variogram-based indicators"
  
  return(output)
}

# Extract the variogram as a data frame from a variogram_sews object
#' 
#' @title extract_variogram() method for variogram_sews objects
#' 
#' @description Extract the empirical variogram from a \code{variogram_sews} 
#'   object
#' 
#' @param x An object produced by \code{variogram_sews}, or
#'   \code{indictest}
#' 
#' @param ... Additional arguments (ignored)
#' 
#' @return A data.frame containing the variogram with the distances 
#'   (column \code{dist}), the empirical semivariance values (\code{gamma}), 
#'   and if object contains more than one matrix, a column named \code{matrixn}. 
#' 
#' @examples 
#' 
#' \donttest{ 
#' vario_indics <- variogram_sews(serengeti)
#' predict(vario_indics)
#' vario_test <- indictest(vario_indics, nulln = 19)
#' predict(vario_test) # same result
#' }
#' 
#' @seealso \code{\link{variogram_sews}}
#' 
#'@export 
extract_variogram <- function(x, ...) { 
  UseMethod("extract_variogram")
}
#'@method extract_variogram variogram_sews_list
#'@export
extract_variogram.variogram_sews_list <- function(x, ...) { 
  all_varios <- Map(function(n, o) { 
    data.frame(matrixn = n, o[["variogram"]]) 
  }, n = seq_along(x), o = x)
  do.call(rbind, all_varios)
}
#'@method extract_variogram variogram_sews_single
#'@export
extract_variogram.variogram_sews_single <- function(x, ...) { 
  x[["variogram"]]
}



# Predict methods 
# --------------------------------------------------
#' @rdname variogram_sews_predict
#' @name variogram_sews_predict
#' 
#' @title predict() method for variogram_sews objects
#' 
#' @description Export the fitted variogram(s)
#' 
#' @param object An object produced by variogram_sews 
#' 
#' @param newdist A vector of distances at which to return the variogram fit 
#'   values (defaults to 128 regularly-spaced values). 
#' 
#' @param ... Additional arguments (ignored)
#' 
#' @return A data.frame with the distances (column \code{dist}), the fitted 
#'   values (\code{gamma}), and if object contains more than one matrix, 
#'   a column \code{matrixn}. 
#' 
#' @examples 
#' 
#' \donttest{ 
#' vario_indics <- variogram_sews(serengeti)
#' predict(vario_indics)
#' vario_test <- indictest(vario_indics, nulln = 19)
#' predict(vario_test) # same result
#' }
#' 
#' @seealso \code{\link{variogram_sews}}
#' 
#'@export 
predict.variogram_sews_list <- function(object, newdist = NULL, ...) { 
  all_preds <- llply(object, predict.variogram_sews_single, newdist = newdist)
  all_preds <- Map(function(n, o) { 
    data.frame(matrixn = n, o)
  }, seq_along(all_preds), all_preds)
  do.call(rbind, all_preds)
}
#'@export 
predict.variogram_sews_single <- function(object, newdist = NULL, ...) { 
  if ( is.null(newdist) ) { 
   newdist <- with(object[["variogram"]], seq(0, max(dist), length.out = 128))
  }
  pars <- object[["value"]]
  if ( object[["pars"]][["model"]] == "sph" ) { 
    y <- spherical_model(newdist, pars[1], pars[2], pars[3])
  } else if ( object[["pars"]][["model"]] == "exp" ) { 
    y <- exponential_model(newdist, pars[1], pars[2], pars[3])
  } else { 
    stop("Unknown model")
  }
  data.frame(dist = newdist, gamma = y)
}

# Plot methods for variogram sews 
#' 
#' @rdname variogram_sews_plot
#' @name variogram_sews_plot
#' 
#' @title Early-warning signals based on variograms 
#' 
#' @description Plot trends of indicators based on variograms
#' 
#' @param x An object produced by \code{\link{variogram_sews}}, or the 
#'   result of applying \code{indictest} on such object.
#' 
#' @param along A vector providing values along which the indicator trends 
#'   will be plotted. If \code{NULL} then the indicator values are plotted
#'   sequentially in their original order. 
#' 
#' @details The \code{plot()} function will display how the estimated 
#'   variogram parameters change along a set of values (passed with argument
#'   \code{along}). If the object passed has been processed through 
#'   \code{indictest}, then the null values are also displayed. 
#'   \code{plot_variogram()} can be used to display the individual variograms
#'   that have been fit to the data. 
#' 
#' @seealso 
#' 
#' \code{\link{variogram_sews}}, \code{\link{indictest}}, 
#'   \code{\link{plot_variogram}}  
#' 
#' @examples 
#'   
#' \donttest{ 
#' serengeti_ews <- variogram_sews(serengeti, model ="exp")
#' 
#' # Display the change in variogram parameters 
#' plot(serengeti_ews, along = serengeti.rain) + 
#'   ggplot2::labs(x = "Rainfall (mm)")
#' 
#' # Visualize the fitted variograms
#' plot_variogram(serengeti_ews, along = serengeti.rain) 
#' 
#'   # Test the trends (nulln should be set to a higher value to obtain 
#'   # meaningful results
#'   serengeti_test <- indictest(serengeti_ews, nulln = 19)
#'   plot(serengeti_test, along = serengeti.rain)
#'   plot_variogram(serengeti_test, along = serengeti.rain)
#' }
#' 
#'@export
plot.variogram_sews <- function(x, along = NULL, ...) { 
  NextMethod("plot")
}
#'@export 
plot.variogram_sews_list <- function(x, along = NULL, ...) { 
  # Adjust x, then pass it on to the method for simple_sews plotting. Here 
  # we use a for loop so that the class attributes is preserved. 
  for (i in seq_along(x) ) { 
    x[["value"]] <- x[["metrics"]]
  }
  
  plot.simple_sews_list(x, along = along, ...)
}

#'@rdname variogram_sews_plot
#'@export
plot_variogram <- function(x, along = NULL, ...) { 
  UseMethod("plot_variogram")
}
#'@rdname variogram_sews_plot
#'@method plot_variogram variogram_sews
#'@export
plot_variogram.variogram_sews <- function(x, along = NULL, ...) { 
  UseMethod("plot_variogram")
}
#'@export
plot_variogram.variogram_sews_single <- function(x, along = NULL, ...) { 
  df <- extract_variogram(x)
  df_pred <- predict(x)
  
  ggobj <- ggplot(df, aes(x = dist, y = gamma)) + 
    geom_point() + 
    geom_line() + 
    geom_line(color = "red", data = df_pred) + 
    theme_spwarnings() + 
    linescale_spwarnings() + 
    labs(x = "Distance", 
         y = "gamma")
  
  return(ggobj)
}
#'@export
plot_variogram.variogram_sews_list <- function(x, along = NULL, ...) { 
  
  # Extract data and add along data
  df <- extract_variogram(x)
  df_pred <- predict(x)
  if ( is.null(along) ) { 
    df[ ,'along'] <- df[ ,"matrixn"]
    df_pred[ ,'along'] <- df_pred[ ,"matrixn"]
  } else { 
    df[ ,"along"] <- along[df[ ,"matrixn"]]
    df_pred[ ,"along"] <- along[df_pred[ ,"matrixn"]]
  }
  
  ggobj <- ggplot(df, aes_string(x = "dist", y = "gamma")) + 
    geom_point() + 
    geom_line() + 
    geom_line(color = "red", data = df_pred) + 
    theme_spwarnings() + 
    linescale_spwarnings() + 
    labs(x = "Distance", 
         y = "gamma") + 
    facet_wrap( ~ along, scales = "free_y")
  
  return(ggobj)
}


# Fit a variogram. We do it by hand and do not use gstat which seems to be 
# a bit unreliable. We initialize the varioram with a segmented regression 
# first, then refine that fit with bounded BFGS. 
fit_variogram <- function(mat, model, nmax, nbins, cutoff) { 
  
  # Handle missing arguments 
  if ( is.null(cutoff) ) { 
    cutoff <- (1/3) * sqrt(nrow(mat)^2 + ncol(mat)^2) 
  }
  
  # Make the model choice case-insensitive
  model <- tolower(model)
  
  # Check if mat is suitable
  check_mat(mat)
  
  # Compute the variogram 
  vario <- variogram_internal(mat, nmax, nbins, cutoff)
  
  # Check if sample has only a single type of values. In that case we return 
  # a default answer with missing values for parameters as the model will not 
  # fit.
  if ( length(unique(vario[ ,"gamma"])) < 2 ) { 
    default_answer <- 
      list(pars = list(nugget = NaN, 
                       psill  = NaN, 
                       range  = NaN), 
           vario = vario)
    return(default_answer)
  }
  
  # Fit variogram. We first use a segmented regression to find a breakpoint in 
  # the empirical variogram, then use that as starting values for correlation 
  # range. Note that the segmented regression can fail, especially when we start giving
  # null data to it, so we wrap it in a try() block.
  brk_model <- try(suppressWarnings({ 
    segmented::segmented(lm(gamma ~ dist, 
                            data = vario))
  }), silent = TRUE)
  
  range0 <- 1/3 * (max(vario[ ,"dist"]) - min(vario[ ,"dist"]))
  psill0 <- 0
  nugget0 <- 0
  if ( ! inherits(brk_model, "try-error") ) { 
    brk_coefs <- coef(brk_model)
    # Extract the break value, if there is one
    nugget0 <- max(0, as.vector(brk_coefs["(Intercept)"]))
    if ( ! is.null(brk_model[["psi"]][2]) ) { 
      range0 <- brk_model[["psi"]][2]
      psill0 <- abs(nugget0 - median(vario[vario$dist > range0, "gamma"]))
    } 
  } 
  
  pars0 <- c(nugget = nugget0, psill  = psill0, range  = range0)
  
  # Objective function
  ssq <- function(pars0) { 
    if ( model == "sph" ) { 
      pred <- spherical_model(vario[ ,"dist"], pars0[1], pars0[2], pars0[3]) 
    } else if ( model == "exp" ) { 
      pred <- exponential_model(vario[ ,"dist"], pars0[1], pars0[2], pars0[3]) 
    }
    sum( ( (vario[ ,"gamma"] - pred) )^2 )
  }
  
  fit <- optim_safe(ssq, pars0,
                    # sensible bounds for nugget, psill, range
                    lower = c(0, 0, 0),
                    upper = c(max(vario[ ,"gamma"]),
                              max(vario[ ,"gamma"]),
                              2*max(vario[ ,"dist"]))
                    )

  vario[ ,"pred"] <- with(as.list(fit$par), 
                          spherical_model(vario[ ,"dist"], 
                                          nugget, psill, range))
#   with(vario, plot(dist, gamma))
#   with(vario, lines(dist, pred, col = "blue"))
#   with(vario, lines(dist, 
#                     spherical_model(dist, pars0[1], pars0[2], pars0[3]), 
#                     col = "red"))
  
  list(pars = as.list(fit$par), 
       vario = vario)
}

# Function that wraps C++ code to compute the variogram
variogram_internal <- function(mat, nmax, nbins, cutoff) { 
  
  v <- variogram_internal_cpp(mat, nmax, nbins, cutoff)
  colnames(v) <- c("dist", "gamma", "np")
  v <- as.data.frame(v)
  
  # Discard all values at small distances which have zero gamma
  line <- 1
  while ( is.na(v[line, "gamma"]) || v[line, "gamma"] == 0 ) { 
    line <- line + 1 
  }
  v[ v[ ,"np"] > 0 & seq.int(nrow(v)) >= line, ]
}

spherical_model <- function(X, nugget, psill, range) { 
  ifelse(X < range, 
         nugget + psill * ( (3/2) * (X/range) - (1/2) * (X/range)^3 ), 
         nugget + psill)
}

exponential_model <- function(X, nugget, psill, range) { 
  nugget + psill * (1 - exp(-X/range))
}

compute_vario_metrics <- function(pars) { 
  structvar <- with(pars, psill / (nugget + psill))
  metrics <- with(pars, c(nugget, psill, range, structvar))
  names(metrics) <- c('nugget', 'psill', 'range', 'structvar')
  return(metrics)
}

# Raw version of the function
#' 
#' @title Variogram parameters
#' 
#' @description Compute the nugget, partial sill, correlation range and 
#'   structural variance on a matrix. 
#' 
#' @param mat A matrix 
#' 
#' @param model The variogram model to use, either "sph" (for a spherical
#'   model) or "exp" (for an exponential model)
#' 
#' @param nmax The maximum number of pairs of cells to use when computing the 
#'   variogram
#' 
#' @param nbins Number of distance bins to use to compute the variogram
#' 
#' @param cutoff Maximum distance to consider in the variogram. If NULL, then 
#'   a distance equal to one third of the diagonal of the matrix is used
#' 
#' @examples 
#' 
#' raw_variogram_metrics(serengeti[[5]])
#' 
#' @seealso variogram_sews, raw_structvar
#' 
#'@export
raw_variogram_metrics <- function(mat, 
                                  model = "sph", 
                                  nmax = 1e5L, 
                                  nbins = 32, 
                                  cutoff = NULL) { 
  if ( ! is.matrix(mat) ) { 
    stop("raw_variogram_metrics only accepts a single matrix as input.")
  }
  
  vario <- fit_variogram(mat, model, nmax, nbins, cutoff)
  compute_vario_metrics(vario[["pars"]])
}

#' @title Structural variance
#' 
#' @description Compute the structural variance on a matrix. 
#' 
#' @param mat A matrix 
#' 
#' @param model The variogram model to use, either "sph" (for a spherical
#'   model) or "exp" (for an exponential model)
#' 
#' @param nmax The maximum number of pairs of cells to use when computing the 
#'   variogram
#' 
#' @param nbins Number of distance bins to use to compute the variogram
#' 
#' @param cutoff Maximum distance to consider in the variogram. If NULL, then 
#'   a distance equal to one third of the diagonal of the matrix is used
#' 
#' @seealso raw_variogram_metrics, variogram_sews
#' 
#' @examples 
#' 
#' raw_structvar(serengeti[[5]])
#' 
#'@export
raw_structvar <- function(mat, 
                          model = "sph", 
                          nmax = 1e5L, 
                          nbins = 32, 
                          cutoff = NULL) { 
  raw_variogram_metrics(mat, model, nmax, nbins, cutoff)["structvar"]
}
