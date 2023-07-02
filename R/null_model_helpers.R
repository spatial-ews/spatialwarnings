# 
# These functions factorize some code between numerical indicators. Most 
# notably it handles creating the null model testing, etc.
# 
# This function will compute a null distribution of indicator values on input, 
# by reshuffling it nulln times. 
# 
compute_indicator_with_null <- function(input, 
                                        nulln, 
                                        indicf, 
                                        null_method, 
                                        null_control) { 
  
  # Create null_control list with default arguments 
  null_control <- null_control_set_args(input, null_control, null_method)
  
  # Compute the observed value
  value  <- indicf(input)
  result <- list(value = value)
  
  if ( nulln > 2 ) { 
    
    null_results <- generate_nulls(input, indicf, nulln, 
                                   null_method, null_control)
    
    # nulldistr is a list so far => combine it to a matrix, with each null 
    # replicate as a row
    null_results[["nulldistr"]] <- do.call(rbind, null_results[["nulldistr"]]) 
    
    # Compute summary stats
    summ_stats <- 
      list(null_mean = apply(null_results[["nulldistr"]], 2, mean), 
           null_sd   = apply(null_results[["nulldistr"]], 2, sd), 
           null_qsup = apply(null_results[["nulldistr"]], 2, safe_quantile,
                             null_control[["qsup"]]), 
           null_qinf = apply(null_results[["nulldistr"]], 2, safe_quantile,
                             null_control[["qinf"]]), 
           z_score   = apply(rbind(value, null_results[["nulldistr"]]), 2, 
                             function(X) (X[1] - mean(X[-1])) / sd(X[-1])),
           pval      = apply(rbind(value, null_results[["nulldistr"]]), 2, 
                             function(X) 1 - rank(X)[1] / length(X)))
    
    # If null_method is a function, we need to convert it to a sensible string 
    # for display in summary(), etc. 
    if ( is.function(null_method) ) { 
      null_method <- "custom function"
    }
    
    # Add to that list what is returned from null + some stats
    result <- list(nulldistr = null_results[["nulldistr"]], 
                   summary_values = summ_stats, 
                   info = list(null_method = null_method, 
                               nulln = nulln, 
                               get_nullmat = null_results[["get_nullmat"]]))
    
  }
  
  return(result)
}

# Returns a list with the values of indicf obtained from random matrices + 
# a function used to generate the nulls
generate_nulls <- function(input, indicf, nulln, null_method, 
                           null_control) {
                             
  if ( ( !is.function(null_method) ) && 
       ( ! null_method %in% c("perm", "intercept", "smooth") ) ) { 
    stop(paste0("Unknown null model method: ", null_method))
  }
  
  # Nulldistr is null as long as we've not computed it yet
  has_computed_nulldistr <- FALSE
  
  
  
  
  # If the user provided a function to generate null matrices
  if ( is.function(null_method) ) { 
    
    # We warn if the function does not return logical values 
    if ( is.logical(input) && ! is.logical(null_method(input)) ) { 
      warning("The original data is logical (TRUE/FALSE), but the null ", 
              "method does not return logical values")
    }
    
    null_mod <- NULL # No model involved when the function is provided
    get_nullmat <- function() { null_method(input) }
  }
  
  
  
  
  # If the null method is set to permutation
  if ( is.character(null_method) && null_method == "perm" ) { 
    # Compute the index on a randomized matrix. Here we use a dedicated 
    # function that can do the shuffling in-place in C++, instead of relying 
    # on R which copies matrices in memory many times. 
    # 
    # shuffle_and_compute will convert all matrices to numeric matrices 
    # internally. We need to explicitely convert back to logical after 
    # shuffling before computing the indicator. 
    if ( is.logical(input) ) { 
      nulldistr <- shuffle_and_compute(input, function(x) indicf(x>0), 
                                       nulln)
    } else { 
      nulldistr <- shuffle_and_compute(input, indicf, nulln)
    }
    has_computed_nulldistr <- TRUE
    null_mod <- NULL # No model involved when we are permuting
    get_nullmat <- function() { 
      if ( is.logical(input) ) { 
        a <- shuffle_matrix(input) > 0.5 
      } else { 
        a <- shuffle_matrix(input) 
      }
      # Pass attributes to the new, random matrix
      attributes(a) <- attributes(input)
      return(a)
    }
  }
  
  
  
  
  # If the null method is to fit an intercept only 
  if ( is.character(null_method) && null_method == "intercept" ) { 
    
    values <- as.vector(input)
#     sub <- select_subset(length(values), null_control[["model_subset"]], 
#                          min(length(values), 512))
    null_mod <- glm(values ~ 1, family = null_control[["family"]])
    get_nullmat <- create_nullmat_generator(input, 
                                            null_mod, 
                                            null_control[["family"]])
  }
  
  
  
  # If the null method is to fit a smoothing spline to the matrix
  if ( is.character(null_method) && null_method == "smooth" ) { 
    
    if ( ! requireNamespace("mgcv", quietly = TRUE) ) { 
      stop("The 'gam' method requires mgcv. Please install it beforehand.")
    }
    
    # We fit a smoothing isotropic spline over the matrix, with automatic
    # parameter selection for the spline. 
    mat_tab <- data.frame(expand.grid(row = seq.int(nrow(input)), 
                                      col = seq.int(ncol(input))), 
                          value = as.vector(input))
#     sub <- select_subset(nrow(mat_tab), 
#                          null_control[["model_subset"]], 
#                          min(nrow(mat_tab), 512))
#     mat_tab <- mat_tab[sub, ]
    null_mod <- mgcv::gam(value ~ s(row, col, bs = "tp"), 
                          data = mat_tab, 
#                           control = mgcv::gam.control(trace = TRUE), 
                          family = null_control[["family"]])
    get_nullmat <- create_nullmat_generator(input, 
                                            null_mod, 
                                            null_control[["family"]])
  }
  
  
  
  # Compute the actual indicator values in the null distribution (if that's 
  # not already done)
  if ( ! has_computed_nulldistr ) { 
    nulldistr <- lapply(seq.int(nulln), function(n) { 
      indicf(get_nullmat())
    })
  }

  list(nulldistr   = nulldistr, 
       get_nullmat = get_nullmat, 
       null_mod    = null_mod)
}

# We use a safe version of quantile that reports as warnings 
# the appearance of NAs in the null distribution.
safe_quantile <- function(nulldistr, p) { 
  if ( any( is.na(nulldistr) ) ) { 
    warning(paste0('Computation of null values produced NAs (', 
                   sum(is.na(nulldistr)), " out of ", length(nulldistr), "). "))
  }
  quantile(nulldistr, p, na.rm = TRUE)
}

# Returns TRUE if fam is a binomial() object
is.binomial <- function(fam) { 
  ( !is.list(fam) && fam == "binomial") || 
    ( inherits(fam, "family") && fam[["family"]] == "binomial" )
}

# Set default arguments for null methods 
null_control_set_args <- function(mat, arglist, null_method) { 
  
  # Choose a sensible default if family is unset. Note that if null_method can 
  # be a function that the user provided. 
  if ( ! "family" %in% names(arglist) ) { 
    
    if ( is.logical(mat) ) { 
      family <- binomial()
      if ( is.character(null_method) && 
           null_method %in% c("smooth", "intercept") ) { 
        warning("indictest: using a binomial() family with default options")
      }
      
    } else { 
      family <- gaussian()
      if ( is.character(null_method) && 
           null_method %in% c("smooth", "intercept") ) { 
        warning("indictest: using a gaussian() family with default options")
      }
    }
    
  } 
  
  args <- list(family = family, 
               qinf = .05, 
               qsup = .95, 
               model_subset = 0.1) # Add other arguments here
  
  for ( i in seq_along(arglist) ) { 
    if ( names(arglist)[i] %in% names(args) ) { 
      args[[ names(arglist)[i] ]] <- arglist[[i]]
    } else { 
      warning("Unknown null model control argument: ", names(arglist)[i])
    }
  }
  
  return(args)
}

# Returns a function that generates matrices, given a null model and the 
# original matrix
create_nullmat_generator <- function(mat, null_mod, family) { 
  function() { 
    newdat <- expand.grid(row = seq.int(nrow(mat)), 
                          col = seq.int(ncol(mat)))
    
    # When the family is binomial, simulate often returns 1/0 instead of 
    # TRUE/FALSE values so we need to convert it back here. 
    sim <- matrix(stats::simulate(null_mod, newdata = newdat)[ ,1], 
                  nrow = nrow(mat), ncol = ncol(mat)) 
    if ( is.binomial(family) ) { 
      sim <- sim > .5
    }
    # Transfer matrix attributes
    attributes(sim) <- attributes(mat)
    return(sim)
  }
}

# # Select a subset of values in N values We take the maximum of at_least points 
# # and the fraction of the N values asked for
# select_subset <- function(N, subset_frac, at_least) { 
#   if ( N <= at_least ) { 
#     return( rep(TRUE, N) )
#   }
#   
#   keep_every <- min(floor(N / at_least), 
#                     floor(1 / subset_frac))
#   seq.int(N) %% keep_every == 0
# }

# 
# # Here newdata is assumed to have columns and x, y
# simulate_newdat <- function(mod, newdata, family) { 
#   
#   # If we passed a vector of values, convert that to a data.frame so it is 
#   # compatible with predict()
#   if ( is.vector(newdata) && !is.data.frame(newdata) ) { 
#     newdata <- data.frame(blank = newdata)
#   }
#   n <- nrow(newdata)
#   
#   # Extract the string from family object if needed
#   fam <- family$family
#   linv <- family$linkinv 
#   
#   if ( fam == "gaussian" ) { 
#     mu_pred <- predict(mod, newdata, type = "response")
#     resp_pred <- rnorm(n, mu_pred, sigma(mod))
#   } else if ( fam == "binomial" ) { 
#     mu_pred <- predict(mod, newdata, type = "response")
#     resp_pred <- rbinom(n, size = 1, prob = mu_pred)
#   } else if ( fam == "poisson" ) { 
#     lambda_pred <- predict(mod, newdata, type = "response")
#     resp_pred <- rpois(n, lambda = lambda_pred)
#   } else { 
#     stop("Family ", fam, " is not handled by this function")
#   }
#   
#   return(resp_pred)
# }
