# 
# This file contains code related to patch size distribution fitting. These 
#   functions can fit Power-law (pl), Truncated Power-law (tpl), Lognormal 
#   (lnorm) and Exponential (exp) distributions using maximum likelihood, as 
#   per Clauset et al. 's (2007) recommandations.
# 
# In addition, it provides the estimation of xmin using the ks-distance for 
# power-laws (same, following Clauset et al. 2007)

# Optimisation global options
ITERLIM <- 10000
GRADTOL <- 1e-10 
STEPTOL <- 1e-10
STEPMAX <- 5

ifnotfinite <- function(x, otherwise = .Machine$double.xmax) { 
  ifelse(is.finite(x), x, sign(x) * otherwise)
}

warnbound <- function() { 
  getOption("spatialwarnings.debug.fit_warn_on_bound", default = FALSE)
}

warnNA <- function() { 
  getOption("spatialwarnings.debug.fit_warn_on_NA", default = FALSE)
}

# This is a safe version of nlm that returns a sensible result (NaNs) when 
# the algorithm fails to converge. This can happen quite often when looking 
# for pathological cases (e.g. fitting distribution based on few points in the 
# tails, etc.). 
optim_safe <- function(f, pars0, 
                       lower = rep(-Inf, length(pars0)), 
                       upper = rep(Inf,  length(pars0)), 
                       fit_on_logscale = FALSE, ...) { 
  
  
  if ( fit_on_logscale ) { 
    optimf <- function(pars) { 
      f(exp(pars))
    }
    pars0 <- log(pars0)
    lower <- suppressWarnings(ifnotfinite(log(lower)))
    upper <- suppressWarnings(ifnotfinite(log(upper)))
  } else { 
    optimf <- f
  }
  
  optiresult <- try({ 
    optim(pars0, optimf, 
          control = list(maxit = ITERLIM), 
          lower = lower, upper = upper, 
          method = "L-BFGS-B", ...)
  }, silent = TRUE)
  
  # Code results above 3 means a true problem, below 3 the 
  # solution is either exact or approximate. 
  # Sometimes L-BFGS-B gets stuck in a very flat area just because our initial 
  # guess was very good, and L-BFGS-B does not like that (it reports 
  # abnormal termination of line search). Or L-BFGS-B will try the bounds of 
  # the parameter space, which will return error. So here we try to use BFGS 
  # which should report success if our initial guess was not too bad. 
  if ( inherits(optiresult, "try-error") || optiresult[["convergence"]] > 3 ) { 
    optiresult_bfgs <- try({ 
      optim(pars0, optimf, 
            control = list(maxit = ITERLIM), 
            method = "BFGS", ...)
    }, silent = TRUE)
    
    # If success, go with BFGS results
    if ( inherits(optiresult_bfgs, "try-error") ) { 
      optiresult <- optiresult_bfgs
    }
  }
  
  # If we could not reach a proper solution, report an error
  if ( inherits(optiresult, "try-error") ) { 
    optiresult <- list(value  = NaN, 
                       par = rep(NaN, length(pars0)), 
                       convergence = 128) 
    return(optiresult)
  }
  
  if ( optiresult[["convergence"]] > 3 ) { 
    warning(paste0('optim returned an error (error code:', 
                   optiresult[["convergence"]], ").\n", 
                   "Make sure the results are reasonable using plot_distr"))
  }
  
  # Convert the estimated pars back from log scale 
  if ( fit_on_logscale ) { 
    optiresult[["par"]] <- exp(optiresult[["par"]])
  }
  
  return(optiresult)
}

# Bounds on parameters, these should be large and no observed distribution should 
# have values beyond them
# Power-laws lambdas
PLMIN <- 1 + sqrt(.Machine$double.eps)
PLMAX <- 10
# Exponential rates
EXPMIN <- sqrt(.Machine$double.eps) # A very close value to, but not, zero
EXPMAX <- 10
# Bounds for truncated power-laws
TPL_EXPOMIN <- -1 # Taken from Clauset's code
TPL_EXPOMAX <- 10
TPL_RATEMIN <- sqrt(.Machine$double.eps)
TPL_RATEMAX <- 10

# Riemann zeta function with xmin taken into account :
# sum( 1/k^-expo ) for i=xmin to i = inf
# This is vectorized over xmins so that we do not sum things several times. 
zeta_w_xmin <- function(expo, xmins) { 
  perm <- order(xmins)
  xmins <- xmins[order(xmins)]
  
  # Compute zeta value
  zetaval <- gsl::zeta(expo)
  
  # Initialize
  output <- rep(NaN, length(xmins))
  current_k <- xmins[1]
  output[perm[1]] <- zetaval - sum_all_one_over_k(from = 1, to = xmins[1], expo)
  
  # If there is only one value, we bail now
  if ( length(xmins) <= 1) { 
    return(output)
  }
  
  for ( i in 2:length(xmins)) { 
    next_k <- xmins[i]
    if (next_k > current_k) { 
      output[perm[i]] <- output[perm[i-1]] - 
                           sum_all_one_over_k(from = current_k, 
                                              to = next_k, expo)
      current_k <- next_k
    } else { 
      output[perm[i]] <- output[perm[i-1]]
    }
  }
  
  return(output)
}




# PL fitting 
# ---------------------------------------

# Normalizing constant for pl with xmin
displnorm <- function(expo, xmin) { 
  # Adjust constant for threshold (note that this has no effect if xmin == 1, 
  #   as expected)
  const <- gsl::zeta(expo)
  const - sum_all_one_over_k(from = 1, to = xmin, expo)
}

# PL: P(x=k)
dpl <- function(x, expo, xmin = 1, log = FALSE) { 
  const <- displnorm(expo, xmin)
  
  # Compute values
  if ( ! log ) { 
    ans <- (1/const) * x^(-expo)
    ans[x < xmin] <- NaN
  } else { 
    if ( const < 0 ) { 
      # Const can be negative as nlm finds its way: the check makes sure 
      # no warning is produced by the log. 
      ans <- NaN
    } else { 
      ans <- -expo * log(x) - log(const)
    }
  }
  
  return(ans)
}

# PL: P(x>=k) 
ippl <- function(x, expo, xmin = 1) { 
  const <- displnorm(expo, xmin)
  
  is_below_xmin <- x < xmin
  
  ps <- zeta_w_xmin(expo, x[!is_below_xmin]) / const 
  
  # Values below threshold are NA'ed
  ans <- NaN*x
  ans[!is_below_xmin] <- ps
  
  return(ans)
}

# PL: Log likelihood
pl_ll <- function(dat, expo, xmin) { 
  sum( dpl(dat, expo, xmin, log = TRUE) )  
}


#' @title Distribution-fitting functions 
#'
#' @description These functions fit parametric distributions to a set of
#'   discrete values. 
#'
#' @param dat The set of values to which the distribution are fit 
#'
#' @param xmin The minimum possible value to consider when fitting the
#'   distribution
#'
#' @return A list containing at list the following components: 
#'
#' \itemize{ 
#'    \item{type: }{The type of distribution fitted (as a character string)}
#'    \item{method: }{The method used for the fit - here, maximum likelihood,
#'                  'll'}
#'    \item{ll: }{The negative log likelihood at the estimated parameter values}
#'    \item{xmin: }{The value of xmin used for the fit}
#'    \item{npars: }{The number of parameters of the distribution}
#'  }
#' 
#' Additionnaly, this list may have one or more of the following parameters 
#'   depending on the type of distribution that has been fitted: 
#'   \itemize{ 
#'     \item{plexpo: }{The exponent of the power-law}
#'     \item{cutoff: }{The rate of truncation, for truncated power law and 
#'                 exponential fits}
#'     \item{meanlog: }{The mean of the lognormal distribution}
#'     \item{sdlog: }{The s.d. of the lognormal distribution}
#'   }
#' 
#' @details These functions will fit distributions to a set of values using 
#'   maximum-likelihood estimation. In the context of the 'spatialwarnings' 
#'   package, they are most-often used to fit parametric distributions on patch
#'   size distributions. As a result, these functions assume that the data 
#'   contains only integer, strictly positive values. The type of distribution
#'   depends on the prefix of the function: 'pl' for power-law, 'tpl' for
#'   truncated power-law, 'lnorm' for lognormal and 'exp' for an exponential
#'   distribution. 
#' 
#' In the context of distribution-fitting, 'xmin' represents the minimum value 
#'   that a distribution can take. It is often used to represent the minimum 
#'   scale at which a power-law model is appropriate (Clauset et al. 2009), and 
#'   can be estimated on an empirical distribution using
#'   \code{\link{xmin_estim}}. Again, please note that the fitting procedure 
#'   assumes here that xmin is equal or grater than one.
#' 
#' Please note that a best effort is made to have the fit converge, but 
#'   it may sometimes fail when the parameters are far from their usual 
#'   range, and numerical issues may occur. It is good practice to make 
#'   sure the fits are sensible when convergence warnings are reported.
#' 
#' For reference, the shape of the distributions is as follow: 
#' 
#' \itemize{
#'   \item{power-law }{\eqn{x^{-a}}{x^(-a)} where a is the power-law exponent}
#'   \item{exponential }{\eqn{exp(-bx)}{exp(-bx)} where b is the truncation rate
#'           of the exponential } 
#'   \item{truncated power-law }{\eqn{x^{-a}exp(-bx)}{x^(-a)exp(-bx)} where a
#'     and b are the exponent of the power law and the rate of truncation}
#' }
#' 
#' The lognormal form follows the \link[=dlnorm]{standard definition}.
#' 
#' The following global options can be used to change the behavior of fitting functions 
#'   and/or produce more verbose output: 
#' \itemize{ 
#'   \item{spatialwarnings.constants.reltol}{the relative tolerance to use to compute 
#'     the power-law normalizing constant
#'     \deqn{sum_{k=1}^{\infty} x^{ak}e^{-bk}}{sum( x^(ak)exp(-bk)) for k in 1:Inf}. 
#'     Increase to increase the precision of this constant, which can be useful in some 
#'     cases, typically with large samples. Default is 1e-8.}
#'   \item{spatialwarnings.constants.maxit}{the maximum number of iterations to compute 
#'     the truncated power-law normalizing constant. Increase if you get a warning 
#'     that the relative tolerance level (defined above) was not reach. Default is 1e8}
#'   \item{spatialwarnings.debug.fit_warn_on_bound}{logical value. Warn if the fit is 
#'     at the boundary of the valid range for distribution parameter}
#'   \item{spatialwarnings.debug.fit_warn_on_NA}{logical value. Warn if the returned fit 
#'     has \code{NA}/\code{NaN} parameters}
#' }
#' 
#' @seealso \code{\link{patchdistr_sews}}, \code{\link{xmin_estim}}
#' 
#' @references
#' 
#' Clauset, Aaron, Cosma Rohilla Shalizi, and M. E. J. Newman. 2009. “Power-Law
#' Distributions in Empirical Data.” SIAM Review 51 (4): 661–703. 
#' https://doi.org/10.1137/070710111.
#' 
#' @examples 
#' 
#' # Fit an exponential model to patch size distribution 
#' exp_fit(patchsizes(forestgap[[8]]))
#'  
#' # Use the estimated parameters as an indicator function
#' \dontrun{
#' 
#' get_truncation <- function(mat) { 
#'    c(exp_cutoff = exp_fit(patchsizes(mat))$cutoff)
#' }
#' trunc_indic <- compute_indicator(forestgap, get_truncation)
#' plot(trunc_indic)
#' plot(indictest(trunc_indic, nulln = 19))
#' 
#' }
#' 
#'@export
pl_fit <- function(dat, xmin = 1) { 
  
  # Cut data to specified range
  dat <- dat[dat >= xmin] 
  
  # Start with the approximation given in Clauset's 
  npts <- length(dat)
  expo_estim <- 1 + npts / (sum(log(dat)) - npts*log(xmin-.5))
  
  negll <- function(expo) {
    result <- - pl_ll(dat, expo, xmin) 
    if ( is.infinite(result) ) { 
      return(NaN)
    } else { 
      return(result)
    }
  }
  
  est <- optim_safe(negll, expo_estim, 
                    lower = PLMIN, upper = PLMAX)
  
  result <- list(type = 'pl',
                 method = 'll', 
                 plexpo = est[["par"]],
                 ll = - est[['value']],
                 xmin = xmin,
                 npars = 1)
  
  if ( warnNA() && is.na(result[["plexpo"]]) ) { 
    warning("Fitting of PL returned NA")
    return(result)
  }
  
  if ( warnbound() && any(abs(result[["plexpo"]] - c(PLMIN, PLMAX)) < 1e-8) ) { 
    warning("Estimated PL exponent is on the boundary of the valid range")
  }
    
  return(result)
}

#' @title Estimate the minimum patch size of a power-law distribution 
#' 
#' @description When fitting a power-law to a discrete distribution, it might 
#'   be worth discarding points below a certain threshold (xmin) to improve 
#'   the fit. This function estimates the optimal xmin based on the 
#'   Kolmogorov-Smirnoff distance between the fit and the empirical 
#'   distribution, as suggested by Clauset et al. (2009). 
#' 
#' @param dat A vector of integer values
#' 
#' @param bounds A vector of two values representing the bounds in which 
#'   the best xmin is searched
#' 
#' @return The estimated xmin as an integer value 
#' 
#' @details The function returns NA if \code{dat} has only three unique values 
#'   or if the power-law fit failed. 
#' 
#' @seealso \code{\link{patchdistr_sews}}
#' 
#' @references 
#' 
#' Clauset, A., Shalizi, C. R., & Newman, M. E. (2009). 
#'   Power-law distributions in empirical data. SIAM review, 51(4), 661-703.
#' 
#' @examples 
#' 
#' \dontrun{ 
#' psd <- patchsizes(forestgap[[5]])
#' xmin_estim(psd)
#' }
#'@export
xmin_estim <- function(dat, bounds = range(dat)) { 
  
  # Create a vector of possible values for xmin
  xmins <- sort(unique(dat))
  
  # We need at least 3 values for a pl fit, so the last value of xmin 
  # needs to have three points after it
  if ( length(xmins) <= 3 ) { 
    return(NaN)
  }
  
  # We build a vector of possible xmins. The last three values are stripped 
  #   away as they won't allow enough data for a fit
  xmins <- head(xmins, length(xmins)-3)
  xmins <- xmins[xmins >= min(bounds) & xmins <= max(bounds)]
  
  # Compute all ks-distances
  kss <- adply(xmins, 1, get_ks_dist, dat = dat)[ ,2]
  
  if ( all(is.nan(kss)) ) { 
    return(NaN)
  }
  
  # Note that sometimes the fit fails, especially when xmin is around the 
  #   distribution tail -> we need to remove some NAs here
  xmin <- xmins[!is.na(kss) & kss == min(kss, na.rm = TRUE)]
  
  # Note that xmin can be NaN
  return(xmin)
}

get_ks_dist <- function(xmin, dat) { 
  # Crop dat to values above xmin and compute cdf
  dat <- dat[dat >= xmin]
  # Compute empirical (inverse) cdf values
  udat <- unique(dat)
  cdf_empirical <- rep(NA, length(dat))
  for ( val in udat ) { 
    cdf_empirical[dat == val] <- mean(dat >= val)
  }
  
  # Fit and retrieve cdf. Note: here we suppress the warnings because finding 
  # xmin requires removing patches below a threshold, which often leads to fit 
  # being done on pathological cases like few unique patch sizes
  fit <- suppressWarnings({ 
    pl_fit(dat, xmin = xmin)
  })
  
  if ( is.na(fit[['plexpo']]) ) { 
    # Note: a warning was already produced in this case as it means that the 
    # fit failed to converge: we do not produce one here again. 
    return(NaN)
  }
  
  cdf_fitted <- ippl(dat, fit[["plexpo"]], fit[["xmin"]])

#   # debug
#   plot(data.frame(dat, rbinom(length(dat), 1, .5)), type = 'n')
#   plot(log10(data.frame(dat, cdf_empirical))) 
#   points(log10(data.frame(dat, cdf_fitted)), col = 'red')
#   browser()
#   zeta.fit(dat, xmin)$exponent
#   fit$expo
  
  # We return the ks distance
  maxks <- max(abs(cdf_empirical - cdf_fitted))
#   cat(xmin, ",", max(dat), "->", maxks,  "\n" )
  return( maxks )
}




# EXP fitting 
# ---------------------------------------
# pexp/dexp is already implemented in R

# EXP: P(x = k)
ddisexp <- function(dat, rate, xmin = 1, log = FALSE) { 
  # sum(P = k) for k = 1 to inf
  if ( log ) { 
    
    const <- log(1 - exp(-rate)) + rate * xmin
    return( ifelse(dat < xmin, NaN, const - rate * dat) )
  
  } else { 
    
    const <- (1 - exp(-rate)) * exp(rate*xmin) 
    return( ifelse(dat < xmin, NaN, const * exp(-rate * dat)) )
    
  }
} 

# EXP: P(x>=k)
# Imported and cleaned up from powerRlaw (def_disexp.R)
ipdisexp <- function(x, rate, xmin) {
  # p >= k
  p <- pexp(x + .5, rate, lower.tail = FALSE) 
  # p >= 1 
  const <-  1 - pexp(xmin + .5, rate) 
  return(p/const)
}

exp_ll <- function(dat, rate, xmin) { 
  sum( ddisexp(dat, rate, xmin, log = TRUE)) 
}

#'@rdname pl_fit
#'@export
exp_fit <- function(dat, xmin = 1) { 
  
  dat <- dat[dat>=xmin]
  
  rate0 <- 1 / mean(dat)
  
  negll <- function(rate) {
    - exp_ll(dat, rate, xmin)
  }
  
  est <- optim_safe(negll, rate0, 
                    fit_on_logscale = TRUE, 
                    lower = EXPMIN, 
                    upper = EXPMAX)
  
  result <- list(type = 'exp',
                 method = 'll', 
                 cutoff = est[['par']], 
                 ll = - est[["value"]],
                 npars = 1)
  
  if ( warnbound() && any(abs(result[["cutoff"]] - c(EXPMIN, EXPMAX)) < 1e-8) ) { 
    warning("Estimated EXP exponent is on the boundary of the valid range")
  }
  
  return(result)
}





# LNORM fitting 
# ---------------------------------------

# LNORM: P(X=k)
ddislnorm <- function(x, meanlog, sdlog, xmin, log = FALSE) { 
  
  p_over_thresh <- plnorm(xmin - .5, meanlog, sdlog, lower.tail = FALSE)
   
  p_equals_k <- plnorm(x-.5, meanlog, sdlog, lower.tail = FALSE) - 
                  plnorm(x+.5, meanlog, sdlog, lower.tail = FALSE)
  
  if ( !log ) { 
    return( ifelse(x<xmin, NaN, p_equals_k / p_over_thresh) )
  } else { 
    return( ifelse(x<xmin, NaN, log(p_equals_k) - log(p_over_thresh)) )
  }
}

# LNORM: P(X>=k)
ipdislnorm <- function(x, meanlog, sdlog, xmin) { 
  px_supto_k <- plnorm(x - .5, meanlog, sdlog, lower.tail = FALSE)
  px_supto_xmin <- plnorm(xmin - .5, meanlog, sdlog, lower.tail = FALSE)
  ifelse(x<xmin, NaN, px_supto_k / px_supto_xmin)
}

# LNORM: LL
lnorm_ll <- function(x, meanlog, sdlog, xmin) { 
  x <- x[x>=xmin]
  sum( ddislnorm(x, meanlog, sdlog, xmin, log = TRUE) ) 
}

# LNORM: fit
#'@rdname pl_fit
#'@export
lnorm_fit <- function(dat, xmin = 1) { 
  
  # Pars[1] holds mean of log-transformed data
  # Pars[2] holds sd 
  pars0 <- c( mean(log(dat)), sd(log(dat)) )
  
  negll <- function(pars) { 
    ll <- - lnorm_ll(dat, pars[1], pars[2], xmin) 
    if ( is.finite(ll) ) ll else 1e10
  }
  
  est <- optim_safe(negll, pars0)
  
  result <- list(type = 'lnorm',
                 method = 'll', 
                 meanlog = est[['par']][1], 
                 sdlog = est[['par']][2], 
                 ll = - est[["value"]], 
                 npars = 2)
  return(result)
}




# TPL fitting 
# ---------------------------------------

tplnorm <- function(expo, rate, xmin) { 
  
  maxit <- getOption("spatialwarnings.constants.maxit", default = 1e8L)
  reltol <- getOption("spatialwarnings.constants.reltol", default = 1e-8)
  
  a <- tplinfsum(expo, rate, xmin, maxit, reltol)
  
  a
}

# P(x=k)
dtpl <- function(x, expo, rate, xmin, log = FALSE) { 
  
  const <- tplnorm(expo, rate, xmin)
  if ( ! log ) { 
    ps <- x^(-expo) * exp(- x * rate) / const
  } else { 
    ps <- - expo * log(x) - rate * x - log(const)
  }
  
  return( ifelse(x < xmin, NaN, ps) )
}

# P(x>=k)
iptpl <- function(x, expo, rate, xmin) { 
  const <- tplnorm(expo, rate, xmin)
  
  # tplsum is vectorized over x
  p_inf_to_k <- tplsum(expo, rate, x, xmin) / const
  
  return( 1 - p_inf_to_k ) 
}

tpl_ll <- function(x, expo, rate, xmin, approximate = FALSE) { 
  x <- x[x>=xmin]
  
  ll <- sum( dtpl(x, expo, rate, xmin, log = TRUE) )
  
  if ( !is.finite(ll) ) { 
    ll <- sign(ll) * .Machine$double.xmax
  }
  return( ll )
} 

#'@rdname pl_fit
#'@export
tpl_fit <- function(dat, xmin = 1) { 
  
  negll <- function(pars) { 
    - tpl_ll(dat, pars[1], pars[2], xmin) * 1000
  }
  
  # Initialize and find minimum
  expo0 <- pl_fit(dat, xmin)[['plexpo']] 
  
  # Do a line search over the cutoff to find a minimum, starting from zero 
  # up to 100
  is <- seq(0, 100, length = 128)
  is <- 10^seq(-7, 2, l = 128)
  lls <- unlist(lapply(is, function(i) { 
    negll(c(expo0, i))
  }))
  llmin <- min(lls[abs(lls) != .Machine$double.xmax])
  expmrate0 <- is[which(lls == llmin)]
  
  # Debug thing lls
#   plot( log(spatialwarnings:::cumpsd(dat)) )
#   vals <- ippl(dat, expo0, xmin = xmin)
#   points(log(dat), log(vals), col = "red")
#   vals <- iptpl(dat, expo0, rate = 0, xmin = xmin)
#   points(log(dat), log(vals), col = "darkgreen")
#   vals <- iptpl(dat, expo0, rate = expmrate0, xmin = xmin)
#   points(log(dat), log(vals), col = "blue")
  
  pars0 <- c(expo0, expmrate0)
  est <- optim_safe(negll, pars0, 
                    lower = c(TPL_EXPOMIN, TPL_RATEMIN), 
                    upper = c(TPL_EXPOMAX, TPL_RATEMAX), 
                    fit_on_logscale = FALSE)
  
  # For very small cutoffs, we may want to fit on log scale to get a descent estimate, 
  # as otherwise there is too little variation in the cutoff for BFGS
  if ( any(is.na(est[["par"]])) ) { 
    est <- optim_safe(negll, pars0, 
                      lower = c(TPL_EXPOMIN, TPL_RATEMIN), 
                      upper = c(TPL_EXPOMAX, TPL_RATEMAX), 
                      fit_on_logscale = TRUE)
  }
  
  result <- list(type = 'tpl',
                 method = 'll', 
                 plexpo = est[['par']][1], 
                 cutoff = est[['par']][2], 
                 ll = - est[["value"]],
                 npars = 2)
  
  if ( warnNA() && is.na(result[["plexpo"]]) ) { 
    warning("Fitting of TPL returned NA")
  }
  
  if ( is.na(result[["plexpo"]]) ) { 
    return(result)
  }
  
  if ( warnbound() && 
       any(abs(result[["plexpo"]] - c(TPL_EXPOMIN, TPL_EXPOMAX)) < 1e-8) ) { 
    warning("Estimated TPL exponent is on the boundary of the valid range")
  }
  
  # We do not warn for reaching minimum exponent because this is quite common when 
  # fitting truncated power laws to have a very small exponent
  if ( warnbound() && 
       any(abs(result[["cutoff"]] - c(TPL_RATEMIN, TPL_RATEMAX)) < 1e-8) ) { 
    warning("Estimated TPL cutoff is on the boundary of the valid range")
  }
  
  return(result)
}
