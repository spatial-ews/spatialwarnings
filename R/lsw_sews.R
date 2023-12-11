# 
# The functions below fit the Lifshitz, Slyozov, Wagner (LSW) distribution to 
# artificially generated patch sizes. Patch size distributions are expected to 
# approach the LSW distribution as cover goes to zero (Siteur et al., 2023).
# 
# Note that the patch sizes represent the patch radii, and not the patch area.
# Patch radius can be estimated by assuming circular patches r=sqrt(a/2pi).
# 
# LSW does not fit well to distributions with a positive skewness, as patch radii 
# greater than 1.5 times the mean patch size have zero probability.
# 
# 
# 
# 
# 
#'@title LSW-based indicators  
#'
#'@description LSW indicators for systems with density-dependent aggregation 
#'
#'
#'@param mat A logical matrix (TRUE/FALSE values) or a list of such matrices
#' 
#'@param wrap Determines whether patches are considered to wrap around the 
#'  matrix when reaching on of its edges
#' 
#'@details 
#' 
#' In systems where a mobile resource or consumer can be fixed in space by a sessile 
#' species, specific patterns are expected to appear. Such systems can include situations 
#' where nutrients available in the environmeent are fixed by a sessile species (e.g. 
#' seagrasses or corals), or where the behavior of herbivores is altered to restrict 
#' them to certain areas (see full theoretical background in Siteur et al. 2023). 
#' 
#' In those systems, as environmental conditions change and the global density of 
#' the sessile species decreases, its spatial structure is expected to change. The area 
#' of patches of the sessile species (as measured by their radii, which assumes circular
#' patches), is expected to go from a log-normal to a  Lifshitz–Slyozov–Wagner (LSW)
#' distribution. Thus, measuring how close the observed distribution of radii are to 
#' those two candidate distributions can constitute an indicator of ecosystem 
#' degradation. 
#' 
#' This function measures this through the relative support based on AIC for the two 
#' distributions (equal to 1 when the empirical distribution is best-approximated by an 
#' LSW, and 0 when it is a \link{}{log-normal distribution}), and the skewness of the
#' observed patch radii, which should approach a value around -0.92 as conditions 
#' worsen.  
#' 
#' @author This code has received contributions from Koen Siteur (\email{k.siteur@gmail.com})
#' 
#' @return 
#' 
#' \code{lsw_sews} returns an object of class \code{simple_sews_single}
#'   (a list) if \code{mat} is a single matrix or an object of class 
#'   \code{simple_sews_list} if \code{mat} is a list. You probably want to use some  
#'   of the methods written for these complicated objects instead of extracting 
#'   values directly (they are displayed using \code{print(<object>)}). 
#' 
#'@seealso dLSW, dda, raw_patch_radii_skewness, raw_lsw_aicw
#'
#'@references 
#'  
#' Siteur, Koen, Quan-Xing Liu, Vivi Rottschäfer, Tjisse van der Heide, Max Rietkerk, 
#' Arjen Doelman, Christoffer Boström, and Johan van de Koppel. 2023. 
#' "Phase-Separation Physics Underlies New Theory for the Resilience of Patchy 
#' Ecosystems." Proceedings of the National Academy of Sciences 120 (2): e2202683120.
#' https://doi.org/10.1073/pnas.2202683120.
#' 
#  # TODO other relevant papers ? 
#'
#'@examples 
#'
#' data(dda)
#' data(dda.pars)
#' 
#' # Compute all indicators at once (skewness and relative AIC support)
#' indics <- lsw_sews(dda)
#' plot(indics, along = dda.pars[ ,"tau"]) 
#' 
#' # Compute individual indicators 
#' 
#' # Skewness of the distribution of patch radii
#' radii_skewness <- compute_indicator(dda, raw_patch_radii_skewness)
#' plot(radii_skewness, along = dda.pars[ ,"tau"])
#' 
#' # Aic weight of LSW distribution relative to a lognormal distribution 
#' lsw_aicw <- compute_indicator(dda, raw_lsw_aicw)
#' plot(lsw_aicw, along = dda.pars[ ,"tau"])
#' 
#'@export
lsw_sews <- function(mat, wrap = FALSE) { 
  compute_indicator(mat, internal_lsw_sews, 
                    wrap = wrap, 
                    taskname = "LSW-based indicators")
}

internal_lsw_sews <- function(mat, wrap) {
  # Convert object to matrix form (or list of matrices) and check if it is
  # suitable for spatialwarnings
  mat <- convert_to_matrix(mat)
  if ( is.list(mat) ) { 
    lapply(mat, check_mat)
  } else { 
    check_mat(mat)
  }
  
  # Compute metrics
  radii <- get_patch_radii(mat, wrap)
  skewness <- cpp_skewness(radii)
  aicw <- get_lsw_aicw(radii)
  
  # For now we just return the numerical values. We will need something a bit more 
  # fancy to have distributions fitting et al. 
  c(cover = mean(mat), 
    skewness = skewness, 
    lsw_prob = aicw)
}

#'@rdname lsw_sews
#'@export
raw_patch_radii_skewness <- function(mat, wrap = FALSE) { 
  radii <- get_patch_radii(mat, wrap)
  c(psd_radii_skewnes = cpp_skewness(radii))
}

#'@rdname lsw_sews
#'@export
raw_lsw_aicw <- function(mat, wrap = FALSE) { 
  # Get patch size distribution. Assume circle patches from the area. 
  radii <- get_patch_radii(mat, wrap)
  aicw <- get_lsw_aicw(radii)
  return( c(lsw_prob = aicw) )
}

get_patch_radii <- function(mat, wrap = FALSE) { 
  # Get patch size distribution. Assume circle patches from the area. 
  psd <- patchsizes(mat, wrap = wrap)
  radii <- sqrt(psd / pi)
}

get_lsw_aicw <- function(radii) {
  fit_lsw <- LSW_fit(radii)
  
  # Fit a lognormal
  fit_lnorm <- fit_lnorm(radii)
  
  # Compute AIC weight of LSW relative to lnorm 
  minaic <- min(fit_lsw[["AIC"]], fit_lnorm[["AIC"]])
  weight_lsw_min   <- exp( - 0.5 * ( fit_lsw[["AIC"]]   - minaic ) ) 
  weight_lnorm_min <- exp( - 0.5 * ( fit_lnorm[["AIC"]] - minaic ) ) 
  
  aicw <- weight_lsw_min / ( weight_lsw_min + weight_lnorm_min )
  
  return(aicw)
}

fit_lnorm <- function(xs) { 
  
  negll <- function(theta) { 
    - sum( dlnorm(xs, meanlog = theta[1], sdlog = theta[2], log = TRUE) )
  }
  
  # TODO: check that this makes sense and I haven't been confused by the logs, also 
  # this is probably useless (MLE = mean of log / sd of log? )
  mle <- optim(par = c(mean(log(xs)), 
                       sd(log(xs))), 
               fn = negll) 
  theta_est <- mle[["par"]]
  
  aic <- 2 - 2 * ( - negll(theta_est) ) 
  list(meanlog = theta_est[1], 
       sdlog   = theta_est[2], 
       AIC     = aic)
}

# The Lifshitz, Slyozov, Wagner (LSW) probability density function. Used for plotting, 
# and for computation of the AIC.
# 
#'@title The Lifshitz-Slyozov-Wagner distribution 
#'
#'@description Density and distribution function for the Lifshitz-Slyozov-Wagner (LSW)
#'  distribution with mean mu. 
#'
#'@param x vector of quantiles 
#'
#'@param mu the mean of the distribution
#'
#'@param log, log.p logical; if \code{TRUE}, probabilities p are given as log(p)
#'
#'@param lower.tail logical; if TRUE (default), probabilities are
#'         \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}
#' 
#'@details 
#' 
#' The LSW distribution is a continuous distribution with density
#' 
#' \deqn{ 
#'   f(x, \mu) = \frac{4x^2}{9\mu^3} ( \frac{3\mu}{3\mu + x} )^{7/3} 
#'                 ( \frac{3\mu}{3\mu - 2x} )^(11/3) e^{\frac{2x}{2x - 3\mu}} 
#' }
#' 
#' where \eqn{\mu} is the mean of the distribution. 
#' 
#' The functions \code{dLSW} gives the probability density, \code{pLSW} gives the 
#' distribution function. \code{qLSW} and \code{rLSW} are not implemented. You 
#' can use \code{LSW_fit} to fit an LSW distribution to a set of observations. 
#' 
#' The length of the results is determined by the length of \code{x}, and \eqn{\mu} can 
#' only be a single value. 
#' 
#' Please note that this distribution has support only on the interval [0,3*mu) 
#'   (TODO: does it??). Probabilities outside this interval are returned as 0. 
#'
#'@references 
#
#TODO: add reference
#
#' Siteur, Koen, Quan-Xing Liu, Vivi Rottschäfer, Tjisse van der Heide, Max Rietkerk, 
#' Arjen Doelman, Christoffer Boström, and Johan van de Koppel. 2023. 
#' "Phase-Separation Physics Underlies New Theory for the Resilience of Patchy 
#' Ecosystems." Proceedings of the National Academy of Sciences 120 (2): e2202683120.
#' https://doi.org/10.1073/pnas.2202683120.
#'
#'@return
#' dLSW gives the density, pLSW gives the distribution function, both as numerical 
#'   vectors determined by the length of x. 
#'
#'@examples
#' 
#' # Plot the density 
#' x <- seq(0, 10, l = 128) 
#' plot(x, dLSW(x, mu = 3), type = "l", col = "black")
#' lines(x, dLSW(x, mu = 5), type = "l", col = "red")
#' lines(x, dLSW(x, mu = 7), type = "l", col = "blue")
#' legend(x = 0, y = max(dLSW(x, mu = 3)), lty = 1, col = c("black", "red", "blue"), 
#'        legend = paste("mu =", c(3, 5, 7)))
#' 
#TODO: what happens to values below 0? We return prob = 0, but is this the right 
#definition?
#'@export
dLSW <- function(x, mu, log = FALSE) { 
  if ( length(mu) > 1 ) { 
    stop("dLSW does not support multiple values for mu")
  }
  
  zero_val <- ifelse(log, -Inf, 0)
  p <- rep(zero_val, length(x))
  x_ok <- which(x<1.5 * mu & x >= 0)
  x <- x[x_ok]
  
  if ( log ) { 
    p[x_ok] <- - 7 * log(x+3*mu) /3+2 * log(x) + 
                  2 * x / (2 * x -  3 * mu) - 
                  ( 11 * log(3 * mu - 2 * x) ) / 3 + 
                  3 * log(mu) + 
                  log(324)
  } else { 
    p[x_ok] <- (4*x^2/(9*mu^3)) * 
                 (3*mu/(3*mu+x))^(7/3) * 
                 (3*mu/(3*mu-2*x))^(11/3) * 
                 exp(2*x/(2*x-3*mu))
  }
  
  p
}

# The LSW distribution, but standard interface
#'@rdname dLSW 
#'
#'@export 
pLSW <- function(x, mu, lower.tail = TRUE, log.p = FALSE) { 
  
  # TODO: is there an analytical expression for the LSW distribution function ? 
  ans <- vapply(x, function(thisr) { 
    integrate(dLSW, lower = 0, upper = thisr, mu = mu, log = FALSE, 
              subdivisions = 512)[["value"]]
  }, numeric(1))
  
  if ( ! lower.tail ) { 
    ans <- 1 - ans
  }
  
  if ( log.p ) { 
    ans <- log(ans)
  }
  
  ans
}
# 
# # The derivative of the log likelihood function to parameter mu. Used to estimate mu
# dLdmu <- function(r, mu){
#   n <- length(r)
#   -3*n+(7/3)*sum(r/(3*mu+r))-(11/3)*sum(2*r/(3*mu-2*r))+3*mu*sum(2*r/(2*r-3*mu)^2)
# }

# Function to fit the LSW probability density function. Returns the estimate of mu and 
# the AIC for model comparison
#'@export
LSW_fit <- function(x) {
  
  negll <- function(mu) { 
    nll <- - dLSW(x, mu, log = TRUE) 
    if ( any( ! is.finite(nll) ) ) { 
      return(Inf)
    }
    sum(nll)
  }
  
  # Here we want to stick to values so that all r < 3mu/2 otherwise the log 
  # likelihood is infinite, which means max(r) < 3mu/2, 
  # i.e  mu > 2*max(r)/3. We find the smallest value that produces a finite negll and 
  # use that as our lower bound for LL minimization later on. 
  boundsearch <- uniroot(f = function(x) is.finite(negll(x)) - 0.5, 
                         lower = min(x), 
                         upper = max(x))
  minbound <- boundsearch[["root"]] + boundsearch[["estim.prec"]]
  
  est <- optim(par = max(x), 
               lower = minbound, 
               upper = max(x), 
               method = "L-BFGS-B", 
               fn = negll)
  
  mu_fit <- est[["par"]]
  AIC_fit <- 2 - 2 * ( - est[["value"]] )
  
  return(list(mu = mu_fit, 
              AIC = AIC_fit))
}



