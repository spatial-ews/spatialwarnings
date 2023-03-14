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
#'@param mat A logical matrix (TRUE/FALSE values) or a list of these
#' 
#'@param wrap Determines whether patches are considered to wrap around the 
#'  matrix when reaching on of its edges
#' 
#'@details 
#' 
#' # TODO
#' 
#' 
#' @return 
#' 
#' \code{lsw_sews} returns an object of class \code{simple_sews_single}
#'   (actually a list) if \code{mat} is a single matrix or an object of class 
#'   \code{simple_sews_list} if \code{mat} is a list. You probably want to use some  
#'   of the methods written for these complicated objects instead of extracting 
#'   values directly (they are displayed using \code{print(<object>)}).
#'
#'@seealso dLSW, dda
#'
#'@references 
#'  
#'  # TODO Siteur 2023, other relevant papers ? 
#'
#'@examples 
#'
#' data(dda)
#' data(dda.pars)
#' 
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
#' 
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
  
  # Compute skewness and warn if positive
  radii <- get_patch_radii(mat, wrap)
  skewness <- cpp_skewness(radii)
  aicw <- get_lsw_aicw(radii)
  
  # For now we just return the skewness and difference in AIC. We will need something 
  # a bit more fancy to have distributions fitting et al. 
  c(cover = mean(mat), 
    skewness = skewness, 
    lsw_aicw = aicw)
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
  return( c(lsw_aicw = aicw) )
}

get_patch_radii <- function(mat, wrap = FALSE) { 
  # Get patch size distribution. Assume circle patches from the area. 
  psd <- patchsizes(mat, wrap = wrap)
  radii <- sqrt(psd / pi)
}

get_lsw_aicw <- function(radii) {
  fit_lsw <- fitLSW(radii)
  
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
#'@param vector of means 
#' 
#'@param log, log.p logical; if \code{TRUE}, probabilities p are given as log(p)
#'
#'@param lower.tail logical; if TRUE (default), probabilities are
#'         \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}
#' 
#'@details 
#' # TODO
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
#'@export
dLSW <- function(x, mu, log = FALSE){
  p <- ifelse(x<1.5*mu,
              (4*x^2/(9*mu^3)) * (3*mu/(3*mu+x))^(7/3) * 
                (3*mu/(3*mu-2*x))^(11/3) * exp(2*x/(2*x-3*mu)),
              0) 
  if ( log ) { 
    # TODO: actually implement this by working out the math on paper for 
    # better numerical stability, but right now it works well
    p <- log(p)
  }
  
  p
}

# The LSW distribution, but standard 
#'@rdname dLSW 
#'
#'@export 
pLSW <- function(r, mu, lower.tail = TRUE, log.p = FALSE) { 
  stop("TODO")
}

# The (approximate) cumulative LSW distribution. Used for plotting.
cLSW <- function(r,mu) {
  # TODO: is there an analytical expression for the LSW distribution ? 
  x <- seq(0, 1.5, 0.01) * mu
  PDF <- dLSW(x, mu)
  CDF <- cumsum(PDF) / sum(PDF)
  return(approx(x, CDF, r)$y)
}

# The (approximate) inverse cumulative LSW distribution. Used for drawing random 
# patch radii
icLSW <- function(X, mu){
  x <- seq(0,1.5,0.001) * mu
  PDF <- dLSW(x,mu)
  CDF <- cumsum(PDF)/sum(PDF)
  return(suppressWarnings((approx(CDF,x,X)$y)))  #I suppress warnings here. THe warning is because the approx function is fed Y values with the same X value. This is because the LSW distribution becomes very flat as it approaches 1.5.
}

# The derivative of the log likelihood function to parameter mu. Used to estimate mu
dLdmu <- function(r, mu){
  n <- length(r)
  -3*n+(7/3)*sum(r/(3*mu+r))-(11/3)*sum(2*r/(3*mu-2*r))+3*mu*sum(2*r/(2*r-3*mu)^2)
}

# Function to fit the LSW probability density function. Returns the estimate of mu and 
# the AIC for model comparison
fitLSW <- function(r){
  
  negll <- function(mu) { 
    nll <- dLSW(r, mu, log = TRUE) 
    if ( any( ! is.finite(nll) ) ) { 
      return(Inf)
    }
    - sum(nll)
  }
  
  # Here we want to stick to values so that all r < 3mu/2 otherwise the log 
  # likelihood is infinite, which means max(r) < 3mu/2, 
  # i.e  mu > 2*max(r)/3. We find the smallest value that produces a finite negll and 
  # use that as our lower bound for LL minimization later on. 
  boundsearch <- uniroot(f = function(x) is.finite(negll(x)) - 0.5, 
                         lower = min(r), 
                         upper = max(r))
  minbound <- boundsearch[["root"]] + boundsearch[["estim.prec"]]
  
  est <- optim(par = max(r), 
               lower = minbound, 
               upper = max(r), 
               method = "L-BFGS-B", 
               fn = negll)
  
  # Show dldmu 
#   mus <- seq(0, max(r), l = 64)
#   mus <- seq(2*max(r)/3 + epsilon, max(r), l = 64)
#   plot(mus, sapply(mus, function(mu) { 
#     dLdmu(r, mu)
#     negll(mu)
#   }))
#   
#   mu_fit <- uniroot(function(mu){ dLdmu(r, mu) },
#                     interval = c(minbound, max(r)))$root
#   AIC_fit <- 2 - 2 * sum(log(dLSW(r, mu_fit)))
#   cat(sprintf("directll mu est: %s / uniroot mu est: %s / ll directll: %s\n", 
#               mu_fit, est[["par"]], -est[["value"]]))
  
  mu_fit <- est[["par"]]
  AIC_fit <- 2 - 2 * ( - est[["value"]] )
  
  return(list(mu = mu_fit, 
              AIC = AIC_fit))
}



