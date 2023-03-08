# The functions below fit the Lifshitz, Slyozov, Wagner (LSW) distribution to 
# artificially generated patch sizes.
# Patch size distributions are expected to approach the LSW distribution as cover goes 
# to zero (Siteur et al., 2023).
# 
# Note that the patch sizes represent the patch radii, and not the patch area.
# Patch radius can be estimated by assuming circular patches r=sqrt(a/2pi).
# 
# LSW does not fit well to distributions with a positive skewness, as patch radii 
# greater than 1.5 times the mean patch size have zero probability.
# 

lsw_sews <- function(mat, wrap = FALSE) { 
  compute_indicator(mat, raw_lsw_sews, 
                    wrap = wrap, 
                    taskname = "LSW-based indicators")
}

#'@export
raw_lsw_sews <- function(mat, wrap) {
  # Convert object to matrix form (or list of matrices) and check if it is
  # suitable for spatialwarnings
  mat <- convert_to_matrix(mat)
  if ( is.list(mat) ) { 
    lapply(mat, check_mat)
  } else { 
    check_mat(mat)
  }
  
  # Get patch size distribution. Assume circle patches from the area. 
  psd <- patchsizes(mat, wrap = wrap)
  radii <- sqrt(psd / pi)
  
  # Compute skewness and warn if positive
  skewness <- cpp_skewness(radii)
  if ( skewness > 0 ) { 
    warning("Positive skewness in radius distribution, LSW fit may fail")
  }
  
  # Fit the LSW distribution 
  fit_lsw <- fitLSW(radii)
  
  # Fit a lognormal
  fit_lnorm <- fit_lnorm(radii)
  
  # For now we just return the skewness and difference in AIC. We will need something 
  # a bit more fancy to have distributions fitting et al. 
  c(cover = mean(mat), 
    skewness = skewness, 
    aic_diff = fit_lsw[["AIC"]] - fit_lnorm[["AIC"]])
}

fit_lnorm <- function(xs) { 
  
  negll <- function(theta) { 
    - sum( dlnorm(xs, mean = theta[1], sd = theta[2]) )
  }
  
  # TODO: check that this makes sense and I haven't been confused by the logs, also 
  # there is probably an explicit formula for the mle of a lognormal 
  # (mean of log / sd of log? )
  mle <- optim(par = c(mean(log(xs)), 
                       sd(log(xs))), 
               fn = negll) 
  theta_est <- mle[["par"]]
  
  aic <- 2 - ( - negll(theta_est) ) 
  list(meanlog = theta_est[1], 
       sdlog   = theta_est[2], 
       AIC     = aic)
}

# The Lifshitz, Slyozov, Wagner (LSW) probability density function. Used for plotting, 
# and for computation of the AIC.
dLSW <- function(r, mu){
  ifelse(r<1.5*mu,
         (4*r^2/(9*mu^3)) * (3*mu/(3*mu+r))^(7/3) * 
          (3*mu/(3*mu-2*r))^(11/3) * exp(2*r/(2*r-3*mu)),
         0)
}

# The (approximate) cumulative LSW distribution. Used for plotting.
cLSW <- function(r,mu){
  x <- seq(0, 1.5, 0.01) * mu
  PDF <- dLSW(x, mu)
  CDF <- cumsum(PDF) / sum(PDF)
  return(approx(x, CDF, r)$y)
}

# The (approximate) inverse cumulative LSW distribution. Used for drawing random 
# patch radii
icLSW <- function(X,mu){
  x <- seq(0,1.5,0.001)*mu
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
#   mu_fit <- uniroot(function(mu){dLdmu(r,mu)},
#                     interval=c(min(r),
#                                max(r)))$root
  mu_fit <- uniroot(function(mu){dLdmu(r,mu)},
                    interval=c(0, max(r)))$root
  AIC_fit <- 2-2*sum(log(dLSW(r,mu_fit)))
  return(list(mu = mu_fit, 
              AIC = AIC_fit))
}



if ( FALSE ) { 
  # generate 100 patch radii with a mu of 1
  patch_radii <- icLSW(runif(100),mu = 1)

  # fit the LSW distribution
  ft <- fitLSW(patch_radii)
  
  # plot of histogram and PDF
  hist(patch_radii,probability = T,xlab = "Patch radius (m)",main="Histogram of patch radii")
  x <- seq(0,1.5,0.02)*ft$mu
  lines(x,dLSW(x,ft$mu),col="red")
  grid()
  
  # plot of cumulative probability distribution and CDF 
  PXltx <- rank(patch_radii,ties.method = "average")/(length(patch_radii)+1)
  plot(patch_radii,PXltx,pch=20,xlab = "Patch radius (m)",ylab="Cumulative probability",main="Cumulative probability distribution of patch radii")
  lines(x,cLSW(x,mu = ft$mu),col="red")
  grid()

  # for some reason 1-CDF is often plotted and a log-log scale is used
  PXltx <- rank(patch_radii,ties.method = "average")/(length(patch_radii)+1)
  plot(patch_radii,1-PXltx,pch=20,log="xy",xlab = "Patch radius (m)",ylab="Cumulative probability",main="Cumulative probability distribution of patch radii")
  lines(x,1-cLSW(x,mu = ft$mu),col="red")
  grid()
}

