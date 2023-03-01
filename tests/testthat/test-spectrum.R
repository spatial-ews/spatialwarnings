# 
# 
# Test the computation of spectra (r-spectrum so far)
# 

context("Test the correct computation of r-spectrum")

test_that("the cpp implementation of the spectrum computations is correct", { 
  
  # Redefine the old functions
  rspectrum_old <- function(mat, step) { 

  nr <- nrow(mat)
  nc <- ncol(mat)
  
  n0x <- floor(nc/2) + 1
  n0y <- floor(nr/2) + 1
  
  # Create distance and angle matrices
  f1 <- t(replicate(nr,seq(1,nc))) - n0x
  f2 <- replicate(nc,seq(1,nr)) - n0y
  DIST <- sqrt(f1^2 + f2^2)
  
  # Calculate DFT
  mi <- 1
  ma <- min(c(n0x,n0y))
  DISTMASK <- DIST>=mi & DIST <= ma
  
  tmp <- fft(mat)
  class(tmp) <- "matrix"
  tmpshift <- myfftshift(tmp)
  tmpshift[n0x, n0y] <- 0
  aspectr2D <- abs(tmpshift)^2 / (n0x*n0y)^4
  
  sig2 <- sum(aspectr2D[DISTMASK]) #Normalisation
  
  aspectr2D <- aspectr2D/sig2 #Normalisation
  
  # Now calculate r-spectrum
  STEP <- 1
  ray <- seq(mi,ma,STEP)
  rspectr <- numeric(length(ray))
  for (i in 1:length(ray))
  {
    m <- DIST >= ray[i] - STEP/2 & DIST < ray[i] + STEP/2
    rspectr[i] <- mean(aspectr2D[m])
  }

  out <- data.frame(dist  = ray, rspec = rspectr)
  
  return(out)
  }
  
  myfftshift <- function(X) {
    nr=dim(X)[1]
    nc=dim(X)[2]
    shiftX = X
    if (nr != nc)
      print("Not a square matrix X")
    else
    {
      n=nc
      shift = floor(n/2)
      for (i in 1:n)
        for (j in 1:n)
        {
          a = (i+shift)%%n
          b = (j+shift)%%n
          if (a==0) a = n
          if (b==0) b = n
          shiftX[a,b] = X[i,j]
        }
    }
    return(shiftX)
  }
  
  for ( i in c(1, 2, 3) ) { 
    testmat <- serengeti[[i]]
    
    # Test if there is more than one value in the matrix, because the behavior 
    # is to return NA whereas the old code just compute things 
    if ( length(unique(as.vector(testmat))) > 1 ) { 
    
    ospec <- rspectrum_old(testmat)
    expect_equal(ospec, 
                 rspectrum(testmat), 
                 tolerance = 1/1000) 
    }
  }
  
  # Test on odd-sized matrix 
  s <- 101
  modd <- matrix(runif(s^2) < 0.5, ncol = s)
  expect_equal(rspectrum(modd), rspectrum_old(modd), 
               tolerance = 1/1000)
  
})

test_that("rspectrum/sdr show realistic values", { 
  
  
  # The expected value of sdr on random matrices is one 
  sdr_random <- mean(replicate(299, { 
    n <- 100
    m <- matrix(rnorm(n^2), nrow = n, ncol = n)
    raw_sdr(m, c(0, 0.2), c(0, 1))
  }))
  
  expect_true({ 
    abs(sdr_random - 1 ) < 0.1
  })
  
  # When there is reddening, the sdr increases 
  sdr_red <- mean(replicate(299, { 
    n <- 100
    m <- matrix(rnorm(n^2), nrow = n, ncol = n)
    m[2:n, 2:n] <- (m[2:n, 2:n] + m[1:(n-1), 1:(n-1)]) / 2
    raw_sdr(m, c(0, 0.2), c(0, 1))
  }))
  
  expect_true({ 
    sdr_red > 1 
  })
  
})
