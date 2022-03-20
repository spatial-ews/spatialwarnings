

context('Test the computation of moran correlation')

test_that("Moran correlation is computed correctly", { 
  
  within_bounds <- function(x) { 
    expect_true( all(x >= -1 & x <= 1) )
  }
  
  # We work on small matrices so that the handling of sides has a big effect
  n <- 5*5
  nreps <- 9999
  
  # Correlation of random matrices is near-zero 
  morans <- replicate(nreps, { 
    raw_moran(matrix(rnorm(n), ncol = sqrt(n)))
  })
  within_bounds(morans)
  # Check that the mean value is equal to the expected value
  expect_true( abs(mean(morans) - ( -1 / ( n - 1 ) )) < 1e-2)
  
  # Correlation of checkerboard is near -1 
  # If n is even, then make it odd, so that the testing matrices are really checkerboards
  ncheck <- ifelse(n %% 2 == 0, (sqrt(n)+1)^2, n)
  c1 <- matrix((1+seq.int(ncheck)) %% 2, ncol = sqrt(ncheck))
  c2 <- matrix(seq.int(ncheck) %% 2, ncol = sqrt(ncheck))
  expect_true( abs(raw_moran(c1) - -1) < 1e-10 )
  expect_true( abs(raw_moran(c2) - -1) < 1e-10 )
  within_bounds(raw_moran(c1))
  within_bounds(raw_moran(c2))
  
  # Correlation of left/right matrix is near +1
  nm <- 200*200
  m1 <- cbind(matrix(1, nrow = sqrt(nm), ncol = floor(sqrt(nm) / 2)), 
              matrix(0, nrow = sqrt(nm), ncol = ceiling(sqrt(nm) / 2)))
  expect_true( abs( raw_moran(m1) - 1 ) < 0.01 )
  within_bounds(raw_moran(m1))
  
  if ( exists("EXTENDED_TESTS") && EXTENDED_TESTS ) { 
    # Test our implementation against the one in package raster
    # Note that raster produces warnings that do not affect us
    if ( suppressWarnings( requireNamespace("raster", quietly = TRUE) ) ) { 
      c <- seq(0.01, 0.99, l = 21)
      w <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), byrow = TRUE, ncol = 3)
      
      for ( i in seq_along(c) ) { 
        m <- matrix(runif(n) < c[i], ncol = sqrt(n))
        diffidx <- abs( raw_moran(m) - raster::Moran(raster::raster(m), w) )
        if ( ! is.nan(diffidx) ) { 
          expect_true( diffidx < 1e-10 )
        }
      }
      expect_true({ 
        abs(raw_moran(c1) - raster::Moran(raster::raster(c1), w)) < 1e-10
      })
      expect_true({ 
        abs(raw_moran(c2) - raster::Moran(raster::raster(c2), w)) < 1e-10
      })
      
    }
  }
  
})

