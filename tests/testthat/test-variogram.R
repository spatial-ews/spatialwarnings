# 
# This file tests the variogram-related low-level functions and indicators
# 

# Make sure it works
context("Computation of variograms")

test_that("Computation of variograms works", { 
  
  # Here so that testthat does not complain about empty tests when missing packages
  expect_true(TRUE)
  
  if ( requireNamespace("gstat", quietly = TRUE) && 
       requireNamespace("sp", quietly = TRUE) ) { 
    if ( exists("EXTENDED_TESTS") && EXTENDED_TESTS) { 
      mats <- list(arizona[[4]], forestgap[[8]], serengeti[[6]]) 
    } else { 
      mats <- list(forestgap[[8]])
    }
    
    for ( example_mat in mats ) { 
      example_mat = serengeti[[6]]
      nbins <- 32
      nmax = 2e6L # 5e5L # prod(dim(example_mat)) ^2
      cutoff = sqrt(ncol(example_mat)^2 + nrow(example_mat)^2) / 10
      
      vario1 <- variogram_internal(example_mat, nmax, nbins, cutoff)
#       pairs(vario1)
      
      # Compare results to gstat
      locations <- expand.grid(seq.int(nrow(example_mat)),
                              seq.int(ncol(example_mat)))
      locations <- locations[sample.int(nrow(locations), 
                                        replace = FALSE, 
                                        size = min(nrow(locations),nmax/10)), ]
      values <- apply(locations, 1, function(X) example_mat[X[1], X[2]])
      locations.gstat <- sp::SpatialPointsDataFrame(locations, 
                                                    data.frame(z = values))
      vario2 <- gstat::variogram(z ~ 1, 
                                 data   = locations.gstat, 
                                 width  = cutoff / nbins, 
                                 cutoff = cutoff)
      plot(vario1[ ,"dist"], vario1[ ,"gamma"] )
      points(vario2[ ,"dist"], vario2[ ,"gamma"], col = "red", pch = 20)
      
      vario1_2 <- approx(vario1[ ,1], vario1[ ,2], vario2[ ,"dist"])$y
      
#       plot(vario1_2 - vario2[ ,"gamma"], vario2[ ,"gamma"])
#       abline(0,1)
      
      # The difference between the two variograms should be zero on average
      mod <- coef(lm(vario1_2 - vario2[ ,'gamma'] ~ 1))
      expect_true(abs(mod[1]) < 0.01) # Intercept ~= 0
      
      # Test if P-value is significant, we expect it not to be
      expect_true(kruskal.test(vario1_2, vario2[ ,"gamma"])$p.value > 0.1)
    }
    
  } # end of requireNamespace()
  
})

