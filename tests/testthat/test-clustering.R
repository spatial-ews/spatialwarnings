# 
# This tests the computation of clustering 
# 

test_that("Computation of clustering are OK", { 

  mat <- matrix(c(1, 0, 0, 
                  0, 0, 0, 
                  0, 0, 0), 
                ncol = 3, byrow = TRUE)
  
  
  # With 4 nbs and wrapping, we have 14 pairs of 0-0, 0 of 1-1
  clusts <- clustering_core(mat, ns = 2, wrap = TRUE, use_8_nb = FALSE)
  expect_true({ 
    all(c(clusts[ ,1] == c(14, 0), 
          clusts[ ,2] == c(8, 1) )) # counts of single cells
  })
  
  # With 8 nbs and wrapping, we have 21 pairs of 0-0, 0 of 1-1
  clusts <- clustering_core(mat, ns = 2, wrap = TRUE, use_8_nb = TRUE)
  expect_true({ 
    all(c(clusts[ ,1] == c(28, 0), 
          clusts[ ,2] == c(8, 1) )) # counts of single cells
  })
  
  # With 4 nbs and no wrapping, we have ten pairs of 0-0, 0 of 1-1
  clusts <- clustering_core(mat, ns = 2, wrap = FALSE, use_8_nb = FALSE)
  expect_true({ 
    all(c(clusts[ ,1] == c(10, 0), 
          clusts[ ,2] == c(8, 1) )) # counts of single cells
  })
  
  # With 8 nbs and no wrapping, we have 13 pairs of 0-0, 0 of 1-1
  clusts <- clustering_core(mat, ns = 2, wrap = FALSE, use_8_nb = TRUE)
  expect_true({ 
    all(c(clusts[ ,1] == c(17, 0), 
          clusts[ ,2] == c(8, 1) )) # counts of single cells
  })
  
  
  
  # Make sure clustering of random matrices is close to one 
  clusts <- replicate(199, { 
    m <- matrix(sample(letters[1:4], size = 100^2, replace = TRUE), 
                nrow = 100, ncol = 100)
    raw_clustering(m) 
  })
  expect_true({ 
    all( abs(apply(clusts, 1, mean) - 1) < 0.10 ) 
  })
  
})

