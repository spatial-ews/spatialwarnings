# 
# Misc tests with working with patch sizes
# 

test_that("Labelling functions & friends work correctly", { 
  
  n <- 512
  m <- matrix(rnorm(n*n) > 0.3, nrow = n, ncol = n)
  
  # Make sure args are dispatched correctly
  X <- patchsizes(list(m, m), nbmask = "moore", wrap = TRUE)
  expect_true({ 
    all(X[[1]] == X[[2]])
  })
  
  # Make sure the type of neighborhood is passed even when using high-level 
  # functions 
  options(spatialwarnings.constants.maxit = 1e9L)
  x1 <- patchdistr_sews(m)
  x2 <- patchdistr_sews(m, nbmask = "moore")
  expect_true({ 
    x1[["npatches"]] > x2[["npatches"]]
  })
  
  options(spatialwarnings.constants.maxit = NULL)
  
  
  
})
