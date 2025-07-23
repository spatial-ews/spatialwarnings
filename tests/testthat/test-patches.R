# 
# Misc tests with working with patch sizes
# 

test_that("Labelling functions & friends work correctly", { 
  
  n <- 512

  # Fix seed to one that does not produce warning in patchdistr_sews()
  oseed <- .Random.seed
  set.seed(123)
  m <- matrix(rnorm(n*n) > 0.3, nrow = n, ncol = n)
  .Random.seed <- oseed

  # Make sure args are dispatched correctly
  X <- patchsizes(list(m, m), nbmask = "moore", wrap = TRUE)
  expect_true({ 
    all(X[[1]] == X[[2]])
  })
  
  
  # Make sure the type of neighborhood is passed even when using high-level 
  # functions. We remove warnings because we fit on random data that can sometimes
  options(spatialwarnings.constants.maxit = 1e9L)
  x1 <- patchdistr_sews(m)
  x2 <- patchdistr_sews(m, nbmask = "moore")
  expect_true({ 
    x1[["npatches"]] > x2[["npatches"]]
  })
  options(spatialwarnings.constants.maxit = NULL)
  
  m <- matrix(c(1, 0, 0, 
                0, 1, 1, 
                0, 0, 0) > 0, 
              ncol = 3, nrow = 3, byrow = TRUE)
  
  expect_true({ 
    length(patchsizes(m)) == 2 && all(patchsizes(m) == c(1, 2))
  })
  
  expect_true({ 
    length(patchsizes(m, nbmask = "moore")) == 1 && patchsizes(m, nbmask = "moore") == 3
  })
  
  m <- matrix(c(1, 0, 0, 
                0, 0, 1, 
                0, 0, 0) > 0, 
              ncol = 3, nrow = 3, byrow = TRUE)
  
  expect_true({ 
    a <- patchdistr_sews(m, nbmask = "von_neumann", wrap = TRUE)
    a[["npatches"]] == 2
  })
  
  expect_true({ 
    a <- patchdistr_sews(m, nbmask = "moore", wrap = TRUE)
    a[["npatches"]] == 1
  })
  
  
  
})
