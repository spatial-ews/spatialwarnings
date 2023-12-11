# 
# Test the implementation of the LSW distribution
# 


test_that("LSW distribution works", { 
  
  cLSW <- function(r, mu) {
    x <- seq(0, 1.5, 0.01) * mu
    PDF <- dLSW(x, mu)
    CDF <- cumsum(PDF) / sum(PDF)
    y <- approx(x, CDF, r)[["y"]]
    
    return(y)
  }

  rs <- seq(0, 20, l = 512)
  mus <- seq(1, 10, l = 5)
  
  # Multiple values for mu is not supported
  expect_error({ 
    dLSW(rs, mus)
  })
  
  for ( mu in mus ) { 
    
    # Make sure log of dLSW is correct 
    logds <- log(dLSW(rs, mu))
    dslog <- dLSW(rs, mu, log = TRUE)
    
    # Compare only values where the log(p) is finite
    expect_true({ 
      all( na.omit(abs(logds[is.finite(logds)] - dslog[is.finite(logds)]))  < 1e-10 )
    })
    
    # Make sure pLSW works
    ps <- pLSW(rs, mu)
    psinv <- pLSW(rs, mu, lower.tail = FALSE)
    expect_true({ 
      all( abs(ps - ( 1 - psinv )) < 1e-10)
    })
    
    # 
    pslog <- pLSW(rs, mu, log.p = TRUE)
    diffs <- suppressWarnings({ 
      ifelse(!(ps > 0), 0, log(ps) - pslog)
    })
    expect_true({ 
      all( abs(diffs) < 1e-10)
    })
    
    # Compare with original functions. There is actually a pretty large difference, 
    # as the original one was quite approximate. We remove NAs that are returned 
    # by cLSW, and only compare values that can be compared.
    orig_ps <- cLSW(rs, mu)
    expect_true({ 
      all( na.omit(abs(ps - orig_ps)) < 0.1 )
    })
    
  }
  
})

test_that("LSW fitting recovers correct values", { 
  
  fits <- vapply(dda, function(m) { 
    LSW_fit(patchsizes(m))[["mu"]]
  }, numeric(1))
  
  # As tau increases, we expect mu to increase. This is a very light test that 
  # will only catch gross errors, but we are not able to produce random samples 
  # from LSW distrib, so we cannot do an in-depth test of the fitting. 
  cor <- cor.test(dda.pars[ ,"tau"], log(fits))[["estimate"]]
  expect_true({ 
    cor > 0
  })
  
})

