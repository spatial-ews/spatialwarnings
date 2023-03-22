# 
# This function makes sure parallelism is enabled 
#

context("Test that parallel computing works")

test_that("Parallelism work", { 
  
  if ( exists("EXTENDED_TESTS") && 
       EXTENDED_TESTS && 
       availableCores() > 1 ) { 
    # We increase the dataset size because otherwise the overhead of setting up workers
    # may make the parallel verison slower overall
    a <- generic_sews( c(forestgap, forestgap) )
    
    plan(sequential)
    b.1 <- system.time( indictest(a, 49) ) 
    
    plan(multisession)
    b.2 <- system.time( indictest(a, 49) ) 
    
    expect_true( b.1["elapsed"] > b.2["elapsed"] )
    
    plan(sequential) # restore plan to no-parallelism
  } else { 
    # This is a dummy test so that testhat does not report an empty test
    expect_true(TRUE)
  }
  
})

