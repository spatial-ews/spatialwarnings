# 
# 
# 

context("Test that matrices can be displayed")

isgg <- function(o) { 
  inherits(o, c("ggplot", "ggplot2::ggplot"))
}

test_that("display_matrix methods work", { 
  
  # Test the matrix methods
  expect_true({ 
    a <- display_matrix(serengeti[[length(serengeti)]])
    isgg(a)
  })
  expect_true({ 
    a <- display_matrix(coarse_grain(serengeti[[length(serengeti)]], 5))
    isgg(a)
  })
  
  # Test the sews methods
  sets <- suppressWarnings({ # spectral_sews will produce warnings
    list(spectral_sews(serengeti[2:3]), 
         compute_indicator(serengeti[2:3], raw_moran))
  })
  
  for ( a in sets ) { 
    expect_true({ 
      all(c(isgg( display_matrix(a, along = NULL) ), 
            isgg( display_matrix(a, along = letters[2:3]) ), 
            isgg( display_matrix(a[[1]]) )))
    })
    
    expect_true({
      all(c(isgg( display_matrix(indictest(a, nulln = 3), along = NULL) ), 
            isgg( display_matrix(indictest(a, nulln = 3), 
                                 along = letters[2:3]) ), 
            isgg( display_matrix(indictest(a, nulln = 3)[[1]]) )))
    })
  }
  
})
