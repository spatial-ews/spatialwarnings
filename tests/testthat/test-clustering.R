#
# This tests the computation of clustering
#

test_that("Neighbor-counting is OK", {

  # Check that total number of pairs matches counting done by hand
  expect_true( total_pair_count(diag(2), TRUE, TRUE) == 16  )
  expect_true( total_pair_count(diag(2), FALSE, TRUE) == 6  )
  expect_true( total_pair_count(diag(2), FALSE, FALSE) == 4 )
  expect_true( total_pair_count(diag(2), TRUE, FALSE) == 8  )
  
  # Check that total number of pairs matches a manual sum after the fact
  for ( wrap in c(TRUE, FALSE) ) { 
    for ( use_8_nb in c(TRUE, FALSE) ) { 
      a <- pair_counts(diag(10), wrap = wrap, use_8_nb = use_8_nb)
      b <- pair_counts(diag(10), prop = FALSE, wrap = wrap, use_8_nb = use_8_nb) 
      b <- b / sum(b, na.rm = TRUE)
      expect_true({ 
        all( abs(a[!is.na(a)] - b[!is.na(b)]) < 1e-8 )
      })
    }
  }
  
})

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
  if ( exists("EXTENDED_TESTS") && EXTENDED_TESTS ) {
    clusts <- replicate(199, {
      m <- matrix(sample(letters[1:4], size = 100^2, replace = TRUE),
                  nrow = 100, ncol = 100)
      raw_clustering(m)
    })
    expect_true({
      all( abs(apply(clusts, 1, mean) - 1) < 0.10 )
    })
  }

})

test_that("Counting of pairs makes sense", {
  n <- 99
  counts <- lapply(seq.int(n), function(n) { 
    m <- matrix(sample(letters[1:4], size = 100^2, replace = TRUE),
                nrow = 100, ncol = 100)
    pair_counts(m, prop = TRUE)
  })
  
  means <- Reduce(`+`, counts) / n
  
  # Sum should be 1
  expect_true({ 
    abs( sum(means, na.rm = TRUE) - 1 ) < 0.00001 
  })
  
  # Diags should be equal
  expect_true({ 
    all(abs(diff(diag(means))) < 0.1)
  })
  
  # Lower tri values should be equal
  expect_true({ 
    all( abs(diff(means[lower.tri(means)])) < 0.01)
  })
  
})
