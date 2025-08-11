# 
# Test for the Kolmogorov complexity
# 
if ( requireNamespace("acss", quietly = TRUE) ) { 

context('Test that the computations for Kolmogorov complexity are correct')

test_that("Kbdm results are correct", { 
  
  skip_on_cran()

  # We can't reproduce exactly the results because we cannot compute 
  # Kbdm using 4x4 submatrices in acss::acss. We just compute the test on 
  # the scale-free patterns which do not show as much noise as the 
  # irregular/turing patterns. 
  dat <- read.csv("./kbdm_data/scale-free/extracted_data.csv", 
                  header = FALSE)

  allfiles <- dir("./kbdm_data/scale-free/", 
                  pattern = "^.*.txt$", full.names = TRUE)
  matrices <- lapply(allfiles, 
                    function(f) { 
                          a <- as.matrix(read.table(f, header = FALSE))
                          attr(a, "dimnames") <- NULL
                          a > mean(a) # binarize as done in Dakos 2017
                          })

  kbdm_vals <- as.data.frame(kbdm_sews(matrices))

  xvals <- sapply(allfiles, function(f) { 

      a <- strsplit(f, split = c("CA_b"))[[1]][2]
      a <- gsub(".txt", "", a)
      a <- as.numeric(a)
      a
    }, USE.NAMES = FALSE)

  # We cannot compare the absolute values, but we compare the ranks 
  ref <- rank(dat[ ,2])
  val <- rank(kbdm_vals[order(xvals),"value"])

  plot(ref, val)
  abline(0, 1, col = "red")
  linmod <- lm(y ~ x, data = data.frame(x = ref, y = val))
  lines(ref, predict(linmod))

  # Slope of relationship should be close to zero
  expect_true({
    abs(1 - coef(linmod)[2]) < .1
  })

  # Both values should have a max that is very close 
  expect_true({ 
    abs(xvals[ref == max(ref)] - xvals[val == max(val)]) < .1
  })

  
})

}
