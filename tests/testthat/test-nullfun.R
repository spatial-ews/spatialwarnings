# 
# 
# 

context("Test that all null model methods work") 

# Here we just run the code to check that it works
test_that("All null model methods work", { 
  
  # We need to increase the max number of iterations otherwise warnings are produced
  options(spatialwarnings.constants.maxit = 10^(8.5))
  options(spatialwarnings.constants.reltol = 1e-8)
  
  # Check that all methods run 
  all_methods <- c("perm", "intercept", "smooth") 
  testmat <- matrix(runif(30*30) > .7, ncol = 30)
  testmat[1:15, ] <- TRUE
  
  dat <- list(testmat, testmat)
  a <- generic_sews(dat)
  b <- patchdistr_sews(dat)
  c <- suppressWarnings( spectral_sews(dat) )
  null_control <- list(family = binomial())
  for ( m in all_methods ) { 
    indictest(a, nulln = 3, null_method = m, 
              null_control = null_control) 
    indictest(b, nulln = 3, null_method = m, 
              null_control = null_control) 
    indictest(c, nulln = 3, null_method = m, 
              null_control = null_control) 
    expect_true(TRUE)
  }
  
  # Check the values returned by the null model 
  ictest <- indictest(a[[2]], 3, null_method = "intercept", 
                      null_control = null_control)
  nullmean <- mean(replicate(99, mean( ictest$get_nullmat() )))
  expect_equal(mean(ictest[["orig_data"]]), nullmean, 
               tol = 0.01)
  
  # Check the values returned by the null model 
  ictest <- indictest(a[[2]], 3, null_method = "smooth", 
                      null_control = null_control)
  nullmean <- mean(replicate(99, mean( ictest$get_nullmat() )))
  expect_equal(mean(ictest[["orig_data"]]), nullmean, 
               tol = 0.01)
  
  # Check that smoothed null model is closer to reality than the intercept
  # model
  ictest <- indictest(a[[2]], 3, null_method = "intercept", 
                      null_control = null_control)
  error_intercept <- replicate(49, { 
    mean( abs(ictest$get_nullmat() - a[[2]][["orig_data"]]) )
  })
  ictest <- indictest(a[[2]], 3, null_method = "smooth", 
                      null_control = null_control)
                      
  error_smooth <- replicate(49, { 
    mean( abs(ictest$get_nullmat() - a[[2]][["orig_data"]]) )
  })
  expect_true({ 
    mean(error_intercept) > mean(error_smooth)
  })
  
  # Check that we warn when the null method is a function and does not return 
  # logical values when the input matrix is logical 
  nullfun <- function(mat) { mat * rnorm(prod(dim(mat))) } 
  expect_warning({ 
    indictest(compute_indicator(serengeti[[1]], raw_cg_moran), 3, 
              null_method = nullfun)
  })
  
  # Check that arguments are properly set 
  expect_warning({ 
    null_control_set_args(serengeti[[1]], list(a = "opl"), "perm")
  })
  expect_true({ 
    a <- null_control_set_args(serengeti[[1]], 
                               list(family = "binomial"), 
                               "intercept")[["family"]]
    is.binomial(a)
  }) 
  expect_true({ 
    a <- null_control_set_args(serengeti[[1]], 
                               list(family = binomial()), 
                               "intercept")[["family"]]
    is.binomial(a)
  }) 
  
  # Model option are honored
  expect_warning({
      null_control_set_args(serengeti[[1]], list(), "intercept")
  }, regexp = "using a binomial")
  expect_warning({
    null_control_set_args(coarse_grain(serengeti[[1]], 5), 
                          list(), "intercept")
  }, regexp = "using a gaussian")
  
  # Set options to their default
  options(spatialwarnings.constants.maxit = NULL)
  options(spatialwarnings.constants.reltol = NULL)
  
  # Test the simulation of new values 
  
#   # Intercept-only models 
#   mu <- 10
#   dat <- rnorm(1000, mu, 2)
#   gmod <- glm(y ~ 1, data = data.frame(y = dat), family = gaussian())
#   newvals <- simulate_newdat(gmod, seq(0, 10000), family = gaussian())
#   expect_equal(mu, mean(newvals), tol = 0.1)
#   
#   lambda <- 10
#   dat <- rpois(1000, lambda)
#   gmod <- glm(y ~ 1, data = data.frame(y = dat), family = poisson())
#   newvals <- simulate_newdat(gmod, seq(0, 10000), family = poisson())
#   expect_equal(lambda, mean(newvals), tol = 0.1)
#   
#   prob <- 0.4
#   dat <- rbinom(1000, size = 1, prob = prob)
#   gmod <- glm(y ~ 1, data = data.frame(y = dat), family =  binomial())
#   newvals <- simulate_newdat(gmod, seq(0, 10000), family = binomial())
#   expect_equal(prob, mean(newvals), tol = 0.1)
#   
#   
#   # Helper function for whats following
#   matricize <- function(tab) { 
#     matrix(tab[ ,3], nrow = max(tab[ ,1]), ncol = max(tab[ ,1]))
#   }
#   
#   
#   # Smooth models: binomial
#   input <- serengeti[[7]]
#   # parameter selection for the spline. 
#   mat_tab <- data.frame(expand.grid(row = seq.int(nrow(input)), 
#                                     col = seq.int(ncol(input))), 
#                         value = as.vector(input))
#   null_mod <- mgcv::gam(value ~ s(row, col, bs = "tp"), 
#                         data = mat_tab, 
#                         family = binomial())
#   ovals <- expand.grid(row = seq.int(nrow(input)), 
#                       col = seq.int(ncol(input)))
#   newvals <- simulate_newdat(null_mod, 
#                             ovals, 
#                             family = binomial())
#   ovals[ ,"y"] <- newvals
#   ovals[ ,"yref"] <- simulate(null_mod)[, 1]
#   ggplot(ovals, aes(x = row, y = col)) + 
#     geom_raster(aes(fill = y > 0.5))
#   
#   m    <- coarse_grain(matricize(ovals[ ,c("row", "col", "y")]),    20)
#   mref <- coarse_grain(matricize(ovals[ ,c("row", "col", "yref")]), 20)
#   expect_true({ 
#     mean( abs(m - mref) ) < 0.05
#   })
#   
#   if ( exists("VISUAL_TESTS") && VISUAL_TESTS ) { 
#     display_matrix(matricize(ovals[ ,c("row", "col", "y")])) + 
#       labs(title = "produced")
#     dev.new()
#     display_matrix(matricize(ovals[ ,c("row", "col", "yref")])) + 
#       labs(title = "ref")
#   
#   }
#   
#   
#   # Smooth models: gaussian
#   input <- matrix(ifelse(serengeti[[7]][], 
#                          rnorm(prod(dim(serengeti[[7]])), mean = 20 , sd = 8), 
#                          rnorm(prod(dim(serengeti[[7]])), mean = 100, sd = 8)), 
#                   nrow = nrow(serengeti[[7]]), 
#                   ncol = ncol(serengeti[[7]]))
#   # parameter selection for the spline. 
#   mat_tab <- data.frame(expand.grid(row = seq.int(nrow(input)), 
#                                     col = seq.int(ncol(input))), 
#                         value = as.vector(input))
#   null_mod <- mgcv::gam(value ~ s(row, col, bs = "tp"), 
#                         data = mat_tab, 
#                         family = gaussian())
#   ovals <- expand.grid(row = seq.int(nrow(input)), 
#                        col = seq.int(ncol(input)))
#   newvals <- simulate_newdat(null_mod, 
#                              ovals, 
#                              family = gaussian())
#   ovals[ ,"y"] <- newvals
#   ovals[ ,"yref"] <- simulate(null_mod)[, 1]
#   
#   m    <- coarse_grain(matricize(ovals[ ,c("row", "col", "y")]),    20)
#   mref <- coarse_grain(matricize(ovals[ ,c("row", "col", "yref")]), 20)
#   expect_true({ 
#     cor(as.vector(m), as.vector(mref)) > 0.95
#   })
#   expect_equal(var(ovals$y), var(ovals$yref), tol = 4)
#   
#   if ( exists("VISUAL_TESTS") && VISUAL_TESTS ) { 
#     display_matrix(matricize(ovals[ ,c("row", "col", "y")])) + 
#       labs(title = "produced")
#     dev.new()
#     display_matrix(matricize(ovals[ ,c("row", "col", "yref")])) + 
#       labs(title = "ref")
#   
#   }
#   
#   
#   
#   
#   # Smooth models: poisson
#   input <- matrix(ifelse(serengeti[[7]][], 
#                          rpois(prod(dim(serengeti[[7]])), lambda = 20 ), 
#                          rpois(prod(dim(serengeti[[7]])), lambda = 100)), 
#                   nrow = nrow(serengeti[[7]]), 
#                   ncol = ncol(serengeti[[7]]))
#   # parameter selection for the spline. 
#   mat_tab <- data.frame(expand.grid(row = seq.int(nrow(input)), 
#                                     col = seq.int(ncol(input))), 
#                         value = as.vector(input))
#   null_mod <- mgcv::gam(value ~ s(row, col, bs = "tp"), 
#                         data = mat_tab, 
#                         family = gaussian())
#   ovals <- expand.grid(row = seq.int(nrow(input)), 
#                        col = seq.int(ncol(input)))
#   newvals <- simulate_newdat(null_mod, 
#                              ovals, 
#                              family = gaussian())
#   ovals[ ,"y"] <- newvals
#   ovals[ ,"yref"] <- simulate(null_mod)[, 1]
#   
#   m    <- coarse_grain(matricize(ovals[ ,c("row", "col", "y")]),    20)
#   mref <- coarse_grain(matricize(ovals[ ,c("row", "col", "yref")]), 20)
#   expect_true({ 
#     cor(as.vector(m), as.vector(mref)) > 0.95
#   })
#   expect_equal(var(ovals$y), var(ovals$yref), tol = 4)
#   
#   if ( exists("VISUAL_TESTS") && VISUAL_TESTS ) { 
#     display_matrix(matricize(ovals[ ,c("row", "col", "y")])) + 
#       labs(title = "produced")
#     dev.new()
#     display_matrix(matricize(ovals[ ,c("row", "col", "yref")])) + 
#       labs(title = "ref")
#   }
  
  
  
  
  
  
})

# 
# 
# img <- serengeti[[length(serengeti)]]
# img_coarse <- coarse_grain(img, 1)
# img_tab <- data.frame(expand.grid(row = seq.int(nrow(img_coarse)), 
#                                   col = seq.int(ncol(img_coarse))), 
#                       value = as.vector(img_coarse))
# 
# 
# # Test if trend 
# mod <- mgcv::gam(value ~ s(row, col, bs = "tp"), data = img_tab, 
#                  family = binomial())
# img_tab[ ,"pred"] <- predict(mod, type = "response")
# img_tab[ ,"sim"]  <- simulate(mod, type = "response") > .5
# 
# 
# 
# plot1 <- ggplot(img_tab) + 
#   geom_raster(aes(x = col, y = row, fill = value)) + 
#   scale_fill_brewer(palette = "Spectral") + 
#   coord_fixed() + 
#   theme_minimal() + 
#   labs(x = "x", y = "y", title = "Observed matrix")
# 
# plot2 <- ggplot(img_tab) + 
#   geom_raster(aes(x = col, y = row, fill = 1 - pred)) + 
#   scale_fill_distiller(palette = "Spectral", 
#                        name = "Cover", 
#                        direction = -1) + 
#   coord_fixed() + 
#   theme_minimal() + 
#   labs(x = "x", y = "y", title = "Fitted model")
# 
# plot3 <- ggplot(img_tab) + 
#   geom_raster(aes(x = col, y = row, fill = sim)) + 
#   scale_fill_brewer(palette = "Spectral") + 
#   coord_fixed() + 
#   theme_minimal() + 
#   labs(x = "x", y = "y", title = "One null matrix")
# 
# 
# 
# library(patchwork)
# plot1 + plot2 + plot3 + 
#   plot_layout(nrow = 1)
# 
# 
# # Create a null matrix 
# library(mgcv)
# library(memoise)
# fit_model <- memoise(function(data) { 
#   mgcv::gam(value ~ s(row, col, bs = "tp"), data = data, 
#             family = binomial())
# })
# fnull <- function(mat) { 
#   mat_tab <- data.frame(expand.grid(row = seq.int(nrow(mat)), 
#                                     col = seq.int(ncol(mat))), 
#                         value = as.vector(mat))
#   mod <- fit_model(data = mat_tab) 
#   mat[ , ] <- simulate(mod) > .5 
#   return(mat)
# }
# 
# indic <- compute_indicator(serengeti, raw_moran)
# test  <- indictest(indic, null_method = fnull, nulln = 99)
# plot(test, along = serengeti.rain)
