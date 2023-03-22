
context('Test that xmin estimations are correct')

if ( exists('EXTENDED_TESTS') && EXTENDED_TESTS ) { 
  
  # Change dir if running tests manually
  if ( file.exists('./tests/testthat') ) { 
    library(testthat)
    setwd('./tests/testthat') 
  }
  
  # Setup pli from Clauzet et al's
  for ( s in dir('./pli-R-v0.0.3-2007-07-25', 
                full.names = TRUE, pattern = '*.R') ) { 
    source(s)
  }
  
  # Compile auxiliary binaries
  system("cd ./pli-R-v0.0.3-2007-07-25/zeta-function/ && make")
  system("cd ./pli-R-v0.0.3-2007-07-25/exponential-integral/ && make")
  system("cd ./pli-R-v0.0.3-2007-07-25/ && \
            gcc discpowerexp.c -lm -o discpowerexp && \
            chmod +x discpowerexp")
  
  test_that("xmins estimation is correct", { 
    
    parms <- expand.grid(expo = c(1.5, 1.3, 1.2), 
                         rate = c(0.001, 0.005, 0.01, 0.1, 0.15, 0.2, 0.25, 0.3, 
                                  0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
    
    estim_xmin <- function(df) { 
      
      pldat <- round(rpowerexp(10000, 1, df[ ,'expo'], df[, 'rate']))
      
      est_xmin <- suppressWarnings( spatialwarnings:::xmin_estim(pldat) )
      # est_xmin <- spatialwarnings:::xmin_estim(pldat) 
      
      # Create pl object and estimate its xmin
      pl_obj <- poweRlaw::displ$new(pldat)
      est_xmin_plpkg <- poweRlaw::estimate_xmin(pl_obj)[["xmin"]]
      
      if ( ! is.na(est_xmin) && 
           ! is.na(est_xmin_plpkg) && 
           length(unique(pldat)) >= 5 ) { 
        # Note: In some pathological cases (few unique patches), there can be 
        # a small difference in xmin, so we use an acceptable error here. 
        fit_is_ok <- ( abs(est_xmin - est_xmin_plpkg) <= 1 ) 
        
        if ( est_xmin != est_xmin_plpkg ) { 
          cat(" Ours: ", est_xmin, ' -> poweRlaw\'s: ', est_xmin_plpkg, " [", 
              length(unique(pldat)), " unique patches]", "\n", sep = "")
        }
        
        # In this case, inspect the fit provided by the poweRlaw package
        if ( exists('GRAPHICAL') && GRAPHICAL && est_xmin != est_xmin_plpkg ) { 
          dev.new()
          par(mfrow = c(1, 2))
          plot(log10(cumpsd(pldat[pldat >= est_xmin])))
          our_fit <- pl_fit(pldat, xmin = est_xmin)
          xs <- unique(round(seq(min(pldat), max(pldat), length.out = 100)))
          lines(log10(xs), 
                log10(ippl(xs, our_fit$plexpo, xmin = est_xmin)))
          title("OUR FIT")
          
          plot(log10(cumpsd(pldat[pldat >= est_xmin_plpkg])))
          plpkg_expo <- poweRlaw::estimate_xmin(pl_obj)$pars
          xs <- unique(round(seq(min(pldat), max(pldat), length.out = 100)))
          lines(log10(xs), 
                log10(ippl(xs, plpkg_expo, xmin = est_xmin_plpkg)))
          title("PWL FIT")
          browser()
        }
      }
      
      data.frame(df, est_xmin = est_xmin, fit_is_ok = fit_is_ok)
    }
    
    xmin_ests <- ddply(parms, ~ expo + rate, estim_xmin)
    
    # Check that fits are overwhelmingly good. Sometimes the poweRlaw code will
    # produce different fits, so we need to be able to accomodate some errors here.
    # We put five as a maximum tolerable number of mismatches (out of 45)
    number_of_ok_fits <- sum(xmin_ests[ ,"fit_is_ok"])
    cat(paste0("Est. xmins OK: ", number_of_ok_fits, "/", nrow(parms), "\n"))
    expect_true({ 
      ( nrow(parms) - number_of_ok_fits ) <= 5
    })
    
  })
  
  # Remove auxiliary binaries now that tests are done
  system("cd ./pli-R-v0.0.3-2007-07-25/zeta-function/ && rm zeta_func zeta_func.o")
  system("cd ./pli-R-v0.0.3-2007-07-25/exponential-integral/ && rm exp_int exp_int.o")
  system("cd ./pli-R-v0.0.3-2007-07-25/ && rm discpowerexp")
  
} else { 
  message('Skipping xmin estimation testing')
}
