# 
# Methods for a simple indicator, i.e. an indicator that produces on or more 
# numeric values for each matrix. 
# 



# as.df methods
# ---------------
#'@method as.data.frame simple_sews_single
#'@export
as.data.frame.simple_sews_single <- function(x, ..., wide = FALSE) { 
  as.data.frame.simple_sews_list( list(x), wide = wide )
}
#'@method as.data.frame simple_sews_list
#'@export
as.data.frame.simple_sews_list <- function(x, ..., wide = FALSE) { 
  
  # Find or create the indicator names
  indicnames <- ifNULLthen(names(x[[1]][['value']]), 
                           paste0("indic_", seq_along(x[[1]][['value']])))
  
  if ( wide ) { 
    output <- Map(function(n, o) { 
        a <- as.data.frame(matrix(o[['value']], nrow = 1))
        names(a) <- indicnames
        data.frame(matrixn = n, a)
      }, seq_along(x), x)
  } else { 
    output <- Map(function(n, o) { 
        data.frame(matrixn = n, indic = indicnames, value = o[['value']])
      }, seq_along(x), x)
  }
  output <- do.call(rbind, output)
  row.names(output) <- NULL
  output
}



# Print methods
# ---------------
#'@method print simple_sews_single
#'@export
print.simple_sews_single <- function(x, ...) { 
  summary.simple_sews_single(x, ...)
}
#'@method print simple_sews_list
#'@export
print.simple_sews_list <- function(x, ...) { 
  summary.simple_sews_list(x, ...)
}



# Summary methods
# ---------------
#'@method summary simple_sews_single
#'@export
summary.simple_sews_single <- function(object, 
                                       indicname = object[["taskname"]], 
                                       ...) { 
  cat('Spatial Early-Warning:', indicname, '\n') 
  cat('\n')
  display_size_info(object)
  cat('\n')
  
  # Format output table
  output <- as.data.frame(object, wide = TRUE)
  names(output)[1] <- c('Mat. #')
  
  print.data.frame(output, row.names = FALSE, digits = DIGITS)
  
  cat('\n')
  cat("The following methods are available: \n")
  cat(list_methods(class(object)), "\n")
}
#'@method summary simple_sews_list
#'@export
summary.simple_sews_list <- function(object, 
                                     indicname = object[[1]][["taskname"]], 
                                     ...) { 
  
  cat('Spatial Early-Warning:', indicname, '\n') 
  cat('\n')
  display_size_info(object)
  cat('\n')
  
  # Format output table
  output <- as.data.frame.simple_sews_list(object, wide = TRUE)
  names(output)[1] <- c('Mat. #')
  
  print.data.frame(output, row.names = FALSE, digits = DIGITS)
  
  cat('\n')
  cat("The following methods are available: \n")
  cat(list_methods(class(object)), "\n")
  
  invisible(output)
}


# Plot methods 
# ------------
# /!\ the doc for this function is documented in ./R/simple_indic_indictest.R
#'@rdname simple_sews_methods
#'@method plot simple_sews_list
#'@export
plot.simple_sews_list <- function(x, along = NULL, ...) { 
  plot.simple_sews_test_list(x, along = along, display_null = FALSE)
}



