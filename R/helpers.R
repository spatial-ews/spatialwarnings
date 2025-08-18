# 
# 

# Function that transforms a set of values into stars for display
pval_stars <- function(value, NA_ret = NA) { 
  
  if ( length(value) > 1 ) { 
    stars <- sapply(value, pval_stars, NA_ret = NA_ret) 
    attributes(stars) <- attributes(value)
    return(stars)
  }
  
  if ( is.na(value) ) { 
    return(NA_ret) 
  }
  
  if (value < 0.001) { 
    return('***' )
  } else if (value < 0.01) { 
    return('** '  )
  } else if (value < 0.05) { 
    return('*  '   )
  } else if (value < 0.1) { 
    return('.  '   )
  } else { 
    return('   ')
  }
}

# Format pvalues to reflect the number of null values instead of printing a 
# 0 when it is below what's detectable. 
format_pvalues <- function(X, nreps) { 
  ifelse(X == 0, 
         paste0("<", formatC(1/nreps, 
                             width = 3, 
                             digits = 0, 
                             format = "e")), 
         formatC(X, format = "f", digits = 3))
}

# Print a line with matrix size information
#' @title Display size information of a sews_result object
#'
#' @param x A \code{sews_result} object (as provided by **_sews functions, such
#'   as \code{generic_sews()}
#'
#' @param ... Other arguments are ignored.
#'
#' @details This function is used to display information about the size of matrices
#'   used by various functions in spatialwarnings. It is mostly destined to internal
#'   use.
#'
#'@export
display_size_info <- function(x, ...) {
  UseMethod('display_size_info')
}

# ... for sews_result_single class
#'@rdname display_size_info
#'@exportS3Method
display_size_info.sews_result_single <- function(x, ...) {
  display_size_info(list(x))
}

# ... for sews_result_list class
#'@rdname display_size_info
#'@exportS3Method
display_size_info.sews_result_list <- function(x, ...) {
  sizes <- sapply(x, function(x) dim(x[["orig_data"]]))
  sizes <- apply(sizes, 1, function(X) length(unique(X)) == 1)
  has_different_sizes <- ! any(sizes)
  if (has_different_sizes) {
    size_text_report <- "variable sizes"
  } else {
    size_text_report <- paste0("size: ", nrow(x[[1]][["orig_data"]]), 'x',
                               ncol(x[[1]][["orig_data"]]))
  }
  cat(' ',
      length(x), ' ',
      ifelse(length(x)>1, 'matrices', 'matrix'), ' ',
      "(", size_text_report,')\n', sep = '')
}

# ...for generic lists
#'@rdname display_size_info
#'@exportS3Method
display_size_info.list <- function(x, ...) {
  display_size_info.sews_result_list(x, ...)
}

ifNULLthen <- function(a, b) { 
  if (is.null(a)) b else a 
}

list_methods <- function(class, 
                         exclude = c("print", "display_size_info", 
                                     "summary")) { 
  all_methods <- lapply(class, function(class) { 
    tab <- attr(methods(class = class), "info")
#     tab[tab[ ,"from"] == "spatialwarnings", "generic"]
    tab[ ,"generic"]
  })
  all_methods <- sort(unique(unlist(all_methods)))
  
  # For _test_ objects, we do not report indictest
  if ( any(grepl("_sews_test_", class)) ) { 
    all_methods <- all_methods[!grepl("indictest", all_methods)]
  }
  
  # Exclude some reported methods 
  all_methods <- all_methods[! all_methods %in% exclude]
  
  return(all_methods)
}

# Wrapper around future_lapply to set future.seed = TRUE so that random number 
# generation in child processes is correct. 
future_lapply_seed <- function(...) { 
  future.apply::future_lapply(..., future.globals = FALSE, future.seed = TRUE)
}
