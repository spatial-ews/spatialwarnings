# 
#
#' @title Custom Spatial Early-Warning signals
#' 
#' @aliases custom_indicator
#' 
#' @description Computation, significance assessment and display of trends 
#'   of a custom, user-defined indicator.
#' 
#' @param fun A function that takes a matrix as input and returns 
#'   a vector of numerical values. If the function returns a named vector, 
#'   then the names will be used in plots and summaries. The function may also
#'   accept extra arguments. 
#' 
#' @param taskname The task name. A character string used used for plots and
#'   textual summaries that describes the indicator (or set of indicators) 
#'   being computed. If a task name cannot be derived from \code{fun}, then 
#'   a default name is used. 
#' 
#' @return 
#' 
#' \code{create_indicator} returns a function that can be used in the same way 
#'   than the other \code{*_sews} functions (e.g. \code{generic_sews}). This 
#'   function as well as \code{compute_indicator} will return
#'   \code{\link[=simple_sews]{simple_sews_*}} objects. 
#' 
#' @details 
#' 
#' \code{spatialwarnings} provides "workflow functions", named \code{*_sews}, 
#'   that assist the user in computing, displaying and assessing the 
#'   significance of indicator values. The functions \code{create_indicator} and 
#'   \code{compute_indicator} provides such workflow for any arbitrary function. 
#'   
#' \code{create_indicator} takes a function `fun` and returns another function 
#'   that can be used as an indicator similar to the \code{*_sews} functions. The 
#'   results of this function can be assessed for significance using 
#'   \code{indictest} and trends can be displayed using 
#'   \code{plot}, \code{summary}, etc. (see Examples). \code{compute_indicator} 
#'   does the same but without needing an intermediate indicator function. 
#' 
#' @seealso \code{\link{simple_sews}}
#' 
#' @examples
#' 
#' # Use the maximum patch size as indicator of degradation
#' maxpatchsize <- function(mat) { 
#'   max(patchsizes(mat))
#' }
#' 
#' # Create the indicator function
#' maxpatch_sews <- create_indicator(maxpatchsize)
#' 
#' # Then work with this function as if it were a function from the *_sews 
#' # family. 
#' mp_indic <- maxpatch_sews(forestgap)
#' summary(mp_indic)
#' 
#' \donttest{ 
#' # Assess significance and display trends
#' mp_test <- indictest(mp_indic, nulln = 49)
#' plot(mp_test)
#' }
#' 
#' 
#' 
#' # Try spatial coefficient of variation as a spatial EWS. This function can 
#' # have arguments. 
#' spatial_cv <- function(mat, subsize) { 
#'   matc <- coarse_grain(mat, subsize)
#'   return( sd(matc) / mean(matc) )
#' }
#' 
#' # Create indicator function
#' cv_sews <- create_indicator(spatial_cv)
#' 
#' # Compute and display trends
#' cv_indic <- cv_sews(serengeti, subsize = 3)
#' plot(cv_indic, along = serengeti.rain)
#' 
#' # We can do the same work in one go using compute_indicator
#' cv_indic2 <- compute_indicator(serengeti, spatial_cv, subsize = 3)
#' plot(cv_indic2, along = serengeti.rain)
#' 
#' \donttest{ 
#' indictest(cv_indic, nulln = 99)
#' }
#'@export
create_indicator <- function(fun, 
                             taskname = as.character(substitute(fun))) { 

  # Test if taskname is derived from an anonymous function
  if ( length(taskname) > 1 ) { 
    warning('A valid name could not be derived from the function passed', 
            'to create_indicator, using a default name')
    taskname <- "Custom indicator"
  }
  
  # Subfunction that works only on a matrix
  get_one_result <- function(mat, ...) { 
    mat <- convert_to_matrix(mat)
    check_mat(mat)
    
    # Handle the arguments passed to the function and store them in the returned
    # object
    fun.args <- as.list(match.call(expand.dots = FALSE))[['...']]
    fun.args <- lapply(fun.args, eval, envir = parent.frame())
    
    # Sometimes the indicator does not take any extra arguments. In that case 
    # ... is NULL, so we need to replace that with an empty list
    if ( is.null(fun.args) ) { 
      fun.args <- list()
    }
    
    result <- list(value     = fun(mat, ...), 
                   orig_data = mat, 
                   fun.args  = fun.args, 
                   taskname = taskname, 
                   indicf = fun)
    
    class(result) <- c('custom_sews_single', 'simple_sews_single',
                       'sews_result_single')
    return(result)
  }
  
  # Actual function produced
  function(mat, ...) { 
    if ( is.list(mat) ) { 
      result <- future_lapply_seed(mat, get_one_result, ...)
      names(result) <- names(mat)
      class(result) <- c('custom_sews_list', 'simple_sews_list',
                         'sews_result_list')
    } else { 
      result <- get_one_result(mat, ...)
    }
    return(result)
  }
  
}


#
#' @rdname create_indicator
#'
#' @param mat A matrix or a list of matrices. 
#' 
#' @param ... Additional arguments to pass to the function \code{fun}
#' 
#'@export
compute_indicator <- function(mat, fun, 
                              taskname = as.character(substitute(fun)), 
                              ...) { 
  indicfun <- create_indicator(fun, taskname = taskname)
  indicfun(mat, ...)
}
