#
# These tags are required to use functions in other packages
# 
#' @import ggplot2
#' @import plyr
#' @import stats
#' @import utils
# 
# These tags are required by Rcpp
# 
#' @useDynLib spatialwarnings
#' @importFrom Rcpp sourceCpp
# 
# 

#
# These tags are here to explicitely show what is imported
#
#' @import stats
#' @import utils

# Spatialwarnings will load the `future` package at boot, but does not use 
# directly its functions. Here we import plan() explicitely to make sure 
# check() does not complain that we depend on a package we don't use. 
#
#'@importFrom future plan 
future::plan
