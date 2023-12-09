#' A list of binary matrices and their associated parameters 
#'
#' @rdname forestgap
#' 
#' @format A list of logical matrices which are the end results of simulations 
#'   from Kubo's Forest Gap model along a gradient of increasing values of 
#'   stress (see references). 
#' 
#' @source Generated using the implementation of Kubo's model in caspr 0.2.0 
#'   \url{https://github.com/fdschneider/caspr}. 
#' 
#' @references
#' 
#' Kubo, T., Iwasa, Y., & Furumoto, N. (1996). Forest spatial dynamics with gap
#'   expansion: Total gap area and gap size distribution. Journal of Theoretical
#'   Biology, 180(3), 229-246. \doi{10.1006/jtbi.1996.0099}
#' 
"forestgap"

#' @rdname forestgap
#' 
#' @format The parameters used for the simulations, as a data frame. 
#' 
#' @details Kubo's forest gap model has three parameters, \eqn{\alpha}{alpha} 
#'   that controls the reproductive rate of trees, \eqn{d}{d} controls the 
#'   non-spatialized mortality and \eqn{\delta}{delta} the increased mortality 
#'   due to the presence of a neighboring gap. 
#' 
"forestgap.pars"



#' @rdname arizona
#' 
#' @title Aerial views of vegetation from Arizona, USA
#' 
#' @format A list of logical matrices which were obtained through the 
#'   classification of aerial images of vegetation taken in Arizona (USA). 
#' 
#' @source Derived from the images provided in the Supplementary Material of 
#'   Rodriguez et al. (2017). 
#' 
#' @references
#' 
#' Rodriguez, F., A. G. Mayor, M. Rietkerk, and S. Bautista. 2017. A null model 
#'   for assessing the cover-independent role of bare soil connectivity as 
#'   indicator of dryland functioning and dynamics. Ecological Indicators.
#' 
"arizona"



#' @rdname serengeti
#' 
#' @title Serengeti dataset
#' 
#' @description Vegetation data along a rainfall gradient in Serengeti national 
#'   park. 
#' 
#' @details 
#' 
#'  The data-set consists of a rectangular area of size 7.5 km x 90 km. 
#'  These data are represented as a list of matrices. Each matrix is a moving 
#'  window of 7.5 km x 7.5 km which moves my 2.5 km along the length of the 
#'  rectangular data-set. 
#'   
#'  Each entry in the matrix is vegetation data at a resolution of 30m as 
#'  classified into binary units with FALSE (grass) and TRUE (forest).
#'  The rainfall data provided here is the average rainfall (mm/yr) of a 
#'  moving window of dimension 7.5km which moves my 2.5 km along the length 
#'  of the rectangular data-set.   
#' 
#' @format A list of logical matrices 
#' 
#' @source Extracted from Eby's et al (2017) supplementary material 
#'   \url{https://github.com/tee-lab/spacetime-csd/}
#' 
#' @references 
#' 
#' Eby, S., Agrawal, A., Majumder, S., Dobson, A.P. & Guttal, V. (2017). 
#' Alternative stable states and spatial indicators of critical slowing down 
#' along a spatial gradient in a savanna ecosystem: Global Ecology 
#' and Biogeography, 26, 638-649
#' 
#' Reed, D. N., Anderson, T. M., Dempewolf, J., Metzger, K., & Serneels, S. (2009). 
#' The spatial distribution of vegetation types in the Serengeti ecosystem: 
#' the influence of rainfall and topographic relief on vegetation patch 
#' characteristics. Journal of Biogeography, 36(4), 770-782.
# 
# The strings below have no effect except forcing R to read this file. 
"serengeti" 

#' @rdname serengeti
#' 
#' @format The annual rainfall corresponding to the matrices in the serengeti 
#'   dataset, in the corresponding order. 
#' 
"serengeti.rain"



# A list of binary matrices for the density-dependent aggregation model 
#' @rdname dda
#' 
#' @title Density-dependent aggregation model
#' 
#' @description \code{dda} is a list of matrices representing results from the 
#'   density-dependent aggregation model (Siteur et al. 2023) in \code{dda}.
#'   \code{dda.pars} is a data frame with the model parameters. Each row of
#'   \code{dda.pars} corresponds to a matrix in \code{dda}, in the same order. All
#'   parameters were maintained constant, except for tau (see model
#'   definition in Siteur et al. 2023).
#' 
#' @source Kindly provided by Koen Siteur
#' 
#' @examples 
#' 
#' ddasews <- lsw_sews(dda)
#' plot(ddasews, along = dda.pars[ ,"tau"]) # tau is the changing parameter
#' 
#' display_matrix(dda[[6]]) 
#' 
#' @references 
#' 
#'   Siteur, Koen, Quan-Xing Liu, Vivi Rottschäfer, Tjisse van der Heide, Max Rietkerk,
#'   Arjen Doelman, Christoffer Boström, and Johan van de Koppel. 2023. "Phase-Separation
#'    Physics Underlies New Theory for the Resilience of Patchy Ecosystems." Proceedings 
#'    of the National Academy of Sciences 120 (2): e2202683120. 
#'    https://doi.org/10.1073/pnas.2202683120.
"dda"

#' @rdname dda
#' 
"dda.pars"
