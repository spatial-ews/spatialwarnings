
Dear CRAN maintainers, 

This is an incremental update to the spatialwarnings package with a few improvements,
along with a few bug fixes. This package version should also fix NOTEs appearing on the CRAN test page (minor errors in the formatting of documentation).

Please note that my email has changed from <a.a.h.genin@uu.nl> to 
<alexandre.genin@sete.cnrs.fr>, as I have been changing institutions. I no longer 
have access to my last email. 

Thanks, 

Alexandre Génin


## Test environments

This package was tested using the following environments: 

 * R-lib Github actions (macOS with R-latest, windows with R-latest, 
    ubuntu-latest with R-devel, ubuntu-latest with R-release, 
    ubuntu-latest with R-oldrel): 
     https://github.com/spatial-ews/spatialwarnings/actions/runs/10738867095
 
 * Local linux computer (Arch Linux, R 4.5.0)
 
 * M1 Mac using https://mac.r-project.org/macbuilder/submit.html (R-release, results at: 
     https://mac.r-project.org/macbuilder/results/1725630635-ffef2c949e143726/)
 
 * Win builder service (R-release, devel, oldrelease)
 
 
## R CMD check results

One remaining NOTE occurred on some platforms: 

The package size is sometimes reported as exceeding 1Mb (on ubuntu platforms 
above), maybe due to the use of Rcpp: 

* checking installed package size ... NOTE
  installed size is  7.7Mb
  sub-directories of 1Mb or more:
    libs   6.4Mb


## Changes in this release

This is a maintenance release to fix an error arising in the package acss. acss is now 
explicitely attached to the search path, rather than just loaded. 

## Package Description

spatialwarnings is a package that assists ecologists in carrying out 
computations of early-warning signals (EWS) of ecosystem degradation.

These EWS are based on the fact that some ecosystems are expected to show 
specific spatial patterns before undergoing non-linear transitions (a wide 
shift in their state despite a small change in external forcings). For example, 
such ecosystems are expected to show an increase in spatial autocorrelation, 
variance and skewness, or, for patchy ecosystems, specific changes in the patch 
size distribution.

This packages assists users with computing these metrics efficiently on matrix 
objects in R, test their significance based on randomizing spatial structure, 
and plot their trends based on ggplot2. A convenient, three-step workflow is 
provided based on summary/plot/etc. generic functions.

Homepage and usage example:

  https://github.com/spatial-ews/spatialwarnings

Reference:
  
  * Génin, A. , Majumder, S. , Sankaran, S. , Danet, A. , Guttal, V. , 
    Schneider, F. D. and Kéfi, S. (2018),
    Monitoring ecosystem degradation using spatial data and the R package 
    'spatialwarnings'. Methods Ecol Evol. 
    doi:10.1111/2041-210X.13058
