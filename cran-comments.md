
Dear CRAN maintainers, 

This is a minor update to the spatialwarnings package with minor improvements, along 
with a few bug fixes. 

This release should fix the problem with leftovers 'assert()s' raised by Prof. Brian 
Ripley. Thanks for pointing those out, I thought they were stripped from CRAN package 
builds. 

Thanks, 

Alexandre Génin



## Test environments

This package was tested using the following environments: 

 * R-lib Github actions (macOS with R-latest, windows with R-latest, 
    ubuntu-latest with R-devel, ubuntu-latest with R-release, 
    ubuntu-latest with R-oldrel): 
     https://travis-ci.org/spatial-ews/spatialwarnings
 
 * Local linux computer (Arch Linux as of 2021-03-10, R 4.1.2)
 
 * M1 Mac using https://mac.r-project.org/macbuilder/submit.html (results at: 
     https://mac.r-project.org/macbuilder/results/1647797175-da4b0c2b39e71edf/)
 
## R CMD check results

One remaining NOTE occurred on some platforms: 

The package size is sometimes reported as exceeding 1Mb (on ubuntu platforms 
above), probably due to the use of Rcpp: 

* checking installed package size ... NOTE
  installed size is  7.9Mb
  sub-directories of 1Mb or more:
    libs   6.6Mb

## Changes in this release

This release has received external contributions from J. Guerber.

Improvements: 

  * In patch size distribution plots, fits are now rescaled when xmin is above 1 
      so they overlay nicely on the observed distribution
  
  * Better handling of errors in segmented::segmented() when fitting variograms 
  
Bug fixes: 
  
  * Fixed an error occurring when using xmin = "estimate" in patchdistr_sews
  
  * Remove calls to assert() that were leftover from testing (CRAN policy)
  
  
  
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
