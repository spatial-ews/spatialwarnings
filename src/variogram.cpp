// 
// This file contains code that will compute the variogram of a matrix, 
// possibly not considering all cells for speed. 
// 

#include <math.h>

#include <Rcpp.h>
using namespace Rcpp; 

#define _DIST 0
#define _SQ 1
#define _NVAL 2

// How much we want to bias the short distances over long distances. A high
// value here will improve the sampling of variance at short distances, but 
// if the value is too high computing times may increase a lot. 
const double SHORT_DIST_BIAS = 3.0; 

//[[Rcpp::export]]
NumericMatrix variogram_internal_cpp(NumericMatrix mat, 
                                     int nmax, 
                                     int bins,
                                     double cutoff) { 
  
  int nr = mat.nrow(); 
  int nc = mat.ncol(); 
  
  // Compute the vector of mindist to maxdist and create variogram matrix. The 
  // latter has three columns in this order: the distance, the squared
  // differences and the number of values in the bin
  double mindist = 1; 
  double maxdist = sqrt( (double)(nr*nr + nc*nc) ); 
  maxdist = maxdist > cutoff ? cutoff : maxdist; 
  double maxrange = (maxdist - mindist);
  NumericMatrix vario(bins, 3); 
  for ( int i=0; i<bins; i++ ) { 
    // We do not consider distance = 0
    double dist = (i+1) * maxrange/bins; 
    vario(i, _DIST) = dist; 
  }
  
  // Cell 1
  int ntot = 0; 
  while ( ntot < nmax ) { 
    
    // Find a candidate point along with its distance
    bool found_coords = false; 
    int i1, j1, i2, j2; 
    double dist; 
    while ( ! found_coords ) { 
      // Get coordinates of cell in matrix
      i1 = (int) floor(R::runif(0, nr));
      j1 = (int) floor(R::runif(0, nc));
      i2 = (int) floor(R::runif(0, nr));
      j2 = (int) floor(R::runif(0, nc));
      
      // Here we use rejection sampling so that the sampled points for the 
      // variogram are biased towards zero distance: this helps having a 
      // variogram that is better-sampled at low distances, which are generally
      // the interesting parts for computing variogram-based EWS metrics. This is
      // quite inefficient and would deserve improvements.
      dist = sqrt( (double)( (i1-i2)*(i1-i2) + (j1-j2)*(j1-j2) ) ); 
      double p_keep = exp( - SHORT_DIST_BIAS * dist / maxrange ); 
//       Rcout << "dist: " << dist << " p_keep: " << p_keep << "\n"; 
      if ( R::runif(0, 1) < p_keep ) { 
        found_coords = true; 
      }
    }
//       Rcout << "cell1: " << cell1 << " i1: " << i1 << " j1: " << j1 << "\n";
    
    if ( dist < maxrange ) { 
      // Compute distance and squared difference
      double sq = pow(mat(i1, j1) - mat(i2, j2), 2); 
      int bin = (int) floor( bins * ( dist / maxrange ) );
//         Rcout << "dist: " << dist << " maxrange: " << maxrange << 
//           " bin: " << bin << "\n"; 
      vario(bin, _SQ) += sq; 
      vario(bin, _NVAL) += 1.0; 
      ntot++; 
    }
    
  }
  
  // Normalize the variogram
  for ( int i=0; i<bins; i++ ) { 
    vario(i, _SQ) = 0.5 * vario(i, _SQ) / vario(i, _NVAL);
  }
  
  return(vario); 
}
