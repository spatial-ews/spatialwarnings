// 
// A small function that computes moran's I index 
// 


#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif

#include <RcppArmadillo.h>
using namespace arma;

//' 
//' @title Spatial correlation at lag 1
//'
//' @description This function computes the Moran's I index of spatial 
//'   correlation at lag 1.
//' 
//' @param mat A matrix containing logical, numeric or integer values
//' 
//' @return The Moran's I numeric value as a numeric number.
//' 
//' @details This function returns the spatial correlation as measured by 
//'   the Moran's I index. If the variance of the matrix is zero, then 
//'   \code{NaN} is returned. This function assumes a 4-way neighborhood (a.k.a. 
//'   von-Neumann neighborhood), and does not wrap around at the sides of the matrix.
//' 
//' @seealso \code{\link{generic_sews}} 
//' 
//' @examples
//' 
//' # Spatial correlation of white noise is close to zero
//' rmat <- matrix(runif(1000) > .5, ncol = 100)
//' raw_moran(rmat) 
//' 
//' # Spatial correlation of a half-ones / half-zeros matrix is close to one. 
//' # This would produce close but inaccurate results in version <3.0.2
//' m <- cbind(matrix(1, nrow = 100, ncol = 50), 
//'            matrix(0, nrow = 100, ncol = 50))
//' 
//' raw_moran(m)
//' 
//'@export
//[[Rcpp::export]]
double raw_moran(arma::mat& mat) { 
  
  uword h = mat.n_rows;
  uword w = mat.n_cols;
  
  // We use vectorize so that mean and var do not do column-wise computations, 
  // or return the covariance matrix
  double m = mean( vectorise(mat) );
  
  // Handle center of matrix (i.e. everything minus first/last row and column)
  double moranI = 0;
  for ( uword j=1; j<(w-1); j++ ) { 
    for ( uword i=1; i<(h-1); i++ ) { 
      moranI += ( mat(i, j) - m ) * (   mat(i, j-1) + mat(i, j+1) 
                                      + mat(i-1, j) + mat(i+1, j) 
                                      - 4*m );
    }
  }
  
  // Handle edges 
  uword c = w - 1; // last column as index
  for ( uword i=1; i<(h-1); i++ ) { 
    // Left edge of matrix (except corners)
    moranI += ( mat(i, 0) - m ) * ( mat(i, 1) + mat(i-1, 0) + mat(i+1, 0) - 3*m );
    // Right edge of matrix (except corners)
    moranI += ( mat(i, c) - m ) * ( mat(i, c-1) + mat(i-1, c) + mat(i+1, c) - 3*m );
  }
  
  uword r = h - 1; // last line as index
  for ( uword j=1; j<(w-1); j++ ) { 
    // Top edge of matrix (except corners)
    moranI += ( mat(0, j) - m ) * ( mat(1,   j) + mat(0, j-1) + mat(0, j+1) - 3*m );
    // Bottom edge of matrix (except corners)
    moranI += ( mat(r, j) - m ) * ( mat(r-1, j) + mat(r, j-1) + mat(r, j+1) - 3*m );
  }
  
  // Handle corners 
  moranI += ( mat(0, 0) - m ) * ( mat(1,   0) + mat(0, 1)   - 2 * m); 
  moranI += ( mat(r, 0) - m ) * ( mat(r-1, 0) + mat(r, 1)   - 2 * m); 
  moranI += ( mat(r, c) - m ) * ( mat(r-1, c) + mat(r, c-1) - 2 * m); 
  moranI += ( mat(0, c) - m ) * ( mat(0, c-1) + mat(1, c)   - 2 * m); 
  
  // Compute 'variance'-like denominator
  double varsum = 0.0; 
  for ( uword i=0; i<h; i++ ) { 
    for ( uword j=0; j<w; j++ ) { 
      varsum += ( mat(i,j) - m ) * ( mat(i,j) - m );
    }
  }
  
  // Compute W, sum of all weights. 4*1 for each cell considered using 4-way 
  // neighborhood
  double W = ( 4.0 * (h - 2) * (w - 2) ) + // Center of matrix -> 4 nbs
               ( 3.0 * (h - 2) * 2 ) + // Left/right edges -> 3 nbs
               ( 3.0 * (w - 2) * 2 ) + // Top/bottom edges -> 3 nbs
               ( 2.0 * 4 ); // Four corners -> 2 nbs
  double N = h * w; 
  
  moranI = ( N / W ) * ( moranI / varsum );
  
  return( moranI );
}
