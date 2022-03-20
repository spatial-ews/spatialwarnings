//
// Function that takes care of coarse-graining
//
// Function to reduce original matrix to submatrix by averaging
// This function is used bt the generic indicators variance, skewness and moran
// correlation. It divides the square matrix into submatrices of mean cover.
//
// Written by Vishwesha Guttal, 7th Nov 2013.
// Modified by Sabiha Majumder
// Converted to c++ by Alexandre Génin
//

#include <RcppArmadillo.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix coarse_grain_cpp(NumericMatrix mat, int subsize) {
  
  // Integer division (round down to nearest integer). We convert to (int) 
  // as mat.n_rows may be a uword.
  int nr = (int)(mat.nrow() / subsize);
  int nc = (int)(mat.ncol() / subsize);
  
  NumericMatrix reduced_matrix = NumericMatrix(nr, nc);
  
  // Fill in values of the submatrix
  for ( int j=0; j<nc; j++ ) {
    for ( int i=0; i<nr; i++ ) {
      
      // Compute mean of the corresponding cells in the original matrix
      double sum = 0;
      for ( int l=(j*subsize); l < (j+1)*subsize; l++ ) {
        for ( int k=(i*subsize); k < (i+1)*subsize; k++ ) {
          sum += mat(k,l);
        }
      }
      
      reduced_matrix(i,j) = sum / ( subsize * subsize );
    }
  }
  
  return reduced_matrix;
}
