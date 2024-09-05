// 
// 
// This file contains a function that will compute the r-spectrum
// 

#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma; 

#define SQ(a) ( (double)(a) * (double)(a) )

// Distance step of the r-spectrum
const double step = 1.0; 

// 
//' @title r-spectrum 
//' 
//' @description Compute the r-spectrum of a matrix 
//' 
//' @param mat A matrix with logical or numeric values 
//' 
//' @return A data.frame with two columns: \code{dist}, the wave number and 
//'   \code{rspec}, the normalized value of the r-spectrum
//' 
//' @details This functions returns a \code{data.frame} with \code{NA}s in the
//'   \code{rspec} column if the input matrix has zero variance. Note that if the matrix 
//'   is not square, then only the largest square matrix fitting in the upper left
//'   corner is used. 
//' 
//' @seealso \code{\link{spectral_sews}}
//' 
//' @examples 
//' 
//' # Spectrum of white noise
//' rmat <- matrix(runif(100*100) > .5, ncol = 100)
//' spec <- rspectrum(rmat) 
//' plot(spec, type = "l")
//' 
//' # Add some spatial correlation and compare the two spectra
//' rmat.cor <- rmat
//' for (i in seq(1, nrow(rmat)-1)) { 
//'   for (j in seq(1, nrow(rmat)-1)) { 
//'     rmat.cor[i,j] <- mean(rmat[(i-1):(i+1), (j-1):(j+1)])
//'   }
//' }
//' spec.cor <- rspectrum(rmat.cor)
//' plot(spec.cor, type = "n")
//' lines(spec, col = "black")
//' lines(spec.cor, col = "blue")
//' 
//' @export
// [[Rcpp::export]]
DataFrame rspectrum(arma::mat mat) { 
  
  // Get number of rows and number of cols to consider. If the matrix has an 
  // odd number of rows/cols, then we will not consider the last column/row so 
  // that the center falls on a cell. This is a minor approximation given the 
  // fact that we should work with matrices of large size in spatialwarnings 
  // (typically > 50 rows/cols)
  int nr = mat.n_rows - ( mat.n_rows % 2 == 0 ? 0 : 1); 
  int nc = mat.n_cols - ( mat.n_cols % 2 == 0 ? 0 : 1); 
  
  // Middle point of matrix. NOTE: this is exact only if nc and nr can be 
  // divided by 2: it is not when nc or nr is odd. This is why we discard the 
  // last row/col above if needed
  int n0x = nc / 2; 
  int n0y = nr / 2; 
  
  // Minimum and maximum distances to consider
  int mi = 1;
  int ma = 1 + (n0x < n0y ? n0x : n0y);
  
  // Initialize the variables that will hold the spectrum values
  arma::vec ray = arma::linspace(mi, ma, ma - mi + 1); // ma - mi + 1 values
  arma::vec rspectr = arma::zeros(ma - mi + 1); 
  
  // We check whether there is more than one value in the supplied matrix
  bool nonzero_variance = false;
  uword i = 1;
  while ( ! nonzero_variance && i < mat.n_elem ) {
    if ( mat(i - 1) != mat(i) ) { 
      nonzero_variance = true;
    }
    i++;
  } 
  
  // We have not found more than one value, or the matrix is too small -> 
  //   return early with NAs
  if ( ( ! nonzero_variance ) || nr < 1 || nc < 1 ) { 
    rspectr.fill(NumericVector::get_na()); 
    return DataFrame::create(_["dist"]  = ray, 
                             _["rspec"] = rspectr); 
  }
  
  // Compute and shift DFT
  arma::cx_mat mat_fft = arma::fft2(mat); 
  
  // Compute aspectr2D and normalize it
  // mat_fft(n0y, n0x) = 0; 
  
  // Compute r-spectrum
  double norm_factor = 0; 
  
  // For all lengths
  for ( uword l=0; l < ray.n_elem; l++ ) { 
    double r = ray(l);
    
    // Go through the matrix and make a sum of relevant fft coefficients
    int total_inmask = 0;
    
    // Note that we do not loop over cells that will not be in the correct 
    // distance range. 
    for ( int j = (n0y - r - 1); j < (n0y + r + 1); j++) { 
      for ( int i = (n0x - r - 1); i < (n0x + r + 1); i++) { 
        
        // Squared distance to center
        double dist = SQ(i - n0x) + SQ(j - n0y);
        
        if ( ( dist > 0 ) && 
             dist >= SQ(r - step/2) && 
             dist  < SQ(r + step/2) ) { 
          
          // We use shifted coordinates to pick the value in mat_fft 
          // instead of making a call to arma::shift() that does a copy of 
          // the matrix
          int shift_i = (i + n0x) % nc;
          int shift_j = (j + n0y) % nr;
          double aspectr2D_ij = SQ( abs(mat_fft(shift_j, shift_i)) ) / 
                                  SQ(SQ( (n0x + 1) * (n0y + 1) ));
          
          rspectr(l) += aspectr2D_ij;
          norm_factor += aspectr2D_ij;
          total_inmask++;
        }
        
      }
    }
    
//     if ( r < 3 ) { 
//       Rcout << "from: " << (r-step/2) << " to: " << (r+step/2) << " -> total_inmask: " << total_inmask << "\n"; 
//       Rcout << "n0x: " << n0x << " n0y: " << n0y << "\n"; 
//     }
    
    rspectr(l) = rspectr(l) / total_inmask;
  }
  
  rspectr = rspectr / norm_factor;
//   Rcout << ray << "\n"; 
//   Rcpp::Rcout << rspectr << "\n\n"; 
  return DataFrame::create(_["dist"]  = ray, 
                           _["rspec"] = rspectr);
}
