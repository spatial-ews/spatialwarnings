// 
// 
// This file contains a function that computes the normalizing constant of 
//   a TPL. The function is vectorized. 
// 

#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif 

#include <RcppArmadillo.h>

using namespace arma; 

//[[Rcpp::export]]
arma::vec tplsum(double expo, double rate, arma::ivec xs, int xmin) { 
  
  arma::vec output(xs.n_elem);
  
  for ( uword i=0; i<xs.n_elem; i++ ) { 
    int x = xs(i);
    double total = 0;
    for (int k=xmin; k<x; k++) { 
      total += pow(k, -expo) * exp(-k*rate);
    }
    output(i) = total;
  }
  
  return(output);
}


// This function will carry out the above computation until "infinity", i.e. 
// when the relative change in the sum value goes below a tolerance threshold
// or reaches a maximum number of iterations. The constant returned is used 
// as a normalization term when computing the probabilities of truncated 
// power laws. 
//[[Rcpp::export]]
double tplinfsum(double expo, 
                 double rate, 
                 double xmin, 
                 arma::uword maxit, 
                 double reltol) { 
  
  double current_term = pow(xmin, -expo) * exp(-xmin*rate);
  double total = current_term;
  double rel_change = 1.0; 
  // We use a double for k just in case we go out of the integer range (this is something
  // R complains about sometimes in kmax)
  double k = xmin + 1.0; 
  uword k_stop = k + maxit; 
  // Rcpp::Rcout << "ct: " << current_term << "\n"; 
  
  while ( reltol < rel_change && k <= k_stop ) { 
  // while ( k <= kmax ) { 
    current_term = pow(k, - expo) * exp(- k * rate);
    rel_change = current_term / total; 
    total += current_term; 
    // if ( (uword) k % (uword) 1e6 == 0 || k > (kmax - 3) ) { 
      // Rcpp::Rcout << "k : " << k << " ct: " << current_term << " relc: " << 
        // rel_change << " total: " << total << "\n"; 
    // }
    k += 1.0; 
  }
  
  // Emit warning if we hit k_stop
  if ( k == k_stop ) { 
    Rcpp::Function warning("warning"); 
    warning("Maximum number of iterations reached in tplinfsum, increase cap with options(spatialwarnings.constants.maxit = <x>"); 
  }
  
  return(total);
}

