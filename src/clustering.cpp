//
// Compute the clustering of a matrix, defined as the density of pairs over the density
// of a state (e.g. rho_aa / rho_a ).
//

#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif

#include <RcppArmadillo.h>
using namespace arma;

// Typedef to be compatible with windows
typedef unsigned short ushort;

//[[Rcpp::export]]
arma::Mat<arma::uword> clustering_core(arma::Mat<unsigned short> m,
                                       arma::uword nstates,
                                       bool wrap,
                                       bool use_8_nb) {

  arma::uvec dpairs(nstates);
  arma::uvec dsingle(nstates);
  dpairs.fill(0);
  dsingle.fill(0);

  uword nr = m.n_rows;
  uword nc = m.n_cols;

  double ncells = (double)nr * nc;

  //
  // Pair distribution for four neighbors (and wrapping, remove the edges if
  // no wrapping)
  //   |   |   |
  // - x - x - x (-) <-- the pairs on the right are the same than on the left
  //   |   |   | ( ) <--
  // - x - x - x (-) <--
  //   |   |   | ( ) <--
  // - x - x - x (-) <--
  //  (|) (|) (|) <- already counted at the top, it is the same pair
  //
  // if we don't wrap, we need to substract all the pairs at the edges we will never
  // consider: there are nrow + ncol
  //
  // Pair distribution for eight neighbors (and wrapping, remove the edges if
  // no wrapping)
  // \ | X | X | /
  // - x - x - x (-) <-- (already counted on the left side of the matrix)
  // X | X | X | (X) <-- (already counted on the left side of the matrix)
  // - x - x - x (-) <-- (already counted on the left side of the matrix)
  // X | X | X | (X) <-- (already counted on the left side of the matrix)
  // - x - x - x (-) <-- (already counted on the left side of the matrix)
  // /(| X | X |) \.
  //              ^-- already counted at the top
  //
  // if we don't wrap, we need to substract all the pairs at the edges and corners
  // we will never consider: there are 4 + nrow + ncol + (nrow-1)*2 + (ncol-1)*2

  // Note: we only considers neighbors that are up/left or upleft, because if we count
  // the ones on the other side we will count pairs twice
  for (uword j = 0; j < m.n_cols; j++) {
    for (uword i = 0; i < m.n_rows; i++) {

      const uword this_state = m(i, j);
      dsingle(this_state)++;

      if (wrap) {
        // left column
        ushort state_left = m(i, (nc + j - 1) % nc);
        ushort state_up = m((nr + i - 1) % nr, j);

        if (state_left == this_state) {
          dpairs(this_state)++;
        }

        if (state_up == this_state) {
          dpairs(this_state)++;
        }

        ushort state_upleft = m((nr + i - 1) % nr, (nc + j - 1) % nc);
        if (use_8_nb && state_upleft == this_state) {
          dpairs(this_state)++;
        }

        ushort state_upright = m((nr + i - 1) % nr, (nc + j - 1) % nc);
        if (use_8_nb && state_upright == this_state) {
          dpairs(this_state)++;
        }

      } else {
        if (j > 0) {
          ushort state_left = m(i, j - 1);
          if (state_left == this_state) {
            dpairs(this_state)++;
          }
        }

        if (i > 0) {
          ushort state_up = m(i - 1, j);
          if (state_up == this_state) {
            dpairs(this_state)++;
          }
        }

        if (use_8_nb && i > 0 && j > 0) {
          ushort state_upleft = m(i - 1, j - 1);
          if (state_upleft == this_state) {
            dpairs(this_state)++;
          }
        }

        if (use_8_nb && i > 0 && j < (nc - 1)) {
          ushort state_upright = m(i, j + 1);
          if (state_upright == this_state) {
            dpairs(this_state)++;
          }
        }
      }
    }
  }

  arma::Mat<uword> counts = arma::join_horiz(dpairs, dsingle);

  return (counts);
}

//[[Rcpp::export]]
arma::Mat<arma::uword> pair_counts_internal(arma::Mat<unsigned short> m,
                                            arma::uword nstates,
                                            bool wrap,
                                            bool use_8_nb) {

  arma::umat pairs(nstates, nstates);
  pairs.fill(0);

  uword nr = m.n_rows;
  uword nc = m.n_cols;

  double ncells = (double)nr * nc;

  //
  // Pair distribution for four neighbors (and wrapping, remove the edges if
  // no wrapping)
  //   |   |   |
  // - x - x - x (-) <-- the pairs on the right are the same than on the left
  //   |   |   | ( ) <--
  // - x - x - x (-) <--
  //   |   |   | ( ) <--
  // - x - x - x (-) <--
  //  (|) (|) (|) <- already counted at the top, it is the same pair
  //
  // if we don't wrap, we need to substract all the pairs at the edges we will never
  // consider: there are nrow + ncol
  //
  // Pair distribution for eight neighbors (and wrapping, remove the edges if
  // no wrapping)
  // \ | X | X | /
  // - x - x - x (-) <-- (already counted on the left side of the matrix)
  // X | X | X | (X) <-- (already counted on the left side of the matrix)
  // - x - x - x (-) <-- (already counted on the left side of the matrix)
  // X | X | X | (X) <-- (already counted on the left side of the matrix)
  // - x - x - x (-) <-- (already counted on the left side of the matrix)
  // /(| X | X |) \.
  //              ^-- already counted at the top
  //
  // if we don't wrap, we need to substract all the pairs at the edges and corners
  // we will never consider: there are 4 + nrow + ncol + (nrow-1)*2 + (ncol-1)*2

  // Note: we only considers neighbors that are up/left or upleft, because if we count
  // the ones on the other side we will count pairs twice
  for (uword j = 0; j < m.n_cols; j++) {
    for (uword i = 0; i < m.n_rows; i++) {

      const uword this_state = m(i, j);

      if (wrap) {
        // left column
        ushort state_left = m(i, (nc + j - 1) % nc);
        ushort state_up = m((nr + i - 1) % nr, j);

        pairs(this_state, state_left)++;
        pairs(this_state, state_up)++;

        if ( use_8_nb ) {
          ushort state_upleft = m((nr + i - 1) % nr, (nc + j - 1) % nc);
          ushort state_upright = m((nr + i - 1) % nr, (nc + j - 1) % nc);
          pairs(this_state, state_upleft)++;
          pairs(this_state, state_upright)++;
        }

      } else { // no wrapping
        if (j > 0) {
          ushort state_left = m(i, j - 1);
          pairs(this_state, state_left)++;
        }

        if (i > 0) {
          ushort state_up = m(i - 1, j);
          pairs(this_state, state_up)++;
        }

        if (use_8_nb && i > 0 && j > 0) {
          ushort state_upleft = m(i - 1, j - 1);
          pairs(this_state, state_upleft)++;
        }

        if (use_8_nb && i > 0 && j < (nc - 1)) {
          ushort state_upright = m(i, j + 1);
          pairs(this_state, state_upright)++;
        }
      }
    }
  }

  return(pairs);
}
