// 
// A function that does the patch labelling
// 

#include <Rcpp.h>
#include <queue>
#include <utility>
#include "./headers.h"

using namespace Rcpp;

// Default value for non-patch cells
#define DEFAULT_VALUE NA_INTEGER

// These are the number rows and columns of the neighbor mask. They are assumed
// to be 3 (only considering the directly-adjacent neighbors) because it speeds
// up computations. 
const int maskrow = 3; 
const int maskcol = 3; 

// [[Rcpp::export]]
IntegerMatrix label_cpp(IntegerMatrix mat, 
                        IntegerMatrix nbmask,
                        bool wrap) { 
  
  int W = mat.ncol();
  int H = mat.nrow();
  
  // Initialize variables 
  IntegerMatrix output(H,W);
  LogicalMatrix is_marked(H,W);
  int patch_id = 1;
  
  std::vector<int> patchsizes;
  IntegerVector patch;
  bool percolation = false;
  
  for (int i=0; i<H; i++) { 
    for (int j=0; j<W; j++) { 
      // We consider the cell (i,j).
      
      // Default value is NA if not in a patch
      if ( ! (mat(i,j) > 0) ) { 
        output(i,j) = DEFAULT_VALUE;
        is_marked(i,j) = 1;
      }
      
      // If it is within a patch and not marked already 
      if ( ! is_marked(i, j) ) { 
        // We flood fill the patch
        std::pair<int,int> xy = std::make_pair(i, j);
        patch = flood_fill(mat, is_marked, output, nbmask, 
                           xy, patch_id, wrap);
        // Update percolation status, patch sizes and increase patch id
        percolation = percolation ? true : patch(0) == 1;
        patchsizes.push_back(patch(1));
        patch_id++;
      }
      output(i,j);
    }
  }
  
  output.attr("psd") = patchsizes;
  output.attr("percolation") = percolation;
  return(output);
}


// This flood fill is implemented with a queue
IntegerVector flood_fill(const IntegerMatrix &mat, 
                         LogicalMatrix &is_marked,
                         IntegerMatrix &output,
                         const IntegerMatrix nbmask,
                         std::pair<int, int> xy,
                         int fillcol,
                         bool wrap) { 
  
  // Create empty queue
  std::queue <std::pair<int, int> > to_fill;
  
  int W = mat.ncol();
  int H = mat.nrow();
  
  // Keep a count of the number of cell counted (== patch size)
  int patch_size = 0;
  
  // Keep a count of the left/right/top/bottom coordinates travelled, 
  int xmin = xy.first;
  int xmax = xy.first;
  int ymin = xy.second;
  int ymax = xy.second;
  int percolation = 0;
  
  // Add our cell to the queue
  to_fill.push(xy);
  is_marked(xy.first, xy.second) = 1; // The cell has been pushed to the queue
  patch_size++;
  
  while ( ! to_fill.empty() ) { 
    xy = to_fill.front();
    to_fill.pop(); 
    int i = xy.first;
    int j = xy.second;
    
    // We paint the pixel in our fill color (it always starts with vegetation)
    output(i, j) = fillcol;
    
    // We update the max coordinates of the patch
    xmin = i < xmin ? i : xmin;
    xmax = i > xmax ? i : xmax;
    ymin = j < ymin ? j : ymin;
    ymax = j > ymax ? j : ymax;
    
    // We consider its neighbors
    std::queue <std::pair<int, int> > nb = get_nb_coords(W, H, xy, nbmask, wrap);
    
    // We add the neighbors to the queue if needed
    while ( ! nb.empty() ) { 
      std::pair<int, int> cur_nb = nb.front(); 
      nb.pop(); 
      int newx = cur_nb.first; 
      int newy = cur_nb.second; 
      if ( !is_marked(newx, newy) && (mat(newx, newy) > 0) ) { 
        // Add the cell to the queue of things to fill
        to_fill.push(cur_nb);
        // We mark this cell to know it has been put in the queue and update
        //   the total number of cells marked (= patch size)
        is_marked(newx, newy) = true;
        patch_size++;
      }
    }
  }
    
  if ( (xmax - xmin + 1) == mat.nrow() || (ymax - ymin + 1) == mat.ncol() ) { 
    percolation = 1;
  }
  
  return( IntegerVector::create(percolation, patch_size) );
}

// A function that returns the coordinates of neighboring cells in a matrix,
//   taking into account the wraparound
std::queue <std::pair<int, int> > get_nb_coords(const int W, // width
                                                const int H, // height
                                                const std::pair<int,int> X,
                                                const IntegerMatrix& nbmask, 
                                                const bool wrap) { 
  
  std::queue <std::pair<int, int> > neighbors_xy; 
  
  for (int nbi=0; nbi<maskrow; nbi++) { 
    for (int nbj=0; nbj<maskcol; nbj++) { 
      // If the neighbor is to be considered, record the necessary shifts
      if ( nbmask(nbi, nbj) > 0 ) { 
        int shift_x = nbi - (maskrow-1)/2;
        int shift_y = nbj - (maskcol-1)/2;
        
        // Actual coordinates of the neighbor
        int nb_x = X.first  + shift_x;
        int nb_y = X.second + shift_y;
        
        // Does the neighbor fall outside the matrix ? 
        bool is_out = (nb_x < 0) || (nb_x >= H) || (nb_y < 0) || (nb_y >= W);
        
        // It is out and we don't wrap around: do not count this neighbor
        if ( (!wrap) && is_out ) { 
          // Nothing
        } else { // The neighbors falls within the field (or we wrap around)
          neighbors_xy.push( std::make_pair((nb_x + H) % H, 
                                            (nb_y + W) % W) );
        }
      }
    }
  }
  
  return neighbors_xy; 
}

