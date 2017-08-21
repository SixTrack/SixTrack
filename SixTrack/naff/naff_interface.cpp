#include "NAFF.h"
#include <iostream>

extern "C" double tunenaff_(double* x,  double* xp, int* maxn, int* x_len, int* xp_len) {

  // Don't mix buffers with Fortran (make sure to flush before calling this code too)
  std::cout << std::flush;
  
  // For debugging of argument passing.
  // std::cout << "**TUNENAFF**" << std::endl << std::flush;
  // std::cout << "maxn   = " << *maxn   << std::endl << std::flush;
  // std::cout << "x_len  = " << *x_len  << std::endl;
  // std::cout << "xp_len = " << *xp_len << std::endl;
  // Flush before these memory accesses...
  // std::cout << "x[0]="  << x[0] << std::endl << std::flush;
  // std::cout << "xp[0]=" << xp[0] << std::endl << std::flush;
  // END for debugging

  // Input sanity checks
  if (maxn <= 0) {
    fprintf(stderr, "CRITICAL ERROR in double tunenaff_(...): maxn = %d <= 0", *maxn);
    exit(EXIT_FAILURE);
  }
  if ((*x_len > 0 and *xp_len > 0) and (*x_len < *maxn or *xp_len < *maxn)) {
    //In some cases, no x_len and xp_len is passed (they are set to 0);
    //then just hope maxn is OK. If not, check the lengths!
    fprintf(stderr, "CRITICAL ERROR in double tunenaff_(...): maxn is bigger than x_len or xp_len.");
    exit(EXIT_FAILURE);
  }
  if (*x_len != *xp_len) {
    fprintf(stderr, "CRITICAL ERROR in double tunenaff_(...): x_len is different than xp_len.");
    exit(EXIT_FAILURE);
  }

  //Copy the data from the FORTRAN arrays and into the vector that will be passed to NAFF
  std::vector<double> data;
  data.reserve(*maxn);
  std::vector<double> data_prime;
  data_prime.reserve(*maxn);
  for( int i = 0; i < *maxn; i++ ) {
    data.push_back(x[i]);
    data_prime.push_back(xp[i]);
  }

  //Call NAFF!
  NAFF naff;
  double tune = naff.get_f1(data,data_prime);

  // FFTW library returns tune from 0-0.5
  if (tune<0.1)
    tune = 1.0-tune;
  
  //More Debugging stuff..
  // std::cout << "tune = " << tune << std::endl;
  
  // Don't mix buffers with Fortran (make sure to flush before calling thus function)
  std::cout << std::flush;
  
  return tune;
}
