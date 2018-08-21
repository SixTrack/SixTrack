#include "NAFF.h"
#include <iostream>

// C++ <-> Fortran interface for NAFF
// K. Sjobak and S. Kostoglou (CERN)
// August 2017

extern "C" double tunenaff(double* x,  double* xp, int maxn, int plane_idx, int norm_flag, double fft_naff) {
  
  // Don't mix buffers with Fortran (make sure to flush before calling this code too)
  std::cout << std::flush;
  
  // For debugging of argument passing.
  /*
  std::cout << "**TUNENAFF**"                 << std::endl << std::flush;
  std::cout << "maxn      = " << maxn        << std::endl << std::flush;
  std::cout << "plane_idx = " << plane_idx   << std::endl << std::flush;
  std::cout << "norm_flag = " << norm_flag   << std::endl << std::flush;
  std::cout << "x[0]      = " << *x          << std::endl << std::flush;
  std::cout << "xp[0]     = " << *xp         << std::endl << std::flush;
  std::cout << "fft_naff     = " << fft_naff         << std::endl << std::flush;*/
  /*
  for( int i = 0; i < maxn; i++ ) {
    std::cout << "i=" << i << std::endl << std::flush;
    std::cout << x[i] << " " << xp[i] << std::endl << std::flush;
  }
  */
  // END debugging of argument passing
  
  // Input sanity checks
  if (maxn <= 0) {
    fprintf(stderr, "CRITICAL ERROR in double tunenaff_(...): maxn = %d <= 0", maxn);
    exit(EXIT_FAILURE);
  }
  
  //Copy the data from the FORTRAN arrays and into the vector that will be passed to NAFF
  std::vector<double> data;
  data.reserve(maxn);
  std::vector<double> data_prime;
  data_prime.reserve(maxn);

  // For physical coordinates only the real x array is used as an input, while for normalized coordinates the complex (x,px) array is passed to the NAFF.
  if ( norm_flag == 0) {
    for( int i = 0; i < maxn; i++ ) {
      //std::cout << i << " " << x[i] << " " << xp[i] << std::endl;
      data.push_back(x[i]);
      data_prime.push_back(0.0);
    }
  }
  else {
    for( int i = 0; i < maxn; i++ ) {
      //std::cout << i << " " << x[i] << " " << xp[i] << std::endl;
      data.push_back(x[i]);
      data_prime.push_back(xp[i]);
    }
  }

  //Call NAFF!
  NAFF naff;

  // Set window to Hann 2nd order for transverse planes, by default Hann 1st order for longitudinal motion 
  if ( plane_idx != 3 )
    naff.set_window_parameter(2, 'h');

  double tune = naff.get_f1(data,data_prime,fft_naff);

  // FFTW library returns tune from 0-0.5
  if ( plane_idx == 3 )
    tune = 1.0-tune;
  
  //More Debugging stuff..
  /*
  std::cout << "tune = " << tune << std::endl;
  std::cout << "Window parameter: " << naff.get_window_parameter()<< std::endl;
  */
  
  // Don't mix buffers with Fortran (make sure to flush before calling thus function)
  std::cout << std::flush;
  
  return tune;
}
