#include "NAFF.h"
#include <iostream>

// C++ <-> Fortran interface for NAFF
// K. Sjobak and S. Kostoglou (CERN)
// August 2017

extern "C" double tunenaff(double** x,  double** xp, int* maxn, int* plane_idx, int* norm_flag) {
  // NOTE: double** x and double** xp are actually pointers to Fortran array descriptors,
  //  and the first 4 or 8 bytes (32- or 64-bit) of these contain the base address of the array.
  //  Therefore, using (*x)[idx] is basically accessing the array as a normal 1D array.
  //
  // The complication is that GFORTRAN and IFORT do not use the same types of array descriptors;
  //  for other compilers I don't know what the situation is.
  //
  // More information:
  //  https://software.intel.com/en-us/node/691959
  //  https://gcc.gnu.org/wiki/ArrayDescriptorUpdate
  //
  // Note that Fortran will pass these for
  //  (1) Allocatable or assumed-shape array
  //  (2) Fortran array pointers
  //  (3) Coarrays
  //  (4) Class objects
  // in the case of an explicit INTERFACE block.
  //
  // A suitable interface block for this function is:
  //    interface
  //       REAL(C_DOUBLE) function tunenaff
  //   &        (x,xp,maxn,plane_idx,norm_flag) BIND(C)
  //       use, intrinsic :: ISO_C_BINDING
  //       IMPLICIT NONE
  //       REAL(C_DOUBLE), dimension(:) :: x,xp
  //       INTEGER(C_INT) :: maxn, plane_idx, norm_flag
  //       end function
  //    end interface
  //
  // It may be more portable to convert the array to a C-style array pointer
  //  using C_LOC() from ISO_C_BINDING, and then pass that.
  // See also: https://stackoverflow.com/a/11935949/6603597
  
  // Don't mix buffers with Fortran (make sure to flush before calling this code too)
  std::cout << std::flush;
  
  // For debugging of argument passing.
  
  /*std::cout << "**TUNENAFF**"                 << std::endl << std::flush;
  std::cout << "maxn      = " << *maxn        << std::endl << std::flush;
  std::cout << "plane_idx = " << *plane_idx   << std::endl << std::flush;
  std::cout << "norm_flag = " << *norm_flag   << std::endl << std::flush;
  //std::cout << "x_len     = " << *x_len       << std::endl << std::flush;
  //std::cout << "xp_len    = " << *xp_len      << std::endl << std::flush;
  std::cout << "x[0]      = " << **x          << std::endl << std::flush;
  std::cout << "xp[0]     = " << **xp         << std::endl << std::flush;
  for( int i = 0; i < *maxn; i++ ) {
    std::cout << "i=" << i << std::endl << std::flush;
    std::cout << (*x)[i] << " " << (*xp)[i] << std::endl << std::flush;
  }*/
  
  // END debugging of argument passing
  
  // Input sanity checks
  if (*maxn <= 0) {
    fprintf(stderr, "CRITICAL ERROR in double tunenaff_(...): maxn = %d <= 0", *maxn);
    exit(EXIT_FAILURE);
  }
  
  //Copy the data from the FORTRAN arrays and into the vector that will be passed to NAFF
  std::vector<double> data;
  data.reserve(*maxn);
  std::vector<double> data_prime;
  data_prime.reserve(*maxn);

  // For physical coordinates only the real x array is used as an input, while for normalized coordinates the complex (x,px) array is passed to the NAFF.
  if ( (*norm_flag) == 0) {
    for( int i = 0; i < *maxn; i++ ) {
      //std::cout << i << " " << (*x)[i] << " " << (*xp)[i] << std::endl;
      data.push_back((*x)[i]);
      data_prime.push_back(0.0);
    }
  }
  else {
    for( int i = 0; i < *maxn; i++ ) {
      //std::cout << i << " " << (*x)[i] << " " << (*xp)[i] << std::endl;
      data.push_back((*x)[i]);
      data_prime.push_back((*xp)[i]);
    }
  }

  //Call NAFF!
  NAFF naff;

  // Set window to Hann 2nd order for transverse planes, by default Hann 1st order for longitudinal motion 
  if ( (*plane_idx) != 3 )
    naff.set_window_parameter(2, 'h');

  double tune = naff.get_f1(data,data_prime);

  // FFTW library returns tune from 0-0.5
  if ( (*plane_idx) == 3 )
    tune = 1.0-tune;
  
  //More Debugging stuff..
  /*std::cout << "tune = " << tune << std::endl;
  std::cout << "Window parameter: " << naff.get_window_parameter()<< std::endl;*/
  
  // Don't mix buffers with Fortran (make sure to flush before calling thus function)
  std::cout << std::flush;
  
  return tune;
}
