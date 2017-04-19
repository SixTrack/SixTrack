/* Initialise the IA32/64 FPU flags from Fortran */
/* An init function which sets FPU flags when needed */

// Bring our own FPU_control.h, as Windows doesn't have it.
// Contains inline ASM macro...
#include "fpu_control.h"

void disable_xp_(void)
{
#ifdef __i386__
  /* Set FPU flags to use double, not double extended,
     with rounding to nearest */
  short unsigned int cw = (_FPU_DEFAULT & ~_FPU_EXTENDED)|_FPU_DOUBLE;
  _FPU_SETCW(cw);
#endif
}

