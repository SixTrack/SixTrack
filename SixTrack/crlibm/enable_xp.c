/* Initialise the IA32/64 FPU flags from Fortran */
/* An init function which sets FPU flags when needed */

#ifdef __APPLE__

void enable_xp_(void)
{
}

#else

#include <fpu_control.h>

void enable_xp_(void)
{
  /* Set FPU flags to use not double but extended,
     with rounding to nearest */
  short unsigned int cw = (_FPU_DEFAULT & ~_FPU_DOUBLE)|_FPU_EXTENDED;
  _FPU_SETCW(cw);
}

#endif
