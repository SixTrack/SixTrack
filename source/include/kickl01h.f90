! start include/kickl01h.f90
#ifndef TILT
dyy1=ekk
dyy2=zero
#else
dyy1=ekk*tiltc(k)
dyy2=ekk*tilts(k)
#endif
! end include/kickl01h.f90
