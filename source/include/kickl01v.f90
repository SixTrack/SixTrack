! start include/kickl01v.f90
#ifndef TILT
dyy1=zero
dyy2=ekk
#else
dyy1=(-one*ekk)*tilts(k)                                         !hr08
dyy2=ekk*tiltc(k)
#endif
! end include/kickl01v.f90
