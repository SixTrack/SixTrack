! start include/kicka01h.f90
#ifndef TILT
mpe=20
dyy1=ekk
dyy2=zero
qu=zero
qv=zero
#else
mpe=20
dyy1=ekk*tiltc(k)
dyy2=ekk*tilts(k)
qu=zero
qv=zero
#endif
! end include/kicka01h.f90
