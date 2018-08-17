!start kicka01v.f90
#ifndef TILT
  mpe=20
  dyy1=zero
  dyy2=ekk
  qu=zero
  qv=zero
#else
  mpe=20
  dyy1=(-one*ekk)*tilts(k)                                         !hr02
  dyy2=ekk*tiltc(k)
  qu=zero
  qv=zero
#endif
!end kicka01v.f90
