! start include/stra06.f90
#ifndef TILT
  strack(i)=smiv(i)*c1m12
#else
  strack(i)=smiv(i)*c1m12
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
! end include/stra06.f90
