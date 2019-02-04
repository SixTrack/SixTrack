! start include/stra09.f90
#ifndef TILT
  strack(i)=smiv(i)*c1m21
#else
  strack(i)=smiv(i)*c1m21
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
! end include/stra09.f90
