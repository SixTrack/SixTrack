! start include/stra03.f90
#ifndef TILT
  strack(i)=smiv(i)*c1m3
#else
  strack(i)=smiv(i)*c1m3
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
! end include/stra03.f90
