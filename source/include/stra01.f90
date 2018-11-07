! start include/stra01.f90
#ifndef TILT
  strack(i)=smiv(i)*c1e3
#else
  strack(i)=smiv(i)*c1e3
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
! end include/stra01.f90
