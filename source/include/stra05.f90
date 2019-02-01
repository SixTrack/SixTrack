! start include/stra05.f90
#ifndef TILT
  strack(i)=smiv(i)*c1m9
#else
  strack(i)=smiv(i)*c1m9
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
! end include/stra05.f90
