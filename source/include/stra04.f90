! start include/stra04.f90
#ifndef TILT
  strack(i)=smiv(1,i)*c1m6
#else
  strack(i)=smiv(1,i)*c1m6
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
! end include/stra04.f90
