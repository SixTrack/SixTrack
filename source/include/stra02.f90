! start include/stra02.f90
#ifndef TILT
  strack(i)=smiv(i)
#else
  strack(i)=smiv(i)
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
! end include/stra02.f90
