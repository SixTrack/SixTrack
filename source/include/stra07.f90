!start stra07.f90
#ifndef TILT
  strack(i)=smiv(1,i)*c1m15
#else
  strack(i)=smiv(1,i)*c1m15
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
!end stra07.f90
