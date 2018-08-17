!start stra03.f90
#ifndef TILT
  strack(i)=smiv(1,i)*c1m3
#else
  strack(i)=smiv(1,i)*c1m3
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif

!end stra03.f90
