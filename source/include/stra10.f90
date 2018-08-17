!start stra10.f90
#ifndef TILT
  strack(i)=smiv(1,i)*c1m24
#else
  strack(i)=smiv(1,i)*c1m24
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif

!end stra10.f90
