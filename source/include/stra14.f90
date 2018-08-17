!start stra14.f90
#ifndef TILT
  strack(i)=dki(ix,2)
#else
  strack(i)=dki(ix,2)
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
!end stra14.f90
