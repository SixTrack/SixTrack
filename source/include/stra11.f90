! start include/stra11.f90
#ifndef TILT
  strack(i)=dki(ix,1)/dki(ix,3)
#else
  strack(i)=dki(ix,1)/dki(ix,3)
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
! end include/stra11.f90
