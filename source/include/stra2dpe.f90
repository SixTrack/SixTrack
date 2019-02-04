! start include/stra2dpe.f90
#ifndef TILT
  strack(i)=zero
  strackx(i)=ed(IX)
  strackz(i)=ek(IX)
#else
  strack(i)=zero
  strackx(i)=ed(IX)*tiltc(i)
  stracks(i)=ed(IX)*tilts(i)
  strackz(i)=ek(IX)*tiltc(i)
  strackc(i)=ek(IX)*tilts(i)
#endif
! end include/stra2dpe.f90
