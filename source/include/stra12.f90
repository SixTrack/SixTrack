#ifndef TILT
  strack(i)=dki(ix,1)
#else
  strack(i)=dki(ix,1)
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
