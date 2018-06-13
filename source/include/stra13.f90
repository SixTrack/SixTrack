#ifndef TILT
  strack(i)=dki(ix,2)/dki(ix,3)
#else
  strack(i)=dki(ix,2)/dki(ix,3)
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
