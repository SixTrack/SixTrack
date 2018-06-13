#ifndef TILT
  strack(i)=smiv(1,i)
#else
  strack(i)=smiv(1,i)
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
