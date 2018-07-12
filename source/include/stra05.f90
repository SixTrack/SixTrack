#ifndef TILT
  strack(i)=smiv(1,i)*c1m9
#else
  strack(i)=smiv(1,i)*c1m9
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
