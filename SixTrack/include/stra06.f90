#ifndef TILT
  strack(i)=smiv(1,i)*c1m12
#else
  strack(i)=smiv(1,i)*c1m12
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
