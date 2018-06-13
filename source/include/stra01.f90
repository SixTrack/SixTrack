#ifndef TILT
  strack(i)=smiv(1,i)*c1e3
#else
  strack(i)=smiv(1,i)*c1e3
  strackc(i)=strack(i)*tiltc(i)
  stracks(i)=strack(i)*tilts(i)
#endif
