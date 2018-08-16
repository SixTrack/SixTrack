#ifndef TILT
  t(i,2)=t(i,2)+qu*t(i,1)
#else
  t(i,2)=t(i,2)+(qu*t(i,1))*tiltc(k)                         !hr08
  t(i,4)=t(i,4)+(qu*t(i,3))*tilts(k)                         !hr08
#endif
