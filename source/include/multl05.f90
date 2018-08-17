! start include/multl05.f90
#ifndef TILT
  t(i,4)=t(i,4)-qu*t(i,3)
#else
  t(i,2)=t(i,2)+(qu*t(i,1))*tilts(k)                         !hr03
  t(i,4)=t(i,4)-(qu*t(i,3))*tiltc(k)                         !hr03
#endif
! end include/multl05.f90
