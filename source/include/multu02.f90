! start include/multu02.f90
#ifndef TILT
  y(j,1)=y(j,1)+qu*x(j,1)
#else
  y(j,1)=y(j,1)+(qu*x(j,1))*tiltc(k)                         !hr03
  y(j,2)=y(j,2)+(qu*x(j,2))*tilts(k)                         !hr03
#endif
! end include/multu02.f90
