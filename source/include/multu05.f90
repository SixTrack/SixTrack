! start include/multu05.f90
#ifndef TILT
  y(j,2)=y(j,2)-qu*x(j,2)
#else
  y(j,1)=y(j,1)+(qu*x(j,1))*tilts(k)                         !hr03
  y(j,2)=y(j,2)-(qu*x(j,2))*tiltc(k)                         !hr03
#endif
! end include/multu05.f90
