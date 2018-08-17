! start include/kickv01h.f90
#ifndef TILT
  yv(1,j)=yv(1,j)+(strack(i)*oidpsv(j))*mtc(j)
#else
  yv(1,j)=yv(1,j)+(strackc(i)*oidpsv(j))*mtc(j)
  yv(2,j)=yv(2,j)+(stracks(i)*oidpsv(j))*mtc(j)
#endif
! end include/kickv01h.f90
