! start include/kickv01v.f90
#ifndef TILT
yv2(j)=yv2(j)+(strack(i)*oidpsv(j))*mtc(j)
#else
yv1(j)=yv1(j)-(stracks(i)*oidpsv(j))*mtc(j)
yv2(j)=yv2(j)+(strackc(i)*oidpsv(j))*mtc(j)
#endif
! end include/kickv01v.f90
