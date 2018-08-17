!start kickvxxv.f90
#ifndef TILT
yv(1,j)=yv(1,j)+((strack(i)*oidpsv(j))*cikve)*mtc(j)
yv(2,j)=yv(2,j)+((strack(i)*oidpsv(j))*crkve)*mtc(j)
#else
yv(1,j)=yv(1,j)+(oidpsv(j)*(strackc(i)*cikve-stracks(i)*crkve))*mtc(j)
yv(2,j)=yv(2,j)+(oidpsv(j)*(strackc(i)*crkve+stracks(i)*cikve))*mtc(j)
#endif
!end kickvxxv.f90
