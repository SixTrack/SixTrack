! start include/kickvxxv.f90
#ifndef TILT
yv1(j)=yv1(j)+((strack(i)*oidpsv(j))*cikve)*mtc(j)
yv2(j)=yv2(j)+((strack(i)*oidpsv(j))*crkve)*mtc(j)
#else
yv1(j)=yv1(j)+(oidpsv(j)*(strackc(i)*cikve-stracks(i)*crkve))*mtc(j)
yv2(j)=yv2(j)+(oidpsv(j)*(strackc(i)*crkve+stracks(i)*cikve))*mtc(j)
#endif
! end include/kickvxxv.f90
