! start include/kickvdpe.f90
#ifndef TILT
  yv1(j)=yv1(j)+((strackx(i)*oidpsv(j))*crkve)*mtc(j)
  yv2(j)=yv2(j)-((strackz(i)*oidpsv(j))*cikve)*mtc(j)
#else
  yv1(j)=yv1(j)+(oidpsv(j)*(strackx(i)*crkve-stracks(i)*cikve))*mtc(j)
  yv2(j)=yv2(j)+(oidpsv(j)*(strackz(i)*cikve+strackc(i)*crkve))*mtc(j)
#endif
! end include/kickvdpe.f90
