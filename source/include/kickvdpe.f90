#ifndef TILT
  yv(1,j)=yv(1,j)+((strackx(i)*oidpsv(j))*crkve)*mtc(j)
  yv(2,j)=yv(2,j)-((strackz(i)*oidpsv(j))*cikve)*mtc(j)
#else
  yv(1,j)=yv(1,j)+(oidpsv(j)*(strackx(i)*crkve-stracks(i)*cikve))*mtc(j)
  yv(2,j)=yv(2,j)+(oidpsv(j)*(strackz(i)*cikve+strackc(i)*crkve))*mtc(j)
#endif
