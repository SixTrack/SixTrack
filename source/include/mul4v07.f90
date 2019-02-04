! start include/mul4v07.f90
#ifndef TILT
  yv1(j)=yv1(j)+yv1j*moidpsv(j)   ! hisix
  yv2(j)=yv2(j)+yv2j*moidpsv(j)   ! hisix
#else
  yv1(j)=yv1(j)+(tiltc(i)*yv1j-tilts(i)*yv2j)*moidpsv(j) ! hisix
  yv2(j)=yv2(j)+(tiltc(i)*yv2j+tilts(i)*yv1j)*moidpsv(j) ! hisix
#endif
! end include/mul4v07.f90
