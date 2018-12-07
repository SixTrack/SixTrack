! start include/mul4v08.f90
#ifndef TILT
  yv1(j)=yv1(j)+bbiv(1,i)*moidpsv(j) ! hisix
  yv2(j)=yv2(j)+aaiv(1,i)*moidpsv(j) ! hisix
#else
  yv1(j)=yv1(j)+(tiltc(i)*bbiv(1,i)-tilts(i)*aaiv(1,i))*moidpsv(j)
  yv2(j)=yv2(j)+(tiltc(i)*aaiv(1,i)+tilts(i)*bbiv(1,i))*moidpsv(j)
#endif
! end include/mul4v08.f90
