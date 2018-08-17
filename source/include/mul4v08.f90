! start include/mul4v08.f90
#ifndef TILT
  yv(1,j)=yv(1,j)+bbiv(1,1,i)*moidpsv(j) ! hisix
  yv(2,j)=yv(2,j)+aaiv(1,1,i)*moidpsv(j) ! hisix
#else
  yv(1,j)=yv(1,j)+(tiltc(i)*bbiv(1,1,i)-tilts(i)*aaiv(1,1,i))*moidpsv(j)
  yv(2,j)=yv(2,j)+(tiltc(i)*aaiv(1,1,i)+tilts(i)*bbiv(1,1,i))*moidpsv(j)
#endif
! end include/mul4v08.f90
