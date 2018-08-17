!start mul4v04.f90
#ifndef TILT
  yv(2,j)=yv(2,j)+strack(i)*(dpsv1(j)+omoidpsv(j)) ! hisix
#else
  yv(1,j)=(yv(1,j)-stracks(i)*(dpsv1(j)+omoidpsv(j)))+((c1e3*dki(ix,2))*moidpsv(j))*tilts(i)
  yv(2,j)=(yv(2,j)+strackc(i)*(dpsv1(j)+omoidpsv(j)))-((c1e3*dki(ix,2))*moidpsv(j))*(one-tiltc(i))
#endif
!end mul4v04.f90
