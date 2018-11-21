! start include/mul4v04.f90
#ifndef TILT
  yv2(j)=yv2(j)+strack(i)*(dpsv1(j)+omoidpsv(j)) ! hisix
#else
  yv1(j)=(yv1(j)-stracks(i)*(dpsv1(j)+omoidpsv(j)))+((c1e3*dki(ix,2))*moidpsv(j))*tilts(i)
  yv2(j)=(yv2(j)+strackc(i)*(dpsv1(j)+omoidpsv(j)))-((c1e3*dki(ix,2))*moidpsv(j))*(one-tiltc(i))
#endif
! end include/mul4v04.f90
