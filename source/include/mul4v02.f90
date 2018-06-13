#ifndef TILT
  yv(1,j)=yv(1,j)-strack(i)*(omoidpsv(j)+dpsv1(j))
#else
  yv(1,j)=(yv(1,j)-strackc(i)*(dpsv1(j) + omoidpsv(j)))+((c1e3*dki(ix,1))*moidpsv(j))*(one-tiltc(i))         ! hisix
  yv(2,j)=(yv(2,j)-stracks(i)*(dpsv1(j) + omoidpsv(j)))+((c1e3*dki(ix,1))*moidpsv(j))*tilts(i)               ! hisix
#endif
