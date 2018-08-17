!start mul4v01.f90
#ifndef TILT
  yv(1,j)=yv(1,j)-((strack(i)*xlvj)*moidpsv(j)+(omoidpsv(j)+dpsv1(j)))*dki(ix,1)
#else
  yv(1,j)=(yv(1,j)-((((strack(i)*xlvj)*moidpsv(j))&
         +(omoidpsv(j)+dpsv1(j)))*dki(ix,1))*tiltc(i))&
         +(((c1e3*dki(ix,1))*moidpsv(j))*(one-tiltc(i)))
  yv(2,j)=(yv(2,j)-((((strack(i)*xlvj)*oidpsv(j))*mtc(j)&
         +omoidpsv(j)+dpsv1(j))*dki(ix,1))*tilts(i))&
         +(((c1e3*dki(ix,1))*oidpsv(j))*tilts(i))*mtc(j)
#endif
!end mul4v01.f90
