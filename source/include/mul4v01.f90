! start include/mul4v01.f90
#ifndef TILT
  yv1(j)=yv1(j)-((strack(i)*xlvj)*moidpsv(j)+(omoidpsv(j)+dpsv1(j)))*dki(ix,1)
#else
  yv1(j)=(yv1(j)-((((strack(i)*xlvj)*moidpsv(j))&
         +(omoidpsv(j)+dpsv1(j)))*dki(ix,1))*tiltc(i))&
         +(((c1e3*dki(ix,1))*moidpsv(j))*(one-tiltc(i)))
  yv2(j)=(yv2(j)-((((strack(i)*xlvj)*oidpsv(j))*mtc(j)&
         +omoidpsv(j)+dpsv1(j))*dki(ix,1))*tilts(i))&
         +(((c1e3*dki(ix,1))*oidpsv(j))*tilts(i))*mtc(j)
#endif
! end include/mul4v01.f90
