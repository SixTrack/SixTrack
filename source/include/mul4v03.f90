!start mul4v03.f90
#ifndef TILT
  yv(2,j)=yv(2,j)-((strack(i)*zlvj)*moidpsv(j)-dpsv1(j)-omoidpsv(j))*dki(ix,2)
#else
  yv(1,j)=(yv(1,j)+(((strack(i)*zlvj)*moidpsv(j)-dpsv1(j)-omoidpsv(j))*dki(ix,2))*tilts(i))&
         +((c1e3*dki(ix,2))*moidpsv(j))*tilts(i)
  yv(2,j)=(yv(2,j)-(((strack(i)*zlvj)*moidpsv(j)-dpsv1(j)-omoidpsv(j))*dki(ix,2))*tiltc(i))&
         -((c1e3*dki(ix,2))*moidpsv(j))*(one-tiltc(i))
#endif
!end mul4v03.f90
