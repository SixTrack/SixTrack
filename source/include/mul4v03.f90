! start include/mul4v03.f90
  yv1(j)=(yv1(j)+(((strack(i)*zlvj)*moidpsv(j)-dpsv1(j)-omoidpsv(j))*dki(ix,2))*tilts(i))&
         +((c1e3*dki(ix,2))*moidpsv(j))*tilts(i)
  yv2(j)=(yv2(j)-(((strack(i)*zlvj)*moidpsv(j)-dpsv1(j)-omoidpsv(j))*dki(ix,2))*tiltc(i))&
         -((c1e3*dki(ix,2))*moidpsv(j))*(one-tiltc(i))
! end include/mul4v03.f90
