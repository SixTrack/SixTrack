! start include/mul4v02.f90
  yv1(j)=(yv1(j)-strackc(i)*(dpsv1(j) + omoidpsv(j)))+((c1e3*dki(ix,1))*moidpsv(j))*(one-tiltc(i))
  yv2(j)=(yv2(j)-stracks(i)*(dpsv1(j) + omoidpsv(j)))+((c1e3*dki(ix,1))*moidpsv(j))*tilts(i)
! end include/mul4v02.f90
