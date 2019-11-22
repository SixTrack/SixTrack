! start include/multl03.f90
  dppi=(c1e3*dki(ix,1))/(one+dpp)
  t(1,2)=(t(1,2)-(dppi*dpp)*tiltc(k))+dppi*(one-tiltc(k))
  t(1,4)=(t(1,4)-(dppi*dpp)*tilts(k))+dppi*tilts(k)
! end include/multl03.f90
