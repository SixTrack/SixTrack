! start include/multl06.f90
#ifndef TILT
  dppi=(c1e3*dki(ix,2))/(one+dpp)                              !hr03
  t(1,4)=t(1,4)+dppi*dpp
#else
  dppi=(c1e3*dki(ix,2))/(one+dpp)                              !hr03
  t(1,2)=(t(1,2)-(dppi*dpp)*tilts(k))+dppi*tilts(k)
  t(1,4)=(t(1,4)+(dppi*dpp)*tiltc(k))-dppi*(one-tiltc(k))
#endif
! end include/multl06.f90
