! start include/multl04.f90
#ifndef TILT
  qu=((dki(ix,2)/dki(ix,3))*dki(ix,2))/(one+dpp)               !hr03
  dppi=(c1e3*dki(ix,2))/(one+dpp)                              !hr03
  t(1,4)=t(1,4)-qu*zl+dppi*dpp
#else
  qu=((dki(ix,2)/dki(ix,3))*dki(ix,2))/(one+dpp)               !hr03
  dppi=(c1e3*dki(ix,2))/(one+dpp)                              !hr03
  t(1,2)=(t(1,2)+(qu*zl-dppi*dpp)*tilts(k))+dppi*tilts(k)
  t(1,4)=(t(1,4)+(dppi*dpp-qu*zl)*tiltc(k))-dppi*(one-tiltc(k))
#endif
! end include/multl04.f90
