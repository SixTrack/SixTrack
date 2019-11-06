! start include/multl04.f90
  qu=((dki(ix,2)/dki(ix,3))*dki(ix,2))/(one+dpp)
  dppi=(c1e3*dki(ix,2))/(one+dpp)
  t(1,2)=(t(1,2)+(qu*zl-dppi*dpp)*tilts(k))+dppi*tilts(k)
  t(1,4)=(t(1,4)+(dppi*dpp-qu*zl)*tiltc(k))-dppi*(one-tiltc(k))
! end include/multl04.f90
