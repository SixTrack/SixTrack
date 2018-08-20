! start include/multu03.f90
#ifndef TILT
  y(1,1)=y(1,1)-((dki(ix,1)*dpp)/(one+dpp))*c1e3               !hr03
#else
  y(1,1)=(y(1,1)-(((dki(ix,1)*dpp)/(one+dpp))*c1e3)*tiltc(k))+((c1e3*dki(ix,1))/(one+dpp))*(one-tiltc(k))
  y(1,2)=(y(1,2)-(((dki(ix,1)*dpp)/(one+dpp))*c1e3)*tilts(k))+((c1e3*dki(ix,1))/(one+dpp))*tilts(k)
#endif
! end include/multu03.f90
