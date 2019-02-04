! start include/multu01.f90
#ifndef TILT
  qu=(((-one*dki(ix,1))/dki(ix,3))*dki(ix,1))/(one+dpp)        !hr03
  y(1,1)=(y(1,1)+qu*xl)-((dpp*c1e3)*dki(ix,1))/(one+dpp)       !hr03
#else
  qu=(((-one*dki(ix,1))/dki(ix,3))*dki(ix,1))/(one+dpp)        !hr03
  y(1,1)=(y(1,1)+(qu*xl-((dpp*c1e3)*dki(ix,1))/(one+dpp))*tiltc(k))+((c1e3*dki(ix,1))/(one+dpp))*(one-tiltc(k))
  y(1,2)=(y(1,2)+(qu*xl-((dpp*c1e3)*dki(ix,1))/(one+dpp))*tilts(k))+((c1e3*dki(ix,1))/(one+dpp))*tilts(k)
#endif
! end include/multu01.f90
