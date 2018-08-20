! start include/multl08.f90
#ifndef TILT
  t(6,2)=t(6,2)-(qu*xl+dppi)/(one+dpp)
#else
  t(6,2)=(t(6,2)-((qu*xl+dppi)/(one+dpp))*tiltc(k))-(dppi/(one+dpp))*(one-tiltc(k))
  t(6,4)=(t(6,4)-((qu*xl+dppi)/(one+dpp))*tilts(k))-(dppi/(one+dpp))*tilts(k)
#endif
! end include/multl08.f90
