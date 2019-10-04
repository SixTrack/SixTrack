! start include/multl10.f90
  t(6,2)=(t(6,2)-((qu*zl+dppi)/(one+dpp))*tilts(k))-(dppi/(one+dpp))*tilts(k)
  t(6,4)=(t(6,4)+((qu*zl+dppi)/(one+dpp))*tiltc(k))+(dppi/(one+dpp))*(one-tiltc(k))
! end include/multl10.f90
