! start include/alignu.f90
#ifndef TILT
  xl=x(1,1)-xs
  zl=x(1,2)-zs
  crkve=xl
  cikve=zl
#else
  xl=(x(1,1)-xs)*tiltc(k)+(x(1,2)-zs)*tilts(k)
  zl=(x(1,2)-zs)*tiltc(k)-(x(1,1)-xs)*tilts(k)                     !hr02
  crkve=xl
  cikve=zl
#endif
! end include/alignu.f90
