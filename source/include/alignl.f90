! start include/alignl.f90
#ifndef TILT
  xl=t(1,1)-xs
  zl=t(1,3)-zs
  crkve=xl
  cikve=zl
#else
  xl=(t(1,1)-xs)*tiltc(k)+(t(1,3)-zs)*tilts(k)
  zl=(t(1,3)-zs)*tiltc(k)-(t(1,1)-xs)*tilts(k)                    !hr02
  crkve=xl
  cikve=zl
#endif
! end include/alignl.f90
