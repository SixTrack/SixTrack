! start include/kicku01h.f90
#ifndef TILT
y(1,1)=y(1,1)+ekk
#else
y(1,1)=y(1,1)+ekk*tiltc(k)
y(1,2)=y(1,2)+ekk*tilts(k)
#endif
! end include/kicku01h.f90
