! start include/kicku01v.f90
#ifndef TILT
y(1,2)=y(1,2)+ekk
#else
y(1,1)=y(1,1)-ekk*tilts(k)
y(1,2)=y(1,2)+ekk*tiltc(k)
#endif
! end include/kicku01v.f90
