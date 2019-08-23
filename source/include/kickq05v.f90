! start include/kickq05v.f90
#ifndef TILT
qu=(four*ekk)*cikve
qv=((-one*four)*ekk)*crkve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(four*ekk)*(tiltck*cikve-tiltsk*crkve)
qv=((-one*four)*ekk)*(tiltck*crkve+tiltsk*cikve)
#endif
! end include/kickq05v.f90
