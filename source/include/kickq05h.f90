! start include/kickq05h.f90
#ifndef TILT
qu=(four*ekk)*crkve
qv=(four*ekk)*cikve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(four*ekk)*(tiltck*crkve+tiltsk*cikve)
qv=(four*ekk)*(tiltck*cikve-tiltsk*crkve)
#endif
! end include/kickq05h.f90
