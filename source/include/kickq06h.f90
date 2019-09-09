! start include/kickq06h.f90
#ifndef TILT
qu=(five*ekk)*crkve
qv=(five*ekk)*cikve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(five*ekk)*(tiltck*crkve+tiltsk*cikve)
qv=(five*ekk)*(tiltck*cikve-tiltsk*crkve)
#endif
! end include/kickq06h.f90
