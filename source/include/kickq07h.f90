! start include/kickq07h.f90
#ifndef TILT
qu=(six*ekk)*crkve
qv=(six*ekk)*cikve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(six*ekk)*(tiltck*crkve+tiltsk*cikve)
qv=(six*ekk)*(tiltck*cikve-tiltsk*crkve)
#endif
! end include/kickq07h.f90
