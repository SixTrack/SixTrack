! start include/kickq09h.f90
#ifndef TILT
qu=(eight*ekk)*crkve
qv=(eight*ekk)*cikve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(eight*ekk)*(tiltck*crkve+tiltsk*cikve)
qv=(eight*ekk)*(tiltck*cikve-tiltsk*crkve)
#endif
! end include/kickq09h.f90
