! start include/kickq08h.f90
#ifndef TILT
qu=(seven*ekk)*crkve
qv=(seven*ekk)*cikve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(seven*ekk)*(tiltck*crkve+tiltsk*cikve)
qv=(seven*ekk)*(tiltck*cikve-tiltsk*crkve)
#endif
! end include/kickq08h.f90
