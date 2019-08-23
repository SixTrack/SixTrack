! start include/kickq03h.f90
#ifndef TILT
qu=(ekk*two)*crkve
qv=(ekk*two)*cikve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(ekk*two)*(tiltck*crkve+tiltsk*cikve)
qv=(ekk*two)*(tiltck*cikve-tiltsk*crkve)
#endif
! end include/kickq03h.f90
