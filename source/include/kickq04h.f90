! start include/kickq04h.f90
#ifndef TILT
qu=(three*ekk)*crkve
qv=(three*ekk)*cikve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(three*ekk)*(tiltck*crkve+tiltsk*cikve)
qv=(three*ekk)*(tiltck*cikve-tiltsk*crkve)
#endif
! end include/kickq04h.f90
