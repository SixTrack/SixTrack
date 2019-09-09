! start include/kickq06v.f90
#ifndef TILT
qu=(five*ekk)*cikve
qv=(-five*ekk)*crkve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(five*ekk)*(tiltck*cikve-tiltsk*crkve)
qv=(-five*ekk)*(tiltck*crkve+tiltsk*cikve)
#endif
! end include/kickq06v.f90
