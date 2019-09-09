! start include/kickq03v.f90
#ifndef TILT
qu=(ekk*two)*cikve
qv=((-one*ekk)*two)*crkve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(ekk*two)*(tiltck*cikve-tiltsk*crkve)
qv=((-one*ekk)*two)*(tiltck*crkve+tiltsk*cikve)
#endif
! end include/kickq03v.f90
