! start include/kickq08v.f90
#ifndef TILT
qu=(seven*ekk)*cikve
qv=(-seven*ekk)*crkve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(seven*ekk)*(tiltck*cikve-tiltsk*crkve)
qv=(-seven*ekk)*(tiltck*crkve+tiltsk*cikve)
#endif
! end include/kickq08v.f90
