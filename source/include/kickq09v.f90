! start include/kickq09v.f90
#ifndef TILT
qu=(eight*ekk)*cikve
qv=(-eight*ekk)*crkve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(eight*ekk)*(tiltck*cikve-tiltsk*crkve)
qv=(-eight*ekk)*(tiltck*crkve+tiltsk*cikve)
#endif
! end include/kickq09v.f90
