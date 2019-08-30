! start include/kickq07v.f90
#ifndef TILT
qu=(six*ekk)*cikve
qv=(-six*ekk)*crkve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(six*ekk)*(tiltck*cikve-tiltsk*crkve)
qv=(-six*ekk)*(tiltck*crkve+tiltsk*cikve)
#endif
! end include/kickq07v.f90
