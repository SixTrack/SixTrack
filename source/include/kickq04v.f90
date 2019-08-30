! start include/kickq04v.f90
#ifndef TILT
qu=(three*ekk)*cikve
qv=((-one*three)*ekk)*crkve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(three*ekk)*(tiltck*cikve-tiltsk*crkve)
qv=((-one*three)*ekk)*(tiltck*crkve+tiltsk*cikve)
#endif
! end include/kickq04v.f90
