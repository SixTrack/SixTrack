! start include/kickq10h.f90
#ifndef TILT
qu=(nine*ekk)*crkve
qv=(nine*ekk)*cikve
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(nine*ekk)*(tiltck*crkve+tiltsk*cikve)
qv=(nine*ekk)*(tiltck*cikve-tiltsk*crkve)
#endif
! end include/kickq10h.f90
