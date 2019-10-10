! start include/kickq05h.f90
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(four*ekk)*(tiltck*crkve+tiltsk*cikve)
qv=(four*ekk)*(tiltck*cikve-tiltsk*crkve)
! end include/kickq05h.f90
