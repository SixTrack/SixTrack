! start include/kickq06h.f90
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(five*ekk)*(tiltck*crkve+tiltsk*cikve)
qv=(five*ekk)*(tiltck*cikve-tiltsk*crkve)
! end include/kickq06h.f90
