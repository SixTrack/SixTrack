! start include/kickq09h.f90
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(eight*ekk)*(tiltck*crkve+tiltsk*cikve)
qv=(eight*ekk)*(tiltck*cikve-tiltsk*crkve)
! end include/kickq09h.f90
