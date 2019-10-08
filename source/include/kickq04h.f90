! start include/kickq04h.f90
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(three*ekk)*(tiltck*crkve+tiltsk*cikve)
qv=(three*ekk)*(tiltck*cikve-tiltsk*crkve)
! end include/kickq04h.f90
