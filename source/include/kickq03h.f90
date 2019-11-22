! start include/kickq03h.f90
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(ekk*two)*(tiltck*crkve+tiltsk*cikve)
qv=(ekk*two)*(tiltck*cikve-tiltsk*crkve)
! end include/kickq03h.f90
