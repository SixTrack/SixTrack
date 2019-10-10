! start include/kickq08v.f90
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(seven*ekk)*(tiltck*cikve-tiltsk*crkve)
qv=(-seven*ekk)*(tiltck*crkve+tiltsk*cikve)
! end include/kickq08v.f90
