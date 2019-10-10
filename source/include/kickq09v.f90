! start include/kickq09v.f90
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(eight*ekk)*(tiltck*cikve-tiltsk*crkve)
qv=(-eight*ekk)*(tiltck*crkve+tiltsk*cikve)
! end include/kickq09v.f90
