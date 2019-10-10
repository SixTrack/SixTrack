! start include/kickq06v.f90
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(five*ekk)*(tiltck*cikve-tiltsk*crkve)
qv=(-five*ekk)*(tiltck*crkve+tiltsk*cikve)
! end include/kickq06v.f90
