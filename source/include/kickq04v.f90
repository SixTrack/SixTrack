! start include/kickq04v.f90
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(three*ekk)*(tiltck*cikve-tiltsk*crkve)
qv=((-one*three)*ekk)*(tiltck*crkve+tiltsk*cikve)
! end include/kickq04v.f90
