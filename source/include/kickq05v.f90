! start include/kickq05v.f90
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(four*ekk)*(tiltck*cikve-tiltsk*crkve)
qv=((-one*four)*ekk)*(tiltck*crkve+tiltsk*cikve)
! end include/kickq05v.f90
