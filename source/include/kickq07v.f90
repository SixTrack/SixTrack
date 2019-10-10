! start include/kickq07v.f90
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(six*ekk)*(tiltck*cikve-tiltsk*crkve)
qv=(-six*ekk)*(tiltck*crkve+tiltsk*cikve)
! end include/kickq07v.f90
