! start include/kickq03v.f90
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(ekk*two)*(tiltck*cikve-tiltsk*crkve)
qv=((-one*ekk)*two)*(tiltck*crkve+tiltsk*cikve)
! end include/kickq03v.f90
