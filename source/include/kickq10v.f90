! start include/kickq10v.f90
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(nine*ekk)*(tiltck*cikve-tiltsk*crkve)
qv=(-nine*ekk)*(tiltck*crkve+tiltsk*cikve)
! end include/kickq10v.f90
