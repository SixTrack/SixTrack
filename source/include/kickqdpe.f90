! start include/kickqdpe.f90
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(ed(IX)*tiltck)/(one+dpp)
qv=((-one*ed(IX))*tiltsk)/(one+dpp)
quz=((-one*ek(IX))*tiltck)/(one+dpp)
qvz=(ek(IX)*tiltsk)/(one+dpp)
! end include/kickqdpe.f90
