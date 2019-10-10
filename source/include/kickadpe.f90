! start include/kickadpe.f90
dyy1=((ed(IX)*tiltc(k))*xl-(ek(IX)*tilts(k))*zl)/(one+dpp)
dyy2=((ek(IX)*tiltc(k))*zl+(ed(IX)*tilts(k))*xl)/(one+dpp)
mpe=20
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(ed(IX)*tiltck)/(one+dpp)
qv=((-one*ed(IX))*tiltsk)/(one+dpp)
quz=((-one*ek(IX))*tiltck)/(one+dpp)
qvz=(ek(IX)*tiltsk)/(one+dpp)
! end include/kickadpe.f90
