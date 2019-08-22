! start include/kickadpe.f90
#ifndef TILT
dyy1=(ed(IX)*xl)/(one+dpp)
dyy2=(ek(IX)*zl)/(one+dpp)
mpe=20
qu=ed(IX)/(one+dpp)
quz=ek(IX)/(one+dpp)
qv=zero
qvz=zero
#else
dyy1=((ed(IX)*tiltc(k))*xl-(ek(IX)*tilts(k))*zl)/(one+dpp)
dyy2=((ek(IX)*tiltc(k))*zl+(ed(IX)*tilts(k))*xl)/(one+dpp)
mpe=20
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(ed(IX)*tiltck)/(one+dpp)
qv=((-one*ed(IX))*tiltsk)/(one+dpp)
quz=((-one*ek(IX))*tiltck)/(one+dpp)
qvz=(ek(IX)*tiltsk)/(one+dpp)
#endif
! end include/kickadpe.f90
