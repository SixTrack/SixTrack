! start include/kickq02v.f90
#ifndef TILT
qu=zero
qv=-one*ekk
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=(-one*ekk)*tiltsk
qv=(-one*ekk)*tiltck
#endif
! end include/kickq02v.f90
