! start include/kickq02h.f90
#ifndef TILT
qu=ekk
qv=zero
#else
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=ekk*tiltck
qv=(-one*ekk)*tiltsk
#endif
! end include/kickq02h.f90
