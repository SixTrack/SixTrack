! start include/kickq02h.f90
#ifndef TILT
qu=ekk
qv=zero
#else
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=ekk*tiltck
qv=(-one*ekk)*tiltsk                                             !hr08
#endif
! end include/kickq02h.f90
