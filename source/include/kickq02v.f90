! start include/kickq02v.f90
#ifndef TILT
qu=zero
qv=-one*ekk                                                      !hr02
#else
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=(-one*ekk)*tiltsk                                             !hr02
qv=(-one*ekk)*tiltck                                             !hr02
#endif
! end include/kickq02v.f90
