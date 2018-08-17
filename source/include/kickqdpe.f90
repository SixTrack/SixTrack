! start include/kickqdpe.f90
#ifndef TILT
qu=ed(IX)/(one+dpp)
quz=ek(IX)/(one+dpp)
qv=zero
qvz=zero
#else
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=(ed(IX)*tiltck)/(one+dpp)                                     !hr02
qv=((-one*ed(IX))*tiltsk)/(one+dpp)                              !hr02
quz=((-one*ek(IX))*tiltck)/(one+dpp)                             !hr02
qvz=(ek(IX)*tiltsk)/(one+dpp)                                    !hr02
#endif
! end include/kickqdpe.f90
