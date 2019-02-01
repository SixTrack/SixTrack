! start include/kickadpe.f90
#ifndef TILT
dyy1=(ed(IX)*xl)/(one+dpp)                                       !hr02
dyy2=(ek(IX)*zl)/(one+dpp)                                       !hr02
mpe=20
qu=ed(IX)/(one+dpp)
quz=ek(IX)/(one+dpp)
qv=zero
qvz=zero
#else
dyy1=((ed(IX)*tiltc(k))*xl-(ek(IX)*tilts(k))*zl)/(one+dpp)       !hr02
dyy2=((ek(IX)*tiltc(k))*zl+(ed(IX)*tilts(k))*xl)/(one+dpp)       !hr02
mpe=20
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=(ed(IX)*tiltck)/(one+dpp)                                     !hr02
qv=((-one*ed(IX))*tiltsk)/(one+dpp)                              !hr02
quz=((-one*ek(IX))*tiltck)/(one+dpp)                             !hr02
qvz=(ek(IX)*tiltsk)/(one+dpp)                                    !hr02
#endif
! end include/kickadpe.f90
