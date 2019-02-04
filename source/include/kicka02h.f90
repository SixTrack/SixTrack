! start include/kicka02h.f90
#ifndef TILT
dyy1=ekk*xl
dyy2=(-one*ekk)*zl                                                     !hr02
mpe=20
qu=ekk
qv=zero
#else
dyy1=ekk*(tiltc(k)*xl+tilts(k)*zl)
dyy2=ekk*(tilts(k)*xl-tiltc(k)*zl)                               !hr08
mpe=20
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=ekk*tiltck
qv=(-one*ekk)*tiltsk                                             !hr02
ab1(2)=qu
ab2(2)=-one*qv                                                   !hr02
#endif
! end include/kicka02h.f90
