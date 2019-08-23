! start include/kicka02h.f90
#ifndef TILT
dyy1=ekk*xl
dyy2=(-one*ekk)*zl
mpe=20
qu=ekk
qv=zero
#else
dyy1=ekk*(tiltc(k)*xl+tilts(k)*zl)
dyy2=ekk*(tilts(k)*xl-tiltc(k)*zl)
mpe=20
tiltck=tiltc(k)**2-tilts(k)**2
tiltsk=(two*tiltc(k))*tilts(k)
qu=ekk*tiltck
qv=(-one*ekk)*tiltsk
ab1(2)=qu
ab2(2)=-one*qv
#endif
! end include/kicka02h.f90
