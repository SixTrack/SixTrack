! start include/kicka02v.f90
#ifndef TILT
  dyy1=ekk*zl
  dyy2=ekk*xl
  mpe=2
  mx=-1
  qu=zero
  qv=-one*ekk                                                      !hr02
  ab2(2)=ekk
#else
  dyy1=ekk*(tiltc(k)*zl-tilts(k)*xl)
  dyy2=ekk*(tiltc(k)*xl+tilts(k)*zl)
  mpe=2
  mx=-1
  tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
  tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
  qu=(-one*ekk)*tiltsk                                             !hr02
  qv=(-one*ekk)*tiltck                                             !hr02
  ab1(2)=qu
  ab2(2)=-one*qv                                                   !hr02
#endif
! end include/kicka02v.f90
