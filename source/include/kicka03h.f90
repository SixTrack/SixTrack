! start include/kicka03h.f90
#ifndef TILT
  mpe=3
  mx=1
  cxzr=xl
  cxzi=zl
  cxzyr=cxzr**2-cxzi**2                                            !hr08
  cxzyi=cxzr*cxzi+cxzi*cxzr
  dyy1=ekk*cxzyr
  dyy2=(-one*ekk)*cxzyi                                            !hr02
  qu=(ekk*two)*xl                                                  !hr02
  qv=(ekk*two)*zl                                                  !hr02
  ab2(2)=-one*qv                                                   !hr02
  ab1(3)=ekk
#else
  mpe=3
  mx=1
  cxzr=xl
  cxzi=zl
  cxzyr=cxzr**2-cxzi**2                                            !hr08
  cxzyi=cxzr*cxzi+cxzi*cxzr
  dyy1=ekk*(tiltc(k)*cxzyr+tilts(k)*cxzyi)
  dyy2=ekk*(tilts(k)*cxzyr-tiltc(k)*cxzyi)                         !hr02
  tiltck=tiltc(k)*tiltc(k)-tilts(k)*tilts(k)
  tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
  qu=(ekk*two)*(tiltck*xl+tiltsk*zl)                               !hr02
  qv=(ekk*two)*(tiltck*zl-tiltsk*xl)                               !hr02
  ab1(2)=qu
  ab2(2)=-one*qv                                                   !hr02
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck=tiltckuk
  ab1(3)=ekk*tiltck
  ab2(3)=ekk*tiltsk
#endif
! end include/kicka03h.f90
