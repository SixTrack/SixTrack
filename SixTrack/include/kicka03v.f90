#ifndef TILT
  mpe=3
  mx=1
  cxzr=xl
  cxzi=zl
  cxzyr=cxzr**2-cxzi**2                                            !hr08
  cxzyi=cxzr*cxzi+cxzi*cxzr
  dyy1=ekk*cxzyi
  dyy2=ekk*cxzyr
  qu=(ekk*two)*zl                                                  !hr02
  qv=((-one*ekk)*two)*xl                                           !hr02
  ab2(2)=-one*qv
  ab2(3)=ekk
#else
  mpe=3
  mx=1
  cxzr=xl
  cxzi=zl
  cxzyr=cxzr**2-cxzi**2                                            !hr08
  cxzyi=cxzr*cxzi+cxzi*cxzr
  dyy1=ekk*(tiltc(k)*cxzyi-tilts(k)*cxzyr)
  dyy2=ekk*(tiltc(k)*cxzyr+tilts(k)*cxzyi)
  tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
  tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
  qu=(ekk*two)*(tiltck*zl-tiltsk*xl)                               !hr02
  qv=((-one*ekk)*two)*(tiltck*xl+tiltsk*zl)                        !hr02
  ab1(2)=qu
  ab2(2)=-one*qv                                                   !hr02
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck=tiltckuk
  ab1(3)=ekk*tiltsk
  ab2(3)=ekk*tiltck
#endif
