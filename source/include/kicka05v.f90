!start kicka05v.f90
#ifndef TILT
  mpe=5
  mx=3
  cxzr=xl
  cxzi=zl
  cxzyr=cxzr**2-cxzi**2                                            !hr08
  cxzyi=cxzr*cxzi+cxzi*cxzr
  ab1(3)=(six*ekk)*cxzyi                                           !hr02
  ab2(3)=(six*ekk)*cxzyr                                           !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  qu=(four*ekk)*cxzyi                                              !hr02
  qv=((-one*four)*ekk)*cxzyr                                       !hr02
  ab2(2)=-one*qv                                                   !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  dyy1=ekk*cxzyi
  dyy2=ekk*cxzyr
  ab1(4)=(four*ekk)*zl                                             !hr08
  ab2(4)=(four*ekk)*xl                                             !hr08
  ab2(5)=ekk
#else
  mpe=5
  mx=3
  cxzr=xl
  cxzi=zl
  cxzyr=cxzr**2-cxzi**2                                            !hr08
  cxzyi=cxzr*cxzi+cxzi*cxzr
  tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
  tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk1=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck1=tiltckuk
  ab1(3)=(six*ekk)*(tiltck1*cxzyi-tiltsk1*cxzyr)                   !hr02
  ab2(3)=(six*ekk)*(tiltck1*cxzyr+tiltsk1*cxzyi)                   !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  qu=(four*ekk)*(tiltck*cxzyi-tiltsk*cxzyr)
  qv=((-one*four)*ekk)*(tiltck*cxzyr+tiltsk*cxzyi)
  ab1(2)=qu
  ab2(2)=-one*qv
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  dyy1=ekk*(tiltc(k)*cxzyi-tilts(k)*cxzyr)
  dyy2=ekk*(tiltc(k)*cxzyr+tilts(k)*cxzyi)
  tiltckuk=tiltck1*tiltc(k)-tiltsk1*tilts(k)
  tiltsk=tiltck1*tilts(k)+tiltsk1*tiltc(k)
  tiltck=tiltckuk
  ab1(4)=(four*ekk)*(tiltck*zl-tiltsk*xl)                          !hr02
  ab2(4)=(four*ekk)*(tiltck*xl+tiltsk*zl)                          !hr02
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck=tiltckuk
  ab1(5)=ekk*tiltsk
  ab2(5)=ekk*tiltck
#endif

!end kicka05v.f90
