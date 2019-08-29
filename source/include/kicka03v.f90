! start include/kicka03v.f90
#ifndef TILT
  mpe=3
  mx=1
  cxzr=xl
  cxzi=zl
  cxzyr=cxzr**2-cxzi**2
  cxzyi=cxzr*cxzi+cxzi*cxzr
  dyy1=ekk*cxzyi
  dyy2=ekk*cxzyr
  qu=(ekk*two)*zl
  qv=((-one*ekk)*two)*xl
  ab2(2)=-one*qv
  ab2(3)=ekk
#else
  mpe=3
  mx=1
  cxzr=xl
  cxzi=zl
  cxzyr=cxzr**2-cxzi**2
  cxzyi=cxzr*cxzi+cxzi*cxzr
  dyy1=ekk*(tiltc(k)*cxzyi-tilts(k)*cxzyr)
  dyy2=ekk*(tiltc(k)*cxzyr+tilts(k)*cxzyi)
  tiltck=tiltc(k)**2-tilts(k)**2
  tiltsk=(two*tiltc(k))*tilts(k)
  qu=(ekk*two)*(tiltck*zl-tiltsk*xl)
  qv=((-one*ekk)*two)*(tiltck*xl+tiltsk*zl)
  ab1(2)=qu
  ab2(2)=-one*qv
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck=tiltckuk
  ab1(3)=ekk*tiltsk
  ab2(3)=ekk*tiltck
#endif
! end include/kicka03v.f90
