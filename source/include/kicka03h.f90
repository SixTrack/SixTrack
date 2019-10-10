! start include/kicka03h.f90
  mpe=3
  mx=1
  cxzr=xl
  cxzi=zl
  cxzyr=cxzr**2-cxzi**2
  cxzyi=cxzr*cxzi+cxzi*cxzr
  dyy1=ekk*(tiltc(k)*cxzyr+tilts(k)*cxzyi)
  dyy2=ekk*(tilts(k)*cxzyr-tiltc(k)*cxzyi)
  tiltck=tiltc(k)*tiltc(k)-tilts(k)*tilts(k)
  tiltsk=(two*tiltc(k))*tilts(k)
  qu=(ekk*two)*(tiltck*xl+tiltsk*zl)
  qv=(ekk*two)*(tiltck*zl-tiltsk*xl)
  ab1(2)=qu
  ab2(2)=-one*qv
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck=tiltckuk
  ab1(3)=ekk*tiltck
  ab2(3)=ekk*tiltsk
! end include/kicka03h.f90
