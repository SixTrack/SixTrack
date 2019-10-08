! start include/kicka02v.f90
  dyy1=ekk*(tiltc(k)*zl-tilts(k)*xl)
  dyy2=ekk*(tiltc(k)*xl+tilts(k)*zl)
  mpe=2
  mx=-1
  tiltck=tiltc(k)**2-tilts(k)**2
  tiltsk=(two*tiltc(k))*tilts(k)
  qu=(-one*ekk)*tiltsk
  qv=(-one*ekk)*tiltck
  ab1(2)=qu
  ab2(2)=-one*qv
! end include/kicka02v.f90
