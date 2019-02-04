! start include/kispa10v.f90
#ifndef TILT
  ekko=ekk
  cxzyr=cxzr**2-cxzi**2                                          !hr08
  cxzyi=cxzr*cxzi+cxzi*cxzr
  ekk=(36.0_fPrec*ekko)*cxzyi                                          !hr03
  call detune(4,ekk,ep,beta,dtu,dtup,dfac)
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  ekk=(126.0_fPrec*ekko)*cxzyi                                         !hr03
  call detune(3,ekk,ep,beta,dtu,dtup,dfac)
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  ekk=(84.0_fPrec*ekko)*cxzyi                                          !hr03
  call detune(2,ekk,ep,beta,dtu,dtup,dfac)
#else
  tiltck=tiltc(k)**2-tilts(k)**2                                 !hr08
  tiltsk=two*tiltc(k)*tilts(k)
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck=tiltckuk
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk4=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck4=tiltckuk
  tiltckuk=tiltck4*tiltc(k)-tiltsk4*tilts(k)
  tiltsk=tiltck4*tilts(k)+tiltsk4*tiltc(k)
  tiltck=tiltckuk
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk6=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck6=tiltckuk
  tiltckuk=tiltck6*tiltc(k)-tiltsk6*tilts(k)
  tiltsk=tiltck6*tilts(k)+tiltsk6*tiltc(k)
  tiltck=tiltckuk
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk8=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck8=tiltckuk
  tiltckuk=tiltck8*tiltc(k)-tiltsk8*tilts(k)
  tiltsk=tiltck8*tilts(k)+tiltsk8*tiltc(k)
  tiltck=tiltckuk
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk10=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck10=tiltckuk
  ekko=ekk
  ekk=(-one*ekko)*tiltsk10                                       !hr03
  call detune(5,ekk,ep,beta,dtu,dtup,dfac)
  cxzyr=cxzr*cxzr-cxzi*cxzi
  cxzyi=cxzr*cxzi+cxzi*cxzr
  ekk=(36.0_fPrec*ekko)*(tiltck8*cxzyi-tiltsk8*cxzyr)                  !hr03
  call detune(4,ekk,ep,beta,dtu,dtup,dfac)
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  ekk=(126.0_fPrec*ekko)*(tiltck6*cxzyi-tiltsk6*cxzyr)                 !hr03
  call detune(3,ekk,ep,beta,dtu,dtup,dfac)
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  ekk=(84.0_fPrec*ekko)*(tiltck4*cxzyi-tiltsk4*cxzyr)                  !hr03
  call detune(2,ekk,ep,beta,dtu,dtup,dfac)
#endif
! end include/kispa10v.f90
