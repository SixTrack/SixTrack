! start include/alignva.f90
#ifndef TILT
  xlv(j)=xv(1,j)-xsiv(1,i)
  zlv(j)=xv(2,j)-zsiv(1,i)
  crkve=xlv(j)
  cikve=zlv(j)
#else
  xlv(j)=(xv(1,j)-xsiv(1,i))*tiltc(i)+(xv(2,j)-zsiv(1,i))*tilts(i)
  zlv(j)=(xv(2,j)-zsiv(1,i))*tiltc(i)-(xv(1,j)-xsiv(1,i))*tilts(i)
  crkve=xlv(j)
  cikve=zlv(j)
#endif
! end include/alignva.f90
