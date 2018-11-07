! start include/alignva.f90
#ifndef TILT
  xlv(j)=xv1(j)-xsiv(i)
  zlv(j)=xv2(j)-zsiv(i)
  crkve=xlv(j)
  cikve=zlv(j)
#else
  xlv(j)=(xv1(j)-xsiv(i))*tiltc(i)+(xv2(j)-zsiv(i))*tilts(i)
  zlv(j)=(xv2(j)-zsiv(i))*tiltc(i)-(xv1(j)-xsiv(i))*tilts(i)
  crkve=xlv(j)
  cikve=zlv(j)
#endif
! end include/alignva.f90
