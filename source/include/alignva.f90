! start include/alignva.f90
#ifndef TILT
  xlv   = xv1(j)-xsiv(i)
  zlv   = xv2(j)-zsiv(i)
  crkve = xlv
  cikve = zlv
#else
  xlv   = (xv1(j)-xsiv(i))*tiltc(i)+(xv2(j)-zsiv(i))*tilts(i)
  zlv   = (xv2(j)-zsiv(i))*tiltc(i)-(xv1(j)-xsiv(i))*tilts(i)
  crkve = xlv
  cikve = zlv
#endif
! end include/alignva.f90
