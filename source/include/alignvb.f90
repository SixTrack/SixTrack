! start include/alignvb.f90
#ifndef TILT
  xlvj=xv1(j)-xsiv(i)
  zlvj=xv2(j)-zsiv(i)
#else
  xlvj=(xv1(j)-xsiv(i))*tiltc(i)+(xv2(j)-zsiv(i))*tilts(i)
  zlvj=(xv2(j)-zsiv(i))*tiltc(i)-(xv1(j)-xsiv(i))*tilts(i)
#endif
! end include/alignvb.f90
