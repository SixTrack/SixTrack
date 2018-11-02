! start include/alignvb.f90
#ifndef TILT
  xlvj=xv1(j)-xsiv(1,i)
  zlvj=xv2(j)-zsiv(1,i)
#else
  xlvj=(xv1(j)-xsiv(1,i))*tiltc(i)+(xv2(j)-zsiv(1,i))*tilts(i)
  zlvj=(xv2(j)-zsiv(1,i))*tiltc(i)-(xv1(j)-xsiv(1,i))*tilts(i)
#endif
! end include/alignvb.f90
