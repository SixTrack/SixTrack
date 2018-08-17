!start alignvb.f90
#ifndef TILT
  xlvj=xv(1,j)-xsiv(1,i)
  zlvj=xv(2,j)-zsiv(1,i)
#else
  xlvj=(xv(1,j)-xsiv(1,i))*tiltc(i)+(xv(2,j)-zsiv(1,i))*tilts(i)
  zlvj=(xv(2,j)-zsiv(1,i))*tiltc(i)-(xv(1,j)-xsiv(1,i))*tilts(i)
#endif

!end alignvb.f90
