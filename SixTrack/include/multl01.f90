#ifndef TILT
  qu=(((-one*dki(ix,1))/dki(ix,3))*dki(ix,1))/(one+dpp)        !hr03
  dppi=(c1e3*dki(ix,1))/(one+dpp)                              !hr03
  t(1,2)=(t(1,2)+qu*xl)-dppi*dpp                               !hr03
#else
  qu=(((-one*dki(ix,1))/dki(ix,3))*dki(ix,1))/(one+dpp)        !hr03
  dppi=(c1e3*dki(ix,1))/(one+dpp)                              !hr03
  t(1,2)=t(1,2)+(qu*xl-dppi*dpp)*tiltc(k)+dppi*(one-tiltc(k))
  t(1,4)=(t(1,4)+(qu*xl-dppi*dpp)*tilts(k))+dppi*tilts(k)
#endif
