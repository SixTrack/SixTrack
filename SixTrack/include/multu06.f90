#ifndef TILT
  y(1,2)=y(1,2)+((dki(ix,2)*dpp)/(one+dpp))*c1e3               !hr03
#else
  y(1,1)=(y(1,1)-(((dki(ix,2)*dpp)/(one+dpp))*c1e3)*tilts(k))+((dki(ix,2)/(one+dpp))*c1e3)*tilts(k)
  y(1,2)=(y(1,2)+(((dki(ix,2)*dpp)/(one+dpp))*c1e3)*tiltc(k))-((dki(ix,2)/(one+dpp))*c1e3)*(one-tiltc(k))
#endif
