#ifndef TILT
  qu=((dki(ix,2)/dki(ix,3))*dki(ix,2))/(one+dpp)               !hr03
  y(1,2)=(y(1,2)-qu*zl)+((dpp*c1e3)*dki(ix,2))/(one+dpp)       !hr03
#else
  qu=((dki(ix,2)/dki(ix,3))*dki(ix,2))/(one+dpp)               !hr03
  y(1,1)=(y(1,1)+(qu*zl-((dpp*c1e3)*dki(ix,2))/(one+dpp))*tilts(k))+((c1e3*dki(ix,2))/(one+dpp))*tilts(k)
  y(1,2)=(y(1,2)+(((dpp*c1e3)*dki(ix,2))/(one+dpp)-qu*zl)*tiltc(k))-((c1e3*dki(ix,2))/(one+dpp))*(one-tiltc(k))
#endif
