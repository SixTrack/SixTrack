#ifndef TILT
!FOX  Y(2)=Y(2)+EKK ;
#else
!FOX  Y(1)=Y(1)-EKK*TILTS(I) ;
!FOX  Y(2)=Y(2)+EKK*TILTC(I) ;
#endif
