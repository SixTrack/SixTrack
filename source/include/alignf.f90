! start include/alignf.f90
#ifndef TILT
!FOX  XL=X(1)-XS ;
!FOX  ZL=X(2)-ZS ;
!FOX  CRKVE=XL ;
!FOX  CIKVE=ZL ;
#else
!FOX  XL=(X(1)-XS)*TILTC(I)+(X(2)-ZS)*TILTS(I) ;
!FOX  ZL=-(X(1)-XS)*TILTS(I)+(X(2)-ZS)*TILTC(I) ;
!FOX  CRKVE=XL ;
!FOX  CIKVE=ZL ;
#endif
! end include/alignf.f90
