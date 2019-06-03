cos_t = cos_mb(temp_angle)
sin_t = -sin_mb(temp_angle)
!FOX  YP(1)=Y(1)*(ONE+DPDA)/MTCDA ;
!FOX  YP(2)=Y(2)*(ONE+DPDA)/MTCDA ;

!FOX  TEMPI(1) = X(1) ;
!FOX  TEMPI(2) = YP(1) ;
!FOX  TEMPI(3) = X(2) ;
!FOX  TEMPI(4) = YP(2) ;

!FOX  X(1) = TEMPI(1)*COS_T - TEMPI(3)*SIN_T ;
!FOX  YP(1) = TEMPI(2)*COS_T - TEMPI(4)*SIN_T ;
!FOX  X(2) = TEMPI(1)*SIN_T + TEMPI(3)*COS_T ;
!FOX  YP(2) = TEMPI(2)*SIN_T + TEMPI(4)*COS_T ;

!FOX  Y(1)=YP(1)*MTCDA/(ONE+DPDA) ;
!FOX  Y(2)=YP(2)*MTCDA/(ONE+DPDA) ;