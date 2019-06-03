cos_t = cos_mb(temp_angle)
sin_t = sin_mb(temp_angle)
tan_t = tan_mb(temp_angle)

print *, "gooiinggg hereee first", temp_angle,  sin_t, cos_t, cos_mb(temp_angle)
!FOX  YP(1)=Y(1)*(ONE+DPDA)/MTCDA ;
!FOX  YP(2)=Y(2)*(ONE+DPDA)/MTCDA ;

!FOX  TEMPI(1) = X(1)*C1M3 ;
!FOX  TEMPI(2) = YP(1)*C1M3 ;
!FOX  TEMPI(3) = X(2)*C1M3 ;
!FOX  TEMPI(4) = YP(2)*C1M3 ;
!FOX  TEMPI(5) = SIGMDA*C1M3 ;
!FOX  TEMPI(6) = ((EJ1-E0)/E0F) ;


!FOX  ZTDA = SQRT((ONE + DPDA)*(ONE + DPDA) 
!FOX  - TEMPI(2)*TEMPI(2) - TEMPI(4)*TEMPI(4)) ;

!FOX  PTTDA = ONE - (TAN_T*TEMPI(4))/ZTDA ;

!FOX  X(1) = X(1) + 
!FOX  C1E3*(TAN_T*TEMPI(3)*TEMPI(2)/(ZTDA*PTTDA)) ;
!FOX  X(2) = C1E3*TEMPI(3)/(COS_T*PTTDA) ;
!FOX  Y(2) = C1E3*(COS_T*TEMPI(4) + SIN_T*ZTDA)/(ONE+DPDA)/MTCDA ;
!FOX  SIGMDA = SIGMDA - C1E3*((TAN_T*TEMPI(3)*
!FOX  (ONE/(E0F/E0)+TEMPI(6))/(ZTDA*PTTDA))*(E0F/E0)) ;
