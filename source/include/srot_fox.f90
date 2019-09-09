! start include/srot_fox.f90
cos_t = cos_mb(temp_angle)
sin_t = -sin_mb(temp_angle)


!FOX  TEMPI(1) = X(1) ;
!FOX  TEMPI(2) = Y(1) ;
!FOX  TEMPI(3) = X(2) ;
!FOX  TEMPI(4) = Y(2) ;

!FOX  X(1) = TEMPI(1)*COS_T - TEMPI(3)*SIN_T ;
!FOX  Y(1) = TEMPI(2)*COS_T - TEMPI(4)*SIN_T ;
!FOX  X(2) = TEMPI(1)*SIN_T + TEMPI(3)*COS_T ;
!FOX  Y(2) = TEMPI(2)*SIN_T + TEMPI(4)*COS_T ;
! end include/srot_fox.f90
