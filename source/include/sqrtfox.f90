#ifdef FAST
!FOX  EL(JX)*(C1E3-RV*(C1E3+(Y(1)*Y(1)+Y(2)*Y(2))*C5M4)) ;
#else
!FOX  EL(JX)*(C1E3-RV*SQRT(C1E6+Y(1)*Y(1)+Y(2)*Y(2))) ;
#endif
