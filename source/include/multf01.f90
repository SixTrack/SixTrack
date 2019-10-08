  if(abs(dki(ix,1)).gt.pieni) then
    if(abs(dki(ix,3)).gt.pieni) then
!FOX  DKIP=DKI(IX,1)*MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=Y(1)-(DKI(IX,1)/DKI(IX,3)*XL+DPDA*C1E3)*
!FOX  TILTC(I)*DKIP
!FOX  +C1E3*DKI(IX,1)*MTCDA/(ONE+DPDA)*(ONE-TILTC(I)) ;
!FOX  Y(2)=Y(2)-(DKI(IX,1)/DKI(IX,3)*XL+DPDA*C1E3)*
!FOX  TILTS(I)*DKIP
!FOX  +C1E3*DKI(IX,1)*MTCDA/(ONE+DPDA)*TILTS(I) ;
    else
!FOX  Y(1)=Y(1)-DKI(IX,1)*DPDA*C1E3*MTCDA/(ONE+DPDA)*TILTC(I)
!FOX  +C1E3*DKI(IX,1)*MTCDA/(ONE+DPDA)*(ONE-TILTC(I)) ;
!FOX  Y(2)=Y(2)-DKI(IX,1)*DPDA*C1E3*MTCDA/(ONE+DPDA)*TILTS(I)
!FOX  +C1E3*DKI(IX,1)*MTCDA/(ONE+DPDA)*TILTS(I) ;
    endif
    if(idp.eq.1.and.iabs(ition).eq.1) then
!FOX  SIGMDA=SIGMDA+RV*DKI(IX,1)*XL ;
    endif
  endif
  if(abs(dki(ix,2)).gt.pieni) then
    if(abs(dki(ix,3)).gt.pieni) then
!FOX  DKIP=DKI(IX,2)*MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=Y(1)+(DKI(IX,2)/DKI(IX,3)*ZL-DPDA*C1E3)*
!FOX  TILTS(I)*DKIP
!FOX  +C1E3*DKI(IX,2)*MTCDA/(ONE+DPDA)*TILTS(I) ;
!FOX  Y(2)=Y(2)-(DKI(IX,2)/DKI(IX,3)*ZL-DPDA*C1E3)*
!FOX  TILTC(I)*DKIP
!FOX  -C1E3*DKI(IX,2)*MTCDA/(ONE+DPDA)*(ONE-TILTC(I)) ;
    else
!FOX  Y(1)=Y(1)-DKI(IX,2)*DPDA*C1E3*MTCDA/(ONE+DPDA)*TILTS(I)
!FOX  +C1E3*DKI(IX,2)*MTCDA/(ONE+DPDA)*TILTS(I) ;
!FOX  Y(2)=Y(2)+DKI(IX,2)*DPDA*C1E3*MTCDA/(ONE+DPDA)*TILTC(I)
!FOX  -C1E3*DKI(IX,2)*MTCDA/(ONE+DPDA)*(ONE-TILTC(I)) ;
    endif
    if(idp.eq.1.and.iabs(ition).eq.1) then
!FOX  SIGMDA=SIGMDA-RV*DKI(IX,2)*ZL ;
    endif
  endif
