!Here it is in angles because the transformation is afterwards
!FOX  YV1J=BBIV(1,I)+BBIV(2,I)*XL+AAIV(2,I)*ZL;
!FOX  YV2J=AAIV(1,I)-BBIV(2,I)*ZL+AAIV(2,I)*XL;
if(dki(ix,3) > pieni .and. abs(dki(ix,1)) > pieni .and. curveff) then !Horizontal curvature effect for quadrupoles

!FOX  YV1J=YV1J -
!FOX  BBIV(2,I)*((DKI(IX,1)/DKI(IX,3))*(XL*XL-0.5*ZL*ZL)*C1M3) ;

!FOX  YV2J=YV2J +
!FOX  BBIV(2,I)*(DKI(IX,1)/DKI(IX,3))*
!FOX  (XL*ZL*C1M3-C1M6*(DKI(IX,1)/DKI(IX,3))*(ZL*ZL*ZL/(6.0)));

endif
!FOX  CRKVE=XL ;
!FOX  CIKVE=ZL ;
