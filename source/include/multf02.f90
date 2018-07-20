!Here it is in angles because the transformation is afterwards
!FOX  YV1J=BBI(I,1)+BBI(I,2)*XL+AAI(I,2)*ZL;
!FOX  YV2J=AAI(I,1)-BBI(I,2)*ZL+AAI(I,2)*XL;
if(dki(I,3).gt. pieni .and. abs(dki(I,1)).gt. pieni .and. curveff.eq.1) then !Horizontal curvature effect for quadrupoles

!FOX  YV1J=YV1J + 
!FOX  BBI(I,2)*((DKI(I,1)/DKI(I,3))*(XL*XL-0.5*ZL*ZL)*C1M3) ; 

!FOX  YV2J=YV2J -
!FOX  BBI(I,2)*(DKI(I,1)/DKI(I,3))*(XL*ZL*C1M3+C1M6*(DKI(I,1))*(ZL*ZL*ZL/(6.0)));

else if(dki(I,3).gt. pieni .and. abs(dki(I,2)).gt. pieni .and. curveff.eq.1) then !Vertical curvature effect for quadrupoles
!FOX  YV1J=YV1J + 
!FOX  BBI(I,2)*(DKI(I,2)/DKI(I,3))*(XL*ZL*C1M3+C1M6*(DKI(I,1))*(XL*XL*XL/(6.0))); 

!FOX  YV2J=YV2J -
!FOX  BBI(I,2)*((DKI(I,2)/DKI(I,3))*(ZL*ZL-0.5*XL*XL)*C1M3); 

endif
!FOX  CRKVE=XL ;
!FOX  CIKVE=ZL ;