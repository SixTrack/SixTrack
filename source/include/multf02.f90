!Here it is in angles because the transformation is afterwards
if(dki(I,3).gt. pieni) then
!FOX  YV1J=BBI(I,1)+BBI(I,2)*XL+AAI(I,2)*ZL + 
!FOX  BBI(I,2)*((DKI(I,1)/DKI(I,3))*(XL*XL-0.5*ZL*ZL)*C1M3) ; 

!FOX  YV2J=AAI(I,1)-BBI(I,2)*ZL+AAI(I,2)*XL-
!FOX  BBI(I,2)*(DKI(I,1)/DKI(I,3))*(XL*ZL*C1M3+C1M6*(DKI(I,1))*(ZL*ZL*ZL/(6.0)));
else ! In case there is no ficitive length

!FOX  YV1J=BBI(I,1)+BBI(I,2)*XL+AAI(I,2)*ZL;
!FOX  YV2J=AAI(I,1)-BBI(I,2)*ZL+AAI(I,2)*XL;
endif

!FOX  CRKVE=XL ;
!FOX  CIKVE=ZL ;