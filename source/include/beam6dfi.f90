! start include/beam6dfi.f90
parbe(ix,4)=(((-one*crad)*ptnfac(ix))*half)*c1m6
!--Hirata's 6D beam-beam kick
dummy=dare(x(1))
!FOX      TRACKI(1)=(X(1)+PARBE(IX,5)-DUMMY)*C1M3 ;
!FOX      YP(1)=Y(1)*(ONE+DPDA)/MTCDA ;
dummy=dare(yp(1))
!FOX      TRACKI(2)=(YP(1)-DUMMY)*C1M3 ;
dummy=dare(x(2))
!FOX      TRACKI(3)=(X(2)+PARBE(IX,6)-DUMMY)*C1M3 ;
!FOX      YP(2)=Y(2)*(ONE+DPDA)/MTCDA ;
dummy=dare(yp(2))
!FOX      TRACKI(4)=(YP(2)-DUMMY)*C1M3 ;
dummy=dare(sigmda)
!FOX      TRACKI(5)=(SIGMDA-DUMMY)*C1M3 ;
dummy=dare(dpda)
!FOX      TRACKI(6)=DPDA-DUMMY ;

!We want to provide a set of canonical variables (15.03.2018)
!FOX      TRACKI(5)=TRACKI(5)/RV ;

call beaminf(tracki,parbe,sigz,bbcu,imbb(i),ix,ibbc)
if(ibeco.eq.1) then
  beamoff1=dare(tracki(1))*c1e3
  beamoff2=dare(tracki(3))*c1e3
  beamoff3=dare(tracki(5))*c1e3 !Added 14.03.2018
  beamoff4=dare(tracki(2))*c1e3
  beamoff5=dare(tracki(4))*c1e3
  beamoff6=dare(tracki(6))
else
  beamoff1=zero
  beamoff2=zero
  beamoff3=zero
  beamoff4=zero
  beamoff5=zero
  beamoff6=zero
endif
dummy=dare(x(1))
!FOX      X(1)=TRACKI(1)*C1E3+DUMMY-BEAMOFF1 ;
dummy=dare(x(2))
!FOX      X(2)=TRACKI(3)*C1E3+DUMMY-BEAMOFF2 ;
dummy=dare(yp(1))
!FOX      YP(1)=TRACKI(2)*C1E3+DUMMY-BEAMOFF4 ;
dummy=dare(yp(2))
!FOX      YP(2)=TRACKI(4)*C1E3+DUMMY-BEAMOFF5 ;
dummy=dare(dpda)

!FOX      DPDA=TRACKI(6)+DUMMY-BEAMOFF6 ;
!FOX      DPDA1=DPDA*C1E3 ;
!FOX      MOIDA=MTCDA/(ONE+DPDA) ;
!FOX      Y(1)=YP(1)*MTCDA/(ONE+DPDA) ;
!FOX      Y(2)=YP(2)*MTCDA/(ONE+DPDA) ;
!FOX      EJF1=E0F*(ONE+DPDA)/(NUCM0/NUCMDA) ;
!FOX      EJ1=SQRT(EJF1*EJF1+NUCMDA*NUCMDA) ;
!FOX      RV=EJ1/E0*E0F/EJF1 ;

!We want to go back to sixtrack variables (15.03.2018)
!FOX      TRACKI(5) = TRACKI(5)*RV ;

dummy=dare(sigmda)  !Added 14.03.2018
!FOX      SIGMDA=TRACKI(5)*C1E3+DUMMY-BEAMOFF3 ;
if(ithick.eq.1) call envada
! end include/beam6dfi.f90
