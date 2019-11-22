! start include/beamcof.f90
!FOX  CRKVEBF=X(1) ;
!FOX  CIKVEBF=X(2) ;
startco=(dare(x(1))-clobeam(1,imbb(i)))+parbe(ix,5)
call dapok(crkvebf,jj,startco)
startco=(dare(x(2))-clobeam(2,imbb(i)))+parbe(ix,6)
call dapok(cikvebf,jj,startco)
if(ibbc.eq.1) then
!FOX  CCCC=CRKVEBF ;
!FOX  CRKVEBF=CCCC*BBCU(IMBB(I),11)+CIKVEBF*BBCU(IMBB(I),12) ;
!FOX  CIKVEBF=-CCCC*BBCU(IMBB(I),12)+CIKVEBF*BBCU(IMBB(I),11) ;
endif
! end include/beamcof.f90
