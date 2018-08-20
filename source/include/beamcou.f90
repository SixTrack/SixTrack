! start include/beamcou.f90
dps0=dps(1)
dps(1)=zero
iqmodc=4
call mydaini(1,2,nd2,nd,nd2,1)
ilinc=2
call mydaini(2,2,nd2,nd,nd2,1)
dps(1)=dps0
! end include/beamcou.f90
