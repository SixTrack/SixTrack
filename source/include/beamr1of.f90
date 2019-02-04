! start include/beamr1of.f90
if(ibbc.eq.0) then
  crk=parbe(ix,5)
  cik=parbe(ix,6)
else
  crk=parbe(ix,5)*bbcu(imbb(i),11) + parbe(ix,6)*bbcu(imbb(i),12)
  cik=parbe(ix,6)*bbcu(imbb(i),11) - parbe(ix,5)*bbcu(imbb(i),12)
endif
rho2b=crk**2+cik**2                                          !hr03
if(rho2b.gt.pieni) &
! end include/beamr1of.f90
