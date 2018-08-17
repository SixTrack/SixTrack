! start include/beamcoo.f90
if(ibbc.eq.0) then
  crkveb(j)=parbe(ix,5)
  cikveb(j)=parbe(ix,6)
else
  crkveb(j)=parbe(ix,5)*bbcu(imbb(i),11)+parbe(ix,6)*bbcu(imbb(i),12)
  cikveb(j)=parbe(ix,6)*bbcu(imbb(i),11)-parbe(ix,5)*bbcu(imbb(i),12)
endif
! end include/beamcoo.f90
