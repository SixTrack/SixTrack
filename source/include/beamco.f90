! start include/beamco.f90
if(ibbc.eq.0) then
    crkveb(j)=(xv1(j)-clobeam(1,imbb(i)))+parbe(ix,5)
    cikveb(j)=(xv2(j)-clobeam(2,imbb(i)))+parbe(ix,6)
  else
    crkveb(j)=((xv1(j)-clobeam(1,imbb(i)))+parbe(ix,5))*bbcu(imbb(i),11)&
             +((xv2(j)-clobeam(2,imbb(i)))+parbe(ix,6))*bbcu(imbb(i),12)
    cikveb(j)=((xv2(j)-clobeam(2,imbb(i)))+parbe(ix,6))*bbcu(imbb(i),11)&
             -((xv1(j)-clobeam(1,imbb(i)))+parbe(ix,5))*bbcu(imbb(i),12)
end if
! end include/beamco.f90
