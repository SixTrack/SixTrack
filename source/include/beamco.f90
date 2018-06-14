if(ibbc.eq.0) then
    crkveb(j)=(xv(1,j)-clobeam(1,imbb(i)))+parbe(ix,5)
    cikveb(j)=(xv(2,j)-clobeam(2,imbb(i)))+parbe(ix,6)
  else
    crkveb(j)=((xv(1,j)-clobeam(1,imbb(i)))+parbe(ix,5))*bbcu(imbb(i),11)&
             +((xv(2,j)-clobeam(2,imbb(i)))+parbe(ix,6))*bbcu(imbb(i),12)
    cikveb(j)=((xv(2,j)-clobeam(2,imbb(i)))+parbe(ix,6))*bbcu(imbb(i),11)&
             -((xv(1,j)-clobeam(1,imbb(i)))+parbe(ix,5))*bbcu(imbb(i),12)
end if
