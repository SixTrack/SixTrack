! start include/beam6d.f90
!--Hirata's 6D beam-beam kick
do j=1,napx
  track6d(1,j)=((xv1(j)+parbe(ix,5)) - clobeam(1,imbb(i)))*c1m3
  track6d(2,j)=(yv1(j)/moidpsv(j)-clobeam(4,imbb(i)))*c1m3
  track6d(3,j)=((xv2(j)+parbe(ix,6)) -  clobeam(2,imbb(i)))*c1m3
  track6d(4,j)=(yv2(j)/moidpsv(j)-clobeam(5,imbb(i)))*c1m3
  track6d(5,j)=(sigmv(j)-clobeam(3,imbb(i)))*c1m3
  track6d(6,j)=dpsv(j)-clobeam(6,imbb(i))

  !We want to provide a set of canonical variables (15.03.2018)
  track6d(5,j) = track6d(5,j)/rvv(j)

end do
call beamint(napx,track6d,parbe,sigz,bbcu,imbb(i),ix,ibtyp,ibbc,mtc)
do j=1,napx
  xv1(j)=(track6d(1,j)*c1e3+clobeam(1,imbb(i)))-beamoff(1,imbb(i))
  xv2(j)=(track6d(3,j)*c1e3+clobeam(2,imbb(i)))-beamoff(2,imbb(i))
  dpsv(j)=(track6d(6,j)+clobeam(6,imbb(i)))-beamoff(6,imbb(i))
enddo
call part_updatePartEnergy(3, .false.)

do j=1,napx
  
  yv1(j)=((track6d(2,j)*c1e3+clobeam(4,imbb(i)))-beamoff(4,imbb(i)))*moidpsv(j)
  yv2(j)=((track6d(4,j)*c1e3+clobeam(5,imbb(i)))-beamoff(5,imbb(i)))*moidpsv(j)
  
  !We want to go back to sixtrack variables (15.03.2018)
  track6d(5,j) = track6d(5,j)*rvv(j)
  sigmv(j)=(track6d(5,j)*c1e3+clobeam(3,imbb(i)))- beamoff(3,imbb(i))
end do

if(ithick == 1) call envarsv(dpsv,moidpsv,rvv,ekv)
! end include/beam6d.f90
