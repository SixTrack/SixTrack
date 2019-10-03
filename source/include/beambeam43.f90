! start include/beambeam43.f90
do j=1,napx
  r2b(j) = two*(sigman2(2,imbb(i))-sigman2(1,imbb(i)))
  rb(j)  = sqrt(r2b(j))
  rkb(j) = (strack(i)*pisqrt)/rb(j)
  if(ibbc == 0) then
    crkveb(j) = (xv1(j)-clobeam(1,imbb(i)))+parbe(ix,5)
    cikveb(j) = (xv2(j)-clobeam(2,imbb(i)))+parbe(ix,6)
  else
    crkveb(j) = ((xv1(j)-clobeam(1,imbb(i)))+parbe(ix,5))*bbcu(imbb(i),11)&
              + ((xv2(j)-clobeam(2,imbb(i)))+parbe(ix,6))*bbcu(imbb(i),12)
    cikveb(j) = ((xv2(j)-clobeam(2,imbb(i)))+parbe(ix,6))*bbcu(imbb(i),11)&
              - ((xv1(j)-clobeam(1,imbb(i)))+parbe(ix,5))*bbcu(imbb(i),12)
  end if
  xrb(j) = abs(crkveb(j))/rb(j)
  zrb(j) = abs(cikveb(j))/rb(j)
  tkb(j) = (crkveb(j)**2/sigman2(1,imbb(i))+cikveb(j)**2/sigman2(2,imbb(i)))*half
  xbb(j) = sigmanq(2,imbb(i))*xrb(j)
  zbb(j) = sigmanq(1,imbb(i))*zrb(j)
end do
if(ibtyp == 0) then
  do j=1,napx
    call errf(zrb(j),xrb(j),crzb(j),crxb(j))
    call errf(zbb(j),xbb(j),cbzb(j),cbxb(j))
  end do
else
  call wzsubv(napx,zrb(1),xrb(1),crzb(1),crxb(1))
  call wzsubv(napx,zbb(1),xbb(1),cbzb(1),cbxb(1))
end if
do j=1,napx
  if(ibbc == 0) then
    yv1(j) = yv1(j)+moidpsv(j)*((rkb(j)*(crzb(j)-exp_mb(-one*tkb(j))*cbzb(j)))*sign(one,crkveb(j))-beamoff(4,imbb(i)))
    yv2(j) = yv2(j)+moidpsv(j)*((rkb(j)*(crxb(j)-exp_mb(-one*tkb(j))*cbxb(j)))*sign(one,cikveb(j))-beamoff(5,imbb(i)))
  else
    cccc = ((rkb(j)*(crzb(j)-exp_mb(-one*tkb(j))*cbzb(j)))*sign(one,crkveb(j))-beamoff(4,imbb(i)))*bbcu(imbb(i),11)&
         - ((rkb(j)*(crxb(j)-exp_mb(-one*tkb(j))*cbxb(j)))*sign(one,cikveb(j))-beamoff(5,imbb(i)))*bbcu(imbb(i),12)
    yv1(j) = yv1(j)+moidpsv(j)*cccc
    cccc = ((rkb(j)*(crzb(j)-exp_mb(-one*tkb(j))*cbzb(j)))*sign(one,crkveb(j))-beamoff(4,imbb(i)))*bbcu(imbb(i),12)&
         + ((rkb(j)*(crxb(j)-exp_mb(-one*tkb(j))*cbxb(j)))*sign(one,cikveb(j))-beamoff(5,imbb(i)))*bbcu(imbb(i),11)
    yv2(j) = yv2(j)+moidpsv(j)*cccc
  end if
end do
! end include/beambeam43.f90
