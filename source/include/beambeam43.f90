! start include/beambeam43.f90
if(ibbc == 0) then
  crkveb(1:napx) = (xv1(1:napx)-clobeam(1,imbb(i))) + parbe(ix,5)
  cikveb(1:napx) = (xv2(1:napx)-clobeam(2,imbb(i))) + parbe(ix,6)
else
  crkveb(1:napx) = ((xv1(1:napx)-clobeam(1,imbb(i))) + parbe(ix,5))*bbcu(imbb(i),11) &
                 + ((xv2(1:napx)-clobeam(2,imbb(i))) + parbe(ix,6))*bbcu(imbb(i),12)
  cikveb(1:napx) = ((xv2(1:napx)-clobeam(2,imbb(i))) + parbe(ix,6))*bbcu(imbb(i),11) &
                 - ((xv1(1:napx)-clobeam(1,imbb(i))) + parbe(ix,5))*bbcu(imbb(i),12)
end if
rb(1:napx)  = sqrt(two*(sigman2(2,imbb(i))-sigman2(1,imbb(i))))
rkb(1:napx) = (strack(i)*pisqrt)/rb(1:napx)
xrb(1:napx) = abs(crkveb(1:napx))/rb(1:napx)
zrb(1:napx) = abs(cikveb(1:napx))/rb(1:napx)
tkb(1:napx) = (crkveb(1:napx)**2/sigman2(1,imbb(i))+cikveb(1:napx)**2/sigman2(2,imbb(i)))*half
xbb(1:napx) = sigmanq(2,imbb(i))*xrb(1:napx)
zbb(1:napx) = sigmanq(1,imbb(i))*zrb(1:napx)
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
    yv1(j) = yv1(j) + moidpsv(j)*( &
        ((rkb(j)*(crzb(j)-exp_mb(-one*tkb(j))*cbzb(j)))*sign(one,crkveb(j))-beamoff(4,imbb(i)))*bbcu(imbb(i),11) &
      - ((rkb(j)*(crxb(j)-exp_mb(-one*tkb(j))*cbxb(j)))*sign(one,cikveb(j))-beamoff(5,imbb(i)))*bbcu(imbb(i),12) &
    )
    yv2(j) = yv2(j) + moidpsv(j)*( &
        ((rkb(j)*(crzb(j)-exp_mb(-one*tkb(j))*cbzb(j)))*sign(one,crkveb(j))-beamoff(4,imbb(i)))*bbcu(imbb(i),12) &
      + ((rkb(j)*(crxb(j)-exp_mb(-one*tkb(j))*cbxb(j)))*sign(one,cikveb(j))-beamoff(5,imbb(i)))*bbcu(imbb(i),11) &
    )
  end if
end do
! end include/beambeam43.f90
