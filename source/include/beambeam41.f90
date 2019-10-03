! start include/beambeam41.f90
if(ibbc == 0) then
  crkveb(1:napx)=(xv1(1:napx)-clobeam(1,imbb(i))) + parbe(ix,5)
  cikveb(1:napx)=(xv2(1:napx)-clobeam(2,imbb(i))) + parbe(ix,6)
else
  crkveb(1:napx) = ((xv1(1:napx)-clobeam(1,imbb(i))) + parbe(ix,5))*bbcu(imbb(i),11) &
                 + ((xv2(1:napx)-clobeam(2,imbb(i))) + parbe(ix,6))*bbcu(imbb(i),12)
  cikveb(1:napx) = ((xv2(1:napx)-clobeam(2,imbb(i))) + parbe(ix,6))*bbcu(imbb(i),11) &
                 - ((xv1(1:napx)-clobeam(1,imbb(i))) + parbe(ix,5))*bbcu(imbb(i),12)
end if
do j=1,napx
  rho2b(j) = crkveb(j)**2+cikveb(j)**2
  if(rho2b(j) <= pieni) cycle
  tkb(j) = rho2b(j)/(two*sigman2(1,imbb(i)))
  if(ibbc == 0) then
    yv1(j) = yv1(j)+moidpsv(j)*(((strack(i)*crkveb(j))/rho2b(j))*(one-exp_mb(-one*tkb(j)))-beamoff(4,imbb(i)))
    yv2(j) = yv2(j)+moidpsv(j)*(((strack(i)*cikveb(j))/rho2b(j))*(one-exp_mb(-one*tkb(j)))-beamoff(5,imbb(i)))
  else
    cccc = (((strack(i)*crkveb(j))/rho2b(j))*(one-exp_mb(-one*tkb(j)))-beamoff(4,imbb(i)))*bbcu(imbb(i),11)&
         - (((strack(i)*cikveb(j))/rho2b(j))*(one-exp_mb(-one*tkb(j)))-beamoff(5,imbb(i)))*bbcu(imbb(i),12)
    yv1(j) = yv1(j)+moidpsv(j)*cccc
    cccc = (((strack(i)*crkveb(j))/rho2b(j))*(one-exp_mb(-one*tkb(j)))-beamoff(4,imbb(i)))*bbcu(imbb(i),12)&
         + (((strack(i)*cikveb(j))/rho2b(j))*(one-exp_mb(-one*tkb(j)))-beamoff(5,imbb(i)))*bbcu(imbb(i),11)
    yv2(j) = yv2(j)+moidpsv(j)*cccc
  end if
end do
! end include/beambeam41.f90
