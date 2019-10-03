! start include/beamco.f90
do j=1,napx
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
  ! start include/beamr1.f90
  rho2b(j)=crkveb(j)**2+cikveb(j)**2
  if(rho2b(j).le.pieni) cycle
  ! end include/beamr1.f90
  ! start include/beamr2.f90
  tkb(j)=rho2b(j)/(two*sigman2(1,imbb(i)))
  ! end include/beamr2.f90
  ! start include/beamr3.f90
  if(ibbc.eq.0) then
    yv1(j)=yv1(j)+moidpsv(j)*(((strack(i)*crkveb(j))/rho2b(j))*(one-exp_mb(-one*tkb(j)))-beamoff(4,imbb(i)))
    yv2(j)=yv2(j)+moidpsv(j)*(((strack(i)*cikveb(j))/rho2b(j))*(one-exp_mb(-one*tkb(j)))-beamoff(5,imbb(i)))
  else
    cccc=(((strack(i)*crkveb(j))/rho2b(j))*(one-exp_mb(-one*tkb(j)))-beamoff(4,imbb(i)))*bbcu(imbb(i),11)&
        -(((strack(i)*cikveb(j))/rho2b(j))*(one-exp_mb(-one*tkb(j)))-beamoff(5,imbb(i)))*bbcu(imbb(i),12)
    yv1(j)=yv1(j)+moidpsv(j)*cccc
    cccc=(((strack(i)*crkveb(j))/rho2b(j))*(one-exp_mb(-one*tkb(j)))-beamoff(4,imbb(i)))*bbcu(imbb(i),12)&
        +(((strack(i)*cikveb(j))/rho2b(j))*(one-exp_mb(-one*tkb(j)))-beamoff(5,imbb(i)))*bbcu(imbb(i),11)
    yv2(j)=yv2(j)+moidpsv(j)*cccc
  end if
end do
! end include/beamr3.f90
