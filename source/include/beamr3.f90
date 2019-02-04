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
! end include/beamr3.f90
