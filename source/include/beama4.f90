! start include/beama4.f90
  if(ibbc.eq.0) then
    yv1(j)=yv1(j)+moidpsv(j)*((rkb(j)*(crzb(j)-exp_mb(-one*tkb(j))*cbzb(j)))*sign(one,crkveb(j))-beamoff(4,imbb(i)))
    yv2(j)=yv2(j)+moidpsv(j)*((rkb(j)*(crxb(j)-exp_mb(-one*tkb(j))*cbxb(j)))*sign(one,cikveb(j))-beamoff(5,imbb(i)))
  else
    cccc=((rkb(j)*(crzb(j)-exp_mb(-one*tkb(j))*cbzb(j)))*sign(one,crkveb(j))-beamoff(4,imbb(i)))*bbcu(imbb(i),11)&
        -((rkb(j)*(crxb(j)-exp_mb(-one*tkb(j))*cbxb(j)))*sign(one,cikveb(j))-beamoff(5,imbb(i)))*bbcu(imbb(i),12)
    yv1(j)=yv1(j)+moidpsv(j)*cccc
    cccc=((rkb(j)*(crzb(j)-exp_mb(-one*tkb(j))*cbzb(j)))*sign(one,crkveb(j))-beamoff(4,imbb(i)))*bbcu(imbb(i),12)&
        +((rkb(j)*(crxb(j)-exp_mb(-one*tkb(j))*cbxb(j)))*sign(one,cikveb(j))-beamoff(5,imbb(i)))*bbcu(imbb(i),11)
    yv2(j)=yv2(j)+moidpsv(j)*cccc
  end if
end do
! end include/beama4.f90
