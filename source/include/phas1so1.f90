! start include/phas1so1.f90
!--solenoid
elseif(kzz.eq.25) then
  do l=1,2
    ll=2*l
    if(abs(t(ll,ll-1)).gt.pieni) then
      phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
    else
      phibf(l)=pi2
    end if
  end do
! end include/phas1so1.f90
