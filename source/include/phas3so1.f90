! start include/phas3so1.f90
do l=1,2
  ll=2*l
  if(abs(t(ll,ll-1)).gt.pieni) then
    dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
  else
    dphi=pi2-phibf(l)
  end if
  phi(l)=phi(l)+dphi/pie
end do
! end include/phas3so1.f90
