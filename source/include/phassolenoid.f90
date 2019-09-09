! start include/phassolenoid.f90
do l=1,2
  ll = 2*l
  if(abs(t(ll,ll-1)) > pieni) then
    phibf(l) = atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
  else
    phibf(l) = pi2
  end if
end do

crkve  = t(i,2) - (t(i,1)*qu)*qv
cikve  = t(i,4) - (t(i,3)*qu)*qv
t(i,2) = crkve*cos_mb(qv) + cikve*sin_mb(qv)
t(i,4) = cikve*cos_mb(qv) - crkve*sin_mb(qv)
crkve  = t(i,1)*cos_mb(qv) + t(i,3)*sin_mb(qv)
cikve  = t(i,3)*cos_mb(qv) - t(i,1)*sin_mb(qv)
t(i,1) = crkve
t(i,3) = cikve

do l=1,2
  ll = 2*l
  if(abs(t(ll,ll-1)) > pieni) then
    dphi = atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
  else
    dphi = pi2-phibf(l)
  end if
  phi(l) = phi(l)+dphi/twopi
end do
! end include/phassolenoid.f90
