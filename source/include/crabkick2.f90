! start include/crabkick.f90
!---------CrabAmp input in MV
!---------ejfv(j) should be in MeV/c --> CrabAmp/c/ejfv(j) is in rad
!---------ejfv(j) should be in MeV ?? --> CrabAmp/ejfv(j) is in rad
!---------CrabFreq input in MHz (ek)
!---------sigmv should be in mm --> sigmv*1e-3/clight*ek*1e6 in rad
crabfreq=ek(ix)*c1e3

do j=1,napx ! loop over particles
  crabamp=ed(ix)*nqq(j)
  kcrab=(((sigmv(j)/(clight*(e0f/e0)))*crabfreq)*two)*pi+crabph(ix)
  yv2(j)=yv2(j) - (crabamp*c1e3)*sin_mb(kcrab)*(moidpsv(j)/e0f)
  ejv(j)=ejv(j) - (((((crabamp*crabfreq)*two)*pi)/clight)*xv2(j))*cos_mb(kcrab)
end do
call part_updatePartEnergy(1,.true.)
if(ithick == 1) call envarsv
! end include/crabkick.f90
