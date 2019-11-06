! start include/acdipkick1.f90
nfree=nturn1(ix)
if(n.gt.nfree) then
  nac=n-nfree
  !---------ACdipAmp input in Tesla*meter converted to KeV/c
  !---------ejfv(j) should be in MeV/c --> ACdipAmp/ejfv(j) is in mrad
  acdipamp=(ed(ix)*clight)*c1m3
  !---------Qd input in tune units
  qd=ek(ix)
  !---------ACphase input in radians
  acphase=acdipph(ix)
  nramp1=nturn2(ix)
  nplato=nturn3(ix)
  nramp2=nturn4(ix)
  do j=1,napx
    acdipamp2=acdipamp*tilts(i)
    acdipamp1=acdipamp*tiltc(i)
    if(nramp1.gt.nac) then
      yv1(j)=yv1(j)+(((acdipamp1*sin_mb((twopi*qd)*real(nac,fPrec)+acphase))*real(nac,fPrec))/real(nramp1,fPrec))/ejfv(j)
      yv2(j)=yv2(j)+(((acdipamp2*sin_mb((twopi*qd)*real(nac,fPrec)+acphase))*real(nac,fPrec))/real(nramp1,fPrec))/ejfv(j)
    endif
    if(nac.ge.nramp1.and.(nramp1+nplato).gt.nac) then
      yv1(j)=yv1(j)+(acdipamp1*sin_mb((twopi*qd)*real(nac,fPrec)+acphase))/ejfv(j)
      yv2(j)=yv2(j)+(acdipamp2*sin_mb((twopi*qd)*real(nac,fPrec)+acphase))/ejfv(j)
    endif
    if(nac.ge.(nramp1+nplato).and.(nramp2+nramp1+nplato).gt.nac)then
      yv1(j)=yv1(j)+((acdipamp1*sin_mb((twopi*qd)*real(nac,fPrec)+acphase))&
            *((-one*real(nac-nramp1-nramp2-nplato,fPrec))/real(nramp2,fPrec)))/ejfv(j)
      yv2(j)=yv2(j)+((acdipamp2*sin_mb((twopi*qd)*real(nac,fPrec)+acphase))&
            *((-one*real(nac-nramp1-nramp2-nplato,fPrec))/real(nramp2,fPrec)))/ejfv(j)
    endif
  enddo
endif
! end include/acdipkick1.f90
