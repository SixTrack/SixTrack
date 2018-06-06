nfree=nturn1(ix)
if(n.gt.nfree) then
  nac=n-nfree
  pi=four*atan_mb(one)
  !---------ACdipAmp input in Tesla*meter converted to KeV/c
  !---------ejfv(j) should be in MeV/c --> ACdipAmp/ejfv(j) is in mrad
  acdipamp=(ed(ix)*clight)*c1m3                                !hr03
  !---------Qd input in tune units
  qd=ek(ix)
  !---------ACphase input in radians
  acphase=acdipph(ix)
  nramp1=nturn2(ix)
  nplato=nturn3(ix)
  nramp2=nturn4(ix)
  do j=1,napx
#ifndef TILT
    if(nramp1.gt.nac) then
      yv(xory,j)=yv(xory,j)+(((acdipamp*sin_mb(((two*pi)*qd)*real(nac,fPrec)+acphase))&
                *real(nac,fPrec))/real(nramp1,fPrec))/ejfv(j)
    endif
    if(nac.ge.nramp1.and.(nramp1+nplato).gt.nac) then
      yv(xory,j)=yv(xory,j)+(acdipamp*sin_mb(((two*pi)*qd)*real(nac,fPrec)+acphase))/ejfv(j)             !hr03
    endif
    if(nac.ge.(nramp1+nplato).and.(nramp2+nramp1+nplato).gt.nac)then
      yv(xory,j)=yv(xory,j)+((acdipamp*sin_mb(((two*pi)*qd)*real(nac,fPrec)+acphase))&
                *((-one*real(nac-nramp1-nramp2-nplato,fPrec))/real(nramp2,fPrec)))/ejfv(j)      !hr03
    endif
#else
    if (xory.eq.1) then
      acdipamp2=acdipamp*tilts(i)
      acdipamp1=acdipamp*tiltc(i)
    else
      acdipamp2=acdipamp*tiltc(i)
      acdipamp1=-acdipamp*tilts(i)
    endif
    if(nramp1.gt.nac) then
      yv(1,j)=yv(1,j)+(((acdipamp1*sin_mb(((two*pi)*qd)*real(nac,fPrec)+acphase))*real(nac,fPrec))/real(nramp1,fPrec))/ejfv(j)
      yv(2,j)=yv(2,j)+(((acdipamp2*sin_mb(((two*pi)*qd)*real(nac,fPrec)+acphase))*real(nac,fPrec))/real(nramp1,fPrec))/ejfv(j)
    endif
    if(nac.ge.nramp1.and.(nramp1+nplato).gt.nac) then
      yv(1,j)=yv(1,j)+(acdipamp1*sin_mb(((two*pi)*qd)*real(nac,fPrec)+acphase))/ejfv(j)             !hr03
      yv(2,j)=yv(2,j)+(acdipamp2*sin_mb(((two*pi)*qd)*real(nac,fPrec)+acphase))/ejfv(j)             !hr03
    endif
    if(nac.ge.(nramp1+nplato).and.(nramp2+nramp1+nplato).gt.nac)then
      yv(1,j)=yv(1,j)+((acdipamp1*sin_mb(((two*pi)*qd)*real(nac,fPrec)+acphase))&
            *((-one*real(nac-nramp1-nramp2-nplato,fPrec))/real(nramp2,fPrec)))/ejfv(j)      !hr03
      yv(2,j)=yv(2,j)+((acdipamp2*sin_mb(((two*pi)*qd)*real(nac,fPrec)+acphase))&
            *((-one*real(nac-nramp1-nramp2-nplato,fPrec))/real(nramp2,fPrec)))/ejfv(j)      !hr03
    endif
#endif
  enddo
endif
