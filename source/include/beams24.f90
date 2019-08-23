! start include/beams24.f90
    endif
    endif
  endif
  goto 290
!--Hirata's 6D beam-beam kick
else if(kzz.eq.20.and.parbe(ix,2).gt.zero) then
  ktrack(i)=44
  parbe(ix,4)=(((-one*crad)*ptnfac(ix))*half)*c1m6
  if(ibeco.eq.1) then
    track6d(1,1)=parbe(ix,5)*c1m3
    track6d(2,1)=zero
    track6d(3,1)=parbe(ix,6)*c1m3
    track6d(4,1)=zero
    track6d(5,1)=zero
    track6d(6,1)=zero
    napx0=napx
    napx=1
    call beamint(napx,track6d,parbe,sigz,bbcu,imbb(i),ix,ibtyp,ibbc,mtc)
    beamoff(1,imbb(i))=track6d(1,1)*c1e3
    beamoff(2,imbb(i))=track6d(3,1)*c1e3
    beamoff(3,imbb(i))=track6d(5,1)*c1e3
    beamoff(4,imbb(i))=track6d(2,1)*c1e3
    beamoff(5,imbb(i))=track6d(4,1)*c1e3
    beamoff(6,imbb(i))=track6d(6,1)
    napx=napx0
  endif
  goto 290
endif
! end include/beams24.f90
