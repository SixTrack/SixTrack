! start include/dalin1.f90
jmel=mel(ix)
if(idp.eq.0.or.ition.eq.0) then
  if(ithick.eq.1) then
    do jb=1,jmel
      jx=mtyp(ix,jb)
      do ip=1,6
        do ien=1,nord+1
          zfeld1(ien)=ald6(jx,1,ip,ien)
        enddo
        if(nvar2.eq.4) then
          call darea6(alda(1,ip),zfeld1,4)
        else if(nvar2.eq.5) then
          call darea6(alda(1,ip),zfeld1,5)
        endif
        do ien=1,nord+1
          zfeld1(ien)=ald6(jx,2,ip,ien)
        enddo
        if(nvar2.eq.4) then
          call darea6(alda(2,ip),zfeld1,4)
        else if(nvar2.eq.5) then
          call darea6(alda(2,ip),zfeld1,5)
        endif
      enddo
! end include/dalin1.f90
