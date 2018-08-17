! start include/dalin4.f90
else
  do jb=1,jmel
    jx=mtyp(ix,jb)
    if(ithick.eq.1) then
      do ip=1,6
        do ien=1,nord+1
          zfeld1(ien)=ald6(jx,1,ip,ien)
          zfeld2(ien)=asd6(jx,1,ip,ien)
        enddo
        if(nvar2.eq.5) then
          call darea6(alda(1,ip),zfeld1,5)
          call darea6(asda(1,ip),zfeld2,5)
        else if(nvar2.eq.6) then
          call darea6(alda(1,ip),zfeld1,6)
          call darea6(asda(1,ip),zfeld2,6)
        endif
        do ien=1,nord+1
          zfeld1(ien)=ald6(jx,2,ip,ien)
          zfeld2(ien)=asd6(jx,2,ip,ien)
        enddo
        if(nvar2.eq.5) then
          call darea6(alda(2,ip),zfeld1,5)
          call darea6(asda(2,ip),zfeld2,5)
        else if(nvar2.eq.6) then
          call darea6(alda(2,ip),zfeld1,6)
          call darea6(asda(2,ip),zfeld2,6)
        endif
      enddo
! end include/dalin4.f90
