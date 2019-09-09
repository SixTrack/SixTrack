! start include/beams1.f90
!start: beam-beam element
      if(nbeam.ge.1) then
        do 15 i=1,nbb
          nbeaux(i)=0
   15   continue
        do i=1,iu
          ix=ic(i)
          if(ix.gt.nblo) then
            ix=ix-nblo
            if(kz(ix).eq.20.and.parbe(ix,2).eq.zero) then
!--round beam
              if(sigman(1,imbb(i)).eq.sigman(2,imbb(i))) then
                if(nbeaux(imbb(i)).eq.2.or.nbeaux(imbb(i)).eq.3) then
                  write(lerr,"(a)") "BEAMBEAM> ERROR At each interaction point the beam must be either "//&
                    "round or elliptical for all particles"
                  call prror
                else
                  nbeaux(imbb(i))=1
                  sigman2(1,imbb(i))=sigman(1,imbb(i))**2
                endif
              endif
!--elliptic beam x>z
              if(sigman(1,imbb(i)).gt.sigman(2,imbb(i))) then
                if(nbeaux(imbb(i)).eq.1.or.nbeaux(imbb(i)).eq.3) then
                  write(lerr,"(a)") "BEAMBEAM> ERROR At each interaction point the beam must be either "//&
                    "round or elliptical for all particles"
                  call prror
                else
                  nbeaux(imbb(i))=2
                  sigman2(1,imbb(i))=sigman(1,imbb(i))**2
                  sigman2(2,imbb(i))=sigman(2,imbb(i))**2
                  sigmanq(1,imbb(i))=sigman(1,imbb(i))/sigman(2,imbb(i))
                  sigmanq(2,imbb(i))=sigman(2,imbb(i))/sigman(1,imbb(i))
                endif
              endif
!--elliptic beam z>x
              if(sigman(1,imbb(i)).lt.sigman(2,imbb(i))) then
                if(nbeaux(imbb(i)).eq.1.or.nbeaux(imbb(i)).eq.2) then
                  write(lerr,"(a)") "BEAMBEAM> ERROR At each interaction point the beam must be either "//&
                    "round or elliptical for all particles"
                  call prror
                else
                  nbeaux(imbb(i))=3
                  sigman2(1,imbb(i))=sigman(1,imbb(i))**2
                  sigman2(2,imbb(i))=sigman(2,imbb(i))**2
                  sigmanq(1,imbb(i))=sigman(1,imbb(i))/sigman(2,imbb(i))
                  sigmanq(2,imbb(i))=sigman(2,imbb(i))/sigman(1,imbb(i))
                endif
              endif
            endif
          endif
        enddo
      endif
!end: beam-beam element
! end include/beams1.f90
