! start include/beams21.f90
!--beam-beam element
        if(kzz.eq.20.and.nbeam.ge.1.and.parbe(ix,2).eq.zero) then
            strack(i)=crad*ptnfac(ix)
            if(abs(strack(i)).le.pieni) then
              ktrack(i)=31
              goto 290
            endif
            if(nbeaux(imbb(i)).eq.1) then
              ktrack(i)=41
              if(ibeco.eq.1) then
                do 42 j=1,napx
! end include/beams21.f90
