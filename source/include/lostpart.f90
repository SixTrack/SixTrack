  do j=1,napx
    llostp(j)=.false.
  end do

  
  !-----------------------------------------------------------------------
  ! check against current aperture marker
  !-----------------------------------------------------------------------

  llost=.false.
  lback=.false.

  if(.not.limifound.or.kape(ix).eq.0) then
    ! limi block not there or aperture type not assigned
    ! general check (set in the ITER block)
    do j=1,napx

      if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll)) then
        llostp(j)=(abs(xv(1,j)).gt.aper(1)).or.(abs(xv(2,j)).gt.aper(2)).or. &
             (xv(1,j).ne.xv(1,j)).or.(xv(2,j).ne.xv(2,j))
        llost=llost.or.llostp(j)
      else if (do_coll) then
        llostp(j)=.false.
      end if
    end do

  else

    ! go through all possible types of aperture
    select case(kape(ix))

    case (-1) ! Transition
      apxx = ape(3,ix)**2.
      apyy = ape(4,ix)**2.
      apxy = apxx * apyy
      do j=1,napx
        if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll)) then
          if(lapeofftlt(ix)) then
            if(lbacktracking) then
              call roffpos(xLast(1,j),xLast(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
            else
              call roffpos(xv(1,j),xv(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
            end if
            llostp(j)=checkTR(xchk(1),xchk(2),ape(1,ix),ape(2,ix),ape(3,ix),ape(4,ix),apxx,apyy,apxy,ape(5,ix),ape(6,ix)).or. &
              isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
          else
            if(lbacktracking) then
              llostp(j)= &
                checkTR(xLast(1,j),xLast(2,j),ape(1,ix),ape(2,ix),ape(3,ix),ape(4,ix),apxx,apyy,apxy,ape(5,ix),ape(6,ix)) .or. &
                isnan_mb(xLast(1,j)).or.isnan_mb(xLast(2,j))
            else
              llostp(j)= &
                checkTR(xv(1,j),xv(2,j),ape(1,ix),ape(2,ix),ape(3,ix),ape(4,ix),apxx,apyy,apxy,ape(5,ix),ape(6,ix))       .or. &
                isnan_mb(xv(1,j)).or.isnan_mb(xv(2,j))
            end if
          end if
          llost=llost.or.llostp(j)

        else if (do_coll) then
          llostp(j)=.false.
        end if
      end do

    case (1) ! circle
      radius2 = ape(3,ix)**2
      do j=1,napx

        if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll)) then

          if(lapeofftlt(ix)) then
            if(lbacktracking) then
              call roffpos(xLast(1,j),xLast(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
            else
              call roffpos(xv(1,j),xv(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
            end if
            llostp(j)=checkCR( xchk(1),xchk(2),radius2 ) .or. &
              isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
          else
            if(lbacktracking) then
              llostp(j)=checkCR( xLast(1,j),xLast(2,j),radius2 ) .or. &
                isnan_mb(xLast(1,j)).or.isnan_mb(xLast(2,j))
            else
              llostp(j)=checkCR( xv(1,j),xv(2,j),radius2 ) .or. &
                isnan_mb(xv(1,j)).or.isnan_mb(xv(2,j))
            end if
          end if
          llost=llost.or.llostp(j)

        else if (do_coll) then
          llostp(j)=.false.
        end if

      end do

    case (2) ! Rectangle
      do j=1,napx
        if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll) ) then

          if(lapeofftlt(ix)) then
            if(lbacktracking) then
              call roffpos(xLast(1,j),xLast(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
            else
              call roffpos(xv(1,j),xv(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
            end if
            llostp(j)=checkRE( xchk(1),xchk(2),ape(1,ix),ape(2,ix) ) .or. &
              isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
          else
            if(lbacktracking) then
              llostp(j)=checkRE( xLast(1,j),xLast(2,j),ape(1,ix),ape(2,ix) ) .or. &
                isnan_mb(xLast(1,j)).or.isnan_mb(xLast(2,j))
            else
              llostp(j)=checkRE( xv(1,j),xv(2,j),ape(1,ix),ape(2,ix) ) .or. &
                isnan_mb(xv(1,j)).or.isnan_mb(xv(2,j))
            end if
          end if
          llost=llost.or.llostp(j)

        else if (do_coll) then
          llostp(j)=.false.
        end if

      end do

    case (3) ! Ellipse
      apxx = ape(3,ix)**2.
      apyy = ape(4,ix)**2.
      apxy = apxx * apyy
      do j=1,napx
        if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll)) then
          if(lapeofftlt(ix)) then
            if(lbacktracking) then
              call roffpos(xLast(1,j),xLast(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
            else
              call roffpos(xv(1,j),xv(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
            end if
            llostp(j)=checkEL( xchk(1),xchk(2),apxx,apyy,apxy ) .or. &
              isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
          else
            if(lbacktracking) then
              llostp(j)=checkEL( xLast(1,j),xLast(2,j),apxx,apyy,apxy ) .or. &
                isnan_mb(xLast(1,j)).or.isnan_mb(xLast(2,j))
            else
              llostp(j)=checkEL( xv(1,j),xv(2,j),apxx,apyy,apxy ) .or. &
                isnan_mb(xv(1,j)).or.isnan_mb(xv(2,j))
            end if
          end if
          llost=llost.or.llostp(j)

        else if (do_coll) then
          llostp(j)=.false.
        end if

      end do

    case (4) ! RectEllipse
      apxx = ape(3,ix)**2.
      apyy = ape(4,ix)**2.
      apxy = apxx * apyy
      do j=1,napx
        if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll)) then
          if(lapeofftlt(ix)) then
            if(lbacktracking) then
              call roffpos(xLast(1,j),xLast(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
            else
              call roffpos(xv(1,j),xv(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
            end if
            llostp(j)=checkRL( xchk(1),xchk(2),ape(1,ix),ape(2,ix),apxx,apyy,apxy ) .or. &
              isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
          else
            if(lbacktracking) then
              llostp(j)=checkRL( xLast(1,j),xLast(2,j),ape(1,ix),ape(2,ix),apxx,apyy,apxy ) .or. &
                isnan_mb(xLast(1,j)).or.isnan_mb(xLast(2,j))
            else
              llostp(j)=checkRL( xv(1,j),xv(2,j),ape(1,ix),ape(2,ix),apxx,apyy,apxy ) .or. &
                isnan_mb(xv(1,j)).or.isnan_mb(xv(2,j))
            end if
          end if
          llost=llost.or.llostp(j)

        else if (do_coll) then
          llostp(j)=.false.
        end if
      end do

    case (5) ! Octagon
      do j=1,napx
        if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll)) then

          if(lapeofftlt(ix)) then
            if(lbacktracking) then
              call roffpos(xLast(1,j),xLast(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
            else
              call roffpos(xv(1,j),xv(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
            end if
            llostp(j)=checkOC(xchk(1),xchk(2),ape(1,ix),ape(2,ix),ape(5,ix),ape(6,ix)).or. &
              isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
          else
            if(lbacktracking) then
              llostp(j)=checkOC(xLast(1,j),xLast(2,j),ape(1,ix),ape(2,ix),ape(5,ix),ape(6,ix)).or. &
                isnan_mb(xLast(1,j)).or.isnan_mb(xLast(2,j))
            else
              llostp(j)=checkOC(xv(1,j),xv(2,j),ape(1,ix),ape(2,ix),ape(5,ix),ape(6,ix)).or. &
                isnan_mb(xv(1,j)).or.isnan_mb(xv(2,j))
            end if
          end if
          llost=llost.or.llostp(j)

        else if (do_coll) then
          llostp(j)=.false.
        end if

      end do

    case (6) ! Racetrack
      !   NB: it follows the MadX definition
      apxy = ape(3,ix)**2.
      do j=1,napx
        if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll)) then
          if(lapeofftlt(ix)) then
            if(lbacktracking) then
              call roffpos(xLast(1,j),xLast(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
            else
              call roffpos(xv(1,j),xv(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
            end if
            llostp(j)=checkRT(xchk(1),xchk(2),ape(1,ix),ape(2,ix),ape(3,ix),apxy).or. &
              isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
          else
            if(lbacktracking) then
              llostp(j)=checkRT(xLast(1,j),xLast(2,j),ape(1,ix),ape(2,ix),ape(3,ix),apxy).or. &
                isnan_mb(xLast(1,j)).or.isnan_mb(xLast(2,j))
            else
              llostp(j)=checkRT(xv(1,j),xv(2,j),ape(1,ix),ape(2,ix),ape(3,ix),apxy).or. &
                isnan_mb(xv(1,j)).or.isnan_mb(xv(2,j))
            end if
          end if
          llost=llost.or.llostp(j)

        else if (do_coll) then
          llostp(j)=.false.
        end if
      end do

    end select
  end if ! if(.not.limifound.or.kape(ix).eq.0)

  !-----------------------------------------------------------------------
  ! dump coordinates in case of losses
  ! if back-tracking is requested, get more detailed point of loss
  ! for the moment, only bi-section method
  !-----------------------------------------------------------------------

  if(llost) then

    if(lbacktracking.and.kape(ix).ne.0.and.iBckTypeLast.ge.0) then
      lback=.true.

      ! Length between elements
      length = dcum(i) - dcum(iLast)

      ! - pay attention to overflow:
      if( length .lt. zero ) then
        length = length+tlen
      end if

      ! - pay attention to too short thick elements
      if( length .le. bktpre ) then
        lback=.false.
      end if

    end if

    ! Number of iterations for bisection method (ln(2x/precision)/ln(2)+1)
    if(lback) then
      niter=nint(inv_ln2*log_mb(two*length/bktpre)+2)
    end if

    do j=1,napx
      if(llostp(j)) then
        ! treat a lost particle

        ! ==============================================================
        ! point of loss
        if(lback) then
          ! A. Mereghetti and P. Garcia Ortega, for the FLUKA Team
          ! last modified: 21-03-2018
          ! back-track particles, in order to better estimate actual loss point

          ylos(1)=yLast(1,j)
          ylos(2)=yLast(2,j)

          ! actual algorithm
          llos = llostp(j)
          step = one

          do jj=1,niter
            ! current step (bisection method):
            if( llos ) then
              step = step - one / (two**(jj))
            else
              step = step + one / (two**(jj))
            end if

            ! - step discretized if last iteration, to compare with BeamLossPattern
            if(jj.eq.niter) then
              slos = int((dcum(iLast)+length*step)/bktpre+one)*bktpre
              step = (slos-dcum(iLast))/length
            end if

            ! - particle coordinates at current step
            select case(iBckTypeLast)
            case (0)
              ! back-track along a drift
              xlos(1) = xLast(1,j)  - yLast(1,j)*((one-step)*length)
              xlos(2) = xLast(2,j)  - yLast(2,j)*((one-step)*length)
              slos    = dcum(iLast) + (step*length)
            end select

            ! - aperture at current step
            call interp_aperture( iLast, ixLast, i, ix, kapert, aprr, slos )

            ! Check aperture
            if( lapeofftlt(ix).or.lapeofftlt(ixLast) ) then
              call roffpos( xlos(1), xlos(2), xchk(1),xchk(2), aprr(7), aprr(8), aprr(9) )
            else
              xchk(1) = xlos(1)
              xchk(2) = xlos(2)
            end if

            select case(kapert)
            case(-1) ! Transition
              apxx = aprr(3)**2.
              apyy = aprr(4)**2.
              apxy = apxx * apyy
              llos=checkTR(xchk(1),xchk(2),aprr(1),aprr(2),aprr(3),aprr(4),apxx,apyy,apxy,aprr(5),aprr(6)).or. &
                isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
            case (1) ! Circle
              radius2 = aprr(3)**2
              llos=checkCR(xchk(1),xchk(2),radius2) .or. &
                isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
            case (2) ! Rectangle
              llos=checkRE(xchk(1),xchk(2),aprr(1),aprr(2)) .or. &
                isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
            case (3) ! Ellipse
              apxx = aprr(3)**2.
              apyy = aprr(4)**2.
              apxy = apxx * apyy
              llos=checkEL( xchk(1),xchk(2),apxx,apyy,apxy )  .or. &
                isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
            case (4) ! RectEllipse
              apxx = aprr(3)**2.
              apyy = aprr(4)**2.
              apxy = apxx * apyy
              llos = checkRL( xchk(1),xchk(2),aprr(1),aprr(2),apxx, apyy, apxy ) .or. &
                isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
            case (5) ! Octagon
              llos=checkOC(xchk(1), xchk(2), aprr(1), aprr(2), aprr(5), aprr(6) ) .or. &
                isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
            case (6) ! RaceTrack
              llos=checkRT( xchk(1), xchk(2), aprr(1), aprr(2), aprr(3), aprr(3)**2. ) .or. &
                isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
            end select
          end do !do jj=1,niter

          ! pay attention to overflow
          if( slos.gt.tlen ) then
            slos=slos-tlen
          end if

        else !if(lback)
          if(lbacktracking) then
            xlos(1) = xLast(1,j)
            xlos(2) = xLast(2,j)
            ylos(1) = yLast(1,j)
            ylos(2) = yLast(2,j)
            slos    = dcum(iLastThick)
          else
            xlos(1) = xv(1,j)
            xlos(2) = xv(2,j)
            ylos(1) = yv(1,j)
            ylos(2) = yv(2,j)
            slos    = dcum(i)
          end if
        end if ! if(lback)

        ! get ready for dumping infos
        if(lbacktracking) then
          ejfvlos = ejfvLast(j)
          ejvlos = ejvLast(j)
          nucmlos = nucmLast(j)
          sigmvlos = sigmv(j)
          dpsvlos = dpsvLast(j)
          naalos = naaLast(j)
          nzzlos = nzzLast(j)
        else
          ejfvlos = ejfv(j)
          ejvlos = ejv(j)
          nucmlos = nucm(j)
          sigmvlos = sigmv(j)
          dpsvlos = dpsv(j)
          naalos = naa(j)
          nzzlos = nzz(j)
        end if

        ! ==============================================================
        ! If lost particles aren't killed, the lost info is dumped only
        ! the first time they hit the aperture. Their secondaries generated
        ! from a lost particles are considered lost as well
        if( apflag ) then
          lparID = .false.
          jjx=1

          !TODO is this really needed?
          if (do_coll) then
            npart_tmp = npart
          else
            npart_tmp = napx
          endif

          do jj=1,npart_tmp
            if(plost(jj).ne.0) then
#ifdef FLUKA
              if( fluka_uid(j).eq.plost(jj).or. fluka_gen(j).eq.plost(jj) ) then
#else
              if ( (     do_coll .and. (  ipart(j) .eq. plost(jj) )) .or. &
                   (.not.do_coll .and. ( nlostp(j) .eq. plost(jj) ))       ) then
#endif
                lparID=.true.
              end if

              jjx=jj+1 !points to the last zero
            end if
          end do

          if(lparID) then
            !old lost particle or secondary, don't print it
            goto 1982
          else
            !new lost particle, store ID and print it
#ifdef FLUKA
            plost(jjx) = fluka_uid(j)
#else
            if (do_coll) then
              plost(jjx) = ipart(j)
            else
              plost(jjx) = j
            endif
#endif
          end if !if(lparID) then
        end if !if( apflag ) then

#ifdef HDF5
        if(h5_useForAPER) then
          call h5_prepareWrite(aper_setLostPart, 1)
          call h5_writeData(aper_setLostPart, 1,  1, n)
          call h5_writeData(aper_setLostPart, 2,  1, i)
          call h5_writeData(aper_setLostPart, 3,  1, ix)
          call h5_writeData(aper_setLostPart, 4,  1, bez(ix))
          call h5_writeData(aper_setLostPart, 5,  1, slos)
          call h5_writeData(aper_setLostPart, 6,  1, xlos(1)*c1m3)
          call h5_writeData(aper_setLostPart, 7,  1, xlos(2)*c1m3)
          call h5_writeData(aper_setLostPart, 8,  1, ylos(1)*c1m3)
          call h5_writeData(aper_setLostPart, 9,  1, ylos(2)*c1m3)
          call h5_writeData(aper_setLostPart, 10, 1, ejfvlos*c1m3)
          call h5_writeData(aper_setLostPart, 11, 1, (ejvlos*(nucm0/nucmlos)-e0)*c1e6)
          call h5_writeData(aper_setLostPart, 12, 1, -c1m3 * (sigmvlos/clight) * (e0/e0f))
          call h5_writeData(aper_setLostPart, 13, 1, naalos)
          call h5_writeData(aper_setLostPart, 14, 1, nzzlos)
#ifdef FLUKA
          call h5_writeData(aper_setLostPart, 15, 1, fluka_uid(j))
          call h5_writeData(aper_setLostPart, 16, 1, fluka_gen(j))
          call h5_writeData(aper_setLostPart, 17, 1, fluka_weight(j))
#endif
          if (do_coll) then
            call h5_writeData(aper_setLostPart, 15, 1, ipart(j))
          endif
#ifndef FLUKA
          if (.not. do_coll) then
            call h5_writeData(aper_setLostPart, 15, 1, nlostp(j))
          endif
#endif
          call h5_finaliseWrite(aper_setLostPart)
        else
  ! END of #ifdef HDF5
#endif

          ! Print to unit 999 (fort.999)
#ifdef FLUKA
          write(losses_unit,'(3(1X,I8),1X,A48,1X,F12.5,2(1X,I8),8(1X,1PE14.7),2(1X,I8))')&
#else
          write(losses_unit,'(3(1X,I8),1X,A48,1X,F12.5,1X,I8,7(1X,1PE14.7),2(1X,I8))')   &
#endif

     &         n, i, ix, bez(ix), slos,                                     &
#ifdef FLUKA
     &         fluka_uid(j), fluka_gen(j), fluka_weight(j),                    &
#else
     &         nlostp(j),                                                      &
#endif

     &         xlos(1)*c1m3, ylos(1)*c1m3, xlos(2)*c1m3, ylos(2)*c1m3,         &
     &         ejfvlos*c1m3, (ejvlos*(nucm0/nucmlos)-e0)*c1e6,                 &
     &         -c1m3 * (sigmvlos/clight) * (e0/e0f),                           &
     &         naalos, nzzlos
#ifdef HDF5
        end if
#endif

#if defined(ROOT)
! root output
        if(root_flag .and. root_ApertureCheck.eq.1) then
          this_name = trim(adjustl(bez(ix))) // C_NULL_CHAR
#if defined(FLUKA)
          call ApertureCheckWriteLossParticleF(n, i, ix, this_name, len_trim(this_name), slos, &
       &  fluka_uid(j), fluka_gen(j), fluka_weight(j), &
       &  xlos(1)*c1m3, ylos(1)*c1m3, xlos(2)*c1m3, ylos(2)*c1m3, ejfvlos*c1m3, (ejvlos-e0)*c1e6, &
       &  -c1m3 * (sigmvlos/clight) * (e0/e0f), naalos, nzzlos)
#else
          call ApertureCheckWriteLossParticle(n, i, ix, this_name, len_trim(this_name), slos, plost(j),&
       &  xlos(1)*c1m3, ylos(1)*c1m3, xlos(2)*c1m3, ylos(2)*c1m3, ejfvlos*c1m3, (ejvlos-e0)*c1e6, &
       &  -c1m3 * (sigmvlos/clight) * (e0/e0f), naalos, nzzlos)
#endif
        end if
#endif

#ifdef FLUKA
        if(nlostp(j).le.aperture_napxStart) then
#else
        if(((nlostp(j).le.aperture_napxStart) .and. do_coll) &
             .or. .not.do_coll) then
#endif
           pstop(nlostp(j))=.true.
           ! Record for postpr
           if(.not.limifound.or.kape(ix).eq.0) then
             aperv(nlostp(j),1) = aper(1)
             aperv(nlostp(j),2) = aper(2)
           else
             aperv(nlostp(j),1) = min(ape(1,ix),ape(3,ix))
             aperv(nlostp(j),2) = min(ape(2,ix),ape(4,ix))
           end if
           ixv(nlostp(j))     = ix
           xvl(1,nlostp(j))   = xlos(1)
           xvl(2,nlostp(j))   = xlos(2)
           yvl(1,nlostp(j))   = ylos(1)
           yvl(2,nlostp(j))   = ylos(2)
           dpsvl(nlostp(j))   = dpsvlos
           ejvl(nlostp(j))    = ejvlos
           sigmvl(nlostp(j))  = sigmvlos
           numxv(nlostp(j))   = numx
           nnumxv(nlostp(j))  = numx

        end if !  (nlostp(j).le.aperture_napxStart) OR
               ! ((nlostp(j).le.aperture_napxStart) .and. do_coll)

1982    continue

      end if ! if(llostp(j))
    end do ! do j=1,napx

    ! flush loss particle file
#ifdef HDF5
    if(.not. h5_useForAPER) then
#endif
       flush(losses_unit)
#ifdef HDF5
    end if
#endif

    call compactArrays

    ! store old particle coordinates
    ! necessary since aperture markers are downstream of lenses...
    if ( lbacktracking ) call aperture_saveLastCoordinates(i,ix,-1)

  end if !if( llost ) then

  !-----------------------------------------------------------------------
  ! closing stuff
  !-----------------------------------------------------------------------

#ifdef FLUKA
  napxo = napx
#endif

  if(napx.eq.0) then
    write(lout,"(a)") ""
    write(lout,"(a)") "************************"
    write(lout,"(a)") "** ALL PARTICLES LOST **"
    write(lout,"(a)") "**   PROGRAM STOPS    **"
    write(lout,"(a)") "************************"
    write(lout,"(a)") ""
#ifdef FLUKA
!skip postpr
    nthinerr = 3000
#else
    nthinerr = 3001
    nnuml=numl
#endif
    return
  end if
