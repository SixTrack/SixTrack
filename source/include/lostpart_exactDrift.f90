! start include/lostpart_exactDrift.f90
  !-----------------------------------------------------------------------
  ! initialise
  llost=.false.
  do j=1,napx
    llostp(j)=.false.
  end do

  !-----------------------------------------------------------------------
  ! check
  if (do_coll) then
    do j=1,napx
      if(part_abs_turn(j).eq.0) then
        llostp(j)=(yv1(j)**2 + yv2(j)**2).gt.c1e6
        llost=llost.or.llostp(j)
      end if
    end do
  else
    do j=1,napx
      llostp(j)=(yv1(j)**2 + yv2(j)**2).gt.c1e6
      llost=llost.or.llostp(j)
    end do
  end if

  !-----------------------------------------------------------------------
  ! dump coordinates in case of losses
  if(llost) then
    ! report losses to user
    call aperture_reportLoss(n,i,ix,.false.)
    if(.not.apflag) call shuffleLostParticles
    ! store old particle coordinates
    ! necessary since aperture markers are downstream of lenses...
    if ( lbacktracking ) call aperture_saveLastCoordinates(i,ix,-1)
  end if

  !-----------------------------------------------------------------------
  ! closing stuff

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
    nnuml = numl
#endif
    return
  end if

  ! stop tracking if no particle survives to this element
  if(nthinerr.ne.0) return
  ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
  ! last modified: 16-07-2018
  if ( lbacktracking ) then
     ! store infos of last aperture marker
     if ( kape(ix).ne.0 ) call aperture_saveLastMarker(i,ix)
  end if
! end include/lostpart_exactDrift.f90
