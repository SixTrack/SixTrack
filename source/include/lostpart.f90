  !-----------------------------------------------------------------------
  ! initialise
  llost=.false.
  do j=1,napx
    llostp(j)=.false.
  end do

  !-----------------------------------------------------------------------
  ! check
  call aperture_checkApeMarker( n, i, ix, llost, nthinerr )
  
!   ! any particle loss?
!   do j=1,napx
!     llost=llost.or.llostp(j)
!     if (llost) exit
!   end do

  !-----------------------------------------------------------------------
  ! dump coordinates in case of losses
  if(llost) then
    ! report losses to user
    call aperture_reportLoss( n, i, ix, llost, nthinerr )
    call compactArrays
    ! store old particle coordinates
    ! necessary since aperture markers are downstream of lenses...
    if ( lbacktracking ) call aperture_saveLastCoordinates(i,ix,-1)
  end if

  !-----------------------------------------------------------------------
  ! closing stuff

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

  ! stop tracking if no particle survives to this element
  if(nthinerr.ne.0) return
  ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
  ! last modified: 16-07-2018
  if ( lbacktracking ) then
     ! store infos of last aperture marker
     if ( kape(ix).ne.0 ) call aperture_saveLastMarker(i,ix)
  end if
