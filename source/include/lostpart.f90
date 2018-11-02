  !-----------------------------------------------------------------------
  ! initialise
  llost=.false.
  do j=1,napx
    llostp(j)=.false.
  end do
  
  !-----------------------------------------------------------------------
  ! check against current aperture marker
  if(.not.limifound.or.kape(ix).eq.0) then
    ! limi block not there or aperture type not assigned
    ! general check (set in the ITER block)
    do j=1,napx
      if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll)) then
        llostp(j)=(abs(xv(1,j)).gt.aper(1)).or.(abs(xv(2,j)).gt.aper(2)).or. &
             (xv(1,j).ne.xv(1,j)).or.(xv(2,j).ne.xv(2,j))
      end if
    end do
  else
    ! go through all possible types of aperture
    call aperture_checkApeMarker(n, i, ix)
  end if

  ! any particle loss?
  do j=1,napx
    if (llostp(j)) then
      llost=.true.
      exit
    end if
  end do
  
  !-----------------------------------------------------------------------
  ! dump coordinates in case of losses
  if(llost) then
    ! report losses to user
    call aperture_reportLoss(n, i, ix)
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
