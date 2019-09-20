! ============================================================================ !
!  Tracking Module
! ~~~~~~~~~~~~~~~~~
!  Module holding various routines used by the main tracking routines
! ============================================================================ !
module tracking

  use floatPrecision

  implicit none

  logical,           public,  save :: tr_is4D   = .false.
  character(len=8),  private, save :: trackMode = " "
  character(len=32), private, save :: trackFmt  = " "

contains

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-20
!  Updated: 2019-09-20
!  Initialisation of variables needed for tracking
! ================================================================================================ !
subroutine trackInit

  use mod_common
  use mathlib_bouncer

  integer oPart, oTurn

  ! Check whether this is 4D or 6D
  tr_is4D = idp == 0 .or. ition == 0

  ! Define the label of the tracking mode
  if(ithick == 1) then
    trackMode = "Thick"
  else
    trackMode = "Thin"
  end if
  if(tr_is4D) then
    trackMode = trim(trackMode)//" 4D"
  else
    trackMode = trim(trackMode)//" 6D"
  end if
  oPart = int(log10_mb(real(napxo, kind=fPrec))) + 1
  oTurn = int(log10_mb(real(numl,  kind=fPrec))) + 1
  write(trackFmt,"(2(a,i0),a)") "(2(a,i",oTurn,"),2(a,i",oPart,"))"

end subroutine trackInit

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Updated: 2019-09-20
! ================================================================================================ !
subroutine trackReport(n)

  use crcoall
  use mod_common, only : numl, napx, napxo

  integer, intent(in) :: n

  write(lout,trackFmt) "TRACKING> "//trim(trackMode)//": Turn ",n," / ",numl,", Particles: ",napx," / ",napxo
  flush(lout)

end subroutine trackReport

! ================================================================================================ !
!  Prepare for Tracking
!  Code merged from old trauthin and trauthick routines
! ================================================================================================ !
subroutine preTracking

  use crcoall
  use mod_common
  use mod_common_track
  use numerical_constants

  use collimation, only : do_coll
  use mod_fluc,    only : fluc_writeFort4

  integer i, ix, jb, jx, kpz
  logical isThick

  isThick = ithick == 1

  if(isThick .and. do_coll) then
    write(lerr,"(a)") "TRACKING> ERROR Collimation is not supported for thick tracking"
    call prror
  end if

  if(mout2 == 1) call fluc_writeFort4

  return

  ! BEGIN Loop over structure elements
  do i=1,iu

    ! Get single and block element index
    ix = ic(i)
    if(ix <= nblo) then
      ! BLOC element
      ktrack(i) = 1
      do jb=1,mel(ix)
        jx        = mtyp(ix,jb)
        strack(i) = strack(i) + el(jx)
      end do
      if(abs(strack(i)) <= pieni) then
        ktrack(i) = 31
      end if
      cycle
    end if

    ix  = ix-nblo     ! Single element index
    kpz = abs(kp(ix)) ! Element type
  end do
  ! END Loop over structure elements

end subroutine preTracking

end module tracking
