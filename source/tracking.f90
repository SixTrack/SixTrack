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
  if(iclo6 > 0) then
    trackMode = trim(trackMode)//" 6D"
  else
    trackMode = trim(trackMode)//" 4D"
  end if
  oPart = int(log10_mb(real(napxo, kind=fPrec))) + 1
  oTurn = int(log10_mb(real(numl,  kind=fPrec))) + 1
  write(trackFmt,"(2(a,i0),a)") "(2(a,i",oTurn,"),2(a,i",oPart,"))"

end subroutine trackInit

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Write a turn report.
!  The isFirst if statement is only computed the first time the routine is called.
! ================================================================================================ !
subroutine trackReport(n)
  
  use crcoall
  use mod_common, only : numl, napx, napxo

  integer, intent(in) :: n

  write(lout,trackFmt) "TRACKING> "//trim(trackMode)//": Turn ",n," / ",numl,", Particles: ",napx," / ",napxo
  flush(lout)

end subroutine trackReport

end module tracking
