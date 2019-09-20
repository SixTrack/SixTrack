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
  use mod_common_main
  use mod_common_track
  use numerical_constants

  use collimation, only : do_coll
  use mod_fluc,    only : fluc_writeFort4, fluc_errAlign
  use cheby,       only : cheby_kz, cheby_ktrack
  use dynk,        only : dynk_enabled, dynk_isused, dynk_pretrack
  use scatter,     only : scatter_elemPointer

  real(kind=fPrec) benkcc, r0, r000, r0a
  integer i, j, ix, jb, jx, kpz, kzz, nmz
  logical isThick

  isThick = ithick == 1

  if(isThick .and. do_coll) then
    write(lerr,"(a)") "TRACKING> ERROR Collimation is not supported for thick tracking"
    call prror
  end if

  if(mout2 == 1) call fluc_writeFort4

  ! BEGIN Loop over structure elements
  do i=1,iu

    ix = ic(i) ! Single element index
    if(ix <= nblo) then
      ! This is a block element
      ktrack(i) = 1
      do jb=1,mel(ix)
        jx        = mtyp(ix,jb)
        strack(i) = strack(i)+el(jx)
      end do
      if(abs(strack(i)) <= pieni) then
        ktrack(i) = 31
      end if
      ! Non-linear/NOT BLOC
      cycle
    end if
    ix = ix-nblo ! Remove the block index offset

    if(mout2 == 1 .and. icextal(i) > 0) then
      write(27,"(a16,2x,1p,2d14.6,d17.9)") bez(ix),&
        fluc_errAlign(1,icextal(i)),fluc_errAlign(2,icextal(i)),fluc_errAlign(3,icextal(i))
    end if

    kpz = abs(kp(ix))
    if(kpz == 6) then
      ktrack(i) = 2
      cycle
    end if

    kzz = kz(ix)
    select case(kzz)

    case(0)
      ktrack(i) = 31

    case(12)
      ! Disabled cavity; enabled cavities have kp=6 and are handled above
      ! Note: kz=-12 are transformed into +12 in daten after reading ENDE.
      ktrack(i) = 31

    case(13,14,17,18,19,21)
      cycle

    case(20) ! Beam-beam element
      call initialise_element(ix,.false.)

    case(15) ! Wire
      ktrack(i) = 45

    case(16) ! AC Dipole
      ktrack(i) = 51

    case(-16) ! AC Dipole
      ktrack(i) = 52

    case(22)
      ktrack(i) = 3

    case(23) ! Crab Cavity
      ktrack(i) = 53

    case(-23) ! Crab Cavity
      ktrack(i) = 54

    case(24) ! DIPEDGE ELEMENT
#include "include/stra2dpe.f90"
      ktrack(i) = 55

    case (25) ! Solenoid
#include "include/solenoid.f90"
      ktrack(i) = 56

    case(26) ! JBG RF CC Multipoles
      ktrack(i) = 57

    case(-26) ! JBG RF CC Multipoles
      ktrack(i) = 58

    case(27)
      ktrack(i) = 59

    case(-27)
      ktrack(i) = 60

    case(28)
      ktrack(i) = 61

    case(-28)
      ktrack(i) = 62

    case(29) ! Electron lens (HEL)
      ktrack(i) = 63

    case(cheby_kz) ! Chebyshev lens
      ktrack(i) = cheby_ktrack

    case(40) ! Scatter
      if(scatter_elemPointer(ix) /= 0) then
        ktrack(i) = 64 ! Scatter thin
      else
        ktrack(i) = 31
      end if

    case(41) ! RF Multipole
      ktrack(i) = 66
    case(43) ! X Rotation
      ktrack(i) = 68
    case(44) ! Y Rotation
      ktrack(i) = 69
    case(45) ! S Rotation
      ktrack(i) = 70

    case(-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10)
      if(abs(smiv(i)) <= pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        if(kzz > 0) then
          ktrack(i) = 10 + kzz
        else
          ktrack(i) = 20 - kzz
        end if
        call setStrack(abs(kzz),i)
      end if

    case (11) ! Multipole block (also in initialise_element)
      r0  = ek(ix)
      nmz = nmu(ix)
      if(abs(r0) <= pieni .or. nmz == 0) then
        if(abs(dki(ix,1)) <= pieni .and. abs(dki(ix,2)) <= pieni) then
          if(dynk_isused(i)) then
            write(lerr,"(a)") "TRACKING> ERROR Element of type 11 (bez = '"//trim(bez(ix))//&
              "') is off in "//trim(fort2)//", but on in DYNK. Not implemented."
            call prror
          end if
          ktrack(i) = 31
        else if(abs(dki(ix,1)) > pieni .and. abs(dki(ix,2)) <= pieni) then
          if(abs(dki(ix,3)) > pieni) then
            ktrack(i) = 33 ! Horizontal Bend with a fictive length
#include "include/stra11.f90"
          else
            ktrack(i) = 35 ! Horizontal Bend without a ficitve length
#include "include/stra12.f90"
          end if
        else if(abs(dki(ix,1)) <= pieni.and.abs(dki(ix,2)) > pieni) then
          if(abs(dki(ix,3)) > pieni) then
            ktrack(i) = 37 ! Vertical bending with fictive length
#include "include/stra13.f90"
          else
            ktrack(i) = 39 ! Vertical bending without fictive length
#include "include/stra14.f90"
          end if
        end if
      else
        ! These are the same as above with the difference that they also will have multipoles associated with them.
        if(abs(dki(ix,1)) <= pieni .and. abs(dki(ix,2)) <= pieni) then
          ktrack(i) = 32
        else if(abs(dki(ix,1)) > pieni .and. abs(dki(ix,2)) <= pieni) then
          if(abs(dki(ix,3)) > pieni) then
            ktrack(i) = 34
#include "include/stra11.f90"
          else
            ktrack(i) = 36
#include "include/stra12.f90"
          end if
        else if(abs(dki(ix,1)) <= pieni .and. abs(dki(ix,2)) > pieni) then
          if(abs(dki(ix,3)) > pieni) then
            ktrack(i) = 38
#include "include/stra13.f90"
          else
            ktrack(i) = 40
#include "include/stra14.f90"
          end if
        end if
      end if
      if(abs(r0) <= pieni .or. nmz == 0) then
        cycle
      end if
      if(mout2 == 1) then
        benkcc = ed(ix)*benkc(irm(ix))
        r0a    = one
        r000   = r0*r00(irm(ix))

        do j=1,mmul
          fake(1,j)=(bbiv(j,i)*r0a)/benkcc
          fake(2,j)=(aaiv(j,i)*r0a)/benkcc
          r0a=r0a*r000
        end do

        write(9,'(a16)') bez(ix)
        write(9,'(1p,3d23.15)') (fake(1,j), j=1,3)
        write(9,'(1p,3d23.15)') (fake(1,j), j=4,6)
        write(9,'(1p,3d23.15)') (fake(1,j), j=7,9)
        write(9,'(1p,3d23.15)') (fake(1,j), j=10,12)
        write(9,'(1p,3d23.15)') (fake(1,j), j=13,15)
        write(9,'(1p,3d23.15)') (fake(1,j), j=16,18)
        write(9,'(1p,2d23.15)') (fake(1,j), j=19,20)
        write(9,'(1p,3d23.15)') (fake(2,j), j=1,3)
        write(9,'(1p,3d23.15)') (fake(2,j), j=4,6)
        write(9,'(1p,3d23.15)') (fake(2,j), j=7,9)
        write(9,'(1p,3d23.15)') (fake(2,j), j=10,12)
        write(9,'(1p,3d23.15)') (fake(2,j), j=13,15)
        write(9,'(1p,3d23.15)') (fake(2,j), j=16,18)
        write(9,'(1p,2d23.15)') (fake(2,j), j=19,20)

        do j=1,20
          fake(1,j)=zero
          fake(2,j)=zero
        end do
      end if

    case default
      ktrack(i) = 31

    end select
  end do
  ! END Loop over structure elements

end subroutine preTracking

subroutine setStrack(kzz, i)

  use crcoall
  use mod_common,       only : tiltc, tilts
  use mod_common_main,  only : smiv
  use mod_common_track, only : strack, strackc, stracks
  use numerical_constants

  integer, intent(in) :: kzz
  integer, intent(in) :: i

  select case(kzz)
  case(1)
    strack(i) = smiv(i)*c1e3
  case(2)
    strack(i) = smiv(i)
  case(3)
    strack(i) = smiv(i)*c1m3
  case(4)
    strack(i) = smiv(i)*c1m6
  case(5)
    strack(i) = smiv(i)*c1m9
  case(6)
    strack(i) = smiv(i)*c1m12
  case(7)
    strack(i) = smiv(i)*c1m15
  case(8)
    strack(i) = smiv(i)*c1m18
  case(9)
    strack(i) = smiv(i)*c1m21
  case(10)
    strack(i) = smiv(i)*c1m24
  case default
    write(lerr,"(a,i0,a)") "TRACKING> ERROR Setting strack for type ",kzz," not possible. This is a bug."
    call prror
  end select

#ifdef TILT
  strackc(i) = strack(i)*tiltc(i)
  stracks(i) = strack(i)*tilts(i)
#endif

end subroutine setStrack

end module tracking
