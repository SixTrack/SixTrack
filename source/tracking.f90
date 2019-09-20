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
  use beam_beam,   only : bb_trackInit
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
  ! if(isThick)    call bb_trackInit

  ! BEGIN Loop over structure elements
  do 290 i=1,iu
    ix=ic(i)
    if(ix <= nblo) then
      !BLOC
      ktrack(i)=1
      do jb=1,mel(ix)
        jx=mtyp(ix,jb)
        strack(i)=strack(i)+el(jx)
      end do
      if(abs(strack(i)).le.pieni) ktrack(i)=31
      !Non-linear/NOT BLOC
      goto 290
    end if
    ix=ix-nblo
    kpz=abs(kp(ix))
    if(kpz.eq.6) then
      ktrack(i)=2
      goto 290
    end if
    kzz=kz(ix)
    if(kzz.eq.0) then
      ktrack(i)=31
      goto 290
    else if(kzz == 12) then
      ! Disabled cavity; enabled cavities have kp=6 and are handled above
      ! Note: kz=-12 are transformed into +12 in daten after reading ENDE.
      ktrack(i)=31
      goto 290
    end if

    !Beam-beam element
    if(kzz.eq.20) then
      call initialise_element(ix,.false.)
      goto 290
    endif

    ! wire
    if(kzz.eq.15) then
      ktrack(i)=45
      goto 290
    endif
    ! acdip1
    if(kzz.eq.16) then
      ktrack(i)=51
      goto 290
    else if(kzz.eq.-16) then
      ktrack(i)=52
      goto 290
    endif
    ! crab
    if(kzz.eq.23) then
      ktrack(i)=53
      goto 290
    else if(kzz.eq.-23) then
      ktrack(i)=54
      goto 290
    endif
    ! JBG RF CC Multipoles
    if(kzz.eq.26) then
      ktrack(i)=57
      goto 290
    else if(kzz.eq.-26) then
      ktrack(i)=58
      goto 290
    endif
    if(kzz.eq.27) then
      ktrack(i)=59
      goto 290
    else if(kzz.eq.-27) then
      ktrack(i)=60
      goto 290
    endif
    if(kzz.eq.28) then
      ktrack(i)=61
      goto 290
    else if(kzz.eq.-28) then
      ktrack(i)=62
      goto 290
    endif
    !electron lens (HEL)
    if(kzz.eq.29) then
      ktrack(i)=63
      goto 290
    endif
    ! Chebyshev lens
    if(kzz.eq.cheby_kz) then
      ktrack(i)=cheby_ktrack
      goto 290
    endif
    ! SCATTER block
    if (kzz.eq.40 .and. scatter_elemPointer(ix).ne.0) then
      ! FOR NOW, ASSUME THIN SCATTER; ktrack(i)=65 RESERVED FOR THICK SCATTER
      ktrack(i)=64
      goto 290
    endif
    if(kzz.eq.22) then
      ktrack(i)=3
      goto 290
    endif
    if(mout2 == 1 .and. icextal(i) > 0) then
      write(27,"(a16,2x,1p,2d14.6,d17.9)") bez(ix),&
        fluc_errAlign(1,icextal(i)),fluc_errAlign(2,icextal(i)),fluc_errAlign(3,icextal(i))
    end if

    select case (kzz)
    case (1)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 11
#include "include/stra01.f90"
      end if
    case (2)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 12
#include "include/stra02.f90"
      end if
    case (3)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 13
#include "include/stra03.f90"
      end if
    case (4)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 14
#include "include/stra04.f90"
      end if
    case (5)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 15
#include "include/stra05.f90"
      end if
    case (6)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 16
#include "include/stra06.f90"
      end if
    case (7)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 17
#include "include/stra07.f90"
      end if
    case (8)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 18
#include "include/stra08.f90"
      end if
    case (9)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 19
#include "include/stra09.f90"
      end if
    case (10)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 20
#include "include/stra10.f90"
      end if
    case (11) ! Multipole block (also in initialise_element)
      r0  = ek(ix)
      nmz = nmu(ix)
      if(abs(r0).le.pieni.or.nmz.eq.0) then
        if(abs(dki(ix,1)).le.pieni.and.abs(dki(ix,2)).le.pieni) then
          if ( dynk_isused(i) ) then
            write(lerr,"(a)") "TRACKING> ERROR Element of type 11 (bez = '"//trim(bez(ix))//&
              "') is off in "//trim(fort2)//", but on in DYNK. Not implemented."
            call prror
          end if
          ktrack(i) = 31
        else if(abs(dki(ix,1)).gt.pieni.and.abs(dki(ix,2)).le.pieni) then
          if(abs(dki(ix,3)).gt.pieni) then
            ktrack(i) = 33 !Horizontal Bend with a fictive length
#include "include/stra11.f90"
          else
            ktrack(i) = 35 !Horizontal Bend without a ficitve length
#include "include/stra12.f90"
          end if
        else if(abs(dki(ix,1)).le.pieni.and.abs(dki(ix,2)).gt.pieni) then
          if(abs(dki(ix,3)).gt.pieni) then
            ktrack(i) = 37 !Vertical bending with fictive length
#include "include/stra13.f90"
          else
            ktrack(i) = 39 !Vertical bending without fictive length
#include "include/stra14.f90"
          end if
        end if
      else
      !These are the same as above with the difference that they also will have multipoles associated with them.
        if(abs(dki(ix,1)).le.pieni.and.abs(dki(ix,2)).le.pieni) then
          ktrack(i) = 32
        else if(abs(dki(ix,1)).gt.pieni.and.abs(dki(ix,2)).le.pieni) then
          if(abs(dki(ix,3)).gt.pieni) then
            ktrack(i) = 34
#include "include/stra11.f90"
          else
            ktrack(i) = 36
#include "include/stra12.f90"
          end if
        else if(abs(dki(ix,1)).le.pieni.and.abs(dki(ix,2)).gt.pieni) then
          if(abs(dki(ix,3)).gt.pieni) then
            ktrack(i) = 38
#include "include/stra13.f90"
          else
            ktrack(i) = 40
#include "include/stra14.f90"
          end if
        end if
      end if
      if(abs(r0).le.pieni.or.nmz.eq.0) goto 290
      if(mout2.eq.1) then
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
    case (12,13,14,15,16,17,18,19,20,21,22,23)
      goto 290
    case (24) ! DIPEDGE ELEMENT
#include "include/stra2dpe.f90"
      ktrack(i) = 55
    case (25) ! Solenoid
#include "include/solenoid.f90"
      ktrack(i) = 56
    case (41) ! RF Multipole
      ktrack(i) = 66
    case (43) ! xrot
      ktrack(i) = 68
    case (44) ! yrot
      ktrack(i) = 69
    case (45) ! srot
      ktrack(i) = 70

    !----------------
    !--Negative KZZ--
    !----------------
    case (-1)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 21
#include "include/stra01.f90"
      end if
    case (-2)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 22
#include "include/stra02.f90"
      end if
    case (-3)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 23
#include "include/stra03.f90"
      end if
    case (-4)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 24
#include "include/stra04.f90"
      end if
    case (-5)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 25
#include "include/stra05.f90"
      end if
    case (-6)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 26
#include "include/stra06.f90"
      end if
    case (-7)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 27
#include "include/stra07.f90"
      end if
    case (-8)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 28
#include "include/stra08.f90"
      end if
    case (-9)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 29
#include "include/stra09.f90"
      end if
    case (-10)
      if(abs(smiv(i)).le.pieni .and. .not.dynk_isused(i)) then
        ktrack(i) = 31
      else
        ktrack(i) = 30
#include "include/stra10.f90"
      end if

    case default
      ktrack(i) = 31
    end select
290 continue
  ! END Loop over structure elements

end subroutine preTracking

end module tracking
