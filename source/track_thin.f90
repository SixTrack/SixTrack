!!--------------------------------------------------------------------------
!!  TRACK THIN LENS PART
!!  F. SCHMIDT
!!  CHANGES FOR COLLIMATION MADE BY G. ROBERT-DEMOLAIZE, October 29th, 2004
!!--------------------------------------------------------------------------
subroutine trauthin(nthinerr)
  ! Updated to Fortran 2015 by V.K.B. Olsen, 19/11/2017
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use scatter, only : scatter_elemPointer
  use dynk, only : dynk_enabled, dynk_isused, dynk_pretrack
  use cheby, only : cheby_kz, cheby_ktrack

  use mod_alloc
  use mod_time
  use mod_units

#ifdef FLUKA
  use mod_fluka
#endif

  use collimation

  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  use mod_fluc, only : fluc_errAlign,fluc_writeFort4
  implicit none
  integer i,ix,j,jb,jj,jx,kpz,kzz,napx0,nmz,nthinerr
  real(kind=fPrec) benkcc,r0,r000,r0a

  real(kind=fPrec), allocatable :: crkveb(:) !(npart)
  real(kind=fPrec), allocatable :: cikveb(:) !(npart)
  real(kind=fPrec), allocatable :: rho2b(:) !(npart)
  real(kind=fPrec), allocatable :: tkb(:) !(npart)
  real(kind=fPrec), allocatable :: r2b(:) !(npart)
  real(kind=fPrec), allocatable :: rb(:) !(npart)
  real(kind=fPrec), allocatable :: rkb(:) !(npart)
  real(kind=fPrec), allocatable :: xrb(:) !(npart)
  real(kind=fPrec), allocatable :: zrb(:) !(npart)
  real(kind=fPrec), allocatable :: xbb(:) !(npart)
  real(kind=fPrec), allocatable :: zbb(:) !(npart)
  real(kind=fPrec), allocatable :: crxb(:) !(npart)
  real(kind=fPrec), allocatable :: crzb(:) !(npart)
  real(kind=fPrec), allocatable :: cbxb(:) !(npart)
  real(kind=fPrec), allocatable :: cbzb(:) !(npart)
  integer :: nbeaux(nbb)
  save

  call alloc(crkveb, npart, zero, "crkveb")
  call alloc(cikveb, npart, zero, "cikveb")
  call alloc(rho2b, npart, zero, "rho2b")
  call alloc(tkb, npart, zero, "tkb")
  call alloc(r2b, npart, zero, "r2b")
  call alloc(rb, npart, zero, "rb")
  call alloc(rkb, npart, zero, "rkb")
  call alloc(xrb, npart, zero, "xrb")
  call alloc(zrb, npart, zero, "zrb")
  call alloc(xbb, npart, zero, "xbb")
  call alloc(zbb, npart, zero, "zbb")
  call alloc(crxb, npart, zero, "crxb")
  call alloc(crzb, npart, zero, "crzb")
  call alloc(cbxb, npart, zero, "cbxb")
  call alloc(cbzb, npart, zero, "cbzb")

  do i=1,nblz
    ktrack(i)=0
    strack(i)=zero
    strackc(i)=zero
    stracks(i)=zero
  end do

  do 290 i=1,iu
    if(mout2.eq.1.and.i.eq.1) call fluc_writeFort4
    ix=ic(i)
    if(ix.gt.nblo) goto 30
    !BLOC
    ktrack(i)=1
    do 20 jb=1,mel(ix)
      jx=mtyp(ix,jb)
      strack(i)=strack(i)+el(jx)
20   continue
    if(abs(strack(i)).le.pieni) ktrack(i)=31
    goto 290
    !Non-linear/NOT BLOC
30   ix=ix-nblo
    kpz=abs(kp(ix))
    if(kpz.eq.6) then
      ktrack(i)=2
      goto 290
    endif
    kzz=kz(ix)
    if(kzz.eq.0) then
      ktrack(i)=31
      goto 290
    else if(kzz.eq.12) then
      !Disabled cavity; enabled cavities have kp=6 and are handled above
      ! Note: kz=-12 are transformed into +12 in daten after reading ENDE.
      ktrack(i)=31
      goto 290
    endif

    !Beam-beam element
    !41 --round beam
    !42 --elliptic beam x>z
    !43--elliptic beam z>x
    !44 -- 6d beam-beam
    if(kzz.eq.20) then
        call initialize_element(ix,.false.)
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
        ktrack(i)=12
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
    case (11) ! Multipole block (also in initialize_element)
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
          fake(1,j)=(bbiv(j,i)*r0a)/benkcc                           !hr01
          fake(2,j)=(aaiv(j,i)*r0a)/benkcc                           !hr01
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
      ktrack(i)=31
    end select

290 continue
  do j=1,napx
    dpsv1(j)=(dpsv(j)*c1e3)/(one+dpsv(j))
  end do

  if (dynk_enabled) call dynk_pretrack
  call time_timeStamp(time_afterPreTrack)

  if ((idp == 0 .or. ition == 0) .and. .not.do_coll) then !4D tracking (not collimat compatible)
    write(lout,"(a)") ""
    write(lout,"(a)") "TRACKING> Calling thin4d subroutine"
    write(lout,"(a)") ""
    call thin4d(nthinerr)
  else !6D tracking
    if(idp == 0 .or. ition == 0) then !Actually 4D, but collimation needs 6D so goto 6D.
      write(lout,"(a)") "TRACKING> WARNING Calling 6D tracking due to collimation! Would normally have called thin4d"
    endif

    hsy(3)=(c1m3*hsy(3))*real(ition,fPrec)
    do jj=1,nele
      if(abs(kz(jj)) == 12) then
        hsyc(jj) = (c1m3*hsyc(jj)) * real(sign(1,kz(jj)),kind=fPrec)
      end if
    end do
    if(abs(phas).ge.pieni) then
      write(lerr,"(a)") "TRACKING> ERROR thin6dua no longer supported. Please use DYNK instead."
      call prror
    else
      write(lout,"(a)") ""
      write(lout,"(a)") "TRACKING> Calling thin6d subroutine"
      write(lout,"(a)") ""
      call thin6d(nthinerr)
    endif !end if(abs(phas).ge.pieni) then
  endif !end if((idp.eq.0.or.ition.eq.0) .and. .not.do_coll) then ... else

  call dealloc(crkveb, "crkveb")
  call dealloc(cikveb, "cikveb")
  call dealloc(rho2b, "rho2b")
  call dealloc(tkb, "tkb")
  call dealloc(r2b, "r2b")
  call dealloc(rb, "rb")
  call dealloc(rkb, "rkb")
  call dealloc(xrb, "xrb")
  call dealloc(zrb, "zrb")
  call dealloc(xbb, "xbb")
  call dealloc(zbb, "zbb")
  call dealloc(crxb, "crxb")
  call dealloc(crzb, "crzb")
  call dealloc(cbxb, "cbxb")
  call dealloc(cbzb, "cbzb")

  return

end subroutine trauthin

!>
!!  TRACK THIN LENS 4D
!!  F. SCHMIDT
!<
subroutine thin4d(nthinerr)
  ! Replaced computed goto with select case. VKBO 27/11/2017

  use floatPrecision
  use string_tools
  use physical_constants
  use numerical_constants
  use mathlib_bouncer
  use mod_particles
  use dynk, only : dynk_enabled, dynk_apply
  use dump, only : dump_linesFirst, dump_lines, ldumpfront
  use collimation, only: do_coll, part_abs_turn
  use aperture

#ifdef FLUKA
  ! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
  ! last modified: 17-07-2013
  ! import mod_fluka
  ! inserted in main code by the 'fluka' compilation flag
  use mod_fluka
#endif

#ifdef ROOT
  use root_output
#endif

  use mod_meta
  use mod_settings
  use postprocessing, only : writebin
  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  use bdex, only : bdex_enable
  use aperture
  use elens
  use cheby, only : cheby_ktrack, cheby_kick
  use utils
  use wire
#ifdef CR
  use checkpoint_restart
#endif
#ifdef BOINC
  use mod_boinc
#endif

  implicit none

  integer i,irrtr,ix,j,k,n,nmz,nthinerr,xory,nac,nfree,nramp1,nplato,nramp2,turnrep,kxxa,nfirst
  real(kind=fPrec) pz,cccc,cikve,crkve,crkveuk,r0,stracki,xlvj,yv1j,yv2j,zlvj,acdipamp,qd,acphase,  &
    acdipamp2,acdipamp1,crabamp,crabfreq,kcrab,RTWO,NNORM,l,cur,dx,dy,tx,ty,embl,chi,xi,yi,dxi,dyi, &
    rrelens,frrelens,xelens,yelens,onedp,fppsig,tan_t,sin_t,cos_t,costh_temp,sinth_temp,pxf,pyf,    &
    r_temp,z_temp,sigf,q_temp,pttemp,xlv,zlv,temp_angle
  logical llost
  real(kind=fPrec) crkveb(npart),cikveb(npart),rho2b(npart),tkb(npart),r2b(npart),rb(npart),        &
    rkb(npart),xrb(npart),zrb(npart),xbb(npart),zbb(npart),crxb(npart),crzb(npart),cbxb(npart),     &
    cbzb(npart)
  real(kind=fPrec) :: krf, x_t, y_t
  complex(kind=fPrec) :: Cp0, Sp1
  complex(kind=fPrec), parameter :: imag=(zero,one)

  save
!-----------------------------------------------------------------------
  nthinerr=0

  ! initialise variables for back-tracking particles
  if (lbacktracking) call aperture_backTrackingInit

#ifdef FLUKA
  napxto = 0
#endif

  ! Determine which turns to print tracking report on
  if(numl > 1000) then
    turnrep = nint(numl/1000.0)
  else
    turnrep = 1
  end if

#ifdef CR
  if(cr_restart) then
    call crstart
    write(crlog,"(2(a,i0))") "TRACKING> Thin 4D restarting on turn ",cr_numl," / ",numl
  end if
  nnuml  = numl
  nfirst = cr_numl
#else
  nfirst = 1
#endif
  do 640 n=nfirst,numl
    if(st_quiet < 3) then
      if(mod(n,turnrep) == 0) then
        call trackReport(n)
      end if
    end if
    meta_nPartTurn = meta_nPartTurn + napx
    numx=n-1

#ifndef FLUKA
    if(mod(numx,nwri) == 0) call writebin(nthinerr)
    if(nthinerr /= 0) return
#endif

#ifdef CR
#ifdef BOINC
    call boinc_turn(n)
#else
    if(mod(numx,numlcp) == 0) call crpoint
#endif
    cr_restart = .false.
    if(st_killswitch) call cr_killSwitch(n)
#endif

    ! A.Mereghetti, for the FLUKA Team
    ! last modified: 03-09-2014
    ! apply dynamic kicks
    ! always in main code
    if ( dynk_enabled ) then
      call dynk_apply(n)
    end if
    call dump_linesFirst(n)

    ! loop over structure elements, single element: name + type + parameter,
    ! structure element = order of single elements/blocks
    do 630 i=1,iu
      ! No if(ktrack(i).eq.1) - a BLOC - is needed in thin tracking,
      ! as no dependency on ix in this case.
      ix=ic(i)-nblo ! ix = index of single element
!Should this be inside "if ktrack .ne. 1"? (time)

      if (ldumpfront) then
        call dump_lines(n,i,ix)
      end if

#ifdef FLUKA
      ! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
      ! last modified: 17-07-2013
      ! is the current entry an instance of a FLUKA element?
      ! inserted in main code by the 'fluka' compilation flag
      if (fluka_enable) then
        if(ktrack(i).ne.1) then ! Skip BLOCs, FLUKA elements must
                                !      be SINGLE ELEMENTs
          if(fluka_type(ix).ne.FLUKA_NONE) then
            if(fluka_type(ix).eq.FLUKA_ELEMENT) then
              call kernel_fluka_element( n, i, ix )
              ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
              ! last modified: 07-03-2018
              ! store old particle coordinates
              if (lbacktracking) call aperture_saveLastCoordinates(i,ix,0)
              goto 620
            else if(fluka_type(ix).eq.FLUKA_ENTRY) then
              fluka_inside = .true.
              call kernel_fluka_entrance( n, i, ix )
              goto 625
            else if(fluka_type(ix).eq.FLUKA_EXIT) then
              fluka_inside = .false.
              call kernel_fluka_exit
              ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
              ! last modified: 07-03-2018
              ! store old particle coordinates
              if (lbacktracking) call aperture_saveLastCoordinates(i,ix,0)
              goto 620
            end if
          end if
        end if
        if(fluka_inside) then
          if(fluka_debug) then
            write(lout,"(a,i0)") "FLUKA> Skipping lattice element at ", i
            write(fluka_log_unit,*) '# Skipping lattice element at ', i
          end if
          goto 630
        end if
      end if
#endif

          if (bdex_enable) then
              write(lerr,"(a)") "BDEX> ERROR BDEX only available for thin6d"
              call prror
          endif

      select case (ktrack(i))
      case (1)
        stracki=strack(i)
        if(iexact) then ! exact drift?
          do j=1,napx
            xv1(j)=xv1(j)*c1m3
            xv2(j)=xv2(j)*c1m3
            yv1(j)=yv1(j)*c1m3
            yv2(j)=yv2(j)*c1m3
            pz=sqrt(one-(yv1(j)**2+yv2(j)**2))
            xv1(j)=xv1(j)+stracki*(yv1(j)/pz)
            xv2(j)=xv2(j)+stracki*(yv2(j)/pz)
            xv1(j)=xv1(j)*c1e3
            xv2(j)=xv2(j)*c1e3
            yv1(j)=yv1(j)*c1e3
            yv2(j)=yv2(j)*c1e3
          enddo
        else
          do j=1,napx
            xv1(j)=xv1(j)+stracki*yv1(j)
            xv2(j)=xv2(j)+stracki*yv2(j)
          end do
        end if
        ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
        ! last modified: 07-03-2018
        ! store old particle coordinates
        if (lbacktracking) call aperture_saveLastCoordinates(i,ix,0)
        goto 630
      case (3)  !Phase Trombone  
        irrtr=imtr(ix)
        do j=1,napx
#include "include/trombone.f90" 
        enddo
      goto 620
      case (2,4,5,6,7,8,9,10)
        goto 630
      case (11) ! HORIZONTAL DIPOLE
        do j=1,napx
#include "include/kickv01h.f90"
        end do
        goto 620
      case (12) ! NORMAL QUADRUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvxxh.f90"
        end do
        goto 620
      case (13) ! NORMAL SEXTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 620
      case (14) ! NORMAL OCTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 620
      case (15) ! NORMAL DECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 620
      case (16) ! NORMAL DODECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 620
      case (17) ! NORMAL 14-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 620
      case (18) ! NORMAL 16-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 620
      case (19) ! NORMAL 18-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 620
      case (20) ! NORMAL 20-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 620
      case (21) ! VERTICAL DIPOLE
        do j=1,napx
#include "include/kickv01v.f90"
        end do
        goto 620
      case (22) ! SKEW QUADRUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvxxv.f90"
        end do
        goto 620
      case (23) ! SKEW SEXTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 620
      case (24) ! SKEW OCTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 620
      case (25) ! SKEW DECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 620
      case (26) ! SKEW DODECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 620
      case (27) ! SKEW 14-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 620
      case (28) ! SKEW 16-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 620
      case (29) ! SKEW 18-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 620
      case (30) ! SKEW 20-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 620
      case (31)
        goto 620
      case (32)
        goto 390
      case (33)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v01.f90"
        end do
        goto 620
      case (34)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v01.f90"
        end do
        goto 390
      case (35)
        do j=1,napx
#include "include/mul4v02.f90"
        end do
        goto 620
      case (36)
        do j=1,napx
#include "include/mul4v02.f90"
        end do
        goto 390
      case (37)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v03.f90"
        end do
        goto 620
      case (38)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v03.f90"
        end do
        goto 390
      case (39)
        do j=1,napx
#include "include/mul4v04.f90"
        end do
        goto 620
      case (40)
        do j=1,napx
#include "include/mul4v04.f90"
        end do
        goto 390
      case (41)
        goto 680
      case (42)
        if(ibtyp.eq.0) then
#include "include/beam11.f90"
#include "include/beama1.f90"
#include "include/beamco.f90"
#include "include/beama2.f90"
#include "include/beam12.f90"
#include "include/beama3.f90"
#include "include/beam13.f90"
#include "include/beama4.f90"
        else if(ibtyp.eq.1) then
#include "include/beam11.f90"
#include "include/beama1.f90"
#include "include/beamco.f90"
#include "include/beama2.f90"
#include "include/beama3.f90"
#include "include/beamwzf1.f90"
#include "include/beama4.f90"
        end if
        goto 620
      case (43)
        if(ibtyp.eq.0) then
#include "include/beam21.f90"
#include "include/beama1.f90"
#include "include/beamco.f90"
#include "include/beama2.f90"
#include "include/beam22.f90"
#include "include/beama3.f90"
#include "include/beam23.f90"
#include "include/beama4.f90"
        else if(ibtyp.eq.1) then
#include "include/beam21.f90"
#include "include/beama1.f90"
#include "include/beamco.f90"
#include "include/beama2.f90"
#include "include/beama3.f90"
#include "include/beamwzf2.f90"
#include "include/beama4.f90"
        end if
        goto 620
      case (44,46,47,48,49,50,57,58,59,60,61,62)
        goto 630
      case (45) ! Wire
#include "include/wirekick.f90"
        goto 620
      case (51)
#include "include/acdipkick1.f90"
        goto 620
      case (52)
#include "include/acdipkick2.f90"
        goto 620
      case (53)
#include "include/crabkick1.f90"
        goto 620
      case (54)
#include "include/crabkick2.f90"
        goto 620
      case (55) ! DIPEDGE ELEMENT
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvdpe.f90"
        end do
        goto 620
      case (56) ! Solenoid
        do j=1,napx
#include "include/kickvso1.f90"
        end do
        goto 620
      case (63) ! Elens
        do j=1,napx
#include "include/kickelens.f90"
        end do
        goto 620
       case (66) ! Rf-multi
#include "include/rfmulti.f90"
        goto 620
      case (cheby_ktrack) ! Chebyshev lens
        call cheby_kick(i,ix,n)
        goto 620
      case (68) ! xrot
        temp_angle = ed(ix)
#include "include/xrot.f90"
        goto 620
      case (69) ! yrot
        temp_angle = ed(ix)
#include "include/yrot.f90"
        goto 620
      case (70) ! srot
        temp_angle = ed(ix)
#include "include/srot.f90"
        goto 620

      end select
      goto 630

390   r0=ek(ix)
      nmz=nmu(ix)
      if(nmz.ge.2) then
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v05.f90"
          do k=3,nmz
#include "include/mul4v06.f90"
          end do
#include "include/mul4v07.f90"
        end do
      else
        do j=1,napx
#include "include/mul4v08.f90"
        end do
      end if
      goto 620


680   continue
      do 690 j=1,napx
#include "include/beamco.f90"
#include "include/beamr1.f90"
      &goto 690
#include "include/beamr2.f90"
#include "include/beamr3.f90"
690   continue
      goto 620

!----------------------------
620 continue

#include "include/lostpart.f90"

625 continue
    if (.not. ldumpfront) then
      call dump_lines(n,i,ix)
    end if

630 continue

#if defined(ROOT)
    if(root_flag .and. root_Collimation.eq.1) then
      call SurvivalRootWrite(n, napx)
    end if
#endif

    if(nthinerr.ne.0) return
    if(ntwin.ne.2) call dist1
#ifndef FLUKA
    if(mod(n,nwr(4)).eq.0) call write6(n)
#endif

#ifdef FLUKA
  ! A.Mereghetti, for the FLUKA Team
  ! last modified: 14-06-2014
  ! increase napxto, to get an estimation of particles*turns
  ! inserted in main code by the 'fluka' compilation flag
  napxto = napxto + napx
#endif

  640 continue

  return

end subroutine thin4d

!>
!!  TRACK THIN LENS 6D
!!  F. SCHMIDT
!<
subroutine thin6d(nthinerr)
  ! Replaced computed gotos with select case, VKBO 27/11/2017
  use floatPrecision
  use string_tools
  use physical_constants
  use numerical_constants
  use mathlib_bouncer
  use mod_particles

  use bdex,       only : bdex_track, bdex_enable, bdex_elementAction
  use scatter,    only : scatter_thin, scatter_debug
  use dynk,       only : dynk_enabled, dynk_apply
  use dump,       only : dump_linesFirst, dump_lines, ldumpfront
  use mod_ffield, only : ffindex,ffield_genAntiQuad,ffield_enterQuad,ffield_exitQuad,ffield_enabled
  use aperture
  use mod_settings
  use mod_meta
  use mod_time

#ifdef FLUKA
  use mod_fluka
#endif
#ifdef ROOT
  use root_output
#endif

  use collimation
  use postprocessing, only : writebin
  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  use aperture
  use elens
  use cheby, only : cheby_ktrack, cheby_kick
  use utils
  use wire
#ifdef CR
  use checkpoint_restart
#endif
#ifdef BOINC
  use mod_boinc
#endif

  implicit none

  integer i,irrtr,ix,j,k,n,nmz,nthinerr,dotrack,xory,nac,nfree,nramp1,nplato,nramp2,turnrep,elemEnd,&
    kxxa,nfirst
  real(kind=fPrec) pz,cccc,cikve,crkve,crkveuk,r0,stracki,xlvj,yv1j,yv2j,zlvj,acdipamp,qd,          &
    acphase,acdipamp2,acdipamp1,crabamp,crabfreq,crabamp2,crabamp3,crabamp4,kcrab,RTWO,NNORM,l,cur, &
    dx,dy,tx,ty,embl,chi,xi,yi,dxi,dyi,rrelens,frrelens,xelens,yelens, onedp,fppsig,costh_temp,     &
    sinth_temp,tan_t,sin_t,cos_t,pxf,pyf,r_temp,z_temp,sigf,q_temp,pttemp,xlv,zlv,temp_angle
  logical llost, doFField
  real(kind=fPrec) crkveb(npart),cikveb(npart),rho2b(npart),tkb(npart),r2b(npart),rb(npart),        &
    rkb(npart),xrb(npart),zrb(npart),xbb(npart),zbb(npart),crxb(npart),crzb(npart),cbxb(npart),     &
    cbzb(npart)
  real(kind=fPrec) :: krf, x_t, y_t
  complex(kind=fPrec) :: Cp0, Sp1
  complex(kind=fPrec), parameter :: imag=(zero,one)
  save

  nthinerr=0

  ! initialise variables for back-tracking particles
  if (lbacktracking) call aperture_backTrackingInit

#ifdef FLUKA
  napxto = 0
#endif

  ! Determine which turns to print tracking report on
  if(numl > 1000) then
    turnrep = nint(numl/1000.0)
  else
    turnrep = 1
  end if

  call ffield_genAntiQuad()

  ! This is the loop over turns: label 660
#ifdef CR
  if(cr_restart) then
    call crstart
    write(crlog,"(2(a,i0))") "TRACKING> Thin 6D restarting on turn ",cr_numl," / ",numl
  end if
  nnuml  = numl
  nfirst = cr_numl
#else
  nfirst = 1
#endif
  do 660 n=nfirst,numl
    if(st_quiet < 3) then
      if(mod(n,turnrep) == 0) then
        call trackReport(n)
      end if
    end if
    meta_nPartTurn = meta_nPartTurn + napx

    if (do_coll) then
      ! This subroutine sets variables iturn and totals
      call collimate_start_turn(n)
    endif

    numx=n-1

#ifndef FLUKA
    if(mod(numx,nwri) == 0) call writebin(nthinerr)
    if(nthinerr /= 0) return
#endif

#ifdef CR
#ifdef BOINC
    call boinc_turn(n)
#else
    if(mod(numx,numlcp) == 0) call crpoint
#endif
    cr_restart = .false.
    if(st_killswitch) call cr_killSwitch(n)
#endif

    ! A.Mereghetti, for the FLUKA Team
    ! last modified: 03-09-2014
    ! apply dynamic kicks
    ! always in main code
    if ( dynk_enabled ) then
      call dynk_apply(n)
    end if

    call dump_linesFirst(n)

    !! This is the loop over each element: label 650
    do 650 i=1,iu !Loop over elements

      if (do_coll) then
        ! This subroutine sets variables myktrack and myix
        call collimate_start_element(i)
      endif

      ! No if(ktrack(i).eq.1) - a BLOC - is needed in thin tracking,
      ! as no dependency on ix in this case.
      ix=ic(i)-nblo

      ! Fringe Fields
      if(ffield_enabled .and. ix > 0) then
        doFField = FFindex(ix) > 0
      else
        doFField = .false.
      end if

#ifdef BEAMGAS
      !YIL Call beamGas subroutine whenever a pressure-element is found
      ! should be faster/safer to first check the turn then do the name search
      if( iturn.eq.1 ) then
        if (bez(myix)(1:5).eq.'PRESS' .or.  bez(myix)(1:5).eq.'press' ) then
          call beamGas(myix, secondary,totals,myenom,ipart)
        end if
      end if
#endif

      if (ldumpfront) then
        call dump_lines(n,i,ix)
      end if

#ifdef FLUKA
      ! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
      ! last modified: 17-07-2013
      ! is the current entry an instance of a FLUKA element?
      ! inserted in main code by the 'fluka' compilation flag
      if (fluka_enable) then
        if(ktrack(i).ne.1) then ! Skip BLOCs, FLUKA elements must
                                !      be SINGLE ELEMENTs
          if(fluka_type(ix).ne.FLUKA_NONE) then
            if(fluka_type(ix).eq.FLUKA_ELEMENT) then
              call kernel_fluka_element( n, i, ix )
              ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
              ! last modified: 07-03-2018
              ! store old particle coordinates
              if (lbacktracking) call aperture_saveLastCoordinates(i,ix,0)
              goto 640
            else if(fluka_type(ix).eq.FLUKA_ENTRY) then
              fluka_inside = .true.
              call kernel_fluka_entrance( n, i, ix )
              goto 645
            else if(fluka_type(ix).eq.FLUKA_EXIT) then
              fluka_inside = .false.
              call kernel_fluka_exit
              ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
              ! last modified: 07-03-2018
              ! store old particle coordinates
              if (lbacktracking) call aperture_saveLastCoordinates(i,ix,0)
              goto 640
            end if
          end if
        end if
        if(fluka_inside) then
          if(fluka_debug) then
            write(lout,"(a,i0)") "FLUKA> Skipping lattice element at ",i
            write(fluka_log_unit,*) '# Skipping lattice element at ', i
          end if
          goto 650
        end if
      end if
#endif

      ! BDEX was in a #ifndef collimat block, and may not be fully collimat-compatible,
      ! so for now the two shall not be mixed.
      ! Also, check that we have a single element, not a BLOC element.
      if(.not.do_coll .and. ix > 0) then
        if(bdex_enable .and. kz(ix) == 0 .and. bdex_elementAction(ix) /= 0) call bdex_track(i,ix,n)
      end if

      if (do_coll) then
        dotrack = myktrack
      else
        dotrack = ktrack(i)
      end if

      select case(dotrack)
      case (1)
        stracki=strack(i)

        if (do_coll) then
          !==========================================
          !Ralph drift length is stracki
          !bez(ix) is name of drift
          totals=totals+stracki

          !________________________________________________________________________
          !++  If we have a collimator then...
          !
          !Feb2006
          !GRD (June 2005) 'COL' option is for RHIC collimators
          !
          !     SR (17-01-2006): Special assignment to the TCS.TCDQ for B1 and B4,
          !     using the new naming as in V6.500.
          !     Note that this must be in the loop "if TCSG"!!
          !
          !     SR, 17-01-2006: Review the TCT assignments because the MADX names
          !     have changes (TCTH.L -> TCTH.4L)
          !
          ! JULY 2008 added changes (V6.503) for names in TCTV -> TCTVA and TCTVB
          ! both namings before and after V6.503 can be used
          !
          elemEnd = len_trim(bez(myix))
          if((    bez(myix)(1:2) == 'TC'  .or. bez(myix)(1:2) == 'tc'   &
            .or.  bez(myix)(1:2) == 'TD'  .or. bez(myix)(1:2) == 'td'   &
            .or.  bez(myix)(1:3) == 'COL' .or. bez(myix)(1:3) == 'col') &
            .and. bez(myix)(elemEnd-2:elemEnd) /= "_AP") then

            call time_startClock(time_clockCOLL)
            call collimate_start_collimator(stracki)

            !++ For known collimators
            if(found) then
              call collimate_do_collimator(stracki)
              call collimate_end_collimator(stracki)
            end if ! end of check for 'found'
            call time_stopClock(time_clockCOLL)
            !------------------------------------------------------------------
            !++  Here leave the known collimator IF loop...
            !_______________________________________________________________________
            !++  If it is just a drift (i.e. did not pass the namecheck)
          else
            ! TODO: Could just as well call normal sixtrack code (below)...
            do j=1,napx
              xv1(j)  = xv1(j) + stracki*yv1(j)
              xv2(j)  = xv2(j) + stracki*yv2(j)
              sigmv(j) = sigmv(j) + stracki*(c1e3-rvv(j)*(c1e3+(yv1(j)*yv1(j)+yv2(j)*yv2(j))*c5m4))
              xj     = (xv1(j)-torbx(ie))/c1e3
              xpj    = (yv1(j)-torbxp(ie))/c1e3
              yj     = (xv2(j)-torby(ie))/c1e3
              ypj    = (yv2(j)-torbyp(ie))/c1e3
              pj     = ejv(j)/c1e3

              if(firstrun) then
                if (iturn.eq.1.and.j.eq.1) then
                  sum_ax(ie)=zero
                  sum_ay(ie)=zero
                end if
              end if

              gammax = (one + talphax(ie)**2)/tbetax(ie)
              gammay = (one + talphay(ie)**2)/tbetay(ie)

              if (part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
                nspx    = sqrt(abs(gammax*(xj)**2 + two*talphax(ie)*xj*xpj + tbetax(ie)*xpj**2)/myemitx0_collgap)
                nspy    = sqrt(abs(gammay*(yj)**2 + two*talphay(ie)*yj*ypj + tbetay(ie)*ypj**2)/myemity0_collgap)
                sum_ax(ie)   = sum_ax(ie) + nspx
                sqsum_ax(ie) = sqsum_ax(ie) + nspx**2
                sum_ay(ie)   = sum_ay(ie) + nspy
                sqsum_ay(ie) = sqsum_ay(ie) + nspy**2
                nampl(ie)    = nampl(ie) + 1
              else
                nspx = zero
                nspy = zero
              end if
              sampl(ie)    = totals
              ename(ie)    = bez(myix)(1:mNameLen)
            end do
          endif
          !GRD END OF THE CHANGES FOR COLLIMATION STUDIES, BACK TO NORMAL SIXTRACK STUFF

        else ! Normal SixTrack drifts
          if(iexact) then
            ! EXACT DRIFT
            do j=1,napx
              xv1(j)=xv1(j)*c1m3
              xv2(j)=xv2(j)*c1m3
              yv1(j)=yv1(j)*c1m3
              yv2(j)=yv2(j)*c1m3
              sigmv(j)=sigmv(j)*c1m3
              pz=sqrt(one-(yv1(j)**2+yv2(j)**2))
              xv1(j)=xv1(j)+stracki*(yv1(j)/pz)
              xv2(j)=xv2(j)+stracki*(yv2(j)/pz)
              sigmv(j)=sigmv(j)+stracki*(one-(rvv(j)/pz))
              xv1(j)=xv1(j)*c1e3
              xv2(j)=xv2(j)*c1e3
              yv1(j)=yv1(j)*c1e3
              yv2(j)=yv2(j)*c1e3
              sigmv(j)=sigmv(j)*c1e3
            enddo
          else
            do j=1,napx
              xv1(j)  = xv1(j) + stracki*yv1(j)
              xv2(j)  = xv2(j) + stracki*yv2(j)
              sigmv(j) = sigmv(j) + stracki*(c1e3-rvv(j)*(c1e3+(yv1(j)**2+yv2(j)**2)*c5m4))
            end do
          end if
        end if
        ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
        ! last modified: 07-03-2018
        ! store old particle coordinates
        ! NB: end up here in case of collimators too, but not in
        !     in case of an aperture marker or other null-length non-active
        !     elements, thanks to trauthin/trauthck and
        !     if(abs(strack(i)).le.pieni) ktrack(i)=31
        if (lbacktracking) call aperture_saveLastCoordinates(i,ix,0)
        goto 650

      case (2)
        do j=1,napx
          ejf0v(j)=ejfv(j)
          if(abs(dppoff).gt.pieni) then
            sigmv(j)=sigmv(j)-sigmoff(i)
          endif
          if(abs(kz(ix)) == 12) then
            ejv(j)=ejv(j)+(ed(ix)*sin_mb(hsyc(ix)*sigmv(j)+phasc(ix)))*nqq(j)
          else
            ejv(j)=ejv(j)+(hsy(1)*sin_mb(hsy(3)*sigmv(j)))*nqq(j)
          endif
          ejfv(j)=sqrt(ejv(j)**2-nucm(j)**2)                               !hr01
          rvv(j)=(ejv(j)*e0f)/(e0*ejfv(j))
          dpsv(j) = (ejfv(j)*(nucm0/nucm(j))-e0f)/e0f
!          dpsv(j)=(ejfv(j)-e0f)/e0f
          oidpsv(j)=one/(one+dpsv(j))
          moidpsv(j)=mtc(j)/(one+dpsv(j))
          omoidpsv(j)=c1e3*((one-mtc(j))*oidpsv(j))
          dpsv1(j)=(dpsv(j)*c1e3)*oidpsv(j)                           !hr01
          yv1(j)=(ejf0v(j)/ejfv(j))*yv1(j)                           !hr01
          yv2(j)=(ejf0v(j)/ejfv(j))*yv2(j)                           !hr01
        end do
        goto 640
      case (3)
        irrtr=imtr(ix)
        do j=1,napx
            !The values are stored in the temp vector which are used for the multiplication.
          temptr(1)=xv1(j)
          temptr(2)=yv1(j)/moidpsv(j)
          temptr(3)=xv2(j)
          temptr(4)=yv2(j)/moidpsv(j)
          temptr(5)=sigmv(j)
          temptr(6)=((mtc(j)*ejv(j)-e0)/e0f)*c1e3*(e0/e0f)
          ! Adding the closed orbit. The previous values are stored in the temptr vector.
          xv1(j)  = cotr(irrtr,1)
          yv1(j)  = cotr(irrtr,2)
          xv2(j)  = cotr(irrtr,3)
          yv2(j)  = cotr(irrtr,4)
          sigmv(j) = cotr(irrtr,5)
          pttemp   = cotr(irrtr,6)

          ! Multiplying the arbitrary matrix to the coordinates.
          do kxxa=1,6
            xv1(j)   =  xv1(j)+temptr(kxxa)*rrtr(irrtr,1,kxxa)
            yv1(j)   =  yv1(j)+temptr(kxxa)*rrtr(irrtr,2,kxxa)
            xv2(j)   =  xv2(j)+temptr(kxxa)*rrtr(irrtr,3,kxxa)
            yv2(j)   =  yv2(j)+temptr(kxxa)*rrtr(irrtr,4,kxxa)
            sigmv(j)  =  sigmv(j)+temptr(kxxa)*rrtr(irrtr,5,kxxa)
            pttemp    =  pttemp+temptr(kxxa)*rrtr(irrtr,6,kxxa)
          enddo
          ! Transforming back to the tracked coordinates of Sixtrack...
          ejv(j)  = (e0f*pttemp/(c1e3*(e0/e0f))+e0)/mtc(j)


          ejfv(j)=sqrt(ejv(j)**2-nucm(j)**2)
          rvv(j)=(ejv(j)*e0f)/(e0*ejfv(j))
          dpsv(j)=(ejfv(j)*(nucm0/nucm(j))-e0f)/e0f
          oidpsv(j)=one/(one+dpsv(j))
          moidpsv(j)=mtc(j)/(one+dpsv(j))
          omoidpsv(j)=c1e3*((one-mtc(j))*oidpsv(j))
          dpsv1(j)=(dpsv(j)*c1e3)*oidpsv(j)


          ! We have to go back to angles after we updated the energy.
          yv1(j) = yv1(j)*mtc(j)/(one+dpsv(j))
          yv2(j) = yv2(j)*mtc(j)/(one+dpsv(j))

          !yv(j,1) = yv(j,1)*moidpsv(j)
          !yv(j,2) = yv(j,2)*moidpsv(j)
        enddo
        goto 640
      case (4,5,6,7,8,9,10)
        goto 650
      case (11) ! HORIZONTAL DIPOLE
        do j=1,napx
#include "include/kickv01h.f90"
        end do
        goto 640
      case (12) ! NORMAL QUADRUPOLE
        if(doFField) then
          if(ic(i) /= ic(i-2) .and. ic(i) /= ic(i-3)) then
            call ffield_enterQuad(i)  !A optimizer!!!
          end if
        end if
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvxxh.f90"
        end do
        if(doFField) then
          if(ic(i) /= ic(i+2) .and. ic(i) /= ic(i+3)) then
            call ffield_exitQuad(i)   !A optimizer!!!
          end if
        end if
        goto 640
      case (13) ! NORMAL SEXTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 640
      case (14) ! NORMAL OCTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 640
      case (15) ! NORMAL DECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 640
      case (16) ! NORMAL DODECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 640
      case (17) ! NORMAL 14-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 640
      case (18) ! NORMAL 16-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 640
      case (19) ! NORMAL 18-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 640
      case (20) ! NORMAL 20-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 640
      case (21) ! VERTICAL DIPOLE
        do j=1,napx
#include "include/kickv01v.f90"
        end do
        goto 640
      case (22) ! SKEW QUADRUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvxxv.f90"
        end do
        goto 640
      case (23) ! SKEW SEXTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 640
      case (24) ! SKEW OCTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 640
      case (25) ! SKEW DECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 640
      case (26) ! SKEW DODECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 640
      case (27) ! SKEW 14-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 640
      case (28) ! SKEW 16-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 640
      case (29) ! SKEW 18-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 640
      case (30) ! SKEW 20-POLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 640
      case (31)
        goto 640
      case (32)
        if(doFField .eqv. .false.) then
          goto 410
        else
          goto 640
        end if
      case (33)
        if(doFField .eqv. .false.) then
          do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v01.f90"
#include "include/mul6v01.f90"
          end do
        end if
        goto 640
      case (34)
        if(doFField .eqv. .false.) then
          do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v01.f90"
#include "include/mul6v01.f90"
          end do
          goto 410
        else
          goto 640
        end if
      case (35)
        if(doFField .eqv. .false.) then
          do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v02.f90"
#include "include/mul6v01.f90"
          end do
        end if
        goto 640
      case (36)
        if(doFField .eqv. .false.) then
          do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v02.f90"
#include "include/mul6v01.f90"
          end do
          goto 410
        else
          goto 640
        end if
      case (37)
        if(doFField .eqv. .false.) then
          do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v03.f90"
#include "include/mul6v02.f90"
          end do
        end if
        goto 640
      case (38)
        if(doFField .eqv. .false.) then
          do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v03.f90"
#include "include/mul6v02.f90"
          end do
          goto 410
        else
          goto 640
        end if
      case (39)
        if(doFField .eqv. .false.) then
          do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v04.f90"
#include "include/mul6v02.f90"
          end do
        end if
        goto 640
      case (40)
        if(doFField .eqv. .false.) then
          do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v04.f90"
#include "include/mul6v02.f90"
          end do
          goto 410
        else
          goto 640
        end if
      case (41) ! 4D BB kick
        do 690 j=1,napx
#include "include/beamco.f90"
#include "include/beamr1.f90"
        &goto 690 !The radius was too small -> Skip
#include "include/beamr2.f90"
#include "include/beamr3.f90"
690     continue
        goto 640
      case (42)
        if(ibtyp.eq.0) then
#include "include/beam11.f90"
#include "include/beama1.f90"
#include "include/beamco.f90"
#include "include/beama2.f90"
#include "include/beam12.f90"
#include "include/beama3.f90"
#include "include/beam13.f90"
#include "include/beama4.f90"
        else if(ibtyp.eq.1) then ! fast kick
#include "include/beam11.f90"
#include "include/beama1.f90"
#include "include/beamco.f90"
#include "include/beama2.f90"
#include "include/beama3.f90"
#include "include/beamwzf1.f90"
#include "include/beama4.f90"
        end if
        goto 640
      case (43)
        if(ibtyp.eq.0) then
#include "include/beam21.f90"
#include "include/beama1.f90"
#include "include/beamco.f90"
#include "include/beama2.f90"
#include "include/beam22.f90"
#include "include/beama3.f90"
#include "include/beam23.f90"
#include "include/beama4.f90"
        else if(ibtyp.eq.1) then
#include "include/beam21.f90"
#include "include/beama1.f90"
#include "include/beamco.f90"
#include "include/beama2.f90"
#include "include/beama3.f90"
#include "include/beamwzf2.f90"
#include "include/beama4.f90"
        end if
        goto 640
      case (44)
#include "include/beam6d.f90"
        goto 640
      case (45) ! Wire
#include "include/wirekick.f90"
        goto 640
      case (46,47,48,49,50)
        goto 650
      case (51)
#include "include/acdipkick1.f90"
        goto 640
      case (52)
#include "include/acdipkick2.f90"
        goto 640
      case (53)
#include "include/crabkick1.f90"
        goto 640
      case (54)
#include "include/crabkick2.f90"
        goto 640
      case (55) ! DIPEDGE ELEMENT
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvdpe.f90"
        end do
        goto 640
      case (56) ! Solenoid
        do j=1,napx
#include "include/kickvso1.f90"
#include "include/kickvso2.f90"
        end do
        goto 640
      case (57) ! JBG RF CC Multipoles
        xory=1
        crabfreq=ek(ix)*c1e3
        do j=1,napx
          crabamp2 = ed(ix)*nqq(j)
          kcrab=(((sigmv(j)/(clight*(e0f/e0)))*crabfreq)*two)*pi + crabph2(ix)
#include "include/alignva.f90"
          yv1(j)=yv1(j) + ((crabamp2*crkve)*moidpsv(j))*cos_mb(kcrab)
          yv2(j)=yv2(j) - ((crabamp2*cikve)*moidpsv(j))*cos_mb(kcrab)
          ejv(j)=ejv(j) - ((((half*(crabamp2))*(crkve**2-cikve**2))*(((crabfreq*two)*pi)/clight))*c1m3)*(sin_mb(kcrab)*e0f)
        end do
        call part_updatePartEnergy(1,.true.)
        goto 640
      case (58) ! JBG RF CC Multipoles
        xory=1
        crabfreq=ek(ix)*c1e3
        do j=1,napx
          crabamp2 = ed(ix)*nqq(j)
          kcrab=(((sigmv(j)/(clight*(e0f/e0)))*crabfreq)*two)*pi + crabph2(ix)
#include "include/alignva.f90"
          yv2(j)=yv2(j) + ((crabamp2*crkve)*moidpsv(j))*cos_mb(kcrab)
          yv1(j)=yv1(j) + ((crabamp2*cikve)*moidpsv(j))*cos_mb(kcrab)
          ejv(j)=ejv(j) - ((((crabamp2)*(cikve*crkve))*(((crabfreq*two)*pi)/clight))*c1m3)*(sin_mb(kcrab)*e0f)
        end do
        call part_updatePartEnergy(1,.true.)
        goto 640
      case (59) ! JBG RF CC Multipoles
        xory=1
        crabfreq=ek(ix)*c1e3
        do j=1,napx
          crabamp3 = ed(ix)*nqq(j)
          kcrab=((sigmv(j)*crabfreq)/(clight*(e0f/e0)))*(two*pi)+crabph3(ix)
#include "include/alignva.f90"
          yv1(j)=yv1(j)+(((crabamp3*moidpsv(j))*c1m3)*(crkve**2-cikve**2))*cos_mb(kcrab)
          yv2(j)=yv2(j)-((two*(((crabamp3*crkve)*cikve)*moidpsv(j)))*c1m3)*cos_mb(kcrab)
          ejv(j)=ejv(j)-(((((one/three)*(crabamp3))*(crkve**3-(three*cikve**2)*crkve))&
                *(((crabfreq*two)*pi)/clight)*c1m6)*sin_mb(kcrab))*e0f
        end do
        call part_updatePartEnergy(1,.true.)
        goto 640
      case (60) ! JBG RF CC Multipoles
        xory=1
        crabfreq=ek(ix)*c1e3
        do j=1,napx
          crabamp3 = ed(ix)*nqq(j)
          kcrab=(((sigmv(j)/(clight*(e0f/e0)))*crabfreq)*two)*pi + crabph3(ix)
#include "include/alignva.f90"
          yv2(j)=yv2(j)-(((crabamp3*moidpsv(j))*c1m3)*((cikve**2)-(crkve**2)))*cos_mb(kcrab)
          yv1(j)=yv1(j)+((two*(crabamp3*(crkve*(cikve*oidpsv(j)))))*c1m3)*cos_mb(kcrab)
          ejv(j)=ejv(j)+(((((one/three)*(crabamp3))*(cikve**3- &
                ((three*crkve**2)*cikve)))*(((crabfreq*two)*pi)/clight))*c1m6)*(sin_mb(kcrab)*e0f)
        end do
        call part_updatePartEnergy(1,.true.)
        goto 640
      case (61) ! JBG RF CC Multipoles
        xory=1
        crabfreq=ek(ix)*c1e3
        do j=1,napx
          crabamp4 = ed(ix)*nqq(j)
          kcrab=(((sigmv(j)/(clight*(e0f/e0)))*crabfreq)*two)*pi + crabph4(ix)
#include "include/alignva.f90"
          yv1(j)=yv1(j) + (((crabamp4*moidpsv(j))*(crkve**3-(three*crkve)*cikve**2))*c1m6)*cos_mb(kcrab)
          yv2(j)=yv2(j) - (((crabamp4*moidpsv(j))*((three*cikve)*crkve**2-cikve**3))*c1m6)*cos_mb(kcrab)
          ejv(j)=ejv(j) - ((((0.25_fPrec*(crabamp4))*(crkve**4-(six*crkve**2)*cikve**2+cikve**4))&
                *(((crabfreq*two)*pi)/clight))*c1m9)*(sin_mb(kcrab)*e0f)
        end do
        call part_updatePartEnergy(1,.true.)
        goto 640
      case (62) ! JBG RF CC Multipoles
        xory=1
        crabfreq=ek(ix)*c1e3
        do j=1,napx
          crabamp4 = ed(ix)*nqq(j)
          kcrab=(((sigmv(j)/(clight*(e0f/e0)))*crabfreq)*two)*pi + crabph4(ix)
#include "include/alignva.f90"
          yv1(j)=yv1(j) - (((crabamp4*moidpsv(j))*(cikve**3-(three*cikve)*crkve**2))*c1m6)*cos_mb(kcrab)
          yv2(j)=yv2(j) - (((crabamp4*moidpsv(j))*((three*crkve)*cikve**2-crkve**3))*c1m6)*cos_mb(kcrab)
          ejv(j)=ejv(j) - ((((crabamp4)*((crkve**3*cikve)-(cikve**3*crkve)))*(((crabfreq*two)*pi)/clight))*c1m9)*(sin_mb(kcrab)*e0f)
        end do
        call part_updatePartEnergy(1,.true.)
        goto 640
      case (63) ! Elens
        do j=1,napx
#include "include/kickelens.f90"
        end do
        goto 640
      case (66) ! Rf-multi
#include "include/rfmulti.f90"
        goto 640

      case (64) ! Scatter (thin)
        !Thin scattering
        ! It is already checked that scatter_elemPointer != 0
        call scatter_thin(i, ix,n)
        goto 640
      case (65) ! Scatter (thick)
        !     TODO
        goto 640
      case (cheby_ktrack) ! Chebyshev lens
        call cheby_kick(i,ix,n)
        goto 640
      case (68) ! xrot
        temp_angle = ed(ix)
#include "include/xrot.f90"
        goto 640
      case (69) ! yrot
        temp_angle = ed(ix)
#include "include/yrot.f90"
        goto 640
      case (70) ! srot
        temp_angle = ed(ix)
#include "include/srot.f90"
        goto 640
      case default
        write(lout,"(3(a,i0),a)") "TRACKING> WARNING Non-handled element in thin6d()!",  &
          " i = ",i,", ix = ",ix,", dotrack = ",dotrack,", bez(ix) = '"//trim(bez(ix))//"' skipped."
      end select
      goto 650

410   r0=ek(ix)
      nmz=nmu(ix)
      if(nmz.ge.2) then
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v05.f90"
          do k=3,nmz
#include "include/mul4v06.f90"
          end do
#include "include/mul4v07.f90"
        end do
      else
        do j=1,napx
#include "include/mul4v08.f90"
        end do
      end if

640   continue ! end of the SELECT CASE over element type (dotrack)

      if (do_coll) then
        call collimate_end_element
      end if

#include "include/lostpart.f90"

645   continue

      if (.not. ldumpfront) then
        call dump_lines(n,i,ix)
      end if

650 continue !END loop over structure elements

    if (do_coll) then
      call collimate_end_turn
    endif

#if defined(ROOT)
    if(root_flag .and. root_Collimation.eq.1) then
      call SurvivalRootWrite(n, napx)
    end if
#endif

    if(nthinerr.ne.0) then
      return
    end if

    if (.not. do_coll) then
      if(ntwin.ne.2) then
        call dist1
      end if
#ifndef FLUKA
      if(mod(n,nwr(4)).eq.0) then
        call write6(n) ! Write to fort.12
      end if
#endif
    endif

#ifdef FLUKA
!   A.Mereghetti, for the FLUKA Team
!   last modified: 14-06-2014
!   increase napxto, to get an estimation of particles*turns
!   inserted in main code by the 'fluka' compilation flag
    napxto = napxto + napx
#endif

    if (do_coll) then
      !GRD HERE WE SET THE FLAG FOR INITIALIZATION TO FALSE AFTER TURN 1
      firstrun = .false.
    endif

660 continue !END loop over turns

end subroutine thin6d

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Write a turn report.
!  The isFirst if statement is only computed the first time the routine is called.
! ================================================================================================ !
subroutine trackReport(n)

  use crcoall
  use parpro,     only : npart
  use mod_common, only : ithick, iclo6, numl, napx, napxo

  implicit none

  integer, intent(in) :: n

  character(len=8)  :: trackMode = " "
  character(len=32) :: trackFmt  = " "
  integer           :: oPart     = 0
  integer           :: oTurn     = 0
  logical           :: isFirst   = .true.

  if(isFirst) then
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
    oPart   = int(log10(real(npart)))+1
    oTurn   = int(log10(real(numl)))+1
    isFirst = .false.
    write(trackFmt,"(2(a,i0),a)") "(2(a,i",oTurn,"),2(a,i",oPart,"))"
  end if

  write(lout,trackFmt) "TRACKING> "//trim(trackMode)//": Turn ",n," / ",numl,", Particles: ",napx," / ",napxo
  flush(lout)

end subroutine trackReport

!-----------------------------------------------------------------------
!
!  F. SCHMIDT
!-----------------------------------------------------------------------
!  3 February 1999
!-----------------------------------------------------------------------
subroutine dist1
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  implicit none
  integer ia,ib2,ib3,ie
  real(kind=fPrec) dam1
  save
!-----------------------------------------------------------------------
  do 20 ia=1,napx,2
    if(.not.pstop(partID(ia)).and..not.pstop(partID(ia)+1).and.(mod(partID(ia),2).ne.0)) then
      ie=ia+1
      dam(ia)=zero
      dam(ie)=zero
      xau(1,1)= xv1(ia)
      xau(1,2)= yv1(ia)
      xau(1,3)= xv2(ia)
      xau(1,4)= yv2(ia)
      xau(1,5)=sigmv(ia)
      xau(1,6)= dpsv(ia)
      xau(2,1)= xv1(ie)
      xau(2,2)= yv1(ie)
      xau(2,3)= xv2(ie)
      xau(2,4)= yv2(ie)
      xau(2,5)=sigmv(ie)
      xau(2,6)= dpsv(ie)
      cloau(1)= clo6v(1)
      cloau(2)=clop6v(1)
      cloau(3)= clo6v(2)
      cloau(4)=clop6v(2)
      cloau(5)= clo6v(3)
      cloau(6)=clop6v(3)
      di0au(1)= di0xs
      di0au(2)=dip0xs
      di0au(3)= di0zs
      di0au(4)=dip0zs
      tau(:,:)=tasau(:,:)

      call distance(xau,cloau,di0au,tau,dam1)
      dam(ia)=dam1
      dam(ie)=dam1
    endif
20 continue
  return
end subroutine dist1

!-----------------------------------------------------------------------
!
!  F. SCHMIDT
!  Write to fort.12
!  TODO: Not using chr_fromreal as it should
!-----------------------------------------------------------------------
!  3 February 1999
!-----------------------------------------------------------------------
subroutine write6(n)

  use crcoall
  use mod_common
  use mod_common_main
  use read_write
  use mod_settings

  implicit none

  integer ia,ig,n

  call writeFort12

  do ia=1,napxo,2
    ig=ia+1
#ifndef CR
#ifndef STF
    flush(91-(ig/2))
#else
    flush(90)
#endif
#endif
    !-- PARTICLES STABLE (Only if QUIET < 2)
    if(.not.pstop(ia).and..not.pstop(ig)) then
      if(st_quiet < 2) write(lout,10000) ia,izu0,dpsv(ia),n
      if(st_quiet < 1) write(lout,10010)                    &
        xv1(ia),yv1(ia),xv2(ia),yv2(ia),sigmv(ia),dpsv(ia), &
        xv1(ig),yv1(ig),xv2(ig),yv2(ig),sigmv(ig),dpsv(ig), &
        e0,ejv(ia),ejv(ig)
    end if
  end do
  return
10000 format(1x/5x,'PARTICLE ',i7,' RANDOM SEED ',i8,' MOMENTUM DEVIATION ',g12.5 /5x,'REVOLUTION ',i8/)
10010 format(10x,f47.33)
end subroutine write6
