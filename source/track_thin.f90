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
  use tracking

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
  use mod_utils
  use wire
#ifdef CR
  use checkpoint_restart
#endif
#ifdef BOINC
  use mod_boinc
#endif

  implicit none

  integer i,irrtr,ix,j,k,n,nmz,nthinerr,xory,nac,nfree,nramp1,nplato,nramp2,kxxa,nfirst
  real(kind=fPrec) pz,cccc,cikve,crkve,crkveuk,r0,stracki,xlvj,yv1j,yv2j,zlvj,acdipamp,qd,acphase,  &
    acdipamp2,acdipamp1,crabamp,crabfreq,kcrab,RTWO,NNORM,l,cur,dx,dy,tx,ty,embl,chi,xi,yi,dxi,dyi, &
    rrelens,frrelens,xelens,yelens,onedp,fppsig,tan_t,sin_t,cos_t,costh_temp,sinth_temp,pxf,pyf,    &
    r_temp,z_temp,sigf,q_temp,pttemp,xlv,zlv,temp_angle
  logical llost
  real(kind=fPrec) crkveb(npart),cikveb(npart),rho2b(npart),tkb(npart),rb(npart),rkb(npart),        &
    xrb(npart),zrb(npart),xbb(npart),zbb(npart),crxb(npart),crzb(npart),cbxb(npart),cbzb(npart)
  real(kind=fPrec) :: krf, x_t, y_t
  complex(kind=fPrec) :: Cp0, Sp1
  complex(kind=fPrec), parameter :: imag=(zero,one)

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
    call trackBeginTurn(n, nthinerr)
    if(nthinerr /= 0) return

    ! loop over structure elements, single element: name + type + parameter,
    ! structure element = order of single elements/blocks
    do 630 i=1,iu
      ! No if(ktrack(i).eq.1) - a BLOC - is needed in thin tracking,
      ! as no dependency on ix in this case.
      ix=ic(i)-nblo ! ix = index of single element
      meta_nPTurnEle = meta_nPTurnEle + napx

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
        if(iexact) then ! EXACT DRIFT
          do j=1,napx
            pz     = sqrt(c1e6 - (yv1(j)**2 + yv2(j)**2))*c1m3 ! pz/p0
            xv1(j) = xv1(j) + stracki*(yv1(j)/pz)
            xv2(j) = xv2(j) + stracki*(yv2(j)/pz)
          end do
        else
          do j=1,napx
            xv1(j) = xv1(j) + stracki*yv1(j)
            xv2(j) = xv2(j) + stracki*yv2(j)
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
          ! The values are stored in the temp vector which are used for the multiplication.
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
          call part_updatePartEnergy(1,.false.)

          ! We have to go back to angles after we updated the energy.
          yv1(j) = yv1(j)*moidpsv(j)
          yv2(j) = yv2(j)*moidpsv(j)
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
#include "include/beambeam41.f90"
        goto 620
      case (42)
#include "include/beambeam42.f90"
        goto 620
      case (43)
#include "include/beambeam43.f90"
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

    if(nthinerr /= 0) return
    if(ntwin /= 2) call trackDistance
#ifndef FLUKA
    if(mod(n,nwr(4)) == 0) call trackPairReport(n)
#else
    ! increase napxto, to get an estimation of particles*turns
    napxto = napxto + napx
#endif
    firstrun = .false.

  640 continue

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
  use tracking

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
  use coll_db
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
  use mod_utils
  use wire
#ifdef CR
  use checkpoint_restart
#endif
#ifdef BOINC
  use mod_boinc
#endif

  implicit none

  integer i,irrtr,ix,j,k,n,nmz,nthinerr,dotrack,xory,nac,nfree,nramp1,nplato,nramp2,kxxa,nfirst
  real(kind=fPrec) pz,cccc,cikve,crkve,crkveuk,r0,stracki,xlvj,yv1j,yv2j,zlvj,acdipamp,qd,          &
    acphase,acdipamp2,acdipamp1,crabamp,crabfreq,crabamp2,crabamp3,crabamp4,kcrab,RTWO,NNORM,l,cur, &
    dx,dy,tx,ty,embl,chi,xi,yi,dxi,dyi,rrelens,frrelens,xelens,yelens, onedp,fppsig,costh_temp,     &
    sinth_temp,tan_t,sin_t,cos_t,pxf,pyf,r_temp,z_temp,sigf,q_temp,pttemp,xlv,zlv,temp_angle
  logical llost, doFField, is_coll
  real(kind=fPrec) crkveb(npart),cikveb(npart),rho2b(npart),tkb(npart),rb(npart),rkb(npart),        &
    xrb(npart),zrb(npart),xbb(npart),zbb(npart),crxb(npart),crzb(npart),cbxb(npart),cbzb(npart)
  real(kind=fPrec) :: krf, x_t, y_t
  complex(kind=fPrec) :: Cp0, Sp1
  complex(kind=fPrec), parameter :: imag=(zero,one)

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
    call trackBeginTurn(n, nthinerr)
    if(nthinerr /= 0) return

    if(do_coll) then
      call coll_startTurn(n)
    end if

    !! This is the loop over each element: label 650
    do 650 i=1,iu !Loop over elements

      ! No if(ktrack(i).eq.1) - a BLOC - is needed in thin tracking,
      ! as no dependency on ix in this case.
      ix=ic(i)-nblo

      if(do_coll) then
        call coll_startElement(i,ix)
      end if
      meta_nPTurnEle = meta_nPTurnEle + napx

      ! Fringe Fields
      if(ffield_enabled .and. ix > 0) then
        doFField = FFindex(ix) > 0
      else
        doFField = .false.
      end if

#ifdef BEAMGAS
      !YIL Call beamGas subroutine whenever a pressure-element is found
      ! should be faster/safer to first check the turn then do the name search
      if(iturn == 1 ) then
        if(bez(c_ix)(1:5).eq.'PRESS' .or.  bez(c_ix)(1:5).eq.'press' ) then
          call beamGas(c_ix,nhit_type,dcum(i),myenom)
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

      ! The below splitting of if-statements is needed to prevent out of bounds error
      ! when building with gfortran/debug
      is_coll = .false.
      if(do_coll) then
        if(cdb_elemMap(myix) > 0) then
          is_coll = .true.
        end if
      end if

      if(is_coll) then
        dotrack = 1
      else
        dotrack = ktrack(i)
      end if

      select case(dotrack)
      case (1)
        stracki = strack(i)
        ! Check if collimation is enabled, and call the collimation code as necessary
        if(do_coll .and. is_coll) then
          ! Collimator is in database, and we're doing collimation
          call collimate_trackThin(stracki,.true.)
        else ! Normal SixTrack drifts
          if(iexact) then ! EXACT DRIFT
            do j=1,napx
              pz       = sqrt(c1e6 - (yv1(j)**2 + yv2(j)**2))*c1m3 ! pz/p0
              xv1(j)   = xv1(j)   + stracki*(yv1(j)/pz)
              xv2(j)   = xv2(j)   + stracki*(yv2(j)/pz)
              sigmv(j) = sigmv(j) + stracki*((one - rvv(j)/pz)*c1e3)
            end do
          else
            do j=1,napx
              xv1(j)   = xv1(j)   + stracki*yv1(j)
              xv2(j)   = xv2(j)   + stracki*yv2(j)
              sigmv(j) = sigmv(j) + stracki*(c1e3-rvv(j)*(c1e3+(yv1(j)**2+yv2(j)**2)*c5m4))
            end do
          end if
          if(do_coll) then
            ! Not a collimator, but collimation still need to perform additional calculations
            call collimate_trackThin(stracki,.false.)
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
          ejfv(j)=sqrt(ejv(j)**2-nucm(j)**2)
          rvv(j)=(ejv(j)*e0f)/(e0*ejfv(j))
          dpsv(j)=(ejfv(j)*(nucm0/nucm(j))-e0f)/e0f
          oidpsv(j)=one/(one+dpsv(j))
          moidpsv(j)=mtc(j)/(one+dpsv(j))
          omoidpsv(j)=c1e3*((one-mtc(j))*oidpsv(j))
          dpsv1(j)=(dpsv(j)*c1e3)*oidpsv(j)
          yv1(j)=(ejf0v(j)/ejfv(j))*yv1(j)
          yv2(j)=(ejf0v(j)/ejfv(j))*yv2(j)
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
#include "include/beambeam41.f90"
        goto 640
      case (42)
#include "include/beambeam42.f90"
        goto 640
      case (43)
#include "include/beambeam43.f90"
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

      if(do_coll) then
        call coll_endElement
      end if

#include "include/lostpart.f90"

645   continue

      if(.not. ldumpfront) then
        call dump_lines(n,i,ix)
      end if

650 continue !END loop over structure elements

    if(do_coll) then
      call coll_endTurn
    end if
#ifdef ROOT
    if(root_flag .and. root_Collimation == 1) then
      call SurvivalRootWrite(n, napx)
    end if
#endif

    if(nthinerr /= 0) return
    if(do_coll .eqv. .false.) then
      if(ntwin /= 2) call trackDistance
#ifndef FLUKA
      if(mod(n,nwr(4)) == 0) call trackPairReport(n)
#endif
    end if
#ifdef FLUKA
    ! increase napxto, to get an estimation of particles*turns
    napxto = napxto + napx
#endif
    firstrun = .false.

660 continue !END loop over turns

end subroutine thin6d
