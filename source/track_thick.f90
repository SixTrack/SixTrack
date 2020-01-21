!-----------------------------------------------------------------------
!
!  TRACK THICK LENS 4D
!
!  F. SCHMIDT
!-----------------------------------------------------------------------
subroutine thck4d(nthinerr)
  ! Replaced computed goto with select case, VKBO 27/11/2017
  use floatPrecision
  use string_tools
  use physical_constants
  use mathlib_bouncer
  use numerical_constants
  use mod_particles
  use bdex, only : bdex_enable
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

  use mod_settings
  use mod_meta
  use postprocessing, only : writebin
  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  use elens, only : elens_ktrack, elens_kick
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

  integer i,idz1,idz2,irrtr,ix,j,jb,jmel,jx,k,n,nmz,nthinerr,nac,nfree,nramp1,nplato,nramp2,   &
    kxxa,nfirst
  real(kind=fPrec) cccc,cikve,crkve,crkveuk,puxve,puxve1,puxve2,puzve1,puzve2,puzve,r0,xlvj,yv1j,   &
    yv2j,zlvj,acdipamp,qd,acphase, acdipamp2,acdipamp1,crabamp,crabfreq,kcrab,RTWO,NNORM,l,cur,dx,  &
    dy,tx,ty,embl,chi,xi,yi,dxi,dyi,rrelens,frrelens,xelens,yelens,onedp,fppsig,costh_temp,         &
    sinth_temp,pxf,pyf,q_temp,xlv,zlv

  logical llost
  real(kind=fPrec) crkveb(npart),cikveb(npart),rho2b(npart),tkb(npart),rb(npart),rkb(npart),        &
    xrb(npart),zrb(npart),xbb(npart),zbb(npart),crxb(npart),crzb(npart),cbxb(npart),cbzb(npart)

#ifdef FLUKA
  logical recompute_linear_matrices
#endif

  idz1=idz(1)
  idz2=idz(2)

#ifdef CR
  if(cr_restart) then
    call crstart
    write(crlog,"(2(a,i0))") "TRACKING> Thick 4D restarting on turn ",cr_numl," / ",numl
  end if
  nnuml  = numl
  nfirst = cr_numl
#else
  nfirst = 1
#endif
  do 490 n=nfirst,numl
    call trackBeginTurn(n, nthinerr)
    if(nthinerr /= 0) return

    do 480 i=1,iu
      if(ktrack(i).eq.1) then
        ix=ic(i)
      else
        ix=ic(i)-nblo
        meta_nPTurnEle = meta_nPTurnEle + napx

        if(ldumpfront) then
          write(lout,"(a)") "TRACKING> DUMP/FRONT not yet supported on thick elements "//&
            "due to lack of test cases. Please contact developers!"
          call prror
        end if

      end if

#ifdef FLUKA
!           A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!           last modified: 17-07-2013
!           is the current entry an instance of a FLUKA element?
!           inserted in main code by the 'fluka' compilation flag
      if (fluka_enable) then
        if(ktrack(i).ne.1) then ! Skip BLOCs, FLUKA elements must
                                !      be SINGLE ELEMENTs
          if(fluka_type(ix).ne.FLUKA_NONE) then
            if(fluka_type(ix).eq.FLUKA_ELEMENT) then
              call kernel_fluka_element( n, i, ix )
!                   re-compute transport matrices of linear elements,
!                      according to momentum of surviving/new particles
              recompute_linear_matrices = .true.
              ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
              ! last modified: 07-03-2018
              ! store old particle coordinates
              if (lbacktracking) call aperture_saveLastCoordinates(i,ix,0)
              goto 470
            else if(fluka_type(ix).eq.FLUKA_ENTRY) then
              fluka_inside = .true.
              call kernel_fluka_entrance( n, i, ix )
              goto 475
            else if(fluka_type(ix).eq.FLUKA_EXIT) then
              fluka_inside = .false.
              call kernel_fluka_exit
              ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
              ! last modified: 07-03-2018
              ! store old particle coordinates
              if (lbacktracking) call aperture_saveLastCoordinates(i,ix,0)
              ! re-compute transport matrices of linear elements,
              !    according to momentum of surviving/new particles
              recompute_linear_matrices = .true.
              goto 470
            end if
          end if
        end if
        if(fluka_inside) then
          if(fluka_debug) then
            write(lout,"(a,i0)") "FLUKA> Skipping lattice element at ",i
            write(fluka_log_unit,*) '# Skipping lattice element at ', i
          end if
          goto 480
        end if
      end if
#endif

            if (bdex_enable) then
               !TODO - if you have a test case, please contact developers!
               write(lout,"(a)") "BDEX> BDEX only available for thin6d"
               call prror
            endif

!----------count=43
      select case (ktrack(i))
      case (1)
        if (lbacktracking) then
          jmel=mel(ix)
          do jb=1,jmel
            jx=mtyp(ix,jb)
            do j=1,napx
#include "include/thcklin.f90"
            end do
          end do
        else
          do j=1,napx
            puxve=xv1(j)
            puzve=yv1(j)
            xv1(j)=bl1v(1,1,j,ix)*puxve+bl1v(2,1,j,ix)*puzve+((real(idz1,fPrec)*bl1v(5,1,j,ix))*dpsv(j))*c1e3
            yv1(j)=bl1v(3,1,j,ix)*puxve+bl1v(4,1,j,ix)*puzve+((real(idz1,fPrec)*bl1v(6,1,j,ix))*dpsv(j))*c1e3
            puxve=xv2(j)
            puzve=yv2(j)
            xv2(j)=bl1v(1,2,j,ix)*puxve+bl1v(2,2,j,ix)*puzve+((real(idz2,fPrec)*bl1v(5,2,j,ix))*dpsv(j))*c1e3
            yv2(j)=bl1v(3,2,j,ix)*puxve+bl1v(4,2,j,ix)*puzve+((real(idz2,fPrec)*bl1v(6,2,j,ix))*dpsv(j))*c1e3
          end do
        end if
        goto 480

      case (2) ! RF Cavity
        goto 480

      case (3) ! Phase Trombone
        irrtr = imtr(ix)
        do j=1,napx
          temptr(1) = xv1(j)
          temptr(2) = yv1(j)
          temptr(3) = xv2(j)
          temptr(4) = yv2(j)

          xv1(j) = cotr(irrtr,1)
          yv1(j) = cotr(irrtr,2)
          xv2(j) = cotr(irrtr,3)
          yv2(j) = cotr(irrtr,4)

          do kxxa=1,4
            xv1(j) = xv1(j) + temptr(kxxa)*rrtr(irrtr,1,kxxa)
            yv1(j) = yv1(j) + temptr(kxxa)*rrtr(irrtr,2,kxxa)
            xv2(j) = xv2(j) + temptr(kxxa)*rrtr(irrtr,3,kxxa)
            yv2(j) = yv2(j) + temptr(kxxa)*rrtr(irrtr,4,kxxa)
          end do
        end do
        goto 470

      case (4,5,6,7,8,9,10)
        goto 480
      case (11) ! HORIZONTAL DIPOLE
        do j=1,napx
#include "include/kickv01h.f90"
        end do
        goto 470
      case (12) ! NORMAL QUADRUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvxxh.f90"
        end do
        goto 470
      case (13) ! NORMAL SEXTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 470
      case (14) ! NORMAL OCTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 470
      case (15) ! NORMAL DECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 470
      case (16) ! NORMAL DODECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 470
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
        goto 470
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
        goto 470
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
        goto 470
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
        goto 470
      case (21) ! VERTICAL DIPOLE
        do j=1,napx
#include "include/kickv01v.f90"
        end do
        goto 470
      case (22) ! SKEW QUADRUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvxxv.f90"
        end do
        goto 470
      case (23) ! SKEW SEXTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 470
      case (24) ! SKEW OCTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 470
      case (25) ! SKEW DECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 470
      case (26) ! SKEW DODECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 470
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
        goto 470
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
        goto 470
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
        goto 470
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
        goto 470
      case (31)
        goto 470
      case (32)
        goto 240
      case (33)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v01.f90"
        end do
        goto 470
      case (34)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v01.f90"
        end do
        goto 240
      case (35)
        do j=1,napx
#include "include/mul4v02.f90"
        end do
        goto 470
      case (36)
        do j=1,napx
#include "include/mul4v02.f90"
        end do
        goto 240
      case (37)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v03.f90"
        end do
        goto 470
      case (38)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v03.f90"
        end do
        goto 240
      case (39)
        do j=1,napx
#include "include/mul4v04.f90"
        end do
        goto 470
      case (40)
        do j=1,napx
#include "include/mul4v04.f90"
        end do
        goto 240
      case (41)
#include "include/beambeam41.f90"
        goto 470
      case (42)
#include "include/beambeam42.f90"
        goto 470
      case (43)
#include "include/beambeam43.f90"
        goto 470
      case (44,46,47,48,49,50,57,58,59,60,61,62)
        goto 480
      case (45) ! Wire
#include "include/wirekick.f90"
        goto 470
      case (51)
#include "include/acdipkick1.f90"
        goto 470
      case (52)
#include "include/acdipkick2.f90"
        goto 470
      case (53)
#include "include/crabkick1.f90"
        goto 470
      case (54)
#include "include/crabkick2.f90"
        goto 470
      case (55) ! DIPEDGE ELEMENT
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvdpe.f90"
        end do
        goto 470
      case (56) ! Solenoid
        do j=1,napx
#include "include/kickvso1.f90"
        end do
        goto 470
      case (elens_ktrack) ! Elens
        call elens_kick(i,ix,n)
        goto 470
      case (cheby_ktrack) ! Chebyshev lens
        call cheby_kick(i,ix,n)
        goto 470
      end select
      goto 480

240   r0=ek(ix)
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
      goto 470

!----------------------------

470   continue

#include "include/lostpart.f90"

#ifdef FLUKA
      ! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
      ! last modified: 17-07-2013
      ! re-compute transport matrices of linear elements
      ! inserted in main code by the 'fluka' compilation flag
      if ( recompute_linear_matrices ) then
        ! after a FLUKA element: additional particles may have been generated
        call envarsv
        recompute_linear_matrices = .false.
      else if ( llost ) then
        ! after any other element: no additional particles, thus update only momentum-dependent matrix elements
        call synuthck
      else
        goto 475
      end if
      ! recompute matrices of BLOCKs
      call blocksv
#endif
#ifndef FLUKA
      if(llost) then
        call synuthck
      end if
#endif

475   continue

      if (.not. ldumpfront) then
        call dump_lines(n,i,ix)
      endif

480 continue

    if(nthinerr /= 0) return
    if(ntwin /= 2) call trackDistance
#ifndef FLUKA
    if(mod(n,nwr(4)) == 0) call trackPairReport(n)
#else
    ! increase napxto, to get an estimation of particles*turns
    napxto = napxto + napx
#endif
    firstrun = .false.

490 continue

end subroutine thck4d

!-----------------------------------------------------------------------
!
!  TRACK THICK LENS 6D
!
!  F. SCHMIDT
!-----------------------------------------------------------------------
subroutine thck6d(nthinerr)
  ! Replaced computed goto with select case, VKBO 28/11/2017
  use floatPrecision
  use string_tools
  use physical_constants
  use mathlib_bouncer
  use numerical_constants
  use mod_particles
  use bdex, only : bdex_enable
  use dump, only : dump_linesFirst, dump_lines, ldumpfront
  use collimation, only: do_coll, part_abs_turn
  use aperture
  use tracking

#ifdef FLUKA
!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 17-07-2013
!     import mod_fluka
!     inserted in main code by the 'fluka' compilation flag
  use mod_fluka
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
  use aperture
  use elens, only : elens_ktrack, elens_kick
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

  integer i,idz1,idz2,irrtr,ix,j,jb,jmel,jx,k,n,nmz,nthinerr,nac,nfree,nramp1,nplato,nramp2,   &
    kxxa,nfirst
  real(kind=fPrec) cccc,cikve,crkve,crkveuk,puxve1,puxve2,puzve1,puzve2,r0,xlvj,yv1j,yv2j,zlvj,     &
    acdipamp,qd,acphase,acdipamp2,acdipamp1,crabamp,crabfreq,kcrab,RTWO,NNORM,l,cur,dx,dy,tx,ty,    &
    embl,chi,xi,yi,dxi,dyi,rrelens,frrelens,xelens,yelens,onedp,fppsig,costh_temp,sinth_temp,pxf,   &
    pyf,r_temp,z_temp,sigf,q_temp,pttemp,xlv,zlv
  logical llost
  real(kind=fPrec) crkveb(npart),cikveb(npart),rho2b(npart),tkb(npart),rb(npart),rkb(npart),        &
    xrb(npart),zrb(npart),xbb(npart),zbb(npart),crxb(npart),crzb(npart),cbxb(npart),cbzb(npart)

#ifdef FLUKA
  logical recompute_linear_matrices
#endif

  idz1=idz(1)
  idz2=idz(2)

#ifdef FLUKA
!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 17-07-2013
!     force re-computation of transport matrices of linear elements
!     inserted in main code by the 'fluka' compilation flag
  recompute_linear_matrices = .false.
#endif

  ! Now the outer loop over turns
#ifdef CR
  if(cr_restart) then
    call crstart
    write(crlog,"(2(a,i0))") "TRACKING> Thick 6D restarting on turn ",cr_numl," / ",numl
  end if
  nnuml  = numl
  nfirst = cr_numl
#else
  nfirst = 1
#endif
  do 510 n=nfirst,numl
    call trackBeginTurn(n, nthinerr)
    if(nthinerr /= 0) return

    do 500 i=1,iu
      if(ktrack(i).eq.1) then
        ix=ic(i)
      else
        ix=ic(i)-nblo
      end if
      meta_nPTurnEle = meta_nPTurnEle + napx

      if (ldumpfront) then
        write(lerr,"(a)") "DUMP> ERROR FRONT not yet supported on thick elements due to lack of test cases. "//&
          "Please contact developers!"
        call prror
      end if

#ifdef FLUKA
!           A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!           last modified: 17-07-2013
!           is the current entry an instance of a FLUKA element?
!           inserted in main code by the 'fluka' compilation flag
      if (fluka_enable) then
        if(ktrack(i).ne.1) then ! Skip BLOCs, FLUKA elements must be SINGLE ELEMENTs
          if(fluka_type(ix).ne.FLUKA_NONE) then
            if(fluka_type(ix).eq.FLUKA_ELEMENT) then
              call kernel_fluka_element( n, i, ix )
              ! Re-compute transport matrices of linear elements, according to momentum of surviving/new particles
              recompute_linear_matrices = .true.
              ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
              ! last modified: 07-03-2018
              ! store old particle coordinates
              if (lbacktracking) call aperture_saveLastCoordinates(i,ix,0)
              goto 490
            else if(fluka_type(ix).eq.FLUKA_ENTRY) then
              fluka_inside = .true.
              call kernel_fluka_entrance( n, i, ix )
              goto 495
            else if(fluka_type(ix).eq.FLUKA_EXIT) then
              fluka_inside = .false.
              call kernel_fluka_exit
              ! Re-compute transport matrices of linear elements, according to momentum of surviving/new particles
              recompute_linear_matrices = .true.
              ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
              ! last modified: 07-03-2018
              ! store old particle coordinates
              if (lbacktracking) call aperture_saveLastCoordinates(i,ix,0)
              goto 490
            end if
          end if
        end if
        if(fluka_inside) then
          if(fluka_debug) then
            write(lout,"(a,i0)") "FLUKA> Skipping lattice element at ",i
            write(fluka_log_unit,*) '# Skipping lattice element at ', i
          end if
          goto 500
        end if
      end if
#endif

            if (bdex_enable) then
               !TODO - if you have a test case, please contact developers!
               write(lout,"(a)") "BDEX> BDEX only available for thin6d"
               call prror
            endif

!----------count 44
!----------count 54! Eric
      select case (ktrack(i))
      case (1)
        jmel=mel(ix)
        do jb=1,jmel
          jx=mtyp(ix,jb)
          do j=1,napx
#include "include/thcklin.f90"
          end do
        end do
        goto 500

      case (2) ! RF Cavity
        if(abs(dppoff) > pieni) then
          sigmv(1:napx) = sigmv(1:napx) - sigmoff(i)
        end if
        if(abs(kz(ix)) == 12) then
          do j=1,napx
            ejv(j) = ejv(j) + (ed(ix)*sin_mb(hsyc(ix)*sigmv(j)+phasc(ix)))*nqq(j)
          end do
        else
          do j=1,napx
            ejv(j) = ejv(j) + (hsy(1)*sin_mb(hsy(3)*sigmv(j)))*nqq(j)
          end do
        end if
        call part_updatePartEnergy(1,.true.)
        goto 490

      case (3) ! Phase Trombone
        irrtr = imtr(ix)
        do j=1,napx
          ! The values are stored in the temp vector which are used for the multiplication.
          temptr(1) = xv1(j)
          temptr(2) = yv1(j)/moidpsv(j)
          temptr(3) = xv2(j)
          temptr(4) = yv2(j)/moidpsv(j)
          temptr(5) = sigmv(j)
          temptr(6) = ((mtc(j)*ejv(j)-e0)/e0f)*c1e3*(e0/e0f)

          ! Adding the closed orbit. The previous values are stored in the temptr vector.
          xv1(j)   = cotr(irrtr,1)
          yv1(j)   = cotr(irrtr,2)
          xv2(j)   = cotr(irrtr,3)
          yv2(j)   = cotr(irrtr,4)
          sigmv(j) = cotr(irrtr,5)
          pttemp   = cotr(irrtr,6)

          ! Multiplying the arbitrary matrix to the coordinates.
          do kxxa=1,6
            xv1(j)   = xv1(j)   + temptr(kxxa)*rrtr(irrtr,1,kxxa)
            yv1(j)   = yv1(j)   + temptr(kxxa)*rrtr(irrtr,2,kxxa)
            xv2(j)   = xv2(j)   + temptr(kxxa)*rrtr(irrtr,3,kxxa)
            yv2(j)   = yv2(j)   + temptr(kxxa)*rrtr(irrtr,4,kxxa)
            sigmv(j) = sigmv(j) + temptr(kxxa)*rrtr(irrtr,5,kxxa)
            pttemp   = pttemp   + temptr(kxxa)*rrtr(irrtr,6,kxxa)
          end do

          ! Transforming back to the tracked coordinates of SixTrack
          ejv(j) = (e0f*pttemp/(c1e3*(e0/e0f))+e0)/mtc(j)
        end do
        call part_updatePartEnergy(1,.false.)

        ! We have to go back to angles after we updated the energy.
        yv1(1:napx) = yv1(1:napx)*moidpsv(1:napx)
        yv2(1:napx) = yv2(1:napx)*moidpsv(1:napx)

        goto 490
      case (4,5,6,7,8,9,10)
        goto 500
      case (11) ! HORIZONTAL DIPOLE
        do j=1,napx
#include "include/kickv01h.f90"
        end do
        goto 490
      case (12) ! NORMAL QUADRUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvxxh.f90"
        end do
        goto 490
      case (13) ! NORMAL SEXTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 490
      case (14) ! NORMAL OCTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 490
      case (15) ! NORMAL DECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 490
      case (16) ! NORMAL DODECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxh.f90"
        end do
        goto 490
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
        goto 490
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
        goto 490
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
        goto 490
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
        goto 490
      case (21) ! VERTICAL DIPOLE
        do j=1,napx
#include "include/kickv01v.f90"
        end do
        goto 490
      case (22) ! SKEW QUADRUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvxxv.f90"
        end do
        goto 490
      case (23) ! SKEW SEXTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 490
      case (24) ! SKEW OCTUPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 490
      case (25) ! SKEW DECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 490
      case (26) ! SKEW DODECAPOLE
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvho.f90"
#include "include/kickvxxv.f90"
        end do
        goto 490
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
        goto 490
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
        goto 490
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
        goto 490
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
        goto 490
      case (31)
        goto 490
      case (32)
        goto 260
      case (33)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v01.f90"
#include "include/mul6v01.f90"
        end do
        goto 490
      case (34)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v01.f90"
#include "include/mul6v01.f90"
        end do
        goto 260
      case (35)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v02.f90"
#include "include/mul6v01.f90"
        end do
        goto 490
      case (36)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v02.f90"
#include "include/mul6v01.f90"
        end do
        goto 260
      case (37)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v03.f90"
#include "include/mul6v02.f90"
        end do
        goto 490
      case (38)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v03.f90"
#include "include/mul6v02.f90"
        end do
        goto 260
      case (39)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v04.f90"
#include "include/mul6v02.f90"
        end do
        goto 490
      case (40)
        do j=1,napx
#include "include/alignvb.f90"
#include "include/mul4v04.f90"
#include "include/mul6v02.f90"
        end do
        goto 260
      case (41)
#include "include/beambeam41.f90"
        goto 490
      case (42)
#include "include/beambeam42.f90"
        goto 490
      case (43)
#include "include/beambeam43.f90"
        goto 490
      case (44)
#include "include/beam6d.f90"
        goto 490
      case (45) ! Wire
#include "include/wirekick.f90"
        goto 490
      case (46,47,48,49,50,57,58,59,60,61,62)
        goto 500
      case (51)
#include "include/acdipkick1.f90"
        goto 490
      case (52)
#include "include/acdipkick2.f90"
        goto 490
      case (53)
#include "include/crabkick1.f90"
        goto 490
      case (54)
#include "include/crabkick2.f90"
        goto 490
      case (55) ! DIPEDGE ELEMENT
        do j=1,napx
#include "include/alignva.f90"
#include "include/kickvdpe.f90"
        end do
        goto 490
      case (56) ! Solenoid
        do j=1,napx
#include "include/kickvso1.f90"
#include "include/kickvso2.f90"
        end do
        goto 490
      case (elens_ktrack) ! Elens
        call elens_kick(i,ix,n)
        goto 490
      case (cheby_ktrack) ! Chebyshev lens
        call cheby_kick(i,ix,n)
        goto 490
      end select
      goto 500

260   r0=ek(ix)
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
      goto 490

490   continue

#include "include/lostpart.f90"

#ifdef FLUKA
      ! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
      ! last modified: 17-07-2013
      ! re-compute transport matrices of linear elements
      ! inserted in main code by the 'fluka' compilation flag
      if ( recompute_linear_matrices ) then
        ! after a FLUKA element: additional particles may have been generated
        call envarsv
        recompute_linear_matrices = .false.
      else if ( llost ) then
        ! after any other element: no additional particles, thus update only momentum-dependent matrix elements
        call synuthck
      end if
#endif
#ifndef FLUKA
      if(llost) then
        call synuthck
      end if
#endif

495   continue

      if (.not. ldumpfront) then
        call dump_lines(n,i,ix)
      end if

500 continue
! End of loop over elements

    if(nthinerr /= 0) return
    if(ntwin /= 2) call trackDistance
#ifndef FLUKA
    if(mod(n,nwr(4)) == 0) call trackPairReport(n)
#else
    ! increase napxto, to get an estimation of particles*turns
    napxto = napxto + napx
#endif
    firstrun = .false.

510 continue
! end loop over turns

end subroutine thck6d

!-----------------------------------------------------------------------
!
!  TRACK THICK LENS PART
!
!
!  F. SCHMIDT
!-----------------------------------------------------------------------
!  3 February 1999
!-----------------------------------------------------------------------
subroutine synuthck
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use numerical_constants
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  implicit none
  integer ih1,ih2,j,kz1,l
  real(kind=fPrec) fokm,fok,fok1,rho,si,co,sm1,sm2,sm3,sm12,sm23,as3,as4,as6,g,gl,rhoc,siq,aek,hi,  &
    fi,hi1,hp,hm,hc,hs,wf,afok,wfa,wfhi,rhoi,fokq

  do j=1,napx
    dpd(j)  = one+dpsv(j)
    dpsq(j) = sqrt(dpd(j))
  end do

  do 160 l=1,il
    if(abs(el(l)).le.pieni) goto 160
    kz1=kz(l)+1
!-----------------------------------------------------------------------
!  DRIFTLENGTH
!-----------------------------------------------------------------------
    if(kz1 == 1) then
      goto 20
!-----------------------------------------------------------------------
!  RECTANGULAR MAGNET
!  HORIZONTAL
!-----------------------------------------------------------------------
    elseif(kz1 == 2 .or. kz1 == 5) then
40    fokm=el(l)*ed(l)
      if(abs(fokm) <= pieni) goto 20
      if(kz1 == 2) then
        ih1=1
        ih2=2
      else
!  RECTANGULAR MAGNET VERTICAL
        ih1=2
        ih2=1
      endif
      do j=1,napx
        fok  = fokm/dpsq(j)
        rho  = (one/ed(l))*dpsq(j)
        fok1 = (tan_mb(fok*half))/rho
        si   = sin_mb(fok)
        co   = cos_mb(fok)
        al(2,ih1,j,l) = rho*si
        al(5,ih1,j,l) = ((-one*dpsv(j))*((rho*(one-co))/dpsq(j)))*c1e3
        al(6,ih1,j,l) = ((-one*dpsv(j))*((two*tan_mb(fok*half))/dpsq(j)))*c1e3

        sm1  = cos_mb(fok)
        sm2  = sin_mb(fok)*rho
        sm3  = -sin_mb(fok)/rho
        sm12 = el(l)-sm1*sm2
        sm23 = sm2*sm3
        as3  = (-one*rvv(j))*(((dpsv(j)*rho)/(two*dpsq(j)))*sm23-(rho*dpsq(j))*(one-sm1))
        as4  = ((-one*rvv(j))*sm23)/c2e3
        as6  = ((-one*rvv(j))*(el(l)+sm1*sm2))/c4e3
        as(1,ih1,j,l) = (el(l)*(one-rvv(j))-rvv(j)*((dpsv(j)**2/(four*dpd(j)))*sm12+dpsv(j)*(el(l)-sm2)))*c1e3
        as(2,ih1,j,l) = (-one*rvv(j))*((dpsv(j)/((two*rho)*dpsq(j)))* sm12-(sm2*dpsq(j))/rho)+fok1*as3
        as(3,ih1,j,l) = as3
        as(4,ih1,j,l) = as4+(two*as6)*fok1
        as(5,ih1,j,l) = ((-one*rvv(j))*sm12)/(c4e3*rho**2)+as6*fok1**2+fok1*as4
        as(6,ih1,j,l) = as6
!--VERTIKAL
        g  = tan_mb(fok*half)/rho
        gl = el(l)*g
        al(1,ih2,j,l) = one-gl
        al(3,ih2,j,l) = (-one*g)*(two-gl)
        al(4,ih2,j,l) = al(1,ih2,j,l)
        as6 = ((-one*rvv(j))*al(2,ih2,j,l))/c2e3
        as(4,ih2,j,l) = ((-one*two)*as6)*fok1
        as(5,ih2,j,l) = (as6*fok1)*fok1
        as(6,ih2,j,l) = as6
      end do
      goto 160
    elseif(kz1 == 4 .or. kz1 == 6) then
!-----------------------------------------------------------------------
!  SEKTORMAGNET
!  HORIZONTAL
!-----------------------------------------------------------------------
60    fokm=el(l)*ed(l)
      if(abs(fokm) <= pieni) goto 20
      if(kz1 == 4) then
        ih1 = 1
        ih2 = 2
      else
!  SECTOR MAGNET VERTICAL
        ih1 = 2
        ih2 = 1
      end if
      do j=1,napx
        fok  = fokm/dpsq(j)
        rho  = (one/ed(l))*dpsq(j)
        si   = sin_mb(fok)
        co   = cos_mb(fok)
        rhoc = (rho*(one-co))/dpsq(j)
        siq  = si/dpsq(j)
        al(1,ih1,j,l) = co
        al(2,ih1,j,l) = rho*si
        al(3,ih1,j,l) = (-one*si)/rho
        al(4,ih1,j,l) = co
        al(5,ih1,j,l) = ((-one*dpsv(j))*rhoc)*c1e3
        al(6,ih1,j,l) = ((-one*dpsv(j))*siq)*c1e3

        sm12 = el(l)-al(1,ih1,j,l)*al(2,ih1,j,l)
        sm23 = al(2,ih1,j,l)*al(3,ih1,j,l)
        as(1,ih1,j,l) = (el(l)*(one-rvv(j))-rvv(j)*((dpsv(j)**2/(four*dpd(j)))*sm12+dpsv(j)*(el(l)-al(2,ih1,j,l))))*c1e3
        as(2,ih1,j,l) = (-one*rvv(j))*((dpsv(j)/(two*rho*dpsq(j)))*sm12-dpd(j)*siq)
        as(3,ih1,j,l) = (-one*rvv(j))*(((dpsv(j)*rho)/(two*dpsq(j)))*sm23-dpd(j)*rhoc)
        as(4,ih1,j,l) = ((-one*rvv(j))*sm23)/c2e3
        as(5,ih1,j,l) = ((-one*rvv(j))*sm12)/((c4e3*rho)*rho)
        as(6,ih1,j,l) = ((-one*rvv(j))*(el(l)+al(1,ih1,j,l)*al(2,ih1,j,l)))/c4e3
!--VERTIKAL
        as(6,ih2,j,l) = ((-one*rvv(j))*al(2,ih2,j,l))/c2e3
      end do
      goto 160
    elseif(kz1 == 3) then
!-----------------------------------------------------------------------
!  QUADRUPOLE
!  FOCUSSING
!-----------------------------------------------------------------------
80   do j=1,napx
        fok = ek(l)*oidpsv(j)
        aek = abs(fok)
        hi  = sqrt(aek)
        fi  = el(l)*hi
        if(fok <= zero) then
          al(1,1,j,l) = cos_mb(fi)
          hi1 = sin_mb(fi)
          if(abs(hi) <= pieni) then
            al(2,1,j,l) = el(l)
          else
            al(2,1,j,l) = hi1/hi
          endif
          al(3,1,j,l) = -hi1*hi
          al(4,1,j,l) = al(1,1,j,l)
          as(1,1,j,l) = el(l)*(one-rvv(j))*c1e3
          as(4,1,j,l) = (((-one*rvv(j))*al(2,1,j,l))*al(3,1,j,l))/c2e3
          as(5,1,j,l) = (((-one*rvv(j))*(el(l)-al(1,1,j,l)*al(2,1,j,l)))*aek)/c4e3
          as(6,1,j,l) = ((-one*rvv(j))*(el(l)+al(1,1,j,l)*al(2,1,j,l)))/c4e3
!--DEFOCUSSING
          hp = exp_mb(fi)
          hm = one/hp
          hc = (hp+hm)*half
          hs = (hp-hm)*half
          al(1,2,j,l) = hc
          if(abs(hi) <= pieni) then
            al(2,2,j,l) = el(l)
          else
            al(2,2,j,l) = hs/hi
          endif
          al(3,2,j,l) = hs*hi
          al(4,2,j,l) = hc
          as(4,2,j,l) = (((-one*rvv(j))*al(2,2,j,l))*al(3,2,j,l))/c2e3
          as(5,2,j,l) = ((rvv(j)*(el(l)-al(1,2,j,l)*al(2,2,j,l)))*aek)/c4e3
          as(6,2,j,l) = ((-one*rvv(j))*(el(l)+al(1,2,j,l)*al(2,2,j,l)))/c4e3
        else
          al(1,2,j,l) = cos_mb(fi)
          hi1 = sin_mb(fi)
          if(abs(hi) <= pieni) then
            al(2,2,j,l) = el(l)
          else
            al(2,2,j,l) = hi1/hi
          endif
          al(3,2,j,l) = (-one*hi1)*hi
          al(4,2,j,l) = al(1,2,j,l)
          as(1,2,j,l) = (el(l)*(one-rvv(j)))*c1e3
          as(4,2,j,l) = (((-one*rvv(j))*al(2,2,j,l))*al(3,2,j,l))/c2e3
          as(5,2,j,l) = (((-one*rvv(j))*(el(l)-al(1,2,j,l)*al(2,2,j,l)))*aek)/c4e3
          as(6,2,j,l) = ((-one*rvv(j))*(el(l)+al(1,2,j,l)*al(2,2,j,l)))/c4e3
!--DEFOCUSSING
          hp = exp_mb(fi)
          hm = one/hp
          hc = (hp+hm)*half
          hs = (hp-hm)*half
          al(1,1,j,l) = hc
          if(abs(hi) <= pieni) then
            al(2,1,j,l) = el(l)
          else
            al(2,1,j,l) = hs/hi
          endif
          al(3,1,j,l) = hs*hi
          al(4,1,j,l) = hc
          as(4,1,j,l) = (((-one*rvv(j))*al(2,1,j,l))*al(3,1,j,l))/c2e3
          as(5,1,j,l) = ((rvv(j)*(el(l)-al(1,1,j,l)*al(2,1,j,l)))*aek)/c4e3
          as(6,1,j,l) = ((-one*rvv(j))*(el(l)+al(1,1,j,l)*al(2,1,j,l)))/c4e3
        endif
      end do
      goto 160
    elseif(kz1 == 7 .or. kz1 == 8) then
!-----------------------------------------------------------------------
!  COMBINED FUNCTION MAGNET HORIZONTAL
!  FOCUSSING
!-----------------------------------------------------------------------
100   if(kz1 == 7) then
        fokq = ek(l)
        ih1  = 1
        ih2  = 2
      else
!  COMBINED FUNCTION MAGNET VERTICAL
        fokq = -ek(l)
        ih1  = 2
        ih2  = 1
      end if
      do j=1,napx
        wf   = ed(l)/dpsq(j)
        fok  = fokq/dpd(j)-wf**2
        afok = abs(fok)
        hi   = sqrt(afok)
        fi   = hi*el(l)
        if(afok <= pieni) then
          as(6,1,j,l) = ((-one*rvv(j))*el(l))/c2e3
          as(6,2,j,l) = as(6,1,j,l)
          as(1,1,j,l) = (el(l)*(one-rvv(j)))*c1e3
        end if
        if(fok < (-one*pieni)) then
          si   = sin_mb(fi)
          co   = cos_mb(fi)
          wfa  = ((wf/afok)*(one-co))/dpsq(j)
          wfhi = ((wf/hi)*si)/dpsq(j)
          al(1,ih1,j,l) = co
          al(2,ih1,j,l) = si/hi
          al(3,ih1,j,l) = (-one*si)*hi
          al(4,ih1,j,l) = co
          al(5,ih1,j,l) = ((-one*wfa)*dpsv(j))*c1e3
          al(6,ih1,j,l) = ((-one*wfhi)*dpsv(j))*c1e3

          sm12 = el(l)-al(1,ih1,j,l)*al(2,ih1,j,l)
          sm23 = al(2,ih1,j,l)*al(3,ih1,j,l)
          as(1,ih1,j,l) = (el(l)*(one-rvv(j))-((rvv(j)*((dpsv(j)**2/(four*dpd(j)))&
            *sm12+dpsv(j)*(el(l)-al(2,ih1,j,l))))/afok)*wf**2)*c1e3
          as(2,ih1,j,l) = (-one*rvv(j))*(((dpsv(j)*wf)/(two*dpsq(j)))*sm12-dpd(j)*wfhi)
          as(3,ih1,j,l) = (-one*rvv(j))*(((((dpsv(j)*half)/afok)/dpd(j))*ed(l))*sm23-dpd(j)*wfa)
          as(4,ih1,j,l) = ((-one*rvv(j))*sm23)/c2e3
          as(5,ih1,j,l) = (((-one*rvv(j))*sm12)*afok)/c4e3
          as(6,ih1,j,l) = ((-one*rvv(j))*(el(l)+al(1,ih1,j,l)*al(2,ih1,j,l)))/c4e3

          aek = abs(ek(l)/dpd(j))
          hi  = sqrt(aek)
          fi  = hi*el(l)
          hp  = exp_mb(fi)
          hm  = one/hp
          hc  = (hp+hm)*half
          hs  = (hp-hm)*half
          al(1,ih2,j,l) = hc
          if(abs(hi) > pieni) al(2,ih2,j,l) = hs/hi
          al(3,ih2,j,l) = hs*hi
          al(4,ih2,j,l) = hc
          as(4,ih2,j,l) = (((-one*rvv(j))*al(2,ih2,j,l))*al(3,ih2,j,l))/c2e3
          as(5,ih2,j,l) = ((rvv(j)*(el(l)-al(1,ih2,j,l)*al(2,ih2,j,l)))*aek)/c4e3
          as(6,ih2,j,l) = ((-one*rvv(j))*(el(l)+al(1,ih2,j,l)*al(2,ih2,j,l)))/c4e3
        end if
!--DEFOCUSSING
        if(fok > pieni) then
          hp = exp_mb(fi)
          hm = one/hp
          hc = (hp+hm)*half
          hs = (hp-hm)*half
          al(1,ih1,j,l) = hc
          al(2,ih1,j,l) = hs/hi
          al(3,ih1,j,l) = hs*hi
          al(4,ih1,j,l) = hc

          wfa  = ((wf/afok)*(one-hc))/dpsq(j)
          wfhi = ((wf/hi)*hs)/dpsq(j)
          al(5,ih1,j,l) = (wfa*dpsv(j))*c1e3
          al(6,ih1,j,l) = ((-one*wfhi)*dpsv(j))*c1e3

          sm12 = el(l)-al(1,ih1,j,l)*al(2,ih1,j,l)
          sm23 = al(2,ih1,j,l)*al(3,ih1,j,l)
          as(1,ih1,j,l) = (((rvv(j)*((dpsv(j)**2/(four*dpd(j)))*sm12&
            +dpsv(j)*(el(l)-al(2,ih1,j,l))))/afok)*wf**2+el(l)*(one-rvv(j)))*c1e3
          as(2,ih1,j,l) = (-one*rvv(j))*(((dpsv(j)*wf)/(two*dpsq(j)))*sm12-dpd(j)*wfhi)
          as(3,ih1,j,l) = rvv(j)*(((((dpsv(j)*half)/afok)/dpd(j))* ed(l))*sm23-dpd(j)*wfa)
          as(4,ih1,j,l) = ((-one*rvv(j))*sm23)/c2e3
          as(5,ih1,j,l) = ((rvv(j)*sm12)*afok)/c4e3
          as(6,ih1,j,l) = ((-one*rvv(j))*(el(l)+al(1,ih1,j,l)*al(2,ih1,j,l)))/c4e3

          aek = abs(ek(l)/dpd(j))
          hi  = sqrt(aek)
          fi  = hi*el(l)
          si  = sin_mb(fi)
          co  = cos_mb(fi)
          al(1,ih2,j,l) = co
          al(2,ih2,j,l) = si/hi
          al(3,ih2,j,l) = (-one*si)*hi
          al(4,ih2,j,l) = co
          as(4,ih2,j,l) = (((-one*rvv(j))*al(2,ih2,j,l))*al(3,ih2,j,l))/c2e3
          as(5,ih2,j,l) = (((-one*rvv(j))*(el(l)-al(1,ih2,j,l)*al(2,ih2,j,l)))*aek)/c4e3
          as(6,ih2,j,l) = ((-one*rvv(j))*(el(l)+al(1,ih2,j,l)*al(2,ih2,j,l)))/c4e3
        end if
      end do
      goto 160
    elseif(kz1 == 9) then
!-----------------------------------------------------------------------
!  EDGE FOCUSSING
!-----------------------------------------------------------------------
140   do j=1,napx
        rhoi = ed(l)/dpsq(j)
        fok  = rhoi*tan_mb((el(l)*rhoi)*half)
        al(3,1,j,l) = fok
        al(3,2,j,l) = -fok
      end do
      goto 160
    else
!Eric
! Is really an error but old code went to 160
      goto 160
    end if
!-----------------------------------------------------------------------
!  DRIFTLENGTH
!-----------------------------------------------------------------------
20  do j=1,napx
      as(6,1,j,l) = ((-one*rvv(j))*el(l))/c2e3
      as(6,2,j,l) = as(6,1,j,l)
      as(1,1,j,l) = (el(l)*(one-rvv(j)))*c1e3
    end do
160 continue
!---------------------------------------  END OF 'ENVARS' (2)
  return
end subroutine synuthck
