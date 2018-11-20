!-----------------------------------------------------------------------
!
!  TRACK THICK LENS PART
!
!  F. SCHMIDT
!-----------------------------------------------------------------------
subroutine trauthck(nthinerr)
  ! Rewritten to remove computed gotos by V.K.B.Olsen on 20/11/2017
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use dynk, only : dynk_enabled, dynk_isused, dynk_pretrack

#ifdef FLUKA
! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
! last modified: 17-07-2013
! import mod_fluka
! inserted in main code by the 'fluka' compilation flag
  use mod_fluka
#endif

  use collimation
  use mod_time

  use crcoall
  use parpro
  use mod_common
  use mod_commonmn
  use mod_commons
  use mod_commont
  use mod_commond
  use mod_fluc, only : fluc_errAlign,fluc_writeFort4
  implicit none

  integer i,ix,j,jb,jj,jx,kpz,kzz,napx0,nbeaux,nmz,nthinerr
  real(kind=fPrec) benkcc,cbxb,cbzb,cikveb,crkveb,crxb,crzb,r0,r000,r0a,r2b,rb,rho2b,rkb,tkb,xbb,xrb,zbb,zrb
  dimension crkveb(npart),cikveb(npart),rho2b(npart),tkb(npart),r2b(npart),rb(npart),rkb(npart),&
  xrb(npart),zrb(npart),xbb(npart),zbb(npart),crxb(npart),crzb(npart),cbxb(npart),cbzb(npart),nbeaux(nbb)
  save

  if (do_coll) then
     write(lout,*) "Error: in trauthck and do_coll is TRUE"
     write(lout,*) "Collimation is not supported for thick tracking"
     call prror(-1)
  endif

  do i=1,npart
    nlostp(i)=i
  end do
  do i=1,nblz
    ktrack(i)=0
    strack(i)=zero
    strackc(i)=zero
    stracks(i)=zero
  end do
#include "include/beams1.f90"
  do 290 i=1,iu
    if(mout2.eq.1.and.i.eq.1) call fluc_writeFort4
    ix=ic(i)
    if(ix.le.nblo) then
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
    else if(kzz.eq.12) then
      ! Disabled cavity; enabled cavities have kp=6 and are handled above
      ! Note: kz=-12 is transformed into +12 in daten after reading ENDE.
      ktrack(i)=31
      goto 290
    end if
#include "include/beams21.f90"
#include "include/beamcoo.f90"
#include "include/beamr1.f90"
     &goto 42
#include "include/beamr2.f90"
#include "include/beamr3o.f90"
#include "include/beams22.f90"
#include "include/beam11.f90"
#include "include/beama1.f90"
#include "include/beamcoo.f90"
#include "include/beama2.f90"
#include "include/beam12.f90"
#include "include/beama3.f90"
#include "include/beam13.f90"
#include "include/beama4o.f90"
    else if(ibtyp.eq.1) then
#include "include/beam11.f90"
#include "include/beama1.f90"
#include "include/beamcoo.f90"
#include "include/beama2.f90"
#include "include/beama3.f90"
#include "include/beamwzf1.f90"
#include "include/beama4o.f90"
#include "include/beams23.f90"
#include "include/beam21.f90"
#include "include/beama1.f90"
#include "include/beamcoo.f90"
#include "include/beama2.f90"
#include "include/beam22.f90"
#include "include/beama3.f90"
#include "include/beam23.f90"
#include "include/beama4o.f90"
    else if(ibtyp.eq.1) then
#include "include/beam21.f90"
#include "include/beama1.f90"
#include "include/beamcoo.f90"
#include "include/beama2.f90"
#include "include/beama3.f90"
#include "include/beamwzf2.f90"
#include "include/beama4o.f90"
#include "include/beams24.f90"
    ! wire
    if(kzz.eq.15) then
      ktrack(i)=45
      goto 290
    endif
    !electron lens (HEL)
    if(kzz.eq.29) then
      ktrack(i)=63
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
        ktrack(i) =31
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
      r0=ek(ix)
      nmz=nmu(ix)
      if(abs(r0).le.pieni.or.nmz.eq.0) then
        if(abs(dki(ix,1)).le.pieni.and.abs(dki(ix,2)).le.pieni) then
          if ( dynk_isused(i) ) then
            write(lout,*) "ERROR: Element of type 11 (bez=",bez(ix),") is off in fort.2, but on in DYNK. Not implemented."
            call prror(-1)
          end if
          ktrack(i)=31
        else if(abs(dki(ix,1)).gt.pieni.and.abs(dki(ix,2)).le.pieni) then
          if(abs(dki(ix,3)).gt.pieni) then
            ktrack(i)=33
#include "include/stra11.f90"
          else
            ktrack(i)=35
#include "include/stra12.f90"
          end if
        else if(abs(dki(ix,1)).le.pieni.and.abs(dki(ix,2)).gt.pieni) then
          if(abs(dki(ix,3)).gt.pieni) then
            ktrack(i)=37
#include "include/stra13.f90"
          else
            ktrack(i)=39
#include "include/stra14.f90"
          end if
        end if
      else
        if(abs(dki(ix,1)).le.pieni.and.abs(dki(ix,2)).le.pieni) then
          ktrack(i)=32
        else if(abs(dki(ix,1)).gt.pieni.and.abs(dki(ix,2)).le.pieni) then
          if(abs(dki(ix,3)).gt.pieni) then
            ktrack(i)=34
#include "include/stra11.f90"
          else
            ktrack(i)=36
#include "include/stra12.f90"
          end if
        else if(abs(dki(ix,1)).le.pieni.and.abs(dki(ix,2)).gt.pieni) then
          if(abs(dki(ix,3)).gt.pieni) then
            ktrack(i)=38
#include "include/stra13.f90"
          else
            ktrack(i)=40
#include "include/stra14.f90"
          end if
        end if
      end if
      if(abs(r0).le.pieni.or.nmz.eq.0) goto 290
      if(mout2.eq.1) then
        benkcc=ed(ix)*benkc(irm(ix))
        r0a=one
        r000=r0*r00(irm(ix))
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
290 continue ! Label is still needed as it is referenced in some of the ca blocks

  do j=1,napx
    dpsv1(j)=(dpsv(j)*c1e3)/(one+dpsv(j))                            !hr01
  end do
  nwri=nwr(3)
  if(nwri.eq.0) nwri=numl+numlr+1
    ! A.Mereghetti, for the FLUKA Team
    ! last modified: 17-07-2013
    ! save original kicks
    ! always in main code
    if (dynk_enabled) call dynk_pretrack
    call time_timeStamp(time_afterPreTrack)

    if(idp.eq.0.or.ition.eq.0) then
      write(lout,*) ''
      write(lout,*) 'Calling thck4d subroutine'
      write(lout,*) ''
      call thck4d(nthinerr)
    else
      hsy(3)=(c1m3*hsy(3))*real(ition,fPrec)                                 !hr01

      do jj=1,nele
        if(kz(jj).eq.12) hsyc(jj)=(c1m3*hsyc(jj))*real(itionc(jj),fPrec)     !hr01
      end do

      if(abs(phas).ge.pieni) then
        write(lout,"(a)") "TRACKING> ERROR thck6dua no longer supported. Please use DYNK instead."
        call prror(-1)
      else
        write(lout,*) ''
        write(lout,*) 'Calling thck6d subroutine'
        write(lout,*) ''
        call thck6d(nthinerr)
      end if
    end if

  return
end subroutine trauthck

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

  use bdex, only : bdex_enable
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

  use mod_settings
  use mod_meta
  use mod_hions
  use postprocessing, only : writebin
  use crcoall
  use parpro
  use mod_common
  use mod_commonmn
  use mod_commons
  use mod_commont
  use mod_commond
  use elens
  use utils
  use wire
#ifdef CR
  use checkpoint_restart
#endif
  implicit none

  integer i,idz1,idz2,irrtr,ix,j,jb,jmel,jx,k,n,nmz,nthinerr,xory,nac,nfree,nramp1,nplato,nramp2,turnrep
  real(kind=fPrec) cccc,cikve,crkve,crkveuk,puxve,puxve1,puxve2,puzve1,puzve2,puzve,r0,xlvj,yv1j,   &
    yv2j,zlvj,acdipamp,qd,acphase, acdipamp2,acdipamp1,crabamp,crabfreq,kcrab,RTWO,NNORM,l,cur,dx,  &
    dy,tx,ty,embl,chi,xi,yi,dxi,dyi,rrelens,frrelens,xelens,yelens,onedp,fppsig,costh_temp,         &
    sinth_temp,pxf,pyf,r_temp,z_temp,sigf,q_temp !solenoid


  logical llost
  real(kind=fPrec) crkveb(npart),cikveb(npart),rho2b(npart),tkb(npart),r2b(npart),rb(npart),        &
    rkb(npart),xrb(npart),zrb(npart),xbb(npart),zbb(npart),crxb(npart),crzb(npart),cbxb(npart),     &
    cbzb(npart)

#ifdef FLUKA
  logical recompute_linear_matrices
#endif

  save

  nthinerr=0
  idz1=idz(1)
  idz2=idz(2)

  ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
  ! last modified: 24-11-2016
  ! initialise variables for back-tracking particles
  if (lbacktracking) call aperture_backTrackingInit

#ifdef FLUKA

!     A.Mereghetti, for the FLUKA Team
!     last modified: 14-06-2014
!     initialise napxto
!     inserted in main code by the 'fluka' compilation flag
  napxto = 0
#endif
  ! Determine which turns to print tracking report on
  if(numl > 1000) then
    turnrep = nint(numl/1000.0)
  else
    turnrep = 1
  end if

#ifdef CR
  if (restart) then
    call crstart
    write(93,*) 'THCK4D ','SIXTRACR restart numlcr',numlcr,'numl',numl
    ! and now reset numl to do only numlmax turns
  end if
  nnuml=min((numlcr/numlmax+1)*numlmax,numl)
  write (93,*) 'numlmax=',numlmax,' DO ',numlcr,nnuml
  ! and reset [n]numxv unless particle is lost
  ! TRYing Eric (and removing postpr fixes).
  if (nnuml.ne.numl) then
    do j=1,napx
      if (numxv(j).eq.numl) numxv(j)=nnuml
      if (nnumxv(j).eq.numl) nnumxv(j)=nnuml
    end do
  end if
  do 490 n=numlcr,nnuml
#endif
#ifndef CR
  do 490 n=1,numl
#endif
    if(st_quiet < 3) then
      if(mod(n,turnrep) == 0) write(lout,"(a,i8,a,i8)") "TRACKING> Thick 4D turn ",n," of ",numl
    end if
    meta_nPartTurn = meta_nPartTurn + napx
#ifdef BOINC
!   call boinc_sixtrack_progress(n,numl)
    call boinc_fraction_done(dble(n)/dble(numl))
    continue
!   call graphic_progress(n,numl)
#endif
  numx=n-1

#ifndef FLUKA
    if(mod(numx,nwri).eq.0) call writebin(nthinerr)
    if(nthinerr.ne.0) return
#endif

#ifdef CR
    !  does not call CRPOINT if restart=.true.
    !  (and note that writebin does nothing if restart=.true.
    if(mod(numx,numlcp).eq.0) call callcrp()
    restart=.false.
#endif

!       A.Mereghetti, for the FLUKA Team
!       last modified: 03-09-2014
!       apply dynamic kicks
!       always in main code
    if ( dynk_enabled ) then
      call dynk_apply(n)
    end if

    call dump_linesFirst(n)

    do 480 i=1,iu
      if(ktrack(i).eq.1) then
        ix=ic(i)
      else
        ix=ic(i)-nblo

        if (ldumpfront) then
          write (lout,*) "DUMP/FRONT not yet supported on thick elements "// &
                         "due to lack of test cases. Please contact developers!"
          call prror(-1)
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
              call kernel_fluka_exit( n, i, ix )
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
            write(lout,*) '[Fluka] Skipping lattice element at ',i
            write(fluka_log_unit,*) '# Skipping lattice element at ', i
          end if
          goto 480
        end if
      end if
#endif

            if (bdex_enable) then
               !TODO - if you have a test case, please contact developers!
               write(lout,*) "BDEX> BDEX only available for thin6d"
               call prror(-1)
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
            xv1(j)=bl1v(1,1,j,ix)*puxve+bl1v(2,1,j,ix)*puzve+((real(idz1,fPrec)*bl1v(5,1,j,ix))*dpsv(j))*c1e3 !hr01
            yv1(j)=bl1v(3,1,j,ix)*puxve+bl1v(4,1,j,ix)*puzve+((real(idz1,fPrec)*bl1v(6,1,j,ix))*dpsv(j))*c1e3 !hr01
            puxve=xv2(j)
            puzve=yv2(j)
            xv2(j)=bl1v(1,2,j,ix)*puxve+bl1v(2,2,j,ix)*puzve+((real(idz2,fPrec)*bl1v(5,2,j,ix))*dpsv(j))*c1e3 !hr01
            yv2(j)=bl1v(3,2,j,ix)*puxve+bl1v(4,2,j,ix)*puzve+((real(idz2,fPrec)*bl1v(6,2,j,ix))*dpsv(j))*c1e3 !hr01
          end do
        end if
        goto 480
      case (2)
        goto 480
      case (3)
        irrtr=imtr(ix)
        do j=1,napx
          temptr(1)=xv1(j)
          temptr(2)=yv1(j)
          temptr(3)=xv2(j)
          temptr(4)=yv2(j)

          xv1(j)  = cotr(irrtr,1)
          yv1(j)  = cotr(irrtr,2)
          xv2(j)  = cotr(irrtr,3)
          yv2(j)  = cotr(irrtr,4)

          do kxxa=1,6
            xv1(j)   =  xv1(j)+temptr(kxxa)*rrtr(irrtr,1,kxxa)
            yv1(j)   =  yv1(j)+temptr(kxxa)*rrtr(irrtr,2,kxxa)
            xv2(j)   =  xv2(j)+temptr(kxxa)*rrtr(irrtr,3,kxxa)
            yv2(j)   =  yv2(j)+temptr(kxxa)*rrtr(irrtr,4,kxxa)
          enddo
        enddo
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
        do 690 j=1,napx
#include "include/beamco.f90"
#include "include/beamr1.f90"
     &goto 690
#include "include/beamr2.f90"
#include "include/beamr3.f90"
690     continue
        goto 470
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
        goto 470
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
      case (63) ! Elens
        do j=1,napx
#include "include/kickelens.f90"
        end do
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
        call envarsv(dpsv,moidpsv,rvv,ekv)
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

490 continue

  return

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
  use bdex, only : bdex_enable
  use dynk, only : dynk_enabled, dynk_apply
  use dump, only : dump_linesFirst, dump_lines, ldumpfront
  use collimation, only: do_coll, part_abs_turn
  use aperture

#ifdef FLUKA
!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 17-07-2013
!     import mod_fluka
!     inserted in main code by the 'fluka' compilation flag
  use mod_fluka
#endif

  use mod_meta
  use mod_settings
  use mod_hions
  use postprocessing, only : writebin
  use crcoall
  use parpro
  use mod_common
  use mod_commonmn
  use mod_commons
  use mod_commont
  use mod_commond
  use aperture
  use elens
  use utils
  use wire
#ifdef CR
  use checkpoint_restart
#endif

  implicit none

  integer i,idz1,idz2,irrtr,ix,j,jb,jmel,jx,k,n,nmz,nthinerr,xory,nac,nfree,nramp1,nplato,nramp2,turnrep
  real(kind=fPrec) cccc,cikve,crkve,crkveuk,puxve1,puxve2,puzve1,puzve2,r0,xlvj,yv1j,yv2j,zlvj,     &
    acdipamp,qd,acphase,acdipamp2,acdipamp1,crabamp,crabfreq,kcrab,RTWO,NNORM,l,cur,dx,dy,tx,ty,    &
    embl,chi,xi,yi,dxi,dyi,rrelens,frrelens,xelens,yelens,onedp,fppsig,costh_temp,sinth_temp,pxf,   &
    pyf,r_temp,z_temp,sigf,q_temp !solenoid
  logical llost
  real(kind=fPrec) crkveb(npart),cikveb(npart),rho2b(npart),tkb(npart),r2b(npart),rb(npart),        &
    rkb(npart),xrb(npart),zrb(npart),xbb(npart),zbb(npart),crxb(npart),crzb(npart),cbxb(npart),     &
    cbzb(npart)

#ifdef FLUKA
  logical recompute_linear_matrices
#endif

  save
#ifdef DEBUG
!-----------------------------------------------------------------------
!===================================================================
! Eric beginthck6dstart
!===================================================================
#endif
  nthinerr=0
  idz1=idz(1)
  idz2=idz(2)

#ifdef FLUKA
!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 17-07-2013
!     force re-computation of transport matrices of linear elements
!     inserted in main code by the 'fluka' compilation flag
  recompute_linear_matrices = .false.
#endif

  ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
  ! last modified: 24-11-2016
  ! initialise variables for back-tracking particles
  if (lbacktracking) call aperture_backTrackingInit

#ifdef FLUKA
!     A.Mereghetti, for the FLUKA Team
!     last modified: 14-06-2014
!     initialise napxto
!     inserted in main code by the 'fluka' compilation flag
  napxto = 0
#endif

  ! Determine which turns to print tracking report on
  if(numl > 1000) then
    turnrep = nint(numl/1000.0)
  else
    turnrep = 1
  end if

! Now the outer loop over turns
#ifdef CR
  if (restart) then
    call crstart
    write(93,*) 'THCK6D ','SIXTRACR restart numlcr',numlcr,'numl',numl
! and now reset numl to do only numlmax turns
  end if
  nnuml=min((numlcr/numlmax+1)*numlmax,numl)
  write (93,*) 'numlmax=',numlmax,' DO ',numlcr,nnuml
! and reset [n]numxv unless particle is lost
! TRYing Eric (and removing postpr fixes).
  if (nnuml.ne.numl) then
    do j=1,napx
      if (numxv(j).eq.numl) numxv(j)=nnuml
      if (nnumxv(j).eq.numl) nnumxv(j)=nnuml
    end do
  end if
  do 510 n=numlcr,nnuml
#endif
#ifndef CR
  do 510 n=1,numl
#endif
    if(st_quiet < 3) then
      if(mod(n,turnrep) == 0) write(lout,"(a,i8,a,i8)") "TRACKING> Thick 6D turn ",n," of ",numl
    end if
    meta_nPartTurn = meta_nPartTurn + napx
! To do a dump and abend
#ifdef BOINC
!   call boinc_sixtrack_progress(n,numl)
    call boinc_fraction_done(dble(n)/dble(numl))
    continue
!   call graphic_progress(n,numl)
#endif
    numx=n-1

#ifndef FLUKA
    if(mod(numx,nwri).eq.0) call writebin(nthinerr)
    if(nthinerr.ne.0) return
#endif

#ifdef CR
!  does not call CRPOINT if restart=.true.
!  (and note that writebin does nothing if restart=.true.
    if(mod(numx,numlcp).eq.0) call callcrp()
    restart=.false.
#endif

!       A.Mereghetti, for the FLUKA Team
!       last modified: 03-09-2014
!       apply dynamic kicks
!       always in main code
    if ( dynk_enabled ) then
      call dynk_apply(n)
    end if
    call dump_linesFirst(n)

#ifdef DEBUG
! Now comes the loop over elements do 500/501
    do 501 i=1,iu
#else
    do 500 i=1,iu
#endif
#ifdef DEBUG
!===================================================================
!===================================================================
! Eric endthck6dstart
! Nothing should be changed in the rest of this loop
!===================================================================
!===================================================================
#endif
      if(ktrack(i).eq.1) then
        ix=ic(i)
      else
        ix=ic(i)-nblo
      end if

      if (ldumpfront) then
        write (lout,*) "DUMP/FRONT not yet supported on thick elements "// &
                       "due to lack of test cases. Please contact developers!"
        call prror(-1)
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
              call kernel_fluka_exit( n, i, ix )
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
            write(lout,*) '[Fluka] Skipping lattice element at ',i
            write(fluka_log_unit,*) '# Skipping lattice element at ', i
          end if
          goto 500
        end if
      end if
#endif

#ifdef DEBUG
!     if (i.ge.673) then
!     call warr('xv12,i,ktrack ',xv1(2),i,ktrack(i),0,0)
!     endif
!     if (i.eq.676) stop
#endif

            if (bdex_enable) then
               !TODO - if you have a test case, please contact developers!
               write(lout,*) "BDEX> BDEX only available for thin6d"
               call prror(-1)
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
      case (2)
        do j=1,napx
          ejf0v(j)=ejfv(j)
          if(abs(dppoff).gt.pieni) sigmv(j)=sigmv(j)-sigmoff(i)
          if(kz(ix).eq.12) then
            ejv(j)=ejv(j)+(ed(ix)*sin_mb(hsyc(ix)*sigmv(j)+phasc(ix)))*nzz(j)
          else
            ejv(j)=ejv(j)+(hsy(1)*sin_mb(hsy(3)*sigmv(j)))*nzz(j)
          end if
          ejfv(j)=sqrt(ejv(j)**2-nucm(j)**2)                             !hr01
          rvv(j)=(ejv(j)*e0f)/(e0*ejfv(j))
          dpsv(j)=(ejfv(j)-e0f)/e0f
          oidpsv(j)=one/(one+dpsv(j))
          moidpsv(j)=mtc(j)/(one+dpsv(j))
          omoidpsv(j)=c1e3*((one-mtc(j))*oidpsv(j))
          dpsv1(j)=(dpsv(j)*c1e3)*oidpsv(j)                          !hr01
          yv1(j)=(ejf0v(j)/ejfv(j))*yv1(j)                         !hr01
          yv2(j)=(ejf0v(j)/ejfv(j))*yv2(j)                         !hr01
        end do
        if(n.eq.1) write(98,'(1p,6(2x,e25.18))') (xv1(j),yv1(j),xv2(j),yv2(j),sigmv(j),dpsv(j),j=1,napx)
#ifdef CR
        ! write(93,*) 'ERIC loop at 40 calling synuthck!!!'
        ! endfile (93,iostat=ierro)
        ! backspace (93,iostat=ierro)
#endif
        call synuthck
        goto 490
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
          dpsv1(j)=(dpsv(j)*c1e3)*oidpsv(j)


          ! We have to go back to angles after we updated the energy.
          yv1(j) = yv1(j)*moidpsv(j)
          yv2(j) = yv2(j)*moidpsv(j)

        enddo
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
        do 690 j=1,napx
#include "include/beamco.f90"
#include "include/beamr1.f90"
     &goto 690
#include "include/beamr2.f90"
#include "include/beamr3.f90"
690     continue
        goto 490
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
        goto 490
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
      case (63) ! Elens
        do j=1,napx
#include "include/kickelens.f90"
        end do
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
        call envarsv(dpsv,moidpsv,rvv,ekv)
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

#ifdef DEBUG
500 continue
    ! if (n.ge.990) then
    !   write(99,*) 'after element i, ktrack ',i,ktrack(i), xv1(1),xv2(1),yv1(1),yv2(1),&
    !     sigmv(1),ejv(1),ejfv(1),rvv(1),dpsv(1),oidpsv(1),dpsv1(1)
    !   endfile (99,iostat=ierro)
    !   backspace (99,iostat=ierro)
    ! end if
501 continue
#endif
#ifndef DEBUG
500 continue
#endif
! End of loop over elements

!===================================================================
!===================================================================
! Eric beginthck6dend
!===================================================================
!===================================================================

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

510 continue
! end loop over turns

!===================================================================
!===================================================================
! Eric endthck6dend
!===================================================================
!===================================================================
  return

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
  use mod_commonmn
  use mod_commons
  use mod_commont
  use mod_commond
  implicit none
  integer ih1,ih2,j,kz1,l
  real(kind=fPrec) fokm
  save
!---------------------------------------  SUBROUTINE 'ENVARS' IN-LINE
#ifdef CR
#ifdef DEBUG
!       write(93,*) 'ERIC synuthck called!!!'
!       write(93,*) 'ERIC il= ',il
!       endfile (93,iostat=ierro)
!       backspace (93,iostat=ierro)
#endif
  sythckcr=.true.
#endif
  do 10 j=1,napx
    dpd(j)=one+dpsv(j)
    dpsq(j)=sqrt(dpd(j))
10 continue
  do 160 l=1,il
    if(abs(el(l)).le.pieni) goto 160
    kz1=kz(l)+1
!       goto(20,40,80,60,40,60,100,100,140),kz1
!       goto 160
!Eric
!-----------------------------------------------------------------------
!  DRIFTLENGTH
!-----------------------------------------------------------------------
    if (kz1.eq.1) then
      goto 20
!-----------------------------------------------------------------------
!  RECTANGULAR MAGNET
!  HORIZONTAL
!-----------------------------------------------------------------------
    elseif (kz1.eq.2.or.kz1.eq.5) then
40     fokm=el(l)*ed(l)
      if(abs(fokm).le.pieni) goto 20
      if(kz1.eq.2) then
        ih1=1
        ih2=2
      else
!  RECTANGULAR MAGNET VERTICAL
        ih1=2
        ih2=1
      endif
      do 50 j=1,napx
        fok(j)=fokm/dpsq(j)
        rho(j)=(one/ed(l))*dpsq(j)
        fok1(j)=(tan_mb(fok(j)*half))/rho(j)
        si(j)=sin_mb(fok(j))
        co(j)=cos_mb(fok(j))
        al(2,ih1,j,l)=rho(j)*si(j)
        al(5,ih1,j,l)=((-one*dpsv(j))*((rho(j)*(one-co(j)))/dpsq(j)))*c1e3
        al(6,ih1,j,l)=((-one*dpsv(j))*((two*tan_mb(fok(j)*half))/dpsq(j)))*c1e3         !hr01
        sm1(j)=cos_mb(fok(j))
        sm2(j)=sin_mb(fok(j))*rho(j)
        sm3(j)=-sin_mb(fok(j))/rho(j)
        sm12(j)=el(l)-sm1(j)*sm2(j)
        sm23(j)=sm2(j)*sm3(j)
        as3(j)=(-one*rvv(j))*(((dpsv(j)*rho(j))/(two*dpsq(j)))*sm23(j)-(rho(j)*dpsq(j))*(one-sm1(j)))
        as4(j)=((-one*rvv(j))*sm23(j))/c2e3                          !hr01
        as6(j)=((-one*rvv(j))*(el(l)+sm1(j)*sm2(j)))/c4e3            !hr01
        as(1,ih1,j,l)=(el(l)*(one-rvv(j))-rvv(j)*((dpsv(j)**2/(four*dpd(j)))*sm12(j)+dpsv(j)*(el(l)-sm2(j))))*c1e3
        as(2,ih1,j,l)=(-one*rvv(j))*((dpsv(j)/((two*rho(j))*dpsq(j)))* sm12(j)-(sm2(j)*dpsq(j))/rho(j))+fok1(j)*as3(j)
        as(3,ih1,j,l)=as3(j)
        as(4,ih1,j,l)=as4(j)+(two*as6(j))*fok1(j)                    !hr01
        as(5,ih1,j,l)=((-one*rvv(j))*sm12(j))/(c4e3*rho(j)**2)+as6(j)*fok1(j)**2+fok1(j)*as4(j)
        as(6,ih1,j,l)=as6(j)
!--VERTIKAL
        g(j)=tan_mb(fok(j)*half)/rho(j)
        gl(j)=el(l)*g(j)
        al(1,ih2,j,l)=one-gl(j)
        al(3,ih2,j,l)=(-one*g(j))*(two-gl(j))                        !hr01
        al(4,ih2,j,l)=al(1,ih2,j,l)
        as6(j)=((-one*rvv(j))*al(2,ih2,j,l))/c2e3                    !hr01
        as(4,ih2,j,l)=((-one*two)*as6(j))*fok1(j)                    !hr01
        as(5,ih2,j,l)=(as6(j)*fok1(j))*fok1(j)                       !hr01
        as(6,ih2,j,l)=as6(j)
50     continue
      goto 160
    elseif (kz1.eq.4.or.kz1.eq.6) then
!-----------------------------------------------------------------------
!  SEKTORMAGNET
!  HORIZONTAL
!-----------------------------------------------------------------------
60     fokm=el(l)*ed(l)
      if(abs(fokm).le.pieni) goto 20
      if(kz1.eq.4) then
        ih1=1
        ih2=2
      else
!  SECTOR MAGNET VERTICAL
        ih1=2
        ih2=1
      endif
      do 70 j=1,napx
        fok(j)=fokm/dpsq(j)
        rho(j)=(one/ed(l))*dpsq(j)
        si(j)=sin_mb(fok(j))
        co(j)=cos_mb(fok(j))
        rhoc(j)=(rho(j)*(one-co(j)))/dpsq(j)                         !hr01
        siq(j)=si(j)/dpsq(j)
        al(1,ih1,j,l)=co(j)
        al(2,ih1,j,l)=rho(j)*si(j)
        al(3,ih1,j,l)=(-one*si(j))/rho(j)                            !hr01
        al(4,ih1,j,l)=co(j)
        al(5,ih1,j,l)=((-one*dpsv(j))*rhoc(j))*c1e3                  !hr01
        al(6,ih1,j,l)=((-one*dpsv(j))*siq(j))*c1e3                   !hr01
        sm12(j)=el(l)-al(1,ih1,j,l)*al(2,ih1,j,l)
        sm23(j)=al(2,ih1,j,l)*al(3,ih1,j,l)
        as(1,ih1,j,l)=(el(l)*(one-rvv(j))-rvv(j)*((dpsv(j)**2/(four*dpd(j)))*sm12(j)+dpsv(j)*(el(l)-al(2,ih1,j,l))))*c1e3
        as(2,ih1,j,l)=(-one*rvv(j))*((dpsv(j)/(two*rho(j)*dpsq(j)))*sm12(j)-dpd(j)*siq(j))
        as(3,ih1,j,l)=(-one*rvv(j))*(((dpsv(j)*rho(j))/(two*dpsq(j)))*sm23(j)-dpd(j)*rhoc(j))
        as(4,ih1,j,l)=((-one*rvv(j))*sm23(j))/c2e3                   !hr01
        as(5,ih1,j,l)=((-one*rvv(j))*sm12(j))/((c4e3*rho(j))*rho(j)) !hr01
        as(6,ih1,j,l)=((-one*rvv(j))*(el(l)+al(1,ih1,j,l)*al(2,ih1,j,l)))/c4e3
!--VERTIKAL
        as(6,ih2,j,l)=((-one*rvv(j))*al(2,ih2,j,l))/c2e3             !hr01
70     continue
      goto 160
    elseif (kz1.eq.3) then
!-----------------------------------------------------------------------
!  QUADRUPOLE
!  FOCUSSING
!-----------------------------------------------------------------------
80   do 90 j=1,napx
        fok(j)=ekv(j,l)*oidpsv(j)
        aek(j)=abs(fok(j))
        hi(j)=sqrt(aek(j))
        fi(j)=el(l)*hi(j)
        if(fok(j).le.zero) then
          al(1,1,j,l)=cos_mb(fi(j))
          hi1(j)=sin_mb(fi(j))
          if(abs(hi(j)).le.pieni) then
            al(2,1,j,l)=el(l)
          else
            al(2,1,j,l)=hi1(j)/hi(j)
          endif
          al(3,1,j,l)=-hi1(j)*hi(j)
          al(4,1,j,l)=al(1,1,j,l)
          as(1,1,j,l)=el(l)*(one-rvv(j))*c1e3
          as(4,1,j,l)=(((-one*rvv(j))*al(2,1,j,l))*al(3,1,j,l))/c2e3
          as(5,1,j,l)=(((-one*rvv(j))*(el(l)-al(1,1,j,l)*al(2,1,j,l)))*aek(j))/c4e3
          as(6,1,j,l)=((-one*rvv(j))*(el(l)+al(1,1,j,l)*al(2,1,j,l)))/c4e3
!--DEFOCUSSING
          hp(j)=exp_mb(fi(j))
          hm(j)=one/hp(j)
          hc(j)=(hp(j)+hm(j))*half
          hs(j)=(hp(j)-hm(j))*half
          al(1,2,j,l)=hc(j)
          if(abs(hi(j)).le.pieni) then
            al(2,2,j,l)=el(l)
          else
            al(2,2,j,l)=hs(j)/hi(j)
          endif
          al(3,2,j,l)=hs(j)*hi(j)
          al(4,2,j,l)=hc(j)
          as(4,2,j,l)=(((-one*rvv(j))*al(2,2,j,l))*al(3,2,j,l))/c2e3
          as(5,2,j,l)=((rvv(j)*(el(l)-al(1,2,j,l)*al(2,2,j,l)))*aek(j))/c4e3
          as(6,2,j,l)=((-one*rvv(j))*(el(l)+al(1,2,j,l)*al(2,2,j,l)))/c4e3
        else
          al(1,2,j,l)=cos_mb(fi(j))
          hi1(j)=sin_mb(fi(j))
          if(abs(hi(j)).le.pieni) then
            al(2,2,j,l)=el(l)
          else
            al(2,2,j,l)=hi1(j)/hi(j)
          endif
          al(3,2,j,l)=(-one*hi1(j))*hi(j)                            !hr01
          al(4,2,j,l)=al(1,2,j,l)
          as(1,2,j,l)=(el(l)*(one-rvv(j)))*c1e3                      !hr01
          as(4,2,j,l)=(((-one*rvv(j))*al(2,2,j,l))*al(3,2,j,l))/c2e3 !hr01
          as(5,2,j,l)=(((-one*rvv(j))*(el(l)-al(1,2,j,l)*al(2,2,j,l)))*aek(j))/c4e3
          as(6,2,j,l)=((-one*rvv(j))*(el(l)+al(1,2,j,l)*al(2,2,j,l)))/c4e3 !hr01
!--DEFOCUSSING
          hp(j)=exp_mb(fi(j))
          hm(j)=one/hp(j)
          hc(j)=(hp(j)+hm(j))*half
          hs(j)=(hp(j)-hm(j))*half
          al(1,1,j,l)=hc(j)
          if(abs(hi(j)).le.pieni) then
            al(2,1,j,l)=el(l)
          else
            al(2,1,j,l)=hs(j)/hi(j)
          endif
          al(3,1,j,l)=hs(j)*hi(j)
          al(4,1,j,l)=hc(j)
          as(4,1,j,l)=(((-one*rvv(j))*al(2,1,j,l))*al(3,1,j,l))/c2e3 !hr01
          as(5,1,j,l)=((rvv(j)*(el(l)-al(1,1,j,l)*al(2,1,j,l)))*aek(j))/c4e3
          as(6,1,j,l)=((-one*rvv(j))*(el(l)+al(1,1,j,l)*al(2,1,j,l)))/c4e3 !hr01
        endif
90     continue
      goto 160
    elseif (kz1.eq.7.or.kz1.eq.8) then
!-----------------------------------------------------------------------
!  COMBINED FUNCTION MAGNET HORIZONTAL
!  FOCUSSING
!-----------------------------------------------------------------------
100     if(kz1.eq.7) then
        do 110 j=1,napx
          fokqv(j)=ekv(j,l)
110       continue
        ih1=1
        ih2=2
      else
!  COMBINED FUNCTION MAGNET VERTICAL
        do 120 j=1,napx
          fokqv(j)=-ekv(j,l)
120       continue
        ih1=2
        ih2=1
      endif
      do 130 j=1,napx
        wf(j)=ed(l)/dpsq(j)
        fok(j)=fokqv(j)/dpd(j)-wf(j)**2                              !hr01
        afok(j)=abs(fok(j))
        hi(j)=sqrt(afok(j))
        fi(j)=hi(j)*el(l)
        if(afok(j).le.pieni) then
          as(6,1,j,l)=((-one*rvv(j))*el(l))/c2e3                     !hr01
          as(6,2,j,l)=as(6,1,j,l)
          as(1,1,j,l)=(el(l)*(one-rvv(j)))*c1e3                      !hr01
        endif
        if(fok(j).lt.(-one*pieni)) then                              !hr06
          si(j)=sin_mb(fi(j))
          co(j)=cos_mb(fi(j))
          wfa(j)=((wf(j)/afok(j))*(one-co(j)))/dpsq(j)               !hr01
          wfhi(j)=((wf(j)/hi(j))*si(j))/dpsq(j)                      !hr01
          al(1,ih1,j,l)=co(j)
          al(2,ih1,j,l)=si(j)/hi(j)
          al(3,ih1,j,l)=(-one*si(j))*hi(j)                           !hr01
          al(4,ih1,j,l)=co(j)
          al(5,ih1,j,l)=((-one*wfa(j))*dpsv(j))*c1e3                 !hr01
          al(6,ih1,j,l)=((-one*wfhi(j))*dpsv(j))*c1e3                !hr01
          sm12(j)=el(l)-al(1,ih1,j,l)*al(2,ih1,j,l)
          sm23(j)=al(2,ih1,j,l)*al(3,ih1,j,l)
          as(1,ih1,j,l)=(el(l)*(one-rvv(j))-((rvv(j)*((dpsv(j)**2/(four*dpd(j)))&
            *sm12(j)+dpsv(j)*(el(l)-al(2,ih1,j,l))))/afok(j))*wf(j)**2)*c1e3
          as(2,ih1,j,l)=(-one*rvv(j))*(((dpsv(j)*wf(j))/(two*dpsq(j)))*sm12(j)-dpd(j)*wfhi(j))
          as(3,ih1,j,l)=(-one*rvv(j))*(((((dpsv(j)*half)/afok(j))/dpd(j))*ed(l))*sm23(j)-dpd(j)*wfa(j))
          as(4,ih1,j,l)=((-one*rvv(j))*sm23(j))/c2e3
          as(5,ih1,j,l)=(((-one*rvv(j))*sm12(j))*afok(j))/c4e3
          as(6,ih1,j,l)=((-one*rvv(j))*(el(l)+al(1,ih1,j,l)*al(2,ih1,j,l)))/c4e3
          aek(j)=abs(ekv(j,l)/dpd(j))
          hi(j)=sqrt(aek(j))
          fi(j)=hi(j)*el(l)
          hp(j)=exp_mb(fi(j))
          hm(j)=one/hp(j)
          hc(j)=(hp(j)+hm(j))*half
          hs(j)=(hp(j)-hm(j))*half
          al(1,ih2,j,l)=hc(j)
          if(abs(hi(j)).gt.pieni) al(2,ih2,j,l)=hs(j)/hi(j)
          al(3,ih2,j,l)=hs(j)*hi(j)
          al(4,ih2,j,l)=hc(j)
          as(4,ih2,j,l)=(((-one*rvv(j))*al(2,ih2,j,l))*al(3,ih2,j,l))/c2e3
          as(5,ih2,j,l)=((rvv(j)*(el(l)-al(1,ih2,j,l)*al(2,ih2,j,l)))*aek(j))/c4e3
          as(6,ih2,j,l)=((-one*rvv(j))*(el(l)+al(1,ih2,j,l)*al(2,ih2,j,l)))/c4e3
        endif
!--DEFOCUSSING
        if(fok(j).gt.pieni) then
          hp(j)=exp_mb(fi(j))
          hm(j)=one/hp(j)
          hc(j)=(hp(j)+hm(j))*half
          hs(j)=(hp(j)-hm(j))*half
          al(1,ih1,j,l)=hc(j)
          al(2,ih1,j,l)=hs(j)/hi(j)
          al(3,ih1,j,l)=hs(j)*hi(j)
          al(4,ih1,j,l)=hc(j)
          wfa(j)=((wf(j)/afok(j))*(one-hc(j)))/dpsq(j)               !hr01
          wfhi(j)=((wf(j)/hi(j))*hs(j))/dpsq(j)                      !hr01
          al(5,ih1,j,l)= (wfa(j)*dpsv(j))*c1e3                       !hr01
          al(6,ih1,j,l)=((-one*wfhi(j))*dpsv(j))*c1e3                !hr01
          sm12(j)=el(l)-al(1,ih1,j,l)*al(2,ih1,j,l)
          sm23(j)=al(2,ih1,j,l)*al(3,ih1,j,l)
          as(1,ih1,j,l)=(((rvv(j)*((dpsv(j)**2/(four*dpd(j)))*sm12(j)&
            +dpsv(j)*(el(l)-al(2,ih1,j,l))))/afok(j))*wf(j)**2+el(l)*(one-rvv(j)))*c1e3
          as(2,ih1,j,l)=(-one*rvv(j))*(((dpsv(j)*wf(j))/(two*dpsq(j)))*sm12(j)-dpd(j)*wfhi(j))
          as(3,ih1,j,l)=rvv(j)*(((((dpsv(j)*half)/afok(j))/dpd(j))* ed(l))*sm23(j)-dpd(j)*wfa(j))
          as(4,ih1,j,l)=((-one*rvv(j))*sm23(j))/c2e3                 !hr01
          as(5,ih1,j,l)=((rvv(j)*sm12(j))*afok(j))/c4e3              !hr01
          as(6,ih1,j,l)=((-one*rvv(j))*(el(l)+al(1,ih1,j,l)*al(2,ih1,j,l)))/c4e3
          aek(j)=abs(ekv(j,l)/dpd(j))
          hi(j)=sqrt(aek(j))
          fi(j)=hi(j)*el(l)
          si(j)=sin_mb(fi(j))
          co(j)=cos_mb(fi(j))
          al(1,ih2,j,l)=co(j)
          al(2,ih2,j,l)=si(j)/hi(j)
          al(3,ih2,j,l)=(-one*si(j))*hi(j)                           !hr01
          al(4,ih2,j,l)=co(j)
          as(4,ih2,j,l)=(((-one*rvv(j))*al(2,ih2,j,l))*al(3,ih2,j,l))/c2e3 !hr01
          as(5,ih2,j,l)=(((-one*rvv(j))*(el(l)-al(1,ih2,j,l)*al(2,ih2,j,l)))*aek(j))/c4e3
          as(6,ih2,j,l)=((-one*rvv(j))*(el(l)+al(1,ih2,j,l)*al(2,ih2,j,l)))/c4e3
        endif
130     continue
      goto 160
    elseif (kz1.eq.9) then
!-----------------------------------------------------------------------
!  EDGE FOCUSSING
!-----------------------------------------------------------------------
140     do 150 j=1,napx
        rhoi(j)=ed(l)/dpsq(j)
        fok(j)=rhoi(j)*tan_mb((el(l)*rhoi(j))*half)                  !hr01
        al(3,1,j,l)=fok(j)
        al(3,2,j,l)=-fok(j)
150     continue
      goto 160
    else
!Eric
! Is really an error but old code went to 160
      goto 160
    endif
!-----------------------------------------------------------------------
!  DRIFTLENGTH
!-----------------------------------------------------------------------
20   do 30 j=1,napx
      as(6,1,j,l)=((-one*rvv(j))*el(l))/c2e3                         !hr01
      as(6,2,j,l)=as(6,1,j,l)
      as(1,1,j,l)=(el(l)*(one-rvv(j)))*c1e3                          !hr01
30   continue
160 continue
!---------------------------------------  END OF 'ENVARS' (2)
  return
end subroutine synuthck
