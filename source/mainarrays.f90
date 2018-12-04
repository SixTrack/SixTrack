! ================================================================================================ !
! Routines for scaling array sizes defined in parpro
! by telling various other module to rescale their arrays;
! It is split off from parpro as a separate module
! in order to avoid circular dependencies.
! ================================================================================================ !

! Allocate arrays scaling with the main memory parameters nele, npart, etc.
subroutine allocate_arrays

  use parpro
  use crcoall

  use mod_common,         only : mod_common_expand_arrays
  use mod_commont,        only : mod_commont_expand_arrays
  use mod_commonmn,       only : mod_commonmn_expand_arrays
  use mod_commond2,       only : mod_commond2_expand_arrays
  use aperture,           only : aperture_expand_arrays
  use elens,              only : elens_allocate_arrays
  use dump,               only : dump_expand_arrays
  use scatter,            only : scatter_expand_arrays
  use bdex,               only : bdex_allocate_arrays
  use dynk,               only : dynk_allocate_arrays
  use wire,               only : wire_expand_arrays
  use mod_hions,          only : hions_allocate_arrays
#ifdef CR
  use checkpoint_restart, only : cr_expand_arrays
#endif
#ifdef FLUKA
  use mod_fluka,          only : fluka_mod_expand_arrays
#endif
  use collimation,        only : collimation_allocate_arrays
  implicit none

  nele  = nele_initial
  nblo  = nblo_initial
  nblz  = nblz_initial
  npart = npart_initial

  call mod_common_expand_arrays(nele,nblo,nblz,npart)
  call mod_commont_expand_arrays(nblz,npart)
  call mod_commonmn_expand_arrays(nblz,npart)
  call mod_commond2_expand_arrays(nele)

  call dump_expand_arrays(nele,nblz)
  call wire_expand_arrays(nele,nblz)
  call scatter_expand_arrays(nele)
  call aperture_expand_arrays(nele,npart)

  call elens_allocate_arrays
  call bdex_allocate_arrays
  call dynk_allocate_arrays
  call hions_allocate_arrays
#ifdef CR
  call cr_expand_arrays(npart)
#endif
  call collimation_allocate_arrays

end subroutine allocate_arrays

! Change the allocation of the arrays scaling with the main memry parameter nele, npart, etc.
subroutine expand_arrays(nele_new, npart_new, nblz_new, nblo_new)

  use mod_alloc, only : alloc_log

  use parpro
  use crcoall

  use mod_common,         only : mod_common_expand_arrays
  use mod_commont,        only : mod_commont_expand_arrays
  use mod_commonmn,       only : mod_commonmn_expand_arrays
  use mod_commond2,       only : mod_commond2_expand_arrays
  use aperture,           only : aperture_expand_arrays
  use elens,              only : elens_expand_arrays
  use dump,               only : dump_expand_arrays
  use scatter,            only : scatter_expand_arrays
  use bdex,               only : bdex_expand_arrays
  use dynk,               only : dynk_expand_arrays
  use wire,               only : wire_expand_arrays
  use mod_hions,          only : hions_expand_arrays
#ifdef CR
  use checkpoint_restart, only : cr_expand_arrays
#endif
#ifdef FLUKA
  use mod_fluka,          only : fluka_mod_expand_arrays
#endif
  use collimation,        only : collimation_expand_arrays
  implicit none

  integer, intent(in) :: nele_new
  integer, intent(in) :: npart_new
  integer, intent(in) :: nblz_new
  integer, intent(in) :: nblo_new

  write(alloc_log,"(a,4(1x,i0))") "ALLOC> Expanding (nele,npart,nblz,nblo):",nele_new,npart_new,nblz_new,nblo_new

  !Call sub-subroutines to actually expand
  call mod_common_expand_arrays(nele_new,nblo_new,nblz_new,npart_new)
  call mod_commont_expand_arrays(nblz_new,npart_new)
  call mod_commonmn_expand_arrays(nblz_new,npart_new)
  call mod_commond2_expand_arrays(nele_new)

  call dump_expand_arrays(nele_new,nblz_new)
  call wire_expand_arrays(nele_new,nblz_new)
  call scatter_expand_arrays(nele_new)
  call aperture_expand_arrays(nele_new, npart_new)

  call elens_expand_arrays(nele_new)
  call bdex_expand_arrays(nele_new)
  call dynk_expand_arrays(nele_new)

  call hions_expand_arrays(npart_new)
#ifdef CR
  call cr_expand_arrays(npart_new)
#endif
#ifdef FLUKA
  call fluka_mod_expand_arrays(npart_new, nele_new)
#endif
  call collimation_expand_arrays(npart_new, nblz_new)

  ! Update array size variables
  nele  = nele_new
  npart = npart_new
  nblz  = nblz_new
  nblo  = nblo_new

end subroutine expand_arrays

! Kicks off the allocation of the thick tracking arrays
subroutine allocate_thickarrays
  use parpro
  use mod_commonmn, only : mod_commonmn_allocate_thickarrays
  use mod_commons,  only : mod_commons_allocate_thickarrays
  implicit none
  call mod_commonmn_allocate_thickarrays
  call mod_commons_allocate_thickarrays
end subroutine allocate_thickarrays

! Kicks off the allocation of the thick tracking arrays
subroutine expand_thickarrays(nele_request, npart_request, nblz_request, nblo_request)
  use parpro
  use mod_commonmn, only : mod_commonmn_expand_thickarrays
  use mod_commons,  only : mod_commons_expand_thickarrays
  use mod_alloc
  implicit none

  integer, intent(in) :: nele_request, npart_request, nblz_request, nblo_request
  integer             :: nele_new, npart_new, nblz_new, nblo_new

  !Decide how much to actually expand -- set nele_new etc. based on some heuristic
  nele_new  = nele_request
  npart_new = npart_request
  nblz_new  = nblz_request
  nblo_new  = nblo_request

  call mod_commonmn_expand_thickarrays(nele_new, npart_new, nblo_new)
  call mod_commons_expand_thickarrays(nele_new, npart_new)

  !Update nele etc.
  nele = nele_new
! npart = npart_new
! nblz = nblz_new

end subroutine expand_thickarrays

! ================================================================================================ !
!  Compact Particle Arrays
!  A. Mereghetti, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-12-04
!  This routine is called to compact all relevant arrays when a particle is lost
! ================================================================================================ !
subroutine compactArrays

  use mod_hions
  use mod_common
  use mod_commont
  use mod_commonmn
  use aperture
  use collimation
#ifdef FLUKA
  use mod_fluka
#endif

  implicit none

  integer j,jj,jj1,ib2,ib3,lnapx

  ! Compact array
  lnapx=napx
  do j=napx,1,-1
    if(llostp(j)) then
      if(j.ne.lnapx) then
        do jj=j,lnapx-1
          jj1=jj+1
          nlostp(jj)=nlostp(jj1)
          xv1(jj)=xv1(jj1)
          xv2(jj)=xv2(jj1)
          yv1(jj)=yv1(jj1)
          yv2(jj)=yv2(jj1)
          dpsv(jj)=dpsv(jj1)
          sigmv(jj)=sigmv(jj1)
          ejfv(jj)=ejfv(jj1)
          ejv(jj)=ejv(jj1)
          rvv(jj)=rvv(jj1)
          ! ph: hisix
          nzz(jj)=nzz(jj1)
          naa(jj)=naa(jj1)
          nucm(jj)=nucm(jj1)
          mtc(jj)=mtc(jj1)
          moidpsv(jj)=moidpsv(jj1)
          omoidpsv(jj)=omoidpsv(jj1)
          ! ph: hisix
          oidpsv(jj)=oidpsv(jj1)
          dpsv1(jj)=dpsv1(jj1)
          clo6v(1,jj)=clo6v(1,jj1)
          clo6v(2,jj)=clo6v(2,jj1)
          clo6v(3,jj)=clo6v(3,jj1)
          clop6v(1,jj)=clop6v(1,jj1)
          clop6v(2,jj)=clop6v(2,jj1)
          clop6v(3,jj)=clop6v(3,jj1)

          !--beam-beam element
          di0xs(jj)=di0xs(jj1)
          dip0xs(jj)=dip0xs(jj1)
          di0zs(jj)=di0zs(jj1)
          dip0zs(jj)=dip0zs(jj1)
          do ib2=1,6
            do ib3=1,6
              tasau(jj,ib2,ib3)=tasau(jj1,ib2,ib3)
            end do
          end do

          ! Backtracking + aperture arrays
          ! These should get reset each time,
          ! but potentially there could be a collimator
          ! losing particles before the next usage
          ! So we compress these for now
          plost(jj) = plost(jj1)
          xLast(1,jj)   =  xLast(1,jj1)   ! position after last thick element [mm] (2,npart)
          xLast(1,jj)   =  xLast(1,jj1)   ! position after last thick element [mm] (2,npart)
          yLast(1,jj)   =  yLast(1,jj1)   ! angles after last thick element [mrad] (2,npart)
          yLast(2,jj)   =  yLast(2,jj1)   ! angles after last thick element [mrad] (2,npart)
          ejfvLast(jj)  =  ejfvLast(jj1)  ! linear momentum [MeV/c] (npart)
          ejvLast(jj)   =  ejvLast(jj1)   ! total energy [MeV] (npart)
          nucmLast(jj)  =  nucmLast(jj1)  ! nuclear mass [GeV/c2] (npart)
          sigmvLast(jj) =  sigmvLast(jj1) ! lag [mm] (npart)
          dpsvLast(jj)  =  dpsvLast(jj1)  ! (npart)
          naaLast(jj)   =  naaLast(jj1)   ! nuclear mass [] (npart)
          nzzLast(jj)   =  nzzLast(jj1)   ! atomic number [] (npart)


          if(do_coll) then
            ! If collimation is enabled,
            ! all the collimation arrays must also be compressed
            xgrd(jj)           = xgrd(jj1)
            ygrd(jj)           = ygrd(jj1)
            xpgrd(jj)          = xpgrd(jj1)
            ypgrd(jj)          = ypgrd(jj1)
            pgrd(jj)           = pgrd(jj1)
            ejfvgrd(jj)        = ejfvgrd(jj1)
            sigmvgrd(jj)       = sigmvgrd(jj1)
            rvvgrd(jj)         = rvvgrd(jj1)
            dpsvgrd(jj)        = dpsvgrd(jj1)
            oidpsvgrd(jj)      = oidpsvgrd(jj1)
            dpsv1grd(jj)       = dpsv1grd(jj1)
            part_hit_pos(jj)   = part_hit_pos(jj1)
            part_hit_turn(jj)  = part_hit_turn(jj1)
            part_abs_pos(jj)   = part_abs_pos(jj1)
            part_abs_turn(jj)  = part_abs_turn(jj1)
            part_select(jj)    = part_select(jj1)
            part_impact(jj)    = part_impact(jj1)
            part_indiv(jj)     = part_indiv(jj1)
            part_linteract(jj) = part_linteract(jj1)
            part_hit_before_pos(jj)  = part_hit_before_pos(jj1)
            part_hit_before_turn(jj) = part_hit_before_turn(jj1)
            secondary(jj)  = secondary(jj1)
            tertiary(jj)   = tertiary(jj1)
            other(jj)      = other(jj1)
            scatterhit(jj) = scatterhit(jj1)
            nabs_type(jj)  = nabs_type(jj1)
            !GRD HERE WE ADD A MARKER FOR THE PARTICLE FORMER NAME
            ipart(jj)      = ipart(jj1)
            flukaname(jj)  = flukaname(jj1)
            do ieff = 1, numeff
              counted_r(jj,ieff) = counted_r(jj1,ieff)
              counted_x(jj,ieff) = counted_x(jj1,ieff)
              counted_y(jj,ieff) = counted_y(jj1,ieff)
            end do
          endif

        end do !do jj=j,lnapx-1

#ifdef FLUKA
        if(fluka_enable) then
          call fluka_lostpart(lnapx, j) ! Inform fluka
        end if
#endif

      end if !if(j.ne.lnapx) then

      lnapx=lnapx-1
    end if !if(llostp(j)) then
  end do !do j=napx,1,-1
  napx=lnapx

end subroutine compactArrays
