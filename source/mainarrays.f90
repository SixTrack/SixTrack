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
  use mod_common_track,        only : mod_commont_expand_arrays
  use mod_common_main,       only : mod_commonmn_expand_arrays
  use mod_commond2,       only : mod_commond2_expand_arrays
  use aperture,           only : aperture_expand_arrays
  use elens,              only : elens_allocate_arrays
  use cheby,              only : cheby_allocate_arrays
  use dump,               only : dump_expand_arrays
  use scatter,            only : scatter_expand_arrays
  use bdex,               only : bdex_allocate_arrays
  use dynk,               only : dynk_allocate_arrays
  use wire,               only : wire_expand_arrays
#ifdef CR
  use checkpoint_restart, only : cr_expand_arrays
#endif
#ifdef FLUKA
  use mod_fluka,          only : fluka_mod_expand_arrays
#endif
  use collimation,        only : collimation_allocate_arrays
  use coll_db,            only : cdb_expand_arrays

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
  call scatter_expand_arrays(nele,npart)
  call aperture_expand_arrays(nele,npart)

  call elens_allocate_arrays
  call cheby_allocate_arrays
  call bdex_allocate_arrays
  call dynk_allocate_arrays
#ifdef CR
  call cr_expand_arrays(npart)
#endif
  call collimation_allocate_arrays
  call cdb_expand_arrays(nele)

end subroutine allocate_arrays

! Change the allocation of the arrays scaling with the main memry parameter nele, npart, etc.
subroutine expand_arrays(nele_new, npart_new, nblz_new, nblo_new)

#ifdef DEBUG
  use mod_alloc, only : alloc_log
#endif

  use parpro
  use crcoall

  use mod_common,         only : mod_common_expand_arrays
  use mod_common_track,   only : mod_commont_expand_arrays
  use mod_common_main,    only : mod_commonmn_expand_arrays
  use mod_commond2,       only : mod_commond2_expand_arrays
  use aperture,           only : aperture_expand_arrays
  use elens,              only : elens_expand_arrays
  use cheby,              only : cheby_expand_arrays
  use dump,               only : dump_expand_arrays
  use scatter,            only : scatter_expand_arrays
  use bdex,               only : bdex_expand_arrays
  use dynk,               only : dynk_expand_arrays
  use wire,               only : wire_expand_arrays
#ifdef CR
  use checkpoint_restart, only : cr_expand_arrays
#endif
#ifdef FLUKA
  use mod_fluka,          only : fluka_mod_expand_arrays
#endif
  use mod_ffield,         only : ffield_mod_expand_arrays
  use collimation,        only : collimation_expand_arrays
  use coll_db,            only : cdb_expand_arrays

  implicit none

  integer, intent(in) :: nele_new
  integer, intent(in) :: npart_new
  integer, intent(in) :: nblz_new
  integer, intent(in) :: nblo_new

#ifdef DEBUG
  write(alloc_log,"(a,4(1x,i0))") "ALLOC> Expanding (nele,npart,nblz,nblo):",nele_new,npart_new,nblz_new,nblo_new
#endif

  !Call sub-subroutines to actually expand
  call mod_common_expand_arrays(nele_new,nblo_new,nblz_new,npart_new)
  call mod_commont_expand_arrays(nblz_new,npart_new)
  call mod_commonmn_expand_arrays(nblz_new,npart_new)
  call mod_commond2_expand_arrays(nele_new)

  call dump_expand_arrays(nele_new,nblz_new)
  call wire_expand_arrays(nele_new,nblz_new)
  call scatter_expand_arrays(nele_new,npart_new)
  call aperture_expand_arrays(nele_new, npart_new)

  call elens_expand_arrays(nele_new)
  call cheby_expand_arrays(nele_new)
  call bdex_expand_arrays(nele_new)
  call dynk_expand_arrays(nele_new)

#ifdef CR
  call cr_expand_arrays(npart_new)
#endif
#ifdef FLUKA
  call fluka_mod_expand_arrays(npart_new, nele_new)
#endif
  call ffield_mod_expand_arrays(npart_new, nele_new)
  call collimation_expand_arrays(npart_new, nblz_new)
  call cdb_expand_arrays(nele_new)

  ! Update array size variables
  nele  = nele_new
  npart = npart_new
  nblz  = nblz_new
  nblo  = nblo_new

end subroutine expand_arrays

! Kicks off the allocation of the thick tracking arrays
subroutine allocate_thickarrays
  use parpro
  use mod_common_main, only : mod_commonmn_expand_thickarrays
  use mod_commons,     only : mod_commons_expand_thickarrays
  implicit none
  call mod_commonmn_expand_thickarrays(nele, npart, nblo)
  call mod_commons_expand_thickarrays(nele, npart)
end subroutine allocate_thickarrays

! Kicks off the allocation of the thick tracking arrays
subroutine expand_thickarrays(nele_request, npart_request, nblz_request, nblo_request)
  use parpro
  use mod_common_main, only : mod_commonmn_expand_thickarrays
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
!  Shuffle Lost Particles
!  A. Mereghetti, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-12-04
!  This routine is called to move all lost particles to the end of particle arrays (after napx)
!  Routine has ben renamed from compactArrays.
! ================================================================================================ !
subroutine shuffleLostParticles

  use floatPrecision
  use mod_common
  use mod_common_track
  use mod_common_main
  use aperture
  use collimation
#ifdef FLUKA
  use mod_fluka
#endif

  implicit none

  integer j, napx_new, tnapx
  logical, allocatable :: tmp_lostP(:)

  napx_new = napx

  allocate(tmp_lostP(npart))
  tmp_lostP(1:npart) = llostp(1:npart)

  tnapx = napx
  do j=napx,1,-1
    if(llostp(j) .eqv. .false.) cycle

    ! Move lost particle to the back
    partID(j:tnapx)    = cshift(partID(j:tnapx),    1)
    parentID(j:tnapx)  = cshift(parentID(j:tnapx),  1)
    tmp_lostP(j:tnapx) = cshift(tmp_lostP(j:tnapx), 1)

    ! Main Particle Arrays
    xv1(j:tnapx)       = cshift(xv1(j:tnapx),       1)
    xv2(j:tnapx)       = cshift(xv2(j:tnapx),       1)
    yv1(j:tnapx)       = cshift(yv1(j:tnapx),       1)
    yv2(j:tnapx)       = cshift(yv2(j:tnapx),       1)
    dpsv(j:tnapx)      = cshift(dpsv(j:tnapx),      1)
    sigmv(j:tnapx)     = cshift(sigmv(j:tnapx),     1)
    ejfv(j:tnapx)      = cshift(ejfv(j:tnapx),      1)
    ejv(j:tnapx)       = cshift(ejv(j:tnapx),       1)
    rvv(j:tnapx)       = cshift(rvv(j:tnapx),       1)

    ! Ion Arrays
    nzz(j:tnapx)       = cshift(nzz(j:tnapx),       1)
    naa(j:tnapx)       = cshift(naa(j:tnapx),       1)
    nqq(j:tnapx)       = cshift(nqq(j:tnapx),       1)
    pdgid(j:tnapx)     = cshift(pdgid(j:tnapx),     1)
    nucm(j:tnapx)      = cshift(nucm(j:tnapx),      1)
    mtc(j:tnapx)       = cshift(mtc(j:tnapx),       1)
    dpsv1(j:tnapx)     = cshift(dpsv1(j:tnapx),     1)
    oidpsv(j:tnapx)    = cshift(oidpsv(j:tnapx),    1)
    moidpsv(j:tnapx)   = cshift(moidpsv(j:tnapx),   1)
    omoidpsv(j:tnapx)  = cshift(omoidpsv(j:tnapx),  1)

    ! Backtracking + Aperture
    plost(j:tnapx)     = cshift(plost(j:tnapx),     1)
    xLast(:,j:tnapx)   = cshift(xLast(:,j:tnapx),   1, 2)
    yLast(:,j:tnapx)   = cshift(yLast(:,j:tnapx),   1, 2)
    ejfvLast(j:tnapx)  = cshift(ejfvLast(j:tnapx),  1)
    ejvLast(j:tnapx)   = cshift(ejvLast(j:tnapx),   1)
    nucmLast(j:tnapx)  = cshift(nucmLast(j:tnapx),  1)
    sigmvLast(j:tnapx) = cshift(sigmvLast(j:tnapx), 1)
    dpsvLast(j:tnapx)  = cshift(dpsvLast(j:tnapx),  1)
    naaLast(j:tnapx)   = cshift(naaLast(j:tnapx),   1)
    nzzLast(j:tnapx)   = cshift(nzzLast(j:tnapx),   1)
    nqqLast(j:tnapx)   = cshift(nqqLast(j:tnapx),   1)
    pdgidLast(j:tnapx) = cshift(pdgidLast(j:tnapx), 1)

    tnapx = tnapx - 1
  end do
  napx_new = tnapx

  ! Collimation
  if(do_coll) then
    tnapx = napx
    do j=napx,1,-1
      if(llostp(j) .eqv. .false.) cycle

      xgrd(j:tnapx)                 = cshift(xgrd(j:tnapx),                 1)
      ygrd(j:tnapx)                 = cshift(ygrd(j:tnapx),                 1)
      xpgrd(j:tnapx)                = cshift(xpgrd(j:tnapx),                1)
      ypgrd(j:tnapx)                = cshift(ypgrd(j:tnapx),                1)
      pgrd(j:tnapx)                 = cshift(pgrd(j:tnapx),                 1)
      ejfvgrd(j:tnapx)              = cshift(ejfvgrd(j:tnapx),              1)
      sigmvgrd(j:tnapx)             = cshift(sigmvgrd(j:tnapx),             1)
      rvvgrd(j:tnapx)               = cshift(rvvgrd(j:tnapx),               1)
      dpsvgrd(j:tnapx)              = cshift(dpsvgrd(j:tnapx),              1)
      oidpsvgrd(j:tnapx)            = cshift(oidpsvgrd(j:tnapx),            1)
      dpsv1grd(j:tnapx)             = cshift(dpsv1grd(j:tnapx),             1)

      part_hit_pos(j:tnapx)         = cshift(part_hit_pos(j:tnapx),         1)
      part_hit_turn(j:tnapx)        = cshift(part_hit_turn(j:tnapx),        1)
      part_abs_pos(j:tnapx)         = cshift(part_abs_pos(j:tnapx),         1)
      part_abs_turn(j:tnapx)        = cshift(part_abs_turn(j:tnapx),        1)
      part_select(j:tnapx)          = cshift(part_select(j:tnapx),          1)
      part_impact(j:tnapx)          = cshift(part_impact(j:tnapx),          1)
      part_indiv(j:tnapx)           = cshift(part_indiv(j:tnapx),           1)
      part_linteract(j:tnapx)       = cshift(part_linteract(j:tnapx),       1)
      part_hit_before_pos(j:tnapx)  = cshift(part_hit_before_pos(j:tnapx),  1)
      part_hit_before_turn(j:tnapx) = cshift(part_hit_before_turn(j:tnapx), 1)

      secondary(j:tnapx)            = cshift(secondary(j:tnapx),            1)
      tertiary(j:tnapx)             = cshift(tertiary(j:tnapx),             1)
      other(j:tnapx)                = cshift(other(j:tnapx),                1)
      scatterhit(j:tnapx)           = cshift(scatterhit(j:tnapx),           1)
      nabs_type(j:tnapx)            = cshift(nabs_type(j:tnapx),            1)
      ipart(j:tnapx)                = cshift(ipart(j:tnapx),                1)
      flukaname(j:tnapx)            = cshift(flukaname(j:tnapx),            1)

      counted_r(j:tnapx,:)          = cshift(counted_r(j:tnapx,:),          1, 1)
      counted_x(j:tnapx,:)          = cshift(counted_x(j:tnapx,:),          1, 1)
      counted_y(j:tnapx,:)          = cshift(counted_y(j:tnapx,:),          1, 1)

      tnapx = tnapx - 1
    end do
  end if

#ifdef FLUKA
  if(fluka_enable) then
    tnapx = napx
    do j=napx,1,-1
      if(llostp(j)) then
        call fluka_shuffleLostParticles(tnapx, j) ! Inform fluka
        tnapx = tnapx - 1
      end if
    end do
  end if
#endif

  napx = napx_new
  call move_alloc(tmp_lostP, llostp)

end subroutine shuffleLostParticles
