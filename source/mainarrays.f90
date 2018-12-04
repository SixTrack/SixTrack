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

  use floatPrecision
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

  integer j, napx_new
  logical, allocatable :: tmp_lostP(:)

  napx_new = napx

  allocate(tmp_lostP(npart))
  tmp_lostP(1:npart) = llostp(1:npart)

  do j=napx,1,-1
    if(llostp(j) .eqv. .false.) cycle

    ! Move lost particle to the back
    nlostp(j:npart)    = cshift(nlostp(j:npart),    1)
    tmp_lostP(j:npart) = cshift(tmp_lostP(j:npart), 1)

    ! Main Particle Arrays
    xv1(j:npart)       = cshift(xv1(j:npart),       1)
    xv2(j:npart)       = cshift(xv2(j:npart),       1)
    yv1(j:npart)       = cshift(yv1(j:npart),       1)
    yv2(j:npart)       = cshift(yv2(j:npart),       1)
    dpsv(j:npart)      = cshift(dpsv(j:npart),      1)
    sigmv(j:npart)     = cshift(sigmv(j:npart),     1)
    ejfv(j:npart)      = cshift(ejfv(j:npart),      1)
    ejv(j:npart)       = cshift(ejv(j:npart),       1)
    rvv(j:npart)       = cshift(rvv(j:npart),       1)

    ! Ion Arrays
    nzz(j:npart)       = cshift(nzz(j:npart),       1)
    naa(j:npart)       = cshift(naa(j:npart),       1)
    nucm(j:npart)      = cshift(nucm(j:npart),      1)
    mtc(j:npart)       = cshift(mtc(j:npart),       1)
    dpsv1(j:npart)     = cshift(dpsv1(j:npart),     1)
    oidpsv(j:npart)    = cshift(oidpsv(j:npart),    1)
    moidpsv(j:npart)   = cshift(moidpsv(j:npart),   1)
    omoidpsv(j:npart)  = cshift(omoidpsv(j:npart),  1)

    ! Beam--Beam
    di0xs(j:npart)     = cshift(di0xs(j:npart),     1)
    dip0xs(j:npart)    = cshift(dip0xs(j:npart),    1)
    di0zs(j:npart)     = cshift(di0zs(j:npart),     1)
    dip0zs(j:npart)    = cshift(dip0zs(j:npart),    1)
    tasau(j:npart,:,:) = cshift(tasau(j:npart,:,:), 1, 1)

    ! Closed Orbit
    clo6v(:,j:npart)   = cshift(clo6v(:,j:npart),   1, 2)
    clop6v(:,j:npart)  = cshift(clop6v(:,j:npart),  1, 2)

    ! Backtracking + Aperture
    plost(j:npart)     = cshift(plost(j:npart),     1)
    xLast(:,j:npart)   = cshift(xLast(:,j:npart),   1, 2)
    yLast(:,j:npart)   = cshift(yLast(:,j:npart),   1, 2)
    ejfvLast(j:npart)  = cshift(ejfvLast(j:npart),  1)
    ejvLast(j:npart)   = cshift(ejvLast(j:npart),   1)
    nucmLast(j:npart)  = cshift(nucmLast(j:npart),  1)
    sigmvLast(j:npart) = cshift(sigmvLast(j:npart), 1)
    dpsvLast(j:npart)  = cshift(dpsvLast(j:npart),  1)
    naaLast(j:npart)   = cshift(naaLast(j:npart),   1)
    nzzLast(j:npart)   = cshift(nzzLast(j:npart),   1)

    napx_new = napx_new - 1
  end do

  ! Collimation
  if(do_coll) then
    do j=napx,1,-1
      if(llostp(j) .eqv. .false.) cycle

      xgrd(j:npart)                 = cshift(xgrd(j:npart),                 1)
      ygrd(j:npart)                 = cshift(ygrd(j:npart),                 1)
      xpgrd(j:npart)                = cshift(xpgrd(j:npart),                1)
      ypgrd(j:npart)                = cshift(ypgrd(j:npart),                1)
      pgrd(j:npart)                 = cshift(pgrd(j:npart),                 1)
      ejfvgrd(j:npart)              = cshift(ejfvgrd(j:npart),              1)
      sigmvgrd(j:npart)             = cshift(sigmvgrd(j:npart),             1)
      rvvgrd(j:npart)               = cshift(rvvgrd(j:npart),               1)
      dpsvgrd(j:npart)              = cshift(dpsvgrd(j:npart),              1)
      oidpsvgrd(j:npart)            = cshift(oidpsvgrd(j:npart),            1)
      dpsv1grd(j:npart)             = cshift(dpsv1grd(j:npart),             1)

      part_hit_pos(j:npart)         = cshift(part_hit_pos(j:npart),         1)
      part_hit_turn(j:npart)        = cshift(part_hit_turn(j:npart),        1)
      part_abs_pos(j:npart)         = cshift(part_abs_pos(j:npart),         1)
      part_abs_turn(j:npart)        = cshift(part_abs_turn(j:npart),        1)
      part_select(j:npart)          = cshift(part_select(j:npart),          1)
      part_impact(j:npart)          = cshift(part_impact(j:npart),          1)
      part_indiv(j:npart)           = cshift(part_indiv(j:npart),           1)
      part_linteract(j:npart)       = cshift(part_linteract(j:npart),       1)
      part_hit_before_pos(j:npart)  = cshift(part_hit_before_pos(j:npart),  1)
      part_hit_before_turn(j:npart) = cshift(part_hit_before_turn(j:npart), 1)

      secondary(j:npart)            = cshift(secondary(j:npart),            1)
      tertiary(j:npart)             = cshift(tertiary(j:npart),             1)
      other(j:npart)                = cshift(other(j:npart),                1)
      scatterhit(j:npart)           = cshift(scatterhit(j:npart),           1)
      nabs_type(j:npart)            = cshift(nabs_type(j:npart),            1)
      ipart(j:npart)                = cshift(ipart(j:npart),                1)
      flukaname(j:npart)            = cshift(flukaname(j:npart),            1)

      counted_r(j:npart,:)          = cshift(counted_r(j:npart,:),          1, 1)
      counted_x(j:npart,:)          = cshift(counted_x(j:npart,:),          1, 1)
      counted_y(j:npart,:)          = cshift(counted_y(j:npart,:),          1, 1)
    end do
  end if

#ifdef FLUKA
  if(fluka_enable) then
    do j=napx,1,-1
      if(llostp(j)) call fluka_lostpart(napx, j) ! Inform fluka
    end do
  end if
#endif

  napx = napx_new
  call move_alloc(tmp_lostP, llostp)

end subroutine compactArrays
