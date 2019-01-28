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
#ifdef COLLIMAT
  use collimation,        only : collimation_allocate_arrays
#endif
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
#ifdef COLLIMAT
  call collimation_allocate_arrays
#endif

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
#ifdef FFIELD
  ! Modification by B.DALENA and T.PUGNAT
  use mod_ffield,         only : ffield_mod_expand_arrays
#endif
#ifdef COLLIMAT
  use collimation,        only : collimation_expand_arrays
#endif
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
#ifdef FFIELD
  call ffield_mod_expand_arrays(npart_new, nele_new)
#endif
#ifdef FFIELD
  call ffield_mod_expand_arrays(npart_new, nele_new)
#endif
#ifdef COLLIMAT
  call collimation_expand_arrays(npart_new, nblz_new)
#endif

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
