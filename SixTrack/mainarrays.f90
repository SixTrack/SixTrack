! ================================================================================================ !
! Module for scaling array sizes defined in parpro
! by telling various other module to rescale their arrays;
! It is split off from parpro as a separate module
! in order to avoid circular dependencies.
! ================================================================================================ !
module parpro_scale
  
  use parpro
  use crcoall
  
  use mod_common,   only : mod_common_allocate_arrays,   mod_common_expand_arrays
  use mod_commonmn, only : mod_commonmn_allocate_arrays, mod_commonmn_expand_arrays
  use mod_commons,  only : mod_commons_allocate_arrays,  mod_commons_expand_arrays
  use mod_commond2, only : mod_commond2_allocate_arrays, mod_commond2_expand_arrays
  use aperture,     only : aperture_allocate_arrays,     aperture_expand_arrays
  use elens,        only : elens_allocate_arrays,        elens_expand_arrays
  use dump,         only : dump_allocate_arrays,         dump_expand_arrays
  use scatter,      only : scatter_allocate_arrays,      scatter_expand_arrays
  use bdex,         only : bdex_allocate_arrays,         bdex_expand_arrays
  use dynk,         only : dynk_allocate_arrays,         dynk_expand_arrays
  use wire,         only : wire_allocate_arrays,         wire_expand_arrays
  use mod_hions,    only : hions_allocate_arrays,    hions_expand_arrays
#ifdef FLUKA
  use mod_fluka,    only : fluka_mod_expand_arrays
#endif

implicit none

contains

! Allocate arrays scaling with the main memory parameters nele, npart, etc.
subroutine allocate_arrays
  
  implicit none
  
  ! Set initial values
  nele = nele_initial
  nblo = nblo_initial
  
  !Call subroutines to actually allocate
  call mod_common_allocate_arrays
  call mod_commonmn_allocate_arrays
  call mod_commons_allocate_arrays
  call mod_commond2_allocate_arrays
  call aperture_allocate_arrays
  call elens_allocate_arrays
  call dump_allocate_arrays
  call scatter_allocate_arrays
  call bdex_allocate_arrays
  call dynk_allocate_arrays
  call wire_allocate_arrays
  call hions_allocate_arrays
  
end subroutine allocate_arrays

! Change the allocation of the arrays scaling with the main memry parameter nele, npart, etc.
subroutine expand_arrays(nele_request, npart_request, nblz_request, nblo_request)
  
  use mod_alloc, only : alloc_log
  
  implicit none
  
  integer, intent(in) :: nele_request, npart_request, nblz_request, nblo_request
  integer             :: nele_new, npart_new, nblz_new, nblo_new

  !Decide how much to actually expand -- set nele_new etc. based on some heuristic
  nele_new = nele_request
  npart_new = npart_request
  nblz_new = nblz_request
  nblo_new = nblo_request

  write(alloc_log,"(a,4(1x,i0))") "ALLOC> Expanding - (nele,npart,nblz,nblo):",nele_new,npart_new,nblz_new,nblo_new
  
  !Call sub-subroutines to actually expand
  call mod_common_expand_arrays(nele_new, nblo_new)
  call mod_commonmn_expand_arrays(nele_new)
  call mod_commons_expand_arrays(nele_new)
  call mod_commond2_expand_arrays(nele_new)
  call aperture_expand_arrays(nele_new)
  call elens_expand_arrays(nele_new)
  call dump_expand_arrays(nele_new)
  call scatter_expand_arrays(nele_new)
  call bdex_expand_arrays(nele_new)
  call dynk_expand_arrays(nele_new)
  call wire_expand_arrays(nele_new)

  call hions_expand_arrays(npart_new)

#ifdef FLUKA
  call fluka_mod_expand_arrays(npart_new, nele_new)
#endif

  ! Update nele etc.
  nele = nele_new
! npart = npart_new
! nblz = nblz_new
  
end subroutine expand_arrays

! Kicks off the allocation of the thick tracking arrays
subroutine allocate_thickarrays
  use mod_commonmn, only : mod_commonmn_allocate_thickarrays
  use mod_commons,  only : mod_commons_allocate_thickarrays
  implicit none
  call mod_commonmn_allocate_thickarrays
  call mod_commons_allocate_thickarrays
end subroutine allocate_thickarrays

! Kicks off the allocation of the thick tracking arrays
subroutine expand_thickarrays(nele_request, npart_request, nblz_request, nblo_request)
  use mod_commonmn, only : mod_commonmn_expand_thickarrays
  use mod_commons,  only : mod_commons_expand_thickarrays
  use mod_alloc
  implicit none
  
  integer, intent(in) :: nele_request, npart_request, nblz_request, nblo_request
  integer             :: nele_new, npart_new, nblz_new, nblo_new
  
  !Decide how much to actually expand -- set nele_new etc. based on some heuristic
  nele_new = nele_request
  npart_new = npart_request
  nblz_new = nblz_request
  nblo_new = nblo_request
  
  call mod_commonmn_expand_thickarrays(nele_new, npart_new, nblo_new)
  call mod_commons_expand_thickarrays(nele_new, npart_new)
  
  !Update nele etc.
  nele = nele_new
! npart = npart_new
! nblz = nblz_new
  
end subroutine expand_thickarrays

end module parpro_scale
