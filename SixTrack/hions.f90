! P. Hermes for Heavy Ion SixTrack (hiSix)
module mod_hions
  
  use floatPrecision
  use parpro
  use, intrinsic :: iso_fortran_env, only : int16, int32, int64
  use mod_alloc
  use numerical_constants, only : zero, one
  
  implicit none
  
  ! initialize the variables for ion tracking
  
  ! Checking for the HION block
  logical, save :: has_hion
  
  ! Rest mass of the reference ion species
  real(kind=fPrec), save :: nucm0,nucmda
  real(kind=fPrec), save :: brhono
  
  ! ien0,ien1: ion energy entering/leaving the collimator
  real(kind=fPrec), save :: ien0, ien1
  
  ! Rest mass of the tracked ion
  real(kind=fPrec), allocatable, save :: nucm(:) !(npart)
  
  ! Relative rigidity offset
  real(kind=fPrec), allocatable, save :: moidpsv(:) !(npart)
  real(kind=fPrec), allocatable, save :: omoidpsv(:) !(npart)
  
  ! Relative mass to charge ratio
  real(kind=fPrec), allocatable, save :: mtc(:) !(npart)
  
  ! Nucleon number of the reference ion species
  integer(kind=int16), save :: aa0
  
  ! Charge multiplicity of the reference ion species
  integer(kind=int16), save :: zz0
  
  integer(kind=int16), save :: nnuc0
  integer(kind=int16), save :: nnuc1
  
  ! Nucleon number of the tracked ion
  integer(kind=int16), allocatable, save :: naa(:) !(npart)
  
  ! Charge multiplicity of the tracked ion
  integer(kind=int16), allocatable, save :: nzz(:) !(npart)
  
  ! SixTrack particle IDs
  integer, allocatable, save :: pids(:) !(npart)
  
!  common/hions/zz0,aa0,nucm0,naa(npart),nzz(npart),nucm(npart),moidpsv(npart),omoidpsv(npart),mtc(npart),nnuc0,nnuc1,brhono,ien0,&
! &ien1,pids(npart)

contains

subroutine mod_hions_allocate_arrays
  
  implicit none
  
  ! Rest mass of the tracked ion
  call alloc(nucm,npart,nucm0,'nucm') !(npart)
  
  ! Relative rigidity offset
  call alloc(moidpsv,npart,one,'moidpsv') !(npart)
  call alloc(omoidpsv,npart,zero,'omoidpsv') !(npart)
  
  ! Relative mass to charge ratio
  call alloc(mtc,npart,one,'mtc') !(npart)
  
  ! Nucleon number of the tracked ion
  call alloc(naa,npart,aa0,'naa') !(npart)
  
  ! Charge multiplicity of the tracked ion
  call alloc(nzz,npart,zz0,'nzz') !(npart)
  print*, "nzzzzz", nzz(1)
  
  ! SixTrack particle IDs
  call alloc(pids,npart,0,'pids') !(npart)
  
end subroutine mod_hions_allocate_arrays

subroutine mod_hions_expand_arrays(npart_new)
  
  implicit none
  
  integer, intent(in) :: npart_new
  
  ! Rest mass of the tracked ion
  call resize(nucm,npart_new,nucm0,'nucm') !(npart)
  
  ! Relative rigidity offset
  call resize(moidpsv,npart_new,one,'moidpsv') !(npart)
  call resize(omoidpsv,npart_new,zero,'omoidpsv') !(npart)
  
  ! Relative mass to charge ratio
  call resize(mtc,npart_new,one,'mtc') !(npart)
  
  ! Nucleon number of the tracked ion
  call resize(naa,npart_new,aa0,'naa') !(npart)
  
  ! Charge multiplicity of the tracked ion
  call resize(nzz,npart_new,zz0,'nzz') !(npart)
  
  ! SixTrack particle IDs
  call resize(pids,npart_new,0,'pids') !(npart)
  
end subroutine mod_hions_expand_arrays

end module mod_hions
