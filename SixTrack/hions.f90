! ================================================================================================ !
!  Heavy Ions Module
!  P. Hermes, J. Molson, T. Persson, V.K. Berglyd Olsen, BE-ABP
!  Last modified: 2018-06-02
! ================================================================================================ !
module mod_hions
  
  use floatPrecision
  use parpro
  use mod_alloc
  use numerical_constants, only : zero, one, c1e3
  
  implicit none
  
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
  integer, save :: aa0
  
  ! Charge multiplicity of the reference ion species
  integer, save :: zz0
  
  integer, save :: nnuc0
  integer, save :: nnuc1
  
  ! Nucleon number of the tracked ion
  integer, allocatable, save :: naa(:) !(npart)
  
  ! Charge multiplicity of the tracked ion
  integer, allocatable, save :: nzz(:) !(npart)
  
  ! SixTrack particle IDs
  integer, allocatable, save :: pids(:) !(npart)
  
contains

subroutine hions_allocate_arrays
  call alloc(nucm,npart,nucm0,'nucm')
  call alloc(moidpsv,npart,one,'moidpsv')
  call alloc(omoidpsv,npart,zero,'omoidpsv')
  call alloc(mtc,npart,one,'mtc')
  call alloc(naa,npart,aa0,'naa')
  call alloc(nzz,npart,zz0,'nzz')
  call alloc(pids,npart,0,'pids')
end subroutine hions_allocate_arrays

subroutine hions_expand_arrays(npart_new)
  integer, intent(in) :: npart_new
  call resize(nucm,npart_new,nucm0,'nucm')
  call resize(moidpsv,npart_new,one,'moidpsv')
  call resize(omoidpsv,npart_new,zero,'omoidpsv')
  call resize(mtc,npart_new,one,'mtc')
  call resize(naa,npart_new,aa0,'naa')
  call resize(nzz,npart_new,zz0,'nzz')
  call resize(pids,npart_new,0,'pids')
end subroutine hions_expand_arrays

subroutine hions_parseInputLine(inLine, iLine, iErr)
  
  use string_tools
  
  implicit none
  
  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr
  
  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit
  logical spErr
  
  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "HIONS> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  
  if(iLine > 1) then
    write(lout,"(a)") "HIONS> WARNING Only expected one input line."
  end if
  
  if(nSplit /= 3) then
    write(lout,"(a,i0)") "HIONS> ERROR Line must have 3 values, got ",nSplit
    iErr = .true.
    return
  end if
  
  call chr_cast(lnSplit(1),aa0,  iErr)
  call chr_cast(lnSplit(2),zz0,  iErr)
  call chr_cast(lnSplit(3),nucm0,iErr)
  
  nucm0 = nucm0*c1e3 ! [GeV/c^2] -> [MeV/c^2]
  
end subroutine hions_parseInputLine

end module mod_hions
