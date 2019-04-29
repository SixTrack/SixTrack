! ================================================================================================ !
!  Heavy Ions Module
!  P. Hermes, J. Molson, T. Persson, V.K. Berglyd Olsen, BE-ABP
!  Last modified: 2018-06-02
! ================================================================================================ !
module mod_hions

  use floatPrecision
  use parpro
  use mod_alloc
  use physical_constants, only : pmap
  use numerical_constants, only : zero, one, c1e3
  use, intrinsic :: iso_fortran_env, only : int16

  implicit none

  ! Checking for the HION block
  logical, save :: has_hion = .false.

  ! Rest mass of the reference ion species
  real(kind=fPrec), save :: nucm0 = pmap
  real(kind=fPrec), save :: nucmda

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
  integer(kind=int16), save :: aa0 = 1

  ! Charge multiplicity of the reference ion species
  integer(kind=int16), save :: zz0 = 1

  integer(kind=int16), save :: nnuc0
  integer(kind=int16), save :: nnuc1

  ! Nucleon number of the tracked ion
  integer(kind=int16), allocatable, save :: naa(:) !(npart)

  ! Charge multiplicity of the tracked ion
  integer(kind=int16), allocatable, save :: nzz(:) !(npart)

  ! SixTrack particle IDs
  integer, allocatable, save :: pids(:) !(npart)

#ifdef CR
  real(kind=fPrec),                 save :: nucmda_cr
  real(kind=fPrec),                 save :: ien0_cr
  real(kind=fPrec),                 save :: ien1_cr
  integer(kind=int16),              save :: nnuc0_cr
  integer(kind=int16),              save :: nnuc1_cr
  real(kind=fPrec),    allocatable, save :: nucm_cr(:)
  real(kind=fPrec),    allocatable, save :: moidpsv_cr(:)
  real(kind=fPrec),    allocatable, save :: omoidpsv_cr(:)
  real(kind=fPrec),    allocatable, save :: mtc_cr(:)
  integer(kind=int16), allocatable, save :: naa_cr(:)
  integer(kind=int16), allocatable, save :: nzz_cr(:)
  integer,             allocatable, save :: pids_cr(:)
#endif

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
  call alloc(nucm,npart_new,nucm0,'nucm')
  call alloc(moidpsv,npart_new,one,'moidpsv')
  call alloc(omoidpsv,npart_new,zero,'omoidpsv')
  call alloc(mtc,npart_new,one,'mtc')
  call alloc(naa,npart_new,aa0,'naa')
  call alloc(nzz,npart_new,zz0,'nzz')
  call alloc(pids,npart_new,0,'pids')
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
    write(lerr,"(a)") "HIONS> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  if(iLine > 1) then
    write(lout,"(a)") "HIONS> WARNING Only expected one input line."
  end if

  if(nSplit /= 3) then
    write(lerr,"(a,i0)") "HIONS> ERROR Line must have 3 values, got ",nSplit
    iErr = .true.
    return
  end if

  call chr_cast(lnSplit(1),aa0,  iErr)
  call chr_cast(lnSplit(2),zz0,  iErr)
  call chr_cast(lnSplit(3),nucm0,iErr)

  nucm0 = nucm0*c1e3 ! [GeV/c^2] -> [MeV/c^2]

end subroutine hions_parseInputLine

subroutine hions_postInput

  use mod_common, only : pma

  if(.not. has_hion) then
    ! If we don't have the HION block, we need to set some variables - default to the proton values
    zz0   = 1
    aa0   = 1
    nucm0 = pma
    write(lout,"(a)")        "HION> No HION block found. Defaulting to the proton values: "
    write(lout,"(a,i0)")     "HION>  * Z = ",zz0
    write(lout,"(a,i0)")     "HION>  * A = ",aa0
    write(lout,"(a,e22.15)") "HION>  * M = ",nucm0
  end if

  ! Init arrays
  mtc(:)      = one
  naa(:)      = aa0
  nzz(:)      = zz0
  nucm(:)     = nucm0
  moidpsv(:)  = one
  omoidpsv(:) = zero

end subroutine hions_postInput

#ifdef CR
subroutine hions_crpoint(fileUnit, writeErr, iErro)

  implicit none

  integer, intent(in)    :: fileUnit
  logical, intent(inout) :: writeErr
  integer, intent(inout) :: iErro

  integer i

  write(fileUnit,err=10,iostat=iErro) nucmda,ien0,ien1,nnuc0,nnuc1
  write(fileUnit,err=10,iostat=iErro) (nucm(i),     i=1, npart)
  write(fileUnit,err=10,iostat=iErro) (moidpsv(i),  i=1, npart)
  write(fileUnit,err=10,iostat=iErro) (omoidpsv(i), i=1, npart)
  write(fileUnit,err=10,iostat=iErro) (mtc(i),      i=1, npart)
  write(fileUnit,err=10,iostat=iErro) (naa(i),      i=1, npart)
  write(fileUnit,err=10,iostat=iErro) (nzz(i),      i=1, npart)
  write(fileUnit,err=10,iostat=iErro) (pids(i),     i=1, npart)
  endfile(fileUnit,iostat=iErro)
  backspace(fileUnit,iostat=iErro)

  return

10 continue
  writeErr = .true.
  write(lout,"(a,i0,a)") "SIXTRACR> ERROR Reading C/R file fort.",fileUnit," in HIONS"
  write(93,  "(a,i0,a)") "SIXTRACR> ERROR Reading C/R file fort.",fileUnit," in HIONS"
  flush(93)

end subroutine hions_crpoint

subroutine hions_crcheck_readdata(fileUnit, readErr)

  use crcoall

  implicit none

  integer, intent(in)  :: fileUnit
  logical, intent(out) :: readErr

  integer i

  call alloc(nucm_cr,    npart,nucm0,"nucm_cr")
  call alloc(moidpsv_cr, npart,one,  "moidpsv_cr")
  call alloc(omoidpsv_cr,npart,zero, "omoidpsv_cr")
  call alloc(mtc_cr,     npart,one,  "mtc_cr")
  call alloc(naa_cr,     npart,aa0,  "naa_cr")
  call alloc(nzz_cr,     npart,zz0,  "nzz_cr")
  call alloc(pids_cr,    npart,0,    "pids_cr")

  read(fileunit,err=10,end=10) nucmda_cr,ien0_cr,ien1_cr,nnuc0_cr,nnuc1_cr
  read(fileunit,err=10,end=10) (nucm_cr(i),     i=1, npart)
  read(fileunit,err=10,end=10) (moidpsv_cr(i),  i=1, npart)
  read(fileunit,err=10,end=10) (omoidpsv_cr(i), i=1, npart)
  read(fileunit,err=10,end=10) (mtc_cr(i),      i=1, npart)
  read(fileunit,err=10,end=10) (naa_cr(i),      i=1, npart)
  read(fileunit,err=10,end=10) (nzz_cr(i),      i=1, npart)
  read(fileunit,err=10,end=10) (pids_cr(i),     i=1, npart)

  readErr = .false.
  return

10 continue
  readErr = .true.
  write(lout,"(a,i0,a)") "SIXTRACR> ERROR Reading C/R file fort.",fileUnit," in HIONS"
  write(93,  "(a,i0,a)") "SIXTRACR> ERROR Reading C/R file fort.",fileUnit," in HIONS"
  flush(93)

end subroutine hions_crcheck_readdata

subroutine hions_crstart

  implicit none

  nucmda = nucmda_cr
  ien0   = ien0_cr
  ien1   = ien1_cr
  nnuc0  = nnuc0_cr
  nnuc1  = nnuc1_cr

  nucm(1:npart)     = nucm_cr(1:npart)
  moidpsv(1:npart)  = moidpsv_cr(1:npart)
  omoidpsv(1:npart) = omoidpsv_cr(1:npart)
  mtc(1:npart)      = mtc_cr(1:npart)
  naa(1:npart)      = naa_cr(1:npart)
  nzz(1:npart)      = nzz_cr(1:npart)
  pids(1:npart)     = pids_cr(1:npart)

  call dealloc(nucm_cr,    "nucm_cr")
  call dealloc(moidpsv_cr, "moidpsv_cr")
  call dealloc(omoidpsv_cr,"omoidpsv_cr")
  call dealloc(mtc_cr,     "mtc_cr")
  call dealloc(naa_cr,     "naa_cr")
  call dealloc(nzz_cr,     "nzz_cr")
  call dealloc(pids_cr,    "pids_cr")

end subroutine hions_crstart
#endif
end module mod_hions
