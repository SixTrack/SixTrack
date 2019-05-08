! ================================================================================================ !
!  Heavy Ions Module
!  P. Hermes, J. Molson, T. Persson, V.K. Berglyd Olsen, BE-ABP
!  Last modified: 2018-06-02
! ================================================================================================ !
module mod_hions

  use floatPrecision
  use, intrinsic :: iso_fortran_env, only : int16

  implicit none

  ! Checking for the HION block
  logical, save :: has_hion = .false.

  ! ien0,ien1: ion energy entering/leaving the collimator
  real(kind=fPrec), save :: ien0, ien1

  integer(kind=int16), save :: nnuc0
  integer(kind=int16), save :: nnuc1

contains

subroutine hions_parseInputLine(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_common

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

end module mod_hions
