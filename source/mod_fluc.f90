! ================================================================================================ !
!  Random Fluctuation Starting Number Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  If besides mean values for the multipole errors (Gaussian) random errors should be considered,
!  this module is used to set the start value for the random generator.
!
!  Moved from main code and updated by V.K. Berglyd Olsen, June 2018
!  Last modified: 2018-06-13
! ================================================================================================ !
module mod_fluc

  use crcoall
  use floatPrecision
  use parpro,         only : nmac
  use mod_common,     only : izu0,mmac,mcut
  use mod_settings,   only : st_debug
  use sixtrack_input, only : sixin_echoVal

  implicit none

contains

subroutine fluc_parseInputLine(inLine, iLine, iErr, mout)

  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr
  integer,          intent(out)   :: mout

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "FLUC> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(iLine > 1) then
    write(lout,"(a)") "FLUC> ERROR This block only takes one line."
    iErr = .true.
    return
  end if

  if(nSplit > 0) call chr_cast(lnSplit(1),izu0,iErr)
  if(nSplit > 1) call chr_cast(lnSplit(2),mmac,iErr)
  if(nSplit > 2) call chr_cast(lnSplit(3),mout,iErr)
  if(nSplit > 3) call chr_cast(lnSplit(4),mcut,iErr)

  if(st_debug) then
    call sixin_echoVal("izu0",izu0,"FLUC",iLine)
    call sixin_echoVal("mmac",mmac,"FLUC",iLine)
    call sixin_echoVal("mout",mout,"FLUC",iLine)
    call sixin_echoVal("mcut",mcut,"FLUC",iLine)
  end if

  ! Process variables
  mcut = iabs(mcut)

  if(mmac > nmac) then
    write(lout,"(a,i0)") "FLUC> ERROR Maximum number of seeds for vectorisation is ",nmac
    iErr = .true.
    return
  end if

end subroutine fluc_parseInputLine

subroutine fluc_readFort8
end subroutine fluc_readFort8

subroutine fluc_readFort16

  implicit none

  logical isOpen
  integer mType, mNum, lineNo

  mType  = 0
  mNum   = 0
  lineNo = 0

  inquire(unit=16, opened=isOpen)
  if(isOpen) close(16)
  open(16,file="fort.16")
  rewind(16)


  close(16)

end subroutine fluc_readFort16

end module mod_fluc
