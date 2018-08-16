! ================================================================================================ !
!  Read Input
! ~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-18
!  Module for parsing special input files like fort.33, fort.13, fort.10, etc.
! ================================================================================================ !

module read_input

  use crcoall
  use string_tools
  use sixtrack_input, only : sixin_echoVal
  use mod_settings
  use mod_units
  use parpro

  implicit none

contains

subroutine readFort33

  use mod_common, only : clo6, clop6

  implicit none

  character(len=:), allocatable :: lnSplit(:)
  character(len=mInputLn)       :: inLine
  integer nSplit, ioStat
  logical spErr, fErr

  call units_openUnit(33,"fort.33",.true.,"r",fErr)

  read(33,"(a)",iostat=ioStat) inLine
  if(ioStat /= 0) then
    write(lout,"(a,i0)") "READ33> ERROR Failed to read line from 'fort.33'. iostat = ",ioStat
    call prror(-1)
  end if

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "READ33> ERROR Failed to parse line from 'fort.33'"
    call prror(-1)
  end if

  if(nSplit > 0) call chr_cast(lnSplit(1),clo6(1), spErr)
  if(nSplit > 1) call chr_cast(lnSplit(2),clop6(1),spErr)
  if(nSplit > 2) call chr_cast(lnSplit(3),clo6(2), spErr)
  if(nSplit > 3) call chr_cast(lnSplit(4),clop6(2),spErr)
  if(nSplit > 4) call chr_cast(lnSplit(5),clo6(3), spErr)
  if(nSplit > 5) call chr_cast(lnSplit(6),clop6(3),spErr)

  if(st_debug) then
    call sixin_echoVal("clo6(1)", clo6(1), "READ33",1)
    call sixin_echoVal("clop6(1)",clop6(1),"READ33",1)
    call sixin_echoVal("clo6(2)", clo6(2), "READ33",1)
    call sixin_echoVal("clop6(2)",clop6(2),"READ33",1)
    call sixin_echoVal("clo6(3)", clo6(3), "READ33",1)
    call sixin_echoVal("clop6(3)",clop6(3),"READ33",1)
  end if

end subroutine readFort33

end module read_input
