! ================================================================================================ !
!  FILE UNITS MODULE
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-12-13
!
!  Module for keeping track of opened file units, their file names, and open the files correctly
!  depending on build flags.
! ================================================================================================ !
module mod_units

  implicit none

  ! Keep track of units
  integer, parameter           :: units_minUnit  = 1                ! First unit to keep track of
  integer, parameter           :: units_maxUnit  = 250              ! Last unit to keep track of
  integer, parameter           :: units_minAuto  = 100              ! First unit available for dynamic allocation
  integer, private, save       :: units_nextUnit = -1               ! Next unit available for dynamic allocation
  integer, private, save       :: units_logUnit  = -1               ! File unit for internal log file
  character(len=14), parameter :: units_logFile  = "file_units.dat" ! File unit for internal log file

  type, private :: unitRecord
    character(len=64), private :: file  = " "     ! The requested file name (not BOINC)
    character(len=3),  private :: mode  = " "     ! Read/write mode
    logical,           private :: taken = .false. ! Whether a unit is known to be taken or not
    logical,           private :: open  = .false. ! Whether file is opened by the module or not
    logical,           private :: fixed = .false. ! Whether the unit was requested as a fixed unit or not
  end type unitRecord

  ! Array to keep track of files
  type(unitRecord), private, save :: units_uList(units_minUnit:units_maxUnit)

contains

subroutine f_initUnits

  ! Grab the first unit for the mod_units log file
  units_logUnit  = units_minAuto
  units_nextUnit = units_nextUnit + 1

  units_uList(units_logUnit)%file  = units_logFile
  units_uList(units_logUnit)%mode  = "w"
  units_uList(units_logUnit)%taken = .true.
  units_uList(units_logUnit)%open  = .true.

  open(units_logUnit,file=units_logFile,form="formatted",status="replace",action="write")
  write(units_logUnit,"(a)") "# File Units Log"
  write(units_logUnit,"(a)") repeat("#",80)
  write(units_logUnit,"(a)") "ASSIGNED Unit ",units_logUnit," to '"//units_logFile//"'"
  flush(units_logUnit)

end subroutine f_initUnits

subroutine f_requestUnit(file,unit)

  character(len=*), intent(in)  :: file
  integer,          intent(out) :: unit

  integer i
  logical isOpen

  unit = -1
  call f_getUnit(file,unit)
  if(unit > 0) then
    write(*,"(a)") "REQUEST> Found #",unit," for file: '"//trim(file)//"'"
  end if

  do i=units_nextUnit:units_maxUnit
    if(units_uList(i)%taken) cycle
    inquire(unit=i, opened=isOpen)
    if(isOpen) then
    end if
  end do


end subroutine f_requestUnit

subroutine f_getUnit(file,unit)

  character(len=*), intent(in)  :: file
  integer,          intent(out) :: unit

  integer i

  unit = -1
  do i=units_minUnit,units_maxUnit
    if(units_uList(i)%file == file) then
      unit = i
      exit
    end if
  end do

end subroutine f_getUnit

subroutine f_open(unit,file,formatted,mode,err,status,recl)

  use crcoall

  implicit none

  integer,                    intent(in)  :: unit
  character(len=*),           intent(in)  :: file
  logical,                    intent(in)  :: formatted
  character(len=*),           intent(in)  :: mode
  logical,                    intent(out) :: err
  character(len=*), optional, intent(in)  :: status
  integer,          optional, intent(in)  :: recl

  ! type(unitSpec),   allocatable :: tmpUnits(:)
  character(len=:), allocatable :: fFileName, fStatus, fAction, fPosition
  character(len=256) :: tmpBoinc
  integer i, fRecl, nUnits, ioStat
  logical fFio, isOpen

  ! The code below breaks CR. Must look into later.
  ! inquire(unit=unit, opened=isOpen)
  ! if(isOpen) then
  !   write(lout,"(a,i0,a)") "UNITS> WARNING Attemting to open already opened unit ",unit," ... ignoring"
  !   return
  ! end if

  ! nUnits      = size(units_uList)
  ! units_nList = units_nList + 1
  ! if(units_nList > nUnits) then
  !   allocate(tmpUnits(units_nList + 10))
  !   tmpUnits(1:units_nList-1) = units_uList(1:units_nList-1)
  !   call move_alloc(tmpUnits,units_uList)
  ! end if

  if(present(recl)) then
    fRecl = recl
  else
    fRecl = 0
  end if

  if(present(status)) then
    fStatus = status
  else
    fStatus = "unknown"
  end if

#ifdef BOINC
  call boincrf(file,tmpBoinc)
  fFileName = trim(tmpBoinc)
#else
  fFileName = file
#endif
#ifdef FIO
  fFio = .true.
#else
  fFio = .false.
#endif
#ifndef NAGFOR
  fRecl = 0
#endif

  if(.not. formatted) fFio = .false.

  select case(mode)
  case("r")
    fAction   = "read"
    fPosition = "asis"
  case("w")
    fAction   = "write"
    fPosition = "asis"
  case("rw")
    fAction   = "readwrite"
    fPosition = "asis"
  case("r-")
    fAction   = "read"
    fPosition = "rewind"
  case("w-")
    fAction   = "write"
    fPosition = "rewind"
  case("rw-")
    fAction   = "readwrite"
    fPosition = "rewind"
  case("w+")
    fAction   = "write"
    fPosition = "append"
  case("rw+")
    fAction   = "readwrite"
    fPosition = "append"
  case default
    fAction   = "read"
    fPosition = "asis"
  end select

  ! units_uList(units_nList)%unit      = unit
  ! units_uList(units_nList)%filename  = file
  ! units_uList(units_nList)%formatted = formatted
  ! units_uList(units_nList)%mode      = mode
  ! units_uList(units_nList)%recl      = fRecl
  ! units_uList(units_nList)%open      = .true.

  err = .false.
  if(formatted) then
    if(fRecl > 0) then
      if(fFio) then
        open(unit,file=fFileName,form="formatted",status=fStatus,iostat=ioStat,&
          action=fAction,position=fPosition,round="nearest",recl=fRecl,err=10)
      else
        open(unit,file=fFileName,form="formatted",status=fStatus,iostat=ioStat,&
          action=fAction,position=fPosition,recl=fRecl,err=10)
      end if
    else
      if(fFio) then
        open(unit,file=fFileName,form="formatted",status=fStatus,iostat=ioStat,&
          action=fAction,position=fPosition,round="nearest",err=10)
      else
        open(unit,file=fFileName,form="formatted",status=fStatus,iostat=ioStat,&
          action=fAction,position=fPosition,err=10)
      end if
    end if
  else
    open(unit,file=fFileName,form="unformatted",status=fStatus,iostat=ioStat,&
      action=fAction,position=fPosition,err=10)
  endif

  if(ioStat /= 0) then
    err = .true.
    write(lout,"(a,i0)") "UNITS> File '"//trim(fFileName)//"' reported iostat = ",ioStat
  end if
  return

10 continue
  err = .true.
  write(lout,"(a)") "UNITS> File '"//trim(fFileName)//"' reported an error"

end subroutine f_open

subroutine f_close(unit)

  implicit none

  integer, intent(in) :: unit

  integer i
  logical isOpen

  inquire(unit=unit, opened=isOpen)
  if(isOpen) then
    flush(unit)
    close(unit)
  end if

  do i=1,units_nList
    if(units_uList(i)%unit == unit) then
      units_uList(i)%open = .false.
    end if
  end do

end subroutine f_close

subroutine f_flush(unit)

  implicit none

  integer, optional, intent(in) :: unit

  integer i
  logical isOpen

  if(present(unit)) then
    inquire(unit=unit, opened=isOpen)
    if(isOpen) flush(unit)
    return
  end if

  do i=1,units_nList
    inquire(unit=units_uList(i)%unit, opened=isOpen)
    if(isOpen) flush(units_uList(i)%unit)
  end do

end subroutine f_flush

end module mod_units
