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
  integer, parameter           :: units_minUnit  = 1                 ! First unit to keep track of
  integer, parameter           :: units_maxUnit  = 250               ! Last unit to keep track of
  integer, parameter           :: units_minAuto  = 100               ! First unit available for dynamic allocation
  integer, private, save       :: units_nextUnit = units_minAuto     ! Next unit available for dynamic allocation
  integer, private, save       :: units_logUnit  = units_maxUnit     ! File unit for internal log file
  character(len=14), parameter :: units_logFile  = "file_units.log"  ! File name for internal log file

  type, private :: unitRecord
    character(len=64), private :: file  = " "     ! The requested file name (not BOINC)
    character(len=3),  private :: mode  = " "     ! Read/write mode
    logical,           private :: taken = .false. ! Whether a unit is known to be taken or not
    logical,           private :: open  = .false. ! Whether file is opened by the module or not
    logical,           private :: fixed = .true.  ! Whether the unit was requested as a fixed unit or not
  end type unitRecord

  ! Array to keep track of files
  type(unitRecord), private, save :: units_uList(units_minUnit:units_maxUnit)

contains

subroutine f_initUnits

  ! All we need to do is open the log file

  units_uList(units_logUnit)%file  = units_logFile
  units_uList(units_logUnit)%mode  = "w"
  units_uList(units_logUnit)%taken = .true.
  units_uList(units_logUnit)%open  = .true.
  units_uList(units_logUnit)%fixed = .true.

  open(units_logUnit,file=units_logFile,form="formatted",status="replace",action="write")
  write(units_logUnit,"(a)") "# File Units Log"
  write(units_logUnit,"(a)") repeat("#",100)
  write(units_logUnit,"(a)") "#   AtTime  Action      Status    Unit  FileName"
  flush(units_logUnit)
  call f_writeLog("ASSIGNED","OK",units_logUnit,units_logFile)

end subroutine f_initUnits

subroutine f_requestUnit(file,unit)

  use, intrinsic :: iso_fortran_env, only : error_unit

  character(len=*), intent(in)  :: file
  integer,          intent(out) :: unit

  integer i
  logical isOpen

  if(len_trim(file) > 64) then
    write(error_unit,"(a,i0)") "UNITS> ERROR Max length of file name in f_requestUnit is 64 characters, got ",len(file)
    call prror
  end if

  unit = -1
  call f_getUnit(file,unit)
  if(unit > 0) then
    call f_writeLog("REQUEST","EXISTS",unit,trim(file))
    return
  end if

  do i=units_nextUnit,units_maxUnit
    if(units_uList(i)%taken) cycle
    inquire(unit=i, opened=isOpen)
    if(isOpen) then
      if(units_uList(i)%taken .eqv. .false.) then
        units_uList(i)%file  = "unknown"
        units_uList(i)%mode  = ""
        units_uList(i)%taken = .true.
        units_uList(i)%open  = .false.
        units_uList(i)%fixed = .false.
        call f_writeLog("REQUEST","TAKEN",i,trim(file))
      end if
    else
      unit = i
      exit
    end if
  end do

  if(unit > 0) then
    call f_writeLog("REQUEST","NEW",unit,trim(file))
    units_uList(unit)%file  = trim(file)
    units_uList(unit)%mode  = ""
    units_uList(unit)%taken = .true.
    units_uList(unit)%open  = .false.
    units_uList(unit)%fixed = .false.
end if

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

  use, intrinsic :: iso_fortran_env, only : error_unit

  implicit none

  integer,                    intent(in)  :: unit
  character(len=*),           intent(in)  :: file
  logical,                    intent(in)  :: formatted
  character(len=*),           intent(in)  :: mode
  logical,                    intent(out) :: err
  character(len=*), optional, intent(in)  :: status
  integer,          optional, intent(in)  :: recl

  ! type(unitSpec),   allocatable :: tmpUnits(:)
  character(len=:), allocatable :: fFileName, fStatus, fAction, fPosition, fMode
  character(len=256) :: tmpBoinc
  integer i, fRecl, nUnits, ioStat, chkUnit
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

  if(len_trim(file) > 64) then
    write(error_unit,"(a,i0)") "UNITS> ERROR Max length of file name in f_open is 64 characters, got ",len(file)
    call prror
  end if

  if(unit < units_minUnit .or. unit > units_maxUnit) then
    write(error_unit,"(3(a,i0),a)") "UNITS> ERROR Unit ",unit," is out of range ",units_minUnit,":",units_maxUnit," in f_open"
    call prror
  end if

#ifdef BOINC
  call boincrf(file,tmpBoinc)
  fFileName = trim(tmpBoinc)
#else
  fFileName = trim(file)
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

  fMode = mode
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
    fMode     = "r"
    fAction   = "read"
    fPosition = "asis"
  end select

  call f_getUnit(trim(file),chkUnit)
  if(chkUnit > 0) then
    ! We already have that file name in the record
    if(chkUnit /= unit) then
      write(error_unit,"(a,i0)") "UNITS> ERROR File '"//trim(file)//"' has already been assigned to unit ",chkUnit
      call prror
    end if
    units_uList(unit)%open  = .true.
  else
    ! The file is opened with a fixed unit, so save the info
    units_uList(unit)%file  = trim(file)
    units_uList(unit)%mode  = trim(fMode)
    units_uList(unit)%taken = .true.
    units_uList(unit)%open  = .true.
    units_uList(unit)%fixed = .true.
  end if

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
    write(error_unit,"(a,i0)") "UNITS> File '"//trim(file)//"' reported iostat = ",ioStat
    call f_writeLog("OPEN","ERROR",unit,file)
  end if

  if(units_uList(unit)%fixed) then
    call f_writeLog("OPEN","FiXED",unit,file)
  else
    call f_writeLog("OPEN","ASSIGNED",unit,file)
  end if
  return

10 continue
  err = .true.
  write(error_unit,"(a)") "UNITS> File '"//trim(file)//"' reported an error"
  call f_writeLog("OPEN","ERROR",unit,file)

end subroutine f_open

subroutine f_close(unit)

  use, intrinsic :: iso_fortran_env, only : error_unit
  
  integer, intent(in) :: unit

  integer i
  logical isOpen

  if(unit < units_minUnit .or. unit > units_maxUnit) then
    write(error_unit,"(3(a,i0),a)") "UNITS> ERROR Unit ",unit," is out of range ",units_minUnit,":",units_maxUnit," in f_close"
    call prror
  end if

  inquire(unit=unit, opened=isOpen)
  if(isOpen) then
    flush(unit)
    close(unit)
    units_uList(unit)%open = .false.
    if(units_uList(unit)%taken) then
      call f_writeLog("CLOSE","CLOSED",unit,units_uList(unit)%file)
    else
      call f_writeLog("CLOSE","CLOSED",unit,"*** Unknown File ***")
    end if
  else
    call f_writeLog("CLOSE","NOTOPEN",unit,units_uList(unit)%file)
  end if

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

  do i=units_minUnit,units_maxUnit
    inquire(unit=i, opened=isOpen)
    if(isOpen) flush(i)
  end do

end subroutine f_flush

subroutine f_writeLog(action,status,unit,file)

  use floatPrecision

  character(len=*), intent(in) :: action
  character(len=*), intent(in) :: status
  integer,          intent(in) :: unit
  character(len=*), intent(in) :: file

  real(kind=fPrec) cpuTime
  character(len=10) wAction
  character(len=8)  wStatus
  character(len=64) wFile
  
  if(units_logUnit <= 0) return ! Only write if we have a log file

  wAction = action
  wStatus = status
  wFile   = file

  call cpu_time(cpuTime)
  write(units_logUnit,"(f10.3,2x,a10,2x,a8,2x,i4,2x,a64)") cpuTime,adjustl(wAction),adjustl(wStatus),unit,adjustl(wFile)
  flush(units_logUnit)

end subroutine f_writeLog

end module mod_units
