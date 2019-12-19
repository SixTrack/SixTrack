! ================================================================================================ !
!  FILE UNITS MODULE
! ===================
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2018-12-13
!  Updated: 2019-05-02
!
!  Module for keeping track of opened file units, their file names, and open the files correctly
!  depending on build flags like BOINC or FIO.
!
!  Note: The maximum number of allowed open units is 256. To ensure we cannot run out of units, the
!  module sets the last assigned unit number to 250. The log file is at 251. this leaves 4 units for
!  debug use in the code, if absolutely necessary.
! ================================================================================================ !
module mod_units

  use parpro, only : mPathName

  implicit none

  logical, public,  save       :: units_beQuiet  = .true.            ! No writing to lout when this flag is .true.

  ! Keep track of units
  integer, parameter           :: units_minUnit  = 1                 ! First unit to keep track of
  integer, parameter           :: units_maxUnit  = 250               ! Last unit to keep track of
  integer, parameter           :: units_minAuto  = 100               ! First unit available for dynamic allocation
  integer, private, save       :: units_nextUnit = units_minAuto     ! Next unit available for dynamic allocation
  integer, private, save       :: units_logUnit  = units_maxUnit+1   ! File unit for internal log file
  character(len=14), parameter :: units_logFile  = "file_units.log"  ! File name for internal log file

  type, private :: unitRecord
    character(len=:), allocatable, private :: file            ! The requested file name (not BOINC)
    character(len=3),              private :: mode  = " "     ! Read/write mode
    logical,                       private :: taken = .false. ! Whether a unit is known to be taken or not
    logical,                       private :: open  = .false. ! Whether file is opened by the module or not
    logical,                       private :: fixed = .true.  ! Whether the unit was requested as a fixed unit or not
  end type unitRecord

  ! Array to keep track of files
  type(unitRecord), private, save :: units_uList(units_minUnit:units_maxUnit+1)

  private :: f_writeLog

contains

! ================================================================================================ !
!  Initialise Units Module
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-12-13
!  This subroutine opens the log file at the highest unit. Make sure this unit is free!
! ================================================================================================ !
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
  write(units_logUnit,"(a)") "#   AtTime  Action    Unit  Status    FileName"
  flush(units_logUnit)
  call f_writeLog("INIT",units_logUnit,"FIXED",units_logFile)

end subroutine f_initUnits

! ================================================================================================ !
!  Request New File Units
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2018-12-13
!  Updated: 2018-12-13
!  Send in a file name, and get a unit back. If the has already been assigned a unit, this is
!  returned. Otherwise, a new unit is selected.
! ================================================================================================ !
subroutine f_requestUnit(file,unit)

  use crcoall

  character(len=*), intent(in)  :: file
  integer,          intent(out) :: unit

  integer i
  logical isOpen

  if(len_trim(file) > mPathName) then
    write(lerr,"(2(a,i0))") "UNITS> ERROR Max length of file path in f_requestUnit is ",mPathName,&
      " characters, got ",len_trim(file)
    call prror
  end if

  unit = -1
  call f_getUnit(file,unit)
  if(unit > 0) then
    call f_writeLog("REQUEST",unit,"EXISTS",trim(file))
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
        call f_writeLog("REQUEST",i,"TAKEN",trim(file))
      end if
    else
      unit = i
      exit
    end if
  end do

  if(unit > 0) then
    call f_writeLog("REQUEST",unit,"NEW",trim(file))
    units_uList(unit)%file  = trim(file)
    units_uList(unit)%mode  = ""
    units_uList(unit)%taken = .true.
    units_uList(unit)%open  = .false.
    units_uList(unit)%fixed = .false.
    units_nextUnit = unit + 1
  else
    write(lerr,"(a,i0)") "UNITS> ERROR Could not find an available file unit within the allowed range."
    call prror
  end if

end subroutine f_requestUnit

! ================================================================================================ !
!  Get Existing File Units
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2018-12-13
!  Updated: 2019-05-02
!  Will search through the record for a filename, and return its unit. -1 if it is not assigned.
! ================================================================================================ !
subroutine f_getUnit(file,unit)

  character(len=*), intent(in)  :: file
  integer,          intent(out) :: unit

  integer i

  unit = -1
  do i=units_minUnit,units_maxUnit
    if(allocated(units_uList(i)%file)) then
      if(units_uList(i)%file == trim(file)) then
        unit = i
        exit
      end if
    end if
  end do

end subroutine f_getUnit

! ================================================================================================ !
!  Open a File
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2018-12-13
!  Updated: 2019-05-02
!  This is a wrapper for Fortran open that also handles all the various build options.
!  The parameters are:
!   - unit      :: The unit number. Either a previously assigned one, or a fixed unit
!   - file      :: The file name. The file name will be checked to ensure that the unit value
!                  matches previous assigned or used unit.
!   - formatted :: .true. for formatted file, .false. for unformatted
!   - mode      :: Short form for action and position keywords.
!                  r, w,and rw correspond to read, write, and readwrite
!                  - corresponds to rewind, + to append, and default is asis
!   - err       :: Optional: If ommitted, errors will cause a call to prror.
!   - status    :: Optional: File status. Defaults to unknown
!   - recl      :: Optional: Record length. Is only used for nagfor
! ================================================================================================ !
subroutine f_open(unit,file,formatted,mode,err,iostat,status,access,recl)

  use crcoall

  integer,                    intent(in)  :: unit
  character(len=*),           intent(in)  :: file
  logical,                    intent(in)  :: formatted
  character(len=*),           intent(in)  :: mode
  logical,          optional, intent(out) :: err
  integer,          optional, intent(out) :: iostat
  character(len=*), optional, intent(in)  :: status
  character(len=*), optional, intent(in)  :: access
  integer,          optional, intent(in)  :: recl

  character(len=:), allocatable :: fFileName, fStatus, fAction, fPosition, fMode, fAccess
#ifdef BOINC
  character(len=mPathName+1) :: tmpBoinc
#endif
  integer fRecl, fStat, chkUnit
  logical fFio

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

  if(present(access)) then
    fAccess = access
  else
    fAccess = "sequential"
  end if

  if(len_trim(file) > mPathName) then
    write(lerr,"(2(a,i0))") "UNITS> ERROR Max length of file path in f_open is ",mPathName," characters, got ",len_trim(file)
    call prror
  end if

  if(unit < units_minUnit .or. unit > units_maxUnit) then
    write(lerr,"(3(a,i0),a)") "UNITS> ERROR Unit ",unit," is out of range ",units_minUnit,":",units_maxUnit," in f_open"
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
    write(lout,"(a)") "UNITS> WARNING Unknown file mode '"//trim(fMode)//"', defaulting to 'r'"
    fMode     = "r"
    fAction   = "read"
    fPosition = "asis"
  end select

  ! Check that we don't have a fixed vs. assigned unit crash
  if(allocated(units_uList(unit)%file)) then
    if(units_uList(unit)%file /= file .and. units_uList(unit)%file /= " ") then
      write(lerr,"(a,i0,a)") "UNITS> ERROR Unit ",unit," has already been assigned to file '"//&
        trim(units_uList(unit)%file)//"' /= '"//trim(file)//"'"
      call prror
    end if
  end if

  call f_getUnit(trim(file),chkUnit)
  if(chkUnit > 0) then
    ! We already have that file name in the record
    if(chkUnit /= unit) then
      write(lerr,"(a,i0)") "UNITS> ERROR File '"//trim(file)//"' has already been assigned to unit ",chkUnit
      call prror
    end if
  else
    ! The file is opened with a fixed unit, so save the info
    units_uList(unit)%file  = trim(file)
    units_uList(unit)%mode  = trim(fMode)
    units_uList(unit)%taken = .true.
    units_uList(unit)%open  = .true.
    units_uList(unit)%fixed = .true.
  end if

  if(formatted) then
    if(fRecl > 0) then
      if(fFio) then
        open(unit,file=fFileName,form="formatted",status=fStatus,iostat=fStat,&
          action=fAction,access=fAccess,position=fPosition,round="nearest",recl=fRecl,err=10)
      else
        open(unit,file=fFileName,form="formatted",status=fStatus,iostat=fStat,&
          action=fAction,access=fAccess,position=fPosition,recl=fRecl,err=10)
      end if
    else
      if(fFio) then
        open(unit,file=fFileName,form="formatted",status=fStatus,iostat=fStat,&
          action=fAction,access=fAccess,position=fPosition,round="nearest",err=10)
      else
        open(unit,file=fFileName,form="formatted",status=fStatus,iostat=fStat,&
          action=fAction,access=fAccess,position=fPosition,err=10)
      end if
    end if
  else
    open(unit,file=fFileName,form="unformatted",status=fStatus,iostat=fStat,&
      action=fAction,access=fAccess,position=fPosition,err=10)
  endif

  if(present(iostat)) then
    iostat = fStat
  end if

  if(fStat /= 0) then
    call f_writeLog("OPEN",unit,"ERROR",file)
    if(present(err)) then
      err = .true.
      if(units_beQuiet .eqv. .false.) then
        write(lout,"(a,i0)") "UNITS> File '"//trim(file)//"' reported iostat = ",fStat
      end if
    else
      write(lerr,"(a,i0)") "UNITS> ERROR File '"//trim(file)//"' reported iostat = ",fStat
      call prror
    end if
  end if

  if(units_uList(unit)%fixed) then
    call f_writeLog("OPEN",unit,"FIXED",file)
  else
    call f_writeLog("OPEN",unit,"ASSIGNED",file)
  end if
  if(present(err)) then
    err = .false.
  end if

  return

10 continue
  call f_writeLog("OPEN",unit,"ERROR",file)
  if(present(iostat)) then
    iostat = fStat
  end if
  if(present(err)) then
    err = .true.
    if(units_beQuiet .eqv. .false.) then
      write(lout,"(a)") "UNITS> Could not open '"//trim(file)//"'"
    end if
  else
    write(lerr,"(a)") "UNITS> ERROR Could not open '"//trim(file)//"'"
    call prror
  end if

end subroutine f_open

! ================================================================================================ !
!  Close File Units
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2018-12-13
!  Updated: 2019-10-10
!  Preferred method for closing file as it keeps the record up to date
! ================================================================================================ !
subroutine f_close(unit, cErr)

  use crcoall

  integer,           intent(in)  :: unit
  logical, optional, intent(out) :: cErr

  integer ioStat
  logical isOpen

  if(unit < units_minUnit .or. unit > units_maxUnit) then
    write(lerr,"(3(a,i0),a)") "UNITS> ERROR Unit ",unit," is out of range ",units_minUnit,":",units_maxUnit," in f_close"
    call prror
  end if

  if(present(cErr)) then
    cErr = .false.
  end if

  inquire(unit=unit, opened=isOpen)
  if(isOpen) then
    if(units_uList(unit)%taken) then
      call f_writeLog("CLOSE",unit,"CLOSED",units_uList(unit)%file)
    else
      call f_writeLog("CLOSE",unit,"CLOSED","*** Unknown File ***")
    end if
    flush(unit)
    close(unit, iostat=ioStat)
    if(ioStat == 0) then
      units_uList(unit)%open = .false.
    else
      if(units_uList(unit)%taken) then
        write(lerr,"(a,i0)") "UNITS> ERROR Failed to close file '"//trim(units_uList(unit)%file)//"', status code ",ioStat
        call f_writeLog("CLOSE",unit,"ERROR",units_uList(unit)%file)
      else
        write(lerr,"(2(a,i0))") "UNITS> ERROR Failed to close file unit ",unit,", status code ",ioStat
        call f_writeLog("CLOSE",unit,"ERROR","*** Unknown File ***")
      end if
      if(present(cErr)) then
        cErr = .true.
      end if
    end if
  else
    if(units_uList(unit)%taken) then
      call f_writeLog("CLOSE",unit,"NOTOPEN",units_uList(unit)%file)
    else
      call f_writeLog("CLOSE",unit,"NOTOPEN","*** Unknown File ***")
    end if
  end if

end subroutine f_close

! ================================================================================================ !
!  Free File Units
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-05-02
!  Updated: 2019-05-02
!  Free a file unit so that it can be assigned again. This implies close.
! ================================================================================================ !
subroutine f_freeUnit(unit)

  use crcoall

  integer, intent(in) :: unit

  logical isOpen

  if(unit < units_minUnit .or. unit > units_maxUnit) then
    write(lerr,"(3(a,i0),a)") "UNITS> ERROR Unit ",unit," is out of range ",units_minUnit,":",units_maxUnit," in f_close"
    call prror
  end if

  inquire(unit=unit, opened=isOpen)
  if(isOpen) then
    call f_close(unit)
  end if

  call f_writeLog("FREE",unit,"FREED",units_uList(unit)%file)

  units_uList(unit)%file  = " "
  units_uList(unit)%mode  = " "
  units_uList(unit)%taken = .false.
  units_uList(unit)%open  = .false.
  units_uList(unit)%fixed = .true.

end subroutine f_freeUnit

! ================================================================================================ !
!  Flush Single or All File Units
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-12-13
! ================================================================================================ !
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

! =================================================================================================
!  Reposition Text File
!
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-24
!  Updated: 2019-12-19
!
!  The routine was moved from Scatter Module, and was originally intended for Checkpoint/Restart,
!  but should also work outside of C/R.
!
!  Parameters:
!   * fileName : The file name of the file to be opened. Needed both to ensure not multiple units
!                is assigned, and to make log messages more readable.
!   * fileUnit : The file unit variable for the file. The unit can be assigned already, but must
!                not be open. This routine opens the file, and returns the unit.
!   * filePos  : The file position variable used to track the latest position in the file. This is
!                the variable the calling code must keep track on and checkpoint.
!   * newPos   : The desired position the file should be truncated to. This is the value the
!                calling code should restore from checkpoint file. It is the checkpointed value of
!                filePos.
! =================================================================================================
subroutine f_positionFile(fileName, fileUnit, filePos, newPos)

  use parpro
  use crcoall

  character(len=*), intent(in)    :: fileName
  integer,          intent(inout) :: fileUnit
  integer,          intent(inout) :: filePos
  integer,          intent(in)    :: newPos

  logical isOpen, fErr
  integer iError
  integer j
  character(len=mInputLn) aRecord

  call f_requestUnit(fileName,fileUnit)
  inquire(unit=fileUnit,opened=isOpen)
  if(isOpen) then
#ifdef CR
    write(crlog,"(a,i0,a)") "UNITS> ERROR Failed while repositioning '"//trim(fileName)//"', "//&
      "unit ",fileUnit," already in use"
    flush(crlog)
#endif
    write(lerr, "(a,i0,a)") "UNITS> ERROR Failed while repositioning '"//trim(fileName)//"', "//&
      "unit ",fileUnit," already in use"
    call prror
  end if

  if(newPos /= -1) then
    call f_open(unit=fileUnit,file=fileName,formatted=.true.,mode="rw",err=fErr,iostat=iError,status="old")
    if(fErr) goto 10
    filePos = 0
    do j=1,newPos
      read(fileUnit,"(a)",end=10,err=10,iostat=iError) aRecord
      filePos = filePos + 1
    end do
    ! Truncate file at requested position
    endfile(fileUnit,iostat=iError)
    call f_close(fileUnit)

    ! Re-open for appending
    call f_open(unit=fileUnit,file=fileName,formatted=.true.,mode="w+",err=fErr,iostat=iError,status="old")
    if(fErr) goto 10
#ifdef CR
    write(crlog,"(2(a,i0))") "UNITS> Sucessfully repositioned '"//trim(fileName)//"': "//&
      "Position: ",filePos,", Position C/R: ",newPos
#endif
    write(lout, "(2(a,i0))") "UNITS> Sucessfully repositioned '"//trim(fileName)//"': "//&
      "Position: ",filePos,", Position C/R: ",newPos
  else
#ifdef CR
    write(crlog,"(a,i0)") "UNITS> Did not attempt repositioning '"//trim(fileName)//"' at line ",newPos
    write(crlog,"(a)")    "UNITS> If anything has been written to the file, it will be correctly truncated"
#endif
    write(lout, "(a,i0)") "UNITS> Did not attempt repositioning '"//trim(fileName)//"' at line ",newPos
    write(lout, "(a)")    "UNITS> If anything has been written to the file, it will be correctly truncated"
  end if

#ifdef CR
  flush(crlog)
#endif
  flush(lout)

  return

10 continue
#ifdef CR
  write(crlog,"(a,i0)")    "UNITS> ERROR While reading '"//trim(fileName)//"', iostat = ",iError
  write(lerr, "(2(a,i0))") "UNITS>       Current position: ",filePos,", requested position: ",newPos
  flush(crlog)
#endif
  write(lerr, "(a,i0)")    "UNITS> ERROR While reading '"//trim(fileName)//"', iostat = ",iError
  write(lerr, "(2(a,i0))") "UNITS>       Current position: ",filePos,", requested position: ",newPos
  call prror

end subroutine f_positionFile

! ================================================================================================ !
!  Internal Log File Writer
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-12-13
! ================================================================================================ !
subroutine f_writeLog(action,unit,status,file)

  use parpro
  use floatPrecision

  character(len=*), intent(in) :: action
  integer,          intent(in) :: unit
  character(len=*), intent(in) :: status
  character(len=*), intent(in) :: file

  real(kind=fPrec)         cpuTime
  character(len=8)         wAction
  character(len=8)         wStatus
  character(len=mFileName) wFile
  integer                  cFile

  if(units_logUnit <= 0) return ! Only write if we have a log file open

  wAction = action
  wStatus = status

  cFile = len_trim(file)
  if(cFile > mFileName) then
    wFile = "[...]"//file(cFile-mFileName+5:cFile)
  else
    wFile = file
  end if

  call cpu_time(cpuTime)
  write(units_logUnit,"(f10.3,2x,a8,2x,i4,2x,a8,2x,a)") cpuTime,adjustl(wAction),unit,adjustl(wStatus),trim(wFile)
  flush(units_logUnit)

end subroutine f_writeLog

end module mod_units
