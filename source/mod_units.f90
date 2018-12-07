module mod_units

  implicit none

  type, private :: unitSpec
    integer,            private :: unit
    character(len=256), private :: filename
    logical,            private :: formatted
    character(len=2),   private :: mode
    integer,            private :: recl
    logical,            private :: open
  end type unitSpec

  type(unitSpec), allocatable, private :: units_uList(:)
  integer,                     private :: units_nList

contains

subroutine units_initUnits
  allocate(units_uList(10))
  units_nList = 0
end subroutine units_initUnits

subroutine units_openUnit(unit,fileName,formatted,mode,err,status,recl)

  use crcoall

  implicit none

  integer,                    intent(in)    :: unit
  character(len=*),           intent(in)    :: fileName
  logical,                    intent(in)    :: formatted
  character(len=*),           intent(in)    :: mode
  logical,                    intent(inout) :: err
  character(len=*), optional, intent(in)    :: status
  integer,          optional, intent(in)    :: recl

  type(unitSpec),   allocatable :: tmpUnits(:)
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

  nUnits      = size(units_uList)
  units_nList = units_nList + 1
  if(units_nList > nUnits) then
    allocate(tmpUnits(units_nList + 10))
    tmpUnits(1:units_nList-1) = units_uList(1:units_nList-1)
    call move_alloc(tmpUnits,units_uList)
  end if

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
  call boincrf(fileName,tmpBoinc)
  fFileName = trim(tmpBoinc)
#else
  fFileName = fileName
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

  units_uList(units_nList)%unit      = unit
  units_uList(units_nList)%filename  = fileName
  units_uList(units_nList)%formatted = formatted
  units_uList(units_nList)%mode      = mode
  units_uList(units_nList)%recl      = fRecl
  units_uList(units_nList)%open      = .true.

  if(formatted) then
    if(fRecl > 0) then
      if(fFio) then
        open(unit,file=fFileName,form="formatted",status=fStatus,iostat=ioStat,&
          action=fAction,position=fPosition,round="nearest",recl=fRecl)
      else
        open(unit,file=fFileName,form="formatted",status=fStatus,iostat=ioStat,&
          action=fAction,position=fPosition,recl=fRecl)
      end if
    else
      if(fFio) then
        open(unit,file=fFileName,form="formatted",status=fStatus,iostat=ioStat,&
          action=fAction,position=fPosition,round="nearest")
      else
        open(unit,file=fFileName,form="formatted",status=fStatus,iostat=ioStat,&
          action=fAction,position=fPosition)
      end if
    end if
  else
    open(unit,file=fFileName,form="unformatted",status=fStatus,iostat=ioStat,&
      action=fAction,position=fPosition)
  endif

  if(ioStat /= 0) then
    err = .true.
    write(lout,"(a,i0)") "UNITS> File '"//trim(fFileName)//"' reported iostat = ",ioStat
  else
    err = .false.
  end if
  return

end subroutine units_openUnit

subroutine units_closeUnits(unit)

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

end subroutine units_closeUnits

subroutine units_flushUnits(unit)

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

end subroutine units_flushUnits

logical function units_isReserved(nUnit)

  integer, intent(in) :: nUnit
  
  integer i
  
  units_isReserved = .false.
  do i=1,units_nList
    if(units_uList(i)%unit == nUnit) then
      units_isReserved = .true.
      return
    end if
  end do

end function units_isReserved

end module mod_units
