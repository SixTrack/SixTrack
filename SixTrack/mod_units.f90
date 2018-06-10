module mod_units
  
  implicit none
  
  type, private :: unitSpec
    integer,            private :: unit
    character(len=256), private :: filename
    logical,            private :: formatted
    character(len=2),   private :: mode
    integer,            private :: recl
  end type unitSpec
  
  type(unitSpec), allocatable, private :: units_uList(:)
  integer,                     private :: units_nList
  
contains

subroutine units_initUnits
  allocate(units_uList(10))
  units_nList = 0
end subroutine units_initUnits

subroutine units_openUnits(unit,fileName,formatted,mode,err,status,recl)
  
  implicit none
  
  integer,                    intent(in)    :: unit
  character(len=*),           intent(in)    :: fileName
  logical,                    intent(in)    :: formatted
  character(len=*),           intent(in)    :: mode
  logical,                    intent(inout) :: err
  character(len=*), optional, intent(in)    :: status
  integer,          optional, intent(in)    :: recl
  
  type(unitSpec),   allocatable :: tmpUnits(:)
  character(len=:), allocatable :: fFileName, fStatus
  character(len=256) :: tmpBoinc
  integer i, fRecl, nUnits
  logical fFio
  
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
  
  units_uList(units_nList)%unit      = unit
  units_uList(units_nList)%filename  = fileName
  units_uList(units_nList)%formatted = formatted
  units_uList(units_nList)%mode      = mode
  units_uList(units_nList)%recl      = fRecl
  
  if(formatted) then
    if(fRecl > 0) then
      if(fFio) then
        open(unit,file=fFileName,form="formatted",status=fStatus,err=10,round="nearest",recl=fRecl)
      else
        open(unit,file=fFileName,form="formatted",status=fStatus,err=10,recl=fRecl)
      end if
    else
      if(fFio) then
        open(unit,file=fFileName,form="formatted",status=fStatus,err=10,round="nearest")
      else
        open(unit,file=fFileName,form="formatted",status=fStatus,err=10)
      end if
    end if
  else
    open(unit,file=fFileName,form="unformatted",status=fStatus,err=10)
  endif
  
  err = .false.
  return
  
10 continue
  err = .true.
  
end subroutine units_openUnits

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

end module mod_units
