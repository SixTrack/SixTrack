module mod_units
  
  implicit none
  
  type, private :: unitSpec
    integer,            private :: unit
    character(len=256), private :: filename
    logical,            private :: formatted
    logical,            private :: boinc
    logical,            private :: fio
    integer,            private :: recl
  end type unitSpec
  
  type(unitSpec), allocatable, private :: units_uList(:)
  integer,                     private :: units_nList
  
contains

subroutine units_initUnits
  allocate(units_uList(10))
  units_nList = 0
end subroutine units_initUnits

subroutine units_openUnits(unit,fileName,formatted,boinc,fio,recl)
  
  implicit none
  
  integer,                    intent(in) :: unit
  character(len=*),           intent(in) :: fileName
  logical,                    intent(in) :: formatted
  logical,          optional, intent(in) :: boinc
  logical,          optional, intent(in) :: fio
  integer,          optional, intent(in) :: recl
  
  type(unitSpec),   allocatable :: tmpUnits(:)
  character(len=:), allocatable :: fFileName
  integer i, fRecl, nUnits
  logical fBoinc, fFio
  
  nUnits      = size(units_uList)
  units_nList = units_nList + 1
  if(units_nList > nUnits) then
    allocate(tmpUnits(units_nList + 10))
    tmpUnits(1:units_nList-1) = units_uList(1:units_nList-1)
    call move_alloc(tmpUnits,units_uList)
  end if
  
  if(present(boinc)) then
    fBoinc = boinc
  else
    fBoinc = .false.
  end if
  
  if(present(fio)) then
    fFio = fio
  else
    fFio = .false.
  end if
  
  if(present(recl)) then
    fRecl = recl
  else
    fRecl = 0
  end if
  
#ifdef BOINC
  if(fBoinc) call boincrf(fileName,fFileName)
#else
  fFileName = fileName
#endif
#ifndef FIO
  fFio = .false.
#endif
#ifndef NAGFOR
  fRecl = 0
#endif
  
  if(.not. formatted) fFio = .false.
  
  units_uList(units_nList)%unit      = unit
  units_uList(units_nList)%filename  = fileName
  units_uList(units_nList)%formatted = formatted
  units_uList(units_nList)%boinc     = fBoinc
  units_uList(units_nList)%fio       = fFio
  units_uList(units_nList)%recl      = fRecl
  
  if(formatted) then
    if(fRecl > 0) then
      if(fFio) then
        open(unit,file=fFileName,form="formatted",status="unknown",round="nearest",recl=fRecl)
      else
        open(unit,file=fFileName,form="formatted",status="unknown",recl=fRecl)
      end if
    else
      if(fFio) then
        open(unit,file=fFileName,form="formatted",status="unknown",round="nearest")
      else
        open(unit,file=fFileName,form="formatted",status="unknown")
      end if
    end if
  else
    open(unit,file=fFileName,form="unformatted",status="unknown")
  endif
  
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
