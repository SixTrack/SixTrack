module mod_units
  
  implicit none
  
  type, public :: unitSpec
    integer,            public :: unit     ! Unit number
    character(len=256), public :: filename ! Standard filename
    integer,            public :: form     ! 0 = unformatted, 1 = formatted
    logical,            public :: sixtrack ! Default open in standard SixTrack
    logical,            public :: sixda    ! Default open in differential algebra version
    logical,            public :: boinc    ! Open in BOINC if enabled
    integer,            public :: recl=0   ! RECL flag, if needed
  end type unitSpec
  
  type(unitSpec), private :: units_uList(200)
  integer,        private :: units_nList
  
contains

subroutine units_initUnits
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
  
  integer i, fRecl
  logical fBoinc, fFio
  
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
  if(fBoinc) call boincrf(unitList(i)%filename, fileName)
#endif
#ifndef FIO
  fFio = .false.
#endif
#ifndef NAGFOR
  fRecl = 0
#endif
  
  if(.not. formatted) fFio = .false.
  
  if(formatted) then
    if(fRecl > 0) then
      if(fFio) then
        open(unit,file=fileName,form="formatted",status="unknown",round="nearest",recl=fRecl)
      else
        open(unit,file=fileName,form="formatted",status="unknown",recl=fRecl)
      end if
    else
      if(fFio) then
        open(unit,file=fileName,form="formatted",status="unknown",round="nearest")
      else
        open(unit,file=fileName,form="formatted",status="unknown")
      end if
    end if
  else
    open(unit,file=fileName,form="unformatted",status="unknown")
  endif
  
end subroutine units_openUnits

end module mod_units
