! =================================================================================================
!  FILE UNITS MODULE
!  V.K. Berglyd Olsen, BE-ABP-HSS, March 2018
!  Last modified: 2018-04-19
! =================================================================================================
module file_units

  use strings

  implicit none

  integer, parameter :: funit_minUnit = 1000
  integer, parameter :: funit_maxUnit = 1999

  integer,      allocatable, private, save :: funit_usedUnits(:)
  type(string), allocatable, private, save :: funit_usedByFile(:)
  integer,                   private, save :: funit_nUnits   = 0
  integer,                   private, save :: funit_nextUnit = funit_minUnit
  integer,                   private, save :: funit_logUnit  = 6

contains

subroutine funit_initUnits

  allocate(funit_usedUnits(10))
  allocate(funit_usedByFile(10))

  funit_logUnit       = funit_minUnit
  funit_nextUnit      = funit_minUnit + 1
  funit_nUnits        = 1
  funit_usedUnits(1)  = funit_logUnit
  funit_usedByFile(1) = "file_units.dat"

  open(funit_logUnit,file="file_units.dat",form="formatted")
  write(funit_logUnit,"(a)")    "# unit  assigned_to"
  write(funit_logUnit,"(i6,a)") funit_usedUnits(1),"  "//funit_usedByFile(1)
  flush(funit_logUnit)

end subroutine funit_initUnits

! Request a new fileunit. The filename parameter is only for internal record keeping. It should
! be a string that describes what file, or set of files, the unit number is used for.
subroutine funit_requestUnit(fileName, fileUnit)

  use crcoall

  implicit none

  character(len=*), intent(in)  :: fileName
  integer,          intent(out) :: fileUnit
  character(len=:), allocatable :: cleanName

  integer,          allocatable :: tmp_usedUnits(:)
  type(string),     allocatable :: tmp_usedByFile(:)

  logical isOpen
  integer i

  cleanName = trim(fileName)

10 continue

  fileUnit       = funit_nextUnit
  funit_nextUnit = funit_nextUnit + 1
  if(funit_nextUnit > funit_maxUnit) goto 30

  if(funit_nUnits + 1 > size(funit_usedUnits)) then
    allocate(tmp_usedUnits(funit_nUnits  + 10))
    allocate(tmp_usedByFile(funit_nUnits + 10))
    tmp_usedUnits(1:funit_nUnits)  = funit_usedUnits(1:funit_nUnits)
    tmp_usedByFile(1:funit_nUnits) = funit_usedByFile(1:funit_nUnits)
    call move_alloc(tmp_usedUnits, funit_usedUnits)
    call move_alloc(tmp_usedByFile,funit_usedByFile)
  end if
  funit_nUnits = funit_nUnits + 1

  inquire(unit=fileUnit, opened=isOpen)
  if(isOpen) then
    funit_usedUnits(funit_nUnits)  = fileUnit
    funit_usedByFile(funit_nUnits) = "Unknown, unit already open."
    write(funit_logUnit,"(i6,a)") funit_usedUnits(funit_nUnits),"  "//funit_usedByFile(funit_nUnits)
    goto 10
  else
    funit_usedUnits(funit_nUnits)  = fileUnit
    funit_usedByFile(funit_nUnits) = cleanName
    write(funit_logUnit,"(i6,a)") funit_usedUnits(funit_nUnits),"  "//funit_usedByFile(funit_nUnits)
    return
  end if
  flush(funit_logUnit)

  return

30 continue
  write(lout,"(a)") "FUNIT> ERROR Failed to find an available file unit for file '"//cleanName//"'"
  call prror(-1)

end subroutine funit_requestUnit

! Lists all assigned file units to lout
subroutine funit_listUnits

  use crcoall

  implicit none

  integer i

  write(lout,"(a)") "FUNIT> Dynamically assigned file units are:"
  do i=1, size(funit_usedUnits)
    write(lout,"(a,i4,a)") "FUNIT>   Unit ",funit_usedUnits(i)," assigned to: "//funit_usedByFile(i)
  end do

end subroutine funit_listUnits

! Flushes all units opened by the module
subroutine funit_flushUnits

  implicit none

  integer chkUnit
  logical isOpen

  do chkUnit=funit_minUnit, funit_nextUnit-1
    inquire(unit=chkUnit, opened=isOpen)
    if(isOpen) flush(chkUnit)
  end do

end subroutine funit_flushUnits

! Closes all units opened by the module (also flushes)
subroutine funit_closeUnits

  implicit none

  integer chkUnit
  logical isOpen

  do chkUnit=funit_minUnit, funit_nextUnit-1
    inquire(unit=chkUnit, opened=isOpen)
    if(isOpen) then
      flush(chkUnit)
      close(chkUnit)
    end if
  end do

end subroutine funit_closeUnits

end module file_units
! =================================================================================================
