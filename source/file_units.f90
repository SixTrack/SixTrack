! =================================================================================================
!  FILE UNITS MODULE
!  V.K. Berglyd Olsen, BE-ABP-HSS, March 2018
!  Last modified: 2018-04-19
! =================================================================================================
module file_units

  use strings

  implicit none

  integer, public :: funit_minUnit, funit_maxUnit
  parameter(funit_minUnit=1000, funit_maxUnit=1999)

  integer,      allocatable, private, save :: funit_usedUnits(:)
  type(string), allocatable, private, save :: funit_usedByFile(:)
  logical,                   private, save :: funit_arrInit  = .true.
  integer,                   private, save :: funit_nUnits   = 0
  integer,                   private, save :: funit_nextUnit = funit_minUnit

contains

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

10  continue

    fileUnit       = funit_nextUnit
    funit_nextUnit = funit_nextUnit + 1

    if(funit_nextUnit > funit_maxUnit) goto 30

    if(funit_arrInit) then
      allocate(funit_usedUnits(1))
      allocate(funit_usedByFile(1))
      funit_arrInit = .false.
    else
      allocate(tmp_usedUnits(funit_nUnits))
      allocate(tmp_usedByFile(funit_nUnits))
      tmp_usedUnits(1:funit_nUnits)  = funit_usedUnits(1:funit_nUnits)
      tmp_usedByFile(1:funit_nUnits) = funit_usedByFile(1:funit_nUnits)
      deallocate(funit_usedUnits)
      deallocate(funit_usedByFile)
      allocate(funit_usedUnits(funit_nUnits+1))
      allocate(funit_usedByFile(funit_nUnits+1))
      funit_usedUnits(1:funit_nUnits)  = tmp_usedUnits(1:funit_nUnits)
      funit_usedByFile(1:funit_nUnits) = tmp_usedByFile(1:funit_nUnits)
      deallocate(tmp_usedUnits)
      deallocate(tmp_usedByFile)
    end if
    funit_nUnits = funit_nUnits + 1

    inquire(unit=fileUnit, opened=isOpen)
    if(isOpen) then
      funit_usedUnits(funit_nUnits)  = fileUnit
      funit_usedByFile(funit_nUnits) = "Unknown, unit already open."
      goto 10
    else
      funit_usedUnits(funit_nUnits)  = fileUnit
      funit_usedByFile(funit_nUnits) = cleanName
      goto 20
    end if

20  return

30  continue
    write(lout,"(a)") "FUNIT> ERROR Failed to find an available file unit for file '"//cleanName//"'"
    stop -1

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

  ! Writes all assigned file units to file_units.dat
  subroutine funit_dumpUnits

    use crcoall

    implicit none

    integer dumpUnit, i
    logical isOpen

    call funit_requestUnit("file_units.dat",dumpUnit)
    inquire(unit=dumpUnit, opened=isOpen)
    if(isOpen) then
      write(lout,*) "ERROR in FUNIT when opening file_units.dat"
      write(lout,*) "Unit ",dumpUnit," was already taken."
      stop -1
    end if

    open(dumpUnit,file="file_units.dat",form="formatted")
    write(dumpUnit,"(a)") "# unit  assigned_to"
    do i=1, size(funit_usedUnits)
      write(dumpUnit,"(i6,a)") funit_usedUnits(i),"  "//funit_usedByFile(i)
    end do
    close(dumpUnit)

  end subroutine funit_dumpUnits

  ! Closes all units opened by the module
  subroutine funit_closeUnits

    implicit none

    integer chkUnit
    logical isOpen

    do chkUnit=funit_minUnit, funit_nextUnit-1
      inquire(unit=chkUnit, opened=isOpen)
      if(isOpen) close(chkUnit)
    end do

  end subroutine funit_closeUnits

end module file_units
! =================================================================================================
