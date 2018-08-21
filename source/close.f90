! This lives in its own file to prevent circular dependencies
subroutine closeUnits

  use file_units,  only : funit_closeUnits, funit_minUnit
  use scatter,     only : scatter_closeUnits
  use dynk,        only : dynk_closeFiles
  use dump,        only : dump_closeUnits
  use bdex,        only : bdex_closeFiles
#ifdef HDF5
  use hdf5_output, only : h5_closeHDF5
#endif

  implicit none

  integer chkUnit
  logical isOpen

  ! Call specific close units routines for modules
  !   that assign units dynamically or via input.
  call funit_closeUnits
  call dump_closeUnits
  call dynk_closeFiles
  call bdex_closeFiles
  call scatter_closeUnits
#ifdef HDF5
  call h5_closeHDF5
#endif

  ! Then iterate through the first 1000 units
  do chkUnit=1, funit_minUnit-1
    ! Do not close C/R files. handled by abend
    if(chkUnit >= 91 .and. chkUnit <= 96) cycle
    inquire(unit=chkUnit, opened=isOpen)
    if(isOpen) close(chkUnit)
  end do

  return

end subroutine closeUnits
