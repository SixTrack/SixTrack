! This lives in its own file to prevent circular dependencies
subroutine closeUnits

  use mod_units,   only : units_maxUnit, f_close
  use dynk,        only : dynk_closeFiles
  use bdex,        only : bdex_closeFiles
#ifdef HDF5
  use hdf5_output, only : h5_closeHDF5
#endif

  implicit none

  integer i
  logical isOpen

  ! Call specific close units routines for modules
  !   that assign units via input or uses other types of files.
  call dynk_closeFiles
  call bdex_closeFiles
#ifdef HDF5
  call h5_closeHDF5
#endif

  ! Then iterate through 1 to units_maxUnit
  do i=1, units_maxUnit
    ! Do not close the following units:
    if(i == 5 .or. i == 6 .or. i == 91 .or. i == 92 .or. i == 93) cycle
    inquire(unit=i, opened=isOpen)
    if(isOpen) call f_close(i)
  end do

end subroutine closeUnits
