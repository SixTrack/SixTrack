! ================================================================================================ !
!  BOINC DUMMY API
! ~~~~~~~~~~~~~~~~~
!  Updated: 2019-05-14 (V.K. Berglyd Olsen)
!  Should match lib/boinc/api/boinc_api_fortran.cpp
!
!  Logicals are returned from BOINC as integer, where non-zero is true and zero is false.
! ================================================================================================ !

subroutine boinc_time_to_checkpoint(doChpoint)
  implicit none
  integer, intent(out) :: doChpoint
  doChpoint = 1 ! Always allow checkpointing when using the dummy API
end subroutine boinc_time_to_checkpoint

subroutine boinc_is_standalone(isStandalone)
  implicit none
  integer, intent(out) :; isStandalone
  isStandalone = 1 ! Always treat dummy API as in standalone mode
end subroutine boinc_is_standalone

subroutine boincrf(fileName, boincName)
  implicit none
  character(*),   intent(in)  :: fileName
  character(256), intent(out) :: boincName
  boincName = fileName
end subroutine boincrf

subroutine boinc_finish(exitStatus)
  implicit none
  integer, intent(in) :: exitStatus
end subroutine boinc_finish

subroutine boinc_fraction_done(fracDone)
  implicit none
  double precision, intent(in) :: fracDone
end subroutine boinc_fraction_done

subroutine boinc_get_fraction_done(fracDone)
  implicit none
  double precision, intent(out) :: fracDone
end subroutine boinc_get_fraction_done

subroutine boinc_init
  implicit none
end subroutine boinc_init

subroutine boinc_checkpoint_completed
  implicit none
end subroutine boinc_checkpoint_completed

subroutine boinc_begin_critical_section
  implicit none
end subroutine boinc_begin_critical_section

subroutine boinc_end_critical_section
  implicit none
end subroutine boinc_end_critical_section

subroutine boinc_zip(mode,zipfile,path,strl1,strl2)
  implicit none
  integer mode,strl1,strl2
  character(*) zipfile,path
end subroutine boinc_zip

subroutine boinc_zipitall
  implicit none
end subroutine boinc_zipitall

subroutine boinc_unzip
  implicit none
end subroutine boinc_unzip
