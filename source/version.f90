module mod_version
  ! Keep data type in sync with 'cr_version' and 'cr_moddate'
  character(len=8),  parameter :: version = "5.5.0"
  integer,           parameter :: numvers = 50500
  character(len=10), parameter :: moddate = "2022-04-04"
  character(len=40), parameter :: git_revision = SIXTRACK_GIT_REVISION ! Git hash set by CMake
end module mod_version
