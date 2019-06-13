module mod_version
  ! Keep data type in sync with 'cr_version' and 'cr_moddate'
  character(len=8),  parameter :: version = "5.2.10"
  integer,           parameter :: numvers = 50210
  character(len=10), parameter :: moddate = "13.06.2019"
  character(len=40), parameter :: git_revision = SIXTRACK_GIT_REVISION ! Git hash set by CMake
end module mod_version
