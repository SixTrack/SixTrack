module mod_version
  ! Keep data type in sync with 'cr_version' and 'cr_moddate'
  character(len=8),  parameter :: version = "5.1.0"
  integer,           parameter :: numvers = 50100
  character(len=10), parameter :: moddate = "11.12.2018"
  character(len=40), parameter :: git_revision = SIXTRACK_GIT_REVISION ! Git hash set by CMake
end module mod_version
