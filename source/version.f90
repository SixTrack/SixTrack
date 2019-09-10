module mod_version
  ! Keep data type in sync with 'cr_version' and 'cr_moddate'
  character(len=8),  parameter :: version = "5.3.3"
  integer,           parameter :: numvers = 50303
  character(len=10), parameter :: moddate = "09.09.2019"
  character(len=40), parameter :: git_revision = SIXTRACK_GIT_REVISION ! Git hash set by CMake
end module mod_version
