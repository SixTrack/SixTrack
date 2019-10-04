program splitSingletrack
  ! ---------------------------------------------------------------
  ! Read a singletrackfile.dat containing multiple particle pairs,
  ! and split it into one file per pair.
  ! K. Sjobak, October 2019
  ! ---------------------------------------------------------------
  use mod_splitSingletrack
  implicit none

  logical                     :: hasInputFile
  character(len=*), parameter :: ifname = 'singletrackfile.dat'

  logical                     :: oldnames = .false.
  logical                     :: getNumPairs = .false.

  integer                     :: cmdarg_i, cmdarg_length, cmdarg_status
  character(len=100)          :: cmdarg_arg

  integer                     :: numPairs

  ! Parse command line arguments
  cmdarg_i = 0
  do
    call get_command_argument(cmdarg_i, cmdarg_arg, cmdarg_length, cmdarg_status)
    if (len_trim(cmdarg_arg) == 0) exit

    if (cmdarg_i .eq. 0) then
      !Skip first argument (command name);
      ! do nothing
    elseif (cmdarg_arg .eq. "--oldnames") then
      oldnames = .true.
    elseif (cmdarg_arg .eq. "--getNumPairs") then
      getnumpairs = .true.
    else if (cmdarg_arg .eq. "-h" .or. cmdarg_arg .eq. "--help") then
      write(*,'(a)') "Usage of 'splitsingletrack'; possible flags"
      write(*,'(a)') "Note that all flags are mutually exclusive."
      write(*,'(a)') "The input file should always be named '"//ifname//"'."
      write(*,'(a)') ""
      write(*,'(a)') "--oldnames    : Split into files named fort.91-i, where 1<=i<=32 is the pair index."
      write(*,'(a)') "--getNumPairs : Print the number of pairs and exit."
      write(*,'(a)') "--h / --help  : Print this help text."
      write(*,'(a)') ""
      write(*,'(a)') "By default it will try to split the file into files on the format '"// &
           ifname//".i', where i>=1 is the pair index"
      stop 0
    else
      write(*,'(a)') "Unrecognized command line argument '" // &
           trim(cmdarg_arg) // "', try -h"
      stop 2
    endif

    cmdarg_i = cmdarg_i + 1

  enddo

  !Check that input file is there
  inquire(file=ifname, exist=hasInputFile)
  if (.not. hasInputFile) then
    write(*, '(a)') "Error in splitSingletrack -- file '"//ifname//"' not found."
    stop 1
  endif

  if (getNumPairs) then
    numPairs = numSTFpairs(ifname)
    write(*,'(i8)') numPairs
    stop !No exit code so that we don't print "STOP 0" etc.
  else
    write(*,'(a)') "Splitting file '"//ifname//"' (run with '-h' for more options)..."
    call convertSTFfile(ifname, oldnames)
    write(*,'(a)') "Done!"
  endif

end program splitSingletrack
