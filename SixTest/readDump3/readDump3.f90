program readDump3
  ! Program to read a dump file of type 2 (Binary version of format 2)
  ! and convert it to a type 2.

  ! For now, no explicit binary->decimal conversion,
  ! since DUMP in SixTrack is also not doing this correctly...

  use, intrinsic :: iso_fortran_env, only : output_unit
  implicit none

  integer :: stat
  logical :: fileExists
  
  INTEGER :: cmdarg_i, cmdarg_length, cmdarg_status
  CHARACTER(len=100) :: cmdarg_arg
  
  character(len=100) :: ifname
  character(len=100) :: ofname

  integer ID, nturn, ktrack
  double precision dcum, x,xp,y,yp,sigma,dEE
  
  ! Read command line arguments: Name of input file, name of output file.
  cmdarg_i = 0
  do
     call get_command_argument(cmdarg_i, cmdarg_arg,cmdarg_length,cmdarg_status)
     if (len_trim(cmdarg_arg)==0) EXIT

     if (cmdarg_i.eq.0) then
        !Skip argument 0, the command name
        cmdarg_i = cmdarg_i+1
        cycle
     endif

     if (cmdarg_status.eq.-1) then
        write(*,*) "Error: Command line argument ", cmdarg_i, " was truncated."
        flush(output_unit)
        stop 2
     endif
     
     if (cmdarg_i.eq.1) then
        ifname = cmdarg_arg
     else if(cmdarg_i.eq.2) then
        ofname = cmdarg_arg
     else
        write(*,*) "Expected 2 arguments: inputfile outputfile"
        flush(output_unit)
        stop 1
     endif
     
     cmdarg_i = cmdarg_i+1
     
  end do

  if (cmdarg_i.ne. 3) then
        write(*,*) "Expected 2 arguments: inputfile outputfile"
        flush(output_unit)
        stop 3
  endif
  ! Check that the input files is there
  INQUIRE(FILE=ifname,EXIST=fileExists)
  if (.not. fileExists) then
     write(*,'(a,a,a)') "Error in readDump3 - input file '"//trim(ifname)//"' was not found"
     flush(output_unit)
     stop 4
  endif

  ! Open files
  open(100, file=ifname,form='UNFORMATTED',status='OLD')
  open(101, file=ofname,form='FORMATTED',status='REPLACE')

  ! Convert (no header)
  do
     read(100,end=5000,err=6000) ID,nturn,dcum, x,xp,y,yp,sigma,dEE, ktrack
     write(101,1986)             ID,nturn,dcum, x,xp,y,yp,sigma,dEE, ktrack
  end do
  
  !End of file
5000 continue
  close(100)
  close(101)

  ! Done :)
  stop

  !Error while reading
6000 continue
  write(*,*) "Error while reading file", ifname
  flush(output_unit)
  stop 5
  
  !Copied from dump_beam_population in sixtrack.s
1985 format (2(1x,I8),1X,F12.5,6(1X,1PE25.18),1X,I8)  !fmt 2 / hiprec
1986 format (2(1x,I8),1X,F12.5,6(1X,1PE16.9),1X,I8)   !fmt 2 / not hiprec
  
end program readDump3
