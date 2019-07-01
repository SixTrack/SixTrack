! =================================================================================================
!  END SIXTRACK
!  Last modified: 2018-11-30
! =================================================================================================

subroutine prror(ier)

  use crcoall
  use mod_common, only : errout_status
#ifdef FLUKA
  use mod_fluka
#endif

  implicit none

  integer, optional, intent(in) :: ier

#ifdef FLUKA
  integer fluka_con
#endif

  if(present(ier)) then
    errout_status = ier
  else
    errout_status = -1
  end if

  write(lout,"(a)") ""
  write(lout,"(a)") "    +++++++++++++++++++++++++++++"
  write(lout,"(a)") "    +      ERROR DETECTED!      +"
  write(lout,"(a)") "    + RUN TERMINATED ABNORMALLY +"
  write(lout,"(a)") "    +++++++++++++++++++++++++++++"
  write(lout,"(a)") ""

#ifdef FLUKA
  call fluka_close
#endif

#ifdef CR
  call abend('                                                  ')
#else
  call closeUnits
  stop 1
#endif

end subroutine prror

#ifdef CR
! =================================================================================================
subroutine abend(cstring)

  use numerical_constants

  use, intrinsic :: iso_fortran_env, only : error_unit

  use crcoall
  use parpro
  use mod_common
  use checkpoint_restart
  use string_tools
  use mod_units
  use mod_version
  use mod_time

  implicit none

  integer i,lstring,j,itot,ttot
  character(len=*) cstring
  character(len=256) filename
  real(kind=fPrec) sumda(60)
  logical fopen, rErr, fErr
  character(len=8192) ch
  character(len=25) ch1
  integer errno,l1,l2
  integer ich

  save

  write(93,*) 'SIXTRACR STOP/ABEND called and closing files'
  endfile(93,iostat=ierro)
  backspace(93,iostat=ierro)

  ! Calling close to be very safe.......96 calls to abend
  ! Easier than adding the call on every abend
  call closeUnits

  ! If fort.10 is inexistent (physics error or some other problem)
  ! we try and write a 0d0 file with a turn number and CPU time
  write(93,"(a)") 'SIXTRACR STOP/ABEND checking fort.10'
  endfile(93,iostat=ierro)
  backspace(93,iostat=ierro)

  call units_openUnit(unit=10,fileName="fort.10",formatted=.true.,mode="rw",err=fErr,recl=8195)
  if(fErr) goto 11

  ! Now we try and read fort.10 i.e. is it empty?
  read(10,'(a1024)',end=11,err=11,iostat=ierro) arecord
  ! Seems to be OK
  goto 12

11 continue
  ! Now we try and write a fort.10
  ! We put some CPU for Igor, a version, and turn number 0
  write(93,*) 'SIXTRACR STOP/ABEND writing a fort.10'
  endfile(93,iostat=ierro)
  backspace(93,iostat=ierro)

  ! The version
  itot=0
  ttot=0
  do i=1,8
    if (version(i:i).ne.' ') then
      if (version(i:i).ne.'.') then
        itot=itot*10+ichar(version(i:i))-ichar('0')
      else
        ttot=ttot*10**2+itot
        itot=0
      end if
    end if
  end do
  ttot=ttot*10**2+itot
  do i=1,60
    sumda(i)=zero
  end do
  sumda(52)=real(ttot,fPrec)
  ! The CPU
  call time_timerCheck(time1)
  trtime=time1-time0
#ifdef CR
  trtime=trtime+crtime3
#endif
  sumda(60)=real(trtime,fPrec)
  ! Note it COULD happen that napxo is 0 for a very very early error and even napx!!!
  if (napxo.eq.0.and.napx.eq.0) napxo=1
  write(93,*) 'SIXTRACR STOP/ABEND writing fort.10, lines',napxo,'/',napx
  endfile(93,iostat=ierro)
  backspace(93,iostat=ierro)
  if (napxo.eq.0.and.napx.eq.0) napxo=1
  if (napxo.eq.0) napxo=napx
  do j=1,napxo
#ifndef CRLIBM
    write(ch,*,iostat=ierro) (sumda(i),i=1,60)
    do ich=8192,1,-1
      if(ch(ich:ich).ne.' ') goto 707
    end do
707 write(10,'(a)',iostat=ierro) ch(:ich)
#else
    l1=1
    do i=1,60
      ! We return the length of the string (always 24)
      call chr_fromReal(sumda(i),ch1,19,2,rErr)
      ch(l1:l1+25)=' '//ch1(1:25)
      l1=l1+26
    end do
    write(10,'(a)',iostat=ierro) ch(1:l1-1)
#endif
    if(ierro.ne.0) then
      write(lout,"(a,i0)") "ABEND> ERROR Problems writing to fort.10. ierro = ",ierro
    end if
  end do

12 continue
  close (10,iostat=ierro)
#ifdef CR
#ifdef DEBUG
  ! call system('../crend   >> crlog')
#endif
  close(91,err=4)
4 continue
  close(94,err=5)
5 continue
  close(95,err=6)
6 continue
  close(96,err=7)
7 continue
  if (lout.eq.92) then
    write(93,*) 'SIXTRACR STOP/ABEND copying fort.92 to fort.6'
    endfile(93,iostat=ierro)
    backspace(93,iostat=ierro)
    rewind 92
3   read(92,'(a1024)',end=1,err=8,iostat=ierro) arecord
    lstring=1024
    do i=1024,2,-1
      lstring=i
      if (arecord(i:i).ne.' ') goto 2
      lstring=lstring-1
    end do
2   write(6,'(a)',iostat=ierro) arecord(1:lstring)
    goto 3
  end if

1 write(6,"(a)",iostat=ierro) 'SIXTRACR> Stop '//cstring
  close(6,iostat=ierro)
  ! and get rid of fort.92
  rewind 92
  endfile(92,iostat=ierro)
  close(92)
  write(93,*) 'SIXTRACR stop '//cstring
  write(93,*)
#ifdef DEBUG
  ! call system('../crend   >> crlog')
#endif
#ifdef BOINC
  do i=2,120
    inquire(i,opened=fopen)
    write(93,*) 'UNIT ',i,' opened ',fopen
  end do
  ! call boinc_zipitall()
  ! call boinc_finish_graphics()
  if(errout_status.ne.0) then
    close(93)
    call boincrf('fort.93',filename)
    call print_lastlines_to_stderr(93,filename)
    call boincrf('fort.6',filename)
    call print_lastlines_to_stderr(6,filename)
  end if
  call boinc_finish(errout_status) !This call does not return
#else
  if(errout_status.ne.0) then
    close(93)
    call print_lastlines_to_stderr(93,"fort.93")
    call print_lastlines_to_stderr(6,"fort.6")

    write(error_unit,'(a,i5)') "Stopping, errout_status=",errout_status
    stop 1
  else
    ! No error
    stop
  end if
#endif

!!!!!! In case of errors when copying fort.92 (lout) -> fort.6 !!!!!!!!!

8 write(93,*) 'SIXTRACR CR ABEND *** ERROR *** reading fort.92, iostat=',ierro
  close(93)
  write(6,*) 'SIXTRACR CR ABEND *** ERROR *** reading fort.92, iostat=',ierro
#ifdef DEBUG
  ! call system('../crend   >> crlog')
#endif
#ifdef BOINC
  do i=2,120
    inquire(i,opened=fopen)
    write(6,*) 'UNIT ',i,' opened ',fopen
  end do
  close(6,err=31)
31 continue
  ! call boinc_zipitall()
  ! call boinc_finish_graphics()
  if(errout_status.ne.0) then
    close(93)
    call boincrf('fort.93',filename)
    call print_lastlines_to_stderr(93,filename)
    call boincrf('fort.6',filename)
    call print_lastlines_to_stderr(6,filename)
  end if
  call boinc_finish(errout_status) !This call does not return
#else
  if(errout_status.ne.0) then
    close(93)
    call print_lastlines_to_stderr(93,"fort.93")
    call print_lastlines_to_stderr(6,"fort.6")

    write(error_unit,'(a,i5)') "Stopping, errout_status=",errout_status
    stop 1
  else
    ! No error
    stop
  end if
#endif
#else
  ! This one should probably remain as write(*,*) or use output_unit
  write(*,*) 'SIXTRACK STOP/ABEND '//cstring
#ifdef DEBUG
  ! call system('../crend   >> crlog')
#endif
  ! No fort.6 and 93 if not CR -> don't do print_lastlines_to_stderr()
  if(errout_status.ne.0) then
    write(error_unit,'(a,i5)') "Stopping, errout_status=",errout_status
    stop 1
  else
    ! No error
    stop
  end if
#endif

end subroutine abend
#endif

! =================================================================================================
!  K. Sjobak, BE-ABP-HSS
!  Last modified: June 2017
!  Subroutine to copy the last nlines lines from a file to stderr
!  It is mainly used just before exiting SixTrack in case there was an error.
!  This is useful since STDERR is often returned from batch systems and BOINC.
! =================================================================================================
subroutine print_lastlines_to_stderr(file_unit, file_name)

  use, intrinsic :: iso_fortran_env, only : error_unit
  use string_tools

  implicit none

  integer,          intent(in) :: file_unit
  character(len=*), intent(in) :: file_name

  integer, parameter :: nlines = 40

  character(len=1024) fileBuff (nlines)
  integer fileBuff_idx
  integer i,j
  integer printLines

  logical lopen
  integer ierro

  ! Clear the buffer
  do i=1,nlines
    fileBuff(i)=''
  end do

  ! Open the file
  inquire(unit=file_unit,opened=lopen)
  if(lopen) then
    write(error_unit,'(a,1x,i5,1x,a)') "Error when opening unit #",file_unit,": The file is already open."
    return
  end if

  open(file_unit,file=file_name,form="formatted",status="old",iostat=ierro)
  if(ierro .ne. 0) then
    write(error_unit,'(a,a,a,1x,i5,1x,a,1x,i5)') &
      "Error when opening file '",trim(file_name),"' on unit #", file_unit,", iostat =",ierro
    return
  end if

  ! Read into fileBuff as a ring buffer.
  fileBuff_idx = 1
  j = 0

1 read(file_unit,'(a1024)',end=3,err=2,iostat=ierro) fileBuff(fileBuff_idx)
  ! write(error_unit,*) fileBuff_idx,":",trim(fileBuff(fileBuff_idx))
  fileBuff_idx = fileBuff_idx+1
  if (fileBuff_idx.ge.nlines) fileBuff_idx = 1
  j = j+1
  goto 1

2 continue ! An error occured
  write(error_unit,'(a,1x,i5)') "An error occured while reading the file, iostat=",ierro


3 continue ! EOF or error; anyway we're done.
  close(file_unit)

  ! Print stuff back out from the buffer
  if (j .lt. nlines) then
    printLines = j
    fileBuff_idx=1
  else
    printLines = nlines
  end if
  write(error_unit,'(a,1x,i5,1x,a,a,a)') "******* Last",printLines,"lines of file '",trim(file_name),"': *******"

  i = fileBuff_idx  ! Position in buffer (we have already incremented i)
  j = 0             ! How many have we printed

10 continue
  if(i.gt.nlines) i=1      ! i is wrapping
    write(error_unit,'(a)') trim(fileBuff(i))
    i = i+1
    j = j+1                   ! j just counts
    if(j.lt.printLines) goto 10
      write(error_unit,'(a,a,a)') "******* Done writing tail of file '",trim(file_name),"' to stderr *******"

end subroutine print_lastlines_to_stderr
