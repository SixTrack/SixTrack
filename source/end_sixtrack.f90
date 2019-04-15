! =================================================================================================
!  END SIXTRACK
!  Last modified: 2018-06-25
! =================================================================================================

subroutine prror(ier)

  use crcoall
  use mod_common, only : errout
#ifdef FLUKA
  use mod_fluka
#endif

  implicit none

  integer, optional, intent(in) :: ier

#ifdef FLUKA
  integer fluka_con
#endif

  if(present(ier)) then
    errout = ier
  else
    errout = -1
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
! ================================================================================================ !
!  End Routine for Checkpoint/Restart Version
!  Changed: 2019-04-15
! ================================================================================================ !
subroutine abend(cstring)

  use, intrinsic :: iso_fortran_env, only : error_unit

  use crcoall
  use parpro
  use mod_common
  use checkpoint_restart
  use numerical_constants
  use string_tools
  use mod_units
  use mod_version
  use mod_time

  implicit none

  character(len=*), intent(in) :: cstring

  real(kind=fPrec) sumda(60)
  integer i, j, k
  logical fOpen, rErr, fErr, fExists
  character(len=25) chBuf
  character(len=mInputLn) inLine, outLine

  write(93,"(a)") "SIXTRACR> STOP/ABEND called and closing files"
  flush(93)

  ! Calling close to be very safe.......96 calls to abend
  ! Easier than adding the call on every abend
  call closeUnits

#ifdef BOINC
  ! If fort.10 is non-existent (physics error or some other problem)
  ! we try and write a 0d0 file with a turn number and CPU time
  write(93,"(a)") "SIXTRACR> STOP/ABEND checking fort.10"
  flush(93)

  call f_open(unit=10,file="fort.10",formatted=.true.,mode="r",err=fErr,status="old",recl=8195)
  if(fErr) goto 11

  ! Now we try and read fort.10 i.e. is it empty?
  read(10,"(a1024)",end=11,err=11,iostat=ierro) inLine
  ! Seems to be OK
  goto 12

11 continue
  ! Now we try and write a fort.10
  ! We put some CPU for Igor, a version, and turn number 0
  ! call f_open(unit=10,file="fort.10",formatted=.true.,mode="w",err=fErr,status="unknown",recl=8195)
  write(93,"(a)") "SIXTRACR> STOP/ABEND writing dummy fort.10"
  flush(93)

  ! Make sure it is closed properly before we re-open for dummy write
  ! inquire(10,opened=fOpen)
  ! if(fOpen) close(10)
  call f_close(10)
  call f_open(unit=10,file="fort.10",formatted=.true.,mode="w",err=fErr,status="unknown",recl=8195)

  sumda(:) = zero
  call time_timerCheck(time1)
  trtime = time1 - time0
  trtime = trtime + crtime3
  sumda(60) = real(trtime,fPrec)  ! Track time
  sumda(52) = real(numvers,fPrec) ! SixTrack version
  outLine = " "
  k = 1
  do i=1,60
    call chr_fromReal(sumda(i),chBuf,19,2,rErr)
    outLine(k:k+25)=" "//chBuf(1:25)
    k = k+26
  end do

  ! Note it COULD happen that napxo is 0 for a very very early error and even napx!!!
  if(napxo == 0 .and. napx == 0) napxo = 1
  if(napxo == 0) napxo = napx
  write(93,"(2(a,i0))") "SIXTRACR> STOP/ABEND writing fort.10, lines ",napxo,"/",napx
  flush(93)
  do j=1,napxo ! Must the dummy file really be napxo times the same dummy line?
    write(10,"(a)",iostat=ierro) outLine(1:60*26)
  end do
  if(ierro /= 0) then
    write(lerr,"(a,i0)") "ABEND> ERROR Problems writing to fort.10. ierro = ",ierro
  end if
#endif

12 continue
  if(lout == 92) then
    write(93,"(a)") "SIXTRACR> STOP/ABEND copying fort.92 to fort.6"
    flush(93)
    rewind(92)
13  continue
    read(92,"(a)",err=14,end=15,iostat=ierro) inLine
    write(6,'(a)',iostat=ierro) trim(inLine)
    goto 13
  end if

14 continue
  write(93,"(a,i0)")   "SIXTRACR> ERROR Reading fort.92 in abend, iostat = ",ierro
  write(lerr,"(a,i0)") "ABEND> ERROR Reading fort.92, iostat = ",ierro

15 continue
  write(6,"(a)",iostat=ierro) "SIXTRACR> Stop "//cstring
  rewind(92) ! Get rid of fort.92
  endfile(92,iostat=ierro)
  call f_close(92)
  write(93,"(a)") "SIXTRACR> Stop "//cstring
  call f_close(93)

#ifdef BOINC
  call copyToStdErr(91,"fort.91")
  call boinc_finish(errout) ! This call does not return
#else
  if(errout /= 0) then
    write(error_unit,"(a,i0)") "ABEND> ERROR Stopping with error ",errout
    stop 1
  else
    stop
  end if
#endif

end subroutine abend
#endif

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Created:   2017-06
!  Rewrittem: 2019-04-15 (VKBO)
!  Updated:   2019-04-15
!  Subroutine to copy the last lines from a file to stderr
!  It is mainly used just before exiting SixTrack in case there was an error.
!  This is useful since STDERR is often returned from batch systems and BOINC.
! =================================================================================================
subroutine copyToStdErr(fUnit,fName)

  use parpro
  use mod_units
  use, intrinsic :: iso_fortran_env, only : error_unit

  implicit none

  integer,          intent(in) :: fUnit
  character(len=*), intent(in) :: fName

  integer, parameter :: maxLines = 10

  integer i, bufIdx, bufMax
  logical isOpen, fErr
  character(len=mInputLn) inLine, inBuf(maxLines)

  inquire(unit=fUnit,opened=isOpen)
  if(isOpen) then
    call f_close(fUnit)
  end if

  fErr = .false.
  call f_open(unit=fUnit,file=fName,formatted=.true.,mode="r",err=fErr,status="old")
  if(fErr) return

  bufIdx = 0
  bufMax = 0
10 continue
  read(fUnit,"(a)",end=20,err=20) inLine
  bufIdx = bufIdx + 1
  if(bufIdx > maxLines) bufIdx = 1
  if(bufIdx > bufMax)   bufMax = bufIdx
  inBuf(bufIdx) = inLine
  goto 10

20 continue
  if(bufIdx > 0) then
    do i=bufIdx+1,bufMax
      write(error_unit,"(a)") trim(inBuf(i))
    end do
    if(bufIdx > 1) then
      do i=1,bufIdx
        write(error_unit,"(a)") trim(inBuf(i))
      end do
    end if
  end if
  call f_close(fUnit)

end subroutine copyToStdErr
