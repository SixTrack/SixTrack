! =================================================================================================
!  END SIXTRACK
!  Last modified: 2018-06-25
! =================================================================================================

subroutine prror

  use crcoall
#ifdef FLUKA
  use mod_fluka
#endif

  implicit none

  ! These should not go to lerr
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
  call abend("ERROR")
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
subroutine abend(endMsg)

  use, intrinsic :: iso_fortran_env, only : error_unit, output_unit

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

  character(len=*), intent(in) :: endMsg

  real(kind=fPrec) sumda(60)
  integer i, j, k
  logical fOpen, rErr, fErr, fExists
  character(len=25) chBuf
  character(len=mInputLn) inLine, outLine

  write(crlog,"(a)") "ABEND_CR> Called"
  write(crlog,"(a)") "ABEND_CR> Closing files"
  flush(crlog)

  ! Calling close to be very safe
  call closeUnits

#ifdef BOINC
  ! If fort.10 is non-existent (physics error or some other problem)
  ! we try and write a 0d0 file with a turn number and CPU time
  write(crlog,"(a)") "ABEND_CR> Checking fort.10"
  flush(crlog)

  ! The safest way to check if a file exists is to try to open it and catch the fail
  call f_requestUnit(fort10,unit10) ! Make sure this is actually set
  call f_open(unit=unit10,file=fort10,formatted=.true.,mode="r",err=fErr,status="old",recl=8195)
  if(fErr) goto 11

  ! Now we try and read fort.10 i.e. is it empty?
  read(unit10,"(a1024)",end=11,err=11,iostat=ierro) inLine
  ! Seems to be OK
  goto 12

11 continue
  ! Now we try and write a fort.10
  write(crlog,"(a)") "ABEND_CR> Writing dummy fort.10"
  flush(crlog)

  ! Make sure it is closed properly before we re-open for dummy write
  call f_close(unit10)
  call f_open(unit=unit10,file=fort10,formatted=.true.,mode="w",err=fErr,status="unknown",recl=8195)

  sumda(:) = zero
  call time_timerCheck(time1)
  trtime = time1 - time0
  trtime = trtime + cr_time
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
  write(crlog,"(2(a,i0))") "ABEND_CR> Writing fort.10, lines ",napxo,"/",napx
  flush(crlog)
  do j=1,napxo,2 ! Must the dummy file really be napxo times the same dummy line?
    write(unit10,"(a)",iostat=ierro) outLine(1:60*26)
  end do
  if(ierro /= 0) then
    write(lerr,"(a,i0)") "ABEND> ERROR Problems writing to fort.10. ierro = ",ierro
  end if
#endif

12 continue
  call cr_copyOut

15 continue
  write(output_unit,"(a)",iostat=ierro) "SIXTRACR> Stop: "//trim(endMsg)
  rewind(cr_outUnit)
  endfile(cr_outUnit,iostat=ierro)

  call copyToStdErr(cr_errUnit,cr_errFile,20)
  call f_close(cr_outUnit)
  write(crlog,"(a)") "ABEND_CR> Stop: "//trim(endMsg)
  call f_close(crlog)

  if(endMsg == "ERROR") then
#ifdef BOINC
    call boinc_finish(1)
#endif
    ! Don't write to stderr, it breaks the error tests.
    write(output_unit,"(a)") "ABEND> Stop: "//trim(endMsg)
    stop 1
  else
#ifdef BOINC
    call boinc_finish(0)
#endif
    stop
  end if

end subroutine abend

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Created:   2017-06
!  Rewritten: 2019-04-15 (VKBO)
!  Updated:   2019-05-03
!  Subroutine to copy the last lines from a file to stderr
!  It is mainly used just before exiting SixTrack in case there was an error.
!  This is useful since STDERR is often returned from batch systems and BOINC.
! =================================================================================================
subroutine copyToStdErr(fUnit,fName,maxLines)

  use crcoall
  use parpro
  use mod_units
  use, intrinsic :: iso_fortran_env, only : error_unit

  implicit none

  integer,          intent(in) :: fUnit
  character(len=*), intent(in) :: fName
  integer,          intent(in) :: maxLines

  integer i, bufIdx, bufMax, lnSize, ioStat, szBuf(maxLines)
  logical isOpen, fErr
  character(len=1024) inLine, inBuf(maxLines)

  inquire(unit=fUnit,opened=isOpen)
  if(isOpen) then
    call f_close(fUnit)
  end if

  fErr = .false.
  call f_open(unit=fUnit,file=fName,formatted=.true.,mode="r-",err=fErr,status="old")
  if(fErr) then
    ! It is crucial that we don't let f_open handle thia error and instead exit SixTrack here
    ! since f_open will call prror on unhandled errors, which calls abend, which calls this
    ! routine again, and sends us into an eternal loop.
    write(error_unit,"(a)") "COPYTOERR> ERROR Critical failure while copying to stderr. Exiting."
    stop 1
  end if

  bufIdx = 0
  bufMax = 0
10 continue
  read(fUnit,"(a1024)",end=20,err=20,iostat=ioStat,size=lnSize,advance="no") inLine
  if(ioStat > 0) goto 20 ! Do not use /= 0
  bufIdx = bufIdx + 1
  if(bufIdx > maxLines) bufIdx = 1
  if(bufIdx > bufMax)   bufMax = bufIdx
  inBuf(bufIdx) = inLine
  szBuf(bufIdx) = lnSize
  write(crlog,"(a)") inLine(1:lnSize)
  goto 10

20 continue
  flush(crlog)
  if(bufIdx > 0) then
    do i=bufIdx+1,bufMax
      write(error_unit,"(a)") inBuf(i)(:szBuf(i))
    end do
    if(bufIdx > 1) then
      do i=1,bufIdx
        write(error_unit,"(a)") inBuf(i)(:szBuf(i))
      end do
    end if
  end if
  call f_close(fUnit)

end subroutine copyToStdErr
#endif
