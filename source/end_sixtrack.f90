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
  use checkpoint_restart
  use mod_units
#ifdef BOINC
  use mod_boinc
#endif

  implicit none

  character(len=*), intent(in) :: endMsg

  write(crlog,"(a)") "ABEND_CR> Called"
  write(crlog,"(a)") "ABEND_CR> Closing files"
  flush(crlog)

#ifdef BOINC
  if(endMsg == "ERROR") then
    ! Write the BOINC summary file with exit code 1
    call boinc_summary(1)
  end if
#endif

  call closeUnits
  call cr_copyOut

  write(output_unit,"(a)") "SIXTRACR> Stop: "//trim(endMsg)
  rewind(cr_outUnit)
  endfile(cr_outUnit)

  call copyToStdErr(cr_errUnit,cr_errFile,20)
  call f_close(cr_outUnit)
  write(crlog,"(a)") "ABEND_CR> Stop: "//trim(endMsg)
  call f_close(crlog)

#ifdef BOINC
  if(endMsg == "ERROR") then
    call boinc_finalise(1)
  else
    call boinc_finalise(0)
  end if
#else
  if(endMsg == "ERROR") then
    ! Don't write to stderr, it breaks the error tests.
    write(output_unit,"(a)") "ABEND> Stop: "//trim(endMsg)
    stop 1
  else
    stop
  end if
#endif

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
