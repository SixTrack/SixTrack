! ================================================================================================ !
!  BOINC HELPER MODULE
! ~~~~~~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-05-15
!  Updated: 2019-05-15
! ================================================================================================ !
module mod_boinc

  use floatPrecision

  implicit none

  integer,            private, save :: boinc_nTurn        = 0
  integer,            private, save :: boinc_cpInterval   = -1
  integer,            private, save :: boinc_progInterval = -1
  integer,            private, save :: boinc_logUnit      = -1
  character(len=12),  private, save :: boinc_logFile      = "cr_boinc.log"
  character(len=256), private, save :: boinc_logBuffer    = " "
  logical,            private, save :: boinc_isStandalone = .false.
  real(kind=fPrec),   private, save :: boinc_lastCPoint   = 0.0
  real(kind=fPrec),   private, save :: boinc_lastProgress = 0.0

contains

! ================================================================================================ !
!  Open the BOINC log file
!  Initialise BOINC
!  Set the checpointing intervals depending on whether BOINC is running in standalone mode or not
! ================================================================================================ !
subroutine boinc_initialise

  use mod_units
  use mod_common

  integer tmpInt

  call f_requestUnit(boinc_logFile, boinc_logUnit)
  call f_open(unit=boinc_logUnit,file=boinc_logFile,formatted=.true.,mode="w+")
  write(boinc_logBuffer,"(a)") "Initialising BOINC"
  call boinc_writeLog

#ifdef API
  call boinc_init
  call boinc_is_standalone(tmpInt)
  boinc_isStandalone = (tmpInt /= 0)
#else
  boinc_isStandalone = .true.
  write(boinc_logBuffer,"(a)") "BOINC API not available"
  call boinc_writeLog
#endif
  write(boinc_logBuffer,"(a,l1)") "Standalone: ",boinc_isStandalone
  call boinc_writeLog

  if(boinc_isStandalone) then
    boinc_cpInterval   = 10
    boinc_progInterval = 1
  else
    boinc_cpInterval   = 120
    boinc_progInterval = 1
  end if

end subroutine boinc_initialise

! ================================================================================================ !
!  Called on every turn in tracking
!  Report progress to BOINC maximum every boinc_progInterval second
!  Checkpoint every boinc_cpInterval seconds, but only if the API allows it
! ================================================================================================ !
subroutine boinc_turn(nTurn)

  use mod_common
  use checkpoint_restart

  integer, intent(in) :: nTurn

  integer          doCheckpoint
  real(kind=fPrec) cpuTime, boincProg

  boinc_nTurn = nTurn
  call cpu_time(cpuTime)

#ifdef API
  if(cpuTime-boinc_lastProgress >= boinc_progInterval) then
    ! Tell BOINC how we're doing
    call boinc_fraction_done(dble(nTurn)/dble(numl))
    boinc_lastProgress = cpuTime
  end if
#endif

  if(cpuTime-boinc_lastCPoint < boinc_cpInterval .and. nTurn > 1) then
    return ! Not checkpointing yet, just return
  end if

  write(boinc_logBuffer,"(a,f0.2,a)") "Last checkpoint ",cpuTime-boinc_lastCPoint," seconds ago"
  call boinc_writeLog

  if(cr_checkp) then
#ifdef API
    call boinc_get_fraction_done(boincProg)
    write(boinc_logBuffer,"(a,f8.3,a)") "Progress: ",100*boincProg," %"
    call boinc_writeLog
    if(boinc_isStandalone) then
      write(boinc_logBuffer,"(a)") "Standalone Mode: Checkpointing permitted"
      call boinc_writeLog
      call crpoint
    else
      call boinc_time_to_checkpoint(doCheckpoint) ! 0 = not allowed, non-0 = allowed
      if(doCheckpoint /= 0) then
        write(boinc_logBuffer,"(a)") "Client Mode: Checkpointing permitted"
        call boinc_writeLog
        call crpoint
        call boinc_checkpoint_completed
      else
        write(boinc_logBuffer,"(a)") "Client Mode: Checkpointing not permitted"
        call boinc_writeLog
      end if
    end if
#else
    write(boinc_logBuffer,"(a,f8.3,a)") "Progress: ",100*dble(nTurn)/dble(numl)," %"
    call boinc_writeLog
    write(boinc_logBuffer,"(a)") "Dummy Mode: Checkpointing permitted"
    call boinc_writeLog
    call crpoint
#endif
  end if

  ! Set the timer for last checkpoint
  boinc_lastCPoint = cpuTime

end subroutine boinc_turn

! ================================================================================================ !
!  Add a final checkpoint after tracking is complete
! ================================================================================================ !
subroutine boinc_post

  use checkpoint_restart

  if(.not. cr_checkp) return

  write(boinc_logBuffer,"(a)") "Final checkpoint after tracking"
  call boinc_writeLog

#ifdef API
  call boinc_begin_critical_section
  call crpoint
  call boinc_end_critical_section
#else
  call crpoint
#endif

end subroutine boinc_post

! ================================================================================================ !
!  Report exit status, close the log file, and tell BOINC to finish
! ================================================================================================ !
subroutine boinc_finalise(exitCode)

  use mod_units

  integer, intent(in) :: exitCode

  call f_requestUnit(boinc_logFile, boinc_logUnit)
  call f_open(unit=boinc_logUnit,file=boinc_logFile,formatted=.true.,mode="w+")

  write(boinc_logBuffer,"(a,i0)") "Exiting with status: ",exitCode
  call boinc_writeLog

  call f_close(boinc_logUnit)

#ifdef API
  ! The API does not return
  call boinc_finish(exitCode)
#else
  if(exitCode == 0) then
    stop
  else
    stop exitCode
  end if
#endif

end subroutine boinc_finalise

! ================================================================================================ !
!  Wrapper for writing the logfile. Just adds a turn number.
! ================================================================================================ !
subroutine boinc_writeLog
  if(boinc_logBuffer /= " ") then
    write(boinc_logUnit,"(a,i8,a)") "[",boinc_nTurn,"] "//trim(boinc_logBuffer)
    flush(boinc_logUnit)
    boinc_logBuffer = " "
  end if
end subroutine boinc_writeLog

end module mod_boinc
