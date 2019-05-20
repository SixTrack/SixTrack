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

  integer,            private, save :: boinc_nTurn        = 0              ! The current turn in tracking
  integer,            private, save :: boinc_logUnit      = -1             ! BOIND C/R debug log unit
  character(len=12),  private, save :: boinc_logFile      = "cr_boinc.log" ! BOIND C/R debug log file
  character(len=256), private, save :: boinc_logBuffer    = " "            ! Log buffer
  logical,            private, save :: boinc_isStandalone = .false.        ! True when the BOINC API is not talking to the manager
  real(kind=fPrec),   private, save :: boinc_cpInterval   = 0.0            ! Number of seconds between checkpoints
  real(kind=fPrec),   private, save :: boinc_progInterval = 0.0            ! Number of seconds between progress updates
  real(kind=fPrec),   private, save :: boinc_lastCPReq    = 0.0            ! CPU time of last checkpointing request
  real(kind=fPrec),   private, save :: boinc_lastProgress = 0.0            ! CPU time of last progress report to the API
  real(kind=fPrec),   private, save :: boinc_trackTotal   = 1.0            ! Max progress count in tracking

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

  ! The logfile is for debugging the checkpoining decisions, and should not be used for other logging
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
    boinc_cpInterval   = 10.0
    boinc_progInterval =  1.0
  else
    boinc_cpInterval   = 61.0 ! BOINC manager defaults to 60 seconds. Setting it slightly higher.
    boinc_progInterval =  1.0
  end if

  if(ipos == 1) then
    boinc_trackTotal = 0.99 ! If posprocesing, tracking progress ends at 99%
  else
    boinc_trackTotal = 1.00
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
    ! We need to re-add BOINC graphics as well here at some point
    call boinc_fraction_done(boinc_trackTotal*dble(nTurn)/dble(numl))
    boinc_lastProgress = cpuTime
  end if
#endif

  if(cpuTime-boinc_lastCPReq < boinc_cpInterval .and. nTurn > 1) then
    return ! Not checkpointing yet, just return
  end if

  write(boinc_logBuffer,"(a,f0.2,a)") "Last checkpoint ",cpuTime-boinc_lastCPReq," seconds ago"
  call boinc_writeLog

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

  ! Set the timer for last checkpoint
  boinc_lastCPReq = cpuTime

end subroutine boinc_turn

! ================================================================================================ !
!  Add a final checkpoint after tracking is complete
! ================================================================================================ !
subroutine boinc_post

  use checkpoint_restart

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
!  PostProcessing Progress
! ================================================================================================ !
subroutine boinc_postpr(nPair)

  use mod_common, only : napxo

  integer, intent(in) :: nPair

  real(kind=fPrec) boincProg

#ifdef API
  call boinc_fraction_done(boinc_trackTotal + dble(nPair)/dble(napxo)/100)
  call boinc_get_fraction_done(boincProg)
  write(boinc_logBuffer,"(a,f8.3,a)") "PostPR Progress: ",100*boincProg," %"
#else
  write(boinc_logBuffer,"(a,f8.3,a)") "PostPR Progress: ",100*boinc_trackTotal + dble(nPair)/dble(napxo)," %"
#endif
  call boinc_writeLog

end subroutine boinc_postpr

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
  call boinc_fraction_done(1.0)
  call boinc_finish(exitCode)
#else
  if(exitCode == 0) then
    stop
  else
    stop 1
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
