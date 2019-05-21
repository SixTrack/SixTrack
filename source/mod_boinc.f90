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

contains

! ================================================================================================ !
!  Open the BOINC log file
!  Initialise BOINC
!  Set the checpointing intervals depending on whether BOINC is running in standalone mode or not
! ================================================================================================ !
subroutine boinc_initialise

  use crcoall
  use mod_units
  use mod_common

  integer tmpInt

  ! The logfile is for debugging the checkpoining decisions, and should not be used for other logging
  call f_requestUnit(boinc_logFile, boinc_logUnit)
  call f_open(unit=boinc_logUnit,file=boinc_logFile,formatted=.true.,mode="w+")
  write(crlog,"(a,l1)") "BOINCAPI> Initialising BOINC"
  write(boinc_logBuffer,"(a)") "Initialising BOINC"
  call boinc_writeLog
  flush(crlog)

  call boinc_init
  call boinc_is_standalone(tmpInt)
  boinc_isStandalone = (tmpInt /= 0)
  write(crlog,"(a,l1)") "BOINCAPI> Standalone: ",boinc_isStandalone
  write(boinc_logBuffer,"(a,l1)") "Standalone: ",boinc_isStandalone
  call boinc_writeLog
  flush(crlog)

  if(boinc_isStandalone) then
    boinc_cpInterval   = 10.0
    boinc_progInterval =  1.0
  else
    boinc_cpInterval   = 61.0 ! BOINC manager defaults to 60 seconds. Setting it slightly higher.
    boinc_progInterval =  1.0
  end if

end subroutine boinc_initialise

! ================================================================================================ !
!  Called on every turn in tracking
!  Report progress to BOINC maximum every boinc_progInterval second
!  Checkpoint every boinc_cpInterval seconds, but only if the API allows it
! ================================================================================================ !
subroutine boinc_turn(nTurn)

  use crcoall
  use mod_common
  use checkpoint_restart

  integer, intent(in) :: nTurn

  integer          doCheckpoint
  double precision cpuTime, boincProg

  boinc_nTurn = nTurn
  call cpu_time(cpuTime)

  ! End of tracking should at 99% complete. Last 1% is for postprocessing.
  boincProg = 0.99*dble(nTurn)/dble(numl)

  if(cpuTime-boinc_lastProgress >= boinc_progInterval) then
    ! Tell BOINC how we're doing, but don't hammer the API if many turns
    ! We need to re-add BOINC graphics as well here at some point
    call boinc_fraction_done(boincProg)
    boinc_lastProgress = cpuTime
    write(boinc_logBuffer,"(a,f8.3,a)") "Progress: ",100*boincProg," %"
    call boinc_writeLog
  end if

  if(cpuTime-boinc_lastCPReq < boinc_cpInterval .and. nTurn > 1) then
    return ! Not checkpointing yet, just return
  end if

  write(boinc_logBuffer,"(a,f0.2,a)") "Last checkpoint ",cpuTime-boinc_lastCPReq," seconds ago"
  call boinc_writeLog

  call boinc_fraction_done(boincProg)
  write(boinc_logBuffer,"(a,f8.3,a)") "Progress: ",100*boincProg," %"
  call boinc_writeLog
  if(boinc_isStandalone) then
    write(crlog,"(a)") "BOINCAPI> Standalone Mode: Checkpointing permitted"
    write(boinc_logBuffer,"(a)") "Standalone Mode: Checkpointing permitted"
    call boinc_writeLog
    flush(crlog)
    call crpoint
  else
    call boinc_time_to_checkpoint(doCheckpoint) ! 0 = not allowed, non-0 = allowed
    if(doCheckpoint /= 0) then
      write(crlog,"(a)") "BOINCAPI> Client Mode: Checkpointing permitted"
      write(boinc_logBuffer,"(a)") "Client Mode: Checkpointing permitted"
      call boinc_writeLog
      flush(crlog)
      call crpoint
      call boinc_checkpoint_completed
    else
      write(crlog,"(a)") "BOINCAPI> Client Mode: Checkpointing not permitted"
      write(boinc_logBuffer,"(a)") "Client Mode: Checkpointing not permitted"
      call boinc_writeLog
      flush(crlog)
    end if
  end if

  ! Set the timer for last checkpoint
  boinc_lastCPReq = cpuTime

end subroutine boinc_turn

! ================================================================================================ !
!  Add a final checkpoint after tracking is complete
! ================================================================================================ !
subroutine boinc_post

  use crcoall
  use checkpoint_restart

  write(crlog,"(a)") "BOINCAPI> Tracking completed. Final checkpoint."
  write(boinc_logBuffer,"(a)") "Tracking completed. Final checkpoint."
  write(boinc_logBuffer,"(a,f8.3,a)") "Progress: ",99.0," %"
  call boinc_writeLog
  flush(crlog)

  call boinc_begin_critical_section
  call crpoint
  call boinc_end_critical_section

end subroutine boinc_post

! ================================================================================================ !
!  Sets the postprocessing progress.
!  These steps are hardcoded, and should only be changed if a time consuming post processing routine
!  is added in main_cr. If so, bump the mSteps parameter and recheck all calls to this routine.
!  Currently there are 5 steps:
!   - Step 0 : Called after tracking, called from main_cr
!   - Step 1 : Called after reports, fort.12 and final state file, called from main_cr
!   - Step 2 : Called after postpr, called from main_cr
!   - Step 3 : Called after FMA, called from main_cr
!   - Step 4 : Called after cleaning, closing and hashing, called from main_cr
!   - Step 5 : Called after generating verification files, called from boinc_done
! ================================================================================================ !
subroutine boinc_postProgress(nStep)
  integer, intent(in) :: nStep
  double precision :: progFrac
  integer, parameter :: mSteps = 5
  if(nStep <= mSteps) then
    progFrac = 0.99 + dble(nStep)/dble(mSteps)/100.0
  else
    progFrac = 1.0
  end if
  call boinc_fraction_done(progFrac)
  write(boinc_logBuffer,"(a,f8.3,a)") "Progress: ",100.0*progFrac," %"
  call boinc_writeLog
end subroutine boinc_postProgress

! ================================================================================================ !
!  Bump progress to 100%, and close the log file
! ================================================================================================ !
subroutine boinc_done

  use crcoall
  use mod_hash
  use mod_time
  use mod_units
  use mod_particles

  character(len=32) :: md5Particles
  character(len=32) :: md5Fort10
  real(kind=fPrec)  :: preTime, trackTime, postTime, totalTime
  integer           :: fUnit

  call part_writeState("boinc_particles.dat",.true.,.true.)
  write(crlog,"(a)") "BOINCAPI> Writing particle final state to file 'boinc_particles.dat'"
  flush(crlog)

  call hash_digestFile("boinc_particles.dat", md5Particles, .true.)
  write(crlog,"(a)") "BOINCAPI> MD5SUM of 'boinc_particles.dat' is "//md5Particles
  flush(crlog)

  call hash_digestFile("fort.10", md5Fort10, .true.)
  write(crlog,"(a)") "BOINCAPI> MD5SUM of 'fort.10' is "//md5Fort10
  flush(crlog)

  call boinc_postProgress(5)
  call f_close(boinc_logUnit)

  ! Write the BOINC summary file for validation
  call time_getSummary(preTime, trackTime, postTime, totalTime)
  call f_requestUnit("boinc_summary.dat",fUnit)
  call f_open(unit=fUnit,file="boinc_summary.dat",formatted=.true.,mode="w",status="replace")
  write(fUnit,"(a,2i0.16,a)") md5Particles//char(10)//md5Fort10//char(10),&
    int(trackTime*1e3),int(totalTime*1e3),char(10)
  call f_close(fUnit)

end subroutine boinc_done

! ================================================================================================ !
!  Report exit status and tell BOINC to finish
! ================================================================================================ !
subroutine boinc_finalise(exitCode)
  integer, intent(in) :: exitCode
  call boinc_fraction_done(1.0)
  call boinc_finish(exitCode) ! The API does not return
end subroutine boinc_finalise

! ================================================================================================ !
!  Wrapper for writing the logfile. Just adds a turn number.
! ================================================================================================ !
subroutine boinc_writeLog
  use mod_time, only : time_lastTick, time_ticToc
  call time_ticToc
  if(boinc_logBuffer /= " ") then
    write(boinc_logUnit,"(a,f10.3,a,i8,a)") "[",time_lastTick,"] [",boinc_nTurn,"] "//trim(boinc_logBuffer)
    flush(boinc_logUnit)
    boinc_logBuffer = " "
  end if
end subroutine boinc_writeLog

end module mod_boinc
