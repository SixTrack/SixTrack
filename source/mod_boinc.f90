! ================================================================================================ !
!  BOINC SERVICE MODULE
! ~~~~~~~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-05-15
!  Updated: 2019-06-05
! ================================================================================================ !
module mod_boinc

  use floatPrecision

  implicit none

  integer,            private, save :: boinc_nTurn        = 0              ! The current turn in tracking
  integer,            private, save :: boinc_logUnit      = -1             ! BOIND C/R debug log unit
  character(len=12),  private, save :: boinc_logFile      = "cr_boinc.log" ! BOIND C/R debug log file
  character(len=256), private, save :: boinc_logBuffer    = " "            ! Log buffer
  character(len=999), private, save :: boinc_sumBuffer    = " "            ! The buffer for the summary file
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
  write(boinc_logBuffer,"(a)") "Initialising BOINC"
  call boinc_writeLog

  call boinc_init
  call boinc_is_standalone(tmpInt)
  boinc_isStandalone = (tmpInt /= 0)

  write(boinc_logBuffer,"(a,l1)") "Standalone: ",boinc_isStandalone
  call boinc_writeLog

  if(boinc_isStandalone) then
    boinc_cpInterval   = 10.0
    boinc_progInterval =  1.0
  else
    boinc_cpInterval   = 61.0 ! BOINC manager defaults to 60 seconds. Setting it slightly higher.
    boinc_progInterval =  1.0
  end if

  call boinc_fraction_done(0.0)
  call boinc_summary(0) ! Make sure we have a summary file as early as possible in case of a crash

end subroutine boinc_initialise

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-05-15
!  Updated: 2019-05-21
!
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

  ! Tracking runs from 1% to 99%. First 1% is for pre tracking, last 1% for post tracking
  boincProg = 0.01 + 0.98*dble(nTurn)/dble(numl)

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
  call boinc_summary(0) ! Make sure the summary file is up to date
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
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-05-15
!  Updated: 2019-05-21
!
!  Called post tracking.
!  Adds a final checkpoint after tracking is complete.
! ================================================================================================ !
subroutine boinc_post

  use crcoall
  use checkpoint_restart

  write(crlog,"(a)") "BOINCAPI> Tracking completed. Final checkpoint."
  write(boinc_logBuffer,"(a)") "Tracking completed. Final checkpoint."
  call boinc_writeLog
  flush(crlog)

  call boinc_begin_critical_section
  call crpoint
  call boinc_end_critical_section

end subroutine boinc_post

! ================================================================================================ !
!  Sets the pre tracking progress.
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-05-21
!  Updated: 2019-05-21
!  These steps are hardcoded, and should only be changed if a time consuming pre tracking routine
!  is added in main_cr. If so, bump the mSteps parameter and recheck all calls to this routine.
!  Currently there are 5 steps:
!   - Step 0 : Implicit, boinc_initialisation
!   - Step 1 : Called after checpoint/restart init, called from main_cr
!   - Step 2 : Called before daten, called from main_cr
!   - Step 3 : Called after daten, called from main_cr
!   - Step 4 : Called before first closed orbit call, called from main_cr
!   - Step 5 : Called after linopt, called from main_cr
!   - Step 6 : Called after closed orbit, called from main_cr
!   - Step 7 : Called after beam distribution, called from main_cr
!   - Step 8 : Called after pre-tracking, called from main_cr
! ================================================================================================ !
subroutine boinc_preProgress(nStep)
  integer, intent(in) :: nStep
  double precision :: progFrac
  integer, parameter :: mSteps = 8
  if(nStep <= mSteps) then
    progFrac = dble(nStep)/dble(mSteps)/100.0
  else
    progFrac = 0.01
  end if
  call boinc_fraction_done(progFrac)
  call boinc_summary(0)
  write(boinc_logBuffer,"(a,f8.3,a)") "Progress: ",100.0*progFrac," %"
  call boinc_writeLog
end subroutine boinc_preProgress

! ================================================================================================ !
!  Sets the post tracking progress.
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-05-21
!  Updated: 2019-05-21
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
  call boinc_summary(0)
  write(boinc_logBuffer,"(a,f8.3,a)") "Progress: ",100.0*progFrac," %"
  call boinc_writeLog
end subroutine boinc_postProgress

! ================================================================================================ !
!  Writes the validation file header, and creates the file. Hash values must be set before this.
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-05-22
!  Updated: 2019-05-22
! ================================================================================================ !
subroutine boinc_summary(exitStatus)

  use mod_time
  use mod_meta
  use mod_units
  use mod_common
  use mod_version

  integer, intent(in) :: exitStatus

  integer          :: fUnit, lenVer
  double precision :: progFrac

  call boinc_get_fraction_done(progFrac)
  call time_ticToc

  lenVer = len_trim(boinc_sumBuffer) - 24
  if(lenVer < 1) lenVer = 7

  write(boinc_sumBuffer( 1:9 ) ,"(i9.9)")   int(time_lastTick*1.0e3)
  write(boinc_sumBuffer(11:16) ,"(i6.6)")   int(progFrac*1.0e5)
  write(boinc_sumBuffer(18:18) ,"(i1.1)")   exitStatus
  write(boinc_sumBuffer(20:24) ,"(i5.5)")   lenVer
  write(boinc_sumBuffer(26:31) ,"(i6.6)")   numvers
  write(boinc_sumBuffer(33:41) ,"(i9.9)")   iu
  write(boinc_sumBuffer(43:51) ,"(i9.9)")   numl
  write(boinc_sumBuffer(53:61) ,"(i9.9)")   napx
  write(boinc_sumBuffer(63:71) ,"(i9.9)")   napxo
  write(boinc_sumBuffer(73:82) ,"(i10.10)") meta_nPartTurn
  write(boinc_sumBuffer(84:101),"(i18.18)") meta_nPTurnEle

  ! Write the BOINC summary file as a binary stream to avoid line endings
  call f_requestUnit("boinc_summary.dat",fUnit)
  call f_open(unit=fUnit,file="boinc_summary.dat",formatted=.false.,access="stream",mode="w",status="replace")
  write(fUnit) trim(boinc_sumBuffer)
  call f_close(fUnit)

end subroutine boinc_summary

! ================================================================================================ !
!  Write the validation file, then bump progress to 100%, and close the log file
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-05-20
!  Updated: 2019-05-22
! ================================================================================================ !
subroutine boinc_done

  use crcoall
  use mod_hash
  use mod_units
  use mod_particles
  use fma,        only : fma_flag, fma_fileName
  use mod_common, only : ipos, unit10, fort10

  character(len=19), parameter :: boincSum = "boinc_particles.dat"
  character(len=32) :: md5Digest
  integer           :: cPos

  boinc_sumBuffer = " "
  cPos = 103

  ! Write our own final state file that does not interfere with the user's sttings in fort.3
  call part_writeState(boincSum,.true.,.true.)
  write(crlog,"(a)") "BOINCAPI> Wrote particle final state to file '"//boincSum//"'"
  flush(crlog)

  ! This file is always hashed and added to the summary
  call hash_digestFile(boincSum, md5Digest, .true.)
  write(crlog,"(a)") "BOINCAPI> MD5SUM '"//boincSum//"': "//md5Digest
  flush(crlog)
  boinc_sumBuffer(cPos:cPos+32) = md5Digest
  cPos = cPos + 33

  ! If postprocessing, try to hash the fort.10
  ! The hash module will handle the case of a missing file
  if(ipos == 1) then
    call f_close(unit10) ! Make sure it is closed before we hash it. This also flushes it.
    call hash_digestFile(fort10, md5Digest, .true.)
    write(crlog,"(a)") "BOINCAPI> MD5SUM '"//trim(fort10)//"': "//md5Digest
    flush(crlog)
    boinc_sumBuffer(cPos:cPos+32) = md5Digest
    cPos = cPos + 33
  end if

  ! If we ran FMA, hash that file too
  if(fma_flag) then
    call hash_digestFile(fma_fileName, md5Digest, .true.)
    write(crlog,"(a)") "BOINCAPI> MD5SUM '"//fma_fileName//"': "//md5Digest
    flush(crlog)
    boinc_sumBuffer(cPos:cPos+32) = md5Digest
    cPos = cPos + 33
  end if

  ! Clean up service module bits, and tell BOINC we're done
  call boinc_postProgress(5)
  call f_close(boinc_logUnit)

  ! Write the BOINC summary file
  call boinc_summary(0)

end subroutine boinc_done

! ================================================================================================ !
!  Report exit status and tell BOINC to finish
! ================================================================================================ !
subroutine boinc_finalise(exitCode)
  integer, intent(in) :: exitCode
  call boinc_fraction_done(1.0) ! Skip to 100%
  call boinc_finish(exitCode)   ! The API does not return
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
