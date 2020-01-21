! ================================================================================================ !
!  SixTrack Timing Module
! ~~~~~~~~~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2018-11-15
!  Updated: 2019-05-29
!
!      This module keeps track of how much time SixTrack spends in the various parts of the code
!  during a simulation run.
!      This is in part done via an array of time stamps which are set when calling time_timeStamp
!  from the relevant places in the code. The total time spent at any given moment is also kept in
!  the variable time_lastTick. This variable is updated when calls are made to this module, but can
!  also be updated by calling time_ticToc. The variable is useful for logging output etc.
!      Since cpu_time is not guaranteed to start at 0 when the executable starts, it will check that
!  the first time stamp is close to zero. If it is not, the offset will be subtracted from all time
!  stamps to ensure they are as accurate as possible.
! ================================================================================================ !
module mod_time

  use floatPrecision

  implicit none

  character(len=12), parameter     :: time_fileName = "sim_time.dat"
  integer,           private, save :: time_fileUnit

  ! Constants for determining timer checkpoints
  ! These are used when calling time_timeStamp() from main_cr and other relevant parts of the code.
  ! To add a new time stamp, add a parameter in the list below in the same order as simulation progress in main_cr,
  ! and if necessary, bump the following parameters so they are also in numerical order.
  ! Then add the corresponding log string to the select case in subroutine time_timeStamp.
  integer, parameter :: time_afterFileUnits      = 1
  integer, parameter :: time_afterDaten          = 2
  integer, parameter :: time_afterCRCheck        = 3
  integer, parameter :: time_afterClosedOrbit    = 4
  integer, parameter :: time_afterBeamDist       = 5
  integer, parameter :: time_afterInitialisation = 6
  integer, parameter :: time_afterPreTrack       = 7
  integer, parameter :: time_afterTracking       = 8
  integer, parameter :: time_afterPostTrack      = 9
  integer, parameter :: time_afterPostProcessing = 10
  integer, parameter :: time_afterFMA            = 11
  integer, parameter :: time_afterAllPostPR      = 12
  integer, parameter :: time_afterHASH           = 13
  integer, parameter :: time_afterZIPF           = 14
  integer, parameter :: time_beforeExit          = 15

  ! Constants for accumulated time
  real(kind=fPrec), public,  save :: time_lastTick = 0.0              ! Latest known cpu-time value. Updated by calling time_ticToc.
  real(kind=fPrec), private, save :: time_timeZero = 0.0              ! Reference time for start of simulation. Usually 0.
  real(kind=fPrec), private, save :: time_timeRecord(time_beforeExit) ! Array of time stamps

  ! Constants for accumulated time
  ! To add a new stop watch, just add a new parameter with an integer value and bump the time_nClocks so the arrays are allocated.
  ! Then it is just a matter of calling time_startClock and time_stopClock around the code to be timed.
  ! The report of the time used must also be added to time_finalise
  integer, parameter :: time_nClocks   = 5  ! Number of stop watches
  integer, parameter :: time_clockCR   = 1
  integer, parameter :: time_clockDUMP = 2
  integer, parameter :: time_clockCOLL = 3
  integer, parameter :: time_clockSCAT = 4
  integer, parameter :: time_clockHDF5 = 5

  real(kind=fPrec), private, save :: time_clockStart(time_nClocks) ! Time at start of stop watch
  real(kind=fPrec), private, save :: time_clockTotal(time_nClocks) ! Total time of stop watch intervals
  integer,          private, save :: time_clockCount(time_nClocks) ! Number of calls to start watch
  integer,          private, save :: time_clockStops(time_nClocks) ! Number of calls to stop watch

#ifdef CR
  ! Additional variables for checkpoint/restart
  real(kind=fPrec), private, save :: time_trackRef = 0.0
  real(kind=fPrec), private, save :: time_trackCR  = 0.0
  real(kind=fPrec), private, save :: time_clockTotal_CR(time_nClocks)
  integer,          private, save :: time_clockCount_CR(time_nClocks)
#endif

contains

subroutine time_initialise

  use crcoall
  use mod_units

  logical fErr

  call cpu_time(time_timeZero)

  fErr = .false.
  call f_requestUnit(time_fileName, time_fileUnit)
  call f_open(unit=time_fileUnit,file=time_fileName,formatted=.true.,mode="w",err=fErr,status="replace")
  if(fErr) then
    write(lerr,"(a,i0)") "TIME> ERROR Opening of '"//time_fileName//"' on unit #",time_fileUnit
    call prror
  end if

  write(time_fileUnit,"(a)") "# SixTrack Simulation Time Data"
  write(time_fileUnit,"(a)") repeat("#",80)
  flush(time_fileUnit)

  call time_writeReal("Internal_ZeroTime",time_timeZero,"s")
  time_lastTick = time_timeZero
  if(time_timeZero > 0.0 .and. time_timeZero < 0.1) then
    ! There is no guarantee that cpu_time is zero at start, but if it is close to 0.0, we will assume that
    ! it was actually the exec start time. If that is the case, it should be within a few ms of 0.0.
    time_timeZero = 0.0
  end if
  call time_writeReal("Stamp_AtStart",time_timeZero,"s")

  time_timeRecord(:) = 0.0
  time_clockStart(:) = 0.0
  time_clockTotal(:) = 0.0
  time_clockCount(:) = 0
  time_clockStops(:) = 0

end subroutine time_initialise

subroutine time_finalise

  use mod_meta
  use mod_units
  use mod_common, only : numl, mbloz
  use numerical_constants, only : zero

  real(kind=fPrec) trackTime, nP, nT, nE, nPT, nPTE

  ! Tracking Averages

  trackTime = time_timeRecord(time_afterTracking) - time_timeRecord(time_afterPreTrack)
  call time_writeReal("Sum_Tracking", trackTime, "s")

  nT   = real(numl,fPrec)
  nE   = real(mbloz,fPrec)
  nPT  = real(meta_nPartTurn,fPrec)
  nPTE = real(meta_nPTurnEle,fPrec)
  if(nT > zero) then
    nP = nPT/nT
  else
    nP = nPT
  end if

  ! Check for zeros just to be safe
  if(nP   > zero) call time_writeReal("Avg_PerParticle",            1.0e3*trackTime/nP,   "ms")
  if(nT   > zero) call time_writeReal("Avg_PerTurn",                1.0e3*trackTime/nT,   "ms")
  if(nE   > zero) call time_writeReal("Avg_PerElement",             1.0e3*trackTime/nE,   "ms")
  if(nPT  > zero) call time_writeReal("Avg_PerParticleTurn",        1.0e6*trackTime/nPT,  "us")
  if(nPTE > zero) call time_writeReal("Avg_PerParticleTurnElement", 1.0e9*trackTime/nPTE, "ns")

  ! Timer Reports
#ifdef CR
  call time_writeReal("Cost_Checkpointing",     time_clockTotal(time_clockCR),   "s", time_clockCount(time_clockCR))
#endif
  call time_writeReal("Cost_DumpModule",        time_clockTotal(time_clockDUMP), "s", time_clockCount(time_clockDUMP))
  call time_writeReal("Cost_CollimationModule", time_clockTotal(time_clockCOLL), "s", time_clockCount(time_clockCOLL))
  call time_writeReal("Cost_ScatterModule",     time_clockTotal(time_clockSCAT), "s", time_clockCount(time_clockSCAT))
  call time_writeReal("Cost_HDF5Module",        time_clockTotal(time_clockHDF5), "s", time_clockCount(time_clockHDF5))

  write(time_fileUnit,"(a)") "# END"
  flush(time_fileUnit)
  call f_close(time_fileUnit)

end subroutine time_finalise

! ================================================================================================ !
!  Return the time used in the three main sections of SixTrack for printing in main_cr
! ================================================================================================ !
subroutine time_getSummary(preTime, trackTime, postTime, totalTime)

  real(kind=fPrec), intent(out) :: preTime
  real(kind=fPrec), intent(out) :: trackTime
  real(kind=fPrec), intent(out) :: postTime
  real(kind=fPrec), intent(out) :: totalTime

  call time_ticToc

  preTime   = time_timeRecord(time_afterPreTrack)
  trackTime = time_timeRecord(time_afterTracking)  - time_timeRecord(time_afterPreTrack)
  postTime  = time_timeRecord(time_afterAllPostPR) - time_timeRecord(time_afterTracking)
  totalTime = time_lastTick

end subroutine time_getSummary

! ================================================================================================ !
!  Update the last clock time. Can be used in log files, etc.
! ================================================================================================ !
subroutine time_ticToc
  call cpu_time(time_lastTick)
#ifdef CR
  time_lastTick = time_lastTick - time_timeZero + time_trackCR
#else
  time_lastTick = time_lastTick - time_timeZero
#endif
end subroutine time_ticToc

! ================================================================================================ !
!  Record simulation time stamps
! ================================================================================================ !
subroutine time_timeStamp(timeStamp)

  integer, intent(in) :: timeStamp

  real(kind=fPrec) timeValue

  call cpu_time(timeValue)
  timeValue = timeValue - time_timeZero

#ifdef CR
  ! Make sure we add checpointed track time to the total time
  if(timeStamp == time_afterPreTrack) then
    time_trackRef = timeValue
  end if
  if(timeStamp > time_afterPreTrack) then
    timeValue = timeValue + time_trackCR
  end if
#endif

  ! Save the time for the current time stamp
  time_timeRecord(timeStamp) = timeValue
  time_lastTick = timeValue

  select case(timeStamp)
  case(time_afterFileUnits)
    call time_writeReal("Stamp_AfterFileUnits",      timeValue, "s")
  case(time_afterDaten)
    call time_writeReal("Stamp_AfterDaten",          timeValue, "s")
  case(time_afterCRCheck)
    call time_writeReal("Stamp_AfterCRCheck",        timeValue, "s")
  case(time_afterClosedOrbit)
    call time_writeReal("Stamp_AfterClosedOrbit",    timeValue, "s")
  case(time_afterBeamDist)
    call time_writeReal("Stamp_AfterBeamDist",       timeValue, "s")
  case(time_afterInitialisation)
    call time_writeReal("Stamp_AfterInitialisation", timeValue, "s")
  case(time_afterPreTrack)
    call time_writeReal("Stamp_AfterPreTrack",       timeValue, "s")
  case(time_afterTracking)
    call time_writeReal("Stamp_AfterTracking",       timeValue, "s")
  case(time_afterPostTrack)
    call time_writeReal("Stamp_AfterPostTrack",      timeValue, "s")
  case(time_afterPostProcessing)
    call time_writeReal("Stamp_AfterPostProcessing", timeValue, "s")
  case(time_afterFMA)
    call time_writeReal("Stamp_AfterFMA",            timeValue, "s")
  case(time_afterAllPostPR)
    call time_writeReal("Stamp_AfterAllPostPR",      timeValue, "s")
  case(time_afterHASH)
    call time_writeReal("Stamp_AfterHASH",           timeValue, "s")
  case(time_afterZIPF)
    call time_writeReal("Stamp_AfterZIPF",           timeValue, "s")
  case(time_beforeExit)
    call time_writeReal("Stamp_BeforeExit",          timeValue, "s")
  end select

end subroutine time_timeStamp

subroutine time_startClock(timerNo)
  integer, intent(in) :: timerNo
  time_clockCount(timerNo) = time_clockCount(timerNo) + 1
  call cpu_time(time_clockStart(timerNo))
end subroutine time_startClock

subroutine time_stopClock(timerNo)
  integer, intent(in) :: timerNo
  real(kind=fPrec) currTime
  call cpu_time(currTime)
  time_clockStops(timerNo) = time_clockStops(timerNo) + 1
  time_clockTotal(timerNo) = time_clockTotal(timerNo) + currTime - time_clockStart(timerNo)
end subroutine time_stopClock

subroutine time_writeReal(timeLabel, timeValue, timeUnit, dataCount)

  character(len=*),  intent(in) :: timeLabel
  real(kind=fPrec),  intent(in) :: timeValue
  character(len=*),  intent(in) :: timeUnit
  integer, optional, intent(in) :: dataCount

  character(len=40) :: endText
  integer iLen

  if(present(dataCount)) then
    iLen = len_trim(timeUnit)
    write(endText,"(a,i0,a)") " "//trim(timeUnit)//repeat(" ",4-iLen)//"[N = ",dataCount,"]"
  else
    endText = " "//timeUnit
  end if

  iLen = len_trim(timeLabel)
  if(iLen > 32) then
    write(time_fileUnit,"(a,f14.6,a)") timeLabel(1:32)//" : ",timeValue,trim(endText)
  else
    write(time_fileUnit,"(a,f14.6,a)") trim(timeLabel)//repeat(" ",32-iLen)//" : ",timeValue,trim(endText)
  end if

  flush(time_fileUnit)

end subroutine time_writeReal

! ================================================================================================ !
!  Extract Integer Value from System Clock
!  Moved from Collimation
!  Last modified: 2019-09-13
! ================================================================================================ !
integer function time_getSysClock()

  use crcoall

  integer mClock, countRate, countMax

  call system_clock(mClock, countRate, countMax)
  if(countMax == 0) then
    write(lout,"(a)") "TIME> System Clock not present or not responding"
  end if
  time_getSysClock = mClock

end function time_getSysClock

#ifdef CR
! ================================================================================================ !
!  CheckPoint/Restart Routines
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-16
! ================================================================================================ !
subroutine time_crcheck(fileUnit, readErr)

  use crcoall

  integer, intent(in)  :: fileUnit
  logical, intent(out) :: readErr

  read(fileUnit, err=10, end=10) time_clockTotal_CR, time_clockCount_CR, time_trackCR

  readErr = .false.
  return

10 continue
  readErr = .true.
  write(lout, "(a,i0,a)") "SIXTRACR> ERROR Reading C/R file unit ",fileUnit," in TIME"
  write(crlog,"(a,i0,a)") "SIXTRACR> ERROR Reading C/R file unit ",fileUnit," in TIME"
  flush(crlog)

end subroutine time_crcheck

subroutine time_crpoint(fileUnit, writeErr)

  use crcoall

  integer, intent(in)  :: fileUnit
  logical, intent(out) :: writeErr

  real(kind=fPrec) timeValue
  call cpu_time(timeValue)
  time_lastTick = timeValue - time_timeZero + time_trackCR
  timeValue     = timeValue - time_trackRef + time_trackCR

  write(fileunit,err=10) time_clockTotal, time_clockCount, timeValue
  flush(fileUnit)

  writeErr = .false.

  return

10 continue
  writeErr = .true.
  write(lout, "(a,i0,a)") "SIXTRACR> ERROR Writing C/R file unit ",fileUnit," in TIME"
  write(crlog,"(a,i0,a)") "SIXTRACR> ERROR Writing C/R file unit ",fileUnit," in TIME"
  flush(crlog)

end subroutine time_crpoint

subroutine time_crstart
  time_clockTotal = time_clockTotal_CR
  time_clockCount = time_clockCount_CR
  call time_ticToc
end subroutine time_crstart
#endif

end module mod_time
