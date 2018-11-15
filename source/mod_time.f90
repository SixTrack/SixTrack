! ================================================================================================ !
!  SixTrack Time Data Module
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-15
!  Records simulation times
! ================================================================================================ !
module mod_time

  use floatPrecision

  implicit none

  character(len=12), parameter     :: time_fileName = "sim_time.dat"
  integer,           private, save :: time_fileUnit

  ! Constants for determining timer checkpoints
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
  integer, parameter :: time_afterZIPF           = 12
  integer, parameter :: time_beforeExit          = 13

  real(kind=fPrec), public, save :: time_timeZero = 0.0
  real(kind=fPrec), public, save :: time_timeRecord(time_beforeExit)

contains

subroutine time_initialise

  use crcoall
  use file_units

  integer ioStat

  call cpu_time(time_timeZero)

  call funit_requestUnit(time_fileName, time_fileUnit)
  open(time_fileUnit,file=time_fileName,status="replace",form="formatted",iostat=ioStat)
  if(ioStat /= 0) then
    write(lout,"(2(a,i0))") "TIME> ERROR Opening of '"//time_fileName//"' on unit #",time_fileUnit," failed with iostat = ",ioStat
    call prror(-1)
  end if

  write(time_fileUnit,"(a)") "# SixTrack Simulation Time Data"
  write(time_fileUnit,"(a)") repeat("#",80)

  call time_writeTime("Internal_ZeroTime",time_timeZero)

  time_timeRecord(:) = 0.0

end subroutine time_initialise

subroutine time_finalise
  close(time_fileUnit)
end subroutine time_finalise

subroutine time_timeStamp(timeStamp)

  integer, intent(in) :: timeStamp

  real(kind=fPrec) timeValue

  call cpu_time(timeValue)
  timeValue = timeValue - time_timeZero
  time_timeRecord(timeStamp) = timeValue

  select case(timeStamp)
  case(time_afterFileUnits)
    call time_writeTime("Stamp_AfterFileUnits",      timeValue)
  case(time_afterDaten)
    call time_writeTime("Stamp_AfterDaten",          timeValue)
  case(time_afterCRCheck)
    call time_writeTime("Stamp_AfterCRCheck",        timeValue)
  case(time_afterClosedOrbit)
    call time_writeTime("Stamp_AfterClosedOrbit",    timeValue)
  case(time_afterBeamDist)
    call time_writeTime("Stamp_AfterBeamDist",       timeValue)
  case(time_afterInitialisation)
    call time_writeTime("Stamp_AfterInitialisation", timeValue)
  case(time_afterPreTrack)
    call time_writeTime("Stamp_AfterPreTrack",       timeValue)
  case(time_afterTracking)
    call time_writeTime("Stamp_AfterTracking",       timeValue)
  case(time_afterPostTrack)
    call time_writeTime("Stamp_AfterPostTrack",      timeValue)
  case(time_afterPostProcessing)
    call time_writeTime("Stamp_AfterPostProcessing", timeValue)
  case(time_afterFMA)
    call time_writeTime("Stamp_AfterFMA",            timeValue)
  case(time_afterZIPF)
    call time_writeTime("Stamp_AfterZIPF",           timeValue)
  case(time_beforeExit)
    call time_writeTime("Stamp_BeforeExit",          timeValue)
  end select

end subroutine time_timeStamp

subroutine time_cpuTime(timeValue)
  real(kind=fPrec), intent(inout) :: timeValue
  call cpu_time(timeValue)
  timeValue = timeValue - time_timeZero
end subroutine time_cpuTime

subroutine time_writeTime(timeLabel, timeValue)
  character(len=*), intent(in) :: timeLabel
  real(kind=fPrec), intent(in) :: timeValue
  integer iLen
  iLen = len_trim(timeLabel)
  if(iLen > 32) then
    write(time_fileUnit,"(a,f13.6,a)") timeLabel(1:32)//" : ",timeValue," sec"
  else
    write(time_fileUnit,"(a,f13.6,a)") trim(timeLabel)//repeat(" ",32-iLen)//" : ",timeValue," sec"
  end if
end subroutine time_writeTime

end module mod_time