! ================================================================================================ !
!  SixTrack Meta Data Module
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-15
!  Records simulation meta data in a file
! ================================================================================================ !
module mod_meta

  use floatPrecision
  use, intrinsic :: iso_fortran_env, only : int64

  implicit none

  character(len=12),   parameter     :: meta_fileName = "sim_meta.dat"
  integer,             private, save :: meta_fileUnit
  logical,             public,  save :: meta_isActive = .false.

  ! Collected MetaData
  integer,             public,  save :: meta_nPartInit = 0   ! Initial number of particles
  integer,             public,  save :: meta_nPartTurn = 0   ! Particle-turns: Counted in tracking routines
  integer(kind=int64), public,  save :: meta_nPTurnEle = 0   ! Particle-turn-eklements: Counted in tracking routines
  integer,             public,  save :: meta_nRestarts = 0   ! Number of C/Rs
  real(kind=fPrec),    public,  save :: meta_sympCheck = 0.0 ! Symplecticity check

  ! Meta Write Interface
  interface meta_write
    module procedure meta_write_char
    module procedure meta_write_real32
    module procedure meta_write_real64
    module procedure meta_write_real128
    module procedure meta_write_int16
    module procedure meta_write_int32
    module procedure meta_write_int64
    module procedure meta_write_log
  end interface meta_write

  private :: meta_write_char
  private :: meta_write_real32
  private :: meta_write_real64
  private :: meta_write_real128
  private :: meta_write_int16
  private :: meta_write_int32
  private :: meta_write_int64
  private :: meta_write_log

#ifdef CR
  integer,             public,  save :: meta_nRestarts_CR = 0
  integer,             private, save :: meta_nPartTurn_CR = 0
  integer(kind=int64), private, save :: meta_nPTurnEle_CR = 0
#endif

contains

subroutine meta_initialise

  use crcoall
  use mod_units

  logical fErr

  fErr = .false.
  call f_requestUnit(meta_fileName, meta_fileUnit)
  call f_open(unit=meta_fileUnit,file=meta_fileName,formatted=.true.,mode="w",err=fErr,status="replace")
  if(fErr) then
    write(lerr,"(a,i0)") "META> ERROR Opening of '"//meta_fileName//"' on unit #",meta_fileUnit
    call prror
  end if

  write(meta_fileUnit,"(a)") "# SixTrack Simulation Meta Data"
  write(meta_fileUnit,"(a)") repeat("#",80)
  flush(meta_fileUnit)

  meta_isActive = .true.

end subroutine meta_initialise

subroutine meta_finalise

  use mod_units
  use mod_alloc
  use mod_common, only : numl

  integer nCRKills1,nCRKills2,tmpUnit
  logical fExist

  nCRKills1 = 0
  nCRKills2 = 0

  call f_requestUnit("crkillswitch.tmp",tmpUnit)
  inquire(file="crkillswitch.tmp",exist=fExist)
  if(fExist) then
    call f_open(unit=tmpUnit,file="crkillswitch.tmp",formatted=.false.,mode="r",access="stream",status="old")
    read(tmpUnit) nCRKills1,nCRKills2
    call f_close(tmpUnit)
  end if

  call meta_write("SymplecticityDeviation",   meta_sympCheck)
  call meta_write("NumParticleTurns",         meta_nPartTurn)
  call meta_write("NumParticleTurnsElement",  meta_nPTurnEle)
  call meta_write("AvgParticlesPerTurn",      real(meta_nPartTurn,fPrec)/numl, "f15.3")
  call meta_write("CR_RestartCount",          meta_nRestarts)
  call meta_write("CR_KillSwitchCount",       nCRKills2)
  call meta_write("PeakDynamicMemAlloc[MiB]", real(maximum_bits,fPrec)/1024/1024/8, "f15.3")
  call meta_write("NumDynamicMemAllocCalls",  alloc_count)

#if defined(GFORTRAN) && defined(MEMUSAGE) && !defined(WIN32)
  call meta_getMemUsage
#endif

  write(meta_fileUnit,"(a)") "# END"
  flush(meta_fileUnit)
  close(meta_fileUnit)

  meta_isActive = .false.

end subroutine meta_finalise

! ================================================================================================ !
!  Interface routines for writing data to sim_meta.dat
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-15
! ================================================================================================ !
subroutine meta_write_char(name, value, fmt)
  character(len=*),           intent(in) :: name
  character(len=*),           intent(in) :: value
  character(len=*), optional, intent(in) :: fmt
  call meta_checkActive
  if(present(fmt)) then
    write(meta_fileUnit,"(a,"//fmt//")") meta_padName(name)//" : "//value
  else
    write(meta_fileUnit,"(a)")           meta_padName(name)//" : "//value
  end if
  flush(meta_fileUnit)
end subroutine meta_write_char

subroutine meta_write_real32(name, value, fmt)
  use, intrinsic :: iso_fortran_env, only : real32
  character(len=*),           intent(in) :: name
  real(kind=real32),          intent(in) :: value
  character(len=*), optional, intent(in) :: fmt
  call meta_checkActive
  if(present(fmt)) then
    write(meta_fileUnit,"(a,"//fmt//")") meta_padName(name)//" : ",value
  else
    write(meta_fileUnit,"(a,es15.7e3)")  meta_padName(name)//" : ",value
  end if
  flush(meta_fileUnit)
end subroutine meta_write_real32

subroutine meta_write_real64(name, value, fmt)
  use, intrinsic :: iso_fortran_env, only : real64
  character(len=*),           intent(in) :: name
  real(kind=real64),          intent(in) :: value
  character(len=*), optional, intent(in) :: fmt
  call meta_checkActive
  if(present(fmt)) then
    write(meta_fileUnit,"(a,"//fmt//")") meta_padName(name)//" : ",value
  else
    write(meta_fileUnit,"(a,es24.16e3)") meta_padName(name)//" : ",value
  end if
  flush(meta_fileUnit)
end subroutine meta_write_real64

subroutine meta_write_real128(name, value, fmt)
  use, intrinsic :: iso_fortran_env, only : real128
  character(len=*),           intent(in) :: name
  real(kind=real128),         intent(in) :: value
  character(len=*), optional, intent(in) :: fmt
  call meta_checkActive
  if(present(fmt)) then
    write(meta_fileUnit,"(a,"//fmt//")") meta_padName(name)//" : ",value
  else
    write(meta_fileUnit,"(a,es41.33e3)") meta_padName(name)//" : ",value
  end if
  flush(meta_fileUnit)
end subroutine meta_write_real128

subroutine meta_write_int16(name, value, fmt)
  use, intrinsic :: iso_fortran_env, only : int16
  character(len=*),           intent(in) :: name
  integer(kind=int16),        intent(in) :: value
  character(len=*), optional, intent(in) :: fmt
  call meta_checkActive
  if(present(fmt)) then
    write(meta_fileUnit,"(a,"//fmt//")") meta_padName(name)//" : ",value
  else
    write(meta_fileUnit,"(a,i6)")        meta_padName(name)//" : ",value
  end if
  flush(meta_fileUnit)
end subroutine meta_write_int16

subroutine meta_write_int32(name, value, fmt)
  use, intrinsic :: iso_fortran_env, only : int32
  character(len=*),           intent(in) :: name
  integer(kind=int32),        intent(in) :: value
  character(len=*), optional, intent(in) :: fmt
  call meta_checkActive
  if(present(fmt)) then
    write(meta_fileUnit,"(a,"//fmt//")") meta_padName(name)//" : ",value
  else
    write(meta_fileUnit,"(a,i11)")       meta_padName(name)//" : ",value
  end if
  flush(meta_fileUnit)
end subroutine meta_write_int32

subroutine meta_write_int64(name, value, fmt)
  use, intrinsic :: iso_fortran_env, only : int64
  character(len=*),           intent(in) :: name
  integer(kind=int64),        intent(in) :: value
  character(len=*), optional, intent(in) :: fmt
  call meta_checkActive
  if(present(fmt)) then
    write(meta_fileUnit,"(a,"//fmt//")") meta_padName(name)//" : ",value
  else
    write(meta_fileUnit,"(a,i20)")       meta_padName(name)//" : ",value
  end if
  flush(meta_fileUnit)
end subroutine meta_write_int64

subroutine meta_write_log(name, value, fmt)
  character(len=*),           intent(in) :: name
  logical,                    intent(in) :: value
  character(len=*), optional, intent(in) :: fmt
  call meta_checkActive
  if(present(fmt)) then
    write(meta_fileUnit,"(a,"//fmt//")") meta_padName(name)//" : ",value
  else
    if(value) then
      write(meta_fileUnit,"(a)") meta_padName(name)//" : true"
    else
      write(meta_fileUnit,"(a)") meta_padName(name)//" : false"
    end if
  end if
  flush(meta_fileUnit)
end subroutine meta_write_log

function meta_padName(inName) result(padName)
  character(len=*), intent(in) :: inName
  character(len=32) padName
  integer i, j, ic
  padName = " "
  j = 0
  do i=1,len(inName)
    ic = ichar(inName(i:i))
    if( ic == 37  .or.  ic == 46  .or. & ! %, .
        ic == 40  .or.  ic == 41  .or. & ! (, )
        ic >= 48  .and. ic <= 57  .or. & ! 0-9
        ic >= 65  .and. ic <= 90  .or. & ! A-Z
        ic == 91  .or.  ic == 93  .or. & ! [, ]
        ic == 45  .or.  ic == 95  .or. & ! -, _
        ic >= 97  .and. ic <= 122 .or. & ! a-z
        ic == 123 .or.  ic == 125      & ! {, }
      ) then
      j = j + 1
      if(j > 32) exit
      padName(j:j) = inName(i:i)
    end if
  end do
end function meta_padName

subroutine meta_checkActive
  use crcoall
  if(meta_isActive .eqv. .false.) then
    write(lerr,"(a)") "META> ERROR Trying to write meta data before initialisation or after finalisation."
    call prror
  end if
end subroutine meta_checkActive

#ifdef CR
! ================================================================================================ !
!  CheckPoint/Restart Routines
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-16
! ================================================================================================ !
subroutine meta_crcheck(fileUnit, readErr)

  use crcoall

  integer, intent(in)  :: fileUnit
  logical, intent(out) :: readErr

  read(fileUnit, err=10, end=10) meta_nRestarts_CR, meta_nPartTurn_CR, meta_nPTurnEle_CR

  readErr = .false.
  return

10 continue
  readErr = .true.
  write(lout, "(a,i0,a)") "SIXTRACR> ERROR Reading C/R file unit ",fileUnit," in META"
  write(crlog,"(a,i0,a)") "SIXTRACR> ERROR Reading C/R file unit ",fileUnit," in META"
  flush(crlog)

end subroutine meta_crcheck

subroutine meta_crpoint(fileUnit, writeErr)

  use crcoall

  integer, intent(in)  :: fileUnit
  logical, intent(out) :: writeErr

  write(fileunit,err=10) meta_nRestarts, meta_nPartTurn, meta_nPTurnEle
  flush(fileUnit)

  writeErr = .false.

  return

10 continue
  writeErr = .true.
  write(lout, "(a,i0,a)") "SIXTRACR> ERROR Writing C/R file unit ",fileUnit," in META"
  write(crlog,"(a,i0,a)") "SIXTRACR> ERROR Writing C/R file unit ",fileUnit," in META"
  flush(crlog)

end subroutine meta_crpoint

subroutine meta_crstart
  meta_nRestarts = meta_nRestarts_CR + 1 ! Restore previous value, and increment
  meta_nPartTurn = meta_nPartTurn_CR
  meta_nPTurnEle = meta_nPTurnEle_CR
end subroutine meta_crstart
#endif

#if defined(GFORTRAN) && defined(MEMUSAGE) && !defined(WIN32)
! ================================================================================================ !
!  Report PID, Peak Virtual Mem Usage and "High Water Mark"
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-04-08
! ================================================================================================ !
subroutine meta_getMemUsage

  use mod_units

  integer cPID
  character(len=30) pPath

  cPID = getpid()
  call meta_write("Exec_PID", cPID)

  write(pPath,"(a,i0,a)") "/proc/",cPID,"/status"
  call execute_command_line("grep VmHWM "//trim(pPath)//" > vmhwm.dat")
  call execute_command_line("grep VmPeak "//trim(pPath)//" > vmpeak.dat")

  call meta_write("Exec_VmHWM[MiB]",  real(meta_extractMemUsage("vmhwm.dat"), kind=fPrec)/1024, "f15.3")
  call meta_write("Exec_VmPeak[MiB]", real(meta_extractMemUsage("vmpeak.dat"),kind=fPrec)/1024, "f15.3")

end subroutine meta_getMemUsage

integer function meta_extractMemUsage(mFile)

  use parpro
  use mod_units
  use string_tools

  character(len=*), intent(in) :: mFile

  character(len=:), allocatable :: lnSplit(:)
  character(len=100) inLine
  logical fErr, sErr
  integer fUnit, nSplit, memKB

  fErr  = .false.
  sErr  = .false.
  memKB = -1024

  call f_requestUnit(mFile,fUnit)
  call f_open(unit=fUnit,file=mFile,formatted=.true.,mode="r",status="old",err=fErr)
  if(fErr) goto 10

  read(fUnit,"(a)",err=10,end=10) inLine
  call chr_split(inLine, lnSplit, nSplit, sErr)
  if(sErr) goto 10
  if(nSplit < 3) goto 10
  call chr_cast(lnSplit(2),memKB,sErr)

10 continue
  call f_close(fUnit)
  meta_extractMemUsage = memKB

end function meta_extractMemUsage
#endif

end module mod_meta
