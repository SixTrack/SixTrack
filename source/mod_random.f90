! ================================================================================================ !
!  Random Number Handling
! ~~~~~~~~~~~~~~~~~~~~~~~~
!  Module for handling random numbers and keeping track of seeds for continuation of random
!  number generatio. The module keeps multiple series of random numbers separate between the
!  different modules of SixTrack.
!
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-23
!  Updated: 2019-09-23
! ================================================================================================ !
module mod_random

  use floatPrecision

  implicit none

  logical, private, save :: rnd_debug       = .false. ! Master debug flag
  integer, private, save :: rnd_masterSeed  = 0       ! The master seed
  integer, private, save :: rnd_genOverride = 0       ! Whether or not to override the default generators for modules
  integer, private, save :: rnd_luxuryLev   = 3       ! The luxury level for the ranlux generator

  ! Series indices for modules (add new series to the end, otherwise the order is disrupted)
  integer, parameter :: rndser_flucErr  = 1 ! Fluc module magnet errors (zfz array)
  integer, parameter :: rndser_scatMain = 2 ! Scatter module scatter generator
  integer, parameter :: rndser_scatPart = 3 ! Scatter module beam sample
  integer, parameter :: rndser_distGen  = 4 ! Dist module generator
  integer, parameter :: rndser_collDist = 5 ! Collimation module dist generator
  integer, parameter :: rnd_nSeries     = 5 ! Number of random number series to store

  integer,           parameter :: rnd_seedStep   = 77377 ! When seeding from a master seed, this is the step between them
  character(len=6),  parameter :: rnd_genName(2) = ["RANECU","RANLUX"]
  character(len=16), parameter :: rnd_seedFile   = "random_seeds.dat"
  integer,       private, save :: rnd_seedUnit   = -1

  type, private :: type_rndRecord
    integer :: generator = 0 ! 1 for ranecu, 2 for ranlux
    integer :: initseed  = 0 ! The initial seed
    integer :: seeds(25) = 0 ! The current state of the seeds
  end type type_rndRecord
  type(type_rndRecord), private, save :: rnd_seriesData(rnd_nSeries)

contains

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-23
!  Updated: 2019-09-23
!  Parse input line for the RANDOM block
! ================================================================================================ !
subroutine rnd_parseInputLine(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common
  use mod_time

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit, i
  logical spErr, cErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "RND> ERROR Failed to parse input line"
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  spErr = .false.
  cErr  = .false.

  select case(lnSplit(1))

  case("DEBUG")
    rnd_debug = .true.
    write(lout,"(a)") "RND> DEBUG mode enabled"

  case("SEED")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "RND> ERROR SEED takes one argument, got ",nSplit-1
      write(lerr,"(a)")    "RND>       SEED seedval|'TIME'"
      iErr = .true.
      return
    end if
    if(chr_toLower(lnSplit(2)) == "time") then
      rnd_masterSeed = time_getSysClock()
    else
      call chr_cast(lnSplit(2), rnd_masterSeed, cErr)
      if(rnd_masterSeed < 0 .or. rnd_masterSeed > 2147483562-rnd_seedStep) then
        write(lerr,"(a,i0)") "RND> ERROR SEED must be an integer between 0 and ",2147483562-rnd_seedStep
      end if
    end if

  end select

end subroutine rnd_parseInputLine

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-23
!  Updated: 2019-09-23
!  Run a self-test
! ================================================================================================ !
subroutine rnd_selfTest

  use mod_units

end subroutine rnd_selfTest

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-23
!  Updated: 2019-09-23
!  Initialise the random number generator by setting all of the storage values
! ================================================================================================ !
subroutine rnd_initSeries

  use mod_units

  integer genRanecu, genRanlux, currSeed

  call f_requestUnit(rnd_seedFile, rnd_seedUnit)
  call f_open(unit=rnd_seedUnit,file=rnd_seedFile,formatted=.true.,mode="w")

  if(rnd_genOverride == 0) then
    genRanecu = 1
    genRanlux = 2
  else
    ! Override the default generators
    genRanecu = rnd_genOverride
    genRanlux = rnd_genOverride
  end if
  currSeed = 0

  write(rnd_seedUnit, "(a,i0)") "# Master Seed         : ",rnd_masterSeed
  write(rnd_seedUnit, "(a,i0)") "# RANLUX Luxury Level : ",rnd_luxuryLev
  write(rnd_seedUnit, "(a)")    "#"
  write(rnd_seedUnit, "(a)")    "# Number Series                  Gen           Seed"

  call rnd_setSeed(rndser_flucErr,  genRanecu, currSeed, .true.)
  call rnd_setSeed(rndser_scatMain, genRanecu, currSeed, .true.)
  call rnd_setSeed(rndser_scatPart, genRanecu, currSeed, .true.)
  call rnd_setSeed(rndser_distGen,  genRanecu, currSeed, .true.)
  call rnd_setSeed(rndser_collDist, genRanlux, currSeed, .true.)

end subroutine rnd_initSeries

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-23
!  Updated: 2019-09-23
!  Set the initial seeds for a random number series
! ================================================================================================ !
subroutine rnd_setSeed(seriesID, genID, seedVal, fromMaster)

  use crcoall
  use mod_ranecu
  use mod_ranlux

  integer, intent(in)    :: seriesID
  integer, intent(in)    :: genID
  integer, intent(inout) :: seedVal
  logical, intent(in)    :: fromMaster

  integer seedArr(25)
  character(len=32) serName

  if(fromMaster) then
    ! Set seed from master seed
    seedVal = rnd_masterSeed + rnd_seedStep*seriesID
  end if

  if(seriesID < 1 .or. seriesID > rnd_nSeries) then
    write(lerr,"(a,i0)") "RND> ERROR SeriesID out of range 1:",rnd_nSeries
    call prror
  end if

  seedArr(:) = 0
  select case(genID)
  case(1) ! Ranecu
    call recuinit(seedVal)
    call recuut(seedArr(1), seedArr(2))
  case(2) ! Ranlux
    call rluxgo(rnd_luxuryLev, seedVal, 0, 0)
    call rluxut(seedArr)
  case default
    write(lerr,"(a,i0)") "RND> ERROR Unknown generator id ",genID
    call prror
  end select

  rnd_seriesData(seriesID)%generator = genID
  rnd_seriesData(seriesID)%initseed  = seedVal
  rnd_seriesData(seriesID)%seeds     = seedArr

  select case(seriesID)
  case(rndser_flucErr)
    serName = "FLUC: Magnet errors"
  case(rndser_scatMain)
    serName = "SCAT: Scattering events"
  case(rndser_scatPart)
    serName = "SCAT: Sample particles"
  case(rndser_distGen)
    serName = "DIST: Distribution generator"
  case(rndser_collDist)
    serName = "COLL: Distribution generator"
  case default
    serName = "Undefined"
  end select
  write(rnd_seedUnit,"(a32,1x,a6,1x,i11)") serName, rnd_genName(genID), seedVal

end subroutine rnd_setSeed

subroutine rnd_uniform(seriesID, rndVec, vLen)

  use crcoall
  use mod_ranecu
  use mod_ranlux

  integer,          intent(in)  :: seriesID
  real(kind=fPrec), intent(out) :: rndVec(*)
  integer,          intent(in)  :: vLen

  integer seedArr(25), tmpArr(25), genID

  seedArr = rnd_seriesData(seriesID)%seeds
  genID   = rnd_seriesData(seriesID)%generator

  select case(genID)
  case(1) ! Raneuc
    call recuut(tmpArr(1), tmpArr(2))
    call recuin(seedArr(1),seedArr(2))
    call ranecuu(rndVec, vLen)
    call recuut(seedArr(1),seedArr(2))
    call recuin(tmpArr(1), tmpArr(2))
    if(rnd_debug) then
      write(lout,"(a,1x,i0,1x,a,2(1x,i0))") "RND> Series:",seriesID," Seeds:",seedArr(1:2)
    end if
  case(2) ! Ranlux
    call rluxut(tmpArr)
    call rluxin(seedArr)
    call ranlux(rndVec, vLen)
    call rluxut(seedArr)
    call rluxin(tmpArr)
    if(rnd_debug) then
      write(lout,"(a,1x,i0,1x,a,25(1x,i0))") "RND> Series:",seriesID," Seeds:",seedArr
    end if
  end select

  rnd_seriesData(seriesID)%seeds = seedArr

end subroutine rnd_uniform

end module mod_random
