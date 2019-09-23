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
  logical, private, save :: rnd_selfTest    = .false. ! Run self test?
  integer, private, save :: rnd_masterSeed  = 0       ! The master seed
  integer, private, save :: rnd_genOverride = 0       ! Whether or not to override the default generators for modules
  integer, private, save :: rnd_luxuryLev   = 3       ! The luxury level for the ranlux generator

  ! Series indices for modules (add new series to the end, otherwise the order is disrupted)
  integer, parameter :: rndser_timeSeed1 = 1 ! Series using ranecu that is always time seeded
  integer, parameter :: rndser_timeSeed2 = 2 ! Series using ranlux that is always time seeded
  integer, parameter :: rndser_selfTest1 = 3 ! Self test of the ranecu generator
  integer, parameter :: rndser_selfTest2 = 4 ! Self test of the ranlux generator
  integer, parameter :: rndser_flucErr   = 5 ! Fluc module magnet errors (zfz array)
  integer, parameter :: rndser_scatMain  = 6 ! Scatter module scatter generator
  integer, parameter :: rndser_scatPart  = 7 ! Scatter module beam sample
  integer, parameter :: rndser_distGen   = 8 ! Dist module generator
  integer, parameter :: rndser_collDist  = 9 ! Collimation module dist generator
  integer, parameter :: rnd_nSeries      = 9 ! Number of random number series to store

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

  case("SELFTEST")
    rnd_selfTest = .true.
    write(lout,"(a)") "RND> Module will run self-test"

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
!  Run a self-test. The test will first generate an array of integer values from 10 to 20, then
!  generate sequences of random numbers of varying length based on these integers. The length is
!  randomised with a time seed, so should be different each time the test is run. However, since
!  the sequences should be able to continue from any interuption, the series generated with a fixed
!  seed should remain the same. In other words, this subroutine should write a file that does not
!  change even if the generators are called with different length of the random number arrays.
! ================================================================================================ !
subroutine rnd_runSelfTest

  use crcoall
  use mod_units
  use mod_ranecu
  use mod_ranlux

  character(len=19), parameter :: fFile = "rnd_selftest.dat"
  integer fUnit, rndInt(25), i, j, k, iLine, rndSeeds(25), currSeed
  real(kind=fPrec) rndReal(20,2), rndOut(250,6)

  call f_requestUnit(fFile,fUnit)
  call f_open(unit=fUnit,file=fFile,formatted=.true.,mode="w")

  write(fUnit,"(a3,4(1x,a13))") "###","RANECU", "RANLUX", "RANECU","RANLUX"
  write(fUnit,"(a3,4(1x,a13))") "###","Uniform","Uniform","Normal","Normal"

  call rnd_uniformInt(rndser_timeSeed1, rndSeeds, 25, 0, 99) ! Some random values to send as "seeds"
  call rnd_uniformInt(rndser_timeSeed1, rndInt,   25, 5, 10) ! The sizes of each sequence
  rndSeeds(25) = 0
  do k=1,3,2

    iLine = 1
    call rnd_setSeed(rndser_selfTest1, 1, currSeed, .true.)
    call rnd_setSeed(rndser_selfTest2, 2, currSeed, .true.)
  
    do i=1,25

      rndReal(:,:) = -1.0_fPrec
      write(lout,"(a,i0,a)") "RND> Running selftest with ",2*rndInt(i)," random numbers"

      ! Interupt the series of random seeds currently in the generators
      call recuin(rndSeeds(1),rndSeeds(2))
      call rluxin(rndSeeds)

      ! Generate numbers from each generator with a random sequence length
      ! These should continue from stored seeds regardless of what was already there
      select case(k)
      case(1)
        call rnd_uniform(rndser_selfTest1, rndReal(:,1), 2*rndInt(i))
        call rnd_uniform(rndser_selfTest2, rndReal(:,2), 2*rndInt(i))
      case(3)
        call rnd_normal (rndser_selfTest1, rndReal(:,1), 2*rndInt(i))
        call rnd_normal (rndser_selfTest2, rndReal(:,2), 2*rndInt(i))
      end select

      do j=1,2*rndInt(i)
        rndOut(iLine,k)   = rndReal(j,1)
        rndOut(iLine,k+1) = rndReal(j,2)
        iLine = iLine + 1
        if(iLine > 250) goto 10
      end do

    end do
10 continue
  end do

  do i=1,250
    write(fUnit,"(i3,4(1x,f13.9))") i,rndOut(i,1:4)
  end do

  ! Reset seeds and close file
  call rnd_setSeed(rndser_selfTest1, 1, currSeed, .true.)
  call rnd_setSeed(rndser_selfTest2, 2, currSeed, .true.)

  call f_close(fUnit)

end subroutine rnd_runSelfTest

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-23
!  Updated: 2019-09-23
!  Initialise the random number generator by setting all of the storage values
! ================================================================================================ !
subroutine rnd_initSeries

  use mod_units
  use mod_time

  integer genRanecu, genRanlux, currSeed, currTime

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

  currTime = time_getSysClock()
  call rnd_setSeed(rndser_timeSeed1, genRanecu, currTime, .false.)
  call rnd_setSeed(rndser_timeSeed2, genRanlux, currTime, .false.)
  call rnd_setSeed(rndser_selfTest1, genRanecu, currSeed, .true.)
  call rnd_setSeed(rndser_selfTest2, genRanlux, currSeed, .true.)
  call rnd_setSeed(rndser_flucErr,   genRanecu, currSeed, .true.)
  call rnd_setSeed(rndser_scatMain,  genRanecu, currSeed, .true.)
  call rnd_setSeed(rndser_scatPart,  genRanecu, currSeed, .true.)
  call rnd_setSeed(rndser_distGen,   genRanecu, currSeed, .true.)
  call rnd_setSeed(rndser_collDist,  genRanlux, currSeed, .true.)

  if(rnd_selfTest) then
    call rnd_runSelfTest
  end if

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
  case(rndser_timeSeed1)
    serName = "RND:  Time seed ranecu"
  case(rndser_timeSeed2)
    serName = "RND:  Time seed ranlux"
  case(rndser_selfTest1)
    serName = "RND:  Self-test ranecu"
  case(rndser_selfTest2)
    serName = "RND:  Self-test ranlux"
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

! ================================================================================================ !
!  GENERATORS
! ================================================================================================ !

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-23
!  Updated: 2019-09-23
!  Uniform distribution between 0 and 1
! ================================================================================================ !
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

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-23
!  Updated: 2019-09-23
!  Uniform distribution between 0 and 1
! ================================================================================================ !
subroutine rnd_uniformInt(seriesID, rndVec, vLen, iStart, iEnd)

  use crcoall

  integer, intent(in)  :: seriesID
  integer, intent(out) :: rndVec(*)
  integer, intent(in)  :: vLen
  integer, intent(in)  :: iStart
  integer, intent(in)  :: iEnd

  real(kind=fPrec) randArr(vLen)
  integer i, rSize

  call rnd_uniform(seriesID, randArr, vLen)

  rSize = iEnd - iStart + 1
  if(rSize < 1) then
    write(lerr,"(a)") "RND> ERROR Uniform int requires end value to be larger than start value"
    call prror
  end if

  do i=1,vLen
    rndVec(i) = int(randArr(i)*real(rSize,fPrec)) + iStart
  end do

end subroutine rnd_uniformInt

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-23
!  Updated: 2019-09-23
!  Normal distribution with Box-Muller using both sine and cosine.
!  At 64-bit, tail truncation occurs at 9.42 sigmas, and 6.66 sigmas at 32-bit.
! ================================================================================================ !
subroutine rnd_normal(seriesID, rndVec, vLen)

  use crcoall
  use mathlib_bouncer
  use numerical_constants, only : twopi, two

  integer,          intent(in)  :: seriesID
  real(kind=fPrec), intent(out) :: rndVec(*)
  integer,          intent(in)  :: vLen

  integer i
  real(kind=fPrec) radVal, rndTmp(vLen), rndAdd(2)

  call rnd_uniform(seriesID, rndTmp, vLen)

  do i=1,vLen-1,2
    radVal      = sqrt((-two)*log_mb(rndTmp(i)))
    rndVec(i)   = radVal * cos_mb(twopi*rndTmp(i+1))
    rndVec(i+1) = radVal * sin_mb(twopi*rndTmp(i+1))
  end do
  if(mod(vLen,2) /= 0) then
    call rnd_uniform(seriesID, rndAdd, 2)
    radVal       = sqrt((-two)*log_mb(rndAdd(1)))
    rndVec(vLen) = radVal * cos_mb(twopi*rndAdd(2))
  end if

end subroutine rnd_normal

end module mod_random
