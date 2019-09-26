! ================================================================================================ !
!  Random Number Handling
! ~~~~~~~~~~~~~~~~~~~~~~~~
!  Module for handling random numbers and keeping track of seeds for continuation of random
!  number generatio. The module keeps multiple series of random numbers separate between the
!  different modules of SixTrack.
!
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-23
!  Updated: 2019-09-26
! ================================================================================================ !
module mod_random

  use floatPrecision
  use numerical_constants, only : zero

  implicit none

  logical, private, save :: rnd_debug       = .false. ! Master debug flag
  logical, private, save :: rnd_selfTest    = .false. ! Run self test?
  integer, private, save :: rnd_masterSeed  = 0       ! The master seed
  integer, private, save :: rnd_genOverride = 0       ! Whether or not to override the default generators for modules
  integer, private, save :: rnd_luxuryLev   = 3       ! The luxury level for the ranlux generator

  ! Series indices for modules (add new series to the end, otherwise the order is disrupted)
  integer, parameter :: rndser_timeSeed1 = 1 ! Series using ranecu that is always time seeded (non-reproducible)
  integer, parameter :: rndser_timeSeed2 = 2 ! Series using ranlux that is always time seeded (non-reproducible)
  integer, parameter :: rndser_genSeq1   = 3 ! General sequence of ranecu random numbers that can be used by any module
  integer, parameter :: rndser_genSeq2   = 4 ! General sequence of ranlux random numbers that can be used by any module
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
    integer          :: generator  = 0       ! 1 for ranecu, 2 for ranlux
    integer          :: initseed   = 0       ! The initial seed
    integer          :: seeds(25)  = 0       ! The current state of the seeds
    real(kind=fPrec) :: normRemain = zero    ! Storage for last generated normal if requesting odd number
    logical          :: hasNormR   = .false. ! A normal distributed remaining number is cached
  end type type_rndRecord
  type(type_rndRecord), private, save :: rnd_seriesData(rnd_nSeries)

  interface rnd_normal
    module procedure rnd_normal_cut
    module procedure rnd_normal_nocut
  end interface rnd_normal

  interface rnd_rayleigh
    module procedure rnd_rayleigh_nocut
    module procedure rnd_rayleigh_maxcut
    module procedure rnd_rayleigh_maxmincut
  end interface rnd_rayleigh

  private :: rnd_normal_cut
  private :: rnd_normal_nocut
  private :: rnd_rayleigh_nocut
  private :: rnd_rayleigh_maxcut
  private :: rnd_rayleigh_maxmincut

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
!  Updated: 2019-09-24
!  Step 1:
!      The test will first generate an array of integer values from 10 to 20, then generate
!  sequences of random numbers of varying length based on these integers. The length is randomised
!  with a time seed, so should be different each time the test is run. However, since the sequences
!  should be able to continue from any interuption, the series generated with a fixed seed should
!  remain the same. In other words, this subroutine should write a file that does not change even if
!  the generators are called with different length of the random number arrays.
!  Step 2:
!      Write a file of 10000 values for each of the generators, testing varying parameters. The
!  results can be plotted with plotRand.py script in the test/random_numbers folder. This generates
!  the numbers in one go, so does not test continuation.
! ================================================================================================ !
subroutine rnd_runSelfTest

  use crcoall
  use mod_units
  use mod_ranecu
  use mod_ranlux

  character(len=19), parameter :: fFile = "rnd_selftest.dat"
  integer fUnit, rndInt(25), i, j, k, iLine, rndSeeds(25), currSeed
  real(kind=fPrec) rndReal(20,2), rndOut(250,8)

  call f_requestUnit(fFile,fUnit)
  call f_open(unit=fUnit,file=fFile,formatted=.true.,mode="w")

  write(fUnit,"(a3,8(1x,a13))") "###","RANECU", "RANLUX", "RANECU","RANLUX","RANECU",  "RANLUX",  "RANECU",    "RANLUX"
  write(fUnit,"(a3,8(1x,a13))") "###","Uniform","Uniform","Normal","Normal","Rayleigh","Rayleigh","Irwin-Hall","Irwin-Hall"

  call rnd_uniformInt(rndser_timeSeed1, rndSeeds, 25,  0, 99) ! Some random values to send as "seeds"
  call rnd_uniformInt(rndser_timeSeed1, rndInt,   25, 10, 20) ! The sizes of each sequence
  rndSeeds(25) = 0
  do k=1,7,2

    iLine = 1
    call rnd_setSeed(rndser_genSeq1, 1, currSeed, .true.)
    call rnd_setSeed(rndser_genSeq2, 2, currSeed, .true.)
  
    do i=1,25

      rndReal(:,:) = -1.0_fPrec
      write(lout,"(a,i0,a)") "RND> Running selftest with ",rndInt(i)," random numbers"

      ! Interupt the series of random seeds currently in the generators
      call recuin(rndSeeds(1),rndSeeds(2))
      call rluxin(rndSeeds)

      ! Generate numbers from each generator with a random sequence length
      ! These should continue from stored seeds regardless of what was already there
      select case(k)
      case(1)
        call rnd_uniform(rndser_genSeq1, rndReal(:,1), rndInt(i))
        call rnd_uniform(rndser_genSeq2, rndReal(:,2), rndInt(i))
      case(3)
        call rnd_normal(rndser_genSeq1, rndReal(:,1), rndInt(i))
        call rnd_normal(rndser_genSeq2, rndReal(:,2), rndInt(i))
      case(5)
        call rnd_rayleigh(rndser_genSeq1, rndReal(:,1), rndInt(i))
        call rnd_rayleigh(rndser_genSeq2, rndReal(:,2), rndInt(i))
      case(7)
        call rnd_irwinHall(rndser_genSeq1, rndReal(:,1), rndInt(i), 4)
        call rnd_irwinHall(rndser_genSeq2, rndReal(:,2), rndInt(i), 4)
      end select

      do j=1,rndInt(i)
        rndOut(iLine,k)   = rndReal(j,1)
        rndOut(iLine,k+1) = rndReal(j,2)
        iLine = iLine + 1
        if(iLine > 250) goto 10
      end do

    end do
10 continue
  end do

  do i=1,250
    write(fUnit,"(i3,8(1x,f13.9))") i,rndOut(i,1:8)
  end do

  ! Reset seeds and close file
  call rnd_setSeed(rndser_genSeq1, 1, currSeed, .true.)
  call rnd_setSeed(rndser_genSeq2, 2, currSeed, .true.)

  call f_close(fUnit)

  ! The call (almost) all possible combinations and write a test file
  do i=1,9
    call rnd_runSelfTestFull(rndser_genSeq1, i, 10000)
    call rnd_runSelfTestFull(rndser_genSeq2, i, 10000)
    call rnd_setSeed(rndser_genSeq1, 1, currSeed, .true.)
    call rnd_setSeed(rndser_genSeq2, 2, currSeed, .true.)
  end do

end subroutine rnd_runSelfTest

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-26
!  Updated: 2019-09-26
!  Writes a series of values to file for a given seriesID and algorithm.
! ================================================================================================ !
subroutine rnd_runSelfTestFull(seriesID, algID, vLen)

  use mod_units
  use string_tools
  use numerical_constants

  integer, intent(in) :: seriesID
  integer, intent(in) :: algID
  integer, intent(in) :: vLen

  integer i, fUnit
  character(len=64) fName
  real(kind=fPrec) :: rndVal(vLen)

  fName = "rnd_sample_"//chr_toLower(trim(rnd_genName(rnd_seriesData(seriesID)%generator)))//"_"
  select case(algID)
  case(1)
    fName = trim(fName)//"uniform.dat"
    call rnd_uniform(seriesID, rndVal, vLen)
  case(2)
    fName = trim(fName)//"normal_nocut.dat"
    call rnd_normal_nocut(seriesID, rndVal, vLen)
  case(3)
    fName = trim(fName)//"normal_cut.dat"
    call rnd_normal_cut(seriesID, rndVal, vLen, two)
  case(4)
    fName = trim(fName)//"rayleigh_nocut.dat"
    call rnd_rayleigh_nocut(seriesID, rndVal, vLen)
  case(5)
    fName = trim(fName)//"rayleigh_maxcut.dat"
    call rnd_rayleigh_maxcut(seriesID, rndVal, vLen, three)
  case(6)
    fName = trim(fName)//"rayleigh_maxmincut.dat"
    call rnd_rayleigh_maxmincut(seriesID, rndVal, vLen, three, one)
  case(7)
    fName = trim(fName)//"irwinhall_2.dat"
    call rnd_irwinHall(seriesID, rndVal, vLen, 2)
  case(8)
    fName = trim(fName)//"irwinhall_4.dat"
    call rnd_irwinHall(seriesID, rndVal, vLen, 4)
  case(9)
    fName = trim(fName)//"irwinhall_6.dat"
    call rnd_irwinHall(seriesID, rndVal, vLen, 6)
  end select

  call f_requestUnit(fName,fUnit)
  call f_open(unit=fUnit,file=fName,formatted=.true.,mode="w")
  do i=1,vLen
    write(fUnit,"(f13.9)") rndVal(i)
  end do
  call f_freeUnit(fUnit)

end subroutine rnd_runSelfTestFull

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
  call rnd_setSeed(rndser_genSeq1, genRanecu, currSeed, .true.)
  call rnd_setSeed(rndser_genSeq2, genRanlux, currSeed, .true.)
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

  rnd_seriesData(seriesID)%generator  = genID
  rnd_seriesData(seriesID)%initseed   = seedVal
  rnd_seriesData(seriesID)%seeds      = seedArr
  rnd_seriesData(seriesID)%normRemain = zero
  rnd_seriesData(seriesID)%hasNormR   = .false.

  select case(seriesID)
  case(rndser_timeSeed1)
    serName = "RND:  Time seed ranecu"
  case(rndser_timeSeed2)
    serName = "RND:  Time seed ranlux"
  case(rndser_genSeq1)
    serName = "RND:  General sequence ranecu"
  case(rndser_genSeq2)
    serName = "RND:  General sequence ranlux"
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

  integer,          intent(in)  :: seriesID, vLen
  real(kind=fPrec), intent(out) :: rndVec(vLen)

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

  integer, intent(in)  :: seriesID, vLen
  integer, intent(out) :: rndVec(vLen)
  integer, intent(in)  :: iStart, iEnd

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
!  Updated: 2019-09-24
!  Normal distribution with Box-Muller.
!  The two numbers generated are independent, so we generate the numbers in pairs. If the last one
!  is not needed, we save it for the next call in order to ensure the sequence can be generated
!  unintrupted. This will also hold after a restart from checkpoint.
!  At 64-bit, tail truncation occurs at 9.42 sigmas, and 6.66 sigmas at 32-bit.
! ================================================================================================ !
subroutine rnd_normal_nocut(seriesID, rndVec, vLen)

  use crcoall
  use mathlib_bouncer
  use numerical_constants, only : twopi, two

  integer,          intent(in)  :: seriesID, vLen
  real(kind=fPrec), intent(out) :: rndVec(vLen)

  integer i, iS, iGen, iRem
  real(kind=fPrec) radVal, remVal, rndTmp(vLen), rndAdd(2)

  if(rnd_seriesData(seriesID)%hasNormR) then
    ! There's already a random number in the cache
    ! Use this one first, and offset the generation of new ones
    rndVec(1) = rnd_seriesData(seriesID)%normRemain
    rnd_seriesData(seriesID)%hasNormR = .false.
    if(vLen < 2) return
    iS   = 1
    iRem = vLen-1
  else
    iS   = 0
    iRem = vLen
  end if

  iGen = iRem-mod(iRem,2) ! Always request an even number of uniform numbers
  call rnd_uniform(seriesID, rndTmp, iGen)
  do i=1,iGen,2
    ! Generate an even number of random numbers using Box-Muller
    radVal         = sqrt((-two)*log_mb(rndTmp(i)))
    rndVec(i+iS)   = radVal * cos_mb(twopi*rndTmp(i+1))
    rndVec(i+iS+1) = radVal * sin_mb(twopi*rndTmp(i+1))
  end do
  if(mod(iRem,2) /= 0) then
    ! If we need one more, generate another pair and save the last for next time
    call rnd_uniform(seriesID, rndAdd, 2)
    radVal       = sqrt((-two)*log_mb(rndAdd(1)))
    rndVec(vLen) = radVal * cos_mb(twopi*rndAdd(2))
    remVal       = radVal * sin_mb(twopi*rndAdd(2))
    rnd_seriesData(seriesID)%normRemain = remVal
    rnd_seriesData(seriesID)%hasNormR   = .true.
  end if

end subroutine rnd_normal_nocut

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-26
!  Updated: 2019-09-26
!  Normal distribution with Box-Muller, and with a maximum cut in the Rayleigh distribution.
!  Since the remaining number, in case of a request of en odd number of values, cannot be saved for
!  next time (due to the cut potentially being changed) we only use the cosine part. This is less
!  efficient than the uncut Box-Muller.
! ================================================================================================ !
subroutine rnd_normal_cut(seriesID, rndVec, vLen, sCut)

  use crcoall
  use mathlib_bouncer
  use numerical_constants, only : twopi, one, two

  integer,          intent(in)  :: seriesID, vLen
  real(kind=fPrec), intent(out) :: rndVec(vLen)
  real(kind=fPrec), intent(in)  :: sCut

  integer i
  real(kind=fPrec) radVal, cFac, rndTmp(vLen*2)

  cFac = one - exp_mb(-(sCut**2/two))
  call rnd_uniform(seriesID, rndTmp, vLen*2)
  do i=1,vLen
    radVal    = sqrt((-two)*log_mb(one-cFac*rndTmp(2*i-1)))
    rndVec(i) = radVal * cos_mb(twopi*rndTmp(2*i))
  end do

end subroutine rnd_normal_cut

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-24
!  Updated: 2019-09-24
!  Rayleigh distribution with radial part of Box-Muller.
! ================================================================================================ !
subroutine rnd_rayleigh_nocut(seriesID, rndVec, vLen)

  use crcoall
  use mathlib_bouncer
  use numerical_constants, only : two

  integer,          intent(in)  :: seriesID, vLen
  real(kind=fPrec), intent(out) :: rndVec(vLen)

  integer i
  real(kind=fPrec) rndTmp(vLen)

  call rnd_uniform(seriesID, rndTmp, vLen)
  do i=1,vLen
    rndVec(i) = sqrt((-two)*log_mb(rndTmp(i)))
  end do

end subroutine rnd_rayleigh_nocut

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-26
!  Updated: 2019-09-26
!  Rayleigh distribution with radial part of Box-Muller and a max cut.
! ================================================================================================ !
subroutine rnd_rayleigh_maxcut(seriesID, rndVec, vLen, sMax)

  use crcoall
  use mathlib_bouncer
  use numerical_constants, only : one, two

  integer,          intent(in)  :: seriesID, vLen
  real(kind=fPrec), intent(out) :: rndVec(vLen)
  real(kind=fPrec), intent(in)  :: sMax

  integer i
  real(kind=fPrec) rndTmp(vLen), cFac

  cFac = one - exp_mb(-(sMax**2/two))
  call rnd_uniform(seriesID, rndTmp, vLen)
  do i=1,vLen
    rndVec(i) = sqrt((-two)*log_mb(one - cFac*rndTmp(i)))
  end do

end subroutine rnd_rayleigh_maxcut

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-26
!  Updated: 2019-09-26
!  Rayleigh distribution with radial part of Box-Muller and a min and max cut.
! ================================================================================================ !
subroutine rnd_rayleigh_maxmincut(seriesID, rndVec, vLen, sMax, sMin)

  use crcoall
  use mathlib_bouncer
  use numerical_constants, only : one, two

  integer,          intent(in)  :: seriesID, vLen
  real(kind=fPrec), intent(out) :: rndVec(vLen)
  real(kind=fPrec), intent(in)  :: sMax, sMin

  integer i
  real(kind=fPrec) rndTmp(vLen), cFac, sMin2

  cFac  = one - exp_mb(-((sMax+sMin)*(sMax-sMin)/two))
  sMin2 = sMin**2
  call rnd_uniform(seriesID, rndTmp, vLen)
  do i=1,vLen
    rndVec(i) = sqrt(sMin2 - two*log_mb(one - cFac*rndTmp(i)))
  end do

end subroutine rnd_rayleigh_maxmincut

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-24
!  Updated: 2019-09-24
!  Irwin-Hall approximation of a normal distribution of order N = 1:10
!  The distribution is centred around zero and extends from - nOrder/2 to + nOrder/2
!  1st order is equivalen to uniform, 2nd order is triangular, and 3rd order and up approaches a
!  Gaussian distribution. First few orders should be faster than Box-Muller.
! ================================================================================================ !
subroutine rnd_irwinHall(seriesID, rndVec, vLen, nOrder)

  use crcoall
  use numerical_constants, only : two

  integer,          intent(in)  :: seriesID, vLen
  real(kind=fPrec), intent(out) :: rndVec(vLen)
  integer,          intent(in)  :: nOrder

  integer i, k
  real(kind=fPrec) rMean
  real(kind=fPrec) rndTmp(nOrder*vLen)
  real(kind=fPrec) rndSum(nOrder,vLen)

  if(nOrder < 1 .or. nOrder > 10) then
    write(lerr,"(a)") "RND> ERROR Irwin-Hall generator must be between order 1 and 10 (inclusive)"
    call prror
  end if

  rMean = real(nOrder,fPrec)/two
  call rnd_uniform(seriesID, rndTmp, vLen*nOrder)

  rndSum = reshape(rndTmp, [nOrder,vLen])
  rndVec(1:vLen) = sum(rndSum,1)-rMean

end subroutine rnd_irwinHall

end module mod_random
