! ================================================================================================ !
!
!  SixTrack - Pythia8 Interface Module
!
!  V.K. Berglyd Olsen, BE-ABP-HSS, 2018
!
! ================================================================================================ !

#ifdef PYTHIA

module mod_pythia

  use floatPrecision
  use numerical_constants

  implicit none

  ! Supported Particles
  integer, parameter :: pythia_partProton      =  2212
  integer, parameter :: pythia_partAntiProton  = -2212
  integer, parameter :: pythia_partNeutron     =  2112
  integer, parameter :: pythia_partAntiNeutron = -2112
  integer, parameter :: pythia_partPionPos     =   211
  integer, parameter :: pythia_partPionNeg     =  -211
  integer, parameter :: pythia_partPionZero    =   111
  integer, parameter :: pythia_partPomeron     =   990
  integer, parameter :: pythia_partGamma       =    22
  integer, parameter :: pythia_partElectron    =    11
  integer, parameter :: pythia_partPositron    =   -11
  integer, parameter :: pythia_partMuonNeg     =    13
  integer, parameter :: pythia_partMuonPos     =   -13

  ! Flags
  logical,            public,  save :: pythia_isActive        = .false.
  logical,            private, save :: pythia_useElastic      = .false.
  logical,            private, save :: pythia_useSDiffractive = .false.

  ! Beam Configuration
  integer,            private, save :: pythia_frameType       = 2
  integer,            private, save :: pythia_beamSpecies(2)  = pythia_partProton
  real(kind=fPrec),   private, save :: pythia_beamEnergy(2)   = zero

  ! Other Settings
  character(len=256), private, save :: pythia_settingsFile    = " "
  logical,            private, save :: pythia_useSettingsFile = .false.
  integer,            private, save :: pythia_rndSeed         = -1

  ! C Interface
  interface

    subroutine pythia_init() bind(C, name="pythiaWrapper_init")
    end subroutine pythia_init

    subroutine pythia_defaults() bind(C, name="pythiaWrapper_defaults")
    end subroutine pythia_defaults

    subroutine pythia_setSeed(rndSeed) bind(C, name="pythiaWrapper_setSeed")
      use, intrinsic :: iso_c_binding
      integer(kind=C_INT), value, intent(in) :: rndSeed
    end subroutine pythia_setSeed

    subroutine pythia_setBeam(frameType,idA,idB,eA,eB) bind(C, name="pythiaWrapper_setBeam")
      use, intrinsic :: iso_c_binding
      integer(kind=C_INT), value, intent(in) :: frameType
      integer(kind=C_INT), value, intent(in) :: idA, idB
      real(kind=C_DOUBLE), value, intent(in) :: eA, eB
    end subroutine pythia_setBeam

    subroutine pythia_readFile(fileName) bind(C, name="pythiaWrapper_readFile")
      use, intrinsic :: iso_c_binding
      character(kind=C_CHAR,len=1), intent(in) :: fileName
    end subroutine pythia_readFile

  end interface

contains

subroutine pythia_parseInputLine(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings,   only : st_debug
  use sixtrack_input, only : sixin_echoVal

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit, iBeam
  logical spErr

  ! Split the input line
  call chr_split(inLine,lnSplit,nSplit,spErr)
  if(spErr) then
    write(lout,"(a)") "PYTHIA> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit == 0) return

  select case(lnSplit(1))

  case("FILE")
    if(nSplit /= 2) then
      write(lout,"(a,i0)") "PYTHIA> ERROR Keyword FILE expected 1 arguments, got ",(nSplit-1)
      iErr = .true.
      return
    end if
    pythia_settingsFile    = trim(lnSplit(2))
    pythia_useSettingsFile = .true.
    write(lout,"(a)") "PYTHIA> Settings will be read from external file '"//trim(pythia_settingsFile)//"'"

  case("SPECIES")
    if(nSplit /= 3) then
      write(lout,"(a,i0)") "PYTHIA> ERROR Keyword BEAM expected 2 arguments, got ",(nSplit-1)
      iErr = .true.
      return
    end if
    do iBeam=1,2
      select case(lnSplit(1+iBeam))
      case("P","PROTON")
        pythia_beamSpecies(iBeam) = pythia_partProton
      case("PBAR","ANTIPROTON")
        pythia_beamSpecies(iBeam) = pythia_partAntiProton
      case("N","NEUTRON")
        pythia_beamSpecies(iBeam) = pythia_partNeutron
      case("NBAR","ANTINEUTRON")
        pythia_beamSpecies(iBeam) = pythia_partAntiNeutron
      case("PI+","PION+")
        pythia_beamSpecies(iBeam) = pythia_partPionPos
      case("PI-","PION-")
        pythia_beamSpecies(iBeam) = pythia_partPionNeg
      case("PI0","PION0")
        pythia_beamSpecies(iBeam) = pythia_partPionZero
      case("PHOTON","GAMMA")
        pythia_beamSpecies(iBeam) = pythia_partGamma
      case("POM","POMERON")
        pythia_beamSpecies(iBeam) = pythia_partPomeron
      case("E","E-","ELECTRON")
        pythia_beamSpecies(iBeam) = pythia_partElectron
      case("E+","POSITRON")
        pythia_beamSpecies(iBeam) = pythia_partPositron
      case("MU-","MUON-")
        pythia_beamSpecies(iBeam) = pythia_partMuonNeg
      case("MU+","MUON+")
        pythia_beamSpecies(iBeam) = pythia_partMuonPos
      case default
        write(lout,"(a)") "PYTHIA> ERROR Unknown beam species '"//trim(lnSplit(1+iBeam))//"'"
        iErr = .true.
        return
      end select
    end do

    if(st_debug) then
      call sixin_echoVal("Species(1)",pythia_beamSpecies(1),"PYTHIA",iLine)
      call sixin_echoVal("Species(2)",pythia_beamSpecies(2),"PYTHIA",iLine)
    end if

  case("ENERGY")
    if(nSplit /= 3) then
      write(lout,"(a,i0)") "PYTHIA> ERROR Keyword ENERGY expected 2 arguments, got ",(nSplit-1)
      iErr = .true.
      return
    end if

    call chr_cast(lnSplit(2),pythia_beamEnergy(1),iErr)
    call chr_cast(lnSplit(3),pythia_beamEnergy(2),iErr)

    if(st_debug) then
      call sixin_echoVal("E(1)",pythia_beamEnergy(1),"PYTHIA",iLine)
      call sixin_echoVal("E(2)",pythia_beamEnergy(2),"PYTHIA",iLine)
    end if

  case("SEED")
    if(nSplit /= 2) then
      write(lout,"(a,i0)") "PYTHIA> ERROR Keyword SEED expected 1 argument, got ",(nSplit-1)
      iErr = .true.
      return
    end if

    call chr_cast(lnSplit(2),pythia_rndSeed,iErr)

    if(st_debug) then
      call sixin_echoVal("seed",pythia_rndSeed,"PYTHIA",iLine)
    end if

  case default
    write(lout,"(a)") "PYTHIA> ERROR Unknown keyword '"//trim(lnSplit(1))//"'."
    iErr = .true.
    return

  end select

end subroutine pythia_parseInputLine

subroutine pythia_inputParsingDone

  use crcoall

  write(lout,"(a)") "PYTHIA> Initialising ..."

  call pythia_defaults

  if(pythia_useSettingsFile) then
    call pythia_readFile(pythia_settingsFile)
  else
    call pythia_setSeed(pythia_rndSeed)
    call pythia_setBeam(pythia_frameType,pythia_beamSpecies(1),pythia_beamSpecies(2),pythia_beamEnergy(1),pythia_beamEnergy(2))
  end if

  call pythia_init

end subroutine pythia_inputParsingDone

end module mod_pythia

#endif
