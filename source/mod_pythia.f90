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

  ! Supported Particles as Defined in Pythia
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

  ! Process Codes as Defined in Pythia
  integer, parameter :: pythia_idNonDiff       = 101
  integer, parameter :: pythia_idElastic       = 102
  integer, parameter :: pythia_idSingleDiffXB  = 103
  integer, parameter :: pythia_idSingleDiffAX  = 104
  integer, parameter :: pythia_idDoubleDiff    = 105
  integer, parameter :: pythia_idCentralDiff   = 106

  ! Flags
  logical,            public,  save :: pythia_isActive        = .false.
  logical,            private, save :: pythia_useElastic      = .false.
  logical,            private, save :: pythia_useSDiffractive = .false.
  logical,            private, save :: pythia_useDDiffractive = .false.
  logical,            private, save :: pythia_useCDiffractive = .false.
  logical,            private, save :: pythia_useNDiffractive = .false.
  logical,            public,  save :: pythia_allowLosses     = .false.
  logical,            private, save :: pythia_useCoulomb      = .false.
  real(kind=fPrec),   private, save :: pythia_elasticTMin     =  5.0e-5_fPrec ! Pythia default value
  real(kind=fPrec),   private, save :: pythia_csElastic       = -1.0_fPrec
  real(kind=fPrec),   private, save :: pythia_csSDiffractive  = -1.0_fPrec
  real(kind=fPrec),   private, save :: pythia_csDDiffractive  = -1.0_fPrec
  real(kind=fPrec),   private, save :: pythia_csCDiffractive  = -1.0_fPrec
  real(kind=fPrec),   private, save :: pythia_csNDiffractive  = -1.0_fPrec

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

    logical(kind=C_BOOL) function pythia_init() bind(C, name="pythiaWrapper_init")
      use, intrinsic :: iso_c_binding
    end function pythia_init

    logical(kind=C_BOOL) function pythia_defaults() bind(C, name="pythiaWrapper_defaults")
      use, intrinsic :: iso_c_binding
    end function pythia_defaults

    subroutine pythia_setProcess(sEL,sSD,sDD,sCD,sND) bind(C, name="pythiaWrapper_setProcess")
      use, intrinsic :: iso_c_binding
      logical(kind=C_BOOL), value, intent(in) :: sEL, sSD, sDD, sCD, sND
    end subroutine pythia_setProcess

    subroutine pythia_setCoulomb(sCMB,tAbsMin) bind(C, name="pythiaWrapper_setCoulomb")
      use, intrinsic :: iso_c_binding
      logical(kind=C_BOOL), value, intent(in) :: sCMB
      real(kind=C_DOUBLE),  value, intent(in) :: tAbsMin
    end subroutine pythia_setCoulomb

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

    subroutine pythia_getCrossSection(sigTot,sigEl) bind(C, name="pythiaWrapper_getCrossSection")
      use, intrinsic :: iso_c_binding
      real(kind=C_DOUBLE), intent(inout) :: sigTot, sigEl
    end subroutine pythia_getCrossSection

    subroutine pythia_getEvent(status,code,t,theta,dEE,dPP) bind(C, name="pythiaWrapper_getEvent")
      use, intrinsic :: iso_c_binding
      logical(kind=C_BOOL), intent(inout) :: status
      integer(kind=C_INT),  intent(inout) :: code
      real(kind=C_DOUBLE),  intent(inout) :: t,theta,dEE,dPP
    end subroutine pythia_getEvent

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
    write(lerr,"(a)") "PYTHIA> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  select case(lnSplit(1))

  case("FILE")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "PYTHIA> ERROR Keyword FILE expected 1 argument, got ",(nSplit-1)
      iErr = .true.
      return
    end if
    pythia_settingsFile    = trim(lnSplit(2))
    pythia_useSettingsFile = .true.
    write(lout,"(a)") "PYTHIA> Settings will be read from external file '"//trim(pythia_settingsFile)//"'"

  case("PROCESS")
    if(nSplit /= 2 .and. nSplit /= 3) then
      write(lerr,"(a,i0)") "PYTHIA> ERROR Keyword PROCESS expected 1 or 2 arguments, got ",(nSplit-1)
      iErr = .true.
      return
    end if
    select case(lnSplit(2))
    case("EL","ELASTIC")
      pythia_useElastic = .true.
      write(lout,"(a)") "PYTHIA> Elastic scattering enabled"
      if(nSplit > 2) then
        call chr_cast(lnSplit(3),pythia_csElastic,iErr)
        write(lout,"(a,f13.8,a)") "PYTHIA> Elastic cross section set to ",pythia_csElastic," mb"
      end if
    case("SD","SINGLEDIFFRACTIVE")
      pythia_useSDiffractive = .true.
      write(lout,"(a)") "PYTHIA> Single diffractive scattering enabled"
      if(nSplit > 2) then
        call chr_cast(lnSplit(3),pythia_csSDiffractive,iErr)
        write(lout,"(a,f13.8,a)") "PYTHIA> Single diffractive cross section set to ",pythia_csSDiffractive," mb"
      end if
    case("DD","DOUBLEDIFFRACTIVE")
      pythia_useDDiffractive = .true.
      write(lout,"(a)") "PYTHIA> Double diffractive scattering enabled"
      if(nSplit > 2) then
        call chr_cast(lnSplit(3),pythia_csDDiffractive,iErr)
        write(lout,"(a,f13.8,a)") "PYTHIA> Double diffractive cross section set to ",pythia_csDDiffractive," mb"
      end if
    case("CD","CENTRALDIFFRACTIVE")
      pythia_useCDiffractive = .true.
      write(lout,"(a)") "PYTHIA> Central diffractive scattering enabled"
      if(nSplit > 2) then
        call chr_cast(lnSplit(3),pythia_csCDiffractive,iErr)
        write(lout,"(a,f13.8,a)") "PYTHIA> Central diffractive cross section set to ",pythia_csCDiffractive," mb"
      end if
    case("ND","NONDIFFRACTIVE")
      pythia_useNDiffractive = .true.
      write(lout,"(a)") "PYTHIA> Non-diffractive scattering enabled"
      if(nSplit > 2) then
        call chr_cast(lnSplit(3),pythia_csNDiffractive,iErr)
        write(lout,"(a,f13.8,a)") "PYTHIA> Non-diffractive cross section set to ",pythia_csNDiffractive," mb"
      end if
    case default
      write(lerr,"(a)") "PYTHIA> ERROR Unknown or unsupported scattering process'"//trim(lnSplit(2))//"'"
      iErr = .true.
      return
    end select

  case("COULOMB")
    if(nSplit /= 2 .and. nSplit /= 3) then
      write(lerr,"(a,i0)") "PYTHIA> ERROR Keyword COULOMB expected 1 or 2 arguments, got ",(nSplit-1)
      iErr = .true.
      return
    end if

    call chr_cast(lnSplit(2),pythia_useCoulomb,iErr)
    if(nSplit == 3) then
      call chr_cast(lnSplit(3),pythia_elasticTMin,iErr)
      if(pythia_elasticTMin < 1.0e-10 .or. pythia_elasticTMin > 1.0e-3) then
        write(lerr,"(a)") "PYTHIA> ERROR Range for COULOMB TMIN is 1e-10 to 1e-3."
        iErr = .true.
        return
      end if
    end if

  case("SPECIES")
    if(nSplit /= 3) then
      write(lerr,"(a,i0)") "PYTHIA> ERROR Keyword SPECIES expected 2 arguments, got ",(nSplit-1)
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
        write(lerr,"(a)") "PYTHIA> ERROR Unknown or unsupported beam species '"//trim(lnSplit(1+iBeam))//"'"
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
      write(lerr,"(a,i0)") "PYTHIA> ERROR Keyword ENERGY expected 2 arguments, got ",(nSplit-1)
      iErr = .true.
      return
    end if

    call chr_cast(lnSplit(2),pythia_beamEnergy(1),iErr)
    call chr_cast(lnSplit(3),pythia_beamEnergy(2),iErr)

    if(st_debug) then
      call sixin_echoVal("E(1)",pythia_beamEnergy(1),"PYTHIA",iLine)
      call sixin_echoVal("E(2)",pythia_beamEnergy(2),"PYTHIA",iLine)
    end if

    pythia_beamEnergy = pythia_beamEnergy*c1m3 ! Pythia expects GeV

  case("SEED")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "PYTHIA> ERROR Keyword SEED expected 1 argument, got ",(nSplit-1)
      iErr = .true.
      return
    end if

    call chr_cast(lnSplit(2),pythia_rndSeed,iErr)

    if(st_debug) then
      call sixin_echoVal("seed",pythia_rndSeed,"PYTHIA",iLine)
    end if

  case default
    write(lerr,"(a)") "PYTHIA> ERROR Unknown keyword '"//trim(lnSplit(1))//"'."
    iErr = .true.
    return

  end select

end subroutine pythia_parseInputLine

subroutine pythia_postInput

  use parpro
  use crcoall
  use, intrinsic :: iso_c_binding

  logical(kind=C_BOOL) pythStat, sEL, sSD, sDD, sCD, sND, sCMB
  integer(kind=C_INT)  rndSeed, frameType, beamSpecies1, beamSpecies2
  real(kind=C_DOUBLE)  beamEnergy1, beamEnergy2, elasticTMin
  real(kind=C_DOUBLE)  sigmaTot, sigmaEl

  rndSeed      = int(pythia_rndSeed,         kind=C_INT)
  frameType    = int(pythia_frameType,       kind=C_INT)
  beamSpecies1 = int(pythia_beamSpecies(1),  kind=C_INT)
  beamSpecies2 = int(pythia_beamSpecies(2),  kind=C_INT)

  beamEnergy1  = real(pythia_beamEnergy(1),      kind=C_DOUBLE)
  beamEnergy2  = real(pythia_beamEnergy(2),      kind=C_DOUBLE)
  elasticTMin  = real(pythia_elasticTMin,        kind=C_DOUBLE)

  sCMB         = logical(pythia_useCoulomb,      kind=C_BOOL)
  sEL          = logical(pythia_useElastic,      kind=C_BOOL)
  sSD          = logical(pythia_useSDiffractive, kind=C_BOOL)
  sDD          = logical(pythia_useDDiffractive, kind=C_BOOL)
  sCD          = logical(pythia_useCDiffractive, kind=C_BOOL)
  sND          = logical(pythia_useNDiffractive, kind=C_BOOL)

  if(pythia_isActive .eqv. .false.) then
    ! No PYTHIA block, so nothing to do.
    return
  end if

  write(lout,"(a)") str_divLine
  write(lout,"(a)") ""
  write(lout,"(a)") "    OOOOOOOOOOOOOOOOOO"
  write(lout,"(a)") "    OO              OO"
  write(lout,"(a)") "    OO    PYTHIA    OO"
  write(lout,"(a)") "    OO              OO"
  write(lout,"(a)") "    OOOOOOOOOOOOOOOOOO"
  write(lout,"(a)") ""

  pythStat = pythia_defaults()
  if(pythStat .eqv. .false.) then
    write(lerr,"(a)") "PYTHIA> ERROR Failed to set default values in libpythia8"
    call prror
  end if

  if(.not. pythia_useElastic .and. pythia_useCoulomb) then
    write(lerr,"(a)") "PYTHIA> ERROR Coulumb corrections to elastic scattering requires elastic scattering to be enabled."
    call prror
  end if
  if(.not. pythia_allowLosses .and. pythia_useDDiffractive) then
    write(lerr,"(a)") "PYTHIA> ERROR Double diffractive scattering requires allowing losses to be enabled."
    call prror
  end if
  if(.not. pythia_allowLosses .and. pythia_useNDiffractive) then
    write(lerr,"(a)") "PYTHIA> ERROR Non-diffractive scattering requires allowing losses to be enabled."
    call prror
  end if

  if(pythia_useSettingsFile) then
    call pythia_readFile(pythia_settingsFile//char(0))
  else
    call pythia_setSeed(rndSeed)
    call pythia_setBeam(frameType,beamSpecies1,beamSpecies2,beamEnergy1,beamEnergy2)
    call pythia_setProcess(sEL,sSD,sDD,sCD,sND)
    call pythia_setCoulomb(sCMB,elasticTMin)
  end if

  pythStat = pythia_init()
  if(pythStat .eqv. .false.) then
    write(lerr,"(a)") "PYTHIA> ERROR Failed to initialise libpythia8"
    call prror
  end if

  call pythia_getCrossSection(sigmaTot, sigmaEl)
  write(lout,"(a)")       "    Cross Sections"
  write(lout,"(a,f12.6)") "    * Total:   ",sigmaTot
  write(lout,"(a,f12.6)") "    * Elastic: ",sigmaEl

  write(lout,"(a)") ""
  write(lout,"(a)") str_divLine

end subroutine pythia_postInput

end module mod_pythia

#endif
