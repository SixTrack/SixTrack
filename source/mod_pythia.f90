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
  use, intrinsic :: iso_c_binding

  implicit none

  logical,            public,  save :: pythia_isActive        = .false.
  logical,            private, save :: pythia_useElastic      = .false.
  logical,            private, save :: pythia_useSDiffractive = .false.

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

  ! Beam Configuration
  integer,            private, save :: pythia_frameType      = 0
  integer,            private, save :: pythia_beamSpecies(2) = pythia_partProton
  real(kind=fPrec),   private, save :: pythia_beamConf(2,3)  = zero     ! E, theta_x, theta_y for beam 1, 2

  ! Other Settings
  character(len=256), private, save :: pythia_settingsFile    = " "
  logical,            private, save :: pythia_useSettingsFile = .false.

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

  case("BEAM")
    if(nSplit < 3 .or. nSplit > 5) then
      write(lout,"(a,i0)") "PYTHIA> ERROR Keyword BEAM expected 2, 3 or 4 arguments, got ",(nSplit-1)
      iErr = .true.
      return
    end if

    if(lnSplit(2) == "CM") then
      pythia_frameType = 1
      iBeam = 1
    else
      call chr_cast(lnSplit(2),iBeam,iErr)
      if(nSplit == 3) then
        pythia_frameType = 2
      else
        pythia_frameType = 3
      end if
    end if

    if(iBeam /= 1 .and. iBeam /= 2) then
      write(lout,"(a)") "PYTHIA> ERROR Keyword BEAM expected parameter 1 to be CM, 1 or 2, got "//trim(lnSplit(2))
      iErr = .true.
      return
    end if

    ! Read E, theta_x and theta_y
    if(nSplit > 2) call chr_cast(lnSplit(3),pythia_beamConf(iBeam,1),iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),pythia_beamConf(iBeam,2),iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),pythia_beamConf(iBeam,3),iErr)

    if(st_debug) then
      if(iBeam == 1) then
        call sixin_echoVal("E(1)",      pythia_beamConf(1,1),"PYTHIA",iLine)
        call sixin_echoVal("theta_x(1)",pythia_beamConf(1,2),"PYTHIA",iLine)
        call sixin_echoVal("theta_y(1)",pythia_beamConf(1,3),"PYTHIA",iLine)
      else
        call sixin_echoVal("E(2)",      pythia_beamConf(2,1),"PYTHIA",iLine)
        call sixin_echoVal("theta_x(2)",pythia_beamConf(2,2),"PYTHIA",iLine)
        call sixin_echoVal("theta_y(2)",pythia_beamConf(2,3),"PYTHIA",iLine)
      end if
    end if

  case default
    write(lout,"(a)") "PYTHIA> ERROR Unknown keyword '"//trim(lnSplit(1))//"'."
    iErr = .true.
    return

  end select

end subroutine pythia_parseInputLine

subroutine pythia_inputParsingDone

  use mathlib_bouncer

  real(kind=fPrec) pxA, pyA, pzA, pxB, pyB, pzB
  real(kind=fPrec) eA, eB, thxA, thxB, thyA, thyB

  ! Set up beams
  call pythia_initBeamType(pythia_frameType, pythia_beamSpecies(1), pythia_beamSpecies(2))

  select case(pythia_frameType)
  case(1)
    call pythia_setBeamCM(pythia_beamConf(1,1))
  case(2)
    call pythia_setBeamEnergy(pythia_beamConf(1,1),pythia_beamConf(2,1))
  case(3)
    eA   = pythia_beamConf(1,1)
    eB   = pythia_beamConf(2,1)
    thxA = pythia_beamConf(1,2)
    thxB = pythia_beamConf(2,2)
    thyA = pythia_beamConf(1,3)
    thyB = pythia_beamConf(2,3)

    pxA  = eA*sin_mb(thxA*c1m3)
    pxB  = eB*sin_mb(thxB*c1m3)
    pyA  = eA*sin_mb(thyA*c1m3)
    pyB  = eB*sin_mb(thyB*c1m3)
    pzA  = sqrt((eA**2 - pxA**2) - pyA**2)
    pzB  = sqrt((eB**2 - pxB**2) - pyB**2)

    call pythia_setBeamMomenta(pxA,pyA,pzA,pxB,pyB,pzB)
  end select

end subroutine pythia_inputParsingDone

! ================================================================================================ !
!  C INTERFACE ROUTINES
! ================================================================================================ !

subroutine pythia_init(partType) bind(C, name="pythiaWrapper_init")
  integer(kind=C_INT), value, intent(in) :: partType
end subroutine pythia_init

subroutine pythia_setBeamType(frameType,idA,idB) bind(C, name="pythiaWrapper_setBeamType")
  integer(kind=C_INT), value, intent(in) :: frameType
  integer(kind=C_INT), value, intent(in) :: idA, idB
end subroutine pythia_setBeamType

subroutine pythia_setBeamCM(eCM) bind(C, name="pythiaWrapper_setBeamCM")
  real(kind=C_DOUBLE), value, intent(in) :: eCM
end subroutine pythia_setBeamCM

subroutine pythia_setBeamEnergy(eA,eB) bind(C, name="pythiaWrapper_setBeamEnergy")
  real(kind=C_DOUBLE), value, intent(in) :: eA, eB
end subroutine pythia_setBeamEnergy

subroutine pythia_setBeamMomenta(pxA,pyA,pzA,pxB,pyB,pzB) bind(C, name="pythiaWrapper_setBeamMomenta")
  real(kind=C_DOUBLE), value, intent(in) :: pxA, pyA, pzA
  real(kind=C_DOUBLE), value, intent(in) :: pxB, pyB, pzB
end subroutine pythia_setBeamMomenta

end module mod_pythia

#endif
