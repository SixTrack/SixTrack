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

  logical,          public,  save :: pythia_isActive        = .false.
  logical,          private, save :: pythia_useElastic      = .false.
  logical,          private, save :: pythia_useSDiffractive = .false.

  ! Beam Configuration
  logical,          private, save :: pythia_centreOfMass    = .true.
  real(kind=fPrec), private, save :: pythia_beamConf(2,4)   = zero      ! E, Px, Py, Pz for beam 1, 2

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

  case("BEAM")
    if(nSplit /= 6) then
      write(lout,"(a,i0)") "PYTHIA> ERROR Keyword BEAM expected 5 arguments, got ",(nSplit-1)
      iErr = .true.
      return
    end if

    if(lnSplit(2) == "CM") then
      pythia_centreOfMass = .true.
      iBeam = 1
    else
      pythia_centreOfMass = .false.
      call chr_cast(lnSplit(2),iBeam,iErr)
    end if

    if(iBeam /= 1 .and. iBeam /= 2) then
      write(lout,"(a)") "PYTHIA> ERROR Keyword BEAM expected parameter 1 to be CM, 1 or 2, got "//trim(lnSplit(2))
      iErr = .true.
      return
    end if

    ! Read E, Px, Py, Pz
    call chr_cast(lnSplit(3),pythia_beamConf(iBeam,1),iErr)
    call chr_cast(lnSplit(4),pythia_beamConf(iBeam,2),iErr)
    call chr_cast(lnSplit(5),pythia_beamConf(iBeam,3),iErr)
    call chr_cast(lnSplit(6),pythia_beamConf(iBeam,4),iErr)

    if(st_debug) then
      if(iBeam == 1) then
        call sixin_echoVal("E(1)", pythia_beamConf(1,1),"PYTHIA",iLine)
        call sixin_echoVal("Px(1)",pythia_beamConf(1,2),"PYTHIA",iLine)
        call sixin_echoVal("Py(1)",pythia_beamConf(1,3),"PYTHIA",iLine)
        call sixin_echoVal("Pz(1)",pythia_beamConf(1,4),"PYTHIA",iLine)
      else
        call sixin_echoVal("E(2)", pythia_beamConf(2,1),"PYTHIA",iLine)
        call sixin_echoVal("Px(2)",pythia_beamConf(2,2),"PYTHIA",iLine)
        call sixin_echoVal("Py(2)",pythia_beamConf(2,3),"PYTHIA",iLine)
        call sixin_echoVal("Pz(2)",pythia_beamConf(2,4),"PYTHIA",iLine)
      end if
    end if

  case default
    write(lout,"(a)") "PYTHIA> ERROR Unknown keyword '"//trim(lnSplit(1))//"'."
    iErr = .true.
    return

  end select

end subroutine pythia_parseInputLine

subroutine pythia_inputParsingDone
end subroutine pythia_inputParsingDone

! ================================================================================================ !
!  C INTERFACE ROUTINES
! ================================================================================================ !

subroutine pythia_init(partType) bind(C, name="pythiaWrapper_init")
  integer(kind=C_INT), value, intent(in) :: partType
end subroutine pythia_init

end module mod_pythia

#endif
