module geant4

  use floatPrecision
  use numerical_constants, only: zero

  implicit none

  real(kind=fPrec) :: g4_recut    = zero
  real(kind=fPrec) :: g4_aecut    = zero
  real(kind=fPrec) :: g4_rcut     = zero
  real(kind=fPrec) :: g4_rangecut = zero
  integer :: g4_physics           = 0
  integer :: g4_keep

  character(len=64) :: phys_str

  logical :: g4_enabled           = .false.
  logical :: g4_debug             = .false.
  logical :: g4_keep_stable       = .false.

contains

subroutine geant4_parseInputLine(inLine,iErr)

  use string_tools
  use crcoall

  implicit none

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable :: lnSplit(:)
  integer nSplit
  logical spErr, cErr

  call chr_split(inLine,lnSplit,nSplit,spErr)
  if(spErr) then
    write(lerr,"(a)") "GEANT4> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit == 0) then
    return
  end if

!Enable/disable debug
  if(lnSplit(1) == 'DEBUG') then
    g4_debug = .true.
    return
  end if

  if(nSplit /= 2) then
    write(lerr,"(a,i0)") "GEANT4> ERROR Expected 2 entries per line, got ",nSplit
    iErr = .true.
    return
  end if

!  For input debugging if needed
!  write(lout,*) '1: ', getfields_fields(1)(1:getfields_lfields(1))
!  write(lout,*) '2: ', getfields_fields(2)(1:getfields_lfields(2))

!relative energy cut
  if(lnSplit(1) == 'RELENERGYCUT') then
    call chr_cast(lnSplit(2),g4_recut,cErr)

!absolute energy cut (GeV)
  else if(lnSplit(1) == 'ABSENERGYCUT') then
    call chr_cast(lnSplit(2),g4_aecut,cErr)

!relative rigidity cut
  else if(lnSplit(1) == 'RELRIGIDITYCUT') then
    call chr_cast(lnSplit(2),g4_rcut,cErr)

!Range cut
  else if(lnSplit(1) == 'RANGECUT') then
    call chr_cast(lnSplit(2),g4_rangecut,cErr)

  else if(lnSplit(1) == 'RETURN') then
    phys_str = trim(lnSplit(2))
    if(phys_str .eq. 'STABLE') then
      g4_keep_stable = .true.
    else if(phys_str .eq. 'IONS') then
      g4_keep = -1
      call g4_keep_id(g4_keep)
    else
      call chr_cast(lnSplit(2),g4_keep,cErr)
      call g4_keep_id(g4_keep)
    end if


!Physics to use number
!FTFP_BERT
!QGSP_BERT
!Anything else? -> error
  else if(lnSplit(1) == 'PHYSICS') then
    phys_str = trim(lnSplit(2))
    if(phys_str .eq. 'FTFP_BERT') then
      g4_physics = 0
    else if(phys_str .eq. 'QGSP_BERT') then
      g4_physics = 1
    else
      write(lout,'(3a)') 'GEANT4> WARNING: Unknown physics model requested: "',phys_str, '" defaulting to FTFP_BERT'
      g4_physics = 0
    end if
  else
    write(lerr,"(2a)") "GEANT4> ERROR unknown keyword", lnSplit(1)
    iErr = .true.
  end if

!Check configuration

!check + enable flags
end subroutine geant4_parseInputLine

subroutine geant4_parseInputDone

  implicit none

  !GEANT4 is enabled
  g4_enabled = .true.
end subroutine geant4_parseInputDone

end module geant4

