module geant4

  use floatPrecision
  use numerical_constants, only: zero
  use, intrinsic :: iso_c_binding

  implicit none

  real(kind=fPrec) :: g4_recut       = zero
  real(kind=fPrec) :: g4_aecut       = zero
  real(kind=fPrec) :: g4_rcut        = zero
  real(kind=fPrec) :: g4_rangecut_mm = zero
  integer :: g4_physics              = 0
  integer :: g4_keep

  character(len=64) :: phys_str

  logical(kind=C_BOOL) :: g4_enabled           = .false.
  logical(kind=C_BOOL) :: g4_debug             = .false.
  logical(kind=C_BOOL) :: g4_keep_stable       = .false.

 interface

!void g4_collimation_init_(double* ReferenceE, int* seed, double* recut, double* aecut, double* rcut,
!int* PhysicsSelect, bool* g4_debug, bool* g4_keep_stable)
  subroutine g4_collimation_init(ReferenceE, seed, recut, aecut, rcut, rangecut_mm, PhysicsSelect, g4_debug, g4_keep_stable) &
& bind(C,name="g4_collimation_init")
    use, intrinsic :: iso_c_binding
    implicit none
    real(kind=C_DOUBLE),          intent(in) :: ReferenceE
    integer(kind=C_INT),          intent(in) :: seed
    real(kind=C_DOUBLE),          intent(in) :: recut
    real(kind=C_DOUBLE),          intent(in) :: aecut
    real(kind=C_DOUBLE),          intent(in) :: rcut
    real(kind=C_DOUBLE),          intent(in) :: rangecut_mm
    integer(kind=C_INT),          intent(in) :: PhysicsSelect
    logical(kind=C_BOOL),         intent(in) :: g4_debug
    logical(kind=C_BOOL),         intent(in) :: g4_keep_stable
  end subroutine g4_collimation_init

!void g4_add_collimator_(char* name, char* material, double* length, double* aperture, double* rotation, double* offset)
  subroutine g4_add_collimator(name, material, length, aperture, rotation, offset) &
& bind(C,name="g4_add_collimator")
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=C_CHAR,len=1), intent(in) :: name
    character(kind=C_CHAR,len=1), intent(in) :: material
    real(kind=C_DOUBLE),          intent(in) :: length
    real(kind=C_DOUBLE),          intent(in) :: aperture
    real(kind=C_DOUBLE),          intent(in) :: rotation
    real(kind=C_DOUBLE),          intent(in) :: offset
  end subroutine g4_add_collimator

!void g4_terminate_()
  subroutine g4_terminate() &
& bind(C,name="g4_terminate")
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine g4_terminate

!void g4_set_collimator_(char* name)
  subroutine g4_set_collimator(name) &
& bind(C,name="g4_set_collimator")
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=C_CHAR,len=1), intent(in) :: name
  end subroutine g4_set_collimator

!void g4_add_particle_(double* x, double* y, double* xp, double* yp, double* e, int32_t* pdgid, int16_t* nzz, int16_t* naa,
! int16_t* nqq, double* mass)
  subroutine g4_add_particle(x, y, xp, yp, e, pdgid, nzz, naa, nqq, mass) &
& bind(C,name="g4_add_particle")
    use, intrinsic :: iso_c_binding
    implicit none
    real(kind=C_DOUBLE),          intent(in) :: x
    real(kind=C_DOUBLE),          intent(in) :: y
    real(kind=C_DOUBLE),          intent(in) :: xp
    real(kind=C_DOUBLE),          intent(in) :: yp
    real(kind=C_DOUBLE),          intent(in) :: e
    integer(kind=C_INT32_T),      intent(in) :: pdgid
    integer(kind=C_INT16_T),      intent(in) :: nzz
    integer(kind=C_INT16_T),      intent(in) :: naa
    integer(kind=C_INT16_T),      intent(in) :: nqq
    real(kind=C_DOUBLE),          intent(in) :: mass
  end subroutine g4_add_particle

!void g4_collimate_()
  subroutine g4_collimate() &
& bind(C,name="g4_collimate")
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine g4_collimate

!void g4_collimate_return_(int* j, double* x, double* y, double* xp, double* yp, double* e, int32_t* pdgid, double* m,
! int16_t* z, int16_t* a, int16_t* q, int *part_hit, int *part_abs, double *part_impact, double *part_indiv, double *part_linteract)
  subroutine g4_collimate_return(j, x, y, xp, yp, e, pdgid, m, z, a, q, part_hit, part_abs, part_impact, part_indiv, &
& part_linteract) &
& bind(C,name="g4_collimate_return")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(kind=C_INT),          intent(in)  :: j
    real(kind=C_DOUBLE),          intent(out) :: x
    real(kind=C_DOUBLE),          intent(out) :: y
    real(kind=C_DOUBLE),          intent(out) :: xp
    real(kind=C_DOUBLE),          intent(out) :: yp
    real(kind=C_DOUBLE),          intent(out) :: e
    integer(kind=C_INT32_T),      intent(out) :: pdgid
    real(kind=C_DOUBLE),          intent(out) :: m
    integer(kind=C_INT16_T),      intent(out) :: z
    integer(kind=C_INT16_T),      intent(out) :: a
    integer(kind=C_INT16_T),      intent(out) :: q
    integer(kind=C_INT),          intent(out) :: part_hit
    integer(kind=C_INT),          intent(out) :: part_abs
    real(kind=C_DOUBLE),          intent(out) :: part_impact
    real(kind=C_DOUBLE),          intent(out) :: part_indiv
    real(kind=C_DOUBLE),          intent(out) :: part_linteract
  end subroutine g4_collimate_return

!void g4_get_particle_count_(int* g4_npart)
  subroutine g4_get_particle_count(g4_npart) &
& bind(C,name="g4_get_particle_count")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(kind=C_INT),          intent(in) :: g4_npart
  end subroutine g4_get_particle_count

!void g4_collimation_clear_()
  subroutine g4_collimation_clear() &
& bind(C,name="g4_collimation_clear")
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine g4_collimation_clear

!void g4_keep_id_(int* id)
  subroutine g4_keep_id(id) &
& bind(C,name="g4_keep_id")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(kind=C_INT),          intent(in) :: id
  end subroutine g4_keep_id

 end interface

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
  else if(lnSplit(1) == 'RANGECUT_MM') then
    call chr_cast(lnSplit(2),g4_rangecut_mm,cErr)

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

