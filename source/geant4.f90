module geant4

  use floatPrecision
  use numerical_constants, only: zero
  use parpro, only : mNameLen
  use, intrinsic :: iso_c_binding

  implicit none

  real(kind=fPrec) :: g4_recut       = zero
  real(kind=fPrec) :: g4_aecut       = zero
  real(kind=fPrec) :: g4_rcut        = zero
  real(kind=fPrec) :: g4_rangecut_mm = zero

! Reference particle velocity
  real(kind=fPrec), save :: g4_v0    = zero
  integer, save          :: g4_keep

  character(len=64) :: g4_phys_str   = ' '
  character(len=64) :: g4_return_str = ' '

  logical(kind=C_BOOL) :: g4_enabled     = .false.
  logical(kind=C_BOOL) :: g4_debug       = .false.
  logical(kind=C_BOOL) :: g4_neutral     = .false. !Keep neutral particles or not
  logical(kind=C_BOOL) :: g4_keep_stable = .false.


! Energy deposition configuration
! Enabled or not?
  logical(kind=C_BOOL) :: g4_edep        = .false.

! Name of an entry (collimator name or ALL)
  character(len=mNameLen) :: g4_edep_colname

! Select either 2d or 3d histogram
  integer :: g4_hist_type

! High resolution energy deposition area
  real(kind=fPrec) :: g4_xstep
  real(kind=fPrec) :: g4_ystep
  real(kind=fPrec) :: g4_zstep
  real(kind=fPrec) :: g4_xmax
  real(kind=fPrec) :: g4_ymax
  real(kind=fPrec) :: g4_zmax

! Outer low resolution regions (if requested)
  real(kind=fPrec) :: g4_xbigstep
  real(kind=fPrec) :: g4_ybigstep
  real(kind=fPrec) :: g4_zbigstep
  real(kind=fPrec) :: g4_xbigmax
  real(kind=fPrec) :: g4_ybigmax
  real(kind=fPrec) :: g4_zbigmax

 interface

!void g4_collimation_init_(double* ReferenceE, int* seed, double* recut, double* aecut, double* rcut,
!int* PhysicsSelect, bool* g4_debug, bool* g4_keep_stable)
  subroutine g4_collimation_init(ReferenceE, seed, recut, aecut, rcut, rangecut_mm, v0, PhysicsSelect, g4_debug, g4_keep_stable, &
& g4_edep, g4_neutral) &
& bind(C,name="g4_collimation_init")
    use, intrinsic :: iso_c_binding
    implicit none
    real(kind=C_DOUBLE),          intent(in) :: ReferenceE
    integer(kind=C_INT),          intent(in) :: seed
    real(kind=C_DOUBLE),          intent(in) :: recut
    real(kind=C_DOUBLE),          intent(in) :: aecut
    real(kind=C_DOUBLE),          intent(in) :: rcut
    real(kind=C_DOUBLE),          intent(in) :: rangecut_mm
    real(kind=C_DOUBLE),          intent(in) :: v0
    character(kind=C_CHAR,len=1), intent(in) :: PhysicsSelect
!    integer(kind=C_INT),          intent(in) :: PhysicsSelect
    logical(kind=C_BOOL),         intent(in) :: g4_debug
    logical(kind=C_BOOL),         intent(in) :: g4_keep_stable
    logical(kind=C_BOOL),         intent(in) :: g4_edep
    logical(kind=C_BOOL),         intent(in) :: g4_neutral
  end subroutine g4_collimation_init

!void g4_add_collimator_(char* name, char* material, double* length, double* aperture, double* rotation, double* offset)
  subroutine g4_add_collimator(name, material, length, aperture, rotation, x_offset, y_offset, onesided) &
& bind(C,name="g4_add_collimator")
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=C_CHAR,len=1), intent(in) :: name
    character(kind=C_CHAR,len=1), intent(in) :: material
    real(kind=C_DOUBLE),          intent(in) :: length
    real(kind=C_DOUBLE),          intent(in) :: aperture
    real(kind=C_DOUBLE),          intent(in) :: rotation
    real(kind=C_DOUBLE),          intent(in) :: x_offset
    real(kind=C_DOUBLE),          intent(in) :: y_offset
    logical(kind=C_BOOL),         intent(in) :: onesided
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
  subroutine g4_add_particle(x, y, xp, yp, e, pdgid, nzz, naa, nqq, mass, sigmv, partID, parentID, weight, spin_x, spin_y, spin_z)&
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
    real(kind=C_DOUBLE),          intent(in) :: sigmv
    integer(kind=C_INT32_T),      intent(in) :: partID
    integer(kind=C_INT32_T),      intent(in) :: parentID
    real(kind=C_DOUBLE),          intent(in) :: weight
    real(kind=C_DOUBLE),          intent(in) :: spin_x
    real(kind=C_DOUBLE),          intent(in) :: spin_y
    real(kind=C_DOUBLE),          intent(in) :: spin_z
  end subroutine g4_add_particle

!void g4_collimate_()
  subroutine g4_collimate() &
& bind(C,name="g4_collimate")
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine g4_collimate

!void g4_collimate_return_(int* j, double* x, double* y, double* xp, double* yp, double* e, int32_t* pdgid, double* m,
! int16_t* z, int16_t* a, int16_t* q, int *part_hit, int *part_abs, double *part_impact, double *part_indiv, double *part_linteract)
  subroutine g4_collimate_return(j, x, y, xp, yp, e, pdgid, m, z, a, q, sigmv, partID, parentID, weight, &
& part_hit, part_abs, part_impact, part_indiv, &
& part_linteract, spin_x, spin_y, spin_z) &
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
    real(kind=C_DOUBLE),          intent(out) :: sigmv
    integer(kind=C_INT32_T),      intent(out) :: partID
    integer(kind=C_INT32_T),      intent(out) :: parentID
    real(kind=C_DOUBLE),          intent(out) :: weight
    integer(kind=C_INT),          intent(out) :: part_hit
    integer(kind=C_INT),          intent(out) :: part_abs
    real(kind=C_DOUBLE),          intent(out) :: part_impact
    real(kind=C_DOUBLE),          intent(out) :: part_indiv
    real(kind=C_DOUBLE),          intent(out) :: part_linteract
    real(kind=C_DOUBLE),          intent(out) :: spin_x
    real(kind=C_DOUBLE),          intent(out) :: spin_y
    real(kind=C_DOUBLE),          intent(out) :: spin_z
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

!extern "C" void g4_add_edep(char* name_in, int* type, double* xstep, double* ystep, double* zstep, double* xmax, double* ymax, double* zmax, double* xbigstep, double* ybigstep, double* zbigstep, double* xbigmax, double* ybigmax, double* zbigmax)
  subroutine g4_add_edep(name_in, edtype, xstep, ystep, zstep, xmax, ymax, zmax, &
  xbigstep, ybigstep, zbigstep, xbigmax, ybigmax, zbigmax) &
& bind(C,name="g4_add_edep")
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=C_CHAR,len=1), intent(in) :: name_in
    integer(kind=C_INT),          intent(in) :: edtype
    real(kind=C_DOUBLE),          intent(in) :: xstep
    real(kind=C_DOUBLE),          intent(in) :: ystep
    real(kind=C_DOUBLE),          intent(in) :: zstep
    real(kind=C_DOUBLE),          intent(in) :: xmax
    real(kind=C_DOUBLE),          intent(in) :: ymax
    real(kind=C_DOUBLE),          intent(in) :: zmax
    real(kind=C_DOUBLE),          intent(in) :: xbigstep
    real(kind=C_DOUBLE),          intent(in) :: ybigstep
    real(kind=C_DOUBLE),          intent(in) :: zbigstep
    real(kind=C_DOUBLE),          intent(in) :: xbigmax
    real(kind=C_DOUBLE),          intent(in) :: ybigmax
    real(kind=C_DOUBLE),          intent(in) :: zbigmax
  end subroutine g4_add_edep

  subroutine g4_get_maximum_particle_id(m_id) bind(C,name="g4_get_maximum_particle_id")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(kind=C_INT32_T),        intent(out) :: m_id
  end subroutine g4_get_maximum_particle_id

  subroutine g4_set_maximum_particle_id(m_id) bind(C,name="g4_set_maximum_particle_id")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(kind=C_INT32_T),        intent(in)  :: m_id
  end subroutine g4_set_maximum_particle_id

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

  if(nSplit /= 2 .and. nSplit /= 4) then
    write(lerr,"(a,i0)") "GEANT4> ERROR Expected 2 or 4 entries per line, got ",nSplit
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
    g4_return_str = trim(lnSplit(2))
    if(g4_return_str .eq. 'STABLE') then
      g4_keep_stable = .true.
    else if(g4_return_str .eq. 'IONS') then
      g4_keep = -1
      call g4_keep_id(g4_keep)
    else
      call chr_cast(lnSplit(2),g4_keep,cErr)
      call g4_keep_id(g4_keep)
    end if

  else if(lnSplit(1) == 'NEUTRAL') then
    g4_neutral = .true.

!Physics list to use string
  else if(lnSplit(1) == 'PHYSICS') then
    g4_phys_str = trim(lnSplit(2)) // C_NULL_CHAR

!Energy deposition configuration
!EDEP COLNAME SETTING VALUE
!HIST2D, HIST3D
!RANGES, STEPSIZE
  else if(lnSplit(1) == 'EDEP') then
    g4_edep = .true.

!Reset these to zero each time
    g4_hist_type = 0

! High resolution energy deposition area
    g4_xstep = zero
    g4_ystep = zero
    g4_zstep = zero
    g4_xmax  = zero
    g4_ymax  = zero
    g4_zmax  = zero

! Outer low resolution regions (if requested)
    g4_xbigstep = zero
    g4_ybigstep = zero
    g4_zbigstep = zero
    g4_xbigmax  = zero
    g4_ybigmax  = zero
    g4_zbigmax  = zero

!   Get the collimator name
!   Use special value 'ALL' for everything
    g4_edep_colname = trim(lnSplit(2)) //C_NULL_CHAR

!   Now look at the configuration variables
    if(lnSplit(3) == 'HIST2D') then
    else if(lnSplit(3) == 'HIST2D') then
      g4_hist_type = 2
    else if(lnSplit(3) == 'HIST3D') then
      g4_hist_type = 3
    else if(lnSplit(3) == 'XSTEP') then
      call chr_cast(lnSplit(4),g4_xstep,cErr)
    else if(lnSplit(3) == 'YSTEP') then
      call chr_cast(lnSplit(4),g4_ystep,cErr)
    else if(lnSplit(3) == 'ZSTEP') then
      call chr_cast(lnSplit(4),g4_zstep,cErr)
    else if(lnSplit(3) == 'XMAX') then
      call chr_cast(lnSplit(4),g4_xmax,cErr)
    else if(lnSplit(3) == 'YMAX') then
      call chr_cast(lnSplit(4),g4_ymax,cErr)
    else if(lnSplit(3) == 'ZMAX') then
      call chr_cast(lnSplit(4),g4_zmax,cErr)
    else if(lnSplit(3) == 'XBIGSTEP') then
      call chr_cast(lnSplit(4),g4_xbigstep,cErr)
    else if(lnSplit(3) == 'YBIGSTEP') then
      call chr_cast(lnSplit(4),g4_ybigstep,cErr)
    else if(lnSplit(3) == 'ZBIGSTEP') then
      call chr_cast(lnSplit(4),g4_zbigstep,cErr)
    else if(lnSplit(3) == 'XBIGMAX') then
      call chr_cast(lnSplit(4),g4_xbigmax,cErr)
    else if(lnSplit(3) == 'YBIGMAX') then
      call chr_cast(lnSplit(4),g4_ybigmax,cErr)
    else if(lnSplit(3) == 'ZBIGMAX') then
      call chr_cast(lnSplit(4),g4_zbigmax,cErr)
    else
      write(lerr,"(2a)") "GEANT4> ERROR unknown keyword", lnSplit(3)
      iErr = .true.
    end if

    call g4_add_edep(g4_edep_colname, g4_hist_type, g4_xstep, g4_ystep, g4_zstep, g4_xmax, g4_ymax, g4_zmax, &
&                    g4_xbigstep, g4_ybigstep, g4_zbigstep, g4_xbigmax, g4_ybigmax, g4_zbigmax)
  else
    write(lerr,"(2a)") "GEANT4> ERROR unknown keyword", lnSplit(1)
    iErr = .true.
  end if

!Check configuration

!check + enable flags
end subroutine geant4_parseInputLine

subroutine geant4_parseInputDone

  use crcoall

#ifdef ROOT
  use root_output
#endif

  implicit none

  logical have_root

  !GEANT4 is enabled
  g4_enabled = .true.

! Is there root support compiled in?
  have_root = .false.
#ifdef ROOT
  have_root = .true.
#endif

  if(g4_edep .eqv. .true.) then
!   Check we have root compiled
    if(have_root .eqv. .true.) then
!     Check we have root enabled
#ifdef ROOT
      if(root_flag .eqv. .false.) then
!     explode and exit
        write(lerr,'(a)') 'GEANT4> ERROR: Energy deposition studies require ROOT output compiled in and enabled!'
        call prror
      end if
#endif
    else
!   No root support compiled in
!   explode and exit
      write(lerr,'(a)') 'GEANT4> ERROR: Energy deposition studies require ROOT output compiled in and enabled!'
      call prror
    end if
  end if

  !If we have no physics string set, use FTFP_BERT as the default
  if(len_trim(g4_phys_str).eq.0) then
    write(lout,'(a)') 'GEANT4> WARNING: No physics model requested, defaulting to FTFP_BERT'
    g4_phys_str = 'FTFP_BERT' // C_NULL_CHAR
  end if

  write(lout,'(2a)') 'GEANT4> Using physics list: ', g4_phys_str

end subroutine geant4_parseInputDone

end module geant4

