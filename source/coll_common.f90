! ============================================================================ !
!  Collimation Common Variables
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ============================================================================ !
module coll_common

  use floatPrecision
  use numerical_constants, only : zero

  implicit none

  integer, parameter :: max_ncoll  = 100
  integer, parameter :: nmat       = 14  ! Materials
  integer, parameter :: nrmat      = 12  ! Materials

  ! Logical Flags
  logical, save :: dowrite_impact    = .false.
  logical, save :: dowrite_dist      = .false.
  logical, save :: dowrite_secondary = .false.
  logical, save :: dowrite_amplitude = .false.
  logical, save :: dowritetracks     = .false.

  ! Various Variables
  integer, save :: rnd_seed   = 0

  ! Collimation Particle Arrays
  real(kind=fPrec), allocatable, save :: rcx(:)
  real(kind=fPrec), allocatable, save :: rcxp(:)
  real(kind=fPrec), allocatable, save :: rcy(:)
  real(kind=fPrec), allocatable, save :: rcyp(:)
  real(kind=fPrec), allocatable, save :: rcp(:)
  real(kind=fPrec), allocatable, save :: rcs(:)

  ! Collimator Materials
  character(4), private, save :: mname(nmat) = &
    ["Be  ","Al  ","Cu  ","W   ","Pb  ","C   ","C2  ","MoGR","CuCD","Mo  ","Glid","Iner","vacu","blac"]

  ! Pencil Beam
  integer,          save :: ipencil       = 0
  integer,          save :: pencil_distr  = 0
  real(kind=fPrec), save :: pencil_offset = zero
  real(kind=fPrec), save :: pencil_rmsx   = zero
  real(kind=fPrec), save :: pencil_rmsy   = zero
  real(kind=fPrec), save :: xp_pencil0    = zero
  real(kind=fPrec), save :: yp_pencil0    = zero
  real(kind=fPrec), allocatable, save :: x_pencil(:)
  real(kind=fPrec), allocatable, save :: y_pencil(:)
  real(kind=fPrec), allocatable, save :: pencil_dx(:)

  ! Output File Names
  character(len=12), parameter :: coll_survivalFile   = "survival.dat"
  character(len=12), parameter :: coll_gapsFile       = "collgaps.dat"
  character(len=10), parameter :: coll_impactFile     = "impact.dat"
  character(len=11), parameter :: coll_tracksFile     = "tracks2.dat"
  character(len=17), parameter :: coll_positionsFile  = "CollPositions.dat"
  character(len=20), parameter :: coll_pencilFile     = "pencilbeam_distr.dat"
  character(len=16), parameter :: coll_ellipseFile    = "coll_ellipse.dat"
  character(len=15), parameter :: coll_allImpactFile  = "all_impacts.dat"
  character(len=19), parameter :: coll_allAbsorbFile  = "all_absorptions.dat"
  character(len=16), parameter :: coll_scatterFile    = "Coll_Scatter.dat"
  character(len=16), parameter :: coll_fstImpactFile  = "FirstImpacts.dat"
  character(len=17), parameter :: coll_flukImpFile    = "FLUKA_impacts.dat"
  character(len=21), parameter :: coll_flukImpAllFile = "FLUKA_impacts_all.dat"
  character(len=13), parameter :: coll_twissLikeFile  = "twisslike.out"
  character(len=17), parameter :: coll_sigmaSetFile   = "sigmasettings.out"
  character(len=16), parameter :: coll_settingsFile   = "collsettings.dat"
  character(len=16), parameter :: coll_jawProfileFile = "jaw_profiles.dat"
  character(len=13), parameter :: coll_ampFile        = "amplitude.dat"
  character(len=17), parameter :: coll_orbitCheckFile = "orbitchecking.dat"
  character(len=16), parameter :: coll_summaryFile    = "coll_summary.dat"
  character(len=14), parameter :: coll_efficFile      = "efficiency.dat"
  character(len=19), parameter :: coll_efficDPFile    = "efficiency_dpop.dat"
  character(len=17), parameter :: coll_effic2DFile    = "efficiency_2d.dat"

  ! Output File Units
  integer, save :: outlun              = -1
  integer, save :: coll_survivalUnit   = -1
  integer, save :: coll_gapsUnit       = -1
  integer, save :: coll_impactUnit     = -1
  integer, save :: coll_tracksUnit     = -1
  integer, save :: coll_positionsUnit  = -1
  integer, save :: coll_pencilUnit     = -1
  integer, save :: coll_ellipseUnit    = -1
  integer, save :: coll_allImpactUnit  = -1
  integer, save :: coll_allAbsorbUnit  = -1
  integer, save :: coll_scatterUnit    = -1
  integer, save :: coll_fstImpactUnit  = -1
  integer, save :: coll_flukImpUnit    = -1
  integer, save :: coll_flukImpAllUnit = -1
  integer, save :: coll_twissLikeUnit  = -1
  integer, save :: coll_sigmaSetUnit   = -1
  integer, save :: coll_settingsUnit   = -1
  integer, save :: coll_jawProfileUnit = -1
  integer, save :: coll_ampUnit        = -1
  integer, save :: coll_orbitCheckUnit = -1
  integer, save :: coll_summaryUnit    = -1
  integer, save :: coll_efficUnit      = -1
  integer, save :: coll_efficDPUnit    = -1
  integer, save :: coll_effic2DUnit    = -1

contains

subroutine coll_expandArrays(npart_new, nblz_new)

  use mod_alloc
  use numerical_constants

  integer, intent(in) :: npart_new
  integer, intent(in) :: nblz_new

  call alloc(rcx,  npart_new, zero, "rcx")
  call alloc(rcxp, npart_new, zero, "rcxp")
  call alloc(rcy,  npart_new, zero, "rcy")
  call alloc(rcyp, npart_new, zero, "rcyp")
  call alloc(rcp,  npart_new, zero, "rcp")
  call alloc(rcs,  npart_new, zero, "rcs")

  call alloc(x_pencil,  max_ncoll, zero, "x_pencil")
  call alloc(y_pencil,  max_ncoll, zero, "y_pencil")
  call alloc(pencil_dx, max_ncoll, zero, "pencil_dx")

end subroutine coll_expandArrays

end module coll_common
