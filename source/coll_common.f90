! ================================================================================================ !
!  Collimation Common Variables
! ================================================================================================ !
module coll_common

  use floatPrecision
  use numerical_constants, only : zero

  implicit none

  integer, parameter :: max_ncoll  = 100

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

! ================================================================================================ !
!  Collimator Material Variables
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Cross section inputs and material property database
! ================================================================================================ !
module coll_materials

  use floatPrecision
  use numerical_constants, only : zero, one, c1e10

  implicit none

  integer, parameter :: nmat  = 14 ! Total number of materials
  integer, parameter :: nrmat = 12 ! Number of real materials

  real(kind=fPrec), public, save :: csect(0:5,nmat)
  real(kind=fPrec), public, save :: xintl(nmat)
  real(kind=fPrec), public, save :: bn(nmat)
  real(kind=fPrec), public, save :: freep(nmat)
  real(kind=fPrec), public, save :: cgen(200,nmat)

  ! pp cross-sections and parameters for energy dependence
  real(kind=fPrec), parameter :: pptref = 0.04_fPrec
  real(kind=fPrec), parameter :: freeco = 1.618_fPrec

  ! Collimator Materials
  character(4), parameter :: colmats(nmat) = &
    ["BE  ","AL  ","CU  ","W   ","PB  ","C   ","C2  ","MoGR","CuCD","Mo  ","Glid","Iner","VA  ","BL  "]

  ! Mean excitation energy (GeV) values added by Claudia for Bethe-Bloch implementation:
  real(kind=fPrec), public, save :: exenergy(nmat) = [ &
    63.7e-9_fPrec, 166.0e-9_fPrec, 322.0e-9_fPrec, 727.0e-9_fPrec, 823.0e-9_fPrec, 78.0e-9_fPrec, 78.0e-9_fPrec, &
    87.1e-9_fPrec, 152.9e-9_fPrec, 424.0e-9_fPrec, 320.8e-9_fPrec, 682.2e-9_fPrec, zero, c1e10 ]

  ! GRD IMPLEMENT CHANGES FROM JBJ, 2/2003 RWA
  real(kind=fPrec), public, save :: anuc(nmat)  = &
    [ 9.01_fPrec,  26.98_fPrec,  63.55_fPrec, 183.85_fPrec, 207.19_fPrec,    12.01_fPrec,  12.01_fPrec,  &
     13.53_fPrec,  25.24_fPrec,  95.96_fPrec,  63.15_fPrec, 166.70_fPrec,     zero,         zero         ]
  real(kind=fPrec), public, save :: zatom(nmat) = &
    [ 4.00_fPrec,  13.00_fPrec,  29.00_fPrec,  74.00_fPrec,  82.00_fPrec,     6.00_fPrec,   6.00_fPrec,  &
      6.65_fPrec,  11.90_fPrec,  42.00_fPrec,  28.80_fPrec,  67.70_fPrec,     zero,         zero         ]
  real(kind=fPrec), public, save :: rho(nmat)   = &
    [ 1.848_fPrec,  2.70_fPrec,   8.96_fPrec,  19.30_fPrec,   11.35_fPrec,    1.67_fPrec,   4.52_fPrec,  &
      2.500_fPrec,  5.40_fPrec,  10.22_fPrec,   8.93_fPrec,   18.00_fPrec,    zero,         zero         ]
  real(kind=fPrec), public, save :: emr(nmat)   = &
    [ 0.22_fPrec,   0.302_fPrec,  0.366_fPrec,  0.520_fPrec,   0.542_fPrec,   0.25_fPrec,   0.25_fPrec,  &
      0.25_fPrec,   0.308_fPrec,  0.481_fPrec,  0.418_fPrec,   0.578_fPrec,   zero,         zero         ]
  real(kind=fPrec), public, save :: hcut(nmat)  = &
    [ 0.02_fPrec,   0.02_fPrec,   0.01_fPrec,   0.01_fPrec,    0.01_fPrec,    0.02_fPrec,    0.02_fPrec, &
      0.02_fPrec,   0.02_fPrec,   0.02_fPrec,   0.02_fPrec,    0.02_fPrec,    zero,          zero        ]
  real(kind=fPrec), public, save :: radl(nmat)  = &
    [ 0.353_fPrec,  0.089_fPrec,  0.0143_fPrec, 0.0035_fPrec,  0.0056_fPrec,  0.2557_fPrec, 0.094_fPrec, &
      0.1193_fPrec, 0.0316_fPrec, 0.0096_fPrec, 0.0144_fPrec,  0.00385_fPrec, 1.0e12_fPrec, 1.0e12_fPrec ]

  ! Nuclear elastic slope from Schiz et al.,PRD 21(3010)1980
  ! MAY06-GRD value for Tungsten (W) not stated. Last 2 ones interpolated
  real(kind=fPrec), public, save :: bnref(nmat) = &
    [74.7_fPrec, 120.3_fPrec, 217.8_fPrec, 440.3_fPrec, 455.3_fPrec, 70.0_fPrec, 70.0_fPrec, &
     76.7_fPrec, 115.0_fPrec, 273.9_fPrec, 208.7_fPrec, 392.1_fPrec, zero,       zero        ]

  ! All cross-sections are in barns,nuclear values from RPP at 20geV
  ! Coulomb is integerated above t=tLcut[Gev2] (+-1% out Gauss mcs)

  ! in Cs and CsRef,1st index: Cross-sections for processes
  ! 0:Total, 1:absorption, 2:nuclear elastic, 3:pp or pn elastic
  ! 4:Single Diffractive pp or pn, 5:Coulomb for t above mcs

  ! Claudia 2013: updated cross section values. Unit: Barn. New 2013:
  real(kind=fPrec), public, save :: csref(0:5,nmat)
  data csref(0,1), csref(1,1), csref(5,1) /0.271_fPrec, 0.192_fPrec, 0.0035e-2_fPrec/
  data csref(0,2), csref(1,2), csref(5,2) /0.643_fPrec, 0.418_fPrec, 0.0340e-2_fPrec/
  data csref(0,3), csref(1,3), csref(5,3) /1.253_fPrec, 0.769_fPrec, 0.1530e-2_fPrec/
  data csref(0,4), csref(1,4), csref(5,4) /2.765_fPrec, 1.591_fPrec, 0.7680e-2_fPrec/
  data csref(0,5), csref(1,5), csref(5,5) /3.016_fPrec, 1.724_fPrec, 0.9070e-2_fPrec/
  data csref(0,6), csref(1,6), csref(5,6) /0.337_fPrec, 0.232_fPrec, 0.0076e-2_fPrec/
  data csref(0,7), csref(1,7), csref(5,7) /0.337_fPrec, 0.232_fPrec, 0.0076e-2_fPrec/
  data csref(0,8), csref(1,8), csref(5,8) /0.362_fPrec, 0.247_fPrec, 0.0094e-2_fPrec/
  data csref(0,9), csref(1,9), csref(5,9) /0.572_fPrec, 0.370_fPrec, 0.0279e-2_fPrec/
  data csref(0,10),csref(1,10),csref(5,10)/1.713_fPrec, 1.023_fPrec, 0.2650e-2_fPrec/
  data csref(0,11),csref(1,11),csref(5,11)/1.246_fPrec, 0.765_fPrec, 0.1390e-2_fPrec/
  data csref(0,12),csref(1,12),csref(5,12)/2.548_fPrec, 1.473_fPrec, 0.5740e-2_fPrec/

  ! Cprob to choose an interaction in iChoix
  real(kind=fPrec), public, save :: cprob(0:5,nmat)
  data cprob(0,1:nmat)/nmat*zero/
  data cprob(5,1:nmat)/nmat*one/

  ! Electron density and plasma energy
  real(kind=fPrec), public, save :: edens(nmat)
  real(kind=fPrec), public, save :: pleng(nmat)

end module coll_materials
