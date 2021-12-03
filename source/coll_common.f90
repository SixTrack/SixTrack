! ================================================================================================ !
!  Collimation Common Variables
! ================================================================================================ !
module coll_common

  use floatPrecision
  use numerical_constants, only : zero

  implicit none

  ! Logical Flags
  logical, save :: coll_debug         = .true.
  logical, save :: dowrite_impact     = .false.
  logical, save :: dowrite_dist       = .false.
  logical, save :: dowrite_secondary  = .false.
  logical, save :: dowrite_amplitude  = .false.
  logical, save :: dowrite_tracks     = .false.
  logical, save :: dowrite_efficiency = .false.
  logical, save :: dowrite_crycoord   = .false.
  logical, save :: coll_hasCrystal    = .false.

  ! Various Variables
  integer, save :: rnd_seed   = 0

  ! Collimation Particle Arrays
  real(kind=fPrec), allocatable, save :: rcx(:)
  real(kind=fPrec), allocatable, save :: rcxp(:)
  real(kind=fPrec), allocatable, save :: rcy(:)
  real(kind=fPrec), allocatable, save :: rcyp(:)
  real(kind=fPrec), allocatable, save :: rcp(:)
  real(kind=fPrec), allocatable, save :: rcs(:)

  ! Process index for interaction with crystals
  integer, allocatable, save :: cry_proc(:)
  integer, allocatable, save :: cry_proc_prev(:)
  integer, allocatable, save :: cry_proc_tmp(:)

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
  real(kind=fPrec), allocatable, save :: xp_pencil(:)
  real(kind=fPrec), allocatable, save :: yp_pencil(:)
  real(kind=fPrec), allocatable, save :: pencil_dx(:)

  ! Other Arrays
  integer,          allocatable, save :: cn_impact(:)
  integer,          allocatable, save :: cn_absorbed(:)
  real(kind=fPrec), allocatable, save :: caverage(:)
  real(kind=fPrec), allocatable, save :: csigma(:)
  real(kind=fPrec), allocatable, save :: gap_rms_error(:)
  real(kind=fPrec), allocatable, save :: csum(:)
  real(kind=fPrec), allocatable, save :: csqsum(:)

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
  character(len=17), parameter :: coll_sigmaSetFile   = "sigmasettings.out"
  character(len=16), parameter :: coll_settingsFile   = "collsettings.dat"
  character(len=16), parameter :: coll_jawProfileFile = "jaw_profiles.dat"
  character(len=13), parameter :: coll_ampFile        = "amplitude.dat"
  character(len=17), parameter :: coll_orbitCheckFile = "orbitchecking.dat"
  character(len=16), parameter :: coll_summaryFile    = "coll_summary.dat"
  character(len=14), parameter :: coll_efficFile      = "efficiency.dat"
  character(len=19), parameter :: coll_efficDPFile    = "efficiency_dpop.dat"
  character(len=17), parameter :: coll_effic2DFile    = "efficiency_2d.dat"
  character(len=16), parameter :: coll_cryEntFile     = "cry_entrance.dat"
  character(len=12), parameter :: coll_cryExitFile    = "cry_exit.dat"
  character(len=19), parameter :: coll_cryInterFile   = "cry_interaction.dat"

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
  integer, save :: coll_sigmaSetUnit   = -1
  integer, save :: coll_settingsUnit   = -1
  integer, save :: coll_jawProfileUnit = -1
  integer, save :: coll_ampUnit        = -1
  integer, save :: coll_orbitCheckUnit = -1
  integer, save :: coll_summaryUnit    = -1
  integer, save :: coll_efficUnit      = -1
  integer, save :: coll_efficDPUnit    = -1
  integer, save :: coll_effic2DUnit    = -1
  integer, save :: coll_cryEntUnit     = -1
  integer, save :: coll_cryExitUnit    = -1
  integer, save :: coll_cryInterUnit   = -1

#ifdef HDF5
  ! Variables to save hdf5 dataset indices
  integer, save :: coll_hdf5_survival
  integer, save :: coll_hdf5_allImpacts
  integer, save :: coll_hdf5_fstImpacts
  integer, save :: coll_hdf5_allAbsorb
  integer, save :: coll_hdf5_collScatter
#endif

#ifdef CR
  ! For resetting file positions
  integer, public,  save :: fort208Pos                = -1
  integer, public,  save :: fort208Pos_CR             =  0
  integer, public,  save :: coll_survivalFilePos      = -1
  integer, public,  save :: coll_survivalFilePos_CR   =  0
  integer, public,  save :: coll_gapsFilePos          = -1
  integer, public,  save :: coll_gapsFilePos_CR       =  0
  integer, public,  save :: coll_settingsFilePos      = -1
  integer, public,  save :: coll_settingsFilePos_CR   =  0
  integer, public,  save :: coll_positionsFilePos     = -1
  integer, public,  save :: coll_positionsFilePos_CR  =  0
  integer, public,  save :: coll_tracksFilePos        = -1
  integer, public,  save :: coll_tracksFilePos_CR     =  0
  integer, public,  save :: coll_pencilFilePos        = -1
  integer, public,  save :: coll_pencilFilePos_CR     =  0
  integer, public,  save :: coll_cryEntFilePos        = -1
  integer, public,  save :: coll_cryEntFilePos_CR     =  0
  integer, public,  save :: coll_cryExitFilePos       = -1
  integer, public,  save :: coll_cryExitFilePos_CR    =  0
  integer, public,  save :: coll_cryInterFilePos      = -1
  integer, public,  save :: coll_cryInterFilePos_CR   =  0
  integer, public,  save :: coll_ellipseFilePos       = -1
  integer, public,  save :: coll_ellipseFilePos_CR    =  0
  integer, public,  save :: coll_allImpactFilePos     = -1
  integer, public,  save :: coll_allImpactFilePos_CR  =  0
  integer, public,  save :: coll_allAbsorbFilePos     = -1
  integer, public,  save :: coll_allAbsorbFilePos_CR  =  0
  integer, public,  save :: coll_scatterFilePos       = -1
  integer, public,  save :: coll_scatterFilePos_CR    =  0
  integer, public,  save :: coll_fstImpactFilePos     = -1
  integer, public,  save :: coll_fstImpactFilePos_CR  =  0
  integer, public,  save :: coll_flukImpFilePos       = -1
  integer, public,  save :: coll_flukImpFilePos_CR    =  0
  integer, public,  save :: coll_flukImpAllFilePos    = -1
  integer, public,  save :: coll_flukImpAllFilePos_CR =  0
  integer, public,  save :: coll_jawProfileFilePos    = -1
  integer, public,  save :: coll_jawProfileFilePos_CR =  0
  integer, public,  save :: coll_efficFilePos         = -1
  integer, public,  save :: coll_efficFilePos_CR      =  0
  integer, public,  save :: coll_efficDPFilePos       = -1
  integer, public,  save :: coll_efficDPFilePos_CR    =  0
  integer, public,  save :: coll_effic2DFilePos       = -1
  integer, public,  save :: coll_effic2DFilePos_CR    =  0
  integer, public,  save :: coll_summaryFilePos       = -1
  integer, public,  save :: coll_summaryFilePos_CR    =  0
  integer, public,  save :: coll_ampFilePos           = -1
  integer, public,  save :: coll_ampFilePos_CR        =  0
  integer, public,  save :: coll_orbitCheckFilePos    = -1
  integer, public,  save :: coll_orbitCheckFilePos_CR =  0
  integer, public,  save :: coll_sigmaSetFilePos      = -1
  integer, public,  save :: coll_sigmaSetFilePos_CR   =  0
  integer, public,  save :: coll_impactFilePos        = -1
  integer, public,  save :: coll_impactFilePos_CR     =  0
  integer, public,  save :: coll_trackoutPos          = -1
  integer, public,  save :: coll_trackoutPos_CR       =  0

! For the RNG
  integer, public,  save :: lux_CR  = 0
  integer, public,  save :: seed_CR = 0
  integer, public,  save :: k1_CR   = 0
  integer, public,  save :: k2_CR   = 0


  integer,          allocatable, save :: cn_impact_cr(:)
  integer,          allocatable, save :: cn_absorbed_cr(:)
  real(kind=fPrec), allocatable, save :: caverage_cr(:)
  real(kind=fPrec), allocatable, save :: csigma_cr(:)
  real(kind=fPrec), allocatable, save :: gap_rms_error_cr(:)
  real(kind=fPrec), allocatable, save :: csum_cr(:)
  real(kind=fPrec), allocatable, save :: csqsum_cr(:)
#endif

contains

subroutine coll_expandArrays(npart_new)

  use mod_alloc
  use numerical_constants

  integer, intent(in) :: npart_new

  call alloc(rcx,  npart_new, zero, "rcx")
  call alloc(rcxp, npart_new, zero, "rcxp")
  call alloc(rcy,  npart_new, zero, "rcy")
  call alloc(rcyp, npart_new, zero, "rcyp")
  call alloc(rcp,  npart_new, zero, "rcp")
  call alloc(rcs,  npart_new, zero, "rcs")

  call alloc(cry_proc, npart_new, -1, "cry_proc")
  call alloc(cry_proc_prev, npart_new, -1, "cry_proc_prev")
  call alloc(cry_proc_tmp, npart_new, -1, "cry_proc_tmp")

end subroutine coll_expandArrays

subroutine coll_expandNColl(nColl)

  use mod_alloc
  use numerical_constants

  integer, intent(in) :: nColl

  call alloc(cn_impact,     nColl, 0,    "cn_impact")
  call alloc(cn_absorbed,   nColl, 0,    "cn_absorbed")
  call alloc(caverage,      nColl, zero, "caverage")
  call alloc(csigma,        nColl, zero, "csigma")
  call alloc(gap_rms_error, nColl, zero, "gap_rms_error")
  call alloc(csum,          nColl, zero, "csum")
  call alloc(csqsum,        nColl, zero, "csqsum")
  call alloc(x_pencil,      nColl, zero, "x_pencil")
  call alloc(y_pencil,      nColl, zero, "y_pencil")
  call alloc(xp_pencil,     nColl, zero, "xp_pencil")
  call alloc(yp_pencil,     nColl, zero, "yp_pencil")
  call alloc(pencil_dx,     nColl, zero, "pencil_dx")

end subroutine coll_expandNColl

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

  integer, parameter :: nmat  = 16 ! Total number of materials
  integer, parameter :: nrmat = 14 ! Number of real materials

  ! pp cross-sections and parameters for energy dependence
  real(kind=fPrec), parameter :: pptref = 0.04_fPrec
  real(kind=fPrec), parameter :: freeco = 1.618_fPrec

  ! Collimator Material Arrays
  character(4),     public, save :: colmats(nmat)   ! Material Names
  real(kind=fPrec), public, save :: exenergy(nmat)  ! Mean excitation energy [GeV]
  real(kind=fPrec), public, save :: anuc(nmat)      ! Atomic mass
  real(kind=fPrec), public, save :: zatom(nmat)     ! Atomic Z
  real(kind=fPrec), public, save :: rho(nmat)       ! Density
  real(kind=fPrec), public, save :: emr(nmat)       ! Nuclear radius
  real(kind=fPrec), public, save :: hcut(nmat)      ! T cut (upper)
  real(kind=fPrec), public, save :: radl(nmat)      ! Remaining length
  real(kind=fPrec), public, save :: bnref(nmat)     ! Nuclear elastic slope from Schiz et al., PRD 21(3010)1980
  real(kind=fPrec), public, save :: csect(0:5,nmat) ! Cross section
  real(kind=fPrec), public, save :: xintl(nmat)     ! Interaction length
  real(kind=fPrec), public, save :: bn(nmat)        ! Nuclear elastic related
  real(kind=fPrec), public, save :: freep(nmat)     ! Number of nucleons involved in single scattering
  real(kind=fPrec), public, save :: cgen(200,nmat)  ! Used by FUNLUX / Rutherford routine

  ! All cross-sections are in barns. Nuclear values from RPP at 20 GeV
  ! Coulomb is integerated above t=tLcut[Gev2] (+-1% out Gauss mcs)

  ! In Cs and CsRef,1st index: Cross-sections for processes
  ! 0:Total, 1:absorption, 2:nuclear elastic, 3:pp or pn elastic
  ! 4:Single Diffractive pp or pn, 5:Coulomb for t above mcs

  ! Claudia 2013: updated cross section values. Unit: Barn. New 2013:
  real(kind=fPrec), public, parameter :: csref(0:5,nmat) = reshape([ &
    [0.271_fPrec, 0.192_fPrec, zero, zero, zero, 0.0035e-2_fPrec], & ! BE
    [0.643_fPrec, 0.418_fPrec, zero, zero, zero, 0.0340e-2_fPrec], & ! AL
    [1.253_fPrec, 0.769_fPrec, zero, zero, zero, 0.1530e-2_fPrec], & ! CU
    [2.765_fPrec, 1.591_fPrec, zero, zero, zero, 0.7680e-2_fPrec], & ! W
    [3.016_fPrec, 1.724_fPrec, zero, zero, zero, 0.9070e-2_fPrec], & ! PB
    [0.337_fPrec, 0.232_fPrec, zero, zero, zero, 0.0076e-2_fPrec], & ! C
    [0.337_fPrec, 0.232_fPrec, zero, zero, zero, 0.0076e-2_fPrec], & ! C2
    [0.664_fPrec, 0.430_fPrec, zero, zero, zero, 0.0390e-2_fPrec], & ! Si
    [1.388_fPrec, 0.844_fPrec, zero, zero, zero, 0.1860e-2_fPrec], & ! Ge
    [0.362_fPrec, 0.247_fPrec, zero, zero, zero, 0.0094e-2_fPrec], & ! MoGR
    [0.572_fPrec, 0.370_fPrec, zero, zero, zero, 0.0279e-2_fPrec], & ! CuCD
    [1.713_fPrec, 1.023_fPrec, zero, zero, zero, 0.2650e-2_fPrec], & ! Mo
    [1.246_fPrec, 0.765_fPrec, zero, zero, zero, 0.1390e-2_fPrec], & ! Glid
    [2.548_fPrec, 1.473_fPrec, zero, zero, zero, 0.5740e-2_fPrec], & ! Iner
    [       zero,        zero, zero, zero, zero,            zero], & ! VA
    [       zero,        zero, zero, zero, zero,            zero]  & ! BL
  ], shape=[6,nmat])

  ! Cprob to choose an interaction in iChoix
  real(kind=fPrec), public, save :: cprob(0:5,nmat) = reshape([ &
    [zero, zero, zero, zero, zero, one], & ! BE
    [zero, zero, zero, zero, zero, one], & ! AL
    [zero, zero, zero, zero, zero, one], & ! CU
    [zero, zero, zero, zero, zero, one], & ! W
    [zero, zero, zero, zero, zero, one], & ! PB
    [zero, zero, zero, zero, zero, one], & ! C
    [zero, zero, zero, zero, zero, one], & ! C2
    [zero, zero, zero, zero, zero, one], & ! Si
    [zero, zero, zero, zero, zero, one], & ! Ge
    [zero, zero, zero, zero, zero, one], & ! MoGR
    [zero, zero, zero, zero, zero, one], & ! CuCD
    [zero, zero, zero, zero, zero, one], & ! Mo
    [zero, zero, zero, zero, zero, one], & ! Glid
    [zero, zero, zero, zero, zero, one], & ! Iner
    [zero, zero, zero, zero, zero, one], & ! VA
    [zero, zero, zero, zero, zero, one]  & ! BL
  ], shape=[6,nmat])

  ! Electron density and plasma energy
  real(kind=fPrec), public, save :: edens(nmat) = zero
  real(kind=fPrec), public, save :: pleng(nmat) = zero

contains

! ================================================================================================ !
!  Initialise Material Database Values
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-19
!  Updated: 2019-09-19
!  Rewritten to make it easier to add new materials.
! ================================================================================================ !
subroutine collmat_init

  use crcoall

  integer iMat

  iMat = 1
  colmats(iMat)  = "BE"
  exenergy(iMat) =  63.7e-9_fPrec
  anuc(iMat)     =   9.01_fPrec
  zatom(iMat)    =   4.00_fPrec
  rho(iMat)      =   1.848_fPrec
  emr(iMat)      =   0.22_fPrec
  hcut(iMat)     =   0.02_fPrec
  radl(iMat)     =   0.353_fPrec
  bnref(iMat)    =  74.7_fPrec

  iMat = iMat + 1
  colmats(iMat)  = "AL"
  exenergy(iMat) = 166.0e-9_fPrec
  anuc(iMat)     =  26.98_fPrec
  zatom(iMat)    =  13.00_fPrec
  rho(iMat)      =   2.70_fPrec
  emr(iMat)      =   0.302_fPrec
  hcut(iMat)     =   0.02_fPrec
  radl(iMat)     =   0.089_fPrec
  bnref(iMat)    = 120.3_fPrec

  ! Copper
  iMat = iMat + 1
  colmats(iMat)  = "CU"
  exenergy(iMat) = 322.0e-9_fPrec
  anuc(iMat)     =  63.55_fPrec
  zatom(iMat)    =  29.00_fPrec
  rho(iMat)      =   8.96_fPrec
  emr(iMat)      =   0.366_fPrec
  hcut(iMat)     =   0.01_fPrec
  radl(iMat)     =   0.0143_fPrec
  bnref(iMat)    = 217.8_fPrec

  ! Tungsten
  iMat = iMat + 1
  colmats(iMat)  = "W"
  exenergy(iMat) = 727.0e-9_fPrec
  anuc(iMat)     = 183.85_fPrec
  zatom(iMat)    =  74.00_fPrec
  rho(iMat)      =  19.30_fPrec
  emr(iMat)      =   0.520_fPrec
  hcut(iMat)     =   0.01_fPrec
  radl(iMat)     =   0.0035_fPrec
  bnref(iMat)    = 440.3_fPrec

  ! Lead
  iMat = iMat + 1
  colmats(iMat)  = "PB"
  exenergy(iMat) = 823.0e-9_fPrec
  anuc(iMat)     = 207.19_fPrec
  zatom(iMat)    =  82.00_fPrec
  rho(iMat)      =  11.35_fPrec
  emr(iMat)      =   0.542_fPrec
  hcut(iMat)     =   0.01_fPrec
  radl(iMat)     =   0.0056_fPrec
  bnref(iMat)    = 455.3_fPrec

  ! Carbon
  iMat = iMat + 1
  colmats(iMat)  = "C"
  exenergy(iMat) =  78.0e-9_fPrec
  anuc(iMat)     =  12.01_fPrec
  zatom(iMat)    =   6.00_fPrec
  rho(iMat)      =   1.67_fPrec
  emr(iMat)      =   0.25_fPrec
  hcut(iMat)     =   0.02_fPrec
  radl(iMat)     =   0.2557_fPrec
  bnref(iMat)    =  70.0_fPrec

  ! Carbon2
  iMat = iMat + 1
  colmats(iMat)  = "C2"
  exenergy(iMat) =  78.0e-9_fPrec
  anuc(iMat)     =  12.01_fPrec
  zatom(iMat)    =   6.00_fPrec
  rho(iMat)      =   4.52_fPrec
  emr(iMat)      =   0.25_fPrec
  hcut(iMat)     =   0.02_fPrec
  radl(iMat)     =   0.094_fPrec
  bnref(iMat)    =  70.0_fPrec

  ! Silicon
  iMat = iMat + 1
  colmats(iMat)  = "Si"
  exenergy(iMat) = 173.0e-9_fPrec
  anuc(iMat)     =  28.08_fPrec
  zatom(iMat)    =  14.00_fPrec
  rho(iMat)      =   2.33_fPrec
  emr(iMat)      =   0.441_fPrec
  hcut(iMat)     =   0.02_fPrec
  radl(iMat)     = one
  bnref(iMat)    = 120.14_fPrec

  ! Germanium
  iMat = iMat + 1
  colmats(iMat)  = "Ge"
  exenergy(iMat) = 350.0e-9_fPrec
  anuc(iMat)     =  72.63_fPrec
  zatom(iMat)    =  32.00_fPrec
  rho(iMat)      =   5.323_fPrec
  emr(iMat)      =   0.605_fPrec
  hcut(iMat)     =   0.02_fPrec
  radl(iMat)     = one
  bnref(iMat)    = 226.35_fPrec

  iMat = iMat + 1
  colmats(iMat)  = "MoGR"
  exenergy(iMat) =  87.1e-9_fPrec
  anuc(iMat)     =  13.53_fPrec
  zatom(iMat)    =   6.65_fPrec
  rho(iMat)      =   2.500_fPrec
  emr(iMat)      =   0.25_fPrec
  hcut(iMat)     =   0.02_fPrec
  radl(iMat)     =   0.1193_fPrec
  bnref(iMat)    =  76.7_fPrec

  iMat = iMat + 1
  colmats(iMat)  = "CuCD"
  exenergy(iMat) = 152.9e-9_fPrec
  anuc(iMat)     =  25.24_fPrec
  zatom(iMat)    =  11.90_fPrec
  rho(iMat)      =   5.40_fPrec
  emr(iMat)      =   0.308_fPrec
  hcut(iMat)     =   0.02_fPrec
  radl(iMat)     =   0.0316_fPrec
  bnref(iMat)    = 115.0_fPrec

  iMat = iMat + 1
  colmats(iMat)  = "Mo"
  exenergy(iMat) = 424.0e-9_fPrec
  anuc(iMat)     =  95.96_fPrec
  zatom(iMat)    =  42.00_fPrec
  rho(iMat)      =  10.22_fPrec
  emr(iMat)      =   0.481_fPrec
  hcut(iMat)     =   0.02_fPrec
  radl(iMat)     =   0.0096_fPrec
  bnref(iMat)    = 273.9_fPrec

  iMat = iMat + 1
  colmats(iMat)  = "Glid"
  exenergy(iMat) = 320.8e-9_fPrec
  anuc(iMat)     =  63.15_fPrec
  zatom(iMat)    =  28.80_fPrec
  rho(iMat)      =   8.93_fPrec
  emr(iMat)      =   0.418_fPrec
  hcut(iMat)     =   0.02_fPrec
  radl(iMat)     =   0.0144_fPrec
  bnref(iMat)    = 208.7_fPrec

  iMat = iMat + 1
  colmats(iMat)  = "Iner"
  exenergy(iMat) = 682.2e-9_fPrec
  anuc(iMat)     = 166.70_fPrec
  zatom(iMat)    =  67.70_fPrec
  rho(iMat)      =  18.00_fPrec
  emr(iMat)      =   0.578_fPrec
  hcut(iMat)     =   0.02_fPrec
  radl(iMat)     =   0.00385_fPrec
  bnref(iMat)    = 392.1_fPrec

  if(iMat > nrmat) then
    write(lerr,"(a)") "COLL> ERROR Variable imat > nrmat in collmat_init. Please increase nrmat."
    call prror
  end if

  ! The following two must always be the last two materials

  ! Vacuum
  iMat = iMat + 1
  colmats(iMat)  = "VA"
  exenergy(iMat) = zero
  anuc(iMat)     = zero
  zatom(iMat)    = zero
  rho(iMat)      = zero
  emr(iMat)      = zero
  hcut(iMat)     = zero
  radl(iMat)     = 1.0e12_fPrec
  bnref(iMat)    = zero

  ! Black Absorber
  iMat = iMat + 1
  colmats(iMat)  = "BL"
  exenergy(iMat) = c1e10
  anuc(iMat)     = zero
  zatom(iMat)    = zero
  rho(iMat)      = zero
  emr(iMat)      = zero
  hcut(iMat)     = zero
  radl(iMat)     = 1.0e12_fPrec
  bnref(iMat)    = zero

  if(iMat > nmat) then
    write(lerr,"(a)") "COLL> ERROR Variable imat > nmat in collmat_init. Please increase nmat."
    call prror
  end if

end subroutine collmat_init

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-16
!  Updated: 2019-09-16
!  Get collimator material number from name (case sensitive)
! ================================================================================================ !
integer function collmat_getCollMatID(matName)

  character(len=*), intent(in) :: matName
  integer i, matID

  matID = -1
  do i=1,nmat
    if(colmats(i) == matName) then
      matID = i
      exit
    end if
  end do
  collmat_getCollMatID = matID

end function collmat_getCollMatID

end module coll_materials
