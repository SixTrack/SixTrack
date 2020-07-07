!
!  This file contains all the main modules holding the SixTrack shared variables
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ================================================================================================ !
!  MAIN ARRAYS SIZE COMMON VARIABLES
!  Last modified: 2018-06-13
! ================================================================================================ !
module parpro

  implicit none

  integer, parameter :: mbea  = 99        ! Maximum number of beam-beam slices
  integer, parameter :: mcor  = 10        ! Maximum number of extra parameters (DIFF block)
  integer, parameter :: mcop  = mcor + 6  ! FOX/DA variable
  integer, parameter :: mpa   = 6         ! Maximum number of trajectories
  integer, parameter :: mran  = 500       ! Maximum number of inserted elements (subroutine ord)
  integer, parameter :: ncom  = 100       ! Maximum number of combinations of elements (COMB block)
  integer, parameter :: ncor1 = 600       ! Maximum number of corrections of closed orbit
  integer, parameter :: nema  = 15        ! Maximum order of the one turn map (DIFF block)
  integer, parameter :: ninv  = 1000      ! Number of invariances (postprocessing)
  integer, parameter :: nlya  = 10000     ! Something something postprocessing
  integer, parameter :: nmon1 = 600       ! Maximum number of monitors (closed orbit)
  integer, parameter :: nper  = 16        ! Maximum number of super periods (BLOC list, line 1)
  integer, parameter :: nplo  = 20000     ! Plotting
  integer, parameter :: npos  = 20000     ! Something something postprocessing
  integer, parameter :: nran  = 2000000   ! Maximum size for scaling nzfz
  integer, parameter :: nrco  = 5         ! Maximum order of compensation (RESO block)
  integer, parameter :: mmul  = 20        ! Maximum order of multipoles
  integer, parameter :: nelb  = 280       ! Maximum elements per BLOC

  ! Maximum length of strings
  integer, parameter :: mNameLen  = 48    ! Maximum length of element names. Keep in sync with MadX (and inside Geant4 code!)
  integer, parameter :: mFileName = 64    ! Maximum length of file names
  integer, parameter :: mPathName = 255   ! Maximum length of path names
  integer, parameter :: mStrLen   = 161   ! Standard string length
  integer, parameter :: mDivLen   = 132   ! Length of lout output lines
  integer, parameter :: mInputLn  = 1600  ! Buffer size for single lines read from input files

  integer :: ntr   =  1   ! Number of phase trombones

  integer :: nzfz  = -1   ! Number of allocated multipole random numbers
  integer :: nele  = -1   ! Number of allocated SINGle elements
  integer :: nblo  = -1   ! Number of allocated BLOCs
  integer :: nblz  = -1   ! Number of allocated STRUcture elements
  integer :: npart = -1   ! Number of allocated particles
  integer :: nbb   = -1   ! Beam-beam lenses

  integer, parameter :: nele_initial  = 500
  integer, parameter :: nblo_initial  = 100
  integer, parameter :: nblz_initial  = 1000
  integer, parameter :: npart_initial = 2
  integer, parameter :: nbb_initial   = 100

  ! Dividing line for output
  character(len=mDivLen), parameter :: str_divLine = repeat("-",mDivLen)

end module parpro

! ================================================================================================ !
!  PARBEAM COMMON VARIABLES
!  Last modified: 2018-06-12
! ================================================================================================ !
module parbeam

  use floatPrecision

  implicit none

  real(kind=fPrec), parameter :: xcut = 7.77_fPrec
  real(kind=fPrec), parameter :: ycut = 7.46_fPrec
  real(kind=fPrec), parameter :: h    =  1.0_fPrec/63.0_fPrec
  integer,          parameter :: nx   = 490
  integer,          parameter :: ny   = 470
  integer,          parameter :: idim = (nx+2)*(ny+2)

  integer,          save :: kstep
  real(kind=fPrec), save :: hrecip
  real(kind=fPrec), save :: wtreal(idim)
  real(kind=fPrec), save :: wtimag(idim)

  ! common /beam_exp/
  integer,          save :: beam_expflag      = 0       ! 0: Old BEAM block, 1: New BEAM::EXPERT
  logical,          save :: beam_expfile_open = .false. ! Have we opened the file 'beam_expert.txt'?

end module parbeam

! ================================================================================================ !
!  Global Settings Module
!  Last modiffied: 2018-06-10
!  Holds global settings values and parameters not directly related to the physics
! ================================================================================================ !
module mod_settings

  implicit none

  ! SETTINGS Block (fort.3)
  logical, save :: st_print        = .false. ! PRINT flag (fort.3)
  integer, save :: st_quiet        = 0       ! QUIET Level 0=verbose, 1=minimal, 2=quiet
  logical, save :: st_debug        = .false. ! Global DEBUG flag
  logical, save :: st_notrack      = .false. ! Flag to disable tracking (exit after initialisation)
  logical, save :: st_partsum      = .false. ! Flag to print final particle summary
  logical, save :: st_writefort12  = .false. ! Flag to write fort.12 after tracking
  logical, save :: st_fStateWrite  = .false. ! Dump particle final state file
  logical, save :: st_fStateText   = .false. ! Dump particle final state file as text
  logical, save :: st_fStateIons   = .false. ! Dump particle final state file with ion columns
  logical, save :: st_iStateWrite  = .false. ! Dump particle initial state file
  logical, save :: st_iStateText   = .false. ! Dump particle initial state file as text
  logical, save :: st_iStateIons   = .false. ! Dump particle initial state file with ion columns

  ! Checpoint/Restart Kills Switch Settings
  logical,              save :: st_killswitch = .false. ! Enables the kill on turn number debug feature
  integer, allocatable, save :: st_killturns(:)         ! List of killswitch turns

end module mod_settings

! ================================================================================================ !
!  THAT BIG COMMON VARIABLES MODULE
!  Last modified: 2019-01-14
! ================================================================================================ !
module mod_common

  use parpro
  use floatPrecision
  use numerical_constants
  use physical_constants, only : pmap
  use, intrinsic :: iso_fortran_env, only : int16, int32

  implicit none

  ! Parameters
  real(kind=fPrec),  parameter :: eps_dcum   = c1m6    ! Tolerance for machine length mismatch [m]

  ! Various Flags and Variables
  character(len=80), save      :: toptit(5)  = " "     ! DANGER: If the len changes, CRCHECK will break
  character(len=80), save      :: sixtit     = " "     ! DANGER: If the len changes, CRCHECK will break
  character(len=80), save      :: commen     = " "     ! DANGER: If the len changes, CRCHECK will break
  logical,           save      :: print_dcum = .false. ! Print dcum. Set in the SETTINGS block
  integer,           save      :: ithick     = 0       ! Thick tracking flag

  ! File Names and Units
  character(len=mFileName), public, save :: fort2   = "fort.2"   ! Name of machine geometry file
  character(len=mFileName), public, save :: fort3   = "fort.3"   ! Name of main input file
  character(len=mFileName), public, save :: fort6   = "fort.6"   ! Name of main output file (stdout)
  character(len=mFileName), public, save :: fort10  = "fort.10"  ! Name of main postprocessing file (text)
  character(len=mFileName), public, save :: fort110 = "fort.110" ! Name of main postprocessing file (binary)
  character(len=mFileName), public, save :: fort208 = "fort.208" ! Collimation/FLUKA

  integer,                  public, save :: unit10  = -1         ! Unit of main postprocessing file (text)
  integer,                  public, save :: unit110 = -1         ! Unit of main postprocessing file (binary)
  integer,                  public, save :: unit208 = -1         ! Collimation/FLUKA

  !  GENERAL VARIABLES
  ! ===================

  ! Error Variables
  integer,          save :: ierro      = 0

  ! Loop Variables
  integer,          save :: iu         = 0    ! Number of entries in the accelerator
  integer,          save :: il         = 0    ! Number of single elements
  integer,          save :: kanf       = 0    ! Structure index where the GO keyword is issued

  ! TRACK Block Variables
  real(kind=fPrec), save :: amp0       = zero    ! End amplitude
  integer,          save :: numl       = 1       ! Number of turns in the forward direction
  integer,          save :: numlr      = 0       ! Number of turns in the backward direction
  integer,          save :: numx       = 0       ! Checkpoint turn (turn-1)
  integer,          save :: napx       = 0       ! Number of amplitude variations
  integer,          save :: ird        = 0       ! Ignored
  integer,          save :: niu(2)     = 0       ! Start and stop structure element for optics calculation
  integer,          save :: numlcp     = 1000    ! How often to write checkpointing files
  integer,          save :: idfor      = 0       ! Add closed orbit to initial coordinates
  integer,          save :: irew       = 0       ! Rewind fort.59-90
  integer,          save :: iclo6      = 0       ! 6D closed orbit flags
  integer,          save :: nde(2)     = 0       ! Number of turns at flat bottom / energy ramping
  integer,          save :: nwr(4)     = 1       ! Writings to singletrackfile.dat
  integer,          save :: ntwin      = 1       ! How to calculate the distance in phase space
  integer,          save :: napxo      = 0       ! Original value of napx
  integer,          save :: napxto     = 0       ! Particles times turns
  integer,          save :: nnuml      = 0       ! Turn number for POSTPR
  logical,          save :: curveff    = .false. ! Enable the curvature effect in a combined function magnet
  logical,          save :: firstrun   = .true.  ! Switched off after the first completed turn (not necessarily turn 1 if C/R)
  logical,          save :: iexact     = .false. ! Exact solution of the equation of motion
  logical,          save :: rdfort13   = .false. ! Wheteher to read distribution from fort.13 or not

  ! INITIAL COORDINATES Block Variables
  real(kind=fPrec), save :: rat        = zero
  integer,          save :: iver       = 0

  ! BLOC DEFINITONS Block
  integer,          save :: mper       = 0    ! Number of super periods
  integer,          save :: msym(nper) = 0
  integer,          save :: mblo       = 0
  integer,          save :: mbloz      = 0

  ! STRUCTURE Block
  logical,          save :: strumcol   = .false. ! Whether the structure block is multicolumn or not

  ! Synchrotron Motion and Differential Algebra Block
  real(kind=fPrec), save :: qs         = zero ! Synchrotron tune [N/turn]
  real(kind=fPrec), save :: pma        = pmap ! Rest mass of the particle in MeV/c^2
  real(kind=fPrec), save :: phas       = zero ! Synchrotron acceleration phase [rad]
  real(kind=fPrec), save :: hsy(3)     = zero ! Cavity [voltage, unused, RF frequency]
  real(kind=fPrec), save :: crad       = zero ! electron radius * electron mass / pma
  real(kind=fPrec), save :: dppoff     = zero ! Offset relative momentum deviation
  real(kind=fPrec), save :: tlen       = one  ! Length of the accelerator [m]
  real(kind=fPrec), save :: mtcda      = one  ! Somthing FOX
  real(kind=fPrec), save :: dpscor     = one  ! Scaling factor for relative momentum deviation
  real(kind=fPrec), save :: sigcor     = one  ! Path length difference
  integer,          save :: iicav      = 0    ! Used between runcav and runda
  integer,          save :: ition      = 0    ! Transition energy switch:
  integer,          save :: idp        = 0    ! Synchrotron motion
  integer,          save :: ixcav      = 0    ! Stores ix, presumably for cavity
  integer,          save :: icode      = 0
  integer,          save :: idam       = 0
  integer,          save :: its6d      = 0

  ! RF Cavities
  integer,          save :: icy        = 0    ! Accelerating cavity: Number of "CAV" locations in STRUCT
  integer,          save :: ncy        = 0    ! Accelerating cavity: Number of "CAV" locations times super periods mper
  integer,          save :: ncy2       = 0    ! Accelerating cavity: Number of cavities (kz = +/- 12) in SING

  ! Organisation of Random Numbers
  integer,          save :: iorg       = 0    ! Organisation of random numbers flag
  integer,          save :: izu0       = 0    ! Start value for the random number generator
  integer,          save :: mcut       = 0    ! Sigma cut for random numbers
  integer,          save :: mout2      = 0    ! Magnet Fluctuations: Write fort.2 > fort.4 (mod_fluc)

  ! Iteration Errors
  real(kind=fPrec), save :: dma        = c1m12 ! Precision of closed orbit displacements
  real(kind=fPrec), save :: dmap       = c1m15 ! Precision of derivative of closed orbit displacements
  real(kind=fPrec), save :: dkq        = c1m10 ! Variations of quadrupole strengths
  real(kind=fPrec), save :: dqq        = c1m10 ! Precision of tunes
  real(kind=fPrec), save :: dsm0       = c1m10 ! Variations of sextupole strengths
  real(kind=fPrec), save :: dech       = c1m10 ! Precision of chromaticity correction
  real(kind=fPrec), save :: de0        = c1m9  ! Variations of momentum spread for chromaticity calculation
  real(kind=fPrec), save :: ded        = c1m9  ! Variations of momentum spread for evaluation of dispersion
  real(kind=fPrec), save :: dsi        = c1m9  ! Precision of desired orbit r.m.s. value
  real(kind=fPrec), save :: aper(2)    = c1e3  ! Precision of aperture limits

  integer,          save :: itco       = 500   ! Number of iterations for closed orbit calculation
  integer,          save :: itcro      = 10    ! Number of iterations for chromaticity correction
  integer,          save :: itqv       = 10    ! Number of iterations for Q adjustment

  ! Closed Orbit
  real(kind=fPrec), save :: di0(2)     = zero
  real(kind=fPrec), save :: dip0(2)    = zero
  real(kind=fPrec), save :: ta(6,6)    = zero

  ! Correction of Closed Orbit
  real(kind=fPrec), save :: sigma0(2)  = zero ! Desired RMS values of the randomly distributed closed orbit
  integer,          save :: iclo       = 0    ! Closed orbit correction switch
  integer,          save :: ncorru     = 0    ! Number of correctors to be used
  integer,          save :: ncorrep    = 0    ! Number of corrections

  integer,          save :: nhmoni     = 0
  integer,          save :: nhcorr     = 0
  integer,          save :: nvmoni     = 0
  integer,          save :: nvcorr     = 0

  real(kind=fPrec), save :: clo6(3)    = zero ! 6D closed orbit correction position
  real(kind=fPrec), save :: clop6(3)   = zero ! 6D closed orbit correction momentum

  ! Tune Variation
  real(kind=fPrec), save :: qw0(3)     = zero ! Qx/Qy/delta_Q
  integer,          save :: iq(3)      = 0    ! Index of elements
  integer,          save :: iqmod      = 0    ! Switch to calculate the tunes
  integer,          save :: iqmod6     = 0    ! Switch to calculate the tunes

  ! Linear Optics
  real(kind=fPrec), save :: eui        = 0    ! Eigenemittance
  real(kind=fPrec), save :: euii       = 0    ! Eigenemittance
  integer,          save :: ilin       = 0    ! Linear optics 4D/6D calculation flag
  integer,          save :: nt         = 0    ! Number of the blocks to which the linear parameter will be printed
  integer,          save :: ntco       = 0    ! Switch to write out linear coupling parameters
  integer,          save :: iprint     = 0    ! Print the linear optics functions at ...
  integer,          save :: nlin       = 0    ! Number of elements for write out

  ! Combination of Elements
  integer,          save :: icomb0(20) = 0    ! Index of single element
  integer,          save :: icoe       = 0    ! Number of combinations (lines in block)

  ! Search for Optimum Places to Compensate Resonances and Sub-Resonance Calculation
  integer,          save :: mesa       = 0    ! Number of positions to be checked
  integer,          save :: mp         = 0    ! Order of the resonance
  integer,          save :: m21        = 0    ! Resonances of order 1
  integer,          save :: m22        = 0    ! Resonances of order 2
  integer,          save :: m23        = 0    ! Resonances of order 3
  integer,          save :: ise1       = 0    ! Distance to a resonance
  integer,          save :: ise2       = 0    ! Distance to a resonance
  integer,          save :: ise3       = 0    ! Distance to a resonance
  integer,          save :: ise        = 0    ! Flag on/off
  integer,          save :: isub       = 0    ! Sub-resonance calculation switch
  integer,          save :: nta        = 0    ! Lowest order of the resonance
  integer,          save :: nte        = 0    ! Highest order of the resonance
  integer,          save :: ipt        = 0    ! Switch to change the nearest distance to the resonance

  integer,          save :: nre        = 0    ! Number of resonances
  integer,          save :: npp        = 0    ! Order of the resonance
  integer,          save :: nrr(5)     = 0    ! Resonances of order n
  integer,          save :: ipr(5)     = 0    ! Distance to the resonance
  integer,          save :: nur        = 0    ! Number of sub-resonances
  integer,          save :: nu(5)      = 0    ! Order of the multipole
  real(kind=fPrec), save :: totl       = zero ! Length of the accelerator [m]
  real(kind=fPrec), save :: qxt        = zero ! Horizontal tune including the integer part
  real(kind=fPrec), save :: qzt        = zero ! Vertical tune including the integer part
  real(kind=fPrec), save :: tam1       = zero ! Horizontal amplitudes [mm]
  real(kind=fPrec), save :: tam2       = zero ! Vertical amplitudes [mm]
  integer,          save :: irmod2     = 0    ! Resonance compensation on/off
  integer,          save :: nch        = 0    ! Switch for the chromaticity correction
  integer,          save :: nqc        = 0    ! Switch for the tune adjustment

  real(kind=fPrec), save :: dtr(10)    = zero
  integer,          save :: ire(12)    = 0

  real(kind=fPrec), save :: rtc(9,18,10,5) = zero
  real(kind=fPrec), save :: rts(9,18,10,5) = zero

  ! Beam-Beam Variables
  real(kind=fPrec), save :: sigz       = zero ! RMS bunch length
  real(kind=fPrec), save :: sige       = zero ! RMS energy spread
  real(kind=fPrec), save :: partnum    = zero ! Number of particles in bunch
  real(kind=fPrec), save :: parbe14    = zero !
  real(kind=fPrec), save :: emitx      = zero ! Horisontal emittance
  real(kind=fPrec), save :: emity      = zero ! Vertical emittance
  real(kind=fPrec), save :: emitz      = zero ! Longitudinal emittance
  real(kind=fPrec), save :: gammar     = one  ! Inverse Lorentz factor
  real(kind=fPrec), save :: brho       = zero ! magnetic rigitidy of beam [Tm]
  integer,          save :: ibb6d      = 0    ! 6D beam-beam switch
  integer,          save :: nbeam      = 0    ! Beam-beam elements flag
  integer,          save :: ibbc       = 0    ! Switch for linear coupling considered in 4D and 6D
  integer,          save :: ibeco      = 1    ! Subtract the closed orbit introduced by BEAM/WIRE
  integer,          save :: ibtyp      = 0    ! Switch to use the fast beam–beam algorithms
  integer,          save :: lhc        = 1    ! Switch for the LHC with its anti-symmetric IR

  ! Decoupling of Motion
  integer,          save :: iskew      = 0    ! Skew switch
  integer,          save :: nskew(6)   = 0    ! Index to skew quadrupole families
  real(kind=fPrec), save :: qwsk(2)    = zero ! Horisontal and vertical tune
  real(kind=fPrec), save :: betx(2)    = zero
  real(kind=fPrec), save :: betz(2)    = zero
  real(kind=fPrec), save :: alfx(2)    = zero
  real(kind=fPrec), save :: alfz(2)    = zero

  ! RF Multipoles
  integer,          save :: iord       = 0
  integer,          save :: nordm      = 0

  ! Reference Particle
  real(kind=fPrec),    save :: e0      = zero ! Reference energy [MeV]
  real(kind=fPrec),    save :: e0f     = zero ! Reference momentum [MeV/c]
  real(kind=fPrec),    save :: nucm0   = pmap ! Reference mass [MeV/c^2]
  real(kind=fPrec),    save :: nucmda  = pmap ! Reference mass [MeV/c^2] (DA)
  real(kind=fPrec),    save :: gamma0  = one  ! Reference beam Lorentz factor
  real(kind=fPrec),    save :: beta0   = zero ! Reference beam relativistic beta
  integer(kind=int16), save :: aa0     = 1    ! Reference nucleon number
  integer(kind=int16), save :: zz0     = 1    ! Reference charge multiplicity
  integer(kind=int16), save :: qq0     = 1    ! Reference charge
  integer(kind=int32), save :: pdgid0  = 2212 ! Reference particle PDG ID
  integer,             save :: pdgyear = 2002 ! Reference particle PDG year

  ! Tracking Particles
  real(kind=fPrec), save :: ej(mpa)    = zero ! Particle energy
  real(kind=fPrec), save :: ejf(mpa)   = zero ! Particle momentum

  ! Other Variables
  integer,          save :: ichromc    = 0
  integer,          save :: ilinc      = 0    ! 2 = Beam-beam closed orbit calc
  integer,          save :: iqmodc     = 0
  real(kind=fPrec), save :: corr(3,3)  = zero
  real(kind=fPrec), save :: chromc(2)  = 9.999999e23_fPrec
  real(kind=fPrec), save :: wxys(3)    = zero
  real(kind=fPrec), save :: clon(6)    = zero

  real(kind=fPrec), save :: aml6(6,6)  = zero ! 6D / FOX
  real(kind=fPrec), save :: edcor(2)   = zero ! Chromaticity-related

  ! Post-Procesing
  real(kind=fPrec), save :: dphix      = zero ! Horisontal angle interval used to stroboscope phase space [rad]
  real(kind=fPrec), save :: dphiz      = zero ! Vertical angle interval used to stroboscope phase space [rad]
  real(kind=fPrec), save :: qx0        = zero ! Horisontal tune
  real(kind=fPrec), save :: qz0        = zero ! Vertical tune
  real(kind=fPrec), save :: dres       = one  ! Maximum distance to the resonance
  real(kind=fPrec), save :: dfft       = one  ! Find resonances with the FFT spectrum below
  real(kind=fPrec), save :: cma1       = one  ! Scale for Lyapunov analysis
  real(kind=fPrec), save :: cma2       = one  ! Scale for Lyapunov analysis

  integer,          save :: nstart     = 0    ! Start turn number for the analysis
  integer,          save :: nstop      = 0    ! Stop turn number for the analysis
  integer,          save :: iskip      = 1    ! Number of samples to skip
  integer,          save :: iconv      = 0    ! Switch: tracking data are not normalised linearly
  integer,          save :: imad       = 0    ! Switch: MAD-X related

  integer,          save :: ipos       = 0    ! Post-processing on/off
  integer,          save :: iav        = 1    ! Averaging interval of the values of the distance in phase space
  integer,          save :: iwg        = 1    ! Switch: weighting of the slope calculation of the distance in phase space
  integer,          save :: ivox       = 1    ! Switch: tune close to an integer
  integer,          save :: ivoz       = 1    ! Switch: tune close to an integer
  integer,          save :: ires       = 1    ! Closest resonances search up to order
  integer,          save :: ifh        = 0    ! FFT analysis tune interval

  integer,          save :: kwtype     = 0    ! Disabled
  integer,          save :: itf        = 0    ! Switch to get PS file of plots
  integer,          save :: icr        = 0    ! Disabled
  integer,          save :: idis       = 0    ! Switch to select plot
  integer,          save :: icow       = 0    ! Switch to select plot
  integer,          save :: istw       = 0    ! Switch to select plot
  integer,          save :: iffw       = 0    ! Switch to select plot
  integer,          save :: nprint     = 1    ! Switch to stop the printing of the post-processing fort.6
  integer,          save :: ndafi      = 1    ! Number of data files to be processed

  !  ALLOCATABLES
  ! ==============

  ! Block Element Indexed (nblo)
  character(len=:), allocatable, save :: bezb(:)        ! Block Elements: Name
  real(kind=fPrec), allocatable, save :: elbe(:)        ! Block Elements: Length
  integer,          allocatable, save :: mtyp(:,:)      ! Block Elements: Indices of single elements (:,nelb)

  real(kind=fPrec), allocatable, save :: bl1(:,:,:)     ! Block Elements: (:,2,6)
  real(kind=fPrec), allocatable, save :: bl2(:,:,:)     ! Block Elements: (:,2,6)
  integer,          allocatable, save :: mel(:)         ! Block Elements: Number of single elements

  ! Single Element Indexed (nele)
  character(len=:), allocatable, save :: bez(:)         ! Single Elements: 1st field, name
  character(len=:), allocatable, save :: bezl(:)        ! Single Elements: Linear optics write out
  real(kind=fPrec), allocatable, save :: ed(:)          ! Single Elements: 3rd field, additional
  real(kind=fPrec), allocatable, save :: ek(:)          ! Single Elements: 4th field, additional
  real(kind=fPrec), allocatable, save :: el(:)          ! Single Elements: 5th field, length [m]
  real(kind=fPrec), allocatable, save :: bbbx(:)        ! Single Elements: 6th field
  real(kind=fPrec), allocatable, save :: bbby(:)        ! Single Elements: 7th field
  real(kind=fPrec), allocatable, save :: bbbs(:)        ! Single Elements: 8th field
  real(kind=fPrec), allocatable, save :: sm(:)          ! Single Elements: Non-linear field, avg strength
  integer,          allocatable, save :: kz(:)          ! Single Elements: 2nd field, type
  integer,          allocatable, save :: kp(:)          ! Single Elements: Additional flag
  real(kind=fPrec), allocatable, save :: dki(:,:)       ! Single Elements: [H bend kick, V bend kick, dipole length] (:,3)

  real(kind=fPrec), allocatable, save :: a(:,:,:)       ! Something Something Matrix (:,2,6)

  real(kind=fPrec), allocatable, save :: hsyc(:)        ! Accelerating cavity: 'Frequency'
  real(kind=fPrec), allocatable, save :: phasc(:)       ! Accelerating cavity: Lag phase

  real(kind=fPrec), allocatable, save :: benkc(:)       ! Multipoles: Bending strength of the dipole [mrad]
  real(kind=fPrec), allocatable, save :: r00(:)         ! Multipoles: Reference radius [mm]
  real(kind=fPrec), allocatable, save :: scalemu(:)     ! Multipoles: Scale (DYNK)
  integer,          allocatable, save :: nmu(:)         ! Multipoles: Max multipole order
  integer,          allocatable, save :: irm(:)         ! Multipoles: Index of the associated element

  real(kind=fPrec), allocatable, save :: freq_rfm(:)    ! RF Multipoles:
  integer,          allocatable, save :: nmu_rf(:)      ! RF Multipoles: Max multipole order
  integer,          allocatable, save :: irm_rf(:)      ! RF Multipoles: Index of the associated element

  real(kind=fPrec), allocatable, save :: xpl(:)         ! Displacement (DISP): Value of horizontal displacement [mm]
  real(kind=fPrec), allocatable, save :: xrms(:)        ! Displacement (DISP): RMS of horizontal displacement [mm]
  real(kind=fPrec), allocatable, save :: zpl(:)         ! Displacement (DISP): Value of vertical displacement [mm]
  real(kind=fPrec), allocatable, save :: zrms(:)        ! Displacement (DISP): RMS of vertical displacement [mm]

  character(len=:), allocatable, save :: bezr(:,:)      ! Organisation of Random Numbers (3,:)

  integer,          allocatable, save :: kpa(:)         ! Tune Variations: Element markers

  integer,          allocatable, save :: ncororb(:)     ! Obrit Correctors: Flag

  real(kind=fPrec), allocatable, save :: crabph(:)      ! Crab Cavities: Phase of the excitation
  real(kind=fPrec), allocatable, save :: crabph2(:)     ! Crab Cavities: Order 2
  real(kind=fPrec), allocatable, save :: crabph3(:)     ! Crab Cavities: Order 3
  real(kind=fPrec), allocatable, save :: crabph4(:)     ! Crab Cavities: Order 4

  real(kind=fPrec), allocatable, save :: ratioe(:)      ! Combination of Elements: Ratio of the magnetic strength
  integer,          allocatable, save :: iratioe(:)     ! Combination of Elements: Index

  integer,          allocatable, save :: isea(:)        ! Compensate Resonance: Element index

  real(kind=fPrec), allocatable, save :: parbe(:,:)     ! Beam-Beam: Input values (:,18)
  real(kind=fPrec), allocatable, save :: ptnfac(:)      ! Beam-Beam: Strength ratio

  real(kind=fPrec), allocatable, save :: acdipph(:)     ! AC Dipole: Phase (el)
  integer,          allocatable, save :: nturn1(:)      ! AC Dipole: Turns to free of excitation at the beginning of the run
  integer,          allocatable, save :: nturn2(:)      ! AC Dipole: Turns to ramp up the excitation amplitude
  integer,          allocatable, save :: nturn3(:)      ! AC Dipole: Turns of constant excitation amplitude
  integer,          allocatable, save :: nturn4(:)      ! AC Dipole: Turns to ramp down the excitation amplitude

  integer,          allocatable, save :: imtr(:)        ! Phase Trombones

  ! Single Element and Multipole Indexed (nele,mmul)
  real(kind=fPrec), allocatable, save :: bk0(:,:)       ! Multipoles: B-value
  real(kind=fPrec), allocatable, save :: ak0(:,:)       ! Multipoles: A-value
  real(kind=fPrec), allocatable, save :: bka(:,:)       ! Multipoles: B-RMS
  real(kind=fPrec), allocatable, save :: aka(:,:)       ! Multipoles: A-RMS
  real(kind=fPrec), allocatable, save :: norrfamp(:,:)  ! RF Multipoles:
  real(kind=fPrec), allocatable, save :: norrfph(:,:)   ! RF Multipoles:
  real(kind=fPrec), allocatable, save :: skrfamp(:,:)   ! RF Multipoles:
  real(kind=fPrec), allocatable, save :: skrfph(:,:)    ! RF Multipoles:

  ! Structure Element Indexed (nblz)
  real(kind=fPrec), allocatable, save :: sigmoff(:)
  real(kind=fPrec), allocatable, save :: tiltc(:)       ! Magnet Tilt: Cos component
  real(kind=fPrec), allocatable, save :: tilts(:)       ! Magnet Tilt: Sin component
  integer,          allocatable, save :: icext(:)       ! Magnet Error: Index (mod_fluc)
  integer,          allocatable, save :: mzu(:)         ! Magnet Error: Random number index
  integer,          allocatable, save :: icextal(:)     ! Magnet Misalignment: index (mod_fluc)

  character(len=:), allocatable, save :: bezs(:)        ! Name of structure element
  real(kind=fPrec), allocatable, save :: elpos(:)       ! Position of structure element
  integer,          allocatable, save :: ic(:)          ! Structure to single/block element map0

  real(kind=fPrec), allocatable, save :: xsi(:)         ! Tracking: Horisontal displacement
  real(kind=fPrec), allocatable, save :: zsi(:)         ! Tracking: Vertical displacement
  real(kind=fPrec), allocatable, save :: smi(:)         ! Tracking: Magnetic kick
  real(kind=fPrec), allocatable, save :: smizf(:)       ! Random Numbers: zfz*ek

  real(kind=fPrec), allocatable, save :: dcum(:)        ! Machine length in m (0:nblz+1)

  integer,          allocatable, save :: imbb(:)        ! Beam-Beam:

  ! Structure Element and Multipole Indexed (nblz,mmul)
  real(kind=fPrec), allocatable, save :: aaiv(:,:)      ! Multipoles:
  real(kind=fPrec), allocatable, save :: bbiv(:,:)      ! Multipoles:
  real(kind=fPrec), allocatable, save :: amultip(:,:)   ! Multipoles:
  real(kind=fPrec), allocatable, save :: bmultip(:,:)   ! Multipoles:

  ! Number of Particles (npart)
  real(kind=fPrec), allocatable, save :: track6d(:,:)   ! Beam-Beam (6,:)

  ! Random Numbers Indexed (nzfz)
  real(kind=fPrec), allocatable, save :: zfz(:)         ! Magnet errors

  ! Number of Multipoles (mmul)
  real(kind=fPrec), allocatable, save :: field_cos(:,:) ! RF Multipoiles (2,:)
  real(kind=fPrec), allocatable, save :: fsddida(:,:)   ! RF Multipoiles (2,:)
  real(kind=fPrec), allocatable, save :: field_sin(:,:) ! RF Multipoiles (2,:)
  real(kind=fPrec), allocatable, save :: fcodda(:,:)    ! RF Multipoiles (2,:)

  ! Number of Monitors (nmon1)
  real(kind=fPrec), allocatable, save :: betam(:,:)     ! Orbit Correction: (:,2)
  real(kind=fPrec), allocatable, save :: pam(:,:)       ! Orbit Correction: (:,2)
  real(kind=fPrec), allocatable, save :: bclorb(:,:)    ! Orbit Correction: (:,2)

  ! Number of Orbit Corrections (ncor1)
  real(kind=fPrec), allocatable, save :: betac(:,:)     ! Orbit Correction: (:,2)
  real(kind=fPrec), allocatable, save :: pac(:,:)       ! Orbit Correction: (:,2)

  ! Number of Combinations of Elements (ncom)
  real(kind=fPrec), allocatable, save :: ratio(:,:)     ! Combination of Elements: Ratio (:,20)
  integer,          allocatable, save :: icomb(:,:)     ! Combination of Elements: Index (:,20)

  ! Number of Beam-Beam Lenses (nbb)
  integer,          allocatable, save :: nbeaux(:)      ! (:)
  real(kind=fPrec), allocatable, save :: sigman(:,:)    ! (2,:)
  real(kind=fPrec), allocatable, save :: sigman2(:,:)   ! sigman^2 (2,:)
  real(kind=fPrec), allocatable, save :: sigmanq(:,:)   ! (2,:)
  real(kind=fPrec), allocatable, save :: clobeam(:,:)   ! (6,:)
  real(kind=fPrec), allocatable, save :: beamoff(:,:)   ! (6,:)
  real(kind=fPrec), allocatable, save :: bbcu(:,:)      ! (:,12)

  ! Number of Phase Trombones (ntr)
  real(kind=fPrec), allocatable, save :: cotr(:,:)      ! (:,6)
  real(kind=fPrec), allocatable, save :: rrtr(:,:,:)    ! (:,6,6)

contains

subroutine mod_common_expand_arrays(nele_new, nblo_new, nblz_new, npart_new, nbb_new)

  use mod_alloc
  use mod_settings
  use numerical_constants, only : zero,one

  implicit none

  integer, intent(in) :: nele_new
  integer, intent(in) :: nblo_new
  integer, intent(in) :: nblz_new
  integer, intent(in) :: npart_new
  integer, intent(in) :: nbb_new

  logical :: firstRun   = .true.
  integer :: nele_prev  = -2
  integer :: nblo_prev  = -2
  integer :: nblz_prev  = -2
  integer :: npart_prev = -2
  integer :: nbb_prev   = -2

  if(nele_new /= nele_prev) then
    call alloc(ed,                   nele_new,       zero,   "ed")
    call alloc(el,                   nele_new,       zero,   "el")
    call alloc(ek,                   nele_new,       zero,   "ek")
    call alloc(sm,                   nele_new,       zero,   "sm")
    call alloc(kz,                   nele_new,       0,      "kz")
    call alloc(kp,                   nele_new,       0,      "kp")
    call alloc(bbbx,                 nele_new,       zero,   "bbbx")
    call alloc(bbby,                 nele_new,       zero,   "bbby")
    call alloc(bbbs,                 nele_new,       zero,   "bbbs")
    call alloc(xpl,                  nele_new,       zero,   "xpl")
    call alloc(zpl,                  nele_new,       zero,   "zpl")
    call alloc(xrms,                 nele_new,       zero,   "xrms")
    call alloc(zrms,                 nele_new,       zero,   "zrms")
    call alloc(a,                    nele_new,2,6,   zero,   "a")
    call alloc(hsyc,                 nele_new,       zero,   "hsyc")
    call alloc(phasc,                nele_new,       zero,   "phasc")
    call alloc(bk0,                  nele_new, mmul, zero,   "bk0")
    call alloc(ak0,                  nele_new, mmul, zero,   "ak0")
    call alloc(bka,                  nele_new, mmul, zero,   "bka")
    call alloc(aka,                  nele_new, mmul, zero,   "aka")
    call alloc(benkc,                nele_new,       zero,   "benkc")
    call alloc(norrfamp,             nele_new, mmul, zero,   "norrfamp")
    call alloc(norrfph,              nele_new, mmul, zero,   "norrfph")
    call alloc(skrfamp,              nele_new, mmul, zero,   "skrfamp")
    call alloc(skrfph,               nele_new, mmul, zero,   "skrfph")
    call alloc(freq_rfm,             nele_new,       zero,   "freq_rfm")
    call alloc(r00,                  nele_new,       zero,   "r00")
    call alloc(scalemu,              nele_new,       one,    "scalemu")
    call alloc(irm,                  nele_new,       0,      "irm")
    call alloc(irm_rf,               nele_new,       0,      "irm_rf")
    call alloc(nmu,                  nele_new,       0,      "nmu")
    call alloc(nmu_rf,               nele_new,       0,      "nmu_rf")
    call alloc(bezr,    mNameLen, 3, nele_new,       " ",    "bezr")
    call alloc(kpa,                  nele_new,       0,      "kpa")
    call alloc(bez,     mNameLen,    nele_new,       " ",    "bez")
    call alloc(bezl,    mNameLen,    nele_new,       " ",    "bezl")
    call alloc(ncororb,              nele_new,       0,      "ncororb")
    call alloc(ratioe,               nele_new,       one,    "ratioe")
    call alloc(iratioe,              nele_new,       0,      "iratioe")
    call alloc(isea,                 nele_new,       0,      "isea")
    call alloc(dki,                  nele_new, 3,    zero,   "dki")
    call alloc(parbe,                nele_new, 18,   zero,   "parbe")
    call alloc(ptnfac,               nele_new,       zero,   "ptnfac")
    call alloc(imtr,                 nele_new,       0,      "imtr")
    call alloc(acdipph,              nele_new,       zero,   "acdipph")
    call alloc(nturn1,               nele_new,       0,      "nturn1")
    call alloc(nturn2,               nele_new,       0,      "nturn2")
    call alloc(nturn3,               nele_new,       0,      "nturn3")
    call alloc(nturn4,               nele_new,       0,      "nturn4")
    call alloc(crabph,               nele_new,       zero,   "crabph")
    call alloc(crabph2,              nele_new,       zero,   "crabph2")
    call alloc(crabph3,              nele_new,       zero,   "crabph3")
    call alloc(crabph4,              nele_new,       zero,   "crabph4")
  end if

  if(nblo_new /= nblo_prev) then
    call alloc(bezb,    mNameLen,    nblo_new,       " ",    "bezb")
    call alloc(elbe,                 nblo_new,       zero,   "elbe")
    call alloc(mel,                  nblo_new,       0,      "mel")
    call alloc(mtyp,                 nblo_new, nelb, 0,      "mtyp")
    call alloc(bl1,                  nblo_new, 2, 6, zero,   "bl1")
    call alloc(bl2,                  nblo_new, 2, 6, zero,   "bl2")
  end if

  if(nblz_new /= nblz_prev) then
    call alloc(ic,                   nblz_new,       0,      "ic")
    call alloc(elpos,                nblz_new+1,     zero,   "elpos", 0)
    call alloc(bezs,      mNameLen,  nblz_new,       " ",    "bezs")
    call alloc(mzu,                  nblz_new,       0,      "mzu")
    call alloc(imbb,                 nblz_new,       0,      "imbb")
    call alloc(icext,                nblz_new,       0,      "icext")
    call alloc(icextal,              nblz_new,       0,      "icextal")
    call alloc(tiltc,                nblz_new,       one,    "tiltc")
    call alloc(tilts,                nblz_new,       zero,   "tilts")
    call alloc(xsi,                  nblz_new,       zero,   "xsi")
    call alloc(zsi,                  nblz_new,       zero,   "zsi")
    call alloc(smi,                  nblz_new,       zero,   "smi")
    call alloc(smizf,                nblz_new,       zero,   "smizf")
    call alloc(aaiv,          mmul,  nblz_new,       zero,   "aaiv")
    call alloc(bbiv,          mmul,  nblz_new,       zero,   "bbiv")
    call alloc(amultip,       mmul,  nblz_new,       zero,   "amultip")
    call alloc(bmultip,       mmul,  nblz_new,       zero,   "bmultip")
    call alloc(dcum,                 nblz_new+1,     zero,   "dcum", 0)
    call alloc(sigmoff,              nblz_new,       zero,   "sigmoff")
  end if

  if(npart_new /= npart_prev) then
    call alloc(track6d, 6,           npart_new,      zero,   "track6d")
  end if

  if(nbb_new /= nbb_prev) then
    call alloc(nbeaux,               nbb,            0,      "nbeaux")
    call alloc(sigman,            2, nbb,            zero,   "sigman")
    call alloc(sigman2,           2, nbb,            zero,   "sigman2")
    call alloc(sigmanq,           2, nbb,            zero,   "sigmanq")
    call alloc(clobeam,           6, nbb,            zero,   "clobeam")
    call alloc(beamoff,           6, nbb,            zero,   "beamoff")
    call alloc(bbcu,                 nbb, 12,        zero,   "bbcu")
  end if

  ! The arrays that don't currently have scalable sizes only need to be allocated once
  if(firstRun) then
    call alloc(betam,                nmon1, 2,       zero,   "betam")
    call alloc(pam,                  nmon1, 2,       zero,   "pam")
    call alloc(bclorb,               nmon1, 2,       zero,   "bclorb")

    call alloc(betac,                ncor1, 2,       zero,   "betac")
    call alloc(pac,                  ncor1, 2,       zero,   "pac")

    call alloc(ratio,                ncom,  20,      zero,   "ratio")
    call alloc(icomb,                ncom,  20,      0,      "icomb")

    call alloc(field_cos,         2, mmul,           zero,   "field_cos")
    call alloc(fsddida,           2, mmul,           zero,   "fsddida")
    call alloc(field_sin,         2, mmul,           zero,   "field_sin")
    call alloc(fcodda,            2, mmul,           zero,   "fcodda")

    call alloc(cotr,                 ntr,  6,        zero,   "cotr")
    call alloc(rrtr,                 ntr,  6,6,      zero,   "rrtr")
  end if

  firstRun   = .false.
  nele_prev  = nele_new
  nblo_prev  = nblo_new
  nblz_prev  = nblz_new
  npart_prev = npart_new

end subroutine mod_common_expand_arrays

end module mod_common

! ================================================================================================ !
!  DA COMMON VARIABLES
!  Last modified: 2019-01-15
! ================================================================================================ !
module mod_common_da

  use parpro
  use floatPrecision
  use numerical_constants

  implicit none

  ! Differential Algebra (DIFF)
  real(kind=fPrec), save :: preda      = c1m38 ! Precision needed by the DA package
  integer,          save :: idial      = 0     ! DIFF block switch
  integer,          save :: nord       = 0     ! Order of the map
  integer,          save :: nvar       = 0     ! Number of the variables
  integer,          save :: nvar2      = 0     ! Number of the variables (not with ncor added)
  integer,          save :: nsix       = 0     ! Switch to calculate a 5x6 instead of a 6x6 map
  integer,          save :: ncor       = 0     ! Number of zero-length elements to be additional parameters
  integer,          save :: ipar(mcor) = 0     ! DIFF variable

  ! Normal Forms (NORM)
  integer,          save :: inorm      = 0     ! NORM block switch
  integer,          save :: nordf      = 0     ! Order of the Normal Form
  integer,          save :: nvarf      = 0     ! Number of variables
  integer,          save :: nord1      = 1     ! 3rd variable in NORM (not in manual)
  integer,          save :: idptr      = 0     ! 4th variable in NORM (not in manual)
  integer,          save :: imod1      = 0     ! Mode
  integer,          save :: imod2      = 0     ! Mode
  integer,          save :: ndimf      = 0     ! NORM vriable

end module mod_common_da

! ================================================================================================ !
!  TRACKING COMMON VARIABLES
!  Last modified: 2018-06-21
! ================================================================================================ !
module mod_common_track

  use parpro
  use floatPrecision
  use numerical_constants

  implicit none

  ! Tracking
  real(kind=fPrec), save :: x(mpa,2) = zero
  real(kind=fPrec), save :: y(mpa,2) = zero
  real(kind=fPrec), save :: amp(2)   = [c1m3,zero]
  real(kind=fPrec), save :: bet0(2)  = zero
  real(kind=fPrec), save :: alf0(2)  = zero
  real(kind=fPrec), save :: clo(2)   = zero
  real(kind=fPrec), save :: clop(2)  = zero
  integer,          save :: nwri     = 0     ! Flag for frequency of calls to writebin. Set by nwr(3) in TRAC

  ! Chromaticity
  real(kind=fPrec), save :: cro(2)   = zero  ! Desired values of the chromaticity
  integer,          save :: crois(2) = 0     ! Index of the elements in the single elements list
  integer,          save :: ichrom   = 0     ! Flag for calculation, 1: "traditional", 2: including beam-beam, 3: both

  ! tas
  real(kind=fPrec), save :: tasm(6,6)

  ! Allocatables
  integer,          allocatable, save :: ktrack(:)  ! (nblz)
  real(kind=fPrec), allocatable, save :: strack(:)  ! (nblz)
  real(kind=fPrec), allocatable, save :: strackc(:) ! (nblz)
  real(kind=fPrec), allocatable, save :: stracks(:) ! (nblz)
  real(kind=fPrec), allocatable, save :: strackx(:) ! (nblz)
  real(kind=fPrec), allocatable, save :: strackz(:) ! (nblz)
  real(kind=fPrec), allocatable, save :: dpsv1(:)   ! (npart)

  ! Linear Optics
  real(kind=fPrec), allocatable, save :: tbetax(:)  ! (nblz)
  real(kind=fPrec), allocatable, save :: tbetay(:)  ! (nblz)
  real(kind=fPrec), allocatable, save :: talphax(:) ! (nblz)
  real(kind=fPrec), allocatable, save :: talphay(:) ! (nblz)
  real(kind=fPrec), allocatable, save :: torbx(:)   ! (nblz)
  real(kind=fPrec), allocatable, save :: torby(:)   ! (nblz)
  real(kind=fPrec), allocatable, save :: torbxp(:)  ! (nblz)
  real(kind=fPrec), allocatable, save :: torbyp(:)  ! (nblz)
  real(kind=fPrec), allocatable, save :: tdispx(:)  ! (nblz)
  real(kind=fPrec), allocatable, save :: tdispy(:)  ! (nblz)
  real(kind=fPrec), allocatable, save :: tdispxp(:) ! (nblz)
  real(kind=fPrec), allocatable, save :: tdispyp(:) ! (nblz)

  ! Substitute variables for x,y and is for DA version
  real(kind=fPrec), save :: xxtr(mpa,2)
  real(kind=fPrec), save :: yytr(mpa,2)

contains

subroutine mod_commont_expand_arrays(nblz_new,npart_new)

  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  integer, intent(in) :: nblz_new
  integer, intent(in) :: npart_new

  call alloc(ktrack,  nblz_new,  0,    "ktrack")
  call alloc(strack,  nblz_new,  zero, "strack")
  call alloc(strackc, nblz_new,  zero, "strackc")
  call alloc(stracks, nblz_new,  zero, "stracks")
  call alloc(strackx, nblz_new,  zero, "strackx")
  call alloc(strackz, nblz_new,  zero, "strackz")

  call alloc(tbetax,  nblz_new,  zero, "tbetax")
  call alloc(tbetay,  nblz_new,  zero, "tbetay")
  call alloc(talphax, nblz_new,  zero, "talphax")
  call alloc(talphay, nblz_new,  zero, "talphay")
  call alloc(torbx,   nblz_new,  zero, "torbx")
  call alloc(torby,   nblz_new,  zero, "torby")
  call alloc(torbxp,  nblz_new,  zero, "torbxp")
  call alloc(torbyp,  nblz_new,  zero, "torbyp")
  call alloc(tdispx,  nblz_new,  zero, "tdispx")
  call alloc(tdispy,  nblz_new,  zero, "tdispy")
  call alloc(tdispxp, nblz_new,  zero, "tdispxp")
  call alloc(tdispyp, nblz_new,  zero, "tdispyp")

  call alloc(dpsv1,   npart_new, zero, "dpsv1")

end subroutine mod_commont_expand_arrays

! Copy from actual variables to temp DA variables
subroutine comt_daStart
  xxtr(1:mpa,1:2) = x(1:mpa,1:2)
  yytr(1:mpa,1:2) = y(1:mpa,1:2)
end subroutine comt_daStart

! Copy from temp DA variables to actual variables
subroutine comt_daEnd
  x(1:mpa,1:2) = xxtr(1:mpa,1:2)
  y(1:mpa,1:2) = yytr(1:mpa,1:2)
end subroutine comt_daEnd

end module mod_common_track

! ================================================================================================ !
!  Main Variables Used for Particle Tracking
!  Last modified: 2019-05-09
! ================================================================================================ !
module mod_common_main

  use parpro
  use floatPrecision
  use numerical_constants
  use, intrinsic :: iso_fortran_env, only : int16, int32

  implicit none

  real(kind=fPrec), save :: qw(2)      = zero
  real(kind=fPrec), save :: qwc(3)     = zero
  real(kind=fPrec), save :: clo0(2)    = zero
  real(kind=fPrec), save :: clop0(2)   = zero
  real(kind=fPrec), save :: ekk(2)     = zero
  real(kind=fPrec), save :: cr(mmul)   = zero
  real(kind=fPrec), save :: ci(mmul)   = zero
  real(kind=fPrec), save :: temptr(6)  = zero
  real(kind=fPrec), save :: clo6v(3)   = zero
  real(kind=fPrec), save :: clop6v(3)  = zero
  real(kind=fPrec), save :: tas(6,6)   = zero
  real(kind=fPrec), save :: tasau(6,6) = zero
  real(kind=fPrec), save :: qwcs(3)    = zero
  real(kind=fPrec), save :: di0xs      = zero
  real(kind=fPrec), save :: di0zs      = zero
  real(kind=fPrec), save :: dip0xs     = zero
  real(kind=fPrec), save :: dip0zs     = zero
  real(kind=fPrec), save :: di0au(4)
  real(kind=fPrec), save :: tau(6,6)
  real(kind=fPrec), save :: wx(3)
  real(kind=fPrec), save :: x1(6)
  real(kind=fPrec), save :: x2(6)
  real(kind=fPrec), save :: fake(2,20)
  real(kind=fPrec), save :: xau(2,6)
  real(kind=fPrec), save :: cloau(6)

  !  Arrays
  ! ========

  ! Number of Structure Elements (nblz)
  real(kind=fPrec), allocatable, save :: zsiv(:)       ! Displacement of elements, including error
  real(kind=fPrec), allocatable, save :: xsiv(:)       ! Displacement of elements, including error
  real(kind=fPrec), allocatable, save :: smiv(:)       ! Magnetic kick, including error

  ! Number of Particles (npart)
  real(kind=fPrec), allocatable, save :: xv1(:)        ! Transverse coordinates: Horisontal position
  real(kind=fPrec), allocatable, save :: yv1(:)        ! Transverse coordinates: Horisontal angle
  real(kind=fPrec), allocatable, save :: xv2(:)        ! Transverse coordinates: Vertical position
  real(kind=fPrec), allocatable, save :: yv2(:)        ! Transverse coordinates: Vertical angle
  real(kind=fPrec), allocatable, save :: sigmv(:)      ! Longitudinal coordinate: Position offset
  real(kind=fPrec), allocatable, save :: dpsv(:)       ! Longitudinal coordinate: Momentum offset from reference momentum e0f
  real(kind=fPrec), allocatable, save :: ejfv(:)       ! Particle momentum
  real(kind=fPrec), allocatable, save :: ejv(:)        ! Particle energy

  real(kind=fPrec), allocatable, save :: oidpsv(:)     ! 1/(1+dpsv)
  real(kind=fPrec), allocatable, save :: moidpsv(:)    ! Relative rigidity offset
  real(kind=fPrec), allocatable, save :: omoidpsv(:)   ! Relative rigidity offset
  real(kind=fPrec), allocatable, save :: rvv(:)        ! Beta_0 / Beta(j)
  real(kind=fPrec), allocatable, save :: ejf0v(:)      ! Temporary array for momentum updates
  real(kind=fPrec), allocatable, save :: dam(:)        ! Distance in phase space
  real(kind=fPrec), allocatable, save :: dpd(:)        ! Thick tracking only: 1+dpsv
  real(kind=fPrec), allocatable, save :: dpsq(:)       ! Thick tracking only: sqrt(1+dpsv)
  real(kind=fPrec), allocatable, save :: ampv(:)       ! Amplitude variations
  real(kind=fPrec), allocatable, save :: nucm(:)       ! Particle mass
  real(kind=fPrec), allocatable, save :: mtc(:)        ! Mass-to-charge ratio

  real(kind=fPrec), allocatable, save :: spin_x(:)     ! x component of the particle spin
  real(kind=fPrec), allocatable, save :: spin_y(:)     ! y component of the particle spin
  real(kind=fPrec), allocatable, save :: spin_z(:)     ! z component of the particle spin

  integer(kind=int16), allocatable, save :: nqq(:)     ! Particle charge
  integer(kind=int16), allocatable, save :: naa(:)     ! Ion atomic mass
  integer(kind=int16), allocatable, save :: nzz(:)     ! Ion atomic number
  integer(kind=int32), allocatable, save :: pdgid(:)   ! Particle PDGid

  integer,          allocatable, save :: numxv(:)      ! Turn in which a particle was lost

! The following variables are int32 for usage with the FLUKA IO TCP/IP communication
! If these are to be changed, remember to also update the FLUKA IO C code (and the CR variables).
! Also update the root output and geant4 interface (all fixed to int32s currently)
  integer(kind=int32), allocatable, save :: partID(:)   ! Particle ID
  integer(kind=int32), allocatable, save :: parentID(:) ! Particle parent ID in case of secondary particles
  integer(kind=int32), save :: MaximumPartID            ! Maximum used particle ID
#if defined(FLUKA) || defined(G4COLLIMATION)
  real(kind=fPrec),    allocatable, save :: partWeight(:) ! Particle weighting for FLUKA and geant4
#endif

  integer,          allocatable, save :: pairID(:,:)   ! The original particle pair ID for a particle
  integer,          allocatable, save :: pairMap(:,:)  ! A reverse map for pairID to index
  logical,          allocatable, save :: pstop(:)      ! Particle lost flag (post-processing)
  logical,          allocatable, save :: llostp(:)     ! Particle lost flag
  real(kind=fPrec), allocatable, save :: aperv(:,:)    ! Aperture at loss
  integer,          allocatable, save :: iv(:)         ! Entry in the sequence where loss occured

  real(kind=fPrec), allocatable, save :: bl1v(:,:,:,:) ! Thick tracking only: Transfer matrix for linear tracking (6,2,npart,nblo)

contains

subroutine mod_commonmn_expand_arrays(nblz_new,npart_new)

  use mod_alloc
  use mod_common, only : nucm0, aa0, zz0, qq0, pdgid0
  use numerical_constants, only : zero, one

  implicit none

  integer, intent(in) :: nblz_new
  integer, intent(in) :: npart_new

  integer :: nblz_prev  = -2
  integer :: npart_prev = -2
  integer npair_new

  npair_new = (npart_new+1)/2

  if(nblz_new /= nblz_prev) then
    call alloc(smiv,     nblz_new,     zero,    "smiv")
    call alloc(zsiv,     nblz_new,     zero,    "zsiv")
    call alloc(xsiv,     nblz_new,     zero,    "xsiv")
  end if

  if(npart_new /= npart_prev) then
    call alloc(xv1,        npart_new, zero,    "xv1")
    call alloc(yv1,        npart_new, zero,    "yv1")
    call alloc(xv2,        npart_new, zero,    "xv2")
    call alloc(yv2,        npart_new, zero,    "yv2")
    call alloc(sigmv,      npart_new, zero,    "sigmv")
    call alloc(dpsv,       npart_new, zero,    "dpsv")
    call alloc(ejv,        npart_new, zero,    "ejv")
    call alloc(ejfv,       npart_new, zero,    "ejfv")
    call alloc(dam,        npart_new, zero,    "dam")
    call alloc(rvv,        npart_new, one,     "rvv")
    call alloc(ejf0v,      npart_new, zero,    "ejf0v")
    call alloc(numxv,      npart_new, 0,       "numxv")
    call alloc(partID,     npart_new, 0,       "partID")
    call alloc(parentID,   npart_new, 0,       "parentID")
#if defined(FLUKA) || defined(G4COLLIMATION)
    call alloc(partWeight, npart_new, one,     "partWeight")
#endif
    call alloc(pairID,  2, npart_new, 0,       "pairID")
    call alloc(pairMap, 2, npair_new, 0,       "pairMap")
    call alloc(pstop,      npart_new, .false., "pstop")
    call alloc(llostp,     npart_new, .false., "llostp")
    call alloc(dpd,        npart_new, zero,    "dpd")
    call alloc(dpsq,       npart_new, zero,    "dpsq")
    call alloc(oidpsv,     npart_new, one,     "oidpsv")
    call alloc(moidpsv,    npart_new, one,     "moidpsv")
    call alloc(omoidpsv,   npart_new, zero,    "omoidpsv")
    call alloc(nucm,       npart_new, zero,    "nucm")
    call alloc(mtc,        npart_new, nucm0,   "mtc")
    call alloc(spin_x,     npart_new, zero,    "spin_x")
    call alloc(spin_y,     npart_new, zero,    "spin_y")
    call alloc(spin_z,     npart_new, zero,    "spin_z")
    call alloc(naa,        npart_new, aa0,     "naa")
    call alloc(nzz,        npart_new, zz0,     "nzz")
    call alloc(nqq,        npart_new, qq0,     "nqq")
    call alloc(pdgid,      npart_new, pdgid0,  "pdgid")
    call alloc(ampv,       npart_new, zero,    "ampv")
    call alloc(aperv,   2, npart_new, zero,    "aperv")
    call alloc(iv,         npart_new, 0,       "iv")
  end if

  nblz_prev  = nblz_new
  npart_prev = npart_new

end subroutine mod_commonmn_expand_arrays

subroutine mod_commonmn_expand_thickarrays(npart_new, nblo_new)

  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  integer,intent(in) :: npart_new, nblo_new

  call alloc(bl1v,6,2,npart_new,nblo_new,zero,"bl1v")

end subroutine mod_commonmn_expand_thickarrays

end module mod_common_main

! ================================================================================================ !
!  SOMETHING-SOMETHING COMMON VARIABLES
!  Last modified: 2018-06-13
! ================================================================================================ !
module mod_commons

  use parpro
  use floatPrecision
  use numerical_constants

  implicit none

  real(kind=fPrec), allocatable, save :: as(:,:,:,:) ! (6,2,npart,nele)
  real(kind=fPrec), allocatable, save :: al(:,:,:,:) ! (6,2,npart,nele)
  real(kind=fPrec), allocatable, save :: at(:,:,:,:) ! (6,2,npart,nele)
  real(kind=fPrec), allocatable, save :: a2(:,:,:,:) ! (6,2,npart,nele)

  real(kind=fPrec), save :: sigm(mpa) = zero
  real(kind=fPrec), save :: dps(mpa)  = zero
  real(kind=fPrec), save :: chi0      = zero  ! Starting phase of the initial coordinate
  real(kind=fPrec), save :: chid      = zero  ! Phase difference between first and second particles
  real(kind=fPrec), save :: exz(2,6)  = zero
  real(kind=fPrec), save :: dp1       = zero
  integer,          save :: idz(2)    = [1,0] ! Coupling on/off
  integer,          save :: itra      = 0     ! Number of particles

contains

subroutine mod_commons_expand_thickarrays(nele_new, npart_new)

  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  integer,intent(in) :: nele_new, npart_new

  call alloc(al,6,2,npart_new,nele_new,zero,"al")
  call alloc(as,6,2,npart_new,nele_new,zero,"as")
  call alloc(at,6,2,2,        nele_new,zero,"at")
  call alloc(a2,6,2,2,        nele_new,zero,"a2")

end subroutine mod_commons_expand_thickarrays

end module mod_commons

! ================================================================================================ !
!  MAIN COMMON VARIABLES
!  Last modified: 2018-06-13
! ================================================================================================ !
module mod_commond2

  use floatPrecision
  use parpro, only : nele, nema
  use crcoall
  use numerical_constants

  implicit none

  !Note: These only to be neccessary for thick4d?
  real(kind=fPrec), allocatable, save :: ald6(:,:,:,:), asd6(:,:,:,:) !(nele,2,6,nema)

contains

subroutine mod_commond2_expand_arrays(nele_new)

  use mod_alloc

  implicit none

  integer, intent(in) :: nele_new

  call alloc(ald6,nele_new,2,6,nema,zero,'ald6')
  call alloc(asd6,nele_new,2,6,nema,zero,'asd6')

end subroutine mod_commond2_expand_arrays

end module mod_commond2

! ================================================================================================ !
!  LIELIB AND DABNEW COMMON VARIABLES
!  Last modified: 2018-06-11
! ================================================================================================ !
module mod_lie_dab

  use floatPrecision

  implicit none

  ! From lielib
  integer, parameter :: ndim  = 3
  integer, parameter :: ndim2 = 6
  integer, parameter :: ntt   = 40
  integer, parameter :: nreso = 20

  integer,          save :: nd,nd2,no,nv
  integer,          save :: ifilt,idpr,iref,itu
  integer,          save :: lienot,iflow,jtune
  integer,          save :: ndc,ndc2,ndpt,ndt
  integer,          save :: nplane(ndim),ista(ndim),idsta(ndim)
  integer,          save :: mx(ndim,nreso),nres
  real(kind=fPrec), save :: epsplane,xplane(ndim)
  real(kind=fPrec), save :: sta(ndim),dsta(ndim),angle(ndim),radn(ndim)
  real(kind=fPrec), save :: ps(ndim),rads(ndim)
  real(kind=fPrec), save :: xintex(0:20)

  ! From dabnew
  integer, save      :: lda = -1
  integer, save      :: lst = -1
  integer, save      :: lea = -1
  integer, save      :: lia = -1
  integer, save      :: lno = -1
  integer, parameter :: lnv = 40

  integer,                       save :: ndat,nda,ndamaxi,lfi
  integer,                       save :: nst,nomax,nvmax,nmmax,nocut
  integer,          allocatable, save :: idano(:),idanv(:),idapo(:)
  integer,          allocatable, save :: idalm(:),idall(:)
  integer,          allocatable, save :: i1(:),i2(:)
  integer,          allocatable, save :: ie1(:),ie2(:),ieo(:),ifi(:)
  integer,          allocatable, save :: ia1(:),ia2(:)
  logical,          allocatable, save :: allvec(:)
  character(len=:), allocatable, save :: daname(:)
  real(kind=fPrec), allocatable, save :: cc(:),facint(:)
  real(kind=fPrec),              save :: eps,epsmac

  ! dascr variables
  integer,           save :: idao
  integer,           save :: iscrda(100)
  integer,           save :: iscrri(100)
  real(kind=fPrec),  save :: rscrri(100)

contains

subroutine mld_allocArrays(da_version)

  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  logical, intent(in) :: da_version

  if(da_version) then
    lda = 10000
    lst = 20050000
    lea = 100000
    lia = 10000000
    lno = 120
  else
    lda = 10000
    lst = 200000
    lea = 500
    lia = 10000
    lno = 120
  end if

  call alloc(idano,lda,0,"idano")
  call alloc(idanv,lda,0,"idanv")
  call alloc(idapo,lda,0,"idapo")
  call alloc(idalm,lda,0,"idalm")
  call alloc(idall,lda,0,"idall")

  call alloc(i1,lst,0,"i1")
  call alloc(i2,lst,0,"i2")

  call alloc(ie1,lea,0,"ie1")
  call alloc(ie2,lea,0,"ie2")
  call alloc(ieo,lea,0,"ieo")
  call alloc(ifi,lea,0,"ifi")

  call alloc(ia1,lia,0,"ia1",0)
  call alloc(ia2,lia,0,"ia2",0)

  call alloc(allvec,   lda,.false.,     "allvec")
  call alloc(daname,10,lda,"          ","daname")

  call alloc(cc,    lst,zero,"cc")
  call alloc(facint,lno,zero,"facint",0)

end subroutine mld_allocArrays

end module mod_lie_dab
