! ================================================================================================ !
! Shared module for general physical constants
! Last modified: 2018-05-25
! See: http://pdg.lbl.gov/2017/reviews/rpp2016-rev-phys-constants.pdf
! Some values are not updated to ensure numerical compatability with older studies
! ================================================================================================ !
module physical_constants

  use floatPrecision

  implicit none

#ifdef FLUKA
  ! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
  ! Last modified: 08-12-2014
  ! Synch masses of proton and electron to values used by FLUKA
! real(kind=fPrec), parameter :: pmap  = 0.938272046e3_fPrec            ! PDG 2014, Fluka-2011-dev
  real(kind=fPrec), parameter :: pmap  = 0.938272310e3_fPrec            ! PDG 20xx, Fluka-2011-pro
  real(kind=fPrec), parameter :: pmae  = 0.510998928_fPrec              ! PDG 2014, Fluka-2011-dev
#else
  ! Proton and Electron Mass (MeV)
  real(kind=fPrec), parameter :: pmap  = 938.271998_fPrec               ! Old
  real(kind=fPrec), parameter :: pmae  = 0.510998902_fPrec              ! PDG 2002
#endif

  ! Proton and Electron Mass (MeV) Latest Values
  real(kind=fPrec), parameter :: pmap_18 = 938.272081_fPrec             ! PDG 2018
  real(kind=fPrec), parameter :: pmae_18 = 0.5109989461_fPrec           ! PDG 2018

  ! Classical electron radius
  real(kind=fPrec), parameter :: crade    = 2.817940285e-15_fPrec       ! Old
  real(kind=fPrec), parameter :: crade_18 = 2.8179403227e-15_fPrec      ! PDG 2018

  ! Speed of light
  real(kind=fPrec), parameter :: clight = 2.99792458e8_fPrec            ! Exact by definition

  ! Avogadro constant
! real(kind=fPrec), parameter :: fnavo  = 6.02214129e23_fPrec           ! Old
! real(kind=fPrec), parameter :: fnavo  = 6.022140857e23_fPrec          ! PDG 2018
  real(kind=fPrec), parameter :: fnavo  = 6.02214076e23_fPrec           ! Exact by definition

  ! Planck constant
! real(kind=fPrec), parameter :: planck = 6.626070040e-34_fPrec
  real(kind=fPrec), parameter :: planck = 6.62607015e-34_fPrec          ! Exact by definition

  !Electric charge
  real(kind=fPrec), parameter :: echarge = 1.602176634e-19_fPrec        ! Exact by definition

  ! Vacuum permittivity
  real(kind=fPrec), parameter :: eps0   = 8.854187817e-12_fPrec

end module physical_constants

! ================================================================================================ !
! Shared module for numerical constants
! Last modified: 2018-05-25
! To ensure that the same values (including round-off) is used throughout the program.
! ================================================================================================ !
module numerical_constants

  use floatPrecision

  implicit none

! real(kind=fPrec), parameter :: pi      = 3.141592653589793238462643383279502884197169399375105820974_fPrec
! real(kind=fPrec), parameter :: inv_ln2 = 1.442695040888963407359924681001892137426645954152985934135_fPrec

  ! We set the parameters of irrational numbers explicitly to their nearest binary representation
  ! without depending on rounding from decimal representation. These use decimal round near.
  ! Small program for checking binary and hex here: https://gist.github.com/vkbo/0e5d0c4b7b1ecb533a0486c79a30f741
  ! Note: Transfer to real128 only works for GNU

#ifdef SINGLE_MATH
  real(kind=fPrec), parameter :: pi      = real(z'40490fdb',fPrec)
  real(kind=fPrec), parameter :: pi2     = real(z'3fc90fdb',fPrec) ! 0.5_fPrec*pi
  real(kind=fPrec), parameter :: twopi   = real(z'40c90fdb',fPrec) ! 2.0_fPrec*pi
  real(kind=fPrec), parameter :: pisqrt  = real(z'3fe2dfc5',fPrec) ! sqrt(pi)
  real(kind=fPrec), parameter :: inv_ln2 = real(z'3fb8aa3b',fPrec) ! 1/log(2)
  real(kind=fPrec), parameter :: rad     = real(z'3c8efa35',fPrec) ! pi/180.0_fPrec
#endif
#ifdef DOUBLE_MATH
  real(kind=fPrec), parameter :: pi      = real(z'400921fb54442d18',fPrec)
  real(kind=fPrec), parameter :: pi2     = real(z'3ff921fb54442d18',fPrec)
  real(kind=fPrec), parameter :: twopi   = real(z'401921fb54442d18',fPrec)
  real(kind=fPrec), parameter :: pisqrt  = real(z'3ffc5bf891b4ef6a',fPrec)
  real(kind=fPrec), parameter :: inv_ln2 = real(z'3ff71547652b82fe',fPrec)
  real(kind=fPrec), parameter :: rad     = real(z'3f91df46a2529d39',fPrec)
#endif
#ifdef QUAD_MATH
#ifdef GFORTRAN
  real(kind=fPrec), parameter :: pi      = real(z'4000921fb54442d18469898cc51701b8',fPrec)
  real(kind=fPrec), parameter :: pi2     = real(z'3fff921fb54442d18469898cc51701b8',fPrec)
  real(kind=fPrec), parameter :: twopi   = real(z'4001921fb54442d18469898cc51701b8',fPrec)
  real(kind=fPrec), parameter :: pisqrt  = real(z'3fffc5bf891b4ef6aa79c3b0520d5db9',fPrec)
  real(kind=fPrec), parameter :: inv_ln2 = real(z'3fff71547652b82fe1777d0ffda0d23a',fPrec)
  real(kind=fPrec), parameter :: rad     = real(z'3ff91df46a2529d3915c1d8becdd290b',fPrec)
#else
  real(kind=fPrec), parameter :: pi      = 3.141592653589793238462643383279502884197169399375105820974_fPrec
  real(kind=fPrec), parameter :: pi2     = 0.5_fPrec*pi
  real(kind=fPrec), parameter :: twopi   = 2.0_fPrec*pi
  real(kind=fPrec), parameter :: pisqrt  = sqrt(pi)
  real(kind=fPrec), parameter :: inv_ln2 = 1.442695040888963407359924681001892137426645954152985934135_fPrec
  real(kind=fPrec), parameter :: rad     = pi/180.0_fPrec
#endif
#endif

  ! The variable pieni is roughly the smallest normal single precision float value.
  ! However, 1e-38 is actually subnormal in 32 bit (0x006ce3ee). Smallest normal is 1.17549435e-38 (0x00800000)
  ! This is probably not a good idea for performance in single prec, but changing it makes a difference in postprocessing it seems.
  real(kind=fPrec), parameter :: pieni  = 1e-38_fPrec

  real(kind=fPrec), parameter :: zero   = 0.0_fPrec
  real(kind=fPrec), parameter :: half   = 0.5_fPrec
  real(kind=fPrec), parameter :: one    = 1.0_fPrec
  real(kind=fPrec), parameter :: two    = 2.0_fPrec
  real(kind=fPrec), parameter :: three  = 3.0_fPrec
  real(kind=fPrec), parameter :: four   = 4.0_fPrec
  real(kind=fPrec), parameter :: five   = 5.0_fPrec
  real(kind=fPrec), parameter :: six    = 6.0_fPrec
  real(kind=fPrec), parameter :: seven  = 7.0_fPrec
  real(kind=fPrec), parameter :: eight  = 8.0_fPrec
  real(kind=fPrec), parameter :: nine   = 9.0_fPrec

  real(kind=fPrec), parameter :: c1e1   = 1.0e1_fPrec
  real(kind=fPrec), parameter :: c1e2   = 1.0e2_fPrec
  real(kind=fPrec), parameter :: c2e2   = 2.0e2_fPrec
  real(kind=fPrec), parameter :: c1e3   = 1.0e3_fPrec
  real(kind=fPrec), parameter :: c2e3   = 2.0e3_fPrec
  real(kind=fPrec), parameter :: c4e3   = 4.0e3_fPrec
  real(kind=fPrec), parameter :: c1e4   = 1.0e4_fPrec
  real(kind=fPrec), parameter :: c1e5   = 1.0e5_fPrec
  real(kind=fPrec), parameter :: c1e6   = 1.0e6_fPrec
  real(kind=fPrec), parameter :: c1e7   = 1.0e7_fPrec
  real(kind=fPrec), parameter :: c1e8   = 1.0e8_fPrec
  real(kind=fPrec), parameter :: c1e9   = 1.0e9_fPrec
  real(kind=fPrec), parameter :: c1e10  = 1.0e10_fPrec
  real(kind=fPrec), parameter :: c1e12  = 1.0e12_fPrec
  real(kind=fPrec), parameter :: c1e13  = 1.0e13_fPrec
  real(kind=fPrec), parameter :: c1e15  = 1.0e15_fPrec
  real(kind=fPrec), parameter :: c1e16  = 1.0e16_fPrec
  real(kind=fPrec), parameter :: c1e27  = 1.0e27_fPrec
  real(kind=fPrec), parameter :: c1e38  = 1.0e38_fPrec

  real(kind=fPrec), parameter :: c180e0 = 1.8e2_fPrec

  real(kind=fPrec), parameter :: c1m1   = 1.0e-1_fPrec
  real(kind=fPrec), parameter :: c2m1   = 2.0e-1_fPrec
  real(kind=fPrec), parameter :: c1m2   = 1.0e-2_fPrec
  real(kind=fPrec), parameter :: c5m2   = 5.0e-2_fPrec
  real(kind=fPrec), parameter :: c1m3   = 1.0e-3_fPrec
  real(kind=fPrec), parameter :: c1m4   = 1.0e-4_fPrec
  real(kind=fPrec), parameter :: c5m4   = 5.0e-4_fPrec
  real(kind=fPrec), parameter :: c1m5   = 1.0e-5_fPrec
  real(kind=fPrec), parameter :: c1m6   = 1.0e-6_fPrec
  real(kind=fPrec), parameter :: c2m6   = 2.0e-6_fPrec
  real(kind=fPrec), parameter :: c1m7   = 1.0e-7_fPrec
  real(kind=fPrec), parameter :: c1m8   = 1.0e-8_fPrec
  real(kind=fPrec), parameter :: c1m9   = 1.0e-9_fPrec
  real(kind=fPrec), parameter :: c2m9   = 2.0e-9_fPrec
  real(kind=fPrec), parameter :: c1m10  = 1.0e-10_fPrec
  real(kind=fPrec), parameter :: c1m12  = 1.0e-12_fPrec
  real(kind=fPrec), parameter :: c1m13  = 1.0e-13_fPrec
  real(kind=fPrec), parameter :: c1m15  = 1.0e-15_fPrec
  real(kind=fPrec), parameter :: c1m18  = 1.0e-18_fPrec
  real(kind=fPrec), parameter :: c1m21  = 1.0e-21_fPrec
  real(kind=fPrec), parameter :: c1m24  = 1.0e-24_fPrec
  real(kind=fPrec), parameter :: c1m27  = 1.0e-27_fPrec
  real(kind=fPrec), parameter :: c1m36  = 1.0e-36_fPrec
  real(kind=fPrec), parameter :: c1m38  = 1.0e-38_fPrec

end module numerical_constants
