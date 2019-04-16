! ================================================================================================ !
! Shared module for general physical constants
! Last modified: 2018-05-25
! See: http://pdg.lbl.gov/2017/reviews/rpp2016-rev-phys-constants.pdf
! Some values are not updated to ensure numerical compatability with older studies
! ================================================================================================ !
module physical_constants

  use floatPrecision

  implicit none

#ifndef FLUKA
  ! Proton mass (MeV)
  real(kind=fPrec), parameter :: pmap   = 938.271998_fPrec              ! old
! real(kind=fPrec), parameter :: pmap   = 938.2720813_fPrec             ! 2017

  ! Electron mass (MeV) from PDG, 2002
  real(kind=fPrec), parameter :: pmae   = 0.510998902_fPrec             ! old
! real(kind=fPrec), parameter :: pmae   = 0.5109989461_fPrec            ! 2017
#else
  ! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
  ! Last modified: 08-12-2014
  ! Synch masses of proton and electron to values used by FLUKA
! real(kind=fPrec), parameter :: pmap  = 0.938272046e3_fPrec            ! PDG 2014, Fluka-2011-dev
  real(kind=fPrec), parameter :: pmap  = 0.938272310e3_fPrec            ! PDG 20xx, Fluka-2011-pro
  real(kind=fPrec), parameter :: pmae  = 0.510998928_fPrec              ! PDG 2014, Fluka-2011-dev
#endif

  ! Classical electron radius
  real(kind=fPrec), parameter :: crade  = 2.817940285e-15_fPrec         ! old
! real(kind=fPrec), parameter :: crade  = 2.8179403227e-15_fPrec        ! 2017

  ! Speed of light
  real(kind=fPrec), parameter :: clight = 2.99792458e8_fPrec            ! Exact

  ! Avogadro constant
  real(kind=fPrec), parameter :: fnavo  = 6.022140857e23_fPrec          ! 2017
! real(kind=fPrec), parameter :: fnavo  = 6.02214129e23_fPrec           ! old

  ! Planck constant
  real(kind=fPrec), parameter :: planck = 6.626070040e-34_fPrec

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

  real(kind=fPrec), parameter :: eulergamma = 0.577215664901532860606512090082402431042159335939923598805_fPrec
  real(kind=fPrec), parameter :: pi         = 3.141592653589793238462643383279502884197169399375105820974_fPrec
  real(kind=fPrec), parameter :: inv_ln2    = 1.442695040888963407359924681001892137426645954152985934135_fPrec

  real(kind=fPrec), parameter :: pi2    = 0.5_fPrec*pi
  real(kind=fPrec), parameter :: twopi  = 2.0_fPrec*pi
  real(kind=fPrec), parameter :: pisqrt = sqrt(pi)
  real(kind=fPrec), parameter :: rad    = pi/180.0_fPrec

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
  real(kind=fPrec), parameter :: c1e3   = 1.0e3_fPrec
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

  real(kind=fPrec), parameter :: c2e3   = 2.0e3_fPrec
  real(kind=fPrec), parameter :: c4e3   = 4.0e3_fPrec
  real(kind=fPrec), parameter :: c180e0 = 180.0_fPrec

  real(kind=fPrec), parameter :: c1m1   = 1.0e-1_fPrec
  real(kind=fPrec), parameter :: c1m2   = 1.0e-2_fPrec
  real(kind=fPrec), parameter :: c1m3   = 1.0e-3_fPrec
  real(kind=fPrec), parameter :: c1m4   = 1.0e-4_fPrec
  real(kind=fPrec), parameter :: c5m4   = 5.0e-4_fPrec
  real(kind=fPrec), parameter :: c1m5   = 1.0e-5_fPrec
  real(kind=fPrec), parameter :: c1m6   = 1.0e-6_fPrec
  real(kind=fPrec), parameter :: c1m7   = 1.0e-7_fPrec
  real(kind=fPrec), parameter :: c1m8   = 1.0e-8_fPrec
  real(kind=fPrec), parameter :: c1m9   = 1.0e-9_fPrec
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
