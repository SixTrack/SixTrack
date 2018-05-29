! M. Fitterer, FNAL, A. Mereghtti, CERN
! last modified: 09-02-2018
! Common block for electron lens definition
module elens
  
  use parpro
  use floatPrecision
  use crcoall
  use mod_alloc
  
  implicit none
  
  ! size of table with elens data
  integer, parameter     :: nelens=125
  ! last elens read
  integer, save          :: melens

  ! index of elens:
  integer,allocatable, save          :: ielens(:) !(nele)

  ! variables to save elens parameters for tracking etc.
  integer, save          :: elens_type(nelens)      ! integer for elens type
                                                    ! 0 : Un-initialized.
                                                    ! 1 : Hollow annular elens, uniform profile
  real(kind=fPrec), save :: elens_theta_r2(nelens)    ! kick strength at R2 [mrad]
  real(kind=fPrec), save :: elens_r2(nelens)          ! outer radius R2 [mm]
  real(kind=fPrec), save :: elens_r1(nelens)          ! inner radius R1 [mm]
  real(kind=fPrec), save :: elens_offset_x(nelens), elens_offset_y(nelens)  ! hor./vert. offset of elens [mm]
  real(kind=fPrec), save :: elens_sig(nelens)         ! sig (Gaussian profile) [mm]
  real(kind=fPrec), save :: elens_geo_norm(nelens)    ! normalisation of f(r)
  real(kind=fPrec), save :: elens_len(nelens)         ! length of eLens (e-beam region) [m]
  real(kind=fPrec), save :: elens_I(nelens)           ! current of e-beam [A]
  real(kind=fPrec), save :: elens_Ek(nelens)          ! kinetic energy of e-beam [keV]
  logical, save          :: elens_lThetaR2(nelens)    ! flag for computing theta@R2
  integer, save          :: elens_iCheby(nelens)      ! mapping to the table with chebyshev coeffs
  real(kind=fPrec), save :: elens_cheby_angle(nelens) ! angle for getting the real bends [deg]
  ! file with chebyshev coefficients
  integer, parameter     :: nelens_cheby_tables=20    ! number of tables with chebyshev coefficients
  integer, parameter     :: elens_cheby_unit=107      ! unit for reading the chebyshev coefficients
  integer, parameter     :: elens_cheby_order=18      ! max order of chebyshev polynomials
  integer, save          :: melens_cheby_tables       ! tables available in memory
  character(len=16), save:: elens_cheby_filename(nelens_cheby_tables) ! names
  real(kind=fPrec), save :: elens_cheby_coeffs(0:elens_cheby_order,0:elens_cheby_order,nelens_cheby_tables)
  real(kind=fPrec), save :: elens_cheby_refCurr(nelens_cheby_tables) ! reference current [A]
  real(kind=fPrec), save :: elens_cheby_refRadius(nelens_cheby_tables) ! reference radius [mm]
  real(kind=fPrec), save :: elens_cheby_refBeta(nelens_cheby_tables) ! reference e-beta []
  
contains

subroutine elens_allocate_arrays
  use crcoall
  implicit none
  integer stat
    call alloc(ielens,nele,0,'ielens')
end subroutine elens_allocate_arrays

subroutine elens_expand_arrays(nele_new)
  implicit none
  integer, intent(in) :: nele_new
  call resize(ielens,nele_new,0,'ielens')
end subroutine elens_expand_arrays

end module elens
