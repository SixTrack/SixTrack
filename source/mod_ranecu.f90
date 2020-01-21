! ================================================================================================ !
! Standard SixTrack RNG
! Last modified: 2018-09-27
!
! References:
!
! [1] F. James, A review of pseudorandom number generators,
!     Computer Physics Communications,
!     Volume 60, Issue 3, 1990, Pages 329-344, ISSN 0010-4655,
!     https://doi.org/10.1016/0010-4655(90)90032-V.
!     http://www.sciencedirect.com/science/article/pii/001046559090032V
!
! [2] P. L'Ecuyer. 1988. Efficient and portable combined random number generators.
!     Commun. ACM 31, 6 (June 1988), 742-751. DOI=10.1145/62959.62969
!     http://doi.acm.org/10.1145/62959.62969
!
! Usage:
!
! recuinit  This routine is meant to be used for initialising the random number series based on user
!           input from input files. It has the additional boundary checks that recuin does not have,
!           and in addition also allows for setting only one seed, to which a value is added to
!           generate the second seed.
! recuin    Takes two seeds as input to continue a random number series already initialised.
! recuut    Returns the two seeds of a random number series so that it can be continued later.
! ================================================================================================ !
module mod_ranecu

  use floatPrecision

  implicit none

  integer, private, save :: iseed1 = 12345
  integer, private, save :: iseed2 = 67890

  ! Constants for converting the integer range 1 to 2147483563 to a real between 0 and 1.0 (non-inclusive)
#ifdef SINGLE_MATH
  real(kind=fPrec), parameter :: rScale = 4.65661287e-10_fPrec ! 0x30000000
#endif
#ifdef DOUBLE_MATH
  real(kind=fPrec), parameter :: rScale = 4.65661305739176810e-10_fPrec ! 0x3e0000000aa00006
#endif
#ifdef QUAD_MATH
  real(kind=fPrec), parameter :: rScale = 4.6566130573917691960466940253828613e-10_fPrec ! 0x3fe0000000aa000070e4004af76831c8
#endif

contains

! Generate LEN random numbers into RVEC with optional sigma cut if normal
!   mode 0 : Generate uniformly distributed random numbers
!   mode 1 : Generate normal distributed random numbers
subroutine ranecu(rvec, len, mode, cut)

  use crcoall
  use mathlib_bouncer
  use numerical_constants

  integer,                    intent(in)  :: len
  real(kind=fPrec),           intent(out) :: rvec(len)
  integer,                    intent(in)  :: mode
  real(kind=fPrec), optional, intent(in)  :: cut

  real(kind=fPrec) rvec0, r(2), docut
  integer i, iz, j, k

  if(present(cut)) then
    docut = cut
  else
    docut = zero
  end if

  i = 1
  rvec0 = zero
  if(mode /= 0 .and. mode /= 1) then
    write(lerr,"(a,i0)") "RANECU> ERROR mode must be 0 (uniform) or 1 (normal), got ", mode
    call prror
  end if

10 continue
  do j = 1,2
    k = iseed1/53668
    iseed1 = 40014*(iseed1-k*53668) - k*12211
    if(iseed1 < 0) iseed1 = iseed1+2147483563
    k = iseed2/52774
    iseed2 = 40692*(iseed2-k*52774) - k*3791
    if(iseed2 < 0) iseed2 = iseed2+2147483399
    iz = iseed1-iseed2
    if(iz < 1) iz = iz+2147483562
    r(j) = real(iz,fPrec)*4.656613e-10_fPrec ! Note: this is still single precision
  end do

  if(mode == 1) then
    ! Convert r(1), r(2) from U(0,1) -> rvec0 as Gaussian with cutoff mcut (#sigmas):
    rvec0 = sqrt((-one*two)*log_mb(r(1))) * cos_mb(twopi*r(2))
  else
    rvec0 = r(1)
  end if

  if(mode == 0 .or. docut <= zero .or. abs(rvec0) <= docut) then
    rvec(i) = rvec0
    i = i + 1
  end if

  if(i <= len) goto 10

end subroutine ranecu

! Uniform-only version of the above
subroutine ranecuu(rvec, len)

  integer,          intent(in)  :: len
  real(kind=fPrec), intent(out) :: rvec(len)

  integer iz, j, k
  do j=1,len
    k = iseed1/53668
    iseed1 = 40014*(iseed1-k*53668) - k*12211
    if(iseed1 < 0) iseed1 = iseed1+2147483563
    k = iseed2/52774
    iseed2 = 40692*(iseed2-k*52774) - k*3791
    if(iseed2 < 0) iseed2 = iseed2+2147483399
    iz = iseed1-iseed2
    if(iz < 1) iz = iz+2147483562
    rvec(j) = real(iz,fPrec)*rScale
  end do

end subroutine ranecuu

! Init the random generator, that is, set seeds with proper value checks as described in [2]
subroutine recuinit(is1,is2)

  use crcoall

  integer,           intent(in) :: is1
  integer, optional, intent(in) :: is2

  if(is1 < 1 .or. is1 > 2147483562) then
    write(lerr,"(a,i0)") "RANECU> ERROR Seed 1 must be an integer in the range 1 to 2147483562, got ", is1
    call prror
  else
    iseed1 = is1
  end if

  if(present(is2)) then
    if(is2 < 1 .or. is2 > 2147483398) then
      write(lerr,"(a,i0)") "RANECU> ERROR Seed 2 must be an integer in the range 1 to 2147483398, got ", is2
      call prror
    else
      iseed2 = is2
    end if
  else
    iseed2 = iseed1 + 19770404
    if(iseed2 > 2147483398) then
      ! Subtract max of seed 1 to ensure lower than max of seed 2
      iseed2 = iseed2 - 2147483562
    end if
  end if

end subroutine recuinit

! Set the seeds
! The recuin routine to be used for continuing a random number stream with seeds extracted using recuut()
subroutine recuin(is1,is2)
  integer, intent(in) :: is1, is2
  iseed1 = is1
  iseed2 = is2
end subroutine recuin

! Get the current seeds
subroutine recuut(is1,is2)
  integer, intent(out) :: is1,is2
  is1 = iseed1
  is2 = iseed2
end subroutine recuut

end module mod_ranecu
