! ================================================================================================ !
! Standard SixTrack RNG
! Last modified: 2018-09-27
!
! References:
!
! F. James, A review of pseudorandom number generators,
!  Computer Physics Communications,
!  Volume 60, Issue 3, 1990, Pages 329-344, ISSN 0010-4655,
!  https://doi.org/10.1016/0010-4655(90)90032-V.
!  http://www.sciencedirect.com/science/article/pii/001046559090032V
!
! P. L'Ecuyer. 1988. Efficient and portable combined random number generators.
!  Commun. ACM 31, 6 (June 1988), 742-751. DOI=10.1145/62959.62969
!  http://doi.acm.org/10.1145/62959.62969
!
! ================================================================================================ !
module mod_ranecu

  implicit none

  integer, private, save :: iseed1 = 12345
  integer, private, save :: iseed2 = 67890

contains

! Generate LEN random numbers into RVEC
!  If mcut ==  0, generate normal distributed random numbers
!  If mcut  >  0, generate normal distributed random numbers r <= mcut
!  If mcut == -1, generate uniformly distributed random numbers
subroutine ranecu(rvec,len,mcut)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  implicit none

  real(kind=fPrec), intent(out), dimension(*) :: rvec
  integer, intent(in) :: len, mcut

  integer i,iz,j,k
  real(kind=fPrec) rvec0
  real(kind=fPrec), dimension(2) :: r

  i=1

10 continue
  do j = 1,2
    k = iseed1/53668
    iseed1 = 40014*(iseed1-k*53668) - k*12211
    if (iseed1 < 0) iseed1 = iseed1+2147483563
    k = iseed2/52774
    iseed2 = 40692*(iseed2-k*52774) - k*3791
    if (iseed2 < 0) iseed2 = iseed2+2147483399
    iz = iseed1-iseed2
    if (iz < 1) iz = iz+2147483562
    r(j) = real(iz,fPrec)*4.656613e-10_fPrec
  end do

  if (mcut >= 0) then ! mcut = -1 => Generate uniform numbers!
    ! Convert r(1), r(2) from U(0,1) -> rvec0 as Gaussian with cutoff mcut (#sigmas):
    rvec0 = sqrt(((-one*two)*log_mb(r(1))))*cos_mb((two*pi)*r(2))
  else if (mcut == -1) then
    rvec0 = r(1)
  end if

  if(abs(rvec0) <= real(mcut,fPrec) .or. mcut == 0 .or. mcut == -1) then
    rvec(i) = rvec0
    i=i+1
  end if
  if(i <= len) goto 10
  return

end subroutine ranecu

! Set the seeds
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
