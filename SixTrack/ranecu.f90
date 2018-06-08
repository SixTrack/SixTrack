! Standard SixTrack RNG
! References:
!
! F. James, A review of pseudorandom number generators,
!  Computer Physics Communications,
!  Volume 60, Issue 3, 1990, Pages 329-344, ISSN 0010-4655,
!  https://doi.org/10.1016/0010-4655(90)90032-V.
!  http://www.sciencedirect.com/science/article/pii/001046559090032V
!
!
! P. L'Ecuyer. 1988. Efficient and portable combined random number generators.
!  Commun. ACM 31, 6 (June 1988), 742-751. DOI=10.1145/62959.62969
!  http://doi.acm.org/10.1145/62959.62969
!
module mod_ranecu
  implicit none

  integer :: iseed1 = 12345
  integer :: iseed2 = 67890

  save iseed1, iseed2
  
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
      real(kind=fPrec) rvec0, pi
      real(kind=fPrec), dimension(2) :: r
      
!-----------------------------------------------------------------------
      pi = four*atan_mb(one)
!     DO 100 I = 1,LEN
      i=1

10    continue
      do j = 1,2
        k = iseed1/53668
        iseed1 = 40014*(iseed1-k*53668) - k*12211
        if (iseed1.lt.0) iseed1 = iseed1+2147483563
        k = iseed2/52774
        iseed2 = 40692*(iseed2-k*52774) - k*3791
        if (iseed2.lt.0) iseed2 = iseed2+2147483399
        iz = iseed1-iseed2
        if (iz.lt.1) iz = iz+2147483562
        r(j) = real(iz,fPrec)*4.656613e-10_fPrec                                     !hr05
     end do

      if (mcut.ge.0) then !mcut = -1 => Generate uniform numbers!
!     Convert r(1), r(2) from U(0,1) -> rvec0 as Gaussian with cutoff mcut (#sigmas):
         rvec0 = sqrt(((-one*two)*log_mb(r(1))))*cos_mb((two*pi)*r(2))      !hr05
      else if (mcut.eq.-1) then
         rvec0 = r(1)
      end if
      
      if(abs(rvec0).le.real(mcut,fPrec).or.mcut.eq.0.or.mcut.eq.-1) then
        rvec(i) = rvec0
        i=i+1
      endif
      if(i.le.len) goto 10
!     RVEC(I) = ((-TWO*LOG(R(1)))**HALF)*COS(TWO*PI*R(2))
! 100 CONTINUE
      return

    end subroutine ranecu

    ! Set the seeds
    subroutine recuin(is1,is2)
      implicit none
      integer, intent(in) :: is1, is2
      
      iseed1 = is1
      iseed2 = is2
    end subroutine recuin

    ! Get the current seeds
    subroutine recuut(is1,is2)
      implicit none
      integer, intent(out) :: is1,is2
      
      is1 = iseed1
      is2 = iseed2
    end subroutine recuut
    
    
  end module mod_ranecu
