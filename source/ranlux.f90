!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module mod_ranlux
  use floatPrecision
  use numerical_constants
  use crcoall

  implicit none

  private

  integer jseed,lp,i,k,inner,izip,izip2,ivec,isk,isd,ilx,iouter

  integer, parameter :: maxlev = 4
  integer, parameter :: lxdflt = 3
  integer, parameter :: igiga = 1000000000
  integer, parameter :: itwo24 = 2**24
  integer, parameter :: jsdflt = 314159265
  integer, parameter :: icons = 2147483563

  integer, save :: ndskip(0:maxlev) = [0, 24, 73, 199, 365]
  integer, save :: isdext(25)
  integer, save :: iseeds(24)
  integer, save :: next(24)

  real(kind=fPrec), parameter :: twop12 = 4096
  real(kind=fPrec), save :: seeds(24)
  real(kind=fPrec), save :: twom12
  real(kind=fPrec), save :: twom24
  real(kind=fPrec), save :: carry = zero
  real(kind=fPrec) :: uni

  integer, save :: luxlev = lxdflt
  logical, save :: notyet = .true.
  integer, save :: in24 = 0
  integer, save :: i24 = 24
  integer, save :: j24 = 10
  integer, save :: kount = 0
  integer, save :: mkount = 0
  integer, save :: inseed
  integer, save :: nskip

!  data notyet, luxlev, in24, kount, mkount /.true., lxdflt, 0,0,0/
!  data in24, kount, mkount / 0,0,0/
!  data i24,j24,carry/24,10,0./

!  save i24, j24, carry, seeds, twom24, twom12, luxlev
!  save nskip, ndskip, in24, next, kount, mkount, inseed

!                               default
!  Luxury Level   0     1     2   *3*    4
!      data ndskip/0,   24,   73,  199,  365 /
!Corresponds to p=24    48    97   223   389
!     time factor 1     2     3     6    10   on slow workstation
!                 1    1.5    2     3     5   on fast mainframe
!

  public ranlux
  public rluxin
  public rluxut
  public rluxat
  public rluxgo
  public rndm4
  public rndm5
  public ran_gauss

contains

subroutine ranlux(rvec,lenv)
!         Subtract-and-borrow random number generator proposed by
!         Marsaglia and Zaman, implemented by F. James with the name
!         RCARRY in 1991, and later improved by Martin Luescher
!         in 1993 to produce "Luxury Pseudorandom Numbers".
!     Fortran 77 coded by F. James, 1993
!
!   LUXURY LEVELS.
!   ------ ------      The available luxury levels are:
!
!  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
!           and Zaman, very long period, but fails many tests.
!  level 1  (p=48): considerable improvement in quality over level 0,
!           now passes the gap test, but still fails spectral test.
!  level 2  (p=97): passes all known tests, but theoretically still
!           defective.
!  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
!           correlations have very small chance of being observed.
!  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
!
!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!  Calling sequences for RANLUX:                                  ++
!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
!!!                   32-bit random floating point numbers between  ++
!!!                   zero (not included) and one (also not incl.). ++
!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
!!!               which is integer between zero and MAXLEV, or if   ++
!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
!!!               should be set to zero unless restarting at a break++
!!!               point given by output of RLUXAT (see RLUXAT).     ++
!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
!!!               which can be used to restart the RANLUX generator ++
!!!               at the current point by calling RLUXGO.  K1 and K2++
!!!               specify how many numbers were generated since the ++
!!!               initialization with LUX and INT.  The restarting  ++
!!!               skips over  K1+K2*E9   numbers, so it can be long.++
!!!   A more efficient but less convenient way of restarting is by: ++
!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
!!!                 32-bit integer seeds, to be used for restarting ++
!!!      ISVEC must be dimensioned 25 in the calling program        ++
!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

  integer :: lenv
  real(kind=fPrec) :: rvec(lenv)

!  NOTYET is .TRUE. if no initialization has been performed yet.
!              Default Initialization by Multiplicative Congruential
      if (notyet) then
         notyet = .false.
         jseed = jsdflt
         inseed = jseed
         write(lout,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',jseed
         luxlev = lxdflt
         nskip = ndskip(luxlev)
         lp = nskip + 24
         in24 = 0
         kount = 0
         mkount = 0
!         WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
!     &        LUXLEV,'      p =',LP
            twom24 = 1.
         do 25 i= 1, 24
            twom24 = twom24 * 0.5
         k = jseed/53668
         jseed = 40014*(jseed-k*53668) -k*12211
         if (jseed .lt. 0)  jseed = jseed+icons
         iseeds(i) = mod(jseed,itwo24)
   25    continue
         twom12 = twom24 * 4096.
         do 50 i= 1,24
         seeds(i) = real(iseeds(i))*twom24
         next(i) = i-1
   50    continue
         next(1) = 24
         i24 = 24
         j24 = 10
         carry = 0.
         if (seeds(24) .eq. 0.) carry = twom24
      endif
!
!          The Generator proper: "Subtract-with-borrow",
!          as proposed by Marsaglia and Zaman,
!          Florida State University, March, 1989
!
      do 100 ivec= 1, lenv
      uni = seeds(j24) - seeds(i24) - carry
      if (uni .lt. 0.)  then
         uni = uni + 1.
         carry = twom24
      else
         carry = 0.
      endif
      seeds(i24) = uni
      i24 = next(i24)
      j24 = next(j24)
      rvec(ivec) = uni
!  small numbers (with less than 12 "significant" bits) are "padded".
      if (uni .lt. twom12)  then
         rvec(ivec) = rvec(ivec) + twom24*seeds(j24)
!        and zero is forbidden in case someone takes a logarithm
         if (rvec(ivec) .eq. 0.)  rvec(ivec) = twom24*twom24
      endif
!        Skipping to luxury.  As proposed by Martin Luscher.
      in24 = in24 + 1
      if (in24 .eq. 24)  then
         in24 = 0
         kount = kount + nskip
         do 90 isk= 1, nskip
         uni = seeds(j24) - seeds(i24) - carry
         if (uni .lt. 0.)  then
            uni = uni + 1.
            carry = twom24
         else
            carry = 0.
         endif
         seeds(i24) = uni
         i24 = next(i24)
         j24 = next(j24)
   90    continue
      endif
  100 continue
      kount = kount + lenv
      if (kount .ge. igiga)  then
         mkount = mkount + 1
         kount = kount - igiga
      endif
      return
end subroutine ranlux

!           Entry to input and float integer seeds from previous run
subroutine rluxin(isdext_tmp)

  implicit none

  integer, intent(in) :: isdext_tmp(25)
  isdext = isdext_tmp
         notyet = .false.
         twom24 = 1.
         do i= 1, 24
           next(i) = i-1
           twom24 = twom24 * 0.5
         end do
         next(1) = 24
         twom12 = twom24 * 4096.
      write(lout,"(a)")      "RANLUX> Full initialization of ranlux with 25 integers: "
      write(lout,"(a,5i12)") "RANLUX> ",isdext
      do 200 i= 1, 24
      seeds(i) = real(isdext(i))*twom24
  200 continue
      carry = 0.
      if (isdext(25) .lt. 0)  carry = twom24
      isd = iabs(isdext(25))
      i24 = mod(isd,100)
      isd = isd/100
      j24 = mod(isd,100)
      isd = isd/100
      in24 = mod(isd,100)
      isd = isd/100
      luxlev = isd
        if (luxlev .le. maxlev) then
          nskip = ndskip(luxlev)
          write(lout,'(A,I2)')' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ', luxlev
        else  if (luxlev .ge. 24) then
          nskip = luxlev - 24
          write(lout,'(A,I5)')' RANLUX P-VALUE SET BY RLUXIN TO:',luxlev
        else
          nskip = ndskip(maxlev)
          write(lout,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',luxlev
          luxlev = maxlev
        endif
      inseed = -1
      return
end subroutine rluxin

! Entry to ouput seeds as integers
subroutine rluxut(isdext_tmp)
  implicit none
  integer, intent(in) :: isdext_tmp(25)
  isdext = isdext_tmp
      do 300 i= 1, 24
         isdext(i) = int(seeds(i)*twop12*twop12)
  300 continue
      isdext(25) = i24 + 100*j24 + 10000*in24 + 1000000*luxlev
      if (carry .gt. 0.)  isdext(25) = -isdext(25)
      return
end subroutine rluxut

! Entry to output the "convenient" restart point
! Note: The first argument was originall called "lout";
! however this conflicts with the variable name used for selecting output unit.
! It was therefore renamed to "lout2".
subroutine rluxat(lout2,inout,k1,k2)
  implicit none
  integer lout2
  integer inout

  integer k1,k2

  lout2 = luxlev
  inout = inseed
  k1 = kount
  k2 = mkount
  return
end subroutine rluxat

! Entry to initialize from one or three integers
subroutine rluxgo(lux,ins,k1,k2)
  implicit none
  integer lux,ins,k1,k2
  if (lux .lt. 0) then
     luxlev = lxdflt
  else if (lux .le. maxlev) then
     luxlev = lux
  else if (lux .lt. 24 .or. lux .gt. 2000) then
     luxlev = maxlev
     write(lout,"(a,i7)") "RANLUX> Illegal luxury rluxgo: ",lux
  else
     luxlev = lux
     do 310 ilx= 0, maxlev
       if (lux .eq. ndskip(ilx)+24)  luxlev = ilx
310       continue
  endif

  if (luxlev .le. maxlev)  then
     nskip = ndskip(luxlev)
     write(lout,"(a,i2,a,i4)") "RANLUX> Luxury level set by rluxgo: ",luxlev," P = ",nskip+24
  else
      nskip = luxlev - 24
      write(lout,"(a,i5)") "RANLUX> P-value set by rluxgo to: ",luxlev
  endif

  in24 = 0

  if (ins .lt. 0)  write(lout,"(a)") "RANLUX> Illegal initialization by RLUXGO, negative input seed"
  if (ins .gt. 0)  then
    jseed = ins
    write(lout,"(a,3i12)") "RANLUX> Initialized by rluxgo from seeds ",jseed,k1,k2
  else
    jseed = jsdflt
    write(lout,"(a)") "RANLUX> Initialized by rluxgo from default seed"
  endif
  inseed = jseed
  notyet = .false.
  twom24 = 1.
  do 325 i= 1, 24
    twom24 = twom24 * 0.5
  k = jseed/53668
  jseed = 40014*(jseed-k*53668) -k*12211
  if (jseed .lt. 0)  jseed = jseed+icons
  iseeds(i) = mod(jseed,itwo24)
325    continue
    twom12 = twom24 * 4096.
   do 350 i= 1,24
   seeds(i) = real(iseeds(i))*twom24
   next(i) = i-1
350    continue
  next(1) = 24
  i24 = 24
  j24 = 10
  carry = 0.
  if (seeds(24) .eq. 0.) carry = twom24
!        If restarting at a break point, skip K1 + IGIGA*K2
!        Note that this is the number of numbers delivered to
!        the user PLUS the number skipped (if luxury .GT. 0).
    kount = k1
    mkount = k2
    if (k1+k2 .ne. 0)  then
      do 500 iouter= 1, k2+1
        inner = igiga
        if (iouter .eq. k2+1)  inner = k1
        do 450 isk= 1, inner
          uni = seeds(j24) - seeds(i24) - carry
          if (uni .lt. 0.)  then
             uni = uni + 1.
             carry = twom24
          else
             carry = 0.
          endif
          seeds(i24) = uni
          i24 = next(i24)
          j24 = next(j24)
450     continue
500   continue
!         Get the right value of IN24 by direct calculation
      in24 = mod(kount, nskip+24)
      if (mkount .gt. 0)  then
         izip = mod(igiga, nskip+24)
         izip2 = mkount*izip + in24
         in24 = mod(izip2, nskip+24)
      endif
!       Now IN24 had better be between zero and 23 inclusive
      if (in24 .gt. 23) then
         write(lout,"(a)")           "RANLUX> ERROR RESTARTING with RLUXGO:"
         write(lout,"(a,3i11,a,i5)") "RANLUX>       The values",ins,k1,k2," cannot occur at luxury level ",luxlev
         in24 = 0
      endif
    endif
    return
end subroutine rluxgo

function rndm4()

  implicit none

  integer len, in
  real(kind=fPrec) rndm4, a
  save IN,a
  parameter ( len =  30000 )
  dimension a(len)
  data in/1/

  if( in.eq.1 ) then
    call ranlux(a,len)
!    call ranecu(a,len,-1)
    rndm4=a(1)
    in=2
  else
    rndm4=a(in)
    in=in+1
    if(in.eq.len+1)in=1
  endif

  return

end function rndm4


!ccccccccccccccccccccccccccccccccccccccc
!-TW-01/2007
! function rndm5(irnd) , irnd = 1 will reset
! inn counter => enables reproducible set of
! random unmbers
!cccccccccccccccccccccccccccccccccc
!
function rndm5(irnd)

  use mathlib_bouncer

  implicit none

  integer len, inn, irnd
  real(kind=fPrec) rndm5,a
  save

  parameter( len =  30000 )
  dimension a(len)
  data inn/1/
!
! reset inn to 1 enable reproducible random numbers
  if( irnd .eq. 1) inn = 1

  if( inn.eq.1 ) then
    call ranlux(a,len)
!     call ranecu(a,len,-1)
    rndm5=a(1)
    inn=2
  else
    rndm5=a(inn)
    inn=inn+1
    if(inn.eq.len+1)inn=1
  end if

  return
end function rndm5

real(kind=fPrec) function ran_gauss(cut)
!*********************************************************************
!
! RAN_GAUSS - will generate a normal distribution from a uniform
!   distribution between [0,1].
!   See "Communications of the ACM", V. 15 (1972), p. 873.
!
! cut - real(kind=fPrec) - cut for distribution in units of sigma
!                the cut must be greater than 0.5
!
!*********************************************************************

  use crcoall
  use parpro
  use mathlib_bouncer
  implicit none

  logical flag
  DATA flag/.TRUE./
  real(kind=fPrec) x, u1, u2, twopi, r,cut

  save

  twopi=eight*atan_mb(one) !Why not 2*pi, where pi is in block "common"?
1 if (flag) then
    r = real(rndm4(),fPrec)
    r = max(r, half**32)
    r = min(r, one-half**32)
    u1 = sqrt(-two*log_mb( r ))
    u2 = real(rndm4(),fPrec)
    x = u1 * cos_mb(twopi*u2)
  else
    x = u1 * sin_mb(twopi*u2)
  endif

  flag = .not. flag

!  cut the distribution if cut > 0.5
  if (cut .gt. half .and. abs(x) .gt. cut) goto 1

  ran_gauss = x
  return
end function ran_gauss

end module mod_ranlux
