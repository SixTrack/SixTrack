! ================================================================================================ !
!  RANLUX
! ~~~~~~~~
!                                 default
!   Luxury Level    0     1     2   *3*    4
!    data ndskip    0    24    73   199   365
! Corresponds to p=24    48    97   223   389
!      time factor  1    2     3     6    10   on slow workstation
!                   1   1.5    2     3     5   on fast mainframe
!
! Reference:
! James, F. ‘RANLUX: A FORTRAN Implementation of the High Quality Pseudorandom Number Generator of
!     Luscher’. Comput.Phys.Commun. 79 (1994): 111–14. https://doi.org/10.1016/0010-4655(94)90233-X.
! ================================================================================================ !
module mod_ranlux

  use floatPrecision
  use numerical_constants

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

  real(kind=fPrec), parameter :: twop12 = 4096.0
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

  public ranlux
  public rluxin
  public rluxut
  public rluxat
  public rluxgo
  public rndm4
  public rndm5
  public coll_rand
  public ran_gauss
  public ran_gauss2

contains

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
subroutine ranlux(rvec,lenv)

  use crcoall

  integer,          intent(in)  :: lenv
  real(kind=fPrec), intent(out) :: rvec(lenv)

  ! NOTYET is .TRUE. if no initialization has been performed yet.
  !             Default Initialization by Multiplicative Congruential
  if(notyet) then
    notyet = .false.
    jseed = jsdflt
    inseed = jseed
    write(lout,"(a,i0)") "RANLUX> Default initialization: ",jseed
    luxlev = lxdflt
    nskip  = ndskip(luxlev)
    lp     = nskip + 24
    in24   = 0
    kount  = 0
    mkount = 0
    twom24 = 1.0
    do i=1,24
      twom24 = twom24 * 0.5
      k = jseed/53668
      jseed = 40014*(jseed-k*53668) -k*12211
      if(jseed < 0) jseed = jseed+icons
      iseeds(i) = mod(jseed,itwo24)
    end do
    twom12 = twom24 * twop12
    do i=1,24
      seeds(i) = real(iseeds(i))*twom24
      next(i) = i-1
    end do
    next(1) = 24
    i24     = 24
    j24     = 10
    carry   = 0.0
    if (seeds(24) == 0.0) carry = twom24
  end if

  ! The Generator proper: "Subtract-with-borrow",
  ! as proposed by Marsaglia and Zaman,
  ! Florida State University, March, 1989

  do ivec=1,lenv
    uni = seeds(j24) - seeds(i24) - carry
    if(uni < 0.0) then
      uni = uni + 1.0
      carry = twom24
    else
      carry = 0.0
    end if
    seeds(i24) = uni
    i24 = next(i24)
    j24 = next(j24)
    rvec(ivec) = uni

    ! small numbers (with less than 12 "significant" bits) are "padded".
    if(uni < twom12) then
      rvec(ivec) = rvec(ivec) + twom24*seeds(j24)
      ! and zero is forbidden in case someone takes a logarithm
      if(rvec(ivec) == 0.0) rvec(ivec) = twom24*twom24
    end if

    ! Skipping to luxury.  As proposed by Martin Luscher.
    in24 = in24 + 1
    if(in24 == 24) then
      in24 = 0
      kount = kount + nskip
      do isk=1,nskip
        uni = seeds(j24) - seeds(i24) - carry
        if(uni < 0.0) then
          uni = uni + 1.0
          carry = twom24
        else
          carry = 0.0
        end if
        seeds(i24) = uni
        i24 = next(i24)
        j24 = next(j24)
      end do
    end if
  end do

  kount = kount + lenv
  if(kount >= igiga) then
    mkount = mkount + 1
    kount = kount - igiga
  end if

end subroutine ranlux

! Entry to input and float integer seeds from previous run
subroutine rluxin(isdext_tmp)

  use crcoall

  integer, intent(in) :: isdext_tmp(25)

  isdext = isdext_tmp
  notyet = .false.
  twom24 = 1.0
  do i=1,24
    next(i) = i-1
    twom24 = twom24 * 0.5
  end do
  next(1) = 24
  twom12 = twom24 * twop12
  write(lout,"(a)")       "RANLUX> Full initialisation of ranlux with 25 integers: "
! write(lout,"(a,25i12)") "RANLUX> ",isdext
  do i=1,24
    seeds(i) = real(isdext(i))*twom24
  end do

  carry = 0.0
  if(isdext(25) < 0) carry = twom24
  isd = iabs(isdext(25))
  i24 = mod(isd,100)
  isd = isd/100
  j24 = mod(isd,100)
  isd = isd/100
  in24 = mod(isd,100)
  isd = isd/100
  luxlev = isd
  if(luxlev <= maxlev) then
    nskip = ndskip(luxlev)
    write(lout,"(a,i0)") "RANLUX> Luxury level set by rluxin to ",luxlev
  else if(luxlev >= 24) then
    nskip = luxlev - 24
    write(lout,"(a,i0)") "RANLUX> P-value set by rluxin to ",luxlev
  else
    nskip = ndskip(maxlev)
    write(lout,"(a,i0)") "RANLUX> Illegal luxury rluxin ",luxlev
    luxlev = maxlev
  end if
  inseed = -1

end subroutine rluxin

! Entry to ouput seeds as integers
subroutine rluxut(isdext_tmp)

  integer, intent(out) :: isdext_tmp(25)

  do i=1,24
    isdext_tmp(i) = int(seeds(i)*twop12*twop12)
  end do
  isdext_tmp(25) = i24 + 100*j24 + 10000*in24 + 1000000*luxlev
  if(carry > 0.0) isdext_tmp(25) = -isdext_tmp(25)

end subroutine rluxut

! Entry to output the "convenient" restart point
subroutine rluxat(lout2,inout,k1,k2)

  integer, intent(out) :: lout2,inout,k1,k2

  lout2 = luxlev
  inout = inseed
  k1    = kount
  k2    = mkount

end subroutine rluxat

! Entry to initialize from one or three integers
subroutine rluxgo(lux,ins,k1,k2)

  use crcoall

  integer, intent(in) :: lux,ins,k1,k2

  if(lux < 0) then
    luxlev = lxdflt
  elseif(lux <= maxlev) then
    luxlev = lux
  elseif(lux < 24 .or. lux > 2000) then
    luxlev = maxlev
    write(lout,"(a,i0)") "RANLUX> Illegal luxury rluxgo ",lux
  else
    luxlev = lux
    do ilx=0,maxlev
      if(lux == ndskip(ilx)+24) luxlev = ilx
    end do
  end if

  if(luxlev <= maxlev) then
    nskip = ndskip(luxlev)
    write(lout,"(a,i0,a,i0)") "RANLUX> Luxury level set by rluxgo: ",luxlev," P = ",nskip+24
  else
    nskip = luxlev - 24
    write(lout,"(a,i0)") "RANLUX> P-value set by rluxgo to ",luxlev
  end if

  in24 = 0

  if(ins < 0) write(lout,"(a)") "RANLUX> Illegal initialisation by RLUXGO, negative input seed"
  if(ins > 0) then
    jseed = ins
    write(lout,"(a,3(1x,i0))") "RANLUX> Initialised by rluxgo from seeds ",jseed,k1,k2
  else
    jseed = jsdflt
    write(lout,"(a)") "RANLUX> Initialised by rluxgo from default seed"
  end if

  inseed = jseed
  notyet = .false.
  twom24 = 1.0
  do i=1,24
    twom24 = twom24 * 0.5
    k = jseed/53668
    jseed = 40014*(jseed-k*53668) - k*12211
    if(jseed < 0) jseed = jseed+icons
    iseeds(i) = mod(jseed,itwo24)
  end do

  twom12 = twom24 * twop12
  do i=1,24
    seeds(i) = real(iseeds(i))*twom24
    next(i) = i-1
  end do

  next(1) = 24
  i24     = 24
  j24     = 10
  carry   = 0.0
  if(seeds(24) == 0.0) carry = twom24

  ! If restarting at a break point, skip K1 + IGIGA*K2
  ! Note that this is the number of numbers delivered to
  ! the user PLUS the number skipped (if luxury .GT. 0).
  kount = k1
  mkount = k2
  if(k1+k2 /= 0) then
    do iouter=1,k2+1
      inner = igiga
      if(iouter == k2+1) inner = k1
      do isk=1,inner
        uni = seeds(j24) - seeds(i24) - carry
        if(uni < 0.) then
          uni = uni + 1.0
          carry = twom24
        else
          carry = 0.0
        end if
        seeds(i24) = uni
        i24 = next(i24)
        j24 = next(j24)
      end do
    end do

    ! Get the right value of IN24 by direct calculation
    in24 = mod(kount, nskip+24)
    if(mkount > 0) then
      izip = mod(igiga, nskip+24)
      izip2 = mkount*izip + in24
      in24 = mod(izip2, nskip+24)
    end if

    ! Now IN24 had better be between zero and 23 inclusive
    if(in24 > 23) then
      write(lerr,"(a)")               "RANLUX> ERROR RESTARTING with RLUXGO:"
      write(lerr,"(a,3(1x,i0),a,i0)") "RANLUX>       The values",ins,k1,k2," cannot occur at luxury level ",luxlev
      in24 = 0
    end if
  end if

end subroutine rluxgo

function rndm4()

  integer len, in
  real(kind=fPrec) rndm4, a

  save in,a

  parameter ( len =  30000 )
  dimension a(len)
  data in/1/

  if(in == 1) then
    call ranlux(a,len)
!    call ranecu(a,len,-1)
    rndm4=a(1)
    in=2
  else
    rndm4=a(in)
    in=in+1
    if(in == len+1)in=1
  endif

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
  if(irnd == 1) inn = 1

  if(inn == 1) then
    call ranlux(a,len)
!     call ranecu(a,len,-1)
    rndm5=a(1)
    inn=2
  else
    rndm5=a(inn)
    inn=inn+1
    if(inn == len+1) inn=1
  end if

end function rndm5

!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! function coll_rand()
! to be used for the collimation code with CR enabled.
! we get a fresh random number every call, so no need to have a 30k buffer here.
! Therefore this is the same as rndm5 with an argument of 1 every time
!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function coll_rand()

  use mathlib_bouncer

  implicit none

  integer len
  real(kind=fPrec) coll_rand,a
  save

  parameter( len =  1 )
  dimension a(len)

  call ranlux(a,len)
  coll_rand=a(1)

end function coll_rand


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
real(kind=fPrec) function ran_gauss(cut)

  use mathlib_bouncer

  real(kind=fPrec), intent(in) :: cut

  logical :: flag = .true.
  real(kind=fPrec) :: x, u1, u2, r

  save

1 if (flag) then
    r  = rndm4()
    r  = max(r, half**32)
    r  = min(r, one-half**32)
    u1 = sqrt(-two*log_mb( r ))
    u2 = rndm4()
    x  = u1 * cos_mb(twopi*u2)
  else
    x  = u1 * sin_mb(twopi*u2)
  end if

  flag = .not. flag

  ! Cut the distribution if cut > 0.5
  if(cut > half .and. abs(x) > cut) goto 1

  ran_gauss = x

end function ran_gauss

!*********************************************************************
!
! ran_gauss2 - will generate a normal distribution from a uniform
!     distribution between [0,1].
!     See "Communications of the ACM", V. 15 (1972), p. 873.
!
!     cut - real(kind=fPrec) - cut for distribution in units of sigma
!     the cut must be greater than 0.5
!
!     changed rndm4 to rndm5(irnd) and defined flag as true
!
!*********************************************************************
real(kind=fPrec) function ran_gauss2(cut)

  use numerical_constants, only : twopi
  use mathlib_bouncer

  logical flag
  real(kind=fPrec) x, u1, u2, r,cut
  save

  flag = .true. !Does this initialize only once, or is it executed every pass?
                !See ran_gauss(cut)

1 if (flag) then
    r = real(rndm5(0),fPrec)
    r = max(r, half**32)
    r = min(r, one-half**32)
    u1 = sqrt(-two*log_mb( r ))
    u2 = real(rndm5(0),fPrec)
    x = u1 * cos_mb(twopi*u2)
  else
    x = u1 * sin_mb(twopi*u2)
  endif

  flag = .not. flag

  ! cut the distribution if cut > 0.5
  if(cut > half .and. abs(x) > cut) goto 1

  ran_gauss2 = x

end function ran_gauss2

end module mod_ranlux
