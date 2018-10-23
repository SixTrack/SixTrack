! =================================================================================================
!  Mathlib Bouncer
!  Last modified: 2018-04-27
! =================================================================================================
module mathlib_bouncer

  use floatPrecision
! use, intrinsic :: ieee_arithmetic
  use, intrinsic :: iso_fortran_env, only : real64, int64
  implicit none

  !These are the actual "bouncy functions" users are supposed to call
  public :: sin_mb, asin_mb, sinh_mb, cos_mb, acos_mb, cosh_mb, tan_mb, atan_mb, atan2_mb, exp_mb, log_mb, log10_mb, isnan_mb

! real(kind=fPrec), parameter :: mb_pi   = 3.1415926535897932d0
! real(kind=fPrec), parameter :: mb_pi2  = 1.5707963267948966d0
  real(kind=fPrec), parameter :: mb_pi   = 3.141592653589793238462643383279502884197169399375105820974_fPrec
  real(kind=fPrec), parameter :: mb_pi2  = 1.570796326794896619231321691639751442098584699687552910487_fPrec
! #ifdef NAGFOR
!   real(kind=fPrec), parameter :: mb_qnan = ieee_value(1.0_fPrec, ieee_quiet_nan)
! #else
!   real(kind=fPrec), parameter :: mb_qnan = transfer(-2251799813685248_int64, 1.0_real64)
! #endif

  !For linking with CRLIBM
#ifdef CRLIBM

#ifdef ROUND_NEAR
  private :: acos_rn, asin_rn, atan2_rn
#endif

#ifdef ROUND_UP
  private :: acos_ru, asin_ru, atan2_ru
#endif

#ifdef ROUND_DOWN
  private :: acos_rd, asin_rd, atan2_rd
#endif

#ifdef ROUND_ZERO
  private :: acos_rz, asin_rz, atan2_rz
#endif

  !Can't declare them private as they have BIND labels,
  ! however they should not be called from outside this module.
  !private :: exp_rn, log_rn, log10_rn, atan_rn, tan_rn, sin_rn, cos_rn, sinh_rn, cosh_rn

  !! Interface definitions for the other functions in crlibm
  interface

#ifdef ROUND_NEAR
     real(kind=c_double) function exp_rn(arg) bind(C,name="exp_rn")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function exp_rn

     real(kind=c_double) function log_rn(arg) bind(C,name="log_rn")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function log_rn

     real(kind=c_double) function log10_rn(arg) bind(C,name="log10_rn")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function log10_rn

     real(kind=c_double) function atan_rn(arg) bind(C,name="atan_rn")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function atan_rn

     real(kind=c_double) function tan_rn(arg) bind(C,name="tan_rn")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function tan_rn

     real(kind=c_double) function sin_rn(arg) bind(C,name="sin_rn")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function sin_rn

     real(kind=c_double) function cos_rn(arg) bind(C,name="cos_rn")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function cos_rn

     real(kind=c_double) function sinh_rn(arg) bind(C,name="sinh_rn")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function sinh_rn

     real(kind=c_double) function cosh_rn(arg) bind(C,name="cosh_rn")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function cosh_rn
#endif

#ifdef ROUND_UP
     real(kind=c_double) function exp_ru(arg) bind(C,name="exp_ru")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function exp_ru

     real(kind=c_double) function log_ru(arg) bind(C,name="log_ru")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function log_ru

     real(kind=c_double) function log10_ru(arg) bind(C,name="log10_ru")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function log10_ru

     real(kind=c_double) function atan_ru(arg) bind(C,name="atan_ru")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function atan_ru

     real(kind=c_double) function tan_ru(arg) bind(C,name="tan_ru")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function tan_ru

     real(kind=c_double) function sin_ru(arg) bind(C,name="sin_ru")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function sin_ru

     real(kind=c_double) function cos_ru(arg) bind(C,name="cos_ru")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function cos_ru

     real(kind=c_double) function sinh_ru(arg) bind(C,name="sinh_ru")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function sinh_ru

     real(kind=c_double) function cosh_ru(arg) bind(C,name="cosh_ru")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function cosh_ru
#endif

#ifdef ROUND_DOWN
     real(kind=c_double) function exp_rd(arg) bind(C,name="exp_rd")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function exp_rd

     real(kind=c_double) function log_rd(arg) bind(C,name="log_rd")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function log_rd

     real(kind=c_double) function log10_rd(arg) bind(C,name="log10_rd")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function log10_rd

     real(kind=c_double) function atan_rd(arg) bind(C,name="atan_rd")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function atan_rd

     real(kind=c_double) function tan_rd(arg) bind(C,name="tan_rd")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function tan_rd

     real(kind=c_double) function sin_rd(arg) bind(C,name="sin_rd")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function sin_rd

     real(kind=c_double) function cos_rd(arg) bind(C,name="cos_rd")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function cos_rd

     real(kind=c_double) function sinh_rd(arg) bind(C,name="sinh_rd")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function sinh_rd

     real(kind=c_double) function cosh_rd(arg) bind(C,name="cosh_rd")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function cosh_rd
#endif

#ifdef ROUND_ZERO
     real(kind=c_double) function exp_rz(arg) bind(C,name="exp_rd")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function exp_rz

     real(kind=c_double) function log_rz(arg) bind(C,name="log_rz")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function log_rz

     real(kind=c_double) function log10_rz(arg) bind(C,name="log10_rz")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function log10_rz

     real(kind=c_double) function atan_rz(arg) bind(C,name="atan_rz")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function atan_rz

     real(kind=c_double) function tan_rz(arg) bind(C,name="tan_rz")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function tan_rz

     real(kind=c_double) function sin_rz(arg) bind(C,name="sin_rz")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function sin_rz

     real(kind=c_double) function cos_rz(arg) bind(C,name="cos_rz")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function cos_rz

     real(kind=c_double) function sinh_rz(arg) bind(C,name="sinh_rz")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function sinh_rz

     real(kind=c_double) function cosh_rz(arg) bind(C,name="cosh_rz")
       use, intrinsic :: iso_c_binding, only : c_double
       implicit none
       real(kind=c_double), intent(in), VALUE :: arg
     end function cosh_rz
#endif

  end interface
#endif

contains

  ! Definition of the MathlibBouncer (_mb) functions

! #ifdef IFORT
!   logical pure elemental function isnan_mb(arg)
!     use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
!     real(kind=fPrec), intent(in) :: arg
!     isnan_mb = ieee_is_nan(arg)
!   end function isnan_mb
! #else
  logical pure elemental function isnan_mb(arg)
    real(kind=fPrec), intent(in) :: arg
    isnan_mb = arg /= arg
  end function isnan_mb
! #endif

  real(kind=fPrec) function sin_mb(arg)
    implicit none
    real(kind=fPrec) arg
    intent(in) arg

#ifndef CRLIBM
    !Input KIND = output KIND
    sin_mb=sin(arg)
#else
#ifdef ROUND_NEAR
    sin_mb=sin_rn(arg)
#endif
#ifdef ROUND_UP
    sin_mb=sin_ru(arg)
#endif
#ifdef ROUND_DOWN
    sin_mb=sin_rd(arg)
#endif
#ifdef ROUND_ZERO
    sin_mb=sin_rz(arg)
#endif
#endif
  end function sin_mb

  real(kind=fPrec) function asin_mb(arg)
    implicit none
    real(kind=fPrec) arg
    intent(in) arg

#ifndef CRLIBM
    !Input KIND = output KIND
    asin_mb=asin(arg)
#else
#ifdef ROUND_NEAR
    asin_mb=asin_rn(arg)
#endif
#ifdef ROUND_UP
    asin_mb=asin_ru(arg)
#endif
#ifdef ROUND_DOWN
    asin_mb=asin_rd(arg)
#endif
#ifdef ROUND_ZERO
    asin_mb=asin_rz(arg)
#endif
#endif
  end function asin_mb

  real(kind=fPrec) function sinh_mb(arg)
    implicit none
    real(kind=fPrec) arg
    intent(in) arg

#ifndef CRLIBM
    !Input KIND = output KIND
    sinh_mb=sinh(arg)
#else
#ifdef ROUND_NEAR
    sinh_mb=sinh_rn(arg)
#endif
#ifdef ROUND_UP
    sinh_mb=sinh_ru(arg)
#endif
#ifdef ROUND_DOWN
    sinh_mb=sinh_rd(arg)
#endif
#ifdef ROUND_ZERO
    sinh_mb=sinh_rz(arg)
#endif
#endif
  end function sinh_mb

  real(kind=fPrec) function cos_mb(arg)
    implicit none
    real(kind=fPrec) arg
    intent(in) arg

#ifndef CRLIBM
    !Input KIND = output KIND
    cos_mb=cos(arg)
#else
#ifdef ROUND_NEAR
    cos_mb=cos_rn(arg)
#endif
#ifdef ROUND_UP
    cos_mb=cos_ru(arg)
#endif
#ifdef ROUND_DOWN
    cos_mb=cos_rd(arg)
#endif
#ifdef ROUND_ZERO
    cos_mb=cos_rz(arg)
#endif
#endif
  end function cos_mb

  real(kind=fPrec) function acos_mb(arg)
    implicit none
    real(kind=fPrec) arg
    intent(in) arg

#ifndef CRLIBM
    !Input KIND = output KIND
    acos_mb=acos(arg)
#else
#ifdef ROUND_NEAR
    acos_mb=acos_rn(arg)
#endif
#ifdef ROUND_UP
    acos_mb=acos_ru(arg)
#endif
#ifdef ROUND_DOWN
    acos_mb=acos_rd(arg)
#endif
#ifdef ROUND_ZERO
    acos_mb=acos_rz(arg)
#endif
#endif
  end function acos_mb

  real(kind=fPrec) function cosh_mb(arg)
    implicit none
    real(kind=fPrec) arg
    intent(in) arg

#ifndef CRLIBM
    !Input KIND = output KIND
    cosh_mb=cosh(arg)
#else
#ifdef ROUND_NEAR
    cosh_mb=cosh_rn(arg)
#endif
#ifdef ROUND_UP
    cosh_mb=cosh_ru(arg)
#endif
#ifdef ROUND_DOWN
    cosh_mb=cosh_rd(arg)
#endif
#ifdef ROUND_ZERO
    cosh_mb=cosh_rz(arg)
#endif
#endif
  end function cosh_mb

  real(kind=fPrec) function tan_mb(arg)
    implicit none
    real(kind=fPrec) arg
    intent(in) arg

#ifndef CRLIBM
    !Input KIND = output KIND
    tan_mb=tan(arg)
#else
#ifdef ROUND_NEAR
    tan_mb=tan_rn(arg)
#endif
#ifdef ROUND_UP
    tan_mb=tan_ru(arg)
#endif
#ifdef ROUND_DOWN
    tan_mb=tan_rd(arg)
#endif
#ifdef ROUND_ZERO
    tan_mb=tan_rz(arg)
#endif
#endif
  end function tan_mb

  real(kind=fPrec) function atan_mb(arg)
    implicit none
    real(kind=fPrec) arg
    intent(in) arg

#ifndef CRLIBM
    !Input KIND = output KIND
    atan_mb=atan(arg)
#else
#ifdef ROUND_NEAR
    atan_mb=atan_rn(arg)
#endif
#ifdef ROUND_UP
    atan_mb=atan_ru(arg)
#endif
#ifdef ROUND_DOWN
    atan_mb=atan_rd(arg)
#endif
#ifdef ROUND_ZERO
    atan_mb=atan_rz(arg)
#endif
#endif
  end function atan_mb

  real(kind=fPrec) function atan2_mb(y,x)
    implicit none
    real(kind=fPrec), intent(in) :: y,x
#ifndef CRLIBM
#ifdef NAGFOR
    ! Nagfor
    if(x == 0.0_fPrec .and. y == 0.0_fPrec) then
      atan2_mb=0.0_fPrec
    else
      atan2_mb=atan2(y,x)
    end if
#else
    atan2_mb=atan2(y,x)
#endif
#else
#ifdef ROUND_NEAR
    atan2_mb=atan2_rn(y,x)
#endif
#ifdef ROUND_UP
    atan2_mb=atan2_ru(y,x)
#endif
#ifdef ROUND_DOWN
    atan2_mb=atan2_rd(y,x)
#endif
#ifdef ROUND_ZERO
    atan2_mb=atan2_rz(y,x)
#endif
#endif
  end function atan2_mb

  real(kind=fPrec) function exp_mb(arg)
    implicit none
    real(kind=fPrec) arg
    intent(in) arg

#ifndef CRLIBM
    !Input KIND = output KIND
    exp_mb=exp(arg)
#else
#ifdef ROUND_NEAR
    exp_mb=exp_rn(arg)
#endif
#ifdef ROUND_UP
    exp_mb=exp_ru(arg)
#endif
#ifdef ROUND_DOWN
    exp_mb=exp_rd(arg)
#endif
#ifdef ROUND_ZERO
    exp_mb=exp_rz(arg)
#endif
#endif
  end function exp_mb

  real(kind=fPrec) function log_mb(arg)
    implicit none
    real(kind=fPrec) arg
    intent(in) arg

#ifndef CRLIBM
    !Input KIND = output KIND
    log_mb=log(arg)
#else
#ifdef ROUND_NEAR
    log_mb=log_rn(arg)
#endif
#ifdef ROUND_UP
    log_mb=log_ru(arg)
#endif
#ifdef ROUND_DOWN
    log_mb=log_rd(arg)
#endif
#ifdef ROUND_ZERO
    log_mb=log_rz(arg)
#endif
#endif
  end function log_mb

  real(kind=fPrec) function log10_mb(arg)
    implicit none
    real(kind=fPrec) arg
    intent(in) arg

#ifndef CRLIBM
    !Input KIND = output KIND
    log10_mb=log10(arg)
#else
#ifdef ROUND_NEAR
    log10_mb=log10_rn(arg)
#endif
#ifdef ROUND_UP
    log10_mb=log10_ru(arg)
#endif
#ifdef ROUND_DOWN
    log10_mb=log10_rd(arg)
#endif
#ifdef ROUND_ZERO
    log10_mb=log10_rz(arg)
#endif
#endif
  end function log10_mb

#ifdef CRLIBM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Functions moved from sixtrack.s, wrapping parts of crlibm !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef ROUND_NEAR
  real(kind=real64) function acos_rn(x)
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    real(kind=real64) x
    if(x /= x) then ! Check if NaN
       acos_rn=x
    elseif (abs(x).eq.0.0d0) then
       acos_rn=mb_pi2
    else
       !       acos_rn=atan_rn(sqrt(1.0d0-x*x)/x)
       ! Try using (1-x)*(1+x) in case x is very small.........
       ! or close to 1.....write a test program!!!
       acos_rn=atan_rn(sqrt((1.0d0-x)*(1.0d0+x))/x)
       if (x.lt.0d0) then
          acos_rn=mb_pi+acos_rn
       endif
    endif
  end function acos_rn

  real(kind=real64) function asin_rn(x)
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    real(kind=real64) x,mb_pi2
    if(x /= x) then ! Check if NaN
       asin_rn=x
       return
    endif
    if (abs(x).eq.1.0d0) then
       asin_rn=sign(mb_pi2,x)
    else
       !       asin_rn=atan_rn(x/sqrt(1.0d0-x*x))
       ! Try using (1-x)*(1+x) in case x is very small.........
       ! or close to 1.....write a test program!!!
       asin_rn=atan_rn(x/sqrt((1.0d0-x)*(1.0d0+x)))
    endif
  end function asin_rn

  real(kind=real64) function atan2_rn(y,x)
    implicit none
    real(kind=real64) x,y
    if (x.eq.0d0) then
      if (y.eq.0d0) then
#ifdef NAGFOR
        atan2_rn=0d0
#else
        ! Let the internal atan2 handle this case according to ieee
        atan2_rn=atan2(y,x)
#endif
      else
        atan2_rn=sign(mb_pi2,y)
      endif
    else
       if (y.eq.0d0) then
        if (x.gt.0d0) then
          ! Let the internal atan2 handle this case according to ieee
          atan2_rn=atan2(y,x)
        else
          atan2_rn=mb_pi
        endif
      else
        atan2_rn=atan_rn(y/x)
        if (x.lt.0d0) then
          atan2_rn=sign(mb_pi,y)+atan2_rn
        endif
       endif
    endif
  end function atan2_rn
#endif

#ifdef ROUND_UP
  real(kind=real64) function acos_ru(x)
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    real(kind=real64) x
    if(x /= x) then ! Check if NaN
       acos_ru=x
    elseif (abs(x).eq.0.0d0) then
       acos_ru=mb_pi2
    else
       !       acos_ru=atan_ru(sqrt(1.0d0-x*x)/x)
       ! Try using (1-x)*(1+x) in case x is very small.........
       ! or close to 1.....write a test program!!!
       acos_ru=atan_ru(sqrt((1.0d0-x)*(1.0d0+x))/x)
       if (x.lt.0d0) then
          acos_ru=mb_pi+acos_ru
       endif
    endif
  end function acos_ru

  real(kind=real64) function asin_ru(x)
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    real(kind=real64) x
    if(x /= x) then ! Check if NaN
       asin_ru=x
       return
    endif
    if (abs(x).eq.1.0d0) then
       asin_ru=sign(mb_pi2,x)
    else
       !       asin_ru=atan_ru(x/sqrt(1.0d0-x*x))
       ! Try using (1-x)*(1+x) in case x is very small.........
       ! or close to 1.....write a test program!!!
       asin_ru=atan_ru(x/sqrt((1.0d0-x)*(1.0d0+x)))
    endif
  end function asin_ru

  real(kind=real64) function atan2_ru(y,x)
    implicit none
    real(kind=real64) x,y
    if (x.eq.0d0) then
      if (y.eq.0d0) then
#ifdef NAGFOR
          atan2_ru=0d0
#else
        ! Let the internal atan2 handle this case according to ieee
        atan2_ru=atan2(y,x)
#endif
      else
        atan2_ru=sign(mb_pi2,y)
      endif
    else
      if (y.eq.0d0) then
        if (x.gt.0d0) then
          ! Let the internal atan2 handle this case according to ieee
          atan2_ru=atan2(y,x)
        else
          atan2_ru=mb_pi
        endif
      else
        atan2_ru=atan_ru(y/x)
        if (x.lt.0d0) then
          atan2_ru=sign(mb_pi,y)+atan2_ru
        endif
      endif
    endif
  end function atan2_ru
#endif

#ifdef ROUND_DOWN
  real(kind=real64) function acos_rd(x)
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    real(kind=real64) x
    if(x /= x) then ! Check if NaN
       acos_rd=x
    elseif (abs(x).eq.0.0d0) then
       acos_rd=mb_pi2
    else
       !       acos_rd=atan_rd(sqrt(1.0d0-x*x)/x)
       ! Try using (1-x)*(1+x) in case x is very small.........
       ! or close to 1.....write a test program!!!
       acos_rd=atan_rd(sqrt((1.0d0-x)*(1.0d0+x))/x)
       if (x.lt.0d0) then
          acos_rd=mb_pi+acos_rd
       endif
    endif
  end function acos_rd

  real(kind=real64) function asin_rd(x)
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    real(kind=real64) x
    if(x /= x) then ! Check if NaN
       asin_rd=x
       return
    endif
    if (abs(x).eq.1.0d0) then
       asin_rd=sign(mb_pi2,x)
    else
       !       asin_rd=atan_rd(x/sqrt(1.0d0-x*x))
       ! Try using (1-x)*(1+x) in case x is very small.........
       ! or close to 1.....write a test program!!!
       asin_rd=atan_rd(x/sqrt((1.0d0-x)*(1.0d0+x)))
    endif
  end function asin_rd

  real(kind=real64) function atan2_rd(y,x)
    implicit none
    real(kind=real64) x,y
    if (x.eq.0d0) then
      if (y.eq.0d0) then
#ifdef NAGFOR
        atan2_rd=0d0
#else
        ! Let the internal atan2 handle this case according to ieee
        atan2_rd=atan2(y,x)
#endif
      else
        atan2_rd=sign(mb_pi2,y)
      endif
    else
      if (y.eq.0d0) then
        if (x.gt.0d0) then
          ! Let the internal atan2 handle this case according to ieee
          atan2_rd=atan2(y,x)
        else
          atan2_rd=mb_pi
        endif
      else
        atan2_rd=atan_rd(y/x)
        if (x.lt.0d0) then
          atan2_rd=sign(mb_pi,y)+atan2_rd
        endif
      endif
    endif
  end function atan2_rd
#endif

#ifdef ROUND_ZERO
  real(kind=real64) function acos_rz(x)
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    real(kind=real64) x
    if(x /= x) then ! Check if NaN
       acos_rz=x
    elseif (abs(x).eq.0.0d0) then
       acos_rz=mb_pi2
    else
       !       acos_rz=atan_rz(sqrt(1.0d0-x*x)/x)
       ! Try using (1-x)*(1+x) in case x is very small.........
       ! or close to 1.....write a test program!!!
       acos_rz=atan_rz(sqrt((1.0d0-x)*(1.0d0+x))/x)
       if (x.lt.0d0) then
          acos_rz=mb_pi+acos_rz
       endif
    endif
  end function acos_rz

  real(kind=real64) function asin_rz(x)
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    real(kind=real64) x
    if(x /= x) then ! Check if NaN
       asin_rz=x
       return
    endif
    if (abs(x).eq.1.0d0) then
       asin_rz=sign(mb_pi2,x)
    else
       !       asin_rz=atan_rz(x/sqrt(1.0d0-x*x))
       ! Try using (1-x)*(1+x) in case x is very small.........
       ! or close to 1.....write a test program!!!
       asin_rz=atan_rz(x/sqrt((1.0d0-x)*(1.0d0+x)))
    endif
  end function asin_rz

  real(kind=real64) function atan2_rz(y,x)
    implicit none
    real(kind=real64) x,y
    if (x.eq.0d0) then
      if (y.eq.0d0) then
#ifdef NAGFOR
        atan2_rz=0d0
#else
        ! Let the internal atan2 handle this case according to ieee
        atan2_rz=atan2(y,x)
#endif
      else
        atan2_rz=sign(mb_pi2,y)
      endif
    else
      if (y.eq.0d0) then
        if (x.gt.0d0) then
          ! Let the internal atan2 handle this case according to ieee
          atan2_rz=atan2(y,x)
        else
          atan2_rz=mb_pi
        endif
      else
        atan2_rz=atan_rz(y/x)
        if (x.lt.0d0) then
          atan2_rz=sign(mb_pi,y)+atan2_rz
        endif
      endif
    endif
  end function atan2_rz
#endif

#endif

end module mathlib_bouncer
