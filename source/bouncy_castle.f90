! =================================================================================================
!  Mathlib Bouncer
!  Last modified: 2018-10-23
! =================================================================================================
module mathlib_bouncer

  use floatPrecision

  implicit none
  private

  public :: sin_mb, asin_mb, sinh_mb, cos_mb, acos_mb, cosh_mb
  public :: tan_mb, atan_mb, atan2_mb, exp_mb, log_mb, log10_mb, isnan_mb

  ! For linking with CRLIBM
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

  ! Interface definitions for the other functions in CRLIBM
  interface

#ifdef ROUND_NEAR
    real(kind=c_double) pure function exp_rn(arg) bind(C,name="exp_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function exp_rn

    real(kind=c_double) pure function log_rn(arg) bind(C,name="log_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function log_rn

    real(kind=c_double) pure function log10_rn(arg) bind(C,name="log10_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function log10_rn

    real(kind=c_double) pure function atan_rn(arg) bind(C,name="atan_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function atan_rn

    real(kind=c_double) pure function tan_rn(arg) bind(C,name="tan_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function tan_rn

    real(kind=c_double) pure function sin_rn(arg) bind(C,name="sin_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function sin_rn

    real(kind=c_double) pure function cos_rn(arg) bind(C,name="cos_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function cos_rn

    real(kind=c_double) pure function sinh_rn(arg) bind(C,name="sinh_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function sinh_rn

    real(kind=c_double) pure function cosh_rn(arg) bind(C,name="cosh_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function cosh_rn
#endif

#ifdef ROUND_UP
    real(kind=c_double) pure function exp_ru(arg) bind(C,name="exp_ru")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function exp_ru

    real(kind=c_double) pure function log_ru(arg) bind(C,name="log_ru")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function log_ru

    real(kind=c_double) pure function log10_ru(arg) bind(C,name="log10_ru")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function log10_ru

    real(kind=c_double) pure function atan_ru(arg) bind(C,name="atan_ru")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function atan_ru

    real(kind=c_double) pure function tan_ru(arg) bind(C,name="tan_ru")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function tan_ru

    real(kind=c_double) pure function sin_ru(arg) bind(C,name="sin_ru")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function sin_ru

    real(kind=c_double) pure function cos_ru(arg) bind(C,name="cos_ru")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function cos_ru

    real(kind=c_double) pure function sinh_ru(arg) bind(C,name="sinh_ru")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function sinh_ru

    real(kind=c_double) pure function cosh_ru(arg) bind(C,name="cosh_ru")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function cosh_ru
#endif

#ifdef ROUND_DOWN
    real(kind=c_double) pure function exp_rd(arg) bind(C,name="exp_rd")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function exp_rd

    real(kind=c_double) pure function log_rd(arg) bind(C,name="log_rd")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function log_rd

    real(kind=c_double) pure function log10_rd(arg) bind(C,name="log10_rd")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function log10_rd

    real(kind=c_double) pure function atan_rd(arg) bind(C,name="atan_rd")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function atan_rd

    real(kind=c_double) pure function tan_rd(arg) bind(C,name="tan_rd")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function tan_rd

    real(kind=c_double) pure function sin_rd(arg) bind(C,name="sin_rd")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function sin_rd

    real(kind=c_double) pure function cos_rd(arg) bind(C,name="cos_rd")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function cos_rd

    real(kind=c_double) pure function sinh_rd(arg) bind(C,name="sinh_rd")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function sinh_rd

    real(kind=c_double) pure function cosh_rd(arg) bind(C,name="cosh_rd")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function cosh_rd
#endif

#ifdef ROUND_ZERO
    real(kind=c_double) pure function exp_rz(arg) bind(C,name="exp_rd")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function exp_rz

    real(kind=c_double) pure function log_rz(arg) bind(C,name="log_rz")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function log_rz

    real(kind=c_double) pure function log10_rz(arg) bind(C,name="log10_rz")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function log10_rz

    real(kind=c_double) pure function atan_rz(arg) bind(C,name="atan_rz")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function atan_rz

    real(kind=c_double) pure function tan_rz(arg) bind(C,name="tan_rz")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function tan_rz

    real(kind=c_double) pure function sin_rz(arg) bind(C,name="sin_rz")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function sin_rz

    real(kind=c_double) pure function cos_rz(arg) bind(C,name="cos_rz")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function cos_rz

    real(kind=c_double) pure function sinh_rz(arg) bind(C,name="sinh_rz")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function sinh_rz

    real(kind=c_double) pure function cosh_rz(arg) bind(C,name="cosh_rz")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function cosh_rz
#endif

  end interface
#endif

contains

  ! ========================================================================== !
  !  Definition of the MathlibBouncer (_mb) callable functions
  ! ========================================================================== !

  logical pure elemental function isnan_mb(arg)
    real(kind=fPrec), intent(in) :: arg
    isnan_mb = arg /= arg
  end function isnan_mb

  real(kind=fPrec) pure function sin_mb(arg)
    real(kind=fPrec), intent(in) :: arg
#ifndef CRLIBM
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

  real(kind=fPrec) pure function asin_mb(arg)
    real(kind=fPrec), intent(in) :: arg
#ifndef CRLIBM
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

  real(kind=fPrec) pure function sinh_mb(arg)
    real(kind=fPrec), intent(in) :: arg
#ifndef CRLIBM
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

  real(kind=fPrec) pure function cos_mb(arg)
    real(kind=fPrec), intent(in) :: arg
#ifndef CRLIBM
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

  real(kind=fPrec) pure function acos_mb(arg)
    real(kind=fPrec), intent(in) :: arg
#ifndef CRLIBM
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

  real(kind=fPrec) pure function cosh_mb(arg)
    real(kind=fPrec), intent(in) :: arg
#ifndef CRLIBM
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

  real(kind=fPrec) pure function tan_mb(arg)
    real(kind=fPrec), intent(in) :: arg
#ifndef CRLIBM
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

  real(kind=fPrec) pure function atan_mb(arg)
    real(kind=fPrec), intent(in) :: arg
#ifndef CRLIBM
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

  real(kind=fPrec) pure function atan2_mb(y,x)
    real(kind=fPrec), intent(in) :: y,x
#ifndef CRLIBM
#ifdef NAGFOR
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

  real(kind=fPrec) pure function exp_mb(arg)
    real(kind=fPrec), intent(in) :: arg
#ifndef CRLIBM
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

  real(kind=fPrec) pure function log_mb(arg)
    real(kind=fPrec), intent(in) :: arg
#ifndef CRLIBM
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

  real(kind=fPrec) pure function log10_mb(arg)
    real(kind=fPrec), intent(in) :: arg
#ifndef CRLIBM
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

  ! ========================================================================== !
  !  ROUND NEAR SPEACIAL FUNCTIONS
  ! ========================================================================== !

#ifdef ROUND_NEAR
  real(kind=real64) pure function acos_rn(x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x
    if(x /= x) then ! Check if NaN
      acos_rn=x
    elseif(abs(x) == 0.0d0) then
      acos_rn=pi2
    else
      ! Try using (1-x)*(1+x) in case x is very small or close to 1. Write a test program?
      acos_rn=atan_rn(sqrt((1.0d0-x)*(1.0d0+x))/x)
      if(x < 0d0) then
        acos_rn=pi+acos_rn
      end if
    end if
  end function acos_rn

  real(kind=real64) pure function asin_rn(x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x
    if(x /= x) then ! Check if NaN
      asin_rn=x
      return
    end if
    if(abs(x) == 1.0d0) then
      asin_rn=sign(pi2,x)
    else
      ! Try using (1-x)*(1+x) in case x is very small or close to 1. Write a test program?
      asin_rn=atan_rn(x/sqrt((1.0d0-x)*(1.0d0+x)))
    end if
  end function asin_rn

  real(kind=real64) pure function atan2_rn(y,x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x,y
    if(x == 0d0) then
      if(y == 0d0) then
#ifdef NAGFOR
        atan2_rn=0d0
#else
        ! Let the internal atan2 handle this case according to ieee
        atan2_rn=atan2(y,x)
#endif
      else
        atan2_rn=sign(pi2,y)
      end if
    else
      if(y == 0d0) then
        if(x > 0d0) then
          ! Let the internal atan2 handle this case according to ieee
          atan2_rn=atan2(y,x)
        else
          atan2_rn=pi
        end if
      else
        atan2_rn=atan_rn(y/x)
        if(x < 0d0) then
          atan2_rn=sign(pi,y)+atan2_rn
        end if
      end if
    end if
  end function atan2_rn
#endif

  ! ========================================================================== !
  !  ROUND UP SPEACIAL FUNCTIONS
  ! ========================================================================== !

#ifdef ROUND_UP
  real(kind=real64) pure function acos_ru(x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x
    if(x /= x) then ! Check if NaN
      acos_ru=x
    elseif(abs(x) == 0.0d0) then
      acos_ru=pi2
    else
      acos_ru=atan_ru(sqrt((1.0d0-x)*(1.0d0+x))/x)
      if(x < 0d0) then
        acos_ru=pi+acos_ru
      end if
    end if
  end function acos_ru

  real(kind=real64) pure function asin_ru(x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x
    if(x /= x) then ! Check if NaN
      asin_ru=x
      return
    end if
    if(abs(x) == 1.0d0) then
      asin_ru=sign(pi2,x)
    else
      asin_ru=atan_ru(x/sqrt((1.0d0-x)*(1.0d0+x)))
    end if
  end function asin_ru

  real(kind=real64) pure function atan2_ru(y,x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x,y
    if(x == 0d0) then
      if(y == 0d0) then
#ifdef NAGFOR
        atan2_ru=0d0
#else
        ! Let the internal atan2 handle this case according to ieee
        atan2_ru=atan2(y,x)
#endif
      else
        atan2_ru=sign(pi2,y)
      end if
    else
      if(y == 0d0) then
        if(x > 0d0) then
          ! Let the internal atan2 handle this case according to ieee
          atan2_ru=atan2(y,x)
        else
          atan2_ru=pi
        end if
      else
        atan2_ru=atan_ru(y/x)
        if(x < 0d0) then
          atan2_ru=sign(pi,y)+atan2_ru
        end if
      end if
    end if
  end function atan2_ru
#endif

  ! ========================================================================== !
  !  ROUND DOWN SPEACIAL FUNCTIONS
  ! ========================================================================== !

#ifdef ROUND_DOWN
  real(kind=real64) pure function acos_rd(x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x
    if(x /= x) then ! Check if NaN
      acos_rd=x
    elseif(abs(x) == 0.0d0) then
      acos_rd=pi2
    else
      acos_rd=atan_rd(sqrt((1.0d0-x)*(1.0d0+x))/x)
      if(x < 0d0) then
        acos_rd=pi+acos_rd
      end if
    end if
  end function acos_rd

  real(kind=real64) pure function asin_rd(x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x
    if(x /= x) then ! Check if NaN
      asin_rd=x
      return
    end if
    if(abs(x) == 1.0d0) then
      asin_rd=sign(pi2,x)
    else
      asin_rd=atan_rd(x/sqrt((1.0d0-x)*(1.0d0+x)))
    end if
  end function asin_rd

  real(kind=real64) pure function atan2_rd(y,x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x,y
    if(x == 0d0) then
      if(y == 0d0) then
#ifdef NAGFOR
        atan2_rd=0d0
#else
        ! Let the internal atan2 handle this case according to ieee
        atan2_rd=atan2(y,x)
#endif
      else
        atan2_rd=sign(pi2,y)
      end if
    else
      if(y == 0d0) then
        if(x > 0d0) then
          ! Let the internal atan2 handle this case according to ieee
          atan2_rd=atan2(y,x)
        else
          atan2_rd=pi
        end if
      else
        atan2_rd=atan_rd(y/x)
        if(x < 0d0) then
          atan2_rd=sign(pi,y)+atan2_rd
        end if
      end if
    end if
  end function atan2_rd
#endif

  ! ========================================================================== !
  !  ROUND ZERO SPEACIAL FUNCTIONS
  ! ========================================================================== !

#ifdef ROUND_ZERO
  real(kind=real64) pure function acos_rz(x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x
    if(x /= x) then ! Check if NaN
      acos_rz=x
    elseif(abs(x) == 0.0d0) then
      acos_rz=pi2
    else
      acos_rz=atan_rz(sqrt((1.0d0-x)*(1.0d0+x))/x)
      if(x < 0d0) then
        acos_rz=pi+acos_rz
      end if
    end if
  end function acos_rz

  real(kind=real64) pure function asin_rz(x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x
    if(x /= x) then ! Check if NaN
      asin_rz=x
      return
    end if
    if(abs(x) == 1.0d0) then
      asin_rz=sign(pi2,x)
    else
      asin_rz=atan_rz(x/sqrt((1.0d0-x)*(1.0d0+x)))
    end if
  end function asin_rz

  real(kind=real64) pure function atan2_rz(y,x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x,y
    if(x == 0d0) then
      if(y == 0d0) then
#ifdef NAGFOR
        atan2_rz=0d0
#else
        ! Let the internal atan2 handle this case according to ieee
        atan2_rz=atan2(y,x)
#endif
      else
        atan2_rz=sign(pi2,y)
      end if
    else
      if(y == 0d0) then
        if(x > 0d0) then
          ! Let the internal atan2 handle this case according to ieee
          atan2_rz=atan2(y,x)
        else
          atan2_rz=pi
        end if
      else
        atan2_rz=atan_rz(y/x)
        if(x < 0d0) then
          atan2_rz=sign(pi,y)+atan2_rz
        end if
      end if
    end if
  end function atan2_rz
#endif

#endif

end module mathlib_bouncer
