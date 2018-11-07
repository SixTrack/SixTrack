! ================================================================================================ !
!  CLOSED ORBIT MODULE
!  Last modified: 2018-11-07
! ================================================================================================ !
module closed_orbit

  use floatPrecision

  implicit none

contains

! ================================================================================================ !
!  CALCULATION OF THE CLOSED ORBIT
!  Last modified: 2018-11-07
! ================================================================================================ !
subroutine calcClosedOrbit(dpp, doWrite)

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use crcoall
  use parpro
  use mod_settings
  use mod_common
  use mod_commons
  use mod_commont

  implicit none

  real(kind=fPrec), intent(in) :: dpp
  logical,          intent(in) :: doWrite

  integer ierr,ii,l,ll
  real(kind=fPrec) am(4,4),cor,dclo(2),dclop(2),dcx,dcxp,dcz,dczp,det,dx(2),dy(2),x0(2),x1(2),y0(2),y1(2)

  ! save ! Saving DPP?

  ierro=0
  do l=1,2
    clo(l)  = dpp*di0(l)
    clop(l) = dpp*dip0(l)
    dx(l)   = c1e6
    dy(l)   = c1e6
  end do

  call envar(dpp)
  call umlauf(dpp,1,ierr)

  ierro=ierr
  if(ierro /= 0) return
  do ii=1,itco
    dcx  = abs(dx(1))
    dcxp = abs(dy(1))
    dcz  = abs(dx(2))
    dczp = abs(dy(2))
    if(dcx <= dma .and. dcz <= dma .and. dcxp <= dmap .and. dczp <= dmap) then
      if(doWrite) then
        goto 50
      else
        return
      end if
    end if

    do l=1,2
      x(1,l) = clo(l)
      y(1,l) = clop(l)
      x0(l)  = x(1,l)
      y0(l)  = y(1,l)
    end do

    call matrix(dpp,am)
    if(ierro.ne.0) return
    do l=1,2
      ll       = 2*l
      x1(l)    = x(1,l)
      y1(l)    = y(1,l)
      det      = (two-am(ll-1,ll-1))-am(ll,ll)
      dx(l)    = x0(l)-x1(l)
      dy(l)    = y0(l)-y1(l)
      dclo(l)  = (dx(l)*(am(ll,ll)-one)-dy(l)*am(ll-1,ll))/det
      dclop(l) = (dy(l)*(am(ll-1,ll-1)-one)-dx(l)*am(ll,ll-1))/det
      clo(l)   = clo(l)+dclo(l)
      clop(l)  = clop(l)+dclop(l)
    end do
  end do

  if(doWrite .eqv. .false.) return

50 continue
  cor = c1e3 * sqrt(dcx**2 + dcz**2)
  if(st_print .and. ncorru /= 1) then
    write(lout,"(a)")               "    Closed Orbit Entry"
    write(lout,"(a,f13.8)")         "      dpp      = ",dpp
    write(lout,"(a,f13.8,a,f13.8)") "      x, xp    = ",clo(1),", ",clop(1)
    write(lout,"(a,f13.8,a,f13.8)") "      y, yp    = ",clo(2),", ",clop(2)
    write(lout,"(a,e13.6,a,i0,a)")  "      accuracy = ",cor," (",ii," iterations)"
    write(lout,"(a)")               ""
  end if

end subroutine calcClosedOrbit

end module closed_orbit
