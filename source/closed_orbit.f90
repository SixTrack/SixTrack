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
subroutine clorb(dpp, doWrite)

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
    if(ierro /= 0) return
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

end subroutine clorb

! ================================================================================================ !
!  CALCULATION OF THE SIX-DIMENSIONAL CLOSED ORBIT
!  Last modified: 2018-11-07
! ================================================================================================ !
subroutine clorda(nn,idummy,am)

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use matrix_inv
  use crcoall
  use parpro
  use mod_common
  use mod_commons
  use mod_commont
  use mod_commond

  implicit none

  integer,          intent(in)  :: nn
  integer,          intent(in)  :: idummy(nn)
  real(kind=fPrec), intent(out) :: am(nn,nn)

  integer i,i4,icheck,ii,j,j4,k,l,ll,nd2,nerror
  real(kind=fPrec) cloc(6),cor,coro,dc(6),dd(6),dlo(6),xx(6)
  character(len=6) chp(3),chd(3)

  save
!-----------------------------------------------------------------------
  nd2=2*ndimf
  write(lout,10010) nd2
  do l=1,nd2
    xx(l)=zero
    cloc(l)=zero
    dd(l)=zero
    dc(l)=zero
    dlo(l)=zero
    do i=1,nd2
      am(l,i)=zero
    end do
  end do
  chp(1)=' CLOX '
  chp(2)=' CLOY '
  chp(3)=' CLOS '
  chd(1)='  D-X '
  chd(2)='  D-Y '
  chd(3)='  D-S '
  cor=zero
  coro=1e38_fPrec
  ii=0
  if(ndimf == 3) then
    do l=1,2
      ll=2*l
      cloc(ll-1)=clo6(l)
      cloc(ll)=clop6(l)
    end do
    cloc(5)=clo6(3)
    cloc(6)=clop6(3)
    if(abs(dppoff) > pieni) cloc(6)=dppoff
  else
    do l=1,ndimf
      ll=2*l
      cloc(ll-1)=clo(l)
      cloc(ll)=clop(l)
    end do
    do l=ndimf+1,3
      ll=2*l
      cloc(ll-1)=zero
      cloc(ll)=zero
    end do
    cloc(6)=dps(1)
  end if
  do 80 ii=1,itco
    do l=1,2
      ll=2*l
      x(1,l)=cloc(ll-1)
      y(1,l)=cloc(ll)
    end do
    sigm(1)=cloc(5)
    dps(1)=cloc(6)
    call umlauda
    do i4=1,nd2
      do j4=1,nd2
        am(i4,j4)=aml6(i4,j4)
      end do
    end do
    call dinv(nd2,am,nd2,idummy,nerror)
    if(nerror /= 0) write(lout,*) ' ATTENTION, MATRIX SINGULAR '
    if(ndimf == 3) then
      do l=1,2
        ll=2*l
        xx(ll-1)=x(1,l)
        xx(ll)=y(1,l)
      end do
      xx(5)=sigm(1)
      xx(6)=dps(1)
    else
      do l=1,ndimf
        ll=2*l
        xx(ll-1)=x(1,l)
        xx(ll)=y(1,l)
      end do
      do l=ndimf+1,3
        ll=2*l
        xx(ll-1)=zero
        xx(ll)=zero
      end do
    end if
    do l=1,nd2
      dd(l)=cloc(l)-xx(l)
      dc(l)=abs(dd(l))
      if(l == 5) dc(5)=dc(5)*c1m2
    end do
    icheck=0
    do l=1,ndimf
      ll=2*l
      if(dc(ll-1) > dma) icheck=1
      if(dc(ll) > dmap) icheck=1
    end do
    if(icheck == 0) goto 90
    do k=1,nd2
      dlo(k)=zero
      do j=1,nd2
        dlo(k)=am(k,j)*dd(j)+dlo(k)
      end do
      if(abs(dppoff) > pieni) dlo(6)=zero
    end do
    write(lout,10020)
    cor=zero
    do l=1,ndimf
      ll=2*l
      write(lout,10060) chp(l),cloc(ll-1),cloc(ll)
      cor=cor+dc(ll-1)**2
    end do
    cor=sqrt(cor)
    if(ii == 1.or.cor < coro) then
      coro=cor
      do l=1,nd2
        cloc(l)=cloc(l)+dlo(l)
      end do
      if(ii /= itco) then
        write(lout,10030)
        do l=1,ndimf
          ll=2*l
          write(lout,10060) chp(l),cloc(ll-1),cloc(ll)
        end do
        write(lout,10080) ii,cor
      end if
    else
      write(lout,10040) nd2,ii
      goto 91
    end if
80   continue
  write(lout,10000) itco
  ii=itco
90   continue
  if(ii /= itco) then
    do k=1,nd2
      dlo(k)=zero
      do j=1,nd2
        dlo(k)=am(k,j)*dd(j)+dlo(k)
      end do
      if(abs(dppoff) > pieni) dlo(6)=zero
    end do
    write(lout,10020)
    cor=zero
    do l=1,ndimf
      ll=2*l
      write(lout,10060) chp(l),cloc(ll-1),cloc(ll)
      cor=cor+dc(ll-1)**2
    end do
    cor=sqrt(cor)
    if(cor < coro) then
      coro=cor
      do l=1,nd2
        cloc(l)=cloc(l)+dlo(l)
      end do
      write(lout,10030)
      do l=1,ndimf
        ll=2*l
        write(lout,10060) chp(l),cloc(ll-1),cloc(ll)
      end do
      write(lout,10080) ii,cor
    else
      write(lout,10040) nd2,ii
      goto 91
    end if
    do l=1,2
      ll=2*l
      x(1,l)=cloc(ll-1)
      y(1,l)=cloc(ll)
    end do
    sigm(1)=cloc(5)
    dps(1)=cloc(6)
    call umlauda
    do i4=1,nd2
      do j4=1,nd2
        am(i4,j4)=aml6(i4,j4)
      end do
    end do
    call dinv(nd2,am,nd2,idummy,nerror)
    if(nerror /= 0) write(lout,*) ' ATTENTION, MATRIX SINGULAR '
    if(ndimf == 3) then
      do l=1,2
        ll=2*l
        xx(ll-1)=x(1,l)
        xx(ll)=y(1,l)
      end do
      xx(5)=sigm(1)
      xx(6)=dps(1)
    else
      do l=1,ndimf
        ll=2*l
        xx(ll-1)=x(1,l)
        xx(ll)=y(1,l)
      end do
      do l=ndimf+1,3
        ll=2*l
        xx(ll-1)=zero
        xx(ll)=zero
      end do
    end if
    do l=1,nd2
      dc(l)=abs(cloc(l)-xx(l))
      if(l == 5) dc(5)=dc(5)*c1m2
    end do
  end if
  write(lout,10050) nd2,ii
  cor=zero
  do l=1,ndimf
    ll=2*l
    write(lout,10070) chp(l),cloc(ll-1),cloc(ll),chd(l),dc(ll-1),dc(ll)
    cor=cor+dc(ll-1)**2
  end do
  cor=sqrt(cor)
  write(lout,10080) ii,cor
91   continue
  if(ndimf == 3) then
    do l=1,2
      ll=2*l
      clo6(l)=cloc(ll-1)
      clop6(l)=cloc(ll)
    end do
    clo6(3)=cloc(5)
    clop6(3)=cloc(6)
  else
    do l=1,ndimf
      ll=2*l
      clo(l)=cloc(ll-1)
      clop(l)=cloc(ll)
    end do
  end if
!-----------------------------------------------------------------------
  return
10000 format(t10,'DA CLOSED ORBIT CALCULATION'/ t10,                    &
  &'MAXIMUM NUMBER OF ITERATIONS ACHIEVED--->',2x,i4/ t10,           &
  &'PROCEDURE MAY NOT HAVE CONVERGED')
10010 format(/131('-')/t10,'ENTERING ',i1,                              &
  &'-D DA CLOSED ORBIT CALCULATION'/)
10020 format(5x,'---- closed orbit before correction----')
10030 format(5x,'---- after DA correction----')
10040 format(/5x,'NO IMPROVEMENT OF ',i1,'-D DA CLOSED ORBIT ',         &
  &'CALCULATION IN ITERATION: ',i4/)
10050 format(t5,'SUCCESSFULL END OF ',i1,'-D DA CLOSED ORBIT ',         &
  &'CALCULATION IN ITERATION: ',i4/)
10060 format(5x,a6,1p,2(1x,g16.9))
10070 format(5x,a6,1p,2(1x,g16.9)/5x,a6,1p,2(1x,g16.9))
10080 format(5x,' ITERAT.=',i3,' ACCURACY=',d13.6/)
end subroutine clorda

end module closed_orbit
