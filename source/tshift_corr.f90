#ifdef SIXDA
! ================================================================================================ !
!  TUNESHIFT CORRECTIONS MODULE
!  Used by SixDA (mainda)
!  Last modified: 2018-06-01
!
!  Moved from main code into module June 2018, by VKBO.
! ================================================================================================ !
module tuneshift_corr

  use floatPrecision

  implicit none

  integer,          private, save :: iamp
  integer,          private, save :: x(10)
  real(kind=fPrec), private, save :: ham(0:3)
  real(kind=fPrec), private, save :: hama(0:4)
  real(kind=fPrec), private, save :: hamp(0:1)

  integer,          private, save :: jeltot1
  integer,          private, save :: maxa
  integer,          private, save :: maxp
  real(kind=fPrec), private, save :: hda1(0:3,3,0:3,0:4000)

  integer,          private, save :: jeltot2
  integer,          private, save :: nordp
  integer,          private, save :: nordm
  integer,          private, save :: norda
  real(kind=fPrec), private, save :: hda2(0:4,5,0:8000)
  real(kind=fPrec), private, save :: hdp(0:1,5,0:8000)

contains

!-----------------------------------------------------------------------
!---- PROGRAM FOR THE TUNESHIFT CORRECTIONS
!----
!----   =========>      FIRST & SECOND ORDER CORRECTIONS      <=========
!----   =========>            AMPLITUDE & CHROMATIC           <=========
!----   =========>                  EFFECTS                   <=========
!----
!-----------------------------------------------------------------------
subroutine coruord

  use numerical_constants
  use mathlib_bouncer
  use mod_commond, only : ncor,namp,nmom,coel
  use crcoall

  implicit none

  integer i,ifail,istate,iter,iuser,iwork,j,jaord,jbound,jcol,jcomp,jconf,jord,jpord,jrow,jsex,jvar,&
    k,kcol,l,liwork,lwork,mcor,n,nclin,ncnln,nconf,ndim2,nout,nrel,nrowa,nrowj,nrowr
  real a,bl,bu,c,cjac,clamda,objf,objgrd,r,work,user
  real(kind=fPrec) ainv,bmat,chia,chib,cmat,cvec,det,detinv,dvec,pi2in,sex,sgn
  external e04udm !,objfun1
  parameter(mcor = 10)
  parameter(ndim2 = 6)
  dimension a(2,10),cjac(1,1),c(1)
  dimension r(10,10),bu(20),bl(20),clamda(20),objgrd(10)
  dimension ainv(2,2),bmat(2,10),cmat(2,10),cvec(2),dvec(2)
  dimension work(450),user(500),sex(10),sgn(10,10)
  dimension istate(20),iwork(40),iuser(3)
  data sgn/100*one/ainv,bmat,cmat,cvec,dvec/48*zero/

  save

  pi2in=one/(eight*atan_mb(one))

  do i=0,3
    do j=1,3
      do k=0,3
        do l=0,4000
          hda1(i,j,k,l)=zero
        end do
      end do
    end do
  end do
  do j=1,10
    sgn(j,j)=-one
  end do

  ! SPECIFIES THE I/O UNITS FOR THE NAG ROUTINES
  nout=26
  call x04abf(1,nout)

  jeltot1=ncor
  jaord=namp
  jpord=nmom

  call readd1(user,jaord,jpord)

  if(jaord.eq.2.or.jpord.eq.3) goto 130

  ! DEFINES THE MATRIX WITH THE LINEAR CONSTRAINTS
  do jrow=1,2
    do jcol=1,jeltot1
      a(jrow,jcol)=hda1(jrow-1,1,1,3**(jcol-1))
    end do
    do jcol=1,jeltot1-2
      bmat(jrow,jcol)=-one*a(jrow,jcol+2)
    end do
    cvec(jrow)=-one*hda1(jrow-1,1,1,0)
  end do

  ! DEFINES THE RELATION BETWEEN THE FIRST TWO SEXTUPOLES AND THE OTHERS
  det=(a(1,1)*a(2,2)-a(2,1)*a(1,2))
  detinv=one/det
  ainv(1,1)=detinv*a(2,2)
  ainv(1,2)=-one*detinv*a(1,2)
  ainv(2,1)=-one*detinv*a(2,1)
  ainv(2,2)=detinv*a(1,1)

  do jrow=1,2
    do jcol=1,jeltot1
      do kcol=1,2
        cmat(jrow,jcol)=cmat(jrow,jcol)+ainv(jrow,kcol) *bmat(kcol,jcol)
      end do
    end do
    do jcol=1,2
      dvec(jrow)=dvec(jrow)+ainv(jrow,jcol)*cvec(jcol)
    end do
  end do

  ! WRITES ON THE EXIT FILE
130 write(lout,10000)
    write(lout,10010)
    write(lout,10020) jeltot1,jaord,jpord
    write(lout,10030)

  nrel=2
  nconf=jeltot1-2
  if(jaord.eq.2.or.jpord.eq.3) then
    nrel=0
    nconf=1
  end if

  ! DEFINES EXTRA PARAMETERS
  do jconf=1,jeltot1
    n=jeltot1
    nclin=nrel
    ncnln=0
    nrowa=2
    nrowj=1
    nrowr=10
    liwork=30
    lwork=450
    do jbound=1,n
      bu(jbound)=c1e10
      bl(jbound)=-c1e10
    end do
    do jbound=1,nclin
      bu(n+jbound)=-hda1(jbound-1,1,1,0)
      bl(n+jbound)=bu(n+jbound)
    end do
    do jvar=1,n
      x(jvar)=zero
    end do
    if(nrel.eq.2) then
      ! DEFINES THE INITIAL GUESS SO THAT IT SATISFIES THE LINEAR CONSTRAINTS
      do jvar=1,2
        do jcol=1,n-2
          x(jvar)=(x(jvar)+cmat(jvar,jcol)*sgn(jcol,jconf)) +dvec(jvar)
          x(jcol+2)=sgn(jcol,jconf)
        end do
      end do
    else if(nrel.eq.0) then
      do jvar=1,n
        x(jvar)=sgn(jvar,jconf)
      end do
    end if
    iuser(1)=jaord
    iuser(2)=jpord
    if(iuser(1).eq.0) iuser(3)=1
    if(iuser(2).eq.0) iuser(3)=iuser(1)+1
    ifail=-1
    call e04uef('MAJOR ITERATION LIMIT = 100')
    call e04ucf(n,nclin,ncnln,nrowa,nrowj,nrowr,a,bl,bu,e04udm,objfun1,iter,istate,c,&
                cjac,clamda,objf,objgrd,r,x,iwork,liwork,work,lwork,iuser,user,ifail)

    if(ifail.ne.0.and.ifail.ne.5) then
      write(lout,10040) ifail
      call prror(-1)
    end if

    do jsex=1,jeltot1
      sex(jsex)=x(jsex)
      write(lout,10050) coel(jsex),sex(jsex)
    end do

    ! COMPUTES THE NEW HAMILTONIAN IN DP/P AFTER THE CORRECTIONS
    do jord=1,maxp
      call hamilton1(1,jord)
      ! WRITES THE VALUES OF THE HAMILTONIAN
      write(lout,10060) jord
      write(lout,10070)
      write(lout,10080) hda1(1,1,jord,0),hda1(0,1,jord,0)
      write(lout,10090)
      write(lout,10080) ham(1),ham(0)
      write(lout,10100)
    end do

    ! COMPUTES THE NEW HAMILTONIAN IN AMP AFTER THE CORRECTIONS
    do jord=2,maxa
      call hamilton1(jord,0)
      ! COMPUTES THE FUNCTION CHI
      if(jord.eq.2) then
        chib=(pi2in/sqrt(three))*sqrt((((two*hda1(0,2,0,0)**2 &
          +hda1(1,2,0,0)**2)+two*hda1(2,2,0,0)**2)+hda1(0,2,0,0) *hda1(1,2,0,0))+hda1(1,2,0,0)*hda1(2,2,0,0))
        chia=(pi2in/sqrt(three))*sqrt((((two*ham(0)**2+ham(1)**2)+two *ham(2)**2)+ham(0)*ham(1))+ham(1)*ham(2))
      else if(jord.eq.3) then
        chib=(pi2in/sqrt(30.d0))*sqrt((((((((27.d0*hda1(3,3,0,0)**2 +5.d0*hda1(2,3,0,0)**2)+5.d0*hda1(1,3,0,0)**2)&
          +27.d0*hda1(0,3,0,0)**2)+(9.d0*hda1(3,3,0,0))*hda1(2,3,0,0))+(9.d0*hda1(1,3,0,0))*hda1(0,3,0,0))+(6.d0*   &
          hda1(2,3,0,0))*hda1(1,3,0,0))+(3.d0*hda1(3,3,0,0))*hda1(1,3,0,0)) +(3.d0*hda1(2,3,0,0))*hda1(0,3,0,0))
        chia=(pi2in/sqrt(30.d0))*sqrt((((((((27.d0*ham(3)**2 +5.d0*ham(2)**2)+5.d0*ham(1)**2)+27.d0*ham(0)**2) &
          +(9.d0*ham(3))*ham(2))+(9.d0*ham(1))*ham(0))+(6.d0*ham(2))*ham(1))+(3.d0*ham(3))*ham(1))+(3.d0*ham(2))*ham(0))
      end if

      ! WRITES THE VALUE OF THE HAMILTONIAN
      write(lout,10110) jord
      write(lout,10120)
      do jcomp=0,jord
        write(lout,10130)jcomp,jord-jcomp,hda1(jcomp,jord,0,0),jcomp,jord-jcomp,ham(jcomp)
      end do
      write(lout,10140) jord-1,chib,jord-1,chia
      write(lout,10100)
    end do
  end do

10000 format(//80('-')//t10,29('O')/t10,2('O'),25x,2('O')/t10,          &
     &'OO  TUNE-SHIFT CORRECTION  OO', /t10,2('O'),25x,2('O')/t10,29('O'&
     &)//80('-')//)
10010 format(//t26,'*** ORDER-BY-ORDER CORRECTIONS ***'//)
10020 format(t10,'NUMBER OF CORRECTOR ELEMENTS ',t48,i8/ t10,           &
     &'TUNE-SHIFT ORDER (AMPLITUDE) ',t48,i8/ t10,                      &
     &'TUNE-SHIFT ORDER (MOMENTUM) ',t48,i8)
10030 format(/,t10,'VALUES OF THE INTEGRATED GRADIENTS OF THE ',        &
     &'CORRECTOR MULTIPOLES',/)
10040 format(//,t10,' ERROR IN ROUTINE E04UCF. IFAIL = ',i8)
10050 format(/,t10,'CORRECTOR ELEMENT  - ',a16,' - ',4x,e21.14)
10060 format(//,t10,'HAMILTONIAN DEPENDENCE OF ORDER ' ,2x,i3,5x,       &
     &'MOMENTUM DEPENDENCE ')
10070 format(/,t10,'BEFORE CORRECTION ')
10080 format(/,'H_1,0    = ',2x,e16.8,7x,'H_0,1    = ', 2x,e16.8)
10090 format(/,t10,'AFTER CORRECTION ')
10100 format(//80('-'))
10110 format(//,t10,'HAMILTONIAN DEPENDENCE OF ORDER ' ,2x,i3,5x,       &
     &'AMPLITUDE DEPENDENCE',/)
10120 format(/,t10,'BEFORE CORRECTION ',20x,'AFTER CORRECTION ')
10130 format(/,t10,'H_',i1,',',i1,'    = ',e16.8,11x, 'H_',i1,',',i1,   &
     &'    = ',e16.8)
10140 format(/,t10,'CHI_',i1,',0  = ',e16.8,11x, 'CHI_',i1,',0  = ',e16.8)

end subroutine coruord

!-----------------------------------------------------------------------
!---- SUBROUTINE TO READ DATA
!-----------------------------------------------------------------------
subroutine readd1(user,jaord,jpord)

  use mathlib_bouncer
  use crcoall

  implicit none

  integer icont,ind,j,j1,j2,j3,j4,j5,j6,jaord,jcomp,jel,jord,jpord,maxcomp,njx,njx1,njz,njz1,nmax,np,ncoef,nord,point,kointer
  real user
  real(kind=fPrec) cc
  dimension ind(10),user(500)
#ifdef CRLIBM
  integer nchars
  parameter (nchars=160)
  character(len=nchars) ch
  character(len=nchars+nchars) ch1
  ! MAXF be kept in sync with value in function fround
  integer maxf,nofields
  parameter (maxf=30)
  parameter (nofields=20)
  character(len=maxf) fields(nofields)
  integer errno,nfields,nunit,lineno,nf
  real(kind=fPrec) fround
  data lineno /0/
#endif

  save

  nmax=40
  maxa=jaord+1
  maxp=jpord
  if(jpord.eq.0) then
    maxp=jaord+1
    maxcomp=maxa
  else if(jaord.eq.0) then
    maxa=jpord
    maxcomp=1
  end if

  rewind 23

  ! Unit 23 is opened with round='nearest' if fio is selected
#ifndef CRLIBM
10  continue
  read(23,'(I6,2X,G21.14,I5,4X,18(2I2,1X))',end=30) ncoef,cc,nord,njx,njx1,njz,njz1,np,(ind(jel),jel=1,jeltot1)
  read(23,*,end=30) cc
#else
10  continue
  read(23,*,end=30) ch
  ch1(:nchars+3)=ch(:nchars)//' / '
  lineno=lineno+1
  call splitfld(errno,23,lineno,nofields,nf,ch1,fields)
  if (nf.gt.0) then
    read (fields(1),*) ncoef
    nf=nf-1
  end if
  if (nf.gt.0) then
    cc=fround(errno,fields,2)
    nf=nf-1
  end if
  if (nf.gt.0) then
    read (fields(3),*) nord
    nf=nf-1
  end if
  if (nf.gt.0) then
    read (fields(4),*) njx
    nf=nf-1
  end if
  if (nf.gt.0) then
    read (fields(5),*) njx1
    nf=nf-1
  end if
  if (nf.gt.0) then
    read (fields(6),*) njz
    nf=nf-1
  end if
  if (nf.gt.0) then
    read (fields(7),*) njz1
    nf=nf-1
  end if
  if (nf.gt.0) then
    read (fields(8),*) np
    nf=nf-1
  end if
  do jel=1,jeltot1
    if (nf.gt.0) then
      read (fields(8+jel),*) ind(jel)
      nf=nf-1
    end if
  end do
  read(23,*,end=30) ch
  lineno=lineno+1
  ch1(:nchars+3)=ch(:nchars)//' / '
  call splitfld(errno,23,lineno,nofields,nf,ch1,fields)
  if (nf.gt.0) then
    cc=fround(errno,fields,1)
    nf=nf-1
  end if
#endif

  ! CODING IND IN BASE 3
  if(njx.eq.njx1.and.njz.eq.njz1.and.(njx+njz).le.maxa.and. np.le.maxp) then
    point=0
    do j=1,jeltot1
      point=point+ind(j)*3**(j-1)
    end do
    if(point.gt.4000) then
      write(lout,10000)
      write(lout,'(a)') "Problem with data in fort.23"
      call prror(-1)
    end if

    ! DATA PROCESSING
    hda1(njx,njx+njz,np,point)=cc+hda1(njx,njx+njz,np,point)
  end if
  goto 10

  ! DEFINES DATA FOR THE ROUTINE OBJFUN
30 continue
  if(jeltot1.eq.1) then
    do jcomp=0,maxcomp
      icont=0
      do jord=0,2
        j1=jord
        icont=icont+1
        kointer=j1
        user(jcomp*nmax+icont)=hda1(jcomp,jaord+1,jpord,kointer)
      end do
    end do
  else if(jeltot1.eq.2) then
    do jcomp=0,maxcomp
      icont=0
      do jord=0,2
        do j1=0,jord
          j2=jord-j1
          icont=icont+1
          kointer=j1+j2*3
          user(jcomp*nmax+icont)=hda1(jcomp,jaord+1,jpord,kointer)
        end do
      end do
    end do
  else if(jeltot1.eq.3) then
    do jcomp=0,maxcomp
      icont=0
      do jord=0,2
        do j1=0,jord
          do j2=0,jord-j1
            j3=jord-j1-j2
            icont=icont+1
            kointer=(j1+j2*3)+j3*3**2
            user(jcomp*nmax+icont)=hda1(jcomp,jaord+1,jpord, kointer)
          end do
        end do
      end do
    end do
  else if(jeltot1.eq.4) then
    do jcomp=0,maxcomp
      icont=0
      do jord=0,2
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              j4=jord-j1-j2-j3
              icont=icont+1
              kointer=((j1+j2*3)+j3*3**2)+j4*3**3
              user(jcomp*nmax+icont)=hda1(jcomp,jaord+1,jpord,kointer)
            end do
          end do
        end do
      end do
    end do
  else if(jeltot1.eq.5) then
    do jcomp=0,maxcomp
      icont=0
      do jord=0,2
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              do j4=0,jord-j1-j2-j3
                j5=jord-j1-j2-j3-j4
                icont=icont+1
                kointer=(((j1+j2*3)+j3*3**2)+j4*3**3)+j5*3**4
                user(jcomp*nmax+icont)=hda1(jcomp,jaord+1, jpord,kointer)
              end do
            end do
          end do
        end do
      end do
    end do
  else if(jeltot1.eq.6) then
    do jcomp=0,maxcomp
      icont=0
      do jord=0,2
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              do j4=0,jord-j1-j2-j3
                do j5=0,jord-j1-j2-j3-j4
                  j6=jord-j1-j2-j3-j4-j5
                  icont=icont+1
                  kointer=((((j1+j2*3)+j3*3**2)+j4*3**3)+j5*3**4)+j6*3**5
                  user(jcomp*nmax+icont)=hda1(jcomp,jaord+1, jpord,kointer)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end if

  return

10000 format(//,t10,' INDEX OUT OF BOUND IN ROUTINE READD ')

end subroutine readd1

!-----------------------------------------------------------------------
!---- COMPUTES THE VALUE OF THE HAMILTONIAN AFTER CORRECTIONS
!-----------------------------------------------------------------------
subroutine hamilton1(ja,jp)

  use numerical_constants
  use mathlib_bouncer

  implicit none

  integer j1,j2,j3,j4,j5,j6,ja,jcomp,jel,jord,jp,l,ncoef,kointer
  real(kind=fPrec) tham
  dimension tham(0:3)
  save

  do jcomp=0,3
    tham(jcomp)=zero
  end do

  if(jp.eq.0) then
    ncoef=ja
  else
    ncoef=1
  end if

  select case (jeltot1)
    case (1)
      do jord=0,2
        j1=jord
        kointer=j1
        do l=0,ncoef
          tham(l)=tham(l)+hda1(l,ja,jp,kointer)*(x(1)**j1)
        end do
      end do

    case (2)
      do jord=0,2
        do j1=0,jord
          j2=jord-j1
          kointer=j1+j2*3
          do l=0,ncoef
            tham(l)=tham(l)+(hda1(l,ja,jp,kointer)*(x(1)**j1))*(x(2)**j2) !hr04
          end do
        end do
      end do

    case (3)
      do jord=0,2
        do j1=0,jord
          do j2=0,jord-j1
            j3=jord-j1-j2
            kointer=(j1+j2*3)+j3*3**2 ! hr04
            do l=0,ncoef
              tham(l)=tham(l)+((hda1(l,ja,jp,kointer)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3) ! hr04
            end do
          end do
        end do
      end do

    case (4)
      do jord=0,2
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              j4=jord-j1-j2-j3
              kointer=((j1+j2*3)+j3*3**2)+j4*3**3 ! hr04
              do l=0,ncoef
                tham(l) = tham(l)+(((hda1(l,ja,jp,kointer)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4) ! hr04
              end do
            end do
          end do
        end do
      end do

    case (5)
      do jord=0,2
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              do j4=0,jord-j1-j2-j3
                j5=jord-j1-j2-j3-j4
                kointer=(((j1+j2*3)+j3*3**2)+j4*3**3)+j5*3**4 ! hr04
                do l=0,ncoef
                  tham(l)=tham(l)+((((hda1(l,ja,jp,kointer)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5) ! hr04
                end do
              end do
            end do
          end do
        end do
      end do

    case (6)
      do jord=0,2
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              do j4=0,jord-j1-j2-j3
                do j5=0,jord-j1-j2-j3-j4
                  j6=jord-j1-j2-j3-j4-j5
                  kointer=((((j1+j2*3)+j3*3**2)+j4*3**3)+j5*3**4)+j6*3**5 ! hr04
                  do l=0,ncoef
                    tham(l)=tham(l)+(((((hda1(l,ja,jp,kointer)*&
                            (x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5))*(x(6)**j6) ! hr04
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

  end select

  do jel=0,ncoef
    ham(jel)=tham(jel)
  end do

  return

end subroutine hamilton1

!-----------------------------------------------------------------------
!---- ROUTINE TO COMPUTE THE VALUE OF THE FUNCTION AND OF ITS
!---- DERIVATIVES
!-----------------------------------------------------------------------
subroutine objfun1(mode,n,x,objf,objgrd,nstate,iuser,user)

  use mathlib_bouncer
  use numerical_constants

  implicit none

  integer icont,iuser,j1,j2,j3,j4,j5,j6,jel,jord,jvar,l,mode,n,nmax,nstate,kointer
  real user
  real(kind=fPrec) fder,fun,objf,objgrd,x
  dimension iuser(*),x(10),objgrd(10),user(*),fun(0:3),fder(0:3,10)
  save

  nmax=40

  do jel=0,3
    do jvar=1,n
      fder(jel,jvar)=zero
    end do
    fun(jel)=zero
  end do

  if(n.eq.1) then
    do l=0,iuser(3)
      icont=0
      do jord=0,2
        j1=jord
        kointer=j1
        icont=icont+1
        fun(l)=fun(l)+user(l*nmax+icont)*(x(1)**j1)
        fder(l,1)=fder(l,1)+(user(l*nmax+icont)*real(j1,fPrec))*(x(1)**(j1-1))
      end do
    end do
  else if(n.eq.2) then
    do l=0,iuser(3)
      icont=0
      do jord=0,2
        do j1=0,jord
          j2=jord-j1
          kointer=j1+j2*3
          icont=icont+1
          fun(l)=fun(l)+(user(l*nmax+icont)*(x(1)**j1))*(x(2)**j2)
          fder(l,1)=fder(l,1)+((user(l*nmax+icont)*real(j1,fPrec))*(x(1)**(j1-1)))*(x(2)**j2)
          fder(l,2)=fder(l,2)+((user(l*nmax+icont)*real(j2,fPrec))*(x(1)**j1))*(x(2)**(j2-1))
        end do
      end do
    end do
  else if(n.eq.3) then
    do l=0,iuser(3)
      icont=0
      do jord=0,2
        do j1=0,jord
          do j2=0,jord-j1
            j3=jord-j1-j2
            kointer=(j1+j2*3)+j3*3**2
            icont=icont+1
            fun(l)=fun(l)+((user(l*nmax+icont)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3)
            fder(l,1)=fder(l,1)+(((user(l*nmax+icont)*real(j1,fPrec))*(x(1)**(j1-1)))*(x(2)**j2))*(x(3)**j3)
            fder(l,2)=fder(l,2)+(((user(l*nmax+icont)*real(j2,fPrec))*(x(1)**j1))*(x(2)**(j2-1)))*(x(3)**j3)
            fder(l,3)=fder(l,3)+(((user(l*nmax+icont)*real(j3,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**(j3-1))
          end do
        end do
      end do
    end do
  else if(n.eq.4) then
    do l=0,iuser(3)
      icont=0
      do jord=0,2
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              j4=jord-j1-j2-j3
              kointer=j1+j2*3+j3*3**2+j4*3**3
              icont=icont+1
              fun(l)=fun(l)+(((user(l*nmax+icont)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4)
              fder(l,1)=fder(l,1)+((((user(l*nmax+icont)*real(j1,fPrec))*(x(1)**(j1-1)))*(x(2)**j2))*(x(3 )**j3))*(x(4)**j4)
              fder(l,2)=fder(l,2)+((((user(l*nmax+icont)*real(j2,fPrec))*(x(1)**j1))*(x(2)**(j2-1)))*(x(3 )**j3))*(x(4)**j4)
              fder(l,3)=fder(l,3)+((((user(l*nmax+icont)*real(j3,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3 )**(j3-1)))*(x(4)**j4)
              fder(l,4)=fder(l,4)+((((user(l*nmax+icont)*real(j4,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3 )**j3))*(x(4)**(j4-1))
            end do
          end do
        end do
      end do
    end do
  else if(n.eq.5) then
    do l=0,iuser(3)
      icont=0
      do jord=0,2
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              do j4=0,jord-j1-j2-j3
                j5=jord-j1-j2-j3-j4
                kointer=(((j1+j2*3)+j3*3**2)+j4*3**3)+j5*3**4
                icont=icont+1
                fun(l)=fun(l)+((((user(l*nmax+icont)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5)
                fder(l,1)=fder(l,1)+(((((user(l*nmax+icont) &
                  *real(j1,fPrec))*(x(1)**(j1-1)))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5)
                fder(l,2)=fder(l,2)+(((((user(l*nmax+icont) &
                  *real(j2,fPrec))*(x(1)**j1))*(x(2)**(j2-1)))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5)
                fder(l,3)=fder(l,3)+(((((user(l*nmax+icont) &
                  *real(j3,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**(j3-1)))*(x(4)**j4))*(x(5)**j5)
                fder(l,4)=fder(l,4)+(((((user(l*nmax+icont) &
                  *real(j4,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**(j4-1)))*(x(5)**j5)
                fder(l,5)=fder(l,5)+(((((user(l*nmax+icont) &
                  *real(j5,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**(j5-1))
              end do
            end do
          end do
        end do
      end do
    end do
  else if(n.eq.6) then
    do l=0,iuser(3)
      icont=0
      do jord=0,2
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              do j4=0,jord-j1-j2-j3
                do j5=0,jord-j1-j2-j3-j4
                  j6=jord-j1-j2-j3-j4-j5
                  kointer=((((j1+j2*3)+j3*3**2)+j4*3**3)+j5*3**4)+j6*3**5
                  icont=icont+1
                  fun(l)=fun(l)+(((((user(l*nmax+icont) *(x(1)**j1)) *(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5))*(x(6)**j6)
                  fder(l,1)=fder(l,1)+((((((user(l*nmax+icont) &
                    *real(j1,fPrec))*(x(1)**(j1-1)))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5))*(x(6)**j6)
                  fder(l,2)=fder(l,2)+((((((user(l*nmax+icont) &
                    *real(j2,fPrec))*(x(1)**j1))*(x(2)**(j2-1)))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5))*(x(6)**j6)
                  fder(l,3)=fder(l,3)+((((((user(l*nmax+icont) &
                    *real(j3,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**(j3-1)))*(x(4)**j4))*(x(5)**j5))*(x(6)**j6)
                  fder(l,4)=fder(l,4)+((((((user(l*nmax+icont) &
                    *real(j4,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**(j4-1)))*(x(5)**j5))*(x(6)**j6)
                  fder(l,5)=fder(l,5)+((((((user(l*nmax+icont) &
                    *real(j5,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**(j5-1)))*(x(6)**j6)
                  fder(l,6)=fder(l,6)+((((((user(l*nmax+icont) &
                    *real(j6,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5))*(x(6)**(j6-1))
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end if

  if(iuser(2).eq.0) then
    if(iuser(1).eq.1) then
      objf=(((2.d0*fun(0)**2+fun(1)**2)+2.d0*fun(2)**2)+fun(0)*fun(1))+fun(2)*fun(1)
      do jvar=1,n
        objgrd(jvar)=((4.d0*fun(0)+fun(1))*fder(0,jvar)+ &
          ((2.d0*fun(1)+fun(0))+fun(2))*fder(1,jvar))+ (4.d0*fun(2)+fun(1))*fder(2,jvar)
      end do
    else
      objf=(((((((27.d0*fun(3)**2+5.d0*fun(2)**2)+5.d0*fun(1)**2)+27.d0*fun(0)**2)+9.d0*fun(3)*fun(2)) &
        +9.d0*fun(1)*fun(0))+6.d0*fun(2)*fun(1))+3.d0*fun(3)*fun(1))+3.d0*fun(2)*fun(0)
      do jvar=1,n
        objgrd(jvar)=((((54.d0*fun(3)+9.d0*fun(2))+3.d0*fun(1)) *fder(3,jvar) &
          +(((10.d0*fun(2)+9.d0*fun(3))+6.d0*fun(1))+3.d0* fun(0))*fder(2,jvar))+(((10.d0*fun(1) &
          +9.d0*fun(0))+6.d0*fun(2))+3.d0*fun(3))*fder(1,jvar))+((54.d0*fun(0)+9.d0*fun(1))+3.d0*fun(2))*fder(0,jvar)
      end do
    end if
  else
    objf=fun(0)**2+fun(1)**2
    do jvar=1,n
      objgrd(jvar)=(two*fun(0))*fder(0,jvar) +(two*fun(1))*fder(1,jvar)
    end do
  end if

  return

end subroutine objfun1

!-----------------------------------------------------------------------
!---- PROGRAM FOR THE TUNESHIFT CORRECTIONS
!----
!----   =========>            CHROMATIC EFFECTS                <========
!----   =========>           GLOBAL MINIMIZATION               <========
!----
!-----------------------------------------------------------------------
subroutine coruglo

  use numerical_constants
  use mathlib_bouncer
  use mod_commond, only : ncor,dpmax,nmom1,nmom2,weig1,weig2,coel
  use crcoall

  implicit none

  integer i,ifail,istate,iter,itype,iuser,iwork,j,jbound,jcol,jcomp,jconf,jord,jrow,jsex,jvar,kcol,l,&
      liwork,lwork,mcor,n,nclin,ncnln,ndim2,nout,nrowa,nrowj,nrowr
  real a,bl,bu,c,cjac,clamda,objf,objgrd,r,user,work
  real(kind=fPrec) ainv,bmat,chia,chib,cmat,cvec,delta,detinv,dvec,pi2in,sex,sgn,value
  external e04udm !,objfun2
  parameter(mcor = 10)
  parameter(ndim2 = 6)

  dimension a(2,10),cjac(1,1),c(1)
  dimension r(10,10),bu(20),bl(20),clamda(20),objgrd(10)
  dimension ainv(2,2),bmat(2,10),cmat(2,10),cvec(2),dvec(2)
  dimension work(450),user(500),sex(10),sgn(10,10)
  dimension istate(20),iwork(40),iuser(2)
  data sgn/100*one/ainv,bmat,cmat,cvec,dvec/48*zero/

  save

  pi2in=one/(eight*atan_mb(one))

  do i=0,4
    do j=1,5
      do l=0,8000
        hda2(i,j,l)=zero
        if(i.le.1) hdp(i,j,l)=zero
      end do
    end do
  end do
  do j=1,10
    sgn(j,j)=-one
  end do

  ! SPECIFIES THE I/O UNITS FOR THE NAG ROUTINES
  nout=26
  call x04abf(1,nout)
  jeltot2=ncor
  delta=dpmax
  nordm=nmom1
  nordp=nmom2
  value=zero
  itype=0
  call readd2(user)

  ! DEFINES THE MATRIX WITH THE LINEAR CONSTRAINTS
  do jrow=1,2
    do jcol=1,jeltot2
      a(jrow,jcol)=hdp(jrow-1,1,(nordp+1)**(jcol-1))
    end do
    do jcol=1,jeltot2-2
      bmat(jrow,jcol)=-one*a(jrow,jcol+2)
    end do
    cvec(jrow)=-one*hdp(jrow-1,1,0)
  end do

  ! DEFINES THE RELATION BETWEEN THE FIRST TWO SEXTUPOLES AND THE OTHERS
  detinv=one/(a(1,1)*a(2,2)-a(2,1)*a(1,2))
  ainv(1,1)=detinv*a(2,2)
  ainv(1,2)=(-one*detinv)*a(1,2)
  ainv(2,1)=(-one*detinv)*a(2,1)
  ainv(2,2)=detinv*a(1,1)

  do jrow=1,2
    do jcol=1,jeltot2
      do kcol=1,2
        cmat(jrow,jcol)=cmat(jrow,jcol)+ainv(jrow,kcol) *bmat(kcol,jcol)
      end do
    end do
    do jcol=1,2
      dvec(jrow)=dvec(jrow)+ainv(jrow,jcol)*cvec(jcol)
    end do
  end do

  ! WRITES ON THE EXIT FILE
  write(lout,10000)
  write(lout,10010)
  write(lout,10020) jeltot2,nordm,nordp,delta,weig1,weig2,value
  write(lout,10030)

  do jconf=1,jeltot2-2
    ! INITIALIZATION
    iuser(1)=nordp
    iuser(2)=nordm
    user(1)=weig1
    user(2)=weig2
    user(3)=value
    user(4)=delta*c1e3
    user(5)=one/(user(4)**(2*nordm+1))

    ! DEFINES EXTRA PARAMETERS
    n=jeltot2
    nclin=2
    ncnln=0
    nrowa=2
    nrowj=1
    nrowr=10
    liwork=30
    lwork=450

    do jbound=1,n
      bu(jbound)=c1e2
      bl(jbound)=-c1e2
    end do

    do jbound=1,nclin
      bu(n+jbound)=-one*hdp(jbound-1,1,0)
      bl(n+jbound)=-one*hdp(jbound-1,1,0)
    end do

    do jvar=1,n
      x(jvar)=zero
    end do

    ! DEFINES THE INITIAL GUESS SO THAT IT SATISFIES THE LINEAR CONSTRAINTS
    do jvar=1,2
      do jcol=1,n-2
        x(jvar)=(x(jvar)+cmat(jvar,jcol)*sgn(jcol,jconf))+dvec(jvar)
        x(jcol+2)=sgn(jcol,jconf)
      end do
    end do

    ifail=-1
    call e04uef('MAJOR ITERATION LIMIT = 100')
    call e04ucf(n,nclin,ncnln,nrowa,nrowj,nrowr,a,bl,bu,e04udm,objfun2,iter,istate,c,cjac,clamda,&
      objf,objgrd,r,x,iwork,liwork,work,lwork,iuser,user,ifail)

    if(ifail.ne.0.and.ifail.ne.5) then
      write(lout,10040) ifail
      call prror(-1)
    end if

    do jsex=1,jeltot2
      sex(jsex)=x(jsex)
      write(lout,10050) coel(jsex),sex(jsex)
    end do

    ! COMPUTES THE NEW HAMILTONIAN IN DP/P AFTER THE CORRECTIONS
    do jord=1,nordp
      iamp=0
      call hamilton2(jord)

      ! WRITES THE VALUES OF THE HAMILTONIAN
      write(lout,10060) jord
      write(lout,10070)
      write(lout,10080) hdp(1,jord,0),hdp(0,jord,0)
      write(lout,10090)
      write(lout,10080) hamp(1),hamp(0)
      write(lout,10100)
    end do

    ! COMPUTES THE NEW HAMILTONIAN IN AMP AFTER THE CORRECTIONS
    do jord=2,norda
      iamp=1
      call hamilton2(jord)

      ! COMPUTES THE FUNCTION CHI
      if(jord.eq.2) then
        chib=(pi2in/sqrt(3.d0))*sqrt((((2.d0*hda2(0,2,0)**2 +hda2(1,2,0)**2)+2.d0*&
          hda2(2,2,0)**2)+hda2(0,2,0) *hda2(1,2,0))+hda2(1,2,0)*hda2(2,2,0))
        chia=(pi2in/sqrt(3.d0))*sqrt((((2.d0*hama(0)**2+hama(1)**2)+2.d0*hama(2)**2)+hama(0)*hama(1))+hama(1)*hama(2))
      else if(jord.eq.3) then
        chib=(pi2in/sqrt(30.d0))*sqrt((((((((27.d0*hda2(3,3,0)**2 +5.d0*hda2 &
          (2,3,0)**2)+5.d0*hda2(1,3,0)**2)+27.d0*hda2(0,3,0)**2)+9.d0*hda2    &
          (3,3,0)*hda2(2,3,0))+9.d0*hda2(1,3,0)*hda2(0,3,0))+6.d0*hda2         &
          (2,3,0)*hda2(1,3,0))+3.d0*hda2(3,3,0)*hda2(1,3,0))+3.d0*hda2         &
          (2,3,0)*hda2(0,3,0))
        chia=(pi2in/sqrt(30.d0))*sqrt((((((((27.d0*hama(3)**2 +5.d0*hama(2)&
          **2)+5.d0*hama(1)**2)+27.d0*hama(0)**2)+9.d0*hama(3)*hama(2))    &
          +9.d0*hama(1)*hama(0))+6.d0*hama(2)*hama(1))+3.d0*hama(3)        &
          *hama(1))+3.d0*hama(2)*hama(0))
      end if

      ! WRITES THE VALUE OF THE HAMILTONIAN
      write(lout,10110) jord
      write(lout,10120)
      do jcomp=0,jord
        write(lout,10130) jcomp, jord-jcomp, hda2(jcomp,jord,0), jcomp, jord-jcomp, hama(jcomp)
      end do
      write(lout,10140) jord-1,chib,jord-1,chia
      write(lout,10100)
    end do
  end do

10000 format(//80('-')//t10,29('O')/t10,2('O'),25x,2('O')/t10,          &
     &'OO  TUNE-SHIFT CORRECTION  OO', /t10,2('O'),25x,2('O')/t10,29('O'&
     &)//80('-')//)
10010 format(//t26,'*** GLOBAL CHROMATIC CORRECTIONS ***'//)
10020 format(t10,'NUMBER OF CORRECTOR ELEMENTS ',t48,i8/ t10,           &
     &'MINIMUM ORDER OF THE MINIMIZATION ',t48,i8/ t10,                 &
     &'MAXIMUM ORDER OF THE MINIMIZATION ',t48,i8/ t10,                 &
     &'MAXIMUM VALUE OF MOMENTUM DEVIATION ',t48,f9.3/ t10,             &
     &'WEIGHT HORIZONTAL TUNE ',t48,f9.3/ t10,'WEIGHT VERTICAL TUNE ',  &
     &t48,f9.3/ t10,'QUADRATIC COEFFICIENT IN DP/P',t48,f9.3)
10030 format(//,t10,'VALUES OF THE INTEGRATED GRADIENTS OF THE ',       &
     &'CORRECTOR MULTIPOLES',//,t10,'GLOBAL MOMENTUM CORRECTIONS')
10040 format(//,t10,' ERROR IN ROUTINE E04UCF. IFAIL = ',i8)
10050 format(//,t10,'CORRECTOR ELEMENT  - ',a16,' - ',4x,e21.14)
10060 format(///,t10,'HAMILTONIAN DEPENDENCE OF ORDER ' ,2x,i3,5x,      &
     &'MOMENTUM DEPENDENCE ')
10070 format(//,t10,'BEFORE CORRECTION ')
10080 format(//,'H_1,0    = ',2x,e16.8,7x,'H_0,1    = ', 2x,e16.8)
10090 format(//,t10,'AFTER CORRECTION ')
10100 format(//80('-'))
10110 format(///,t10,'HAMILTONIAN DEPENDENCE OF ORDER ' ,2x,i3,5x,      &
     &'AMPLITUDE DEPENDENCE',/)
10120 format(//,t10,'BEFORE CORRECTION ',20x,'AFTER CORRECTION ')
10130 format(//,t10,'H_',i1,',',i1,'    = ',e16.8,11x, 'H_',i1,',',i1,  &
     &'    = ',e16.8)
10140 format(//,t10,'CHI_',i1,',0  = ',e16.8,11x, 'CHI_',i1,',0  = ',e16&
     &.8)

end subroutine coruglo

!-----------------------------------------------------------------------
!---- SUBROUTINE TO READ DATA
!-----------------------------------------------------------------------
subroutine readd2(user)

  use mathlib_bouncer
  use crcoall

  implicit none

  integer icont,ind,j,j1,j2,j3,j4,j5,j6,jcomp,jel,jord,jp,ncoef,njx,njx1,njz,njz1,nor,np,point,kointer
  real user
  real(kind=fPrec) cc
  dimension ind(10),user(500)

#ifdef CRLIBM
  integer nchars
  parameter (nchars=160)
  character(len=nchars) ch
  character(len=nchars+nchars) ch1
  integer nofields
  parameter (nofields=20)
  character(len=nchars) fields(nofields)
  integer errno,nfields,nunit,lineno,maxf,nf
  real(kind=fPrec) fround
  data lineno /0/
#endif

  save

  rewind 23

  ! Unit 23 is opened round='nearest' if fio is selected
#ifndef CRLIBM
10 continue
  read(23,*,end=40) ncoef,cc,nor,njx,njx1,njz,njz1,np,(ind(jel),jel=1,jeltot2)
  read(23,*,end=40) cc
#else
10 continue
  nunit=23
  read(23,*,end=40) ch
  lineno=lineno+1
  ch1(:nchars+3)=ch(:nchars)//' / '
  call splitfld(errno,23,lineno,nofields,nf,ch1,fields)
  if (nf.gt.0) then
    read (fields(1),*) ncoef
    nf=nf-1
  end if
  if (nf.gt.0) then
    cc=fround(errno,fields,2)
    nf=nf-1
  end if
  if (nf.gt.0) then
    read (fields(3),*) nor
    nf=nf-1
  end if
  if (nf.gt.0) then
    read (fields(4),*) njx
    nf=nf-1
  end if
  if (nf.gt.0) then
    read (fields(5),*) njx1
    nf=nf-1
  end if
  if (nf.gt.0) then
    read (fields(6),*) njz
    nf=nf-1
  end if
  if (nf.gt.0) then
    read (fields(7),*) njz1
    nf=nf-1
  end if
  if (nf.gt.0) then
    read (fields(8),*) np
    nf=nf-1
  end if
  do jel=1,jeltot2
    if (nf.gt.0) then
      read (fields(8+jel),*) ind(jel)
      nf=nf-1
    end if
  end do
  read(23,*,end=40) ch
  lineno=lineno+1
  ch1(:nchars+3)=ch(:nchars)//' / '
  call splitfld(errno,23,lineno,nofields,nf,ch1,fields)
  if (nf.gt.0) then
    cc=fround(errno,fields,1)
    nf=nf-1
  end if
#endif

  ! CODING IND IN BASE NORDP+1
  if(njx.eq.njx1.and.njz.eq.njz1.and.(njx+njz).eq.1.and. np.gt.0.and.np.le.nordp) then
    point=0
    do j=1,jeltot2
      point=point+ind(j)*(nordp+1)**(j-1)
    end do
    if(point.gt.8000) then
      write(lout,10000)
      call prror(-1)
    end if

    ! DATA PROCESSING
    hdp(njx,np,point)=cc+hdp(njx,np,point)
  else if(njx.eq.njx1.and.njz.eq.njz1.and.(njx+njz).lt.10.and. np.eq.0) then
    point=0
    norda=njx+njz
    do j=1,jeltot2
      point=point+ind(j)*(nordp+1)**(j-1)
    end do
    if(point.gt.8000) then
      write(lout,10000)
      call prror(-1)
    end if

    ! DATA PROCESSING
    hda2(njx,njx+njz,point)=cc+hda2(njx,njx+njz,point)
  end if

  goto 10

  ! DEFINES DATA FOR THE ROUTINE OBJFUN
40 continue
  icont=5
  if(jeltot2.eq.1) then
    do jcomp=0,1
      do jp=nordm,nordp
        do jord=0,jp
          j1=jord
          icont=icont+1
          kointer=j1
          user(icont)=hdp(jcomp,jp,kointer)
        end do
      end do
    end do
  else if(jeltot2.eq.2) then
    do jcomp=0,1
      do jp=nordm,nordp
        do jord=0,jp
          do j1=0,jord
            j2=jord-j1
            icont=icont+1
            kointer=j1+j2*(nordp+1)
            user(icont)=hdp(jcomp,jp,kointer)
          end do
        end do
      end do
    end do
  else if(jeltot2.eq.3) then
    do jcomp=0,1
      do jp=nordm,nordp
        do jord=0,jp
          do j1=0,jord
            do j2=0,jord-j1
              j3=jord-j1-j2
              icont=icont+1
              kointer=(j1+j2*(nordp+1))+j3*(nordp+1)**2
              user(icont)=hdp(jcomp,jp,kointer)
            end do
          end do
        end do
      end do
    end do
  else if(jeltot2.eq.4) then
    do jcomp=0,1
      do jp=nordm,nordp
        do jord=0,jp
          do j1=0,jord
            do j2=0,jord-j1
              do j3=0,jord-j1-j2
                j4=jord-j1-j2-j3
                icont=icont+1
                kointer=((j1+j2*(nordp+1))+j3*(nordp+1)**2)+j4*(nordp+1)**3
                user(icont)=hdp(jcomp,jp,kointer)
              end do
            end do
          end do
        end do
      end do
    end do
  else if(jeltot2.eq.5) then
    do jcomp=0,1
      do jp=nordm,nordp
        do jord=0,jp
          do j1=0,jord
            do j2=0,jord-j1
              do j3=0,jord-j1-j2
                do j4=0,jord-j1-j2-j3
                  j5=jord-j1-j2-j3-j4
                  icont=icont+1
                  kointer=(((j1+j2*(nordp+1))+j3*(nordp+1)**2)+j4*(nordp+1)**3)+j5*(nordp+1)**4
                  user(icont)=hdp(jcomp,jp,kointer)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  else if(jeltot2.eq.6) then
    do jcomp=0,1
      do jp=nordm,nordp
        do jord=0,jp
          do j1=0,jord
            do j2=0,jord-j1
              do j3=0,jord-j1-j2
                do j4=0,jord-j1-j2-j3
                  do j5=0,jord-j1-j2-j3-j4
                    j6=jord-j1-j2-j3-j4-j6
                    icont=icont+1
                    kointer=((((j1+j2*(nordp+1))+j3*(nordp+1)**2)+j4*(nordp+1)**3)+j5*(nordp+1)**4)+j6*(nordp+1)**5
                    user(icont)=hdp(jcomp,jp,kointer)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end if

  return

10000 format(//,t10,' INDEX OUT OF BOUND IN ROUTINE READD ')

end subroutine readd2

!-----------------------------------------------------------------------
!---- COMPUTES THE VALUE OF THE HAMILTONIAN AFTER CORRECTIONS
!-----------------------------------------------------------------------
subroutine hamilton2(jp)

  use numerical_constants
  use mathlib_bouncer

  implicit none

  integer j,j1,j2,j3,j4,j5,j6,jel,jord,jp,l,kointer
  real(kind=fPrec) thama,thamp

  dimension thamp(0:1),thama(0:4)
  save

  if(iamp.eq.0) then
    thamp(0)=zero
    thamp(1)=zero
    if(jeltot2.eq.1) then
      do jord=0,jp
        j1=jord
        kointer=j1
        do l=0,1
          thamp(l)=thamp(l)+hdp(l,jp,kointer)*(x(1)**j1)
        end do
      end do
    else if(jeltot2.eq.2) then
      do jord=0,jp
        do j1=0,jord
          j2=jord-j1
          kointer=j1+j2*(nordp+1)
          do l=0,1
            thamp(l)=thamp(l)+hdp(l,jp,kointer)*(x(1)**j1)*(x(2)**j2)
          end do
        end do
      end do
    else if(jeltot2.eq.3) then
      do jord=0,jp
        do j1=0,jord
          do j2=0,jord-j1
            j3=jord-j1-j2
            kointer=(j1+j2*(nordp+1))+j3*(nordp+1)**2
            do l=0,1
              thamp(l)=thamp(l)+((hdp(l,jp,kointer)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3)
            end do
          end do
        end do
      end do
    else if(jeltot2.eq.4) then
      do jord=0,jp
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              j4=jord-j1-j2-j3
              kointer=((j1+j2*(nordp+1))+j3*(nordp+1)**2)+j4*(nordp+1)**3
              do l=0,1
                thamp(l)=thamp(l)+(((hdp(l,jp,kointer)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4)
              end do
            end do
          end do
        end do
      end do
    else if(jeltot2.eq.5) then
      do jord=0,jp
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              do j4=0,jord-j1-j2-j3
                j5=jord-j1-j2-j3-j4
                kointer=(((j1+j2*(nordp+1))+j3*(nordp+1)**2)+j4*(nordp+1)**3)+j5*(nordp+1)**4
                do l=0,1
                  thamp(l)=thamp(l)+((((hdp(l,jp,kointer)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5)
                end do
              end do
            end do
          end do
        end do
      end do
    else if(jeltot2.eq.6) then
      do jord=0,jp
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              do j4=0,jord-j1-j2-j3
                do j5=0,jord-j1-j2-j3-j4
                  j6=jord-j1-j2-j3-j4-j5
                  kointer=((((j1+j2*(nordp+1))+j3*(nordp+1)**2)+j4*(nordp+1)**3)+j5*(nordp+1)**4)+j6*(nordp+1)**5
                  do l=0,1
                    thamp(l)=thamp(l)+(((((hdp(l,jp,kointer)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5))*(x(6)**j6)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end if
    do jel=0,1
      hamp(jel)=thamp(jel)
    end do
  else
    do j=0,4
      thama(j)=0.d0
    end do
    if(jeltot2.eq.1) then
      do jord=0,jp
        j1=jord
        kointer=j1
        do l=0,jp
          thama(l)=thama(l)+hda2(l,jp,kointer)*(x(1)**j1)
        end do
      end do
    else if(jeltot2.eq.2) then
      do jord=0,jp
        do j1=0,jord
          j2=jord-j1
          kointer=j1+j2*(nordp+1)
          do l=0,jp
            thama(l)=thama(l)+hda2(l,jp,kointer)*(x(1)**j1)*(x(2)**j2)
          end do
        end do
      end do
    else if(jeltot2.eq.3) then
      do jord=0,jp
        do j1=0,jord
          do j2=0,jord-j1
            j3=jord-j1-j2
            kointer=(j1+j2*(nordp+1))+j3*(nordp+1)**2
            do l=0,jp
              thama(l)=thama(l)+((hda2(l,jp,kointer)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3)
            end do
          end do
        end do
      end do
    else if(jeltot2.eq.4) then
      do jord=0,jp
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              j4=jord-j1-j2-j3
              kointer=((j1+j2*(nordp+1))+j3*(nordp+1)**2)+j4*(nordp+1)**3
              do l=0,jp
                thama(l)=thama(l)+(((hda2(l,jp,kointer)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4)
              end do
            end do
          end do
        end do
      end do
    else if(jeltot2.eq.5) then
      do jord=0,jp
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              do j4=0,jord-j1-j2-j3
                j5=jord-j1-j2-j3-j4
                kointer=(((j1+j2*(nordp+1))+j3*(nordp+1)**2)+j4*(nordp+1)**3)+j5*(nordp+1)**4
                do l=0,jp
                  thama(l)=thama(l)+((((hda2(l,jp,kointer)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5)
                end do
              end do
            end do
          end do
        end do
      end do
    else if(jeltot2.eq.6) then
      do jord=0,jp
        do j1=0,jord
          do j2=0,jord-j1
            do j3=0,jord-j1-j2
              do j4=0,jord-j1-j2-j3
                do j5=0,jord-j1-j2-j3-j4
                  j6=jord-j1-j2-j3-j4-j5
                  kointer=((((j1+j2*(nordp+1))+j3*(nordp+1)**2)+j4*(nordp+1)**3)+j5*(nordp+1)**4)+j6*(nordp+1)**5
                  do l=0,jp
                    thama(l)=thama(l)+(((((hda2(l,jp,kointer)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5))*(x(6)**j6)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end if
    do jel=0,jp
      hama(jel)=thama(jel)
    end do
  end if

  return

end subroutine hamilton2

!-----------------------------------------------------------------------
!---- ROUTINE TO COMPUTE THE VALUE OF THE FUNCTION AND OF ITS
!---- DERIVATIVES
!-----------------------------------------------------------------------
subroutine objfun2(mode,n,x,objf,objgrd,nstate,iuser,user)

  use numerical_constants
  use mathlib_bouncer

  implicit none

  integer icont,iuser,j1,j2,j3,j4,j5,j6,jel,jord,jp,jvar,kord,l,mode,n,nstate,kointer
  real user
  real(kind=fPrec) add1,add2,fder,fun,objf,objgrd,sgn,tunedx,tunedy,tunex,tuney,weight,x
  dimension x(10),objgrd(10),user(*),fun(0:1,10),fder(0:1,10,10)
  dimension iuser(*),tunedx(10),tunedy(10)

  save

  do jel=0,1
    do jord=1,iuser(1)
      do jvar=1,n
        fder(jel,jord,jvar)=0.d0
      end do
      fun(jel,jord)=0.d0
    end do
  end do

  icont=5
  if(n.eq.1) then
    do l=0,1
      do jp=iuser(2),iuser(1)
        do jord=0,jp
          j1=jord
          kointer=j1
          icont=icont+1
          fun(l,jp)=fun(l,jp)+user(icont)*(x(1)**j1)
          fder(l,jp,1)=fder(l,jp,1)+(user(icont)*real(j1,fPrec))*(x(1)**(j1-1))
        end do
      end do
    end do
  else if(n.eq.2) then
    do l=0,1
      do jp=iuser(2),iuser(1)
        do jord=0,jp
          do j1=0,jord
            j2=jord-j1
            kointer=j1+j2*(iuser(1)+1)
            icont=icont+1
            fun(l,jp)=fun(l,jp)+(user(icont)*(x(1)**j1))*(x(2)**j2)
            fder(l,jp,1)=fder(l,jp,1)+((user(icont)*real(j1,fPrec))*(x(1)**(j1-1)))*(x(2)**j2)
            fder(l,jp,2)=fder(l,jp,2)+((user(icont)*real(j2,fPrec))*(x(1)**j1))*(x(2)**(j2-1))
          end do
        end do
      end do
    end do
  else if(n.eq.3) then
    do l=0,1
      do jp=iuser(2),iuser(1)
        do jord=0,jp
          do j1=0,jord
            do j2=0,jord-j1
              j3=jord-j1-j2
              kointer=(j1+j2*(iuser(1)+1))+j3*(iuser(1)+1)**2
              icont=icont+1
              fun(l,jp)=fun(l,jp)+((user(icont)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3)
              fder(l,jp,1)=fder(l,jp,1)+(((user(icont)*real(j1,fPrec))*(x(1)**(j1-1)))*(x(2)**j2))*(x(3)**j3)
              fder(l,jp,2)=fder(l,jp,2)+(((user(icont)*real(j2,fPrec))*(x(1)**j1))*(x(2)**(j2-1)))*(x(3)**j3)
              fder(l,jp,3)=fder(l,jp,3)+(((user(icont)*real(j3,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**(j3-1))
            end do
          end do
        end do
      end do
    end do
  else if(n.eq.4) then
    do l=0,1
      do jp=iuser(2),iuser(1)
        do jord=0,jp
          do j1=0,jord
            do j2=0,jord-j1
              do j3=0,jord-j1-j2
                j4=jord-j1-j2-j3
                kointer=((j1+j2*(iuser(1)+1))+j3*(iuser(1)+1)**2)+j4*(iuser(1)+1)**3
                icont=icont+1
                fun(l,jp)=fun(l,jp)+(((user(icont)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4)
                fder(l,jp,1)=fder(l,jp,1)+((((user(icont)*real(j1,fPrec))*(x(1)**(j1-1)))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4)
                fder(l,jp,2)=fder(l,jp,2)+((((user(icont)*real(j2,fPrec))*(x(1)**j1))*(x(2)**(j2-1)))*(x(3)**j3))*(x(4)**j4)
                fder(l,jp,3)=fder(l,jp,3)+((((user(icont)*real(j3,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**(j3-1)))*(x(4)**j4)
                fder(l,jp,4)=fder(l,jp,4)+((((user(icont)*real(j4,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**(j4-1))
              end do
            end do
          end do
        end do
      end do
    end do
  else if(n.eq.5) then
    do l=0,1
      do jp=iuser(2),iuser(1)
        do jord=0,jp
          do j1=0,jord
            do j2=0,jord-j1
              do j3=0,jord-j1-j2
                do j4=0,jord-j1-j2-j3
                  j5=jord-j1-j2-j3-j4
                  kointer=(((j1+j2*(iuser(1)+1))+j3*(iuser(1)+1)**2)+j4*(iuser(1)+1)**3)+j5*(iuser(1)+1)**4
                  icont=icont+1
                  fun(l,jp)=fun(l,jp)+((((user(icont)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5)
                  fder(l,jp,1)=fder(l,jp,1)+(((((user(icont)* &
                    real(j1,fPrec))*(x(1)**(j1-1)))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5)
                  fder(l,jp,2)=fder(l,jp,2)+(((((user(icont)* &
                    real(j2,fPrec))*(x(1)**j1))*(x(2)**(j2-1)))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5)
                  fder(l,jp,3)=fder(l,jp,3)+(((((user(icont)* &
                    real(j3,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**(j3-1)))*(x(4)**j4))*(x(5)**j5)
                  fder(l,jp,4)=fder(l,jp,4)+(((((user(icont)* &
                    real(j4,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**(j4-1)))*(x(5)**j5)
                  fder(l,jp,5)=fder(l,jp,5)+(((((user(icont)* &
                    real(j5,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**(j5-1))
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  else if(n.eq.6) then
    do l=0,1
      do jp=iuser(2),iuser(1)
        do jord=0,jp
          do j1=0,jord
            do j2=0,jord-j1
              do j3=0,jord-j1-j2
                do j4=0,jord-j1-j2-j3
                  do j5=0,jord-j1-j2-j3-j4
                    j6=jord-j1-j2-j3-j4-j5
                    kointer=((((j1+j2*(iuser(1)+1))+j3*(iuser(1)+1)**2)+j4*(iuser(1)+1)**3) +j5*(iuser(1)+1)**4)+j6*(iuser(1)+1)**4
                    icont=icont+1
                    fun(l,jp)=fun(l,jp)+(((((user(icont)*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5))*(x(6)**j6)
                    fder(l,jp,1)=fder(l,jp,1)+((((((user(icont)* &
                      real(j1,fPrec))*(x(1)**(j1-1)))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5))*(x(6)**j6)
                    fder(l,jp,2)=fder(l,jp,2)+((((((user(icont)* &
                      real(j2,fPrec))*(x(1)**j1))*(x(2)**(j2-1)))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5))*(x(6)**j6)
                    fder(l,jp,3)=fder(l,jp,3)+((((((user(icont)* &
                      real(j3,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**(j3-1)))*(x(4)**j4))*(x(5)**j5))*(x(6)**j6)
                    fder(l,jp,4)=fder(l,jp,4)+((((((user(icont)* &
                      real(j4,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**(j4-1)))*(x(5)**j5))*(x(6)**j6)
                    fder(l,jp,5)=fder(l,jp,5)+((((((user(icont)* &
                      real(j5,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**(j5-1)))*(x(6)**j6)
                    fder(l,jp,6)=fder(l,jp,6)+((((((user(icont)* &
                      real(j6,fPrec))*(x(1)**j1))*(x(2)**j2))*(x(3)**j3))*(x(4)**j4))*(x(5)**j5))*(x(6)**(j6-1))
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end if

  tunex=zero
  tuney=zero
  do jvar=1,n
    tunedx(jvar)=zero
    tunedy(jvar)=zero
  end do

  do jord=iuser(2),iuser(1)
    do kord=iuser(2),iuser(1)
      sgn=one
      if(((jord+kord)/2)*2.ne.(jord+kord)) sgn=-one
      add1=one/real((jord+kord)+1,fPrec)
      add2=user(3)/real((jord+kord)+3,fPrec)
      weight=((user(4)**((jord+kord)+1))*(one+sgn))*(add1+add2)*user(5)
      tunex=tunex+(fun(0,jord)*fun(0,kord))*weight
      tuney=tuney+(fun(1,jord)*fun(1,kord))*weight
      do jvar=1,n
        tunedx(jvar)=tunedx(jvar)+(fun(0,jord)*fder(0,kord, jvar))*weight
        tunedy(jvar)=tunedy(jvar)+(fun(1,jord)*fder(1,kord, jvar))*weight
      end do
    end do
  end do

  objf=user(1)**2*tunex+user(2)**2*tuney
  do jvar=1,n
    objgrd(jvar)=(two*user(1)**2)*tunedx(jvar)+ (two*user(2)**2)*tunedy(jvar)
  end do

  return

end subroutine objfun2

! #ifndef NAGLIB
! subroutine e04ucf(n,nclin,ncnln,lda,ldcj,ldr,a,bl,bu,confun,objfun,iter,ierroe,c,cjac,clamda,objf,  &
!                   objgrd,r,x,iwork,liwork,work,lwork,iuser,user,ifail)
!   implicit none
!   integer n,nclin,ncnln,lda,ldcj,ldr,iter,ierroe(*),istate(n+nclin+ncnln),liwork,iwork(liwork),lwork,iuser(*),ifail,x(n)
!   real a(lda,*),bl(n+nclin+ncnln),bu(n+nclin+ncnln),c(*),cjac(ldcj,*),clamda(n+nclin+ncnln),objf,&
!     objgrd(n),r(ldr,n),work(lwork),user(*)
!   external confun,objfun
!   save
!   return
! end subroutine e04ucf

! subroutine e04uef(c1)
!   implicit none
!   character c1
!   save
!   return
! end subroutine e04uef

! subroutine e04udm(c1)
!   implicit none
!   character c1
!   save
!   return
! end subroutine e04udm

! subroutine x04abf(n1,n2)
!   implicit none
!   integer n1,n2
!   save
!   return
! end subroutine x04abf
! #endif

end module tuneshift_corr
#endif