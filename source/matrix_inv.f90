module matrix_inv

  implicit none

contains

! ================================================================================================ !
! REPLACES A BY ITS INVERSE.
! (PARAMETERS AS FOR DEQINV.)
! CALLS ... DFACT, DFINV, F010PR, ABEND.
! ================================================================================================ !
subroutine dinv(n,a,idim,ir,ifail)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

  implicit none

  integer idim,ifail,jfail,k,kprnt,n
  integer ir
  real(kind=fPrec) t1,t2,t3
  real(kind=fPrec) a,det,temp,s,c11,c12,c13,c21,c22,c23,c31,c32,c33
  character(len=6) name
  dimension ir(n),a(idim,n)
  data name/'DINV'/,kprnt/0/
  save

  ! TEST FOR PARAMETER ERRORS.
  if((n.lt.1).or.(n.gt.idim)) goto 7

  ! TEST FOR N.LE.3.
  if(n.gt.3) goto 6
  ifail=0
  if(n.lt.3) goto 4

  ! N=3 CASE.
  ! COMPUTE COFACTORS.
  c11=a(2,2)*a(3,3)-a(2,3)*a(3,2)
  c12=a(2,3)*a(3,1)-a(2,1)*a(3,3)
  c13=a(2,1)*a(3,2)-a(2,2)*a(3,1)
  c21=a(3,2)*a(1,3)-a(3,3)*a(1,2)
  c22=a(3,3)*a(1,1)-a(3,1)*a(1,3)
  c23=a(3,1)*a(1,2)-a(3,2)*a(1,1)
  c31=a(1,2)*a(2,3)-a(1,3)*a(2,2)
  c32=a(1,3)*a(2,1)-a(1,1)*a(2,3)
  c33=a(1,1)*a(2,2)-a(1,2)*a(2,1)
  t1=abs(a(1,1))
  t2=abs(a(2,1))
  t3=abs(a(3,1))

  ! (SET TEMP=PIVOT AND DET=PIVOT*DET.)
  if(t1.ge.t2) goto 1
  if(t3.ge.t2) goto 2
  ! (PIVOT IS A21)
  temp=a(2,1)
  det=c13*c32-c12*c33
  goto 3
1 if(t3.ge.t1) goto 2
  ! (PIVOT IS A11)
  temp=a(1,1)
  det=c22*c33-c23*c32
  goto 3
  ! (PIVOT IS A31)
2 temp=a(3,1)
  det=c23*c12-c22*c13

  ! SET ELEMENTS OF INVERSE IN A.
3 if(det.eq.zero) goto 8
  s=temp/det
  a(1,1)=s*c11
  a(1,2)=s*c21
  a(1,3)=s*c31
  a(2,1)=s*c12
  a(2,2)=s*c22
  a(2,3)=s*c32
  a(3,1)=s*c13
  a(3,2)=s*c23
  a(3,3)=s*c33
  return

4 if(n.lt.2) goto 5
  ! N=2 CASE BY CRAMERS RULE.
  det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
  if(det.eq.zero) goto 8
  s=one/det
  c11   =s*a(2,2)
  a(1,2)=(-one*s)*a(1,2)
  a(2,1)=(-one*s)*a(2,1)
  a(2,2)=s*a(1,1)
  a(1,1)=c11
  return

  !  N=1 CASE.
5 if(a(1,1).eq.zero) goto 8
  a(1,1)=one/a(1,1)
  return

  ! N.GT.3 CASES.  FACTORIZE MATRIX AND INVERT.

6 call dfact(n,a,idim,ir,ifail,det,jfail)
  if(ifail.ne.0) return
  call dfinv(n,a,idim,ir)
  return

  ! ERROR EXITS.
7 ifail=+1
  call f010pr(name,n,idim,k,kprnt)
  return

8 ifail=-1
  return

end subroutine dinv

! ================================================================================================ !
subroutine kermtr(ercode,log,mflag,rflag)

  use floatPrecision
  use mathlib_bouncer
  use crcoall

  implicit none

  integer i,kounte,l,lgfile,limitm,limitr,log,logf
  parameter(kounte = 27)
  character(len=6)    ercode,   code(kounte)
  logical             mflag,    rflag
  integer             kntm(kounte),       kntr(kounte)

  data      logf      /  0  /
  data      code(1), kntm(1), kntr(1)  / 'C204.1', 255, 255 /
  data      code(2), kntm(2), kntr(2)  / 'C204.2', 255, 255 /
  data      code(3), kntm(3), kntr(3)  / 'C204.3', 255, 255 /
  data      code(4), kntm(4), kntr(4)  / 'C205.1', 255, 255 /
  data      code(5), kntm(5), kntr(5)  / 'C205.2', 255, 255 /
  data      code(6), kntm(6), kntr(6)  / 'C305.1', 255, 255 /
  data      code(7), kntm(7), kntr(7)  / 'C308.1', 255, 255 /
  data      code(8), kntm(8), kntr(8)  / 'C312.1', 255, 255 /
  data      code(9), kntm(9), kntr(9)  / 'C313.1', 255, 255 /
  data      code(10),kntm(10),kntr(10) / 'C336.1', 255, 255 /
  data      code(11),kntm(11),kntr(11) / 'C337.1', 255, 255 /
  data      code(12),kntm(12),kntr(12) / 'C341.1', 255, 255 /
  data      code(13),kntm(13),kntr(13) / 'D103.1', 255, 255 /
  data      code(14),kntm(14),kntr(14) / 'D106.1', 255, 255 /
  data      code(15),kntm(15),kntr(15) / 'D209.1', 255, 255 /
  data      code(16),kntm(16),kntr(16) / 'D509.1', 255, 255 /
  data      code(17),kntm(17),kntr(17) / 'E100.1', 255, 255 /
  data      code(18),kntm(18),kntr(18) / 'E104.1', 255, 255 /
  data      code(19),kntm(19),kntr(19) / 'E105.1', 255, 255 /
  data      code(20),kntm(20),kntr(20) / 'E208.1', 255, 255 /
  data      code(21),kntm(21),kntr(21) / 'E208.2', 255, 255 /
  data      code(22),kntm(22),kntr(22) / 'F010.1', 255,   0 /
  data      code(23),kntm(23),kntr(23) / 'F011.1', 255,   0 /
  data      code(24),kntm(24),kntr(24) / 'F012.1', 255,   0 /
  data      code(25),kntm(25),kntr(25) / 'F406.1', 255,   0 /
  data      code(26),kntm(26),kntr(26) / 'G100.1', 255, 255 /
  data      code(27),kntm(27),kntr(27) / 'G100.2', 255, 255 /
  save

  log  =  logf
  do i=1, kounte
    if(ercode .eq. code(i)) goto 21
  end do
  write(lout,1000)  ercode
  write(lout,"(a)") "KERNLIB> Library Error"
  call prror
  return

21 continue
  rflag  =  kntr(i) .ge. 1
  if(rflag  .and.  (kntr(i) .lt. 255))  kntr(i)  =  kntr(i) - 1
  mflag  =  kntm(i) .ge. 1
  if(mflag  .and.  (kntm(i) .lt. 255))  kntm(i)  =  kntm(i) - 1
  if(.not. rflag)  then
    if(logf .lt. 1)  then
      write(lout,1001)  code(i)
    else
      write(logf,1001)  code(i)
    end if
  end if
  if(mflag .and. rflag)  then
    if(logf .lt. 1)  then
      write(lout,1002)  code(i)
    else
      write(logf,1002)  code(i)
    end if
  end if
  return

1000  format('KERNLIB LIBRARY ERROR. ERROR CODE ',a6,' NOT RECOGNIZED BY KERMTR ERROR MONITOR. RUN ABORTED.')
1001  format(/' ***** RUN TERMINATED BY CERN LIBRARY ERROR CONDITION ',a6)
1002  format(/' ***** CERN LIBRARY ERROR CONDITION ',a6)

end subroutine kermtr

! This routine is unused
! ================================================================================================ !
!     REPLACES A BY ITS INVERSE.
!     (PARAMETERS AS FOR REQINV.)
!     CALLS ... RFACT, RFINV, F010PR, ABEND.
! ================================================================================================ !
subroutine rinv(n,a,idim,ir,ifail)

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer

  implicit none
  integer idim,ifail,ir,jfail,k,kprnt,n
  real(kind=fPrec) t1,t2,t3,a,det,temp,s,c11,c12,c13,c21,c22,c23,c31,c32,c33

  character(len=6) name
  dimension ir(n),a(idim,n)
  data name/'RINV'/,kprnt/0/
  save

  ! TEST FOR PARAMETER ERRORS.
  if((n.lt.1).or.(n.gt.idim)) goto 7

  ! TEST FOR N.LE.3.
  if(n.gt.3) goto 6
  ifail=0
  if(n.lt.3) goto 4

  ! N=3 CASE.
  ! COMPUTE COFACTORS.
  c11=a(2,2)*a(3,3)-a(2,3)*a(3,2)
  c12=a(2,3)*a(3,1)-a(2,1)*a(3,3)
  c13=a(2,1)*a(3,2)-a(2,2)*a(3,1)
  c21=a(3,2)*a(1,3)-a(3,3)*a(1,2)
  c22=a(3,3)*a(1,1)-a(3,1)*a(1,3)
  c23=a(3,1)*a(1,2)-a(3,2)*a(1,1)
  c31=a(1,2)*a(2,3)-a(1,3)*a(2,2)
  c32=a(1,3)*a(2,1)-a(1,1)*a(2,3)
  c33=a(1,1)*a(2,2)-a(1,2)*a(2,1)
  t1=abs(a(1,1))
  t2=abs(a(2,1))
  t3=abs(a(3,1))

  ! (SET TEMP=PIVOT AND DET=PIVOT*DET.)
  if(t1.ge.t2) goto 1
  if(t3.ge.t2) goto 2
  ! (PIVOT IS A21)
  temp=a(2,1)
  det=c13*c32-c12*c33
  goto 3
1 if(t3.ge.t1) goto 2
  ! (PIVOT IS A11)
  temp=a(1,1)
  det=c22*c33-c23*c32
  goto 3
  ! (PIVOT IS A31)
2 temp=a(3,1)
  det=c23*c12-c22*c13
  ! SET ELEMENTS OF INVERSE IN A.
3 if(det.eq.zero) goto 8
  s=temp/det
  a(1,1)=s*c11
  a(1,2)=s*c21
  a(1,3)=s*c31
  a(2,1)=s*c12
  a(2,2)=s*c22
  a(2,3)=s*c32
  a(3,1)=s*c13
  a(3,2)=s*c23
  a(3,3)=s*c33
  return

4 if(n.lt.2) goto 5

  ! N=2 CASE BY CRAMERS RULE.
  det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
  if(det.eq.zero) goto 8
  s=one/det
  c11   =s*a(2,2)
  a(1,2)=(-one*s)*a(1,2)
  a(2,1)=(-one*s)*a(2,1)
  a(2,2)=s*a(1,1)
  a(1,1)=c11
  return

  ! N=1 CASE.
5 if(a(1,1).eq.zero) goto 8
  a(1,1)=one/a(1,1)
  return

  ! N.GT.3 CASES.  FACTORIZE MATRIX AND INVERT.
6 call rfact(n,a,idim,ir,ifail,det,jfail)
  if(ifail.ne.0) return
  call rfinv(n,a,idim,ir)
  return

  ! ERROR EXITS.
7 ifail=+1
  call f010pr(name,n,idim,k,kprnt)
  return

8 ifail=-1
  return

end subroutine rinv

! ================================================================================================ !
!  PRINT ROUTINE FOR PARAMETER ERRORS IN MATRIX SUBROUTINES $EQINV,
!  $EQN, $INV (WHERE $ IS A LETTER SPECIFYING THE ARITHMETIC TYPE).
!  NAME         (CHARACTER(6)) NAME OF THE CALLING ROUTINE.
!  N,IDIM,K     PARAMETERS OF THE CALLING ROUTINE (WITH K=0 IF K IS
!               NOT TO BE PRINTED).
!  KPRNT        PRINT FLAG FOR K (K IS NOT PRINTED IF KPRNT=0).
! ================================================================================================ !
subroutine f010pr(name,n,idim,k,kprnt)

  use floatPrecision
  use mathlib_bouncer
  use crcoall

  implicit none

  integer idim,k,kprnt,lgfile,n
  character(len=6) name
  logical mflag,rflag
  save

  call kermtr('F010.1',lgfile,mflag,rflag)
  if(mflag) then
    if(lgfile.eq.0)  then
      if(kprnt.eq.0) write(lout,2000) name,n,idim
      if(kprnt.ne.0) write(lout,2001) name,n,idim,k
    else
      if(kprnt.eq.0) write(lgfile,2000) name,n,idim
      if(kprnt.ne.0) write(lgfile,2001) name,n,idim,k
    end if
  end if
  if(.not. rflag) then
    write(lerr,"(a)") "KERNLIB> ERROR F010PR: "//name
    call prror
  end if
  return

2000 format(7x,'subroutine ',a6,' ... parameter error (n.lt.1 or n.gt.idim).',6x,'n =',i4,6x,'idim =',i4,'.')
2001 format(7x,'subroutine ',a6,' ... parameter error (n.lt.1 or n.gt.idim or k.lt.1).',6x,'n =',i4,6x,'idim =',i4,6x,'k =',i4,'.')
end subroutine f010pr

subroutine rfact(n,a,idim,ir,ifail,det,jfail)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

  implicit none

  integer i,idim,ifail,imposs,ir,j,jfail,jm1,jover,jp1,jrange,junder,k,l,n,normal,nxch

  real(kind=fPrec) a,det,g1,g2,p,q,t,tf,x,y
  real(kind=fPrec) s11,s12

  character(len=6) hname

  dimension ir(*),a(idim,*)
  !      data      g1, g2              /  1.e-37,  1.e37  /
  data      g1, g2              /  1.0e-37_fPrec,  1.0e37_fPrec  /
  data      hname               /  ' RFACT'  /
  data      normal, imposs      /  0, -1  /
  data      jrange, jover, junder  /  0, +1, -1  /
  save

  if(idim .ge. n  .and.  n .gt. 0)  goto 110
    call tmprnt(hname,n,idim,0)
    return
110 ifail  =  normal
    jfail  =  jrange
    nxch   =  0
    det    =  one
    do 144    j  =  1, n
120   k = j
      p  =  pivotf(a(j,j))
      if(j .eq. n)  goto 122
      jp1  =  j+1
      do 121    i  =  jp1, n
        q = pivotf(a(i,j))
        if(q .le. p)  goto 121
        k  =  i
        p  =  q
121   continue
      if(k .ne. j)  goto 123
122   if(p .gt. zero)  goto 130
      det    =  zero
      ifail  =  imposs
      jfail  =  jrange
      return
123   do 124    l  =  1, n
        tf      =  a(j,l)
        a(j,l)  =  a(k,l)
        a(k,l)  =  tf
124   continue
      nxch      =  nxch + 1
      ir(nxch)  =  ipairf(j,k)
130   det     =  det * a(j,j)
      a(j,j)  =  one / a(j,j)
      t  =  sizef(det)
      if(t .lt. g1)  then
        det    =  zero
        if(jfail .eq. jrange)  jfail  =  junder
      elseif(t .gt. g2)  then
        det    =  one
        if(jfail .eq. jrange)  jfail  =  jover
      endif
      if(j .eq. n)  goto 144
      jm1  =  j-1
      jp1  =  j+1
      do 143   k  =  jp1, n
        s11  =  -one*a(j,k)
        s12  =  -one*a(k,j+1)
        if(j .eq. 1)  goto 142
        do 141  i  =  1, jm1
          s11  =  dotf(a(i,k),a(j,i),s11)
          s12  =  dotf(a(i,j+1),a(k,i),s12)
141     continue
142     a(j,k)   =  (-one*s11) * a(j,j)
        a(k,j+1) =   -one*dotf(a(j,j+1),a(k,j),s12)
143   continue
144 continue
150 if(mod(nxch,2) .ne. 0)  det  =  -one*det
    if(jfail .ne. jrange)   det  =  zero
    ir(n)  =  nxch
    return

  contains

  real(kind=fPrec) function dotf(x,y,s11)
    real(kind=fPrec), intent(in) :: x, y, s11
    dotf = x * y + s11
  end function dotf

  integer function ipairf(j,k)
    integer, intent(in) :: j, k
    ipairf = j*2**12 + k
  end function ipairf

  real(kind=fPrec) function pivotf(x)
    real(kind=fPrec), intent(in) :: x
    pivotf = abs(x)
  end function pivotf

  real(kind=fPrec) function sizef(x)
    real(kind=fPrec), intent(in) :: x
    sizef = abs(x)
  end function sizef

end subroutine rfact

subroutine dfact(n,a,idim,ir,ifail,det,jfail)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

  implicit none

  integer i,idim,ifail,imposs,ir,j,jfail,jm1,jover,jp1,jrange,junder,k,l,n,normal,nxch
  real(kind=fPrec) g1,g2,p,q,t
  real(kind=fPrec) a,det,s11,s12,x,y,tf
  character(len=6)         hname
  dimension ir(*),a(idim,*)

  !      data      g1, g2              /  1.e-37,  1.e37  /
  data      g1, g2              /  1.0e-37_fPrec,  1.0e37_fPrec  /
  data      hname               /  ' DFACT'  /
  data      normal, imposs      /  0, -1  /
  data      jrange, jover, junder  /  0, +1, -1  /
  save

  if(idim .ge. n  .and.  n .gt. 0) goto 110
  call tmprnt(hname,n,idim,0)
  return

110 continue
  ifail = normal
  jfail = jrange
  nxch  = 0
  det   = one

  do 144    j  =  1, n
120 k  =  j
  p  =  pivotf(a(j,j))
  if(j .eq. n)  goto 122
  jp1  =  j+1
  do 121    i  =  jp1, n
    q  =  pivotf(a(i,j))
    if(q .le. p)  goto 121
    k  =  i
    p  =  q
121 continue
    if(k .ne. j)  goto 123
122 if(p .gt. 0.)  goto 130
    det    =  zero
    ifail  =  imposs
    jfail  =  jrange
    return
123 do 124    l  =  1, n
      tf      =  a(j,l)
      a(j,l)  =  a(k,l)
      a(k,l)  =  tf
124 continue
    nxch      =  nxch + 1
    ir(nxch)  =  ipairf(j,k)
130 det     =  det * a(j,j)
    a(j,j)  =  one / a(j,j)
    t  =  sizef(det)
    if(t .lt. g1)  then
      det    =  zero
      if(jfail .eq. jrange)  jfail  =  junder
    elseif(t .gt. g2)  then
      det    =  one
      if(jfail .eq. jrange)  jfail  =  jover
    endif
    if(j .eq. n)  goto 144
    jm1  =  j-1
    jp1  =  j+1
    do 143   k  =  jp1, n
      s11  =  -one*a(j,k)
      s12  =  -one*a(k,j+1)
      if(j .eq. 1)  goto 142
      do 141  i  =  1, jm1
        s11  =  dotf(a(i,k),a(j,i),s11)
        s12  =  dotf(a(i,j+1),a(k,i),s12)
141   continue
142   a(j,k)    =  (-one*s11) * a(j,j)
      a(k,j+1)  =  -one*dotf(a(j,j+1),a(k,j),s12)
143 continue
144 continue
150 if(mod(nxch,2) .ne. 0)  det  =  -one*det
  if(jfail .ne. jrange)   det  =  zero
  ir(n)  =  nxch
  return

  contains

  real(kind=fPrec) function dotf(x,y,s11)
    real(kind=fPrec), intent(in) :: x, y, s11
    dotf = x * y + s11
  end function dotf

  integer function ipairf(j,k)
    integer, intent(in) :: j, k
    ipairf = j*2**12 + k
  end function ipairf

  real(kind=fPrec) function pivotf(x)
    real(kind=fPrec), intent(in) :: x
    pivotf = abs(x)
  end function pivotf

  real(kind=fPrec) function sizef(x)
    real(kind=fPrec), intent(in) :: x
    sizef = abs(x)
  end function sizef

end subroutine dfact

  ! This subroutine is not used
subroutine rfeqn(n,a,idim,ir,k,b)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

  implicit none

  integer i,idim,ij,im1,ir,j,k,l,m,n,nm1,nmi,nmjp1,nxch
  real(kind=fPrec) a,b,te,x,y
  real(kind=fPrec) s21,s22
  character(len=6) hname
  dimension ir(*),a(idim,*),b(idim,*)
  data      hname               /  ' RFEQN'  /
  save

  if(idim .ge. n  .and.  n .gt. 0  .and.  k .gt. 0) goto 210
  call tmprnt(hname,n,idim,k)
  return

210 continue
  nxch = ir(n)
  if(nxch .eq. 0) goto 220
  do m=1,nxch
    ij = ir(m)
    i  = ij / 4096
    j  = mod(ij,4096)
    do l=1,k
      te     = b(i,l)
      b(i,l) = b(j,l)
      b(j,l) = te
    end do
  end do

220 continue
  do l=1,k
    b(1,l) = a(1,1)*b(1,l)
  end do
  if(n .eq. 1) goto 299
  do l=1,k
    do i=2,n
      im1 = i-1
      s21 = -one*b(i,l)
      do j=1,im1
        s21 = dotf(a(i,j),b(j,l),s21)
      end do
      b(i,l) = (-one*a(i,i))*s21
    end do
    nm1 = n-1
    do i=1,nm1
      nmi = n-i
      s22 = -one*b(nmi,l)
      do j=1,i
        nmjp1 = n - j+1
        s22   = dotf(a(nmi,nmjp1),b(nmjp1,l),s22)
      end do
      b(nmi,l) = -one*s22
    end do
  end do

299 continue
  return

  contains

  real(kind=fPrec) function dotf(x,y,s21)
    real(kind=fPrec), intent(in) :: x, y, s21
    dotf = x * y + s21
  end function dotf

end subroutine rfeqn

subroutine rfinv(n,a,idim,ir)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

  implicit none

  integer i,idim,ij,im2,ir,j,k,m,n,nm1,nmi,nxch
  real(kind=fPrec) a,ti,x,y
  real(kind=fPrec) s31,s32,s33,s34
  character(len=6) hname
  dimension ir(*),a(idim,*)
  data      hname               /  ' RFINV'  /
  save

  if(idim .ge. n .and. n .gt. 0) goto 310
  call tmprnt(hname,n,idim,0)
  return
310 if(n .eq. 1) return

  a(2,1)  = (-one*a(2,2)) * dotf(a(1,1),a(2,1),zero)
  a(1,2)  =  -one*a(1,2)
  if(n .eq. 2) goto 330
  do i=3,n
    im2 = i-2
    do j=1,im2
      s31  =  zero
      s32  =  a(j,i)
      do k=j,im2
        s31  =  dotf(a(k,j),a(i,k),s31)
        s32  =  dotf(a(j,k+1),a(k+1,i),s32)
      end do
      a(i,j) = (-one*a(i,i)) * dotf(a(i-1,j),a(i,i-1),s31)
      a(j,i) =  -one*s32
    end do
    a(i,i-1) = (-one*a(i,i))*dotf(a(i-1,i-1),a(i,i-1),zero)
    a(i-1,i) =  -one*a(i-1,i)
  end do

330 continue
  nm1 = n-1
  do i=1,nm1
    nmi  =  n-i
    do j=1,i
      s33 = a(i,j)
      do k=1,nmi
        s33 = dotf(a(i+k,j),a(i,i+k),s33)
      end do
      a(i,j) = s33
    end do
    do j=1,nmi
      s34 = zero
      do k=j,nmi
        s34 = dotf(a(i+k,i+j),a(i,i+k),s34)
      end do
      a(i,i+j) = s34
    end do
  end do
  nxch = ir(n)
  if(nxch .eq. 0)  return
  do m=1,nxch
    k   =  nxch - m+1
    ij  =  ir(k)
    i   =  ij / 4096
    j   =  mod(ij,4096)
    do k=1,n
      ti     = a(k,i)
      a(k,i) = a(k,j)
      a(k,j) = ti
    end do
  end do

  return

  contains

  real(kind=fPrec) function dotf(x,y,s31)
    real(kind=fPrec), intent(in) :: x, y, s31
    dotf = x * y + s31
  end function dotf

end subroutine rfinv

subroutine dfinv(n,a,idim,ir)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

  implicit none

  integer i,idim,ij,im2,ir,j,k,m,n,nm1,nmi,nxch
  real(kind=fPrec) a,s31,s32,s33,s34,ti,x,y
  character(len=6) hname
  dimension ir(*),a(idim,*)
  data      hname               /  ' DFINV'  /
  save

  if(idim .ge. n  .and.  n .gt. 0) goto 310
  call tmprnt(hname,n,idim,0)
  return
310 if(n .eq. 1) return

  a(2,1) = (-one*a(2,2)) * dotf(a(1,1),a(2,1),zero)
  a(1,2) = -one*a(1,2)
  if(n .eq. 2) goto 330
  do i=3,n
    im2 = i-2
    do j=1,im2
      s31 = zero
      s32 = a(j,i)
      do k=j,im2
        s31 = dotf(a(k,j),a(i,k),s31)
        s32 = dotf(a(j,k+1),a(k+1,i),s32)
      end do
      a(i,j) = (-one*a(i,i)) * dotf(a(i-1,j),a(i,i-1),s31)
      a(j,i) = -one*s32
    end do
    a(i,i-1) = (-one*a(i,i)) * dotf(a(i-1,i-1),a(i,i-1),zero)
    a(i-1,i) = -one*a(i-1,i)
  end do

330 continue
  nm1 = n-1
  do i=1,nm1
    nmi = n-i
    do j=1,i
      s33 = a(i,j)
      do k=1,nmi
        s33 = dotf(a(i+k,j),a(i,i+k),s33)
      end do
      a(i,j)  =  s33
    end do
    do j=1,nmi
      s34 = zero
      do k = j,nmi
        s34 = dotf(a(i+k,i+j),a(i,i+k),s34)
      end do
      a(i,i+j) = s34
    end do
  end do
  nxch = ir(n)
  if(nxch .eq. 0) return
  do m=1,nxch
    k   =  nxch - m+1
    ij  =  ir(k)
    i   =  ij / 4096
    j   =  mod(ij,4096)
    do k=1,n
      ti      =  a(k,i)
      a(k,i)  =  a(k,j)
      a(k,j)  =  ti
    end do
  end do

  return

  contains

  real(kind=fPrec) function dotf(x,y,s31)
    real(kind=fPrec), intent(in) :: x, y, s31
    dotf = x * y + s31
  end function dotf

end subroutine dfinv

subroutine tmprnt(name,n,idim,k)

  use floatPrecision
  use mathlib_bouncer
  use crcoall

  implicit none

  integer idim,k,lgfile,n
  character(len=6) name
  logical mflag,rflag
  save

  if(name(2:2) .eq. 'S') then
    call kermtr('F012.1',lgfile,mflag,rflag)
  else
    call kermtr('F011.1',lgfile,mflag,rflag)
  end if
  if(mflag) then
    if(lgfile .eq. 0) then
      if(name(3:6) .eq. 'FEQN') then
        write(lout,1002) name, n, idim, k
      else
        write(lout,1001) name, n, idim
      endif
    else
      if(name(3:6) .eq. 'FEQN') then
        write(lgfile,1002) name, n, idim, k
      else
        write(lgfile,1001) name, n, idim
      endif
    endif
  endif
  if(.not. rflag) then
    write(lerr,"(a)") "KERNLIB> ERROR TMPRNT: "//name
    call prror
  endif
  return
1001 format(7x,' parameter error in subroutine ',a6,' (n.lt.1 or idim.lt.n).',5x,'n =',i4,5x,'idim =',i4,'.')
1002 format(7x,' parameter error in subroutine ',a6,' (n.lt.1 or idim.lt.n or k.lt.1).',5x,'n =',i4,5x,'idim =',i4,5x,'k =',i4,'.')
end subroutine tmprnt

! This subrtoutine is not used
subroutine dfeqn(n,a,idim,ir,k,b)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

  implicit none

  integer i,idim,ij,im1,ir,j,k,l,m,n,nm1,nmi,nmjp1,nxch
  real(kind=fPrec) a,b,x,y,te
  real(kind=fPrec) s21,s22
  character(len=6) hname
  dimension ir(*),a(idim,*),b(idim,*)
  data      hname               /  ' DFEQN'  /
  save

  if(idim .ge. n  .and.  n .gt. 0  .and.  k .gt. 0) goto 210
  call tmprnt(hname,n,idim,k)
  return

210 continue
  nxch = ir(n)
  if(nxch .eq. 0) goto 220
  do m=1,nxch
    ij = ir(m)
    i  = ij / 4096
    j  = mod(ij,4096)
    do l=1,k
      te     = b(i,l)
      b(i,l) = b(j,l)
      b(j,l) = te
    end do
  end do
220 continue
  do l=1,k
    b(1,l) = a(1,1)*b(1,l)
  end do
  if(n .eq. 1) goto 299
  do l=1,k
    do i=2,n
      im1 = i-1
      s21 = - b(i,l)
      do j=1,im1
        s21 = dotf(a(i,j),b(j,l),s21)
      end do
      b(i,l) = (-one*a(i,i))*s21
    end do
    nm1 = n-1
    do i=1,nm1
      nmi = n-i
      s22 = -one*b(nmi,l)
      do j=1,i
        nmjp1 = n - j+1
        s22   = dotf(a(nmi,nmjp1),b(nmjp1,l),s22)
      end do
      b(nmi,l) = -one*s22
    end do
  end do

299 continue
  return

  contains

  real(kind=fPrec) function dotf(x,y,s21)
    real(kind=fPrec), intent(in) :: x, y, s21
    dotf = x * y + s21
  end function dotf

end subroutine dfeqn

end module matrix_inv
