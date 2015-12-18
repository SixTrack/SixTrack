+cd crlibco
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
+cd crcoall
      integer lout
      common /crflags/lout
+cd dascr
      integer idao,is,iscrri
      double precision rs
      common/dascr/is(100),rs(100),iscrri(100),idao
+cd alloc
      logical allvec(lda)
      common /alloc/ allvec
+cd alloctot
      integer ndat
      common /alloctot/ ndat
+cd daname
      character daname(lda)*10
      common / daname / daname
+cd dabinc
      integer i1,i2,ia1,ia2,idall,idalm,idano,idanv,idapo,ie1,ie2,ieo,  &
     &ifi,lfi,nda,ndamaxi,nmmax,nocut,nomax,nst,nvmax
      double precision cc,eps,epsmac,facint
+if small
      parameter(lda=10000,lst=200000,lea=500,                           &
     &lia=10000,lno=120,lnv=40)
      common /fordes/ nda,ndamaxi
      common / da / cc(lst),eps,epsmac
      common / dai / i1(lst),i2(lst),                                   &
     &ie1(lea),ie2(lea),ieo(lea),ia1(0:lia),ia2(0:lia),ifi(lea),        &
     &idano(lda),idanv(lda),idapo(lda),idalm(lda),idall(lda),           &
     &nst,nomax,nvmax,nmmax,nocut,lfi
      common /factor/ facint(0:lno)
+ei
+if big
      parameter(lda=10000,lst=20050000,lea=100000,                      &
     &lia=5000000,lno=120,lnv=40)
      common /fordes/ nda,ndamaxi
      common / da / cc(lst),eps,epsmac
      common / dai / i1(lst),i2(lst),                                   &
     &ie1(lea),ie2(lea),ieo(lea),ia1(0:lia),ia2(0:lia),ifi(lea),        &
     &idano(lda),idanv(lda),idapo(lda),idalm(lda),idall(lda),           &
     &nst,nomax,nvmax,nmmax,nocut,lfi
      common /factor/ facint(0:lno)
+ei
+if ctrack
      parameter(lda=1000,lst=3000000,lea=20000,                         &
     &lia=20000,lno=120,lnv=40)
      common /fordes/ nda,ndamaxi
      common / da / cc(lst),eps,epsmac
      common / dai / i1(lst),i2(lst),                                   &
     &ie1(lea),ie2(lea),ieo(lea),ia1(0:lia),ia2(0:lia),ifi(lea),        &
     &idano(lda),idanv(lda),idapo(lda),idalm(lda),idall(lda),           &
     &nst,nomax,nvmax,nmmax,nocut,lfi
      common /factor/ facint(0:lno)
+ei
+dk daini
!******************************************************************************
!                                                                             *
!                                                                             *
!                                                                             *
!               DIFFERENTIAL ALGEBRA PACKAGE OF M. BERZ                       *
!                     ****************************                            *
!                                                                             *
!                                                                             *
!                                                                             *
!                                                                             *
!         VERSION FOR MACHINE IN LINE THAT IS NOT COMMENTED OFF               *
!        TO CREATE DIFFERENT VERSIONS, USE THE PROGRAM 'VERSION'              *
!                                                                             *
!                                                                             *
!                                                                             *
!                                                                             *
!        THIS PACKAGE WAS INITIALLY WRITTEN BY PROF. M. BERZ WHILE AT         *
!        THE LAWRENCE BERKELEY LABORATORY.                                    *
!        IT HAS BEEN EXTENSIVELY MODIFIED BY THE MEMBERS OF THE ESG GROUP.    *
!        THEREFORE PROF. BERZ SHOULD NOT BE HELD RESPONSIBLE FOR ANY BUGS.    *
!                                                                             *
!                  NEW RULES OF THE GAME (EXHAUSTIVE)                         *
!                 **********************************                          *
!                         THERE ARE NONE                                      *
!                                                                             *
!******************************************************************************
!
!
!     THIS FILE CONTAINS ROUTINES TO PERFORM DIFFERENTIAL ALGEBRA (DA)
!     AS AN OPTION, ALSO COMPONENTWISE ALGEBRA (CA) CAN BE PERFORMED.
!     A DESCRIPTION OF THE INTERNAL ARRAYS USED BY THE ROUTINES CAN
!     BE FOUND IN BLOCKDATA DABLD.
!
!
!     SHORT REFERENCE CHART
!     *********************
!
!     THE PARAMETERS USED BELOW HAVE THE FOLLOWING MEANING:
!
!     A,B:                NAME OF INPUT DA VECTORS   (INTEGER)
!     C:                  NAME OF OUTPUT DA VECTOR   (INTEGER)
!     X,Y:                NAME OF INPUT DA MATRIX    (INTEGER(...))
!     Z:                  NAME OF OUTPUT DA MATRIX   (INTEGER(...))
!
!     F:                  NAME OF A DA FUNCTION      (CHARACTER*4)
!     G:                  NAME OF EXTERNAL FUNCTION  (DOUBLE PRECISION)
!     JJ:                 ARRAY OF EXPONENTS         (INTEGER(20))
!     O:                  ORDER                      (INTEGER)
!     N:                  NUMBER OF VARIABLES        (INTEGER)
!     I,J,K:              INTEGER NUMBER             (INTEGER
!     R,RA,RB:            REAL NUMBERS               (DOUBLE PRECISION)
!     H:                  ARRAY OF LENGTH LH         (DOUBLE PRECISION)
!     U:                  OUTPUT UNIT NUMBER         (INTEGER)
!     T:                  COMMENT TEXT               (CHARACTER*10)
!
!
!               SUBROUTINES AND THEIR CALLING PARAMETERS
!               ****************************************
!
!     DAINI(O,N,U):       INITIALIZES CONTROL ARRAYS AND SETS MAX. ORDER O AND
!                         MAX. NUMBER OF VARIABLES N. MUST BE CALLED BEFORE ANY
!                         OTHER DA ROUTINE CAN BE USED.
!
!     DAALL(A,I,T,O,N):   ALLOCATES SPACE FOR I VECTORS A. T: CHARACTER NAME
!     DADAL(A,I):         DEALLOCATES THE I VECTORS A.
!!     DAVAR(A,R,I):       MAKES A INDEPENDENT VARIABLE # I WITH INITIAL VALUE R
!!     DACON(A,R):         SETS A TO CONSTANT R
!     DANOT(O):           SETS NEW TRUNCATION ORDER O FOR DA OPERATIONS
!     DAEPS(R):           SETS NEW ZERO TOLERANCE EPSILON
!
!!     DAPEK(A,JJ,R):      RETURNS COEF R OF MONOMIAL WITH EXPONENTS JJ OF A
!!     DAPOK(A,JJ,R):      SETS COEF OF MONOMIAL WITH EXPONENTS JJ OF A TO R
!
!!     DACOP(A,C):         PERFORMS C = A
!!     DAADD(A,B,C):       PERFORMS C = A + B
!!    DASUB(A,B,C):       PERFORMS C = A - B
!!     DAMUL(A,B,C):       PERFORMS C = A * B
!!     DADIV(A,B,C):       PERFORMS C = A / B
!!     DASQR(A,C):         PERFORMS C = A^2           (SQUARE OF A)
!
!!     DACAD(A,RA,C):      PERFORMS C = A + RA
!!     DACSU(A,RA,C):      PERFORMS C = A - RA
!!     DASUC(A,RA,C):      PERFORMS C = RA - A
!!     DACMU(A,RA,C):      PERFORMS C = A * RA
!!    DACDI(A,RA,C):      PERFORMS C = A / RA
!!     DADIC(A,RA,C):      PERFORMS C = RA / A
!!     DACMA(A,B,RB,C):    PERFORMS C = A + RB*B
!!DAMULIN(A,B,RA,C,D,RB,C):    PERFORMS C = A*B*RA + C*D*RB
!!     DALIN(A,RA,B,RB,C): PERFORMS C = A*RA + B*RB
!!     DAFUN(F,A,C):       PERFORMS C = F(A)          (DA FUNCTION)
!
!!     DAABS(A,R):         PERFORMS R = |A|           (NORM OF A)
!!     DACOM(A,B,R):       PERFORMS R = |A-B|         (NORM OF A-B)
!!     DAPOS(A,C):         PERFORMS C(I) = |A(I)|     (MAKE SIGNS POSITIVE)
!
!!     DACCT(X,I,Y,J,Z,K)  CONCATENATES Z = X o Y;   I,J,K: # OF VECTORS IN X,Y,
!!     DAINV(X,I,Z,K)      INVERTS Z = X^-1;           I,J: # OF VECTORS IN X,Y
!!     DAPIN(X,I,Z,K,JJ)   PARTIALLY INVERTS Z = X^-1; I,J: # OF VECTORS IN X,Y,
!                         JJ: ARRAY; NONZERO ENTRIES DENOTE TO BE INVERTED LINES
!
!!     DADER(I,A,C):       PERFORMS C = DA/DI (DERIV. WITH RESPECT TO VARIABLE I
!!     DAPOI(A,B,C,I):     PERFORMS C = [A,B] (POISSON BRACKET, 2*I: # PHASEVARS
!!     DACFU(A,G,C):       MULTIPLIES COEFFICIENTS WITH FUNCTION G(JJ)
!
!     DAPRI(A,U):         PRINTS DA VECTOR A TO UNIT U
!     DAREA(A,U):         READS DA VECTOR A FROM UNIT U
!     DADEB(U,T,I):       DEBUGGER, DUMPS TO U. T: MEMO, I=0: RETURN, I=1:STOP
!!     DARAN(A,R,seed):         FILLS A WITH RANDOM NUMBERS. R: FILLFACTOR
!     DANUM(O,N,I):       COMPUTES NUMBER OF MONOMIALS IN N VAR THROUGH ORDER O
!
!
!     ADDITIONAL ROUTINES THE USER DOES NOT NEED TO CALL:
!
!     DAINF: RETURNS INFOS ABOUT A DA VECTOR PREVIOUSLY DECLARED
!     DAPAC: PACKS DA VECTORS
!     DACHK: CHECKS IF DA VECTORS HAVE COMPATIBLE ATTRIBUTES
!     DCODE: TRANSFORMS DIGITS IN A CERTAIN BASE TO A DECIMAL INTEGER
!     NCODE: EXTRACTS DIGITS IN A CERTAIN BASE FROM A DECIMAL INTEGER
!
!
!     FURTHER WISHES
!     **************
!
!     - CHECK DAREA AND DAPRI FOR CA VECTORS
!     - MAKE DAFUN USE DASQR
!
!
!      BLOCKDATA DABLD
!     ***************
!
!
!     PARAMETERS:
!
!     LDA: MAXIMUM NUMBER OF DA-VECTORS;    CAN BE CHANGED QUITE ARBITRARILY
!     LST: LENGTH OF MAIN STORAGE STACK;    CAN BE CHANGED QUITE ARBITRARILY
!     LEA: MAXIMUM NUMBER OF MONOMIALS;     CAN BE INCREASED FOR LARGE NO,NV
!     LIA: DIMENSION OF IA1,IA2;            CAN BE INCREASED FOR LARGE NO,NV
!     LNO: MAXIMUM ORDER;                   CAN BE INCREASED TO ABOUT 1000
!     LNV: MAXIMUM NUMBER OF VARIABLES;     CAN BE INCREASED TO ABOUT 1000
!
!     ALL THE CHANGES IN THE VALUES OF PARAMETERS HAVE TO BE MADE BY GLOBAL
!     SUBSTITUTIONS IN ALL SUBROUTINES.
!
!     DANAME:   NAME OF DA VECTOR
!
!     CC:       STACK OF DOUBLE PRECISON COEFFICIENTS
!     I1:       FIRST CHARACTERISTIC INTEGER (CF DAINI)
!     I2:       SECOND CHARACTERISTIC INTEGER (CF DAINI)
!
!     IE1:      CHARACTERISTIC INTEGER 1 OF UNPACKED REPRESENTATION (CF DAINI)
!     IE2:      CHARACTERISTIC INTEGER 2 OF UNPACKED REPRESENTATION (CF DAINI)
!     IEO:      ORDER OF ENTRY IN UNPACKED REPRESENTATION
!     IA1:      REVERSE TO IE1 (CF DAINI)
!     IA2:      REVERSE TO IE2 (CF DAINI)
!
!     IDANO:    ORDER OF DA VECTOR; IN CA, NUMBER OF COMPONENTS
!     IDANV:    NUMBER OF VARIABLES; IF 0, INDICATES CA VECTOR
!     IDAPO:    FIRST ADDRESS IN STACK
!     IDALM:    NUMBER OF RESERVED STACK POSITIONS
!     IDALL:    NUMBER OF MOMENTARILY REQUIRED STACK POSITIONS
!
!     NDA:      NUMBER OF DA VECTORS MOMENTARILY DEFINED
!     NST:      NUMBER OF STACK POSITIONS MOMENTARILY ALLOCATED
!     NOMAX:    MAXIMUM REQUESTED ORDER  (CF DAINI)
!     NVMAX:    MAXIMUM REQUESTED NUMBER OF VARIABLES (CF DAINI)
!     NMMAX:    MAXIMUM NUMBER OF MONOMIALS FOR NOMAX, NVMAX (CF DAINI)
!     NOCUT:    MOMENTARY TRUNCATION ORDER
!     EPS:      TRUNCATION ACCURACY (CAN BE SET BY USER)
!     EPSMAC:   MANTISSA LENGTH OF MACHINE (PESSIMISTIC ESTIMATE)
!
!-----------------------------------------------------------------------------1
!
      subroutine daini(no,nv,iunit)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,iall,ibase,ic1,ic2,icmax,io1,io2,iout,iunit,j,jd,jj,jjj,&
     &jjjj,jl,js,k,lda,lea,lia,lno,lnv,lst,n,nn,no,nv
!     *****************************
!
!     THIS SUBROUTINE SETS UP THE MAJOR ORDERING AND ADDRESSING ARRAYS IN
!     COMMON BLOCK DAINI. IF IUNIT > 0, THE ARRAYS WILL BE PRINTED TO UNIT
!     NUMBER IUNIT. AN EXAMPLE FOR THE ARRAYS GENERATED BY DAINI CAN BE
!     FOUND AFTER THE ROUTINE.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!-----------------------------------------------------------------------------9
!      COMMON / DASCR /  IS(20), RS(20)                                        1
+ca dascr
!-----------------------------------------------------------------------------2

+ca alloc

!
      character aa*10
      dimension n(lnv+1),k(0:lnv),j(lnv),jj(lnv)
+if boinc
      character*256 filename
+ei
!
      if(eps.le.0.d0) eps=1.d-38
!      if(EPS.le.0.d0) eps=1.d-90
      epsmac=1.d-7
      if(nv.eq.0) return
      ndamaxi=0
!
      do 5 i=1, lda
        allvec(i) = .false.
5     continue

!*****************************************

!     INITIALIZING VARIABLES IN COMMON / DA /
!     ***************************************
!
      nda   = 0
      nst   = 0
      nomax = no
      nvmax = nv
      call danum(no,nv,nmmax)
      nocut = no
      lfi   = 0
!
      do 10 i=0,lia
      ia1(i) = 0
      ia2(i) = 0
  10  continue
!
! Eric
      do i=1,lst
        i1(i)=0
        i2(i)=0
      enddo

      idao=0 
      do i=1,100
        is(i) = 0
        iscrri(i)=0
        rs(i)=0d0
      enddo
!
      if(nv.gt.lnv.or.no.gt.lno) then
+if cr
         write(lout,*)'ERROR IN SUBROUTINE DAINI, NO, NV = ',no,nv
+ei
+if .not.cr
         write(*,*)'ERROR IN SUBROUTINE DAINI, NO, NV = ',no,nv
+ei
         call dadeb(31,'ERR DAINI ',1)
      endif
!
      ibase = no+1
      js    = nv/2
!hr10 if(float(ibase)**((nv+1)/2).gt.float(lia)) then
      if(real(ibase)**((nv+1)/2).gt.real(lia)) then                      !hr10
+if cr
         write(lout,*)'ERROR, NO = ',no,', NV = ',nv,' TOO LARGE FOR',
+ei
+if .not.cr
         write(*,*)'ERROR, NO = ',no,', NV = ',nv,' TOO LARGE FOR',     &
+ei
     &' LIA = ',lia
         call dadeb(31,'ERR DAINI ',1)
      endif
!
      icmax = 0
      nn    = 0
      k(0)  = 0
!
      do 100 io2=0,no
!     ***************
!
      n(1)  = io2
      jl    = 0
      jd    = 1
!
  50  jl    = jl + jd
!

!old
!      IF(JL.EQ.0) THEN
!old
!     modified according to Wu Ying
!
      if(jl.le.0) then
!
         goto 100
      elseif(jd.eq.1) then
         j(jl) = 0
      else
         j(jl) = j(jl) + 1
      endif
!
      k(jl)    = k(jl-1)*ibase + j(jl)
      n(jl+1)  = n(jl) - j(jl)
!
      if(j(jl).gt.n(jl)) then
         jd    = -1
         goto 50
      elseif(jl.lt.js) then
         jd    = 1
         goto 50
      else
         j(jl) = n(jl)
         k(jl) = k(jl-1)*ibase + j(jl)

         ic2   = k(jl)
         icmax = max(icmax,ic2)
         k(jl) = 0
!
         ia2(ic2) = nn
!
         do 80 io1=0,no-io2
!        ******************
!
         n(js+1) = io1
         jd      = 1
!
  70     jl      = jl + jd
!
         if(jl.eq.js) then
            goto 80
         elseif(jd.eq.1) then
            j(jl) = 0
         else
            j(jl) = j(jl) + 1
         endif
!
         k(jl)    = k(jl-1)*ibase + j(jl)
         n(jl+1)  = n(jl) - j(jl)
!
         if(j(jl).gt.n(jl)) then
            jd    = -1
            goto 70
         elseif(jl.lt.nv) then
            jd    = 1
            goto 70
         else
            jd    = -1
            j(jl) = n(jl)
            k(jl) = k(jl-1)*ibase + j(jl)
            ic1   = k(jl)
            icmax = max(icmax,ic1)
            nn = nn + 1
!
            ie1(nn) = ic1
            ie2(nn) = ic2
            i1 (nn) = ic1
            i2 (nn) = ic2
            if(ic2.eq.0) ia1(ic1) = nn
            ieo(nn) = io1 + io2
!
            goto 70
         endif
!
   80    continue
!
         jd = -1
         goto 50
      endif
!
  100 continue
  110 continue
!
      if(nn.gt.lea) then
+if cr
         write(lout,*)'ERROR IN DAINI, NN = ',nn,' EXCEEDS LEA'
+ei
+if .not.cr
         write(*,*)'ERROR IN DAINI, NN = ',nn,' EXCEEDS LEA'
+ei
         call dadeb(31,'ERR DAINI ',1)
      endif
!
!     ALLOCATING SCRATCH VARIABLES
!     ****************************
!
      iall = 0
      call daall(iall,1,'$$UNPACK$$',nomax,nvmax)
!
      do 150 i=0,nomax
      aa = '$$MUL   $$'
      write(aa(6:10),'(I5)') i
      iall = 0
!      CALL DAALL(IALL,1,AA,I,NVMAX)
      call daall(iall,1,aa,nomax,nvmax)
 150  continue
!
      idall(1) = nmmax
!
!     DOUBLE CHECKING ARRAYS IE1,IE2,IA1,IA2
!     **************************************
!
      do 300 i=1,nmmax
!
      jjj = ia1(ie1(i)) + ia2(ie2(i))
      if(jjj.ne.i) then
+if cr
         write(lout,*)                                                  &
     &'ERROR IN DAINI IN ARRAYS IE1,IE2,IA1,IA2 AT I = ',i
+ei
+if .not.cr
         write(*,*)                                                     &
     &'ERROR IN DAINI IN ARRAYS IE1,IE2,IA1,IA2 AT I = ',i
+ei
         call dadeb(31,'ERR DAINI ',1)
      endif
!
 300  continue
!
      if(iunit.eq.0) return
!
+if cr
      write(lout,*)'ARRAY SETUP DONE, BEGIN PRINTING'
+ei
+if .not.cr
      write(*,*)'ARRAY SETUP DONE, BEGIN PRINTING'
+ei
!
      iout = 32
+if boinc
      call boincrf('DAINI.DAT',filename)
+if fio
      open(iout,file=filename,status='new',round='nearest')
+ei
+if .not.fio
      open(iout,file=filename,status='new')
+ei
+ei
+if .not.boinc
+if fio
      open(iout,file='DAINI.DAT',status='NEW',round='nearest')
+ei
+if .not.fio
      open(iout,file='DAINI.DAT',status='NEW')
+ei
+ei
!CRAY OPEN(IOUT,FILE='DAINI',STATUS='UNKNOWN',FORM='FORMATTED')          *CRAY
!CRAY REWIND IOUT                                                        *CRAY
!
      write(iout,'(/A/A/)') ' ARRAYS I1 THROUGH I20, IE1,IE2,IEO',      &
     &' **********************************'
      do 200 i=1,nmmax
      call dancd(ie1(i),ie2(i),jj)
      write(iout,'(1X,I5,2X,4(5I2,1X),3I6)') i,(jj(jjjj),jjjj=1,lnv),   &
     &ie1(i),ie2(i),ieo(i)
 200  continue
!
      write(iout,'(/A/A/)') ' ARRAYS IA1,IA2',' **************'
      do 210 i=0,icmax
      write(iout,'(3I10)') i,ia1(i),ia2(i)
 210  continue
!
      return
      end
+dk daexter
      subroutine daexter
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,lda,lea,lia,lno,lnv,lst
!     *****************************
!
!-----------------------------------------------------------------------------1
+ca dabinc
!-----------------------------------------------------------------------------9
+ca alloc

!
      do 5 i=1, lda
        allvec(i)=.false.
5     continue

      return
      end

+dk dallsta
      subroutine dallsta(ldanow)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,lda,ldanow,lea,lia,lno,lnv,lst
!     *****************************
!
!-----------------------------------------------------------------------------1
+ca dabinc
!-----------------------------------------------------------------------------9
+ca alloc

!
      ldanow=0
      do 5 i=1, lda
        if(allvec(i)) ldanow=ldanow+1
5     continue

+if cr
      write(lout,*) ' ALLOCATED ',ldanow
+ei
+if .not.cr
      write(*,*) ' ALLOCATED ',ldanow
+ei

      return
      end

!
! EXAMPLE: ARRAYS I1 THROUGH I20, IE1,IE2,IEO (NOMAX=3,NVMAX=4)
! *************************************************************
!     I   I1               THROUGH               I20     IE1   IE2   IEO
!     1   0 0 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      0     0     0
!     2   1 0 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      1     0     1
!     3   0 1 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      4     0     1
!     4   2 0 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      2     0     2
!     5   1 1 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      5     0     2
!     6   0 2 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      8     0     2
!     7   3 0 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      3     0     3
!     8   2 1 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      6     0     3
!     9   1 2 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      9     0     3
!    10   0 3 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0     12     0     3
!    11   0 0 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      0     1     1
!    12   1 0 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      1     1     2
!    13   0 1 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      4     1     2
!    14   2 0 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      2     1     3
!    15   1 1 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      5     1     3
!    16   0 2 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      8     1     3
!    17   0 0 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      0     4     1
!    18   1 0 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      1     4     2
!    19   0 1 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      4     4     2
!    20   2 0 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      2     4     3
!    21   1 1 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      5     4     3
!    22   0 2 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      8     4     3
!    23   0 0 0 0 0  0 0 0 0 0  2 0 0 0 0  0 0 0 0 0      0     2     2
!    24   1 0 0 0 0  0 0 0 0 0  2 0 0 0 0  0 0 0 0 0      1     2     3
!    25   0 1 0 0 0  0 0 0 0 0  2 0 0 0 0  0 0 0 0 0      4     2     3
!    26   0 0 0 0 0  0 0 0 0 0  1 1 0 0 0  0 0 0 0 0      0     5     2
!    27   1 0 0 0 0  0 0 0 0 0  1 1 0 0 0  0 0 0 0 0      1     5     3
!    28   0 1 0 0 0  0 0 0 0 0  1 1 0 0 0  0 0 0 0 0      4     5     3
!    29   0 0 0 0 0  0 0 0 0 0  0 2 0 0 0  0 0 0 0 0      0     8     2
!    30   1 0 0 0 0  0 0 0 0 0  0 2 0 0 0  0 0 0 0 0      1     8     3
!    31   0 1 0 0 0  0 0 0 0 0  0 2 0 0 0  0 0 0 0 0      4     8     3
!    32   0 0 0 0 0  0 0 0 0 0  3 0 0 0 0  0 0 0 0 0      0     3     3
!    33   0 0 0 0 0  0 0 0 0 0  2 1 0 0 0  0 0 0 0 0      0     6     3
!    34   0 0 0 0 0  0 0 0 0 0  1 2 0 0 0  0 0 0 0 0      0     9     3
!    35   0 0 0 0 0  0 0 0 0 0  0 3 0 0 0  0 0 0 0 0      0    12     3
!
!    ARRAYS IA1,IA2
!    **************
!    I        IA1       IA2
!    0         1         0   IE1,IE2 AND IA1,IA2 ALLOW THE EASY COMPUTATION
!    1         2        10   OF THE ADDRESS OF THE PRODUCT OF TWO MONOMIALS.
!    2         4        22   LET IX AND IY BE THE POSITIONS OF THE TWO
!    3         7        31   FACTORS. THEN THE POSITION IZ OF THE PRODUCT OF
!    4         3        16   THE TWO FACTORS IS GIVEN BY
!    5         5        25
!    6         8        32   IZ = IA1(IE1(IX)+IE1(IY)) + IA2(IE2(IX)+IE2(IY))
!    7         0         0
!    8         6        28
!    9         9        33   THE OTHER VARIABLES SET BY DAINI WOULD HAVE THE
!   10         0         0   VALUES
!   11         0         0
!   12        10        34   NOMAX = 3,  NVMAX = 4, NMMAX = 35
!
+dk daallno
      subroutine daallno(ic,l,ccc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ind,l,lda,lea,lia,lno,lnv,lst,ndanum,no,nv
      double precision x
!     ********************************
!
!     THIS SUBROUTINE ALLOCATES STORAGE FOR A DA VECTOR WITH
!     ORDER NOmax AND NUMBER OF VARIABLES NVmax
!
!-----------------------------------------------------------------------------1
+ca dabinc
+ca alloc

!-----------------------------------------------------------------------------9
+ca daname
!-----------------------------------------------------------------------------3
!
      integer ic(*)
      logical incnda
      character c*10,ccc*10
!
      no=nomax
      nv=nvmax
      ind = 1
      do 10 i=1,l
        if(ic(i).gt.0.and.ic(i).le.nda) then
!         DANAME(IC(I)) = C
!         IF(IDANO(IC(I)).EQ.NO.AND.IDANV(IC(I)).EQ.NV) THEN
        else
          if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
+if cr
             write(lout,*)'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',
+ei
+if .not.cr
             write(*,*)'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',    &
+ei
     &no,nv,' NOMAX, NVMAX = ',nomax,nvmax
            call dadeb(31,'ERR DAALL ',1)
          endif
!
20        if (allvec(ind)) then
           ind = ind + 1
           goto 20
          endif

          incnda = .false.
          if (ind .gt. nda) then
             incnda = .true.
             nda = nda + 1
             if(nda.gt.lda) then
+if cr
         write(lout,*)                                                  &
     &'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
+ei
+if .not.cr
         write(*,*)                                                     &
     &'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
+ei
             call dadeb(31,'ERR DAALL ',1)
             endif
          endif

          allvec(ind) = .true.

          ic(i) = ind
!
          if(nv.ne.0) then
            call danum(no,nv,ndanum)
          else
            ndanum = no
          endif

          c = ccc
          if(l.ne.1) write(c(6:10),'(I5)') i
          daname(ind) = c

           if (incnda) then
            if(ind.gt.nomax+2) then
             idano(ind) = nomax
             idanv(ind) = nvmax
             idapo(ind) = nst + 1
             idalm(ind) = nmmax
             idall(ind) = 0
             nst = nst + nmmax
            else
             idano(ind) = no
             idanv(ind) = nv
             idapo(ind) = nst + 1
             idalm(ind) = ndanum
             idall(ind) = 0
             nst = nst + ndanum
            endif
          endif
!
          if(nst.gt.lst) then
            x=-1.d0
+if cr
            write(lout,*)'ERROR IN DAALL, STACK EXHAUSTED '
+ei
+if .not.cr
            write(*,*)'ERROR IN DAALL, STACK EXHAUSTED '
+ei
+if cr
            write(lout,*) ' NST,LST '
+ei
+if .not.cr
            write(*,*) ' NST,LST '
+ei
+if cr
            write(lout,*)  nst,lst
+ei
+if .not.cr
            write(*,*)  nst,lst
+ei
+if cr
            write(lout,*) ' NDA,NDANUM,NDA*NDANUM '
+ei
+if .not.cr
            write(*,*) ' NDA,NDANUM,NDA*NDANUM '
+ei
+if cr
            write(lout,*)  nda,ndanum,nda*ndanum
+ei
+if .not.cr
            write(*,*)  nda,ndanum,nda*ndanum
+ei
!            X=DSQRT(X)
            call dadeb(31,'ERR DAALL ',1)
          endif
!
          if(nv.eq.0.or.nomax.eq.1) then
            call daclr(ic(i))
            idall(ic(i)) = idalm(ic(i))
          endif
        endif
  10  continue
!
      if(nda.gt.ndamaxi) ndamaxi=nda

      return
      end
+dk daall
      subroutine daall(ic,l,ccc,no,nv)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ind,l,lda,lea,lia,lno,lnv,lst,ndanum,no,nv
      double precision x
!     ********************************
!
!     THIS SUBROUTINE ALLOCATES STORAGE FOR A DA VECTOR WITH
!     ORDER NO AND NUMBER OF VARIABLES NV
!
!-----------------------------------------------------------------------------1
+ca dabinc
+ca alloc
+ca alloctot

!-----------------------------------------------------------------------------9
+ca daname
!-----------------------------------------------------------------------------3
!
      integer ic(*)
      logical incnda
      character c*10,ccc*10
!
      ind = 1

      do 10 i=1,l
        if(ic(i).gt.0.and.ic(i).le.nda) then
!         DANAME(IC(I)) = C
!         IF(IDANO(IC(I)).EQ.NO.AND.IDANV(IC(I)).EQ.NV) THEN
        else
          if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
+if cr
             write(lout,*)'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',
+ei
+if .not.cr
             write(*,*)'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',    &
+ei
     &no,nv,' NOMAX, NVMAX = ',nomax,nvmax
            call dadeb(31,'ERR DAALL ',1)
          endif
!
20        if (allvec(ind)) then
           ind = ind + 1
           goto 20
          endif

          incnda = .false.
          if (ind .gt. nda) then
             incnda = .true.
             nda = nda + 1
             ndat=nda
             if(nda.gt.lda) then
+if cr
        write(lout,*)                                                   &
     &'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
+ei
+if .not.cr
        write(*,*)                                                      &
     &'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
+ei
             call dadeb(31,'ERR DAALL ',1)
             endif
          endif

          allvec(ind) = .true.

          ic(i) = ind
!
          if(nv.ne.0) then
            call danum(no,nv,ndanum)
          else
            ndanum = no
          endif

          c = ccc
          if(l.ne.1) write(c(6:10),'(I5)') i

          daname(ind) = c

           if (incnda) then
            if(ind.gt.nomax+2) then
             idano(ind) = nomax
             idanv(ind) = nvmax
             idapo(ind) = nst + 1
             idalm(ind) = nmmax
             idall(ind) = 0
             nst = nst + nmmax
            else
             idano(ind) = no
             idanv(ind) = nv
             idapo(ind) = nst + 1
             idalm(ind) = ndanum
             idall(ind) = 0
             nst = nst + ndanum
            endif
          endif
!
          if(nst.gt.lst) then
            x=-1.d0
+if cr
            write(lout,*)'ERROR IN DAALL, STACK EXHAUSTED '
+ei
+if .not.cr
            write(*,*)'ERROR IN DAALL, STACK EXHAUSTED '
+ei
+if cr
            write(lout,*) ' NST,LST '
+ei
+if .not.cr
            write(*,*) ' NST,LST '
+ei
+if cr
            write(lout,*)  nst,lst
+ei
+if .not.cr
            write(*,*)  nst,lst
+ei
+if cr
            write(lout,*) ' NDA,NDANUM,NDA*NDANUM '
+ei
+if .not.cr
            write(*,*) ' NDA,NDANUM,NDA*NDANUM '
+ei
+if cr
            write(lout,*)  nda,ndanum,nda*ndanum
+ei
+if .not.cr
            write(*,*)  nda,ndanum,nda*ndanum
+ei
!            X=DSQRT(X)
            call dadeb(31,'ERR DAALL ',1)
          endif
!
!          IF(NV.EQ.0) THEN
          if(nv.eq.0.or.nomax.eq.1) then
            call daclr(ic(i))
            idall(ic(i)) = idalm(ic(i))
          endif
        endif
  10  continue
!
      if(nda.gt.ndamaxi) ndamaxi=nda

      return
      end
!
+dk dadal
      subroutine dadal(idal,l)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,l,lda,lea,lia,lno,lnv,lst
!     ************************
!
!     THIS SUBROUTINE DEALLOCATES THE VECTORS IDAL
!
!-----------------------------------------------------------------------------1
+ca dabinc
+ca daname
!-----------------------------------------------------------------------------3
+ca alloc

      integer idal(*)

      do 10 i=l,1,-1
        if(idal(i).le.nomax+2.or.idal(i).gt.nda) then
+if cr
          write(lout,*)                                                 &
     &'ERROR IN ROUTINE DADAL, IDAL(I),NDA = ',idal(i),nda
+ei
+if .not.cr
          write(*,*)                                                    &
     &'ERROR IN ROUTINE DADAL, IDAL(I),NDA = ',idal(i),nda
+ei
          call dadeb(31,'ERR DADAL ',1)
        endif
        if(idal(i).eq.nda) then
!       deallocate
          nst = idapo(nda) - 1
          nda = nda - 1
!        else
!        write(6,'(a10)')daname(i)
!        write(6,*)' etienne',idal(i),nda
!        write(6,*) sqrt(-1.d0)
        endif

        allvec(idal(i)) = .false.

!        IDANO(IDAL(I)) = 0
!        IDANV(IDAL(I)) = 0
!        IDAPO(IDAL(I)) = 0
!        IDALM(IDAL(I)) = 0
         idall(idal(i)) = 0

        idal(i) = 0
  10  continue

      return
      end

+dk davar
      subroutine davar(ina,ckon,i)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ibase,ic1,ic2,illa,ilma,ina,inoa,inva,ipoa,lda,lea,lia, &
     &lno,lnv,lst
      double precision ckon
!     ****************************
!
!     THIS SUBROUTINE DECLARES THE DA VECTOR
!     AS THE INDEPENDENT VARIABLE NUMBER I.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!
!
      if(i.gt.inva) then
+if cr
         write(lout,*)'ERROR IN DAVAR, I = ',i,' EXCEEDS INVA = ',inva
+ei
+if .not.cr
         write(*,*)'ERROR IN DAVAR, I = ',i,' EXCEEDS INVA = ',inva
+ei
         call dadeb(31,'ERR DAVAR ',1)
      endif
!
      if(nomax.eq.1) then
         if(i.gt.inva) then
+if cr
            write(lout,*)                                               &
     &'ERROR IN DAVAR, I = ',i,' EXCEEDS INVA = ',inva
+ei
+if .not.cr
            write(*,*)                                                  &
     &'ERROR IN DAVAR, I = ',i,' EXCEEDS INVA = ',inva
+ei
!           CALL DADEB(31,'ERR DAVAR3',1)
         endif
         call daclr(ina)
         cc(ipoa) = ckon
         cc(ipoa+i) = 1d0
         return
      endif

      ibase = nomax+1
!
      if(i.gt.(nvmax+1)/2) then
        ic1 = 0
!hr10   ic2 = ibase**(i-(nvmax+1)/2-1)
        ic2 = ibase**((i-(nvmax+1)/2)-1)                                 !hr10
      else
        ic1 = ibase**(i-1)
        ic2 = 0
      endif
!
      if(abs(ckon).gt.eps) then
         idall(ina) = 2
         cc(ipoa) = ckon
         i1(ipoa) = 0
         i2(ipoa) = 0
!
         cc(ipoa+1) = 1.d0
         i1(ipoa+1) = ic1
         i2(ipoa+1) = ic2
      else
         idall(ina) = 1
         cc(ipoa) = 1.d0
         i1(ipoa) = ic1
         i2(ipoa) = ic2
      endif
!
      return
      end
!
+dk dacon
      subroutine dacon(ina,ckon)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illa,ilma,ina,inoa,inva,ipoa,lda,lea,lia,lno,lnv,lst
      double precision ckon
!     **************************
!
!     THIS SUBROUTINE SETS THE VECTOR C TO THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!
!
      if(nomax.eq.1) then
         call daclr(ina)
         cc(ipoa) = ckon
         return
      endif

      idall(ina) = 1
      cc(ipoa) = ckon
      i1(ipoa) = 0
      i2(ipoa) = 0
      if(abs(ckon).lt.eps) idall(ina) = 0
!
      return
      end
!
+dk danot
      subroutine danot(not)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer lda,lea,lia,lno,lnv,lst,not
!     *********************
!
!     THIS SUBROUTINE RESETS THE TRUNCATION ORDER NOCUT TO A NEW VALUE
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      if(not.gt.nomax) then
+if cr
         write(lout,*)'ERROR, NOCUT = ',nocut,' EXCEEDS NOMAX = ',nomax
+ei
+if .not.cr
         write(*,*)'ERROR, NOCUT = ',nocut,' EXCEEDS NOMAX = ',nomax
+ei
         call dadeb(31,'ERR DANOT ',1)
      endif
!
      nocut = not
+if debug
!     call warr('nocut',0d0,0,0,0,0)
+ei
!
      return
      end

+dk getdanot
      subroutine getdanot(not)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer lda,lea,lia,lno,lnv,lst,not
!     *********************
!
!     THIS SUBROUTINE RESETS THE TRUNCATION ORDER NOCUT TO A NEW VALUE
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      if(not.gt.nomax) then
+if cr
         write(lout,*)'ERROR, NOCUT = ',nocut,' EXCEEDS NOMAX = ',nomax
+ei
+if .not.cr
         write(*,*)'ERROR, NOCUT = ',nocut,' EXCEEDS NOMAX = ',nomax
+ei
         call dadeb(31,'ERR DANOT ',1)
      endif
!
      not=nocut
!
      return
      end

+dk daeps
      subroutine daeps(deps)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer lda,lea,lia,lno,lnv,lst
      double precision deps
!     **********************
!
!     THIS SUBROUTINE RESETS THE TRUNCATION ORDER NOCUT TO A NEW VALUE
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      if(deps.ge.0.d0) then
        eps = deps
      else
        deps=eps
      endif
!
      return
      end
!
+dk dapek
      subroutine dapek(ina,jj,cjj)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ibase,ic,ic1,ic2,icu,icz,ii1,ikk,illa,ilma,ina,         &
     &inoa,inva,ipek,ipoa,iu,iz,jj,jj1,lda,lea,lia,lno,lnv,lst,mchk
      double precision cjj
!     ****************************
!
!     THIS SUBROUTINE DETERMINES THE COEFFICIENT OF THE ARRAY
!     OF EXPONENTS JJ AND RETURNS IT IN CJJ
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      dimension jj(lnv)
+if debug
!Eric
!     integer umcalls,dapcalls,dokcalls,dumpl
!     common /mycalls/ umcalls,dapcalls,dokcalls,dumpl
+ei
+if debug
!     dapcalls=dapcalls+1
!     if (dapcalls.ge.606380)                                           &
!    &call warr('dapek1',0d0,dapcalls,ina,jj(1),1)
+ei
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!
+if debug
!     if (dapcalls.ge.606380) then
!     call warr('ina',0d0,ina,0,0,0)
!     call warr('inoa',0d0,inoa,0,0,0)
!     call warr('inva',0d0,inva,0,0,0)
!     call warr('ipoa',0d0,ipoa,0,0,0)
!     call warr('ilma',0d0,ilma,0,0,0)
!     call warr('illa',0d0,illa,0,0,0)
!     endif
+ei
!
      if(illa.eq.0) then   ! etienne shit
!hr10   cjj = 0
        cjj = 0d0                                                        !hr10
       return
      endif
      jj1 = 1
      if(inva.eq.0.or.nomax.eq.1) then
         if(inva.ne.0.and.nomax.eq.1) then
            if(illa.ge.2) then
               do 115 i=1,illa - 1
                  if(jj(i).eq.1) jj1 = i + 1
 115           continue
            else
               jj1 = jj(1) + 1
            endif
         else
            jj1 = jj(1)
         endif
         if(jj1.lt.1.or.jj1.gt.illa) then
+if cr
            write(lout,*)                                               &
     &'ERROR IN DAPEK, INDEX OUTSIDE RANGE, JJ(1) = ',jj1
+ei
+if .not.cr
            write(*,*)                                                  &
     &'ERROR IN DAPEK, INDEX OUTSIDE RANGE, JJ(1) = ',jj1
+ei
!           CALL DADEB(31,'ERR DAPEK1',1)
         endif
         ipek = ipoa + jj1 - 1
         cjj = cc(ipek)
+if debug
!     if (dapcalls.ge.606380)                                           &
!    &call warr('dapek2',cjj,2,ipek,0,0)
!     if (dapcalls.ge.606381) then
!       call dumpda('in dapek',606381,2)
!     call abend('in dapek 606381                                   ')
!     endif
+ei
         return
      endif

      ii1 = (nvmax+1)/2
      ibase = nomax+1
!
!     DETERMINE INDEX TO BE SEARCHED FOR
!     **********************************
!
      call dadcd(jj,ic1,ic2)
!
!ETIENNE
      if(ic1.gt.lia.or.ic2.gt.lia) then
+if cr
       write(lout,*) 'DISASTER IN DAPEK, INA= ',ina
       write(lout,*) ic1,ic2
       write(lout,*) (jj(ikk),ikk=1,lnv)
+ei
+if .not.cr
       write(*,*) 'DISASTER IN DAPEK, INA= ',ina
       write(*,*) ic1,ic2
       write(*,*) (jj(ikk),ikk=1,lnv)
+ei
      endif
!ETIENNE
      ic = ia1(ic1) + ia2(ic2)
!
!     DETERMINE IF MONOMIAL TO BE POKED CONFORMS WITH INOA, INVA,NOCUT
!     ****************************************************************
!
!      IF(ICO.GT.INOA.OR.ICV.GT.INVA.OR.ICO.GT.NOCUT) THEN
!         CJJ = 0
!         RETURN
!      ENDIF
!
!     DETERMINE IF MONOMIAL IS INSIDE FIRST AND LAST MONOMIALS OF A
!     *************************************************************
!
      iu = ipoa
      iz = ipoa + illa - 1
      icu = ia1(i1(iu))+ia2(i2(iu))
      icz = ia1(i1(iz))+ia2(i2(iz))
!
      if(illa.eq.0) then
!hr10    cjj = 0
         cjj = 0d0                                                       !hr10
         return
      elseif(ic.eq.icu) then
         cjj = cc(iu)
+if debug
!     if (dapcalls.ge.606380)                                           &
!    &call warr('dapek3',cjj,3,iu,0,0)
!     if (dapcalls.ge.606400) then
!       call dumpda('in dapek',606400,3)
!     call abend('                                                  ')
!     endif
+ei
         return
      elseif(ic.eq.icz) then
         cjj = cc(iz)
+if debug
!     if (dapcalls.ge.606380)                                           &
!    &call warr('dapek4',cjj,4,iz,0,0)
!     if (dapcalls.ge.606400) then
!       call dumpda('in dapek',606400,4)
!     call abend('                                                  ')
!     endif
+ei
         return
      elseif(ic.lt.icu.or.ic.gt.icz) then
!hr10    cjj = 0
         cjj = 0d0                                                       !hr10
         return
      endif
!
!     SEARCHING PROPER MONOMIAL
!     *************************
!
 10   continue
      if(iz-iu.le.1) then
         cjj = 0
         return
      endif
      i = (iu+iz)/2
!
!     if(ia1(i1(i))+ia2(i2(i)) - ic) 20,30,40
!hr10 mchk=ia1(i1(i))+ia2(i2(i)) - ic
      mchk=(ia1(i1(i))+ia2(i2(i))) - ic                                  !hr10
      if(mchk.lt.0) goto 20
      if(mchk.eq.0) goto 30
      if(mchk.gt.0) goto 40
 20   iu = i
      goto 10
 30   cjj = cc(i)
+if debug
!     if (dapcalls.ge.606380)                                           &
!    &call warr('dapek5',cjj,5,i,0,0)
!     if (dapcalls.ge.606400) then
!       call dumpda('in dapek',606400,5)
!     call abend('                                                  ')
!     endif
+ei
      return
 40   iz = i
      goto 10
!
      end
!
+dk dapok
      subroutine dapok(ina,jj,cjj)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,ic1,ic2,icu,icz,ii,illa,ilma,ina,inoa,inva,          &
     &ipoa,ipok,iu,iz,jj,jj1,lda,lea,lia,lno,lnv,lst,mchk
      double precision cjj
!     ****************************
!
!     THIS SUBROUTINE SETS THE COEFFICIENT OF THE ARRAY
!     OF EXPONENTS JJ TO THE VALUE CJJ
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      dimension jj(lnv)
!
+if debug
!Eric
!     integer umcalls,dapcalls,dokcalls,dumpl
!     common /mycalls/ umcalls,dapcalls,dokcalls,dumpl
+ei
+if debug
!     dokcalls=dokcalls+1
+ei
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!
!
+if debug
!      if (dokcalls.ge.445959) then
!      call wda('dapokcalls',cjj,dokcalls,0,0,0)
!      endif
+ei
+if debug
!      if (dokcalls.eq.445999) then
!      call dumpda('dapok666',999,8)
!      read (666) 
!      endif
+ei
      jj1 = 1
      if(inva.eq.0.or.nomax.eq.1) then
         if(inva.ne.0.and.nomax.eq.1) then
            if(illa.ge.2) then
               do 115 i=1,illa - 1
                  if(jj(i).eq.1) jj1 = i + 1
 115           continue
            else
               jj1 = jj(1) + 1
            endif
         else
            jj1 = jj(1)
         endif
         if(jj1.lt.1.or.jj1.gt.illa) then
+if cr
            write(lout,*)                                               &
     &'ERROR IN DAPOK, INDEX OUTSIDE RANGE, JJ(1) = ',jj1
+ei
+if .not.cr
            write(*,*)                                                  &
     &'ERROR IN DAPOK, INDEX OUTSIDE RANGE, JJ(1) = ',jj1
+ei
!           CALL DADEB(31,'ERR DAPOK1',1)
         endif
         ipok = ipoa + jj1 - 1
         cc(ipok) = cjj
+if debug
!      if (dokcalls.ge.445959) then
!      call wda('dapok',cjj,ipok,ipoa,jj1,0)
!      endif
+ei
+if debug
!      if (dokcalls.eq.445999) then
!      call dumpda('dapok666',999,9)
!      read (666) 
!      endif
+ei
         return
      endif

!     DETERMINE INDEX TO BE SEARCHED FOR
!     **********************************
!
      call dadcd(jj,ic1,ic2)
!
      ic = ia1(ic1) + ia2(ic2)
!
!     DETERMINE IF MONOMIAL TO BE POKED CONFORMS WITH INOA, INVA,NOCUT
!     ****************************************************************
!
!      IF(ICO.GT.INOA.OR.ICV.GT.INVA) THEN
!         write(6,*)'ERROR IN DAPOK, MONOMIAL NOT ALLOWED FOR ',A
!         CALL DADEB(31,'ERR DAPOK ',1)
!      ENDIF
!      IF(ICO.GT.NOCUT) RETURN
!
      if(illa.ne.0) then ! etienne shit
      iu = ipoa
      iz = ipoa + illa - 1
!
!     DETERMINE IF MONOMIAL IS INSIDE FIRST AND LAST MONOMIALS OF A
!     *************************************************************
!
      icu = ia1(i1(iu))+ia2(i2(iu))
      icz = ia1(i1(iz))+ia2(i2(iz))
      endif

      if(illa.eq.0) then
         i = ipoa
         goto 100
      elseif(ic.eq.icu) then
         cc(iu) = cjj
         i = iu
         goto 200
      elseif(ic.eq.icz) then
         cc(iz) = cjj
         i = iz
         goto 200
      elseif(ic.lt.icu) then
         i = iu
         goto 100
      elseif(ic.gt.icz) then
         i = iz + 1
         goto 100
      endif
!
!
!     SEARCHING PLACE TO POKE INTO OR BEFORE WHICH TO POKE
!     ****************************************************
!
      iu = ipoa
      iz = ipoa + illa
!
 10   continue
      if(iz-iu.le.1) then
         i = iz
         goto 100
      endif
      i = (iu+iz)/2
!
!      if(ia1(i1(i))+ia2(i2(i)) - ic) 20,30,40
!hr10 mchk=ia1(i1(i))+ia2(i2(i)) - ic
      mchk=(ia1(i1(i))+ia2(i2(i))) - ic                                  !hr10
      if(mchk.lt.0) goto 20
      if(mchk.eq.0) goto 30
      if(mchk.gt.0) goto 40
 20   iu = i
      goto 10
 30   cc(i) = cjj
      goto 200
 40   iz = i
      goto 10
!
!     INSERTING THE MONOMIAL, MOVING THE REST
!     ***************************************
!
 100  continue
!
+if debug
!      if (dokcalls.ge.445959) then
!      call wda('eps',eps,0,0,0,0)
!      call wda('cjj',cjj,0,0,0,0)
!      endif
+ei
      if(abs(cjj).lt.eps) return
!
      do 110 ii=ipoa+illa,i+1,-1
      cc(ii) = cc(ii-1)
      i1(ii) = i1(ii-1)
      i2(ii) = i2(ii-1)
 110  continue
!
      cc(i) = cjj
      i1(i) = ic1
      i2(i) = ic2
!
      idall(ina) = illa + 1
      if(idall(ina).gt.idalm(ina)) then
+if cr
         write(lout,*)'ERROR IN DAPOK '
+ei
+if .not.cr
         write(*,*)'ERROR IN DAPOK '
+ei
         call dadeb(31,'ERR DAPOK ',1)
      endif
!
      return
!
!     CASE OF CJJ = 0 WHICH MEANS MOVING THE REST
!     *********************************************
!
 200  continue
+if debug
!      if (dokcalls.ge.445959) then
!      call wda('eps',eps,1,1,1,1)
!      call wda('cjj',cjj,1,1,1,1)
!      endif
+ei
      if(abs(cjj).lt.eps) then
         do 210 ii=i,ipoa+illa-2
         cc(ii) = cc(ii+1)
         i1(ii) = i1(ii+1)
         i2(ii) = i2(ii+1)
 210     continue
         idall(ina) = illa - 1
      endif
      return
!
      end
!
+dk daclr
      subroutine daclr(inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,illc,ilmc,inc,inoc,invc,ipoc,lda,lea,lia,lno,lnv,lst
!     *********************
!
!     THIS SUBROUTINE SETS ALL THE STACK SPACE RESERVED FOR VARIABLE
!     C TO ZERO
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
      do 100 i=ipoc,ipoc+ilmc-1
!
      cc(i) = 0.d0
!
 100  continue
!
      return
      end
!
+dk dacop
      subroutine dacop(ina,inb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ia,ib,iif,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,&
     &ipoa,ipob,lda,lea,lia,lno,lnv,lst
!     *************************
!
!     THIS SUBROUTINE COPIES THE DA VECTOR A TO THE DA VECTOR B
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
+if debug
!     integer umcalls,dapcalls,dokcalls,dumpl
!     common /mycalls/ umcalls,dapcalls,dokcalls,dumpl
+ei
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)
!
      ib = ipob - 1
!
      iif = 0
      if(nomax.eq.1.or.inva.eq.0) iif = 1
+if debug
!     if (dokcalls.ge.445959) then
!       call wda('dacopiif',0d0,iif,ina,inb,nocut)
!ERIC
!       if (ina.eq.105.and.inb.eq.11) then
!         call dumpda('dacopiif',1,0)
!         read (555)
!       endif
!     endif
+ei

      do 100 ia = ipoa,ipoa+illa-1
!
      if(iif.eq.0) then
        if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 100
      endif
      ib = ib + 1
+if debug
!     if (dokcalls.ge.445959) then
!       call wda('daibb',cc(ib),ib,0,0,0)
!       call wda('daiba',cc(ia),ia,0,0,0)
!     endif
+ei
      cc(ib) = cc(ia)
      i1(ib) = i1(ia)
      i2(ib) = i2(ia)
+if debug
!     if (dokcalls.ge.445959) then
!       call wda('daibb2',cc(ib),ib,0,0,0)
!       call wda('daiba2',cc(ia),ia,0,0,0)
!     endif
+ei
!
 100  continue
!
!hr10 idall(inb) = ib - ipob + 1
      idall(inb) = (ib - ipob) + 1                                       !hr10
      if(idall(inb).gt.idalm(inb)) then
+if cr
         write(lout,*)'ERROR IN DACOP'
+ei
+if .not.cr
         write(*,*)'ERROR IN DACOP'
+ei
         call dadeb(31,'ERR DACOP ',1)
      endif
!
      return
      end

!
+dk datrashn
      subroutine datrashn(idif,ina,inbb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,idif,illa,ilma,ina,inb,inbb,inoa,inva,ipoa,lda,lea,  &
     &lia,lno,lnv,lst
      double precision rr
!     *************************
!
!     THIS SUBROUTINE COPIES THE DA VECTOR A TO THE DA VECTOR B
!
!-----------------------------------------------------------------------------1
+ca dabinc
      integer jd(lnv)
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      inb=0

      if(inbb.eq.ina) then
         call daall(inb,1,'$$DAADD $$',inoa,inva)
        else
         inb=inbb
      endif

      call daclr(inb)
!
!
      do 100 ia = ipoa,ipoa+illa-1
!
      if(nomax.ne.1) then
        call dancd(i1(ia),i2(ia),jd)
      else
        do i=1,lnv
          jd(i)=0
        enddo
        if(ia.ne.ipoa) then
          jd(ia-ipoa+1)=1
        endif
      endif

      call dapek(ina,jd,rr)
      jd(idif)=0
      if(abs(rr).gt.0.d0) call dapok(inb,jd,rr)
!
 100  continue
!
!
      if(inbb.eq.ina) then
         call dacop(inb,inbb)
         call dadal(inb,1)
      endif

      return
      end
!
+dk daadd
      subroutine daadd(ina,inb,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer idaadd,illc,ilmc,ina,inb,inc,inoc,invc,ipoc
!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA ADDITION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C.
!
      if(ina.ne.inc.and.inb.ne.inc) then
         call dalin(ina,+1.d0,inb,+1.d0,inc)
      else
         idaadd = 0
         call dainf(inc,inoc,invc,ipoc,ilmc,illc)
         call daall(idaadd,1,'$$DAADD $$',inoc,invc)
         call dalin(ina,+1.d0,inb,+1.d0,idaadd)
         call dacop(idaadd,inc)
         call dadal(idaadd,1)
      endif
!
      return
      end
!
+dk dasub
      subroutine dasub(ina,inb,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer idasub,illc,ilmc,ina,inb,inc,inoc,invc,ipoc
!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA SUBTRACTION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C.
!
      if(ina.ne.inc.and.inb.ne.inc) then
         call dalin(ina,+1.d0,inb,-1.d0,inc)
      else
         idasub = -1
         call dainf(inc,inoc,invc,ipoc,ilmc,illc)
         call daall(idasub,1,'$$DASUB $$',inoc,invc)
         call dalin(ina,+1.d0,inb,-1.d0,idasub)
         call dacop(idasub,inc)
         call dadal(idasub,1)
      endif
!
      return
      end
!
+dk damulin
      subroutine damulin(ina,inb,coe1,inc,ind,coe2,ine)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ina,inb,inc,incc,ind,ine,inoc,invc,lda,lea,lia,lno,lnv,   &
     &lst
      double precision coe1,coe2
!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA MULTIPLICATION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C. AS TEMPORARY STORAGE, THE STACK SPACE
!     OF THE (NOMAX+2) SCRATCH VARIABLES ALLOCATED BY DAINI IS USED.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
!

      call daall(incc,1,'$$DAJUNK$$',inoc,invc)
      call damul(ina,inb,incc)
      call damul(inc,ind,ine)
      call dalin(incc,coe1,ine,coe2,ine )
      call dadal(incc,1)

      return
      end

!
!
! ANFANG UNTERPROGRAMM
+dk daexx
      subroutine daexx(ina,inb,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illc,ilmc,ina,inaa,inb,inbb,inc,inoc,invc,ipoc,lda,lea,   &
     &lia,lno,lnv,lst
+ca dabinc
!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INE WITH THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
+if cr
       write(lout,*) "daexx"
+ei
+if .not.cr
       write(*,*) "daexx"
+ei
      if(ina.ne.inc.and.inb.ne.inc) then
         call daexxt(ina,inb,inc)
      else
         inaa = 0
         inbb = 0
         call dainf(inc,inoc,invc,ipoc,ilmc,illc)
         call daall(inaa,1,'$$DAADD $$',inoc,invc)
         call daall(inbb,1,'$$DAADD $$',inoc,invc)
         call dacop(ina,inaa)
         call dacop(inb,inbb)
         call daexxt(inaa,inbb,inc)
         call dadal(inaa,1)
         call dadal(inbb,1)
      endif

      return
      end

! ANFANG UNTERPROGRAMM
+dk daexxt
      subroutine daexxt(ina,inb,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer idaexx,illa,illb,illc,ilma,ilmb,ilmc,ina,inb,inc,inoa,    &
     &inob,inoc,inva,invb,invc,ipoa,ipob,ipoc,lda,lea,lia,lno,lnv,lst
!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INA WITH INB
!
!-----------------------------------------------------------------------------1
+ca dabinc
!-----------------------------------------------------------------------------9
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
      idaexx = 0
      call daall(idaexx,1,'$$DAEXX $$',inoc,invc)
      call dafun('LOG   ',ina,inc)
      call damul(inb,inc,idaexx)
      call dafun('EXP   ',idaexx,inc)
      call dadal(idaexx,1)
!
      return
      end

! ANFANG UNTERPROGRAMM
+dk dacex
      subroutine dacex(ina,ckon,inb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illc,ilmc,ina,inb,inc,incc,inoc,invc,ipoc,lda,lea,lia,lno,&
     &lnv,lst
      double precision ckon
+ca dabinc
!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INE WITH THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
+if cr
        write(lout,*) "dacex"
+ei
+if .not.cr
        write(*,*) "dacex"
+ei
      if(ina.eq.inb) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall(incc,1,'$$DAJUNK$$',inoc,invc)
        call dacext(ina,ckon,incc)
        call dacop(incc,inb)
        call dadal(incc,1)
      else
        call dacext(ina,ckon,inb)
      endif

      return
      end
+dk dacext
      subroutine dacext(ina,ckon,inb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer idacex,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,   &
     &ipoa,ipob,lda,lea,lia,lno,lnv,lst
      double precision ckon
!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES THE CONSTANT CKON WITH INA
!
!-----------------------------------------------------------------------------1
+ca dabinc
!-----------------------------------------------------------------------------9
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
      if(ckon.le.0) then
+if cr
         write(lout,*)'ERROR IN DACEX, CKON NOT POSITIVE'
+ei
+if .not.cr
         write(*,*)'ERROR IN DACEX, CKON NOT POSITIVE'
+ei
!        CALL DADEB(31,'ERR DACEX1',1)
      endif
!
      idacex = 0
      call daall(idacex,1,'$$DACEX $$',inob,invb)
+if crlibm
      ckon = log_rn(ckon)
+ei
+if .not.crlibm
      ckon = log(ckon)
+ei
      call dacmu(ina,ckon,idacex)
      call dafun('EXP   ',idacex,inb)
      call dadal(idacex,1)
!
      return
      end

+dk daexc
      subroutine daexc(ina,ckon,inb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illc,ilmc,ina,inb,inc,incc,inoc,invc,ipoc,lda,lea,lia,lno,&
     &lnv,lst
      double precision ckon
+ca dabinc
!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INE WITH THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
!        write(6,*) "daexc"

      if(ina.eq.inb) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall(incc,1,'$$DAJUNK$$',inoc,invc)
        call daexct(ina,ckon,incc)
        call dacop(incc,inb)
        call dadal(incc,1)
      else
        call daexct(ina,ckon,inb)
      endif

      return
      end

+dk daexct
      subroutine daexct(ina,ckon,inb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,idaexc,illa,illb,ilma,ilmb,ina,inb,inoa,inob,        &
     &inva,invb,ipoa,ipob,lda,lea,lia,lno,lnv,lst
      double precision ckon,xic
+ca dabinc
!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INE WITH THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
      idaexc = 0
      call daall(idaexc,1,'$$DAEXC $$',inob,invb)

      if(ckon.lt.0.d0) then
        call dafun('LOG   ',ina,inb)
        call dacmu(inb,ckon,idaexc)
        call dafun('EXP   ',idaexc,inb)
      else
!hr10   xic=abs(ckon-dble(idint(ckon)))
        xic=abs(ckon-dble(int(ckon)))                                    !hr10
        if(xic.gt.eps) then
          call dafun('LOG   ',ina,inb)
          call dacmu(inb,ckon,idaexc)
          call dafun('EXP   ',idaexc,inb)
        else
          ic=idint(ckon)
          call dacon(idaexc,1.d0)
          do i=1,ic
            call damul(idaexc,ina,idaexc)
          enddo
          call dacop(idaexc,inb)
        endif
      endif
      call dadal(idaexc,1)
!
      return
      end

+dk damul
      subroutine damul(ina,inb,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illc,ilmc,ina,inb,inc,incc,inoc,invc,ipoc,lda,lea,lia,lno,&
     &lnv,lst
!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA MULTIPLICATION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C. AS TEMPORARY STORAGE, THE STACK SPACE
!     OF THE (NOMAX+2) SCRATCH VARIABLES ALLOCATED BY DAINI IS USED.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
!

      if(ina.eq.inc.or.inb.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall(incc,1,'$$DAJUNK$$',inoc,invc)
        call damult(ina,inb,incc)
        call dacop(incc,inc)
        call dadal(incc,1)
      else
        call damult(ina,inb,inc)
      endif

      return
      end

+dk damult
      subroutine damult(ina,inb,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,i1ia,i2ia,ia,ib,ic,illa,illb,illc,ilma,ilmb,ilmc,ina,   &
     &inb,inc,inoa,inob,inoc,inva,invb,invc,ioffb,ipno,ipoa,ipob,ipoc,  &
     &ipos,lda,lea,lia,lno,lnv,lst,minv,noff,noib,nom
      double precision ccia,ccipoa,ccipob
!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA MULTIPLICATION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C. AS TEMPORARY STORAGE, THE STACK SPACE
!     OF THE (NOMAX+2) SCRATCH VARIABLES ALLOCATED BY DAINI IS USED.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      dimension ipno(0:lno),noff(0:lno)
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!      CALL DACHK(INA,INOA,INVA, INB,INOB,INVB, INC,INOC,INVC)
!
!     CASE OF FIRST ORDER ONLY
!     ************************

      if(nomax.eq.1) then
         minv = min(inva,invb,invc)
         ccipoa = cc(ipoa)
         ccipob = cc(ipob)
         cc(ipoc) = ccipoa*ccipob
         do 20 i=1,minv
  20     cc(ipoc+i) = ccipoa*cc(ipob+i) + ccipob*cc(ipoa+i)
         do 30 i=ipoc+minv+1,ipoc+invc
  30     cc(i) = 0d0
         return
      endif
!
!     GENERAL CASE
!     ************
!
      do 10 i=0,nomax
      noff(i) = idapo(i+2)
  10  ipno(i) = 0
!
      call daclr(1)
!
!     RE-SORTING THE VECTOR B INTO PIECES THAT ARE OF ONLY ONE ORDER
!     *************************************************************
!
      do 50 ib=ipob,ipob+illb-1
!
      noib = ieo(ia1(i1(ib))+ia2(i2(ib)))
      ipos = ipno(noib) + 1
      ipno(noib) = ipos
      inob = noff(noib) + ipos
!
      cc(inob) = cc(ib)
      i1(inob) = i1(ib)
      i2(inob) = i2(ib)
!
  50  continue
!
      do 60 i=0,nomax
  60  idall(i+2) = ipno(i)
!
!     PERFORMING ACTUAL MULTIPLICATION
!     ********************************
!
      nom = min(nocut,inoc)
!
      do 100 ia=ipoa,ipoa+illa-1
!
      i1ia = i1(ia)
      i2ia = i2(ia)
      ccia = cc(ia)
!
      do 100 noib = 0,nom-ieo(ia1(i1(ia))+ia2(i2(ia)))
!
      ioffb = noff(noib)
!
      do 100 ib = ioffb+1,ioffb+ipno(noib)
!
      ic = ia2(i2ia+i2(ib)) + ia1(i1ia + i1(ib))
      cc(ic) = cc(ic) + ccia*cc(ib)
!
 100  continue
!
      call dapac(inc)
!
      return
      end
!
+dk dadiv
      subroutine dadiv(ina,inb,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer idadiv,illc,ilmc,ina,inb,inc,inoc,invc,ipoc
!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA DIVISION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C.
!
      idadiv = 0
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
      call daall(idadiv,1,'$$DADIV $$',inoc,invc)
      call dafun('INV ',inb,idadiv)
      call damul(ina,idadiv,inc)
      call dadal(idadiv,1)
!
      return
      end
!
+dk dasqr
      subroutine dasqr(ina,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illc,ilmc,ina,inc,incc,inoc,invc,ipoc,lda,lea,lia,lno,lnv,&
     &lst
!     *************************
!
!     THIS SUBROUTINE SQUARES THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall(incc,1,'$$DAJUNK$$',inoc,invc)
        call dasqrt(ina,incc)
        call dacop(incc,inc)
        call dadal(incc,1)
      else
        call dasqrt(ina,inc)
      endif

      return
      end
+dk dasqrt
      subroutine dasqrt(ina,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,i1ia,i2ia,ia,ib,ib1,ic,illa,illc,ilma,ilmc,ina,inc,inoa,&
     &inoc,inva,invc,ioffa,ioffb,ipno,ipoa,ipoc,ipos,lda,lea,lia,lno,   &
     &lnv,lst,minv,noff,noia,noib,nom
      double precision ccia,ccipoa
!     *************************
!
!     THIS SUBROUTINE SQUARES THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      dimension ipno(0:lno),noff(0:lno)
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!      CALL DACHK(INA,INOA,INVA,'          ',-1,-1,INC,INOC,INVC)
!
      if(inva+invc.eq.0) then
         do 5 i=0,illa-1
!hr10  5      cc(ipoc+i) = cc(ipoa+i) * cc(ipoa+i)
  5      cc(ipoc+i) = cc(ipoa+i)**2                                      !hr10
         idall(inc) = idall(ina)
         if(idall(inc).gt.idalm(inc)) then
+if cr
            write(lout,*)'ERROR IN DASQR '
+ei
+if .not.cr
            write(*,*)'ERROR IN DASQR '
+ei
            call dadeb(31,'ERR DASQR ',1)
         endif
         return
      endif
!
!     CASE OF FIRST ORDER ONLY
!     ************************

      if(nomax.eq.1) then
         minv = min(inva,invc)
         ccipoa = cc(ipoa)
!hr10    cc(ipoc) = ccipoa*ccipoa
         cc(ipoc) = ccipoa**2                                            !hr10
         do 20 i=1,minv
  20     cc(ipoc+i) = 2d0*ccipoa*cc(ipoa+i)
         do 30 i=ipoc+minv+1,ipoc+invc
  30     cc(i) = 0d0
         return
      endif
!
!     GENERAL CASE
!     ************
!
      do 10 i=0,nomax
      noff(i) = idapo(i+2)
  10  ipno(i) = 0
!
      call daclr(1)
!
!     RESORTING THE VECTOR A INTO PIECES THAT ARE OF ONLY ONE ORDER
!     *************************************************************
!
      do 50 ia=ipoa,ipoa+illa-1
!
      noia = ieo(ia1(i1(ia))+ia2(i2(ia)))
      ipos = ipno(noia) + 1
      ipno(noia) = ipos
      inoa = noff(noia) + ipos
!
      cc(inoa) = cc(ia)
      i1(inoa) = i1(ia)
      i2(inoa) = i2(ia)
!
  50  continue
!
      do 60 i=0,nomax
  60  idall(i+2) = ipno(i)
!
!     PERFORMING ACTUAL MULTIPLICATION
!     ********************************
!
      nom = min(nocut,inoc)
!
      do 100 noia = 0,nom/2
!
      ioffa = noff(noia)
!
      do 100 ia=ioffa+1,ioffa+ipno(noia)
!
      i1ia = i1(ia)
      i2ia = i2(ia)
      ccia = cc(ia)
!
      ic = ia2(i2ia+i2ia) + ia1(i1ia+i1ia)
!hr10 cc(ic) = cc(ic) + ccia*ccia
      cc(ic) = cc(ic) + ccia**2                                          !hr10
      ccia = ccia + ccia
!
      do 100 noib = noia,nom-noia
!
      ioffb = noff(noib)
      if(noib.eq.noia) then
         ib1 = ia + 1
      else
         ib1 = ioffb + 1
      endif
!
      do 100 ib = ib1,ioffb+ipno(noib)
!
      ic = ia2(i2ia+i2(ib)) + ia1(i1ia + i1(ib))
      cc(ic) = cc(ic) + ccia*cc(ib)
!
 100  continue
!
      call dapac(inc)
!
      return
      end
!
+dk dacad
      subroutine dacad(ina,ckon,inb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,ipoa,ipob,&
     &lda,lea,lia,lno,lnv,lst
      double precision ckon,const
!     ******************************
!
!     THIS SUBROUTINE ADDS THE CONSTANT CKON TO THE VECTOR A
!
!-----------------------------------------------------------------------------1
+ca dabinc
      integer jj(lnv)
      data jj / lnv*0 /
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
!
      call dacop(ina,inb)
!
      if(nomax.eq.1) then
         cc(ipob) = cc(ipob) + ckon
         return
      endif
!
      call dapek(inb,jj,const)
      call dapok(inb,jj,const+ckon)
!
      return
      end
!
+dk dacsu
      subroutine dacsu(ina,ckon,inb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,ipoa,ipob,&
     &lda,lea,lia,lno,lnv,lst
      double precision ckon,const
!     ******************************
!
!     THIS SUBROUTINE SUBTRACTS THE CONSTANT CKON FROM THE VECTOR A
!
!-----------------------------------------------------------------------------1
+if debug
!     integer umcalls,dapcalls,dokcalls,dumpl
!     common /mycalls/ umcalls,dapcalls,dokcalls,dumpl
+ei
+ca dabinc
      integer jj(lnv)
      data jj / lnv*0 /
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
+if debug
!     if (umcalls.eq.8) then
!       if (dumpl.ne.0) then
! write the i's
!     call warr('ina',0d0,ina,0,0,0)
!     call warr('inoa',0d0,inoa,0,0,0)
!     call warr('inva',0d0,inva,0,0,0)
!     call warr('ipoa',0d0,ipoa,0,0,0)
!     call warr('ilma',0d0,ilma,0,0,0)
!     call warr('illa',0d0,illa,0,0,0)
!     call warr('inb',0d0,inb,0,0,0)
!     call warr('inob',0d0,inob,0,0,0)
!     call warr('invb',0d0,invb,0,0,0)
!     call warr('ipob',0d0,ipob,0,0,0)
!     call warr('ilmb',0d0,ilmb,0,0,0)
!     call warr('illb',0d0,illb,0,0,0)
!     call wda('bdacsu',0d0,0,0,0,0)
!       endif
!     endif
+ei
!
!
      call dacop(ina,inb)
!
+if debug
!     if (umcalls.eq.8) then
!       if (dumpl.ne.0) then
! write nomax
!     call wda('bnomax',0d0,nomax,0,0,0)
!       endif
!     endif
+ei

      if(nomax.eq.1) then
+if debug
!     if (umcalls.eq.8) then
!       call wda('bnomax1',cc(ipob),ipob,0,0,0)
!     endif
+ei
         cc(ipob) = cc(ipob) - ckon
+if debug
!ERIC THIS IS IT!
!     if (umcalls.eq.8) then
!       call wda('anomaxck',ckon,nomax,0,0,0)
!       call wda('anomax',cc(ipob),ipob,0,0,0)
!       if (dumpl.ne.0) then
!         call dumpda('adacsux',1,0)
!       read (444)
!       endif
!     endif
+ei
         return
+if debug
!     if (umcalls.eq.8) then
!       call wda('dacsu',cc(ipob),ipob,0,0,0)
!       if (dumpl.ne.0) then
!         call dumpda('adacsu',2,0)
!       read (444)
!       endif
!     endif
+ei
      endif
!
      call dapek(inb,jj,const)
+if debug
!     if (dokcalls.ge.445959) then
!       call wda('dacsu',const,inb,jj,0,0)
!     endif
+ei
      call dapok(inb,jj,const-ckon)
+if debug
!     if (dokcalls.ge.445959) then
!       call wda('dacsucc',const-ckon,inb,jj,0,0)
!       call wda('dacsuck',ckon,inb,jj,0,0)
!     endif
+ei
+if debug
!     if (umcalls.eq.8) then
!       if (dumpl.ne.0) then
!         call dumpda('adacsu',2,0)
!       read (444)
!       endif
!     endif
+ei
!
      return
      end
!
+dk dasuc
      subroutine dasuc(ina,ckon,inb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,ipoa,ipob,&
     &lda,lea,lia,lno,lnv,lst
      double precision ckon
!     ******************************
!
!     THIS SUBROUTINE SUBTRACTS THE VECTOR INA FROM THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
!
      call dacsu(ina,ckon,inb)
      call dacmu(inb,-1.d0,inb)
!
      return
      end
!
+dk dacmu
      subroutine dacmu(ina,ckon,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illc,ilmc,ina,inc,incc,inoc,invc,ipoc,lda,lea,lia,lno,lnv,&
     &lst
      double precision ckon
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1
+ca dabinc
+if debug
!     integer umcalls,dapcalls,dokcalls,dumpl
!     common /mycalls/ umcalls,dapcalls,dokcalls,dumpl
+ei
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall(incc,1,'$$DAJUNK$$',inoc,invc)
        call dacmut(ina,ckon,incc)
        call dacop(incc,inc)
+if debug
!     if (dokcalls.ge.445959) then
!       call wda('dacmuz',ckon,ina,inc,incc,0)
!     endif
+ei
        call dadal(incc,1)
      else
        call dacmut(ina,ckon,inc)
+if debug
!     if (dokcalls.ge.445959) then
!       call wda('dacmunz',ckon,ina,inc,0,0)
!     endif
+ei
      endif
      return
      end

+dk dacmut
      subroutine dacmut(ina,ckon,inb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,ib,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,  &
     &ipoa,ipob,lda,lea,lia,lno,lnv,lst,minv
      double precision ckon
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1
+ca dabinc
+if debug
!     integer umcalls,dapcalls,dokcalls,dumpl
!     common /mycalls/ umcalls,dapcalls,dokcalls,dumpl
+ei
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)
!
      if(nomax.eq.1) then
         minv = min(inva,invb)
         do 20 i=0,minv
  20     cc(ipob+i) = cc(ipoa+i) * ckon
         do 30 i=ipob+minv+1,ipob+invb
  30     cc(i) = 0d0
+if debug
!     if (dokcalls.ge.445959) then
!       call wda('dacmutz',ckon,ipoa,ipob,0,0)
!     endif
+ei
         return
      endif
!
      if(abs(ckon).lt.eps) then
         idall(inb) = 0
+if debug
!     if (dokcalls.ge.445959) then
!       call wda('dacmuteps',eps,inb,0,0,0)
!       call wda('dacmutck',ckon,inb,0,0,0)
!     endif
+ei
         return
      endif
!
      ib = ipob - 1
!
      do 100 ia=ipoa,ipoa+illa-1
!
      if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 100
      ib = ib + 1
      cc(ib) = cc(ia)*ckon
      i1(ib) = i1(ia)
      i2(ib) = i2(ia)
!
 100  continue
+if debug
!     if (dokcalls.ge.445959) then
!       call wda('dacmut100',ckon,ipoa,illa,0,0)
!     endif
+ei
!
!hr10 idall(inb) = ib-ipob+1
      idall(inb) = (ib-ipob)+1                                           !hr10
      if(idall(inb).gt.idalm(inb)) then
+if cr
         write(lout,*)'ERROR IN DACMU '
+ei
+if .not.cr
         write(*,*)'ERROR IN DACMU '
+ei
         call dadeb(31,'ERR DACMU ',1)
      endif
!
      return
      end
!
+dk dacdi
      subroutine dacdi(ina,ckon,inb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,ipoa,ipob,&
     &lda,lea,lia,lno,lnv,lst
      double precision ckon
!     ******************************
!
!     THIS SUBROUTINE DIVIDES THE VECTOR INA BY THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      if(ckon.eq.0.d0) then
+if cr
         write(lout,*)'ERROR IN DACDI, CKON IS ZERO'
+ei
+if .not.cr
         write(*,*)'ERROR IN DACDI, CKON IS ZERO'
+ei
         call dadeb(31,'ERR DACDI ',1)
      endif
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
!
      call dacmu(ina,1.d0/ckon,inb)
!
      return
      end
!
+dk dadic
      subroutine dadic(ina,ckon,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer idadic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,   &
     &ipoa,ipoc,lda,lea,lia,lno,lnv,lst
      double precision ckon,zero
      parameter(zero=0d0)
!     ******************************
!
!     THIS SUBROUTINE DIVIDES THE CONSTANT CKON BY THE VECTOR INA
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!
      if(abs(ckon).lt.eps) then
         call dacon(inc,zero)
         return
      endif

      idadic = 0
      call daall(idadic,1,'$$DADIC $$',inoc,invc)

      if(ckon.eq.0.d0) then
+if cr
         write(lout,*)'ERROR IN DACDI and DADIC, CKON IS ZERO'
+ei
+if .not.cr
         write(*,*)'ERROR IN DACDI and DADIC, CKON IS ZERO'
+ei
         call dadeb(31,'ERR DACDI ',1)
      endif
      call dacdi(ina,ckon,idadic)
      call dafun('INV ',idadic,inc)
      call dadal(idadic,1)
!
      return
      end
!
+dk dacma
      subroutine dacma(ina,inb,bfac,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer idacma,illc,ilmc,ina,inb,inc,inoc,invc,ipoc,lda,lea,lia,  &
     &lno,lnv,lst
      double precision bfac
!     **********************************
!
!     THIS SUBROUTINE PERFORMS THE OPERATIONS C = A + B*BFAC, WHERE A,B,C ARE
!     DA VECTORS AND BFAC IS A DOUBLE PRECISION. A AND C CAN BE IDENTICAL.
!     CAN LATER BE REPLACED BY SOMETHING LIKE DAADD WITH MINOR CHANGES.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
      idacma = 0
      call daall(idacma,1,'$$DACMA $$',inoc,invc)
      call dalin(ina,+1.d0,inb,bfac,idacma)
      call dacop(idacma,inc)
      call dadal(idacma,1)
!
      return
      end
!
+dk dalin
      subroutine dalin(ina,afac,inb,bfac,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illc,ilmc,ina,inb,inc,incc,inoc,invc,ipoc,lda,lea,lia,lno,&
     &lnv,lst
      double precision afac,bfac
!     ***************************************
!
!     THIS SUBROUTINE COMPUTES THE LINEAR COMBINATION
!     C = AFAC*A + BFAC*B. IT IS ALSO USED TO ADD AND SUBTRACT.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!

      if(ina.eq.inc.or.inb.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall(incc,1,'$$DAJUNK$$',inoc,invc)
        call dalint(ina,afac,inb,bfac,incc)
        call dacop(incc,inc)
        call dadal(incc,1)
      else
        call dalint(ina,afac,inb,bfac,inc)
      endif

      return
      end

+dk dalint
      subroutine dalint(ina,afac,inb,bfac,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,iamax,ib,ibmax,ic,icmax,illa,illb,illc,ilma,ilmb,    &
     &ilmc,ina,inb,inc,inoa,inob,inoc,inva,invb,invc,ipoa,ipob,ipoc,is, &
     &ismax,ismin,ja,jb,lda,lea,lia,lno,lnv,lst,minv,mchk
      double precision afac,bfac,ccc,copf
!     ***************************************
!
!     THIS SUBROUTINE COMPUTES THE LINEAR COMBINATION
!     C = AFAC*A + BFAC*B. IT IS ALSO USED TO ADD AND SUBTRACT.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!      CALL DACHK(INA,INOA,INVA, INB,INOB,INVB, INC,INOC,INVC)
!
!
      if(nomax.eq.1) then
         minv = min(inva,invb,invc)
         do 7 i=0,minv
 7       cc(ipoc+i) = cc(ipoa+i) * afac + cc(ipob+i) * bfac
         do 8 i=ipoc+minv+1,ipoc+invc
 8       cc(i) = 0d0
         return
      endif

      ia = ipoa
      ib = ipob
      ic = ipoc - 1
!hr10 iamax = ipoa+illa-1
      iamax = (ipoa+illa)-1                                              !hr10
!hr10 ibmax = ipob+illb-1
      ibmax = (ipob+illb)-1                                              !hr10
!hr10 icmax = ipoc+ilmc-1
      icmax = (ipoc+ilmc)-1                                              !hr10
      ja = ia1(i1(ia)) + ia2(i2(ia))
      jb = ia1(i1(ib)) + ia2(i2(ib))
!
      if(ia.gt.iamax) then
         ismin = ib
         ismax = ibmax
         copf  = bfac
         goto 50
      endif
      if(ib.gt.ibmax) then
         ismin = ia
         ismax = iamax
         copf  = afac
         goto 50
      endif
!
!     COMPARING
!     *********
!
  10  continue
!      if(ja-jb) 30,20,40
      mchk=ja-jb
      if(mchk.lt.0) goto 30
      if(mchk.eq.0) goto 20
      if(mchk.gt.0) goto 40
!
!     ADDING TWO TERMS
!     ****************
!
  20  continue
      ccc = cc(ia)*afac + cc(ib)*bfac
      if(abs(ccc).lt.eps) goto 25
      if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 25
      ic = ic + 1
      cc(ic) = ccc
      i1(ic) = i1(ia)
      i2(ic) = i2(ia)
  25  continue
      ia = ia + 1
      ib = ib + 1
      if(ia.gt.iamax) then
         ismin = ib
         ismax = ibmax
         copf  = bfac
         goto 50
      endif
      if(ib.gt.ibmax) then
         ismin = ia
         ismax = iamax
         copf  = afac
         goto 50
      endif
      ja = ia1(i1(ia)) + ia2(i2(ia))
      jb = ia1(i1(ib)) + ia2(i2(ib))
      goto 10
!
!     STORING TERM A
!     **************
!
  30  continue
      if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 35
      ccc = cc(ia)*afac
      if(abs(ccc).lt.eps) goto 35
      ic = ic + 1
      cc(ic) = ccc
      i1(ic) = i1(ia)
      i2(ic) = i2(ia)
  35  continue
      ia = ia + 1
      if(ia.gt.iamax) then
         ismin = ib
         ismax = ibmax
         copf  = bfac
         goto 50
      endif
      ja = ia1(i1(ia)) + ia2(i2(ia))
      goto 10
!
!     STORING TERM B
!     **************
!
  40  continue
      if(ieo(ia1(i1(ib))+ia2(i2(ib))).gt.nocut) goto 45
      ccc = cc(ib)*bfac
      if(abs(ccc).lt.eps) goto 45
      ic = ic + 1
      cc(ic) = ccc
      i1(ic) = i1(ib)
      i2(ic) = i2(ib)
  45  continue
      ib = ib + 1
      if(ib.gt.ibmax) then
         ismin = ia
         ismax = iamax
         copf  = afac
         goto 50
      endif
      jb = ia1(i1(ib)) + ia2(i2(ib))
      goto 10
!
!     COPYING THE REST
!     ****************
!
  50  continue
      do 60 is=ismin,ismax
      if(ieo(ia1(i1(is))+ia2(i2(is))).gt.nocut) goto 60
      ccc = cc(is)*copf
      if(abs(ccc).lt.eps) goto 60
      ic = ic + 1
      cc(ic) = ccc
      i1(ic) = i1(is)
      i2(ic) = i2(is)
  60  continue
!
!hr10 idall(inc) = ic - ipoc + 1
      idall(inc) = (ic - ipoc) + 1                                       !hr10
!
      if(idall(inc).gt.idalm(inc)) then
+if cr
         write(lout,*)'ERROR IN DALIN, RESULT HAS TOO MANY TERMS '
+ei
+if .not.cr
         write(*,*)'ERROR IN DALIN, RESULT HAS TOO MANY TERMS '
+ei
         call dadeb(31,'ERR DALIN ',1)
      endif
!
      return
      end
!
+dk dafun
      subroutine dafun(cf,ina,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illc,ilmc,ina,inc,incc,inoc,invc,ipoc,lda,lea,lia,lno,lnv,&
     &lst
!     ****************************
!
!     THIS SUBROUTINE COMPUTES THE FUNCTION CF OF THE DA VECTOR A
!     AND STORES THE RESULT IN C.
!     AT PRESENT, SOME FUNCTIONS CAN BE COMPUTED ONLY TO FIFTH ORDER.
!     THIS HAS TO BE FIXED IN THE FUTURE.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      character cf*4
+if debug
!     integer umcalls,dapcalls,dokcalls,dumpl
!     common /mycalls/ umcalls,dapcalls,dokcalls,dumpl
+ei

      if(ina.eq.inc) then
       call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
       call daall(incc,1,'$$DAJUNK$$',inoc,invc)
       call dafunt(cf,ina,incc)
       call dacop(incc,inc)
       call dadal(incc,1)
      else
       call dafunt(cf,ina,inc)
      endif
+if debug
!     if (umcalls.eq.8) then
!       call wda('dafun',0d0,0,0,0,0)
!       if (dumpl.ne.0) then
!         if (dumpl.eq.1) then
!           dumpl=dumpl+1
!         else
!           call dumpda('adafun',dumpl,0)
!           read (444)
!         endif
!       endif
!     endif
+ei

      return
      end

+dk dafunt
      subroutine dafunt(cf,ina,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,illa,illc,ilma,ilmc,ina,inc,ind,inoa,inoc,inon,         &
     &inva,invc,ipoa,ipoc,ipow,iscr,jj,lda,lea,lfun,lia,lno,lnv,lst,    &
     &no
      double precision a0,a1,a2,a3,a4,a5,ca,e1,                         &
     &e2,ea,era,p,ra,rpi4,sa,scr,                                       &
     &t,xf
!     ****************************
!
!     THIS SUBROUTINE COMPUTES THE FUNCTION CF OF THE DA VECTOR A
!     AND STORES THE RESULT IN C.
!     AT PRESENT, SOME FUNCTIONS CAN BE COMPUTED ONLY TO FIFTH ORDER.
!     THIS HAS TO BE FIXED IN THE FUTURE.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      character cf*4,cfh*4,abcs*26,abcc*26
      dimension xf(0:lno),jj(lnv)
!
      data jj /lnv*0/
      data abcs /'abcdefghijklmnopqrstuvwxyz'/
      data abcc /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
!
      if(cf(1:1).eq.' ') then
         cfh(1:3) = cf(2:4)
         cfh(1:4) = ' '
         cf = cfh
      endif
!
      do 5 i=1,4
      ind = index(abcs,cf(i:i))
      if(ind.ne.0) cf(i:i) = abcc(ind:ind)
  5   continue
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!     CASE OF NV = 0 WHICH MEANS COORDINATEWISE OPERATION
!     ***************************************************
!
!     CASE OF NV > 0 WHICH MEANS DIFFERENTIAL ALGEBRAIC OPERATION
!     ***********************************************************
!
      if(cf.eq.'SQR ') then
         call dasqr(ina,inc)
         return
      endif
!
!     ALLOCATE VARIABLES, PICK ZEROTH ORDER TERM
!     ******************************************
!
      ipow = 0
      inon = 0
      iscr = 0
!
      call daall(ipow,1,'$$DAFUN1$$',inoc,invc)
      call daall(inon,1,'$$DAFUN2$$',inoc,invc)
      call daall(iscr,1,'$$DAFUN3$$',inoc,invc)
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      call dapek(ina,jj,a0)
!
      no = min(nocut,inoa,inoc)
!
!     BRANCHING TO DIFFERENT FUNCTIONS
!     ********************************
!
      if(cf.eq.'INV ') then
!        1/(A0+P) = 1/A0*(1-(P/A0)+(P/A0)**2-...)
!hr10    if(a0.eq.0) then
         if(a0.eq.0d0) then                                              !hr10
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
         xf(0) = 1.d0/a0
         do 601 i=1,no
!hr10  601    xf(i) = -xf(i-1)/a0
  601    xf(i) = (-1d0*xf(i-1))/a0
!
      elseif(cf.eq.'SQRT') then
!        SQRT(A0+P) = SQRT(A0)*(1+1/2(P/A0)-1/8*(P/A0)**2+...)
!hr10    if(a0.le.0) then
         if(a0.le.0d0) then                                              !hr10
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
!hr10    ra = dsqrt(a0)
         ra = sqrt(a0)                                                   !hr10
         xf(0) = ra
         do 602 i=1,no
!hr10  602    xf(i) = -xf(i-1)/a0/dble(2*i)*dble(2*i-3)
  602    xf(i) = (((-1d0*xf(i-1))/a0)/dble(2*i))*dble(2*i-3)             !hr10
!
      elseif(cf.eq.'ISRT') then
!        1/SQRT(A0+P) = 1/SQRT(A0)*(1-1/2(P/A0)+3/8*(P/A0)**2-...)
!hr10    if(a0.le.0) then
         if(a0.le.0d0) then                                              !hr10
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
!hr10    era = 1.d0/dsqrt(a0)
         era = 1.d0/sqrt(a0)                                             !hr10
         xf(0) = era
         do 603 i=1,no
!hr10  603    xf(i) = -xf(i-1)/a0/dble(2*i)*dble(2*i-1)
  603    xf(i) = (((-1d0*xf(i-1))/a0)/dble(2*i))*dble(2*i-1)             !hr10
!
      elseif(cf.eq.'EXP ') then
!        EXP(A0+P) = EXP(A0)*(1+P+P**2/2!+...)
+if crlibm
         ea  = exp_rn(a0)
+ei
+if .not.crlibm
         ea  = exp(a0)
+ei
         xf(0) = ea
         do 604 i=1,no
  604    xf(i) = xf(i-1)/dble(i)
!
      elseif(cf.eq.'LOG ') then
!        LOG(A0+P) = LOG(A0) + (P/A0) - 1/2*(P/A0)**2 + 1/3*(P/A0)**3 - ...)
!hr10    if(a0.le.0) then
         if(a0.le.0d0) then                                              !hr10
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
+if crlibm
         ea  = log_rn(a0)
+ei
+if .not.crlibm
         ea  = log(a0)
+ei
         xf(0) = ea
         xf(1) = 1.d0/a0
         do 605 i=2,no
!hr10  605    xf(i) = -xf(i-1)/a0/dble(i)*dble(i-1)
  605    xf(i) = (((-1d0*xf(i-1))/a0)/dble(i))*dble(i-1)                 !hr10
!
      elseif(cf.eq.'SIN ') then
!        SIN(A0+P) = SIN(A0)*(1-P**2/2!+P**4/4!) + COS(A0)*(P-P**3/3!+P**5/5!)
+if crlibm
         sa  = sin_rn(a0)
+ei
+if .not.crlibm
         sa  = sin(a0)
+ei
+if crlibm
         ca  = cos_rn(a0)
+ei
+if .not.crlibm
         ca  = cos(a0)
+ei
         xf(0) = sa
         xf(1) = ca
         do 606 i=2,no
!hr10  606    xf(i) = -xf(i-2)/dble(i*(i-1))
  606    xf(i) = (-1d0*xf(i-2))/dble(i*(i-1))                            !hr10
!
      elseif(cf.eq.'COS ') then
!        COS(A0+P) = COS(A0)*(1-P**2/2!+P**4/4!) - SIN(A0)*(P-P**3/3!+P**5/5!)
+if crlibm
         sa  = sin_rn(a0)
+ei
+if .not.crlibm
         sa  = sin(a0)
+ei
+if crlibm
         ca  = cos_rn(a0)
+ei
+if .not.crlibm
         ca  = cos(a0)
+ei
         xf(0) = ca
!hr10    xf(1) = -sa
         xf(1) = -1d0*sa                                                 !hr10
         do 607 i=2,no
!hr10  607    xf(i) = -xf(i-2)/dble(i*(i-1))
  607    xf(i) = (-1d0*xf(i-2))/dble(i*(i-1))                            !hr10
!
      elseif(cf.eq.'SIRX') then
!        SIN(SQRT(P))/SQRT(P) = 1 - P/3! + P**2/5! - P**3/7! + ...
!hr10    if(a0.ne.0) then
         if(a0.ne.0d0) then                                              !hr10
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
         xf(0)=1.d0
         do 608 i=1,no
!hr10  608    xf(i) = -xf(i-1)/dble(2*i*(2*i+1))
  608    xf(i) = (-1d0*xf(i-1))/dble((2*i)*(2*i+1))                      !hr10
!
      elseif(cf.eq.'CORX') then
!        COS(SQRT(P)) = 1 - P/2! + P**2/4! - P**3/6! + ...
!hr10    if(a0.ne.0) then
         if(a0.ne.0d0) then                                              !hr10
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
         xf(0)=1.d0
         do 609 i=1,no
!hr10  609    xf(i) = -xf(i-1)/dble(2*i*(2*i-1))
  609    xf(i) = (-1d0*xf(i-1))/dble((2*i)*(2*i-1))                      !hr10
!
      elseif(cf.eq.'SIDX') then
!        SIN(P)/P = 1 - P**2/3! + P**4/5! - P**6/7! + ...
!hr10    if(a0.ne.0) then
         if(a0.ne.0d0) then                                              !hr10
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
         xf(0)=1.d0
         xf(1)=0.d0
         do 610 i=2,no
!hr10  610    xf(i) = -xf(i-2)/dble(i*(i+1))
  610    xf(i) = (-1d0*xf(i-2))/dble(i*(i+1))                            !hr10
!
      elseif(cf.eq.'TAN ') then
+if crlibm
         if(abs(cos_rn(a0)).lt.epsmac) then
+ei
+if .not.crlibm
         if(abs(cos(a0)).lt.epsmac) then
+ei
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
+if crlibm
         sa  = sin_rn(a0)
+ei
+if .not.crlibm
         sa  = sin(a0)
+ei
+if crlibm
         ca  = cos_rn(a0)
+ei
+if .not.crlibm
         ca  = cos(a0)
+ei
         xf(0) = sa/ca
!hr10    xf(1) = 1.d0/ca/ca
         xf(1) = (1.d0/ca)/ca                                            !hr10
!hr10    xf(2) = 2.d0*sa/ca/ca/ca/2.d0
         xf(2) = ((((2.d0*sa)/ca)/ca)/ca)/2.d0
!hr10    xf(3) = (2.d0*ca*ca+6.d0*sa*sa)/ca/ca/ca/ca/6.d0
         xf(3) = (((((2.d0*ca**2+6.d0*sa**2)/ca)/ca)/ca)/ca)/6.d0        !hr10
!hr10    xf(4) = (16*sa+8.d0*sa*sa*sa)/ca/ca/ca/ca/ca/24.d0
         xf(4) = ((((((16.d0*sa+8.d0*sa**3)/ca)/ca)/ca)/ca)/ca)/24.d0    !hr10
!hr10    xf(5) = (16.d0*ca*ca+24.d0*ca*ca*sa*sa+80.d0*sa*sa+            &
!hr10&40.d0*sa*sa*sa*sa)/ca/ca/ca/ca/ca/ca/120.d0
         xf(5) = (((((((((16.d0*ca**2+(24.d0*ca**2)*sa**2)+80.d0*sa**2)+&!hr10
     &40.d0*sa**4)/ca)/ca)/ca)/ca)/ca)/ca)/120.d0                        !hr10
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
            call dadeb(31,'ERR DAFUN ',1)
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
            stop 
+ei
         endif
      elseif(cf.eq.'COT ') then
+if crlibm
         if(abs(sin_rn(a0)).lt.epsmac) then
+ei
+if .not.crlibm
         if(abs(sin(a0)).lt.epsmac) then
+ei
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
+if crlibm
         sa  = sin_rn(a0)
+ei
+if .not.crlibm
         sa  = sin(a0)
+ei
+if crlibm
         ca  = cos_rn(a0)
+ei
+if .not.crlibm
         ca  = cos(a0)
+ei
         xf(0) = ca/sa
!hr10    xf(1) = -1.d0/sa/sa
         xf(1) = (-1.d0/sa)/sa
!hr10    xf(2) = 2.d0*ca/sa/sa/sa/2.d0
         xf(2) = ((((2.d0*ca)/sa)/sa)/sa)/2.d0                           !hr10
!hr10    xf(3) = -(2.d0*sa*sa+6.d0*ca*ca)/sa/sa/sa/sa/6.d0
         xf(3) = (((((-1d0*(2.d0*sa**2+6.d0*ca**2))/sa)/sa)/sa)/sa)/6.d0 !hr10
!hr10    xf(4) = (16*ca+8.d0*ca*ca*ca)/sa/sa/sa/sa/sa/24.d0
         xf(4) = ((((((16d0*ca+8.d0*ca**3)/sa)/sa)/sa)/sa)/sa)/24.d0     !hr10
!hr10    xf(5) = -(16.d0*sa*sa+24.d0*sa*sa*ca*ca+80.d0*ca*ca+           &
!hr10&40.d0*ca*ca*ca*ca)/sa/sa/sa/sa/sa/sa/120.d0
         xf(5) = (((((((-1d0*(((16.d0*sa**2+(24.d0*sa**2)*ca**2)+       &!hr10
     &80.d0*ca**2)+ 40.d0*ca**4))/sa)/sa)/sa)/sa)/sa)/sa)/120.d0         !hr10
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
            call dadeb(31,'ERR DAFUN ',1)
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
            stop 
+ei
         endif
      elseif(cf.eq.'ASIN') then
         if((1.d0-abs(a0)).lt.0.d0) then
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
+if crlibm
         xf(0) = asin_rn(a0)
+ei
+if .not.crlibm
         xf(0) = asin(a0)
+ei
!hr10 This code is not tested so leave **(-0.5d0) as it is.
!hr10 lf95 opt 1 gives a different result to opt 0 so should be changed to SQRT.
!hr10    xf(1) = (1.d0-a0*a0)**(-0.5d0)
!        xf(1) = (1.d0-a0**2)**(-0.5d0)                                  !hr10
         xf(1) = sqrt(1.d0-a0*a0)                                        !eric
!hr10    xf(2) = a0*xf(1)**3.d0/2.d0
!        xf(2) = (a0*xf(1)**3.d0)/2.d0                                   !hr10
         xf(2) = (a0*(xf(1)*xf(1)*xf(1)))/2.d0                            !eric
!hr10    xf(3) = (1+2.d0*a0*a0)*xf(1)**5.d0/6.d0
!        xf(3) = ((1.d0+2.d0*a0**2)*xf(1)**5.d0)/6.d0                    !hr10
         xf(3) = ((1.d0+2.d0*(a0*a0))*                                  &
     &           (xf(1)*xf(1)*xf(1)*xf(1)*xf(1)))/6.d0                   !eric
!hr10    xf(4) = (9.d0*a0+6.d0*a0*a0*a0)*xf(1)**7.d0/24.d0
!        xf(4) = ((9.d0*a0+6.d0*a0**3)*xf(1)**7.d0)/24.d0                !hr10
         xf(4) = ((9.d0*a0+6.d0*(a0*a0*a0))*                            &
     &           (xf(1)*xf(1)*xf(1)*xf(1)*xf(1)*xf(1)*xf(1)))/24.d0      !eric
!hr10    xf(5) = (9.d0+72.d0*a0*a0+24.d0*a0*a0*a0*a0)*xf(1)**9.d0/120.d0
!        xf(5) = ((9.d0+72.d0*a0**2+24.d0*a0**4)*xf(1)**9.d0)/120.d0     !hr10
         xf(5) = ((9.d0+72.d0*(a0*a0)+24.d0*(a0*a0*a0*a0))*             &
     &   (xf(1)*xf(1)*xf(1)*xf(1)*xf(1)*xf(1)*xf(1)*xf(1)*xf(1)))/120.d0 !eric
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
            stop 
+ei
         endif
      elseif(cf.eq.'ACOS')then
         if((1.d0-abs(a0)).lt.0.d0) then
            call dadeb(31,'ERR DAFUN ',1)
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            lfun = 0
            return
         endif
+if crlibm
         xf(0) =  acos_rn(a0)
+ei
+if .not.crlibm
         xf(0) =  acos(a0)
+ei
!hr10 This code is not tested so leave **(-0.5d0) as it is.
!hr10 lf95 opt 1 gives a different result to opt 0 so should be changed to SQRT.
!hr10    scr =  (1.d0-a0*a0)**(-0.5d0)
         scr =  (1.d0-a0**2)**(-0.5d0)                                   !hr10
         xf(1) =  -1d0*scr
!hr10    xf(2) = -a0*scr**3.d0/2.d0
         xf(2) = ((-1d0*a0)*scr**3.d0)/2.d0                              !hr10
!hr10    xf(3) = -(1+2.d0*a0*a0)*scr**5.d0/6.d0
         xf(3) = ((-1d0*(1.d0+2.d0*a0**2))*scr**5.d0)/6.d0               !hr10
!hr10    xf(4) = -(9.d0*a0+6.d0*a0*a0*a0)*scr**7.d0/24.d0
         xf(4) = ((-1d0*(9.d0*a0+6.d0*a0**3))*scr**7.d0)/24.d0           !hr10
!hr10    xf(5) = -(9.d0+72.d0*a0*a0+24.d0*a0*a0*a0*a0)*scr**9.d0/120.d0
         xf(5) =((-1d0*(9.d0+72.d0*a0**2+24.d0*a0**4))*scr**9.d0)/120.d0 !hr10
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'ATAN') then
!        ATAN(A0+P) = ATAN(A0)+1/(1+A0**2)*P-A0/(1+A0**2)**2*P**2+....)
+if crlibm
         xf(0) = atan_rn(a0)
+ei
+if .not.crlibm
         xf(0) = atan(a0)
+ei
         xf(1) = 1.d0/(1.d0+a0*a0)
         xf(2) = -a0*(xf(1)*xf(1))
         xf(3) = (a0*a0-1.d0/3.d0)*xf(1)**3
         xf(4) = (a0-a0*a0*a0)*xf(1)**4
         xf(5) = (1.d0/5.d0+a0**4-2.d0*a0*a0)*xf(1)**5
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'ACOT') then
+if crlibm
         xf(0) = 2.d0*atan_rn(1.d0)-atan_rn(a0)
+ei
+if .not.crlibm
         xf(0) = 2.d0*datan(1.d0)-atan(a0)
+ei
         scr = 1.d0/(1.d0+a0*a0)
         xf(1) = -scr
         xf(2) = a0*(scr*scr)
         xf(3) = -(a0*a0-1.d0/3.d0)*scr**3
         xf(4) = -(a0-a0*a0*a0)*scr**4
         xf(5) = -(1.d0/5.d0+a0**4-2.d0*a0*a0)*scr**5
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'SINH') then
+if crlibm
         sa  = sinh_rn(a0)
+ei
+if .not.crlibm
         sa  = sinh(a0)
+ei
+if crlibm
         ca  = cosh_rn(a0)
+ei
+if .not.crlibm
         ca  = cosh(a0)
+ei
         xf(0) = sa
         xf(1) = ca
         xf(2) = sa/2.d0
         xf(3) = ca/6.d0
         xf(4) = sa/24.d0
         xf(5) = ca/120.d0
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'COSH') then
+if crlibm
         sa  = sinh_rn(a0)
+ei
+if .not.crlibm
         sa  = sinh(a0)
+ei
+if crlibm
         ca  = cosh_rn(a0)
+ei
+if .not.crlibm
         ca  = cosh(a0)
+ei
         xf(0) = ca
         xf(1) = sa
         xf(2) = ca/2.d0
         xf(3) = sa/6.d0
         xf(4) = ca/24.d0
         xf(5) = sa/120.d0
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'TANH') then
+if crlibm
         sa  = sinh_rn(a0)
+ei
+if .not.crlibm
         sa  = sinh(a0)
+ei
+if crlibm
         ca  = cosh_rn(a0)
+ei
+if .not.crlibm
         ca  = cosh(a0)
+ei
         xf(0) = sa/ca
         xf(1) = 1.d0/ca/ca
         xf(2) = -2.d0*sa/ca/ca/ca/2.d0
         xf(3) = (-2.d0*ca*ca+6.d0*sa*sa)/ca/ca/ca/ca/6.d0
         xf(4) = (16*sa-8.d0*sa*sa*sa)/ca/ca/ca/ca/ca/24.d0
         xf(5) = (16.d0*ca*ca-24.d0*ca*ca*sa*sa-80.d0*sa*sa+            &
     &40.d0*sa*sa*sa*sa)/ca/ca/ca/ca/ca/ca/120.d0
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'COTH') then
+if crlibm
         if(abs(sinh_rn(a0)).lt.epsmac) then
+ei
+if .not.crlibm
         if(abs(sinh(a0)).lt.epsmac) then
+ei
            lfun = 0
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            call dadeb(31,'ERR DAFUN ',1)
            return
         endif
+if crlibm
         sa  = sinh_rn(a0)
+ei
+if .not.crlibm
         sa  = sinh(a0)
+ei
+if crlibm
         ca  = cosh_rn(a0)
+ei
+if .not.crlibm
         ca  = cosh(a0)
+ei
         xf(0) = ca/sa
         xf(1) = -1.d0/sa/sa
         xf(2) =  2.d0*ca/sa/sa/sa/2.d0
         xf(3) = (2.d0*sa*sa-6.d0*ca*ca)/sa/sa/sa/sa/6.d0
         xf(4) = (16*ca+8.d0*ca*ca*ca)/sa/sa/sa/sa/sa/24.d0
         xf(5) = (16.d0*sa*sa+24.d0*sa*sa*ca*ca-80.d0*ca*ca-            &
     &40.d0*ca*ca*ca*ca)/sa/sa/sa/sa/sa/sa/120.d0
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'ASNH') then
+if crlibm
         xf(0) = log_rn(a0+dsqrt(a0*a0+1.d0))
+ei
+if .not.crlibm
         xf(0) = log(a0+dsqrt(a0*a0+1.d0))
+ei
         xf(1) = (1.d0+a0*a0)**(-0.5d0)
         xf(2) = -a0*xf(1)**3.d0/2.d0
         xf(3) = (2.d0*a0*a0-1.d0)*xf(1)**5.d0/6.d0
         xf(4) = (9.d0*a0-6.d0*a0*a0*a0)*xf(1)**7.d0/24.d0
         xf(5) = (9.d0-72.d0*a0*a0+24.d0*a0*a0*a0*a0)*xf(1)**9.d0/120.d0
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'ACSH') then
         if((1.d0-a0).ge.0.d0) then
            lfun = 0
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            call dadeb(31,'ERR DAFUN ',1)
            return
         endif
+if crlibm
         xf(0) = log_rn(a0+dsqrt(a0*a0-1.d0))
+ei
+if .not.crlibm
         xf(0) = log(a0+dsqrt(a0*a0-1.d0))
+ei
         xf(1) = (a0*a0-1.d0)**(-0.5d0)
         xf(2) = -a0*xf(1)**3.d0/2.d0
         xf(3) = (2.d0*a0*a0+1.d0)*xf(1)**5.d0/6.d0
         xf(4) = (-9.d0*a0-6.d0*a0*a0*a0)*xf(1)**7.d0/24.d0
         xf(5) = (9.d0+72.d0*a0*a0+24.d0*a0*a0*a0*a0)*xf(1)**9.d0/120.d0
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'ATNH') then
         if((abs(a0)-1.d0).ge.0.d0) then
            lfun = 0
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            call dadeb(31,'ERR DAFUN ',1)
            return
         endif
+if crlibm
         xf(0) =  0.5d0*log_rn((1+a0)/(1-a0))
+ei
+if .not.crlibm
         xf(0) =  0.5d0*log((1+a0)/(1-a0))
+ei
         xf(1) =  1.d0/(1.d0-a0*a0)
         xf(2) =  a0*(xf(1)*xf(1))
         xf(3) = (a0*a0+1.d0/3.d0)*xf(1)**3
         xf(4) = (a0+a0*a0*a0)*xf(1)**4
         xf(5) = (1.d0/5.d0+a0**4+2.d0*a0*a0)*xf(1)**5
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'ACTH') then
         if(1.d0-abs(a0).ge.0.d0) then
            lfun = 0
+if cr
            write(lout,1000) cf,ina,a0
+ei
+if .not.cr
            write(*,1000) cf,ina,a0
+ei
            call dadeb(31,'ERR DAFUN ',1)
            return
         endif
+if crlibm
         xf(0) =  0.5d0*log_rn((a0+1)/(a0-1))
+ei
+if .not.crlibm
         xf(0) =  0.5d0*log((a0+1)/(a0-1))
+ei
         scr =  1.d0/(-1.d0+a0*a0)
         xf(1) = -scr
         xf(2) =  a0*(scr*scr)
         xf(3) = (-a0*a0-1.d0/3.d0)*scr**3.d0
         xf(4) = (a0+a0*a0*a0)*scr**4.d0
         xf(5) = (-1.d0/5.d0-a0**4-2.d0*a0*a0)*scr**5.d0
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
            call dadeb(31,'ERR DAFUN ',1)
         endif
!      ELSEIF(CF.EQ.'ABF ') THEN
!
!     DIESE FUNKTION BESCHREIBT DEN FELDABFALL BEI IONENOPTISCHEN ELEMENTEN
!     ABF=1/(1+EXP(A0+X))
!        =1/(1+EXP(A0)*(1-EXP(A0)/(1+EXP(A0))*X+....)
!         XF(0) = 1.D0/(1+DEXP(A0))
!         E1  = DEXP(A0)*X1
!         E1  = DEXP(A0)*Xf(0)
!         E2  = E1 * E1
!         E3  = E2 * E1
!         E4  = E3 * E1
!         E5  = E4 * E1
!         XF(1) = X1*(-E1)
!         XF(2) = X1*(-0.5D0* E1 + E2)
!         XF(3) = X1*(-E1/6.D0 + E2 - E3)
!         XF(4) = X1*(-E1/24.D0 + E2*7.D0/12.D0 - E3*3.D0/2.D0 + E4)
!         XF(5) = X1*(-E1/120.D0 + E2/4.D0 - E3*5.D0/4.D0 +
!     *         E4*2.D0 - E5)
!         IF(NO.GT.5) THEN
!            write(6,*)'ERROR IN DAFUN, ',CF, ' ONLY UP TO NO = 5'
!            CALL DADEB(31,'ERR DAFUN ',1)
!         ENDIF
!      ELSEIF(CF.EQ.'GAUS') THEN
!
!     DIESE FUNKTION BESCHREIBT DIE ENTWICKLUNG VON EXP(-X*X)
!
!         XF(0) = DEXP(-A0*A0)
!         XF(1) = -2.D0*A0*X1
!         XF(2) = (-1.D0+2.D0*A0*A0)*X1
!         XF(3) = (12.D0*A0-8.D0*A0*A0*A0)/6.D0*X1
!         XF(4) = (16.D0*A0*A0*A0*A0-48.D0*A0*A0+12.D0)/24.D0*X1
!         XF(5) = (-32.D0*A0*A0*A0*A0*A0+160.D0*A0*A0*A0-120.D0*A0)/
!     *           120.D0*X1
!         IF(NO.GT.5) THEN
!            write(6,*)'ERROR IN DAFUN, ',CF, ' ONLY UP TO NO = 5'
!            CALL DADEB(31,'ERR DAFUN ',1)
!         ENDIF
      elseif(cf.eq.'ERF ') then
!
!    ERF(X) STELLT DAS INTEGRAL VON 0 BIS X VON [ 2/SQRT(PI) * EXP(-X*X) ]
!    DAR
!
+if crlibm
         e1 = exp_rn(-a0*a0)
+ei
+if .not.crlibm
         e1 = exp(-a0*a0)
+ei
         a1 = .254829592d0
         a2 = -.284496736d0
         a3 = 1.421413741d0
         a4 = -1.453152027d0
         a5 = 1.061405429d0
         p  = .3275911d0
+if crlibm
         rpi4 = sqrt(atan_rn(1.d0))
+ei
+if .not.crlibm
         rpi4 = sqrt(atan(1.d0))
+ei
         t  = 1.d0/(1.d0+p*a0)
         e2 = 1.d0-t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))*e1
         xf(0)= e2
         xf(1) = e1/rpi4
         xf(2) = -a0*e1/rpi4
         xf(3) = (-2.d0+4.d0*a0*a0)/6.d0*e1/rpi4
         xf(4) = (12.d0*a0-8.d0*a0*a0*a0)/24.d0*e1/rpi4
         xf(5) = (16.d0*a0*a0*a0*a0-48.d0*a0*a0+12.d0)/120.d0*e1/rpi4
         if(no.gt.5) then
+if cr
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
+if .not.cr
            write(*,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
+ei
            call dadeb(31,'ERR DAFUN ',1)
         endif
      else
+if cr
         write(lout,*)'ERROR, UNSOPPORTED FUNCTION ',cf
+ei
+if .not.cr
         write(*,*)'ERROR, UNSOPPORTED FUNCTION ',cf
+ei
      endif
!
      call dacon(inc,xf(0))
      call dacop(ina,inon)
      call dapok(inon,jj,0.d0)
      call dacon(ipow,1.d0)
!
      do 800 i=1,min(no,nocut)
!
      call damul(inon,ipow,iscr)
      call dacop(iscr,ipow)
      call dacma(inc,ipow,xf(i),inc)
!
 800  continue
!
 1000 format('ERROR IN DAFUN, ',a4,' DOES NOT EXIST FOR VECTOR ',i10,   &
     &'CONST TERM  = ',e12.5)
!
      call dadal(iscr,1)
      call dadal(inon,1)
      call dadal(ipow,1)
!
      return
      end
!

+dk daabs
      subroutine daabs(ina,anorm)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,illa,ilma,ina,inoa,inva,ipoa,lda,lea,lia,lno,lnv,lst
      double precision anorm
!     ***************************
!
!     THIS SUBROUTINE COMPUTES THE NORM OF THE DA VECTOR A
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!
      anorm = 0.d0
      do 100 i=ipoa,ipoa+illa-1
      anorm = anorm + abs(cc(i))
 100  continue
!
      return
      end
!
+dk dacom
      subroutine dacom(ina,inb,dnorm)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer idacom,illc,ilmc,ina,inb,inc,inoc,invc,ipoc
      double precision dnorm
!     *******************************
!
!     THIS SUBROUTINE COMPARES TWO DA VECTORS BY RETURNING THE NORM
!     OF THE DIFFERENCE
!
      idacom = 0
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
      call daall(idacom,1,'$$DACOM $$',inoc,invc)
      call dasub(ina,inb,idacom)
      call daabs(idacom,dnorm)
      call dadal(idacom,1)
!
      return
      end
!

+dk dapos
      subroutine dapos(ina,inb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ia,ib,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,    &
     &ipoa,ipob,lda,lea,lia,lno,lnv,lst
!     *************************
!
!     THIS SUBROUTINE MAKES THE SIGNS OF ALL THE COEFFICIENTS OF A POSITIVE
!
!-----------------------------------------------------------------------------1
+ca dabinc

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)
!
      ib = ipob - 1
!
      do 100 ia = ipoa,ipoa+illa-1
!
      if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 100
      ib     = ib + 1
      cc(ib) = abs(cc(ia))
      i1(ib) = i1(ia)
      i2(ib) = i2(ia)
!
 100  continue
!
      idall(inb) = ib - ipob + 1
      if(idall(inb).gt.idalm(inb)) then
+if cr
         write(lout,*)'ERROR IN DAPOS '
+ei
+if .not.cr
         write(*,*)'ERROR IN DAPOS '
+ei
         call dadeb(31,'ERR DAPOS ',1)
      endif
!
      return
      end
!
+dk dacct
      subroutine dacct(ma,ia,mb,ib,mc,ic)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,ib,ic,ij,illc,ilmc,inoc,invc,ipoc,lda,lea,lia,lno,   &
     &lnv,lst
!     ***********************************
!
!     THIS SUBROUTINE PERFORMS A CONCATENATION MA = MB o MC
!     WHERE MA, MB AND MC ARE MATRICES CONSISTING OF IA, IB AND IC
!     DA VECTORS, RESPECTIVELY.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      integer mon(lnv),mb(*),mc(*),ma(*)

      if(ma(1).eq.mc(1).or.mb(1).eq.mc(1)) then
        call dainf(mc(1),inoc,invc,ipoc,ilmc,illc)
        do 101 ij=1,ic
 101    mon(ij)=0
        call daall(mon,ic,'$$DAJUNK$$',inoc,invc)
        call dacctt(ma,ia,mb,ib,mon,ic)
        do 9 i=1,ic
 9      call dacop(mon(i),mc(i))
        call dadal(mon,ic)
      else
        call dacctt(ma,ia,mb,ib,mc,ic)
      endif

      return
      end

+dk dacctt
      subroutine dacctt(mb,ib,mc,ic,ma,ia)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,ib,ic,iia,iib,iic,illa,illb,illc,ilma,ilmb,ilmc,inoa,&
     &inob,inoc,inva,invb,invc,ipoa,ipob,ipoc,iv,jl,jv,lda,lea,lia,lno, &
     &lnv,lst
      double precision ccf
!     ***********************************
!
!     THIS SUBROUTINE PERFORMS A CONCATENATION MA = MB o MC
!     WHERE MA, MB AND MC ARE MATRICES CONSISTING OF IA, IB AND IC
!     DA VECTORS, RESPECTIVELY.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
!      INTEGER MON(LNO+1),ICC(LNV),MB(*),MC(*),MA(*)
!ETIENNE
      integer mon(lno+1),icc(lno),mb(*),mc(*),ma(*)
!ETIENNE
!
!     CONSISTENCY CHECKS
!     ******************
!
      iia = ma(1)
      iib = mb(1)
      iic = mc(1)
      call dainf(iia,inoa,inva,ipoa,ilma,illa)
      call dainf(iib,inob,invb,ipob,ilmb,illb)
      call dainf(iic,inoc,invc,ipoc,ilmc,illc)
!
      call damch(ma,ia)
      call damch(mb,ib)
!
      if(ia.ne.ib) then
+if cr
         write(lout,*)'ERROR IN DACCT, IA .NE. IB'
+ei
+if .not.cr
         write(*,*)'ERROR IN DACCT, IA .NE. IB'
+ei
         call dadeb(31,'ERR DACCT1',1)
      elseif(ic.ne.invb) then
+if cr
         write(lout,*)'ERROR IN DACCT, IC.NE.INVB'
+ei
+if .not.cr
         write(*,*)'ERROR IN DACCT, IC.NE.INVB'
+ei
         call dadeb(31,'ERR DACCT2',1)
      endif
!
!     ALLOCATING LOCAL VECTORS AND CALLING MTREE
!     ******************************************
!
      do 50 i=1,ib
  50  icc(i) = 0
!
      do 60 i=1,nomax+1
  60  mon(i) = 0
!
      call daall(icc,ib,'$$DACCT $$',nomax,nvmax)
      call daall(mon,nomax+1,'$$DAMON $$',inoc,invc)
!
      call mtree(mb,ib,icc,ib)
!
!     PERFORMING CONCATENATION
!     ************************
!
      do 80 i=1,ia
      call dacon(ma(i),cc(idapo(icc(i))))
  80  continue
!
      call dacon(mon(1),1.d0)
!
      do 100 i=1,idall(icc(1))-1
!
      jl = i1(idapo(icc(1))+i)
      jv = i2(idapo(icc(1))+i)
!
      call damul(mon(jl),mc(jv),mon(jl+1))
!
      do 100 iv=1,ia
!
      ccf = cc(idapo(icc(iv))+i)
      if(abs(ccf).gt.eps) call dacma(ma(iv),mon(jl+1),ccf,ma(iv))
!
 100  continue
!
      call dadal(mon,nomax+1)
      call dadal(icc,ib)
!
      return
      end
!
+dk mtree
      subroutine mtree(mb,ib,mc,ic)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ib,ib1,ibi,ic,ic1,ic2,icc,ichk,ii,iib,iic,illb,illc,    &
     &ilmb,ilmc,inob,inoc,invb,invc,ipob,ipoc,j,jl,jnon,lda,lea,lia,lno,&
     &lnv,lst,nterm,ntermf
      double precision apek,bbijj,chkjj
!     *****************************
!
!     THIS SUBROUTINE IS USED FOR CONCATENATION AND TRACKING OF VECTORS
!     THROUGH A DA MAP. IT COMPUTES THE TREE THAT HAS TO BE TRANSVERSED
!     MB IS THE DA MATRIX WITH IA TERMS. THE OUTPUT MC IS A CA MATRIX WHICH
!     CONTAINS COEFFICIENTS AND CONTROL INTEGERS USED FOR THE TRAVERSAL.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      integer jj(lnv),jv(0:lno),mb(*),mc(*)
!
!     CONSISTENCY CHECKS
!     ******************
!
      iib = mb(1)
      iic = mc(1)
      call dainf(iib,inob,invb,ipob,ilmb,illb)
      call dainf(iic,inoc,invc,ipoc,ilmc,illc)
!
      call damch(mb,ib)
      call damch(mc,ic)
!
      if(ib.ne.ic) then
+if cr
         write(lout,*)'ERROR IN MTREE, IB .NE. IC'
+ei
+if .not.cr
         write(*,*)'ERROR IN MTREE, IB .NE. IC'
+ei
         call dadeb(31,'ERR MTREE1',1)
      endif
!
!     ALLOCATING LOCAL VECTORS
!     ************************
!
      ichk = 0
      call daall(ichk,1,'$$MTREE $$',nomax,nvmax)
!
!     FIND ALL THE ENTRIES TO BE LOOKED FOR
!     *************************************
!
      call daclr(1)
!
      cc(1) = 1.d0
!
      do 100 i=1,ib
      if(nomax.eq.1) then
      do 91 ib1 = 2,7
      cc(ib1) = 1d0
   91 continue
      else
      do  90 ibi = idapo(mb(i)),idapo(mb(i))+idall(mb(i))-1
      icc = ia1(i1(ibi)) + ia2(i2(ibi))
      if(ieo(icc).gt.inob) goto 90
      cc(icc) = 1.d0
   90 continue
      endif
  100 continue
!
      do 150 ii=1,inob
!
!     SEARCHING FOR FATHER FOR EACH TERM
!
      do 140 i=1,nmmax
      if(cc(i).lt.0.5d0) goto 140
!
      jnon = 0
      call dancd(i1(i),i2(i),jj)
      do 130 j=1,invb
      if(jj(j).eq.0) goto 130
      jnon = j
      jj(j) = jj(j) - 1
      call dadcd(jj,ic1,ic2)
      apek = cc(ia1(ic1)+ia2(ic2))
      jj(j) = jj(j) + 1
      if(apek.ge.0.5d0) goto 140
  130 continue
!
      if(jnon.eq.0) goto 140
!
!     TERM IS AN ORPHAN, SO CREATE FOSTER FATHER
!
      jj(jnon) = jj(jnon) - 1
      call dadcd(jj,ic1,ic2)
      cc(ia1(ic1)+ia2(ic2)) = 1.d0
!
  140 continue
  150 continue
!
      call dapac(ichk)
!ETIENNE      CALL DAPRI(ICHK,32)
!
!     SETTING UP TREE STRUCTURE
!     *************************
!
      ntermf = idall(ichk)
!
!     ZEROTH ORDER TERMS
!     ******************
!
      do 160 i=1,lnv
 160  jj(i) = 0
!
      do 170 i=1,ib
      call dapek(mb(i),jj,bbijj)
      i1(idapo(mc(i))) = 0
      i2(idapo(mc(i))) = 0
      cc(idapo(mc(i))) = bbijj
 170  continue
!
      call dapek(ichk,jj,chkjj)
      if(chkjj.gt.0.5d0) then
         call dapok(ichk,jj,-1.d0)
      else
+if cr
         write(lout,*)                                                  &
     &'ERROR IN MTREE, ZEROTH ORDER TERM OF ICHK IS ZERO'
+ei
+if .not.cr
         write(*,*)'ERROR IN MTREE, ZEROTH ORDER TERM OF ICHK IS ZERO'
+ei
         call dadeb(31,'ERR MTREE2',1)
      endif
!
      nterm = 1
!
!     HIGHER ORDER TERMS
!     ******************
!
      do 180 jl=1,inob
 180  jv(jl) = 0
!
      jl = 0
      chkjj = 1.d0
!
 200  continue
      if(jl.eq.0.and.chkjj.le.0.5d0) goto 250
      if(jl.lt.inob.and.chkjj.gt.0.5d0) then
         jl = jl + 1
         jj(1) = jj(1) + 1
         jv(jl) = 1
      elseif(jv(jl).eq.invb) then
         jj(jv(jl)) = jj(jv(jl)) - 1
         jv(jl) = 0
         jl = jl - 1
         chkjj = 0.d0
         goto 200
      else
         jj(jv(jl)) = jj(jv(jl)) - 1
         jv(jl) = jv(jl) + 1
         jj(jv(jl)) = jj(jv(jl)) + 1
      endif
!
      call dapek(ichk,jj,chkjj)
!
      if(chkjj.le.0.5d0) goto 200
!
      nterm = nterm + 1
      if(nterm.gt.idalm(mc(1))) then
+if cr
         write(lout,*)'ERROR IN MTREE, NTERM TOO LARGE'
+ei
+if .not.cr
         write(*,*)'ERROR IN MTREE, NTERM TOO LARGE'
+ei
         call dadeb(31,'ERR MTREE3',1)
      endif
!
      call dapok(ichk,jj,-1.d0)
!
!     write(6,*)'JL,JV = ',JL,JV(JL)
      do 210 i=1,ib
      call dapek(mb(i),jj,bbijj)
      i1(idapo(mc(i))+nterm-1) = jl
      i2(idapo(mc(i))+nterm-1) = jv(jl)
      cc(idapo(mc(i))+nterm-1) = bbijj
 210  continue
!
      goto 200
!
 250  continue
!
      do 260 i=1,ib
 260  idall(mc(i)) = nterm
!
!     PERFORMING CROSS CHECKS
!     ***********************
!
      if(nterm.ne.ntermf.or.nterm.ne.idall(ichk)) then
+if cr
         write(lout,*)'ERROR IN MTREE, NTERM, NTERMF, IDALL(ICHK) =  '
+ei
+if .not.cr
         write(*,*)'ERROR IN MTREE, NTERM, NTERMF, IDALL(ICHK) =  '     &
+ei
     &,nterm,ntermf,idall(ichk)
         call dadeb(31,'ERR MTREE4',1)
      endif
!
      do 270 i=idapo(ichk),idapo(ichk)+nterm-1
      if(abs(cc(i)+1.d0).gt.epsmac) then
+if cr
         write(lout,*)'ERROR IN MTREE, NOT ALL TERMS IN ICHK ARE -1'
+ei
+if .not.cr
         write(*,*)'ERROR IN MTREE, NOT ALL TERMS IN ICHK ARE -1'
+ei
         call dadeb(31,'ERR MTREE5',1)
      endif
 270  continue
!
      call dadal(ichk,1)
!
      return
      end
!
+dk ppushpr
      subroutine ppushpri(mc,ic,mf,jc,line)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,iv,jc,jl,jv,lda,lea,lia,lno,lnv,lst,mc,mf
+ca dabinc
      dimension mc(*)
      character*20 line
      if(mf.le.0) return
      write(mf,*) 0,0,jc+1,0,line
      do 10 i=1,ic
        jc=1+jc
        write(mf,*) jc,jl,jv,cc(idapo(mc(i)))
  10  continue
!     xf(i) = cc(idapo(mc(i)))
!      xm(1) = 1.d0
      do 100 i=1,idall(mc(1))-1
      jl = i1(idapo(mc(1))+i)
      jv = i2(idapo(mc(1))+i)
!      xx = xm(jl)*xi(jv)
!      xm(jl+1) = xx
      do 100 iv=1,ic
        jc=1+jc
        write(mf,*) jc,jl,jv,cc(idapo(mc(iv))+i)
!        xf(iv) = xf(iv) + cc(idapo(mc(iv))+i) * xx
 100  continue
      return
      end
!

+dk ppush
      subroutine ppush(mc,ic,xi,xf)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,iv,jl,jv,lda,lea,lia,lno,lnv,lst,mc
      double precision xf,xi,xm,xt,xx
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE MATRIX WHOSE TREE IS STORED IN CA VECTOR MC
!     TO THE COORDINATES IN XI AND STORES THE RESULT IN XF
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      dimension mc(*),xf(*),xi(*),xm(lno+1) ,xt(lno)
!
      do i=1,ic
        xt(i)=xi(i)
      enddo
      do 10 i=1,ic
  10  xf(i) = cc(idapo(mc(i)))
!
      xm(1) = 1.d0
!
      do 100 i=1,idall(mc(1))-1
!
      jl = i1(idapo(mc(1))+i)
      jv = i2(idapo(mc(1))+i)
      xx = xm(jl)*xt(jv)
      xm(jl+1) = xx
!
      do 100 iv=1,ic
      xf(iv) = xf(iv) + cc(idapo(mc(iv))+i) * xx
 100  continue
!
      return
      end

+dk ppush1
      subroutine ppush1(mc,xi,xf)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,jl,jv,lda,lea,lia,lno,lnv,lst,mc
      double precision xf,xi,xm,xt,xx
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE MATRIX WHOSE TREE IS STORED IN CA VECTOR MC
!     TO THE COORDINATES IN XI AND STORES THE RESULT IN XF
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      dimension xi(*),xm(lno+1) ,xt(lno)
!
      do i=1,nvmax
        xt(i)=xi(i)
      enddo

      xf = cc(idapo(mc))
!
      xm(1) = 1.d0
!
      do 100 i=1,idall(mc)-1
!
      jl = i1(idapo(mc)+i)
      jv = i2(idapo(mc)+i)
      xx = xm(jl)*xt(jv)
      xm(jl+1) = xx
!
      xf = xf + cc(idapo(mc)+i) * xx
 100  continue
!
      return
      end

+dk dainv
      subroutine dainv(ma,ia,mb,ib)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,ib,ij,illb,ilmb,inob,invb,ipob,lda,lea,lia,lno,lnv,  &
     &lst
      double precision x
!     *****************************
!
!     THIS SUBROUTINE INVERTS THE MATRIX MA WITH IA DA VECTORS AND
!     STORES THE RESULT IN MI
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      integer jj(lnv),ml(lnv),ma(*),mb(*)
!
      dimension x(lnv)
!

      do 1 i=1,lnv
 1    jj(i)=0

      if(ma(1).eq.mb(1)) then
        call dainf(mb(1),inob,invb,ipob,ilmb,illb)
        do 9 i=1,ia
 9      call dapok(ma(i),jj,0.d0)
        do 101 ij=1,ib
 101    ml(ij)=0
        call daall(ml,ib,'$$DAJUNK$$',inob,invb)
        call dainvt(ma,ia,ml,ib)
        do 90 i=1,ib
 90     call dacop(ml(i),mb(i))
        call dadal(ml,ib)
      else
        do 99 i=1,ia
        call dapek(ma(i),jj,x(i))
 99     call dapok(ma(i),jj,0.d0)
        call dainvt(ma,ia,mb,ib)
        do 999 i=1,ia
 999    call dapok(ma(i),jj,x(i))
      endif

      return
      end

+dk dainvt
      subroutine dainvt(ma,ia,mb,ib)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,ib,ie,ier,illa,illb,ilma,ilmb,inoa,inob,inva,invb,   &
     &ipoa,ipob,j,k,lda,lea,lia,lno,lnv,lst,nocut0
      double precision aa,ai,amjj,amsjj,prod
!     *****************************
!
!     THIS SUBROUTINE INVERTS THE MATRIX MA WITH IA DA VECTORS AND
!     STORES THE RESULT IN MI
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      integer jj(lnv),ms(lnv),ml(lnv),ma(*),mb(*)
!
      dimension aa(lnv,lnv),ai(lnv,lnv)
!
      call dainf(ma(1),inoa,inva,ipoa,ilma,illa)
      call dainf(mb(1),inob,invb,ipob,ilmb,illb)
!
!     CONSISTENCY CHECKS
!     ******************
!
      call damch(ma,ia)
      call damch(mb,ib)
!etienne
      do 959 ie=1,ib
 959  call dacon(mb(ie),0.d0)
!etienne
!
      if(ia.ne.ib) then
+if cr
         write(lout,*)'ERROR IN DAINV, IA .NE. IB'
+ei
+if .not.cr
         write(*,*)'ERROR IN DAINV, IA .NE. IB'
+ei
         call dadeb(31,'ERR DAINV1',1)
      elseif(ia.ne.inva.or.ib.ne.invb) then
+if cr
         write(lout,*)'ERROR IN DAINV, IA.NE.INVA.OR.IB.NE.INVB'
+ei
+if .not.cr
         write(*,*)'ERROR IN DAINV, IA.NE.INVA.OR.IB.NE.INVB'
+ei
         call dadeb(31,'ERR DAINV2',1)
      endif
!
!     ALLOCATING LOCAL VECTORS
!     ************************
!
      do 10 i=1,ia
      ms(i) = 0
      ml(i) = 0
  10  continue
!
      call daall(ms,ia,'$$INV   $$',inoa,inva)
      call daall(ml,ia,'$$INVL  $$',inoa,inva)
!
!     EXTRACTING LINEAR MATRIX, GENERATING NONLINEAR PART OF A
!     ********************************************************
!
      do 115 i=1,ib
      do 110 j=1,ib
      do 105 k=1,ib
 105  jj(k) = 0
      jj(j) = 1
      call dapek(ma(i),jj,amjj)
      if(abs(amjj).gt.eps) call dapok(ma(i),jj,0.d0)
 110  aa(i,j) = amjj
      call dacmu(ma(i),-1.d0,ma(i))
 115  continue
!
!     INVERTING LINEAR MATRIX, CHECKING RESULT AND STORING IN ML
!     **********************************************************
!
      call matinv(aa,ai,ia,lnv,ier)
!
      if(ier.eq.132) then
+if cr
         write(lout,*)'ERROR IN ROUTINE DAINV'
+ei
+if .not.cr
         write(*,*)'ERROR IN ROUTINE DAINV'
+ei
         call dadeb(31,'ERR DAINV3',1)
      endif
!
      ier = 0
      do 140 i=1,ib
      do 140 j=1,ib
      prod = 0.d0
      do 120 k=1,ib
      jj(k) = 0
 120  prod = prod + aa(i,k)*ai(k,j)
      if(i.eq.j) prod = prod - 1.d0
      if(abs(prod).gt.100*epsmac) then
+if cr
         write(lout,*)                                                  &
     &'ERROR IN DAINV, INVERSION DID NOT WORK,I,J,PROD = ',
+ei
+if .not.cr
         write(*,*)'ERROR IN DAINV, INVERSION DID NOT WORK,I,J,PROD = ',&
+ei
     &i,j,prod,epsmac,eps
         ier = 1
!ETIENNE
      return
!ETIENNE
      endif
      jj(j) = 1
      call dapok(mb(i),jj,ai(i,j))
      call dapok(ml(i),jj,ai(i,j))
 140  continue
!
      if(ier.eq.1) call dadeb(31,'ERR DAINV4',1)
!
!     ITERATIVELY COMPUTING DIFFERENT PARTS OF THE INVERSE
!     ****************************************************
!
!     MB (OF ORDER I) = A1^-1 o [ E - ANL (NONLINEAR) o MB (OF ORDER I) ]
!
      nocut0 = nocut
!
      do 250 i=2,nocut
!
      nocut = i
!
      call dacct(ma,ia,mb,ib,ms,ia)
      do 240 j=1,ib
      do 230 k=1,ib
 230  jj(k) = 0
      jj(j) = 1
      call dapek(ms(j),jj,amsjj)
      call dapok(ms(j),jj,amsjj+1.d0)
 240  continue
!
      call dacct(ml,ia,ms,ia,mb,ib)
!
 250  continue
!
      nocut = nocut0
!
!     FLIPPING BACK SIGN OF A, FILLING UP FIRST ORDER PART AGAIN
!     **********************************************************
!
      do 320 i=1,ib
      call dacmu(ma(i),-1.d0,ma(i))
      do 320 j=1,ib
      do 310 k=1,ib
 310  jj(k) = 0
      jj(j) = 1
      call dapok(ma(i),jj,aa(i,j))
 320  continue
!
      call dadal(ml,ia)
      call dadal(ms,ia)
!
      return
      end
!
+dk matinv
      subroutine matinv(a,ai,n,nmx,ier)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ier,indx,j,n,nmax,nmx
      double precision a,ai,aw,d
!     *********************************
!
!     THIS SUBROUTINE INVERTS THE MATRIX A AND STORES THE RESULT IN AI
!     INPUT  A   - SAVED
!            N   - ORDER OF MATRIX < 100
!     OUTPUT AI  - A INVERSE
!            IER - 0 NO ERROR
!                  132 ZERO DETERMINANT
!
      parameter (nmax=400)
      dimension a(nmx,nmx),ai(nmx,nmx),aw(nmax,nmax),indx(nmax)

      do 12 i=1,n
         do 11 j=1,n
            aw(i,j) = a(i,j)
11       ai(i,j) = 0.0
12    ai(i,i) = 1.d0

      call ludcmp(aw,n,nmax,indx,d,ier)
      if (ier .eq. 132) return
      do 13 j=1,n
13    call lubksb(aw,n,nmax,indx,ai(1,j),nmx)
!
      return
      end
!
+dk ludcmp
      subroutine ludcmp(a,n,np,indx,d,ier)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ier,imax,indx,j,k,n,nmax,np
      double precision a,aamax,d,dum,sum,tiny,vv
!     ************************************
!
!     THIS SUBROUTINE DECOMPOSES A MATRIX INTO LU FORMAT
!     INPUT A: NXN MATRIX - WILL BE OVERWRITTEN BY THE LU DECOMP.
!           NP: PHYSICAL DIMENSION OF A
!           INDX: ROW PERMUTATION VECTOR
!           D: EVEN OR ODD ROW INTERCHANGES
!
!     REFERENCE: NUMERICAL RECIPIES BY PRESS ET AL (CAMBRIDGE) PG. 35
!
      parameter (nmax = 400, tiny = 1.0e-20)
      dimension a(np,np), indx(np), vv(nmax)
      ier=0.
      d=1.d0
      do 12 i=1,n
         aamax=0.d0
         do 11 j=1,n
            if(abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11       continue
         if(aamax.eq.0.d0) then
            ier=132
            return
         endif
         vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
         if(j.gt.1) then
            do 14 i=1,j-1
               sum=a(i,j)
               if(i.gt.1) then
                  do 13 k=1,i-1
                     sum=sum-a(i,k)*a(k,j)
13                continue
                  a(i,j)=sum
               endif
14          continue
         endif
         aamax=0.d0
         do 16 i=j,n
            sum=a(i,j)
            if (j.gt.1) then
               do 15 k=1,j-1
                  sum=sum-a(i,k)*a(k,j)
15             continue
               a(i,j)=sum
            endif
            dum=vv(i)*abs(sum)
            if(dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
16       continue
         if (j.ne.imax) then
            do 17 k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
17          continue
            d=-d
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(j.ne.n) then
            if(a(j,j).eq.0.d0) a(j,j)=tiny
            dum=1./a(j,j)
            do 18 i=j+1,n
               a(i,j)=a(i,j)*dum
18          continue
         endif
19    continue
      if(a(n,n).eq.0.d0) a(n,n)=tiny
      return
      end
!
+dk lubksb
      subroutine lubksb(a,n,np,indx,b,nmx)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ii,indx,j,ll,n,nmx,np
      double precision a,b,sum
!     ************************************
!
!     THIS SUBROUTINE SOLVES SET OF LINEAR EQUATIONS AX=B,
!     INPUT A: NXN MATRIX IN lu FORM GIVEN BY LUDCMP
!           NP: PHYSICAL DIMENSION OF A
!           INDX: ROW PERMUTATION VECTOR
!           D: EVEN OR ODD ROW INTERCHANGES
!           B: RHS OF LINEAR EQUATION - WILL BE OVERWRITTEN BY X
!
!     REFERENCE: NUMERICAL RECIPIES BY PRESS ET AL (CAMBRIDGE) PG. 36
!
      dimension a(np,np), indx(np), b(nmx)
      ii = 0
      do 12 i=1,n
         ll = indx(i)
         sum = b(ll)
         b(ll) = b(i)
         if(ii.ne.0) then
            do 11 j=ii,i-1
               sum = sum-a(i,j)*b(j)
11          continue
         else if (sum.ne.0.d0) then
            ii = i
         endif
         b(i)=sum
12    continue
      do 14 i=n,1,-1
         sum=b(i)
         if(i.lt.n) then
            do 13 j=i+1,n
               sum = sum-a(i,j)*b(j)
13          continue
         endif

         b(i)=sum/a(i,i)

14    continue
      return
      end
!

+dk dapin
      subroutine dapin(ma,ia,mb,ib,jx)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,ib,ij,illb,ilmb,inob,invb,ipob,lda,lea,lia,lno,lnv,  &
     &lst
      double precision x
!     *****************************
!
!     THIS SUBROUTINE INVERTS THE MATRIX MA WITH IA DA VECTORS AND
!     STORES THE RESULT IN MI
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      integer jj(lnv),ml(lnv),ma(*),mb(*),jx(*)
!
      dimension x(lnv)
!

      do 1 i=1,lnv
 1    jj(i)=0

      if(ma(1).eq.mb(1)) then
        call dainf(mb(1),inob,invb,ipob,ilmb,illb)
        do 9 i=1,ia
 9      call dapok(ma(i),jj,0.d0)
        do 101 ij=1,ib
 101    ml(ij)=0
        call daall(ml,ib,'$$DAJUNK$$',inob,invb)
        call dapint(ma,ia,ml,ib,jx)
        do 90 i=1,ib
 90     call dacop(ml(i),mb(i))
        call dadal(ml,ib)
      else
        do 99 i=1,ia
        call dapek(ma(i),jj,x(i))
 99     call dapok(ma(i),jj,0.d0)
        call dapint(ma,ia,mb,ib,jx)
        do 999 i=1,ia
 999    call dapok(ma(i),jj,x(i))
      endif

      return
      end

+dk dapint
      subroutine dapint(ma,ia,mb,ib,jind)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,ib,illa,ilma,inoa,inva,ipoa,k,lda,lea,lia,lno,lnv,   &
     &lst
!     **********************************
!
!     THIS SUBROUTINE PERFORMS A PARTIAL INVERSION OF THE ROWS MARKED WITH
!     NONZERO ENTRIES IN JJ OF THE MATRIX A. THE RESULT IS STORED IN B.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      integer jj(lnv),jind(*),ma(*),mb(*),mn(lnv),mi(lnv),me(lnv)
!
      call dainf(ma(1),inoa,inva,ipoa,ilma,illa)
!

      do 5 i=1,ia
      mn(i) = 0
      mi(i) = 0
      me(i) = 0
  5   continue
!
      call daall(mn,ia,'$$PIN1  $$',inoa,inva)
      call daall(mi,ia,'$$PIN2  $$',inoa,inva)
      call daall(me,ia,'$$PIN3  $$',inoa,inva)
!
      do 10 i=1,ia
      do  8 k=1,nvmax
  8   jj(k) = 0
      jj(i) = 1
  10  call dapok(me(i),jj,1.d0)
!
      do 20 i=1,ia
      call dacop(ma(i),mn(i))
      if(jind(i).eq.0) call dacop(me(i),mn(i))
  20  continue
!
      call dainv(mn,ia,mi,ia)
!
      do 30 i=1,ia
      if(jind(i).eq.0) call dacop(ma(i),me(i))
  30  continue
!
      call dacct(me,ia,mi,ia,mb,ib)
!
      call dadal(me,ia)
      call dadal(mi,ia)
      call dadal(mn,ia)
!
      return
      end
!
+dk dader
      subroutine dader(idif,ina,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer idif,illc,ilmc,ina,inc,incc,inoc,invc,ipoc,lda,lea,lia,   &
     &lno,lnv,lst
!     ******************************
!
!     THIS SUBROUTINE COMPUTES THE DERIVATIVE WITH RESPECT TO VARIABLE I
!     OF THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall(incc,1,'$$DAJUNK$$',inoc,invc)
        call dadert(idif,ina,incc)
        call dacop(incc,inc)
        call dadal(incc,1)
      else
        call dadert(idif,ina,inc)
      endif

      return
      end

+dk dadert
      subroutine dadert(idif,ina,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ibase,ic,ider1,ider1s,ider2,ider2s,idif,iee,ifac,illa,  &
     &illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,ipoa,ipoc,jj,lda,      &
     &lea,lia,lno,lnv,lst
      double precision rr,x,xdivi
!     ******************************
!
!     THIS SUBROUTINE COMPUTES THE DERIVATIVE WITH RESPECT TO VARIABLE I
!     OF THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1
+ca dabinc
      integer jd(lnv)
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!
      if(nomax.eq.1) then
!         PRINT*,'ERROR, DADER CALLED WITH NOMAX = 1'
!        CALL DADEB(31,'ERR DADER1',1)
+if cr
!     call abend('666                                               ')
+ei
+if .not.cr
!        stop 666
+ei
      do i=1,lnv
        jd(i)=0
        enddo
        jd(idif)=1
        call dapek(ina,jd,rr)
        call dacon(inc,rr)
        return
      endif
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      ibase = nomax + 1
!
      if(idif.gt.(nvmax+1)/2) then
         ider1  = 0
         ider1s = 0
         ider2  = idif-(nvmax+1)/2
         ider2s = 1
         do 10 jj=1,ider2-1
  10     ider2s = ider2s*ibase
         xdivi  = ider2s*ibase
      else
         ider1  = idif
         ider1s = 1
         do 20 jj=1,ider1-1
  20     ider1s = ider1s*ibase
         ider2  = 0
         ider2s = 0
         xdivi  = ider1s*ibase
      endif
!
      ibase = nomax+1
!
      ic = ipoc-1
!
      do 100 i=ipoa,ipoa+illa-1
!
      if(ider1.eq.0) then
         iee = i2(i)
      else
         iee = i1(i)
      endif
!
      x = iee/xdivi
      ifac = int(ibase*(x-int(x+epsmac)+epsmac))
!
      if(ifac.eq.0) goto 100
!
      ic = ic + 1
      cc(ic) = cc(i)*ifac
      i1(ic) = i1(i) - ider1s
      i2(ic) = i2(i) - ider2s
!
 100  continue
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
+if cr
         write(lout,*)'ERROR IN DADER '
+ei
+if .not.cr
         write(*,*)'ERROR IN DADER '
+ei
         call dadeb(31,'ERR DADER2',1)
      endif
!
      return
      end
!
+dk dapoi
      subroutine dapoi(ina,inb,inc,n)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ina,inb,inc,lda,lea,lia,lno,lnv,lst,n
!     *******************************
!
!     THIS SUBROUTINE COMPUTES THE POISSON BRACKET OF THE VECTORS A AND
!     B AND STORES THE RESULT IN C. N IS THE DEGREE OF FREEDOM OF THE SYSTEM.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      integer is(4)
!
      is(1) = 0
      is(2) = 0
      is(3) = 0
      is(4) = 0
      call daall(is,4,'$$DAPOI $$',nomax,nvmax)
!
!
      do 100 i=1,n
!
      call dader(2*i-1,ina,is(1))
      call dader(2*i,  inb,is(2))
      call damul(is(1),is(2),is(3))
      call daadd(is(4),is(3),is(1))
      call dacop(is(1),is(4))
!
      call dader(2*i,  ina,is(1))
      call dader(2*i-1,inb,is(2))
      call damul(is(1),is(2),is(3))
      call dasub(is(4),is(3),is(1))
      call dacop(is(1),is(4))
!
 100  continue

      call dacop(is(4),inc)
!
      call dadal(is,4)
!
      return
      end
!
+dk dacfur
      subroutine dacfur(ina,fun,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illc,ilmc,ina,inc,incc,inoc,invc,ipoc,lda,lea,lia,lno,lnv,&
     &lst
      double complex fun
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall(incc,1,'$$DAJUNK$$',inoc,invc)
        call dacfurt(ina,fun,incc)
        call dacop(incc,inc)
        call dadal(incc,1)
      else
        call dacfurt(ina,fun,inc)
      endif

      return
      end
+dk dacfurt
      subroutine dacfurt(ina,fun,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,  &
     &ipoa,ipoc,j,lda,lea,lia,lno,lnv,lst
      double precision cfac,rr
      double complex fun
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      dimension j(lnv)
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!
      if(nomax.eq.1) then
        do i=1,lnv
        j(i)=0
      enddo
      call dapek(ina,j,rr)
      cfac = dreal(fun(j))
      rr=cfac*rr
      call dapok(inc,j,rr)
      do i=1,lnv
        j(i)=1
        call dapek(ina,j,rr)
        cfac = dreal(fun(j))
        rr=cfac*rr
        call dapok(inc,j,rr)
        j(i)=0
      enddo
!         PRINT*,'ERROR, DACFU CALLED WITH NOMAX = 1'
!        CALL DADEB(31,'ERR DACFU ',1)
+if cr
!     call abend('667                                               ')
+ei
+if .not.cr
!        stop 667
+ei
      return
      endif
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      ic = ipoc - 1
!
      do ia=ipoa,ipoa+illa-1
!
      call dancd(i1(ia),i2(ia),j)
      cfac = dreal(fun(j))
!      IF(dABS(CFAC).LT.EPS) GOTO 100
!      IF(dABS(CFAC*CC(IA)).LT.EPS) GOTO 100
       if(abs(cfac*cc(ia)).lt.eps.or.abs(cc(ia)).lt.eps) goto 100
!
      ic = ic + 1
      cc(ic) = cc(ia)*cfac
      i1(ic) = i1(ia)
      i2(ic) = i2(ia)
!
 100  continue
      enddo
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
+if cr
         write(lout,*)'ERROR IN DACFU '
+ei
+if .not.cr
         write(*,*)'ERROR IN DACFU '
+ei
         call dadeb(31,'ERR DACFU ',1)
      endif
!
      return
      end
!
+dk dacfu
      subroutine dacfu(ina,fun,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illc,ilmc,ina,inc,incc,inoc,invc,ipoc,lda,lea,lia,lno,lnv,&
     &lst
      double precision fun
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE PRECISION FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall(incc,1,'$$DAJUNK$$',inoc,invc)
        call dacfut(ina,fun,incc)
        call dacop(incc,inc)
        call dadal(incc,1)
      else
        call dacfut(ina,fun,inc)
      endif

      return
      end
+dk dacfui
      subroutine dacfui(ina,fun,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illc,ilmc,ina,inc,incc,inoc,invc,ipoc,lda,lea,lia,lno,lnv,&
     &lst
      double complex fun
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall(incc,1,'$$DAJUNK$$',inoc,invc)
        call dacfuit(ina,fun,incc)
        call dacop(incc,inc)
        call dadal(incc,1)
      else
        call dacfuit(ina,fun,inc)
      endif

      return
      end
+dk dacfuit
      subroutine dacfuit(ina,fun,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,  &
     &ipoa,ipoc,j,lda,lea,lia,lno,lnv,lst
      double precision cfac,rr
      double complex fun
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      dimension j(lnv)
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!
      if(nomax.eq.1) then
        do i=1,lnv
        j(i)=0
      enddo
      call dapek(ina,j,rr)
      cfac = dimag(fun(j))
      rr=cfac*rr
      call dapok(inc,j,rr)
      do i=1,lnv
        j(i)=1
        call dapek(ina,j,rr)
        cfac = dimag(fun(j))
        rr=cfac*rr
        call dapok(inc,j,rr)
        j(i)=0
      enddo
!         PRINT*,'ERROR, DACFU CALLED WITH NOMAX = 1'
!        CALL DADEB(31,'ERR DACFU ',1)
+if cr
!     call abend('667                                               ')
+ei
+if .not.cr
!        stop 667
+ei
      return
      endif
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      ic = ipoc - 1
!
      do ia=ipoa,ipoa+illa-1
!
      call dancd(i1(ia),i2(ia),j)
      cfac = dimag(fun(j))
!      IF(dABS(CFAC).LT.EPS) GOTO 100
!      IF(dABS(CFAC*CC(IA)).LT.EPS) GOTO 100
       if(abs(cfac*cc(ia)).lt.eps.or.abs(cc(ia)).lt.eps) goto 100
!
      ic = ic + 1
      cc(ic) = cc(ia)*cfac
      i1(ic) = i1(ia)
      i2(ic) = i2(ia)
!
 100  continue
      enddo
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
+if cr
         write(lout,*)'ERROR IN DACFU '
+ei
+if .not.cr
         write(*,*)'ERROR IN DACFU '
+ei
         call dadeb(31,'ERR DACFU ',1)
      endif
!
      return
      end
!
+dk dacfut
      subroutine dacfut(ina,fun,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,  &
     &ipoa,ipoc,j,lda,lea,lia,lno,lnv,lst
      double precision cfac,fun,rr
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE PRECISION FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      dimension j(lnv)
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!
      if(nomax.eq.1) then
        do i=1,lnv
        j(i)=0
      enddo
      call dapek(ina,j,rr)
      cfac = fun(j)
      rr=cfac*rr
      call dapok(inc,j,rr)
      do i=1,lnv
        j(i)=1
        call dapek(ina,j,rr)
        cfac = fun(j)
        rr=cfac*rr
        call dapok(inc,j,rr)
        j(i)=0
      enddo
!         PRINT*,'ERROR, DACFU CALLED WITH NOMAX = 1'
!        CALL DADEB(31,'ERR DACFU ',1)
+if cr
!     call abend('667                                               ')
+ei
+if .not.cr
!        stop 667
+ei
      return
      endif
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      ic = ipoc - 1
!
      do 100 ia=ipoa,ipoa+illa-1
!
      call dancd(i1(ia),i2(ia),j)
      cfac = fun(j)
!      IF(dABS(CFAC).LT.EPS) GOTO 100
!      IF(dABS(CFAC*CC(IA)).LT.EPS) GOTO 100
       if(abs(cfac*cc(ia)).lt.eps.or.abs(cc(ia)).lt.eps) goto 100
!
      ic = ic + 1
      cc(ic) = cc(ia)*cfac
      i1(ic) = i1(ia)
      i2(ic) = i2(ia)
!
 100  continue
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
+if cr
         write(lout,*)'ERROR IN DACFU '
+ei
+if .not.cr
         write(*,*)'ERROR IN DACFU '
+ei
         call dadeb(31,'ERR DACFU ',1)
      endif
!
      return
      end
!

+dk dapri
      subroutine dapri(ina,iunit)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ii,iii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,iunit,j,k, &
     &lda,lea,lia,lno,lnv,lst
!     ***************************
!       Frank
!     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!-----------------------------------------------------------------------------9
      dimension j(lnv)
+ca daname
!-----------------------------------------------------------------------------3
!
!
      if(ina.lt.1.or.ina.gt.nda) then
+if cr
         write(lout,*)'ERROR IN DAPRI, INA = ',ina
+ei
+if .not.cr
         write(*,*)'ERROR IN DAPRI, INA = ',ina
+ei
!        X = SQRT(-ONE)
!        PRINT*,X
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
         stop 
+ei
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
!
      write(iunit,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)')                       &
     &daname(ina),', NO =',inoa,', NV =',inva,', INA =',ina,            &
     &'***********'//'**********************************'
!
      iout = 0
      ioa = 0

      if(inva.eq.0) then
         write(iunit,'(A)')                                             &
     &'    I  VALUE  '
+if crlibm
+if .not.lf95
                                                  call enable_xp()
+ei
+ei
         do 80 i = ipoa,ipoa+illa-1
         write(iunit,'(I6,2X,G20.14)') i-ipoa, cc(i)
!Eric
 80      write(111) cc(i)
+if crlibm
+if .not.lf95
                                                  call disable_xp()
+ei
+ei
      elseif(nomax.eq.1) then
         if(illa.ne.0) write(iunit,'(A)')                               &
     &'    I  COEFFICIENT          ORDER   EXPONENTS'
         if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
+if crlibm
+if .not.lf95
                                                  call enable_xp()
+ei
+ei
         do 90 i=1,illa
           do k=1,inva
             j(k)=0
           enddo
             iout=iout+1
           if(i.ne.1) then
             j(i-1)=1
             ioa=1
           endif
         write(iunit,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')                 &
     &iout,cc(ipoa+i-1),ioa,(j(iii),iii=1,nvmax)
         write(111) cc(ipoa+i-1)
!Eric
!        write(iunit,*) cc(ipoa+i-1)
         write(iunit,'(G20.14)') cc(ipoa+i-1)
         write(111) cc(ipoa+i-1)
 90     continue
+if crlibm
+if .not.lf95
                                                  call disable_xp()
+ei
+ei
      else
         if(illa.ne.0) write(iunit,'(A)')                               &
     &'    I  COEFFICIENT          ORDER   EXPONENTS'
         if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
      do 100 ioa = 0,inoa
        do 100 ii=ipoa,ipoa+illa-1
          if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
          call dancd(i1(ii),i2(ii),j)
!ETIENNE
          if(abs(cc(ii)).gt.eps) then
!ETIENNE
          iout = iout+1
+if crlibm
+if .not.lf95
                                                  call enable_xp()
+ei
+ei
          write(iunit,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')                &
     &iout,cc(ii),ioa,(j(iii),iii=1,nvmax)
!Eric
         write(111) cc(ii)
!ETIENNE
!Eric
!         write(iunit,* ) cc(ii)
          write(iunit,'(G20.14)') cc(ii)
          write(111) cc(ii)
          endif
!ETIENNE
!
+if crlibm
+if .not.lf95
                                                  call disable_xp()
+ei
+ei
 100  continue
!
      endif

      write(iunit,'(A)') '                                      '
!
!Eric
      write(111) 0d0
      return
      end

+dk dapri77
      subroutine dapri77(ina,iunit)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,iunit,j,lda,   &
     &lea,lia,lno,lnv,lst
      character c10*10,k10*10
!     ***************************
!       Etienne
!     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
!
!-----------------------------------------------------------------------------1
+ca dabinc
      dimension j(lnv)
+ca daname
!-----------------------------------------------------------------------------3
!
        if(iunit.eq.0) return
!
      if(ina.lt.1.or.ina.gt.nda) then
+if cr
         write(lout,*)'ERROR IN DAPRI, INA = ',ina
+ei
+if .not.cr
         write(*,*)'ERROR IN DAPRI, INA = ',ina
+ei
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
         stop 
+ei
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
!
!      WRITE(IUNIT,*) INA, ' in dapri ', DANAME(INA)
!      WRITE(6,*) INA, ' in dapri ', DANAME(INA)
! 611  WRITE(6,*) ' MORE '
!        READ(5,*) MORE
!        IF(MORE.GT.0) THEN
!        WRITE(6,*) MORE,' ',DANAME(MORE)
!        GOTO 611
!        ENDIF
      write(iunit,'(/1X,A10,A6,I5,A6,I5,A7,I5/1X,A/)')                  &
     &daname(ina),', NO =',inoa,', NV =',inva,', INA =',ina,            &
     &'***********'//'**********************************'
!
      if(illa.ne.0) write(iunit,'(A)')                                  &
     &'    I  COEFFICIENT          ORDER   EXPONENTS'
      if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
!
      c10='      NO ='
      k10='      NV ='

      write(iunit,'(A10,I6,A10,I6)') c10,inoa,k10,inva

      iout = 0
!
!      DO 100 IOA = 0,INOA
       do 100 ioa = 0,nocut
       do 100 ii=ipoa,ipoa+illa-1
         if(nomax.ne.1) then
         if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
      endif
!ETIENNE
      if(abs(cc(ii)).gt.eps) then
!ETIENNE

      if(nomax.ne.1) then
        call dancd(i1(ii),i2(ii),j)
        iout = iout+1
      else
        if(ii.eq.ipoa.and.ioa.eq.1) goto 100
        if(ii.gt.ipoa.and.ioa.eq.0) goto 100
        do i=1,lnv
          j(i)=0
        enddo
        if(ii.ne.ipoa) j(ii-ipoa)=1
        iout = iout+1
      endif
!

+if crlibm
+if .not.lf95
                                                  call enable_xp()
+ei
+ei
!      WRITE(IUNIT,*) IOA,CC(II),(J(I),I=1,INVA)
      if(abs(cc(ii)).gt.eps) then
      if(eps.gt.1.e-37) then
       write(iunit,501) ioa,cc(ii),(j(i),i=1,inva)
      else
       write(iunit,503) ioa,cc(ii),(j(i),i=1,inva)
      endif
      endif
 501  format(' ', i3,1x,g23.16,1x,100(1x,i2))
 503  format(' ', i3,1x,g23.16,1x,100(1x,i2))
 502  format(' ', i5,1x,g23.16,1x,100(1x,i2))
+if crlibm
+if .not.lf95
                                                  call disable_xp()
+ei
+ei

      endif
!ETIENNE
!
 100  continue
!
      do 111 i=1,lnv
 111  j(i)=0

      if(iout.eq.0) iout=1
+if crlibm
+if .not.lf95
                                                  call enable_xp()
+ei
+ei

      write(iunit,502) -iout,0.d0,(j(i),i=1,inva)
!
+if crlibm
+if .not.lf95
                                                  call disable_xp()
+ei
+ei
      return
      end

+dk dashift
      subroutine dashift(ina,inc,ishift)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,j,lda,         &
     &lea,lia,lno,lnv,lst
+ca dabinc
!-----------------------------------------------------------------------------9
      dimension j(lnv)
+ca daname
!-----------------------------------------------------------------------------3
!
      integer inb,ishift,ich,ik,jd(lnv),inc
!-----------------------------------------------------------------------------3
!
!

      inb=0
      if(ina.lt.1.or.ina.gt.nda) then
+if cr
         write(lout,*)'ERROR IN DAPRI, INA = ',ina
+ei
+if .not.cr
         write(*,*)'ERROR IN DAPRI, INA = ',ina
+ei
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
         stop 
+ei
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
      call daall(inb,1,'$$DAJUNK$$',inoa,inva)

!      WRITE(IUNIT,*) INA, ' in dapri ', DANAME(INA)
!      WRITE(6,*) INA, ' in dapri ', DANAME(INA)
! 611  WRITE(6,*) ' MORE '
!        READ(5,*) MORE
!        IF(MORE.GT.0) THEN
!        WRITE(6,*) MORE,' ',DANAME(MORE)
!        GOTO 611
!        ENDIF
      iout = 0
!
!      DO 100 IOA = 0,INOA
       do 100 ioa = 0,nocut
       do 100 ii=ipoa,ipoa+illa-1
         if(nomax.ne.1) then
         if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
      endif
!ETIENNE
      if(abs(cc(ii)).gt.eps) then
!ETIENNE

      if(nomax.ne.1) then
        call dancd(i1(ii),i2(ii),j)
        iout = iout+1
      else
        if(ii.eq.ipoa.and.ioa.eq.1) goto 100
        if(ii.gt.ipoa.and.ioa.eq.0) goto 100
        do i=1,lnv
          j(i)=0
        enddo
        if(ii.ne.ipoa) j(ii-ipoa)=1
        iout = iout+1
      endif
!

!      WRITE(IUNIT,*) IOA,CC(II),(J(I),I=1,INVA)
      if(abs(cc(ii)).gt.eps) then
      if(eps.gt.1.e-37) then
+if crlibm
!                                                 call enable_xp()
+ei
!       write(iunit,501) ioa,cc(ii),(j(i),i=1,inva)
+if crlibm
!                                                 call disable_xp()
+ei
!      write(111) cc(ii)
       ich=1
       do ik=1,ishift
         if(j(ik).ne.0) ich=0
       enddo
       if(ich.eq.1) then
         do ik=1,lnv
           jd(ik)=0
         enddo
         do ik=ishift+1,lnv
           jd(ik-ishift)=j(ik)  !%%%%%%etienne
         enddo
       endif
       call dapok(inb,jd,cc(ii))
      else
+if crlibm
!                                                 call enable_xp()
+ei
!       write(iunit,503) ioa,cc(ii),(j(i),i=1,inva)
+if crlibm
!                                                 call disable_xp()
+ei
!       write(111) c(ii)
        ich=1
        do ik=1,ishift
          if(j(ik).ne.0) ich=0
        enddo
        if(ich.eq.1) then
          do ik=1,lnv
            jd(ik)=0
          enddo
          do ik=ishift+1,lnv
            jd(ik-ishift)=j(ik) !%%%%%%etienne
          enddo
        endif
        call dapok(inb,jd,cc(ii))
      endif
      endif
 501  format(' ', i3,1x,g23.16,1x,100(1x,i2))
 503  format(' ', i3,1x,g23.16,1x,100(1x,i2))
 502  format(' ', i5,1x,g23.16,1x,100(1x,i2))

      endif
!ETIENNE
!
 100  continue
!
      do 111 i=1,lnv
 111  j(i)=0

      call dacop(inb,inc)
      call dadal(inb,1)
!
      return
      end

+dk darea
      subroutine darea(ina,iunit)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,iche,ii,ii1,ii2,iin,illa,ilma,ina,inoa,inva,io,io1,  &
     &ipoa,iunit,iwarin,iwarno,iwarnv,j,lda,lea,lia,lno,lnv,lst,nno
      double precision c
!       Frank
+ca dabinc
!-----------------------------------------------------------------------------9
+ca daname
!-----------------------------------------------------------------------------3
!
      character c10*10
      dimension j(lnv)
!
      if(ina.lt.1.or.ina.gt.nda) then
+if cr
         write(lout,*)'ERROR IN DAREA, INA = ',ina
+ei
+if .not.cr
         write(*,*)'ERROR IN DAREA, INA = ',ina
+ei
!        X = SQRT(-ONE)
!        PRINT*,X
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
         stop 
+ei
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
!
      do 5 i=1,lnv
   5  j(i) = 0
!
      call daclr(1)
!
      ic = 0
!
      iwarno = 0
      iwarnv = 0
      iwarin = 0
!
      read(iunit,'(A10)') c10
      read(iunit,'(18X,I4)') nno
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
!
!
      iin = 0
!
  10  continue
+if debug
!     c=0.d0
!     call wda('dar1c',c,1,0,0,0)
+ei
      iin = iin + 1
+if crlibm
+if .not.lf95
                                                  call enable_xp()
+ei
+ei
      read(iunit,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')                     &
     &ii,c,io,(j(i),i=1,inva)
+if crlibm
+if .not.lf95
                                                  call disable_xp()
+ei
+ei
!Eric
      read(111) c
+if debug
!     call wda('dar2c',c,2,0,0,0)
+ei
!
      if(ii.eq.0) goto 20
!ETIENNE
+if debug
!     call wda('dar3c',c,3,0,0,0)
+ei
+if crlibm
+if .not.lf95
                                                  call enable_xp()
+ei
+ei
!Eric
      read(iunit,'(G20.14)') c
!Eric
      read(111) c
+if crlibm
+if .not.lf95
                                                  call disable_xp()
+ei
+ei
+if debug
!     call wda('dar4c',c,4,0,0,0)
+ei
!ETIENNE
      if(ii.ne.iin) then
         iwarin = 1
      endif
      io1 = 0
      do 15 i=1,inva
  15  io1 = io1 + j(i)
!
      if(io1.ne.io) then
         iwarnv = 1
         goto 10
      endif
      if(io.gt.inoa) then
!        IF(IWARNO.EQ.0) PRINT*,'WARNING IN DAREA, FILE ',
!    *              'CONTAINS HIGHER ORDERS THAN VECTOR '
         iwarno = 1
         goto 10
      endif
!
      if(nomax.ne.1) then
        ic = ic + 1
        call dadcd(j,ii1,ii2)
        ic = ia1(ii1) + ia2(ii2)
+if debug
!     call wda('dar5c',c,5,0,0,0)
+ei
        cc(ic) = c
+if debug
!     call wda('dar6c',c,6,0,0,0)
+ei
        goto 10
      else
        iche=0
        do i=1,inva
          if(j(i).eq.1) iche=i
        enddo
+if debug
!     call wda('dar7c',c,7,0,0,0)
+ei
        cc(ipoa+iche)=c
+if debug
!     call wda('dar8c',c,8,0,0,0)
+ei
        goto 10
      endif

!
  20  continue
!
      if(nomax.ne.1) call dapac(ina)
!
+if debug
!     call wda('dar9c',c,9,0,0,0)
+ei
      return
      end
!FF
!
+dk darea77
      subroutine darea77(ina,iunit)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,iche,ii,ii1,ii2,iin,illa,ilma,ina,inoa,inva,ipoa,    &
     &iunit,j,k,lda,lea,lia,lno,lnv,lst,nojoh,nvjoh
      double precision c
!     ***************************
!     Etienne
!     THIS SUBROUTINE READS THE DA VECTOR INA FROM UNIT IUNIT.
!
!-----------------------------------------------------------------------------1
+ca dabinc
+ca daname
!-----------------------------------------------------------------------------3
!
      character c10*10,k10*10
      dimension j(lnv)
!
      if(ina.lt.1.or.ina.gt.nda) then
+if cr
         write(lout,*)'ERROR IN DAPRI, INA = ',ina
+ei
+if .not.cr
         write(*,*)'ERROR IN DAPRI, INA = ',ina
+ei
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
         stop 
+ei
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
!
      do 5 i=1,lnv
   5  j(i) = 0
!
      call daclr(1)
!
!
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10,I6,A10,I6)') c10,nojoh,k10,nvjoh
!
      iin = 0
!
  10  continue
      iin = iin + 1
      read(iunit,*) ii,c,(j(k),k=1,nvjoh)
      if(ii.lt.0) goto 20

      do 15 i=inva+1,nvjoh
        if(j(i).ne.0) goto 10
  15  continue
      iche=0
      do 16 i=1,inva
        iche=iche+j(i)
  16  continue
      if(iche.gt.nomax) goto 10
      if(nomax.ne.1) then
        call dadcd(j,ii1,ii2)
        ic = ia1(ii1) + ia2(ii2)
        cc(ic) = c
      else
        iche=0
        do i=1,inva
          if(j(i).eq.1) iche=i
        enddo
        cc(ipoa+iche)=c
      endif
      goto 10
!
  20  continue
!
      if(nomax.ne.1) call dapac(ina)
!
      return
      end

+dk dadeb
      subroutine dadeb(iunit,c,istop)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer istop,iunit,lda,lea,lia,lno,lnv,lst
!     *******************************
!
!     THIS SUBROUTINE SERVES AS A DEBUGGING TOOL. IT PRINTS ALL
!     NONZERO INFORMATION IN THE COMMON BLOCKS AND ALL DA  VECTORS.
!
!-----------------------------------------------------------------------------1
+ca dabinc
+ca daname
!-----------------------------------------------------------------------------3
!
      character c*10
!
!etienne

+if cr
      write(lout,*) '  ',c
+ei
+if .not.cr
      write(*,*) '  ',c
+ei
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
      stop 
+ei
      end
!
!
!
+dk danum
      subroutine danum(no,nv,numda)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,mm,no,numda,nv
!     *****************************
!
!     THIS SUBROUTINE COMPUTES THE NUMBER OF MONOMIALS OF
!     ORDER NO AND NUMBER OF VARIABLES NV
!
      numda = 1
      mm = max(nv,no)
!
      do 5 i=1,min(nv,no)
  5   numda = (numda*(mm+i))/i
!
      return
      end
!
+dk dainf
      subroutine dainf(inc,inoc,invc,ipoc,ilmc,illc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illc,ilmc,inc,inoc,invc,ipoc,lda,lea,lia,lno,lnv,lst
!     **********************************************
!
!     THIS SUBROUTINE SEARCHES THE NUMBER OF DA VECTOR C
!     AND RETURS THE INFORMATION IN COMMON DA
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      if(inc.ge.1.and.inc.le.nda) then
         inoc = idano(inc)
         invc = idanv(inc)
         ipoc = idapo(inc)
         ilmc = idalm(inc)
         illc = idall(inc)
         return
      endif
!
+if cr
      write(lout,*) 'ERROR IN DAINF, DA VECTOR ',inc,' NOT FOUND '
+ei
+if .not.cr
      write(*,*) 'ERROR IN DAINF, DA VECTOR ',inc,' NOT FOUND '
+ei
      call dadeb(31,'ERR DAINF ',1)
!
      return
      end
!
+dk dapac
      subroutine dapac(inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,illc,ilmc,inc,inoc,invc,ipoc,lda,lea,lia,lno,lnv,lst
      double precision ccc
!     ************************
!
!     THIS SUBROUTINE PACKS THE INFORMATION IN THE SCRATCH VECTOR 1
!     INTO THE VECTOR INC. IF LF = 1, THE FILTERING (CF DAMUF) IS
!     PERFORMED.
!     INVERSE IS DAUNP.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
      ic = ipoc - 1
!
         do 100 i=1,nmmax
         ccc = cc(i)
         if(abs(ccc).lt.eps) goto 100
         ic = ic + 1
         cc(ic) = ccc
         i1(ic) = ie1(i)
         i2(ic) = ie2(i)
 100     continue
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
+if cr
         write(lout,*)'ERROR IN DAPAC '
+ei
+if .not.cr
         write(*,*)'ERROR IN DAPAC '
+ei
         call dadeb(31,'ERR DAPAC ',1)
      endif
!
      return
      end
!
!
+dk dachk
      subroutine dachk(ina,inoa,inva, inb,inob,invb, inc,inoc,invc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ierr,ina,inb,inc,inoa,inob,inoc,inva,invb,invc,invsum,lsw
!     *************************************************************
!
!     THIS SUBROUTINE CHECKS IF THE VECTORS A, B AND C
!     HAVE COMPATIBLE ATTRIBUTES
!
      parameter(lsw=1)
!
      if(lsw.eq.1) return
!
      ierr = 0
!
!     CASE OF A UNARY OPERATION
!     *************************
!
      if(inob.eq.-1.and.invb.eq.-1) then
         invsum = inva + invc
         if(invsum.eq.0) then
            if(inoa.gt.inoc) ierr = 1
         elseif(invsum.eq.1) then
            ierr = 1
         else
            if(inoa.gt.inoc.or.inva.gt.invc) ierr = 1
         endif
         if(ierr.eq.1) then
+if cr
            write(lout,*)'ERROR IN DACHK, ',ina,' AND ',inc,
+ei
+if .not.cr
            write(*,*)'ERROR IN DACHK, ',ina,' AND ',inc,               &
+ei
     &' ARE INCOMPATIBLE',inoa,inva,inoc,invc
            call dadeb(31,'ERR DACHK1',1)
         endif
!
!     CASE OF A BINARY OPERATION
!     **************************
!
      else
         invsum = inva + invb + invc
         if(invsum.eq.0) then
            if(inoa.gt.inoc.or.inob.gt.inoc) ierr = 1
         elseif(invsum.eq.1.or.invsum.eq.2) then
            ierr = 1
         else
            if(inoa.gt.inoc.or.inob.gt.inoc.or.                         &
     &inva.gt.invc.or.invb.gt.invc) ierr = 1
         endif
         if(ierr.eq.1) then
+if cr
            write(lout,*)'ERROR IN DACHK, ',ina,',',inb,' AND ',inc,
+ei
+if .not.cr
            write(*,*)'ERROR IN DACHK, ',ina,',',inb,' AND ',inc,       &
+ei
     &' ARE INCOMPATIBLE'
            call dadeb(31,'ERR DACHK2',1)
         endif
      endif
!
      return
      end
!
+dk damch
      subroutine damch(iaa,ia)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,iaa,illa,ilma,ino1,inoi,inv1,invi,ipoa
!     ************************
!
!     THIS SUBROUTINE CHECKS IF THE IA VECTORS IN THE MATRIX IA HAVE
!     IDENTICAL ATTRIBUTES.
!
      dimension iaa(*)
!
      call dainf(iaa(1),ino1,inv1,ipoa,ilma,illa)
!
      do 10 i=2,ia
      call dainf(iaa(i),inoi,invi,ipoa,ilma,illa)
      if(ino1.ne.inoi.or.inv1.ne.invi) then
+if cr
         write(lout,*)'ERROR IN DAMCH, VECTORS ',iaa(1),' AND ',iaa(i),
+ei
+if .not.cr
         write(*,*)'ERROR IN DAMCH, VECTORS ',iaa(1),' AND ',iaa(i),    &
+ei
     &' ARE INCOMPATIBLE '
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
         stop 
+ei
      endif
  10  continue
!
      return
      end
!
+dk dadcd
      subroutine dadcd(jj,ic1,ic2)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ibase,ic1,ic2,isplit,lda,lea,lia,lno,lnv,lst
!     ****************************
!
!     THIS SUBROUTINE CODES THE EXPONENTS IN JJ INTO THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      integer jj(lnv)
      ibase = nomax + 1
      isplit = (nvmax+1)/2
      ic1 = 0
      ic2 = 0
!
      do 10 i=nvmax,isplit+1,-1
 10   ic2 = ic2*ibase + jj(i)
!
      do 20 i=isplit,1,-1
 20   ic1 = ic1*ibase + jj(i)
!
      return
      end
!
+dk dancd
      subroutine dancd(ic1,ic2,jj)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ibase,ic,ic1,ic2,isplit,lda,lea,lia,lno,lnv,lst
      double precision x
!     ****************************
!
!     THIS SUBROUTINE ENCODES THE EXPONENTS IN JJ FROM THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      integer jj(*)
      ibase = nomax + 1
      isplit = (nvmax+1)/2
!
      ic = ic1
      do 10 i=1,isplit
      x  = dble(ic)/dble(ibase)
      ic = int(x+epsmac)
  10  jj(i) = nint(ibase*(x-ic))
!
      ic = ic2
      do 20 i=isplit+1,nvmax
      x  = dble(ic)/dble(ibase)
      ic = int(x+epsmac)
  20  jj(i) = nint(ibase*(x-ic))
!
      do 30 i=nvmax+1,lnv
  30  jj(i) = 0
!
      return
      end

!ETIENNE
+dk datra
      subroutine datra(idif,ina,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ibase,ic,ider1,ider1s,ider2,ider2s,idif,iee,ifac,illa,  &
     &illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,ipoa,ipoc,jj,lda,      &
     &lea,lia,lno,lnv,lst
      double precision x,xdivi
!     ******************************
!
!     THIS SUBROUTINE COMPUTES THE PSEUDO DERIVATIVE WITH RESPECT TO VARIABLE I
!     OF THE VECTOR A AND STORES THE RESULT IN C.
!
!     dx^n/dx= x^(n-1)
!-----------------------------------------------------------------------------1
+ca dabinc
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!
!       CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!

      if(nomax.eq.1) then
        call dader(idif,ina,inc)
        return
      endif
      ibase = nomax + 1
!
      if(idif.gt.(nvmax+1)/2) then
         ider1  = 0
         ider1s = 0
         ider2  = idif-(nvmax+1)/2
         ider2s = 1
         do 10 jj=1,ider2-1
  10     ider2s = ider2s*ibase
         xdivi  = ider2s*ibase
      else
         ider1  = idif
         ider1s = 1
         do 20 jj=1,ider1-1
  20     ider1s = ider1s*ibase
         ider2  = 0
         ider2s = 0
         xdivi  = ider1s*ibase
      endif
!
      ibase = nomax+1
!
      ic = ipoc-1
!
      do 100 i=ipoa,ipoa+illa-1
!
      if(ider1.eq.0) then
         iee = i2(i)
      else
         iee = i1(i)
      endif
!
      x = iee/xdivi
      ifac = int(ibase*(x-int(x+epsmac)+epsmac))
!
      if(ifac.eq.0) goto 100
!
!etienne      IFAC = INT(IBASE*(X-INT(X)+1.D-8))
!
!etienne      IF(IFAC.EQ.0) GOTO 100
!
      ic = ic + 1
      cc(ic) = cc(i)
      i1(ic) = i1(i) - ider1s
      i2(ic) = i2(i) - ider2s
!
 100  continue
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
+if cr
         write(lout,*)'ERROR IN DADTRA'
+ei
+if .not.cr
         write(*,*)'ERROR IN DADTRA'
+ei
         call dadeb(111,'ERR DADTRA',1)
      endif
!
      return
      end

+dk etred
      subroutine etred(no1,nv1,ic1,ic2,no2,nv2,i11,i21)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,i11,i21,ic,ic1,ic2,lda,lea,lia,lno,lnv,lst,no1,no2,nv1, &
     &nv2
!     ****************************
!
!     THIS SUBROUTINE CODES THE EXPONENTS IN JJ INTO THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      integer jj(lnv)

      if(nv1.gt.lnv.or.nv2.gt.lnv) then
+if cr
      write(lout,*) ' ERROR IN RECODING '
+ei
+if .not.cr
      write(*,*) ' ERROR IN RECODING '
+ei
+if cr
      call abend('123                                               ')
+ei
+if .not.cr
      stop 123
+ei
      endif
      call dehash(no1,nv1,ic1,ic2,jj)
      ic=0

      do 1 i=nv1+1,nv2
      if(jj(i).gt.0) then
      i11=-1
      i21=0
      return
      endif
 1    continue

      do 2 i=1,nv2
      ic=ic+jj(i)
 2    continue

      if(ic.gt.no2) then
      i11=-1
      i21=0
      return
      endif

      call hash(no2,nv2,jj,i11,i21)

!
      return
      end

+dk hash
      subroutine hash(no1,nv1,jj,ic1,ic2)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ibase,ic1,ic2,isplit,lda,lea,lia,lno,lnv,lst,no1,nv1
!     ****************************
!
!     THIS SUBROUTINE CODES THE EXPONENTS IN JJ INTO THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1
+ca dabinc
      integer jj(*)

      ibase = no1 + 1
      isplit = (nv1+1)/2
      ic1 = 0
      ic2 = 0
!
      do 10 i=nv1,isplit+1,-1
 10   ic2 = ic2*ibase + jj(i)
!
      do 20 i=isplit,1,-1
 20   ic1 = ic1*ibase + jj(i)
!
      return
      end
!
+dk dehash
      subroutine dehash(no1,nv1,ic1,ic2,jj)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ibase,ic,ic1,ic2,isplit,lda,lea,lia,lno,lnv,lst,        &
     &no1,nv1
      double precision x
!     ****************************
!
!     THIS SUBROUTINE ENCODES THE EXPONENTS IN JJ FROM THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      integer jj(*)

      epsmac=1.e-7
      ibase = no1 + 1
      isplit = (nv1+1)/2
!
      ic = ic1
      do 10 i=1,isplit
      x  = dble(ic)/dble(ibase)
      ic = int(x+epsmac)
  10  jj(i) = nint(ibase*(x-ic))
!
      ic = ic2
      do 20 i=isplit+1,nv1
      x  = dble(ic)/dble(ibase)
      ic = int(x+epsmac)
  20  jj(i) = nint(ibase*(x-ic))
!
      return
      end

+dk daswap
      subroutine daswap(j1,j2,inb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ia,ic,ic1,ic2,illb,ilmb,inb,inob,invb,ipob,j1,j2,jj,k1,k2,&
     &lda,lea,lia,lno,lnv,lst
!     *************************
!
!     SWAP A DA VECTOR
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      dimension jj(lnv)

      call dainf(inb,inob,invb,ipob,ilmb,illb)

      call daclr(1)
!

      do 100 ia = ipob,ipob+illb-1

      call dehash(nomax,nvmax,i1(ia),i2(ia),jj)
      k1=jj(j1)
      k2=jj(j2)
      jj(j1)=k2
      jj(j2)=k1
      call hash(nomax,nvmax,jj,ic1,ic2)

      ic=ia1(ic1)+ia2(ic2)

      cc(ic) = cc(ia)
      i1(ic) = ic1
      i2(ic) = ic2
!
 100  continue
!
!
      call dapac(inb)
      return
      end

+dk dagauss
      subroutine dagauss(ina,inb,nd2,anorm)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,ib,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,  &
     &ipoa,ipob,ja,jb,lda,lea,lia,lno,lnv,lst,nd2
      double precision anorm,gau
!     ***************************
!
!     THIS SUBROUTINE COMPUTES THE NORM OF THE DA VECTOR A
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      dimension ja(lnv),jb(lnv)

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
      anorm = 0.d0

      do 100 ia=ipoa,ipoa+illa-1
      do 101 ib=ipob,ipob+illb-1
      call dancd(i1(ia),i2(ia),ja)
      call dancd(i1(ib),i2(ib),jb)
      gau=1.d0
      do 102 i=1,nd2
 102  gau= facint(ja(i)+jb(i))*gau
      anorm = anorm + cc(ia)*cc(ib)*gau
 101  continue
 100  continue

!
      return
      end

+dk daran
      subroutine daran(ina,cm,xran)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,illa,ilma,ina,inoa,inva,ipoa,lda,lea,lia,lno,lnv,       &
     &lst
      double precision bran,cm,xran
!     ************************
!
!     THIS SUBROUTINE FILLS THE DA VECTOR A WITH RANDOM ENTRIES.
!     FOR CM > 0, THE VECTOR IS FILLED WITH REALS,
!     FOR CM < 0, THE VECTOR IS FILLED WITH SINGLE DIGIT INTEGERS
!     ABS(CM) IS THE FILLING FACTOR
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!
      if(inva.eq.0.or.nomax.eq.1) then
         do 10 i=ipoa,ipoa+ilma-1
         if(cm.gt.0d0) then
            cc(i) = bran(xran)
            if(cc(i).gt.cm) cc(i) = 0d0
         elseif(cm.lt.0d0) then
            cc(i) = int(1+10*bran(xran))
            if(cc(i).gt.-1d1*cm) cc(i) = 0d0
         endif
  10     continue
         idall(ina) = idalm(ina)
         return
      endif
!
      if(inoa.ne.nomax.or.inva.ne.nvmax) then
+if cr
         write(lout,*)'ERROR IN DARAN, ONLY VECTORS WITH NO = NOMAX AND'
+ei
+if .not.cr
         write(*,*)'ERROR IN DARAN, ONLY VECTORS WITH NO = NOMAX AND'   &
+ei
     &//' NV = NVMAX ALLOWED'
         call dadeb(31,'ERR DARAN1',1)
      endif
!
      call daclr(1)
!
      do 100 i=1,nmmax
      if(cm.gt.0.d0) then
         cc(i) = bran(xran)
         if(cc(i).gt.cm) cc(i) = 0.d0
      elseif(cm.lt.0.d0) then
         cc(i) = int(1+10*bran(xran))
         if(cc(i).gt.-10.d0*cm) cc(i) = 0.d0
      else
+if cr
         write(lout,*)'ERROR IN ROUTINE DARAN'
+ei
+if .not.cr
         write(*,*)'ERROR IN ROUTINE DARAN'
+ei
         call dadeb(31,'ERR DARAN2',1)
      endif
 100  continue
!
      call dapac(ina)
!
      return
      end
!
+dk bran
      double precision function bran(xran)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      double precision xran
!     ************************************
!
!     VERY SIMPLE RANDOM NUMBER GENERATOR
!
      xran = xran + 10.d0
      if(xran.gt.1.d4) xran = xran - 9999.12345
+if crlibm
      bran = abs(sin_rn(xran))
+ei
+if .not.crlibm
      bran = abs(sin(xran))
+ei
      bran = 10.0D0*bran
      bran = bran - int(bran)
!      IF(BRAN.LT. .1D0) BRAN = BRAN + .1D0
!
      return
      end
!
!

+dk danorm2
      subroutine danorm2(ina,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illc,ilmc,ina,inc,incc,inoc,invc,ipoc,lda,lea,lia,lno,lnv,&
     &lst
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall(incc,1,'$$DAJUNK$$',inoc,invc)
        call danorm2t(ina,incc)
        call dacop(incc,inc)
        call dadal(incc,1)
      else
        call danorm2t(ina,inc)
      endif

      return
      end

+dk danorm2t
      subroutine danorm2t(ina,inb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ia,ib,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,    &
     &ipoa,ipob,lda,lea,lia,lno,lnv,lst
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)
!
!
      ib = ipob - 1
!
      do 100 ia=ipoa,ipoa+illa-1
!
      if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 100
      ib = ib + 1
      cc(ib) = cc(ia)**2
      i1(ib) = i1(ia)
      i2(ib) = i2(ia)
!
 100  continue
!
      idall(inb) = ib-ipob+1
      if(idall(inb).gt.idalm(inb)) then
+if cr
         write(lout,*)'ERROR IN DANORM'
+ei
+if .not.cr
         write(*,*)'ERROR IN DANORM'
+ei
         call dadeb(31,'ERR DANOR1',1)
      endif
!
      return
      end

+dk danormr
      subroutine danormr(ina,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer illc,ilmc,ina,inc,incc,inoc,invc,ipoc,lda,lea,lia,lno,lnv,&
     &lst
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall(incc,1,'$$DAJUNK$$',inoc,invc)
        call danormrt(ina,incc)
        call dacop(incc,inc)
        call dadal(incc,1)
      else
        call danormrt(ina,inc)
      endif

      return
      end

+dk danormrt
      subroutine danormrt(ina,inb)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ia,ib,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,    &
     &ipoa,ipob,lda,lea,lia,lno,lnv,lst
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)
!
!
      ib = ipob - 1
!
      do 100 ia=ipoa,ipoa+illa-1
!
      if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 100
      ib = ib + 1
      cc(ib) = dsqrt(cc(ia))
      i1(ib) = i1(ia)
      i2(ib) = i2(ia)
!
 100  continue
!
      idall(inb) = ib-ipob+1
      if(idall(inb).gt.idalm(inb)) then
+if cr
         write(lout,*)'ERROR IN DANORM '
+ei
+if .not.cr
         write(*,*)'ERROR IN DANORM '
+ei
         call dadeb(31,'ERR DANOR2',1)
      endif
!
      return
      end
+dk dakey
      subroutine dakey(c)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      character c*(*)
!

      return
!
      end
! ANFANG UNTERPROGRAMM
+dk dapri6
      subroutine dapri6(ina,result,ien,i56)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,i56,ien,ihp,ii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,   &
     &j,lda,lea,lia,lno,lnv,lst
      double precision result
!     *************************************
!
!     THIS SUBROUTINE IS FOR REDUCED STORAGE DA VERSION JULY 91
!     RESULT CONTAINS THE (IEN-1) TH DERIVATIVE TO ENERGY AFTER EXECUTION
!     I56 SAYS WHETHER THE 5TH OR THE 6TH COORDINATE IS THE ENERGY
!     AND MUST HAVE THE VALUE 5 OR 6 ACCORDINGLY
!-----------------------------------------------------------------------------1
+ca dabinc
+ca daname
!-----------------------------------------------------------------------------3
      dimension j(lnv)
      result=0.
      if(ina.lt.1.or.ina.gt.nda) then
+if cr
        write(lout,*)'ERROR IN DAPRI6, INA = ',ina
+ei
+if .not.cr
        write(*,*)'ERROR IN DAPRI6, INA = ',ina
+ei
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
        stop 
+ei
      endif
      inoa=idano(ina)
      inva=idanv(ina)
      ipoa=idapo(ina)
      ilma=idalm(ina)
      illa=idall(ina)
      iout=0
      if(nomax.eq.1) then
        do 90 i=1,illa
          if(ien.eq.1) then
            if(i-1.ne.0) goto 90
            result=cc(ipoa+i-1)
            return
          endif
          if(ien.eq.2) then
            if(i-1.ne.i56) goto 90
            result=cc(ipoa+i-1)
            return
          endif
 90     continue
      else
        do 100 ioa=0,inoa
        do 100 ii=ipoa,ipoa+illa-1
          if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
          iout=iout+1
          call dancd(i1(ii),i2(ii),j)
          if(i56.eq.6) then
            do ihp=1,5
              if(j(ihp).ne.0) goto 100
            enddo
            if(j(6).eq.(ien-1)) then
              result=cc(ii)
              return
            endif
          else if(i56.eq.5) then
            do ihp=1,4
              if(j(ihp).ne.0) goto 100
            enddo
            if(j(5).eq.(ien-1)) then
              result=cc(ii)
              return
            endif
          else if(i56.eq.4) then
            do ihp=1,3
              if(j(ihp).ne.0) goto 100
            enddo
            if(j(4).eq.(ien-1)) then
              result=cc(ii)
              return
            endif
          endif
 100    continue
      endif
      return
      end
! ANFANG UNTERPROGRAMM
+dk darea6
      subroutine darea6(ina,zfeld,i56)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,i56,ic,ii1,ii2,iin,illa,ilma,ina,inoa,inva,io,io1,ip,   &
     &ipoa,iwarin,iwarno,iwarnv,j,lda,lea,lia,lno,lnv,lst
      double precision zfeld
!     *************************************
!
!     THIS SUBROUTINE IS FOR REDUCED STORAGE DA VERSION JULY 91
!     I56 SAYS WHETHER THE 5TH OR THE 6TH COORDINATE IS THE ENERGY
!     AND MUST HAVE THE VALUE 5 OR 6 ACCORDINGLY
!
!-----------------------------------------------------------------------------1
+ca dabinc
+ca daname
      dimension zfeld(100)
!-----------------------------------------------------------------------------3
      dimension j(lnv)
      if(ina.lt.1.or.ina.gt.nda) then
+if cr
        write(lout,*)'ERROR IN DAREA6, INA = ',ina
+ei
+if .not.cr
        write(*,*)'ERROR IN DAREA6, INA = ',ina
+ei
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
        stop 
+ei
      endif
      inoa=idano(ina)
      inva=idanv(ina)
      ipoa=idapo(ina)
      ilma=idalm(ina)
      illa=idall(ina)
      do 5 i=1,lnv
        j(i)=0
   5  continue
      call daclr(1)
      ic=0
      iwarno=0
      iwarnv=0
      iwarin=0
      iin=0
      if(nomax.eq.1) then
        do 90 i=1,illa
          if (i-1.eq.0) then
            cc(ipoa+i-1)=zfeld(1)
          else if (i-1.eq.i56) then
            cc(ipoa+i-1)=zfeld(2)
          endif
  90    continue
        return
      endif
      do 8 ip=1,inva
        j(ip)=0
   8  continue
      io=0
  10  continue
      iin=iin+1
      io1=0
      do 15 i=1,inva
        io1=io1+j(i)
  15  continue
      if(io1.ne.io) then
+if cr
        if(iwarnv.eq.0) write(lout,*)'WARNING IN DAREA6, FILE ',
+ei
+if .not.cr
        if(iwarnv.eq.0) write(*,*)'WARNING IN DAREA6, FILE ',           &
+ei
     &'CONTAINS MORE VARIABLES THAN VECTOR'
        iwarnv = 1
        goto 10
      endif
      if(io.gt.inoa) then
        iwarno = 1
        goto 10
      endif
      ic = ic + 1
      call dadcd(j,ii1,ii2)
      ic = ia1(ii1) + ia2(ii2)
      cc(ic) = zfeld(io+1)
      j(i56)=j(i56)+1
      io=io+1
      if (io.gt.inoa) goto 20
      goto 10
  20  continue
      call dapac(ina)
      return
      end
! ANFANG FUNKTION
+dk dare
      double precision function dare(ina)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ii,illa,ilma,ina,inoa,inva,ioa,ipoa,j,jj,lda,lea,lia,lno, &
     &lnv,lst
!     ***********************************
!     NEW VERSION OF DARE, AUGUST 1992
!     SUPPOSED TO TREAT THE 0TH COMPONENT ACCURATELY
!
!     30.10 1997 E.Mcintosh & F.Schmidt
!
!-----------------------------------------------------------------------------1
+ca dabinc
      dimension j(lnv)
!-----------------------------------------------------------------------------9
!
!      CALL DAINF(INA,INOA,INVA,IPOA,ILMA,ILLA)
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)

!FRS 30.10.1997
      if(nomax.eq.1) then
        dare = cc(ipoa)
        return
      endif
!FRS 30.10.1997
!FRS March 1997
!      IF(NOMAX.EQ.1) goto 110
!FRS March 1997

      ioa = 0
      do 100 ii=ipoa,ipoa+illa-1
      if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
      call dancd(i1(ii),i2(ii),j)
      do 110 jj=1,inva
      if(j(jj).ne.0) goto 100
 110  continue
      dare = cc(ipoa)
      return
 100  continue
      dare = 0d0
      return
      end
! ANFANG UNTERPROGRAMM
+dk daprimax
      subroutine daprimax(ina,iunit)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ii,iii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,iunit,j,   &
     &lda,lea,lia,lno,lnv,lst
!     ***************************
!
!     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!-----------------------------------------------------------------------------9
+ca daname
!-----------------------------------------------------------------------------3
!
      dimension j(lnv)
!
      if(ina.lt.1.or.ina.gt.nda) then
+if cr
         write(lout,*)'ERROR IN DAPRI, INA = ',ina
+ei
+if .not.cr
         write(*,*)'ERROR IN DAPRI, INA = ',ina
+ei
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
         stop 
+ei
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
!
      iout = 0
      ioa = 0

      if(inva.eq.0) then
+if crlibm
+if .not.lf95
                                                  call enable_xp()
+ei
+ei
         do 80 i = ipoa,ipoa+illa-1
 80      write(iunit,'(I6,2X,G20.14)') i-ipoa, cc(i)
+if crlibm
+if .not.lf95
                                                  call disable_xp()
+ei
+ei
      elseif(nomax.eq.1) then
         do 90 i=1,illa
             iout=iout+1
           if(i.ne.1) then
             j(i-1)=1
             ioa=1
           endif
+if crlibm
+if .not.lf95
                                                  call enable_xp()
+ei
+ei
         write(iunit,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')                 &
     &iout,cc(ipoa+i-1),ioa,(j(iii),iii=1,nvmax)
         write(iunit,*) cc(ipoa+i-1)
+if crlibm
+if .not.lf95
                                                  call disable_xp()
+ei
+ei
 90      continue
      else
        iout = 0
        do 100 ioa = 0,inoa
        do 100 ii=ipoa,ipoa+illa-1
          if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
          call dancd(i1(ii),i2(ii),j)
!ETIENNE
          if(abs(cc(ii)).gt.eps) then
!ETIENNE
          iout = iout+1
+if crlibm
+if .not.lf95
                                                  call enable_xp()
+ei
+ei
          write(iunit,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')                &
     &iout,cc(ii),ioa,(j(iii),iii=1,nvmax)
!ETIENNE
          write(iunit,*) cc(ii)
+if crlibm
+if .not.lf95
                                                  call disable_xp()
+ei
+ei
          endif
!ETIENNE
!
 100  continue
      endif
!

!     WRITE(IUNIT,'(A)') '                                      '
!
      return
      end
!FF

!  unknown stuff
+dk damono
      subroutine damono(ina,jd,cfac,istart,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ia,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,    &
     &ipoa,ipoc,istart,jd,lda,lea,lia,lno,lnv,lst
      double precision cfac
!     *****************************
!
!     THIS SUBROUTINE RETURNS THE MONOMIALS ONE BY ONE
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      dimension jd(*)
!
      if(ina.eq.inc) then
+if cr
       write(lout,*) ' USE DIFFERENT POWER SERIES IN DAMONO '
+ei
+if .not.cr
       write(*,*) ' USE DIFFERENT POWER SERIES IN DAMONO '
+ei
+if cr
      call abend('999                                               ')
+ei
+if .not.cr
       stop 999
+ei
      endif
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
      if(istart.eq.0) then
        istart=illa
        return
      endif
!
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      ic = ipoc - 1
!
      ia=ipoa+istart-1
!
      call dancd(i1(ia),i2(ia),jd)

!
      ic = ic + 1
      cc(ic) = cc(ia)
       cfac=cc(ia)
      i1(ic) = i1(ia)
      i2(ic) = i2(ia)
 100  continue
!
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
+if cr
         write(lout,*)'ERROR IN DAMONO'
+ei
+if .not.cr
         write(*,*)'ERROR IN DAMONO'
+ei
         call dadeb(31,'ERR DAMONO',1)
      endif
!
      return
      end
!
!

+dk dacycle
      subroutine dacycle(ina,ipresent,value,j,illa)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ii,illa,ilma,ina,inoa,inva,iout,ipoa,ipresent,j,lda,lea,&
     &lia,lno,lnv,lst
      double precision value
!     ***************************
!
!     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
!
!-----------------------------------------------------------------------------1
+ca dabinc
      dimension j(lnv)
+ca daname
!-----------------------------------------------------------------------------3
!
!
      if(ina.lt.1.or.ina.gt.nda) then
+if cr
        write(lout,*)'ERROR IN DAPRI, INA = ',ina
+ei
+if .not.cr
        write(*,*)'ERROR IN DAPRI, INA = ',ina
+ei
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
        stop 
+ei
      endif
!
      if(ina.eq.0) then
        value=0.d0
        illa=0
        do 111 i=1,lnv
 111    j(i)=0
        return
      endif

      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
      iout = 0
      ipresent=1+ipresent
        if(ipresent.gt.illa) then
          ipresent=1
        endif
        ii=ipresent+ipoa-1
        call dancd(i1(ii),i2(ii),j)
        value=cc(ii)
      return

      end
+dk daorder
      subroutine daorder(ina,iunit,jx,invo,nchop)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,ii,ii1,ii2,iin,illa,ilma,ina,inoa,inva,invo,io,io1,  &
     &ipoa,iunit,iwarin,iwarno,iwarnv,j,jh,jt,jx,lda,lea,lia,lno,lnv,   &
     &lst,nchop
      double precision c
!     ***************************
!
!     THIS SUBROUTINE READS THE DA VECTOR INA FROM UNIT IUNIT.
!
!-----------------------------------------------------------------------------1
+ca dabinc
+ca daname
!-----------------------------------------------------------------------------3
!
      character c10*10
      dimension j(lnv),jx(lnv),jt(lnv)
!
      if(ina.lt.1.or.ina.gt.nda) then
+if cr
         write(lout,*)'ERROR IN DAPRI, INA = ',ina
+ei
+if .not.cr
         write(*,*)'ERROR IN DAPRI, INA = ',ina
+ei
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
         stop 
+ei
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
!
      do 5 i=1,lnv
      jt(i)=0
   5  j(i) = 0
!
      call daclr(1)
!
      ic = 0
!
      iwarno = 0
      iwarnv = 0
      iwarin = 0
!
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
!
      iin = 0
!
  10  continue
      iin = iin + 1
+if crlibm
+if .not.lf95
                                                  call enable_xp()
+ei
+ei
      read(iunit,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')                     &
     &ii,c,io,(jt(i),i=1,invo)
+if crlibm
+if .not.lf95
                                                  call disable_xp()
      read(111) c
+ei
+ei
!
      if(ii.eq.0) goto 20
!etienne
!Eric
!     read(iunit,*) c
+if crlibm
+if .not.lf95
                                                  call enable_xp()
+ei
+ei
      read(iunit,'(G20.14)') c
+if crlibm
+if .not.lf95
                                                  call disable_xp()
+ei
+ei
      read(111) c
      do 999 jh=1,invo
 999  j(jh)=jt(jx(jh))
      do 998 jh=nchop+1,inva
 998  j(jh)=0
!etienne
      io1 = 0
      do 15 i=1,inva
  15  io1 = io1 + j(i)
!
      ic = ic + 1
      call dadcd(j,ii1,ii2)
      ic = ia1(ii1) + ia2(ii2)
      cc(ic) = c
      goto 10
!
  20  continue
!
      call dapac(ina)
!
      return
      end
!
!ETIENNE
+dk datrash
      subroutine datrash(idif,ina,inc)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ibase,ic,ider1,ider1s,ider2,ider2s,idif,ikil1,ikil2,    &
     &illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,ipoa,ipoc,jj,     &
     &lda,lea,lia,lno,lnv,lst
      double precision xdivi
!     ******************************
!
!     THIS SUBROUTINE COMPUTES THE DERIVATIVE WITH RESPECT TO VARIABLE I
!     OF THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1
+ca dabinc
!
      integer jx(lnv)

!      call daclr(1)
!      call dacop(ina,1)
!      call dapac(ina)

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      ibase = nomax + 1
!
      if(idif.gt.(nvmax+1)/2) then
         ider1  = 0
         ider1s = 0
         ider2  = idif-(nvmax+1)/2
         ider2s = 1
         do 10 jj=1,ider2-1
  10     ider2s = ider2s*ibase
         xdivi  = ider2s*ibase
      else
         ider1  = idif
         ider1s = 1
         do 20 jj=1,ider1-1
  20     ider1s = ider1s*ibase
         ider2  = 0
         ider2s = 0
         xdivi  = ider1s*ibase
      endif
!
      ibase = nomax+1
!
      ic = ipoc-1
!
      do 100 i=ipoa,ipoa+illa-1
!
      call dancd(i1(i),i2(i),jx)

      ikil1=0
      ikil2=0
      if(idif.gt.(nvmax+1)/2) then
       ikil2=jx(idif)
      else
       ikil1=jx(idif)
      endif
!
!      X = IEE/XDIVI
!etienne      IFAC = INT(IBASE*(X-INT(X)+1.D-8))
!
!etienne      IF(IFAC.EQ.0) GOTO 100
!
      ic = ic + 1
      cc(ic) = cc(i)
      i1(ic) = i1(i) - ikil1*ider1s
      i2(ic) = i2(i) - ikil2*ider2s
!
 100  continue
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
+if cr
         write(lout,*)'ERROR IN DATRASH '
+ei
+if .not.cr
         write(*,*)'ERROR IN DATRASH '
+ei
         call dadeb(111,'ERR DATRAS',1)
      endif
!
      return
      end
+dk dumps
+if debug
!DUMPS
      subroutine dumpda(dumpname,n,i)
      implicit none
      integer i,lda,lea,lia,lno,lnv,lst
+ca dabinc
      integer n
      character*(*) dumpname
      character*10 mydump
      mydump=dumpname
      write(99) mydump,n,i
      write(99) cc
      endfile 99
      backspace 99
      end
      subroutine wda(vname,value,i,j,k,l)
      implicit none
      integer i,lda,lea,lia,lno,lnv,lst
+ca dabinc
      integer n
      character*(*) vname
      double precision value
      integer j,k,l
      character*(16) myname,ccname
      myname=vname
      ccname='cc(50)'
      write(100) myname,value,i,j,k,l
      write(100) ccname,cc(50),50,0,0,0
      ccname='cc(64)'
      write(100) ccname,cc(64),64,0,0,0
      end
!DUMPS
+ei
