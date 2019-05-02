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
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : idao,iscrda,iscrri,rscrri,allvec,eps,epsmac,nda,ndamaxi,nst,nomax,nvmax,  &
        nmmax,nocut,lfi,idall,i1,i2,ie1,ie2,ieo,ia1,ia2,lda,lst,lea,lia,lno,lnv
      use crcoall
      use mod_units
      implicit none
      integer i,ibase,ic1,ic2,icmax,io1,io2,iout,iunit,j,jd,jj,jjj,jjjj,jl,js,k,n,nn,no,nv
      integer iall(1)
!     *****************************
!
!     THIS SUBROUTINE SETS UP THE MAJOR ORDERING AND ADDRESSING ARRAYS IN
!     COMMON BLOCK DAINI. IF IUNIT > 0, THE ARRAYS WILL BE PRINTED TO UNIT
!     NUMBER IUNIT. AN EXAMPLE FOR THE ARRAYS GENERATED BY DAINI CAN BE
!     FOUND AFTER THE ROUTINE.
!
!-----------------------------------------------------------------------------1

!-----------------------------------------------------------------------------9
!      COMMON / DASCR /  IS(20), RS(20)                                        1
!-----------------------------------------------------------------------------2


      character(len=10) aa
      dimension n(lnv+1),k(0:lnv),j(lnv),jj(lnv)
#ifdef BOINC
      character(len=256) filename
#endif

      if(eps.le.zero) eps=1.e-38_fPrec ! Why is this not pieni?
!      if(EPS.le.0.d0) eps=1.d-90
      epsmac=c1m7
      if(nv.eq.0) return
      ndamaxi=0

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

      do 10 i=0,lia
      ia1(i) = 0
      ia2(i) = 0
  10  continue

! Eric
      do i=1,lst
        i1(i)=0
        i2(i)=0
      enddo

      idao=0
      do i=1,100
        iscrda(i)=0
        iscrri(i)=0
        rscrri(i)=zero
      enddo

      if(nv.gt.lnv.or.no.gt.lno) then
         write(lout,*)'ERROR IN SUBROUTINE DAINI, NO, NV = ',no,nv
         call dadeb(31,'ERR DAINI ',1)
      endif

      ibase = no+1
      js    = nv/2
      if(real(ibase,fPrec)**((nv+1)/2).gt.real(lia,fPrec)) then          !hr10
         write(lout,*) 'ERROR, NO = ',no,', NV = ',nv,' TOO LARGE FOR',' LIA = ',lia
         call dadeb(31,'ERR DAINI ',1)
      endif

      icmax = 0
      nn    = 0
      k(0)  = 0

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

         ia2(ic2) = nn

         do 80 io1=0,no-io2
!        ******************

         n(js+1) = io1
         jd      = 1

  70     jl      = jl + jd

         if(jl.eq.js) then
            goto 80
         elseif(jd.eq.1) then
            j(jl) = 0
         else
            j(jl) = j(jl) + 1
         endif

         k(jl)    = k(jl-1)*ibase + j(jl)
         n(jl+1)  = n(jl) - j(jl)

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

            ie1(nn) = ic1
            ie2(nn) = ic2
            i1 (nn) = ic1
            i2 (nn) = ic2
            if(ic2.eq.0) ia1(ic1) = nn
            ieo(nn) = io1 + io2

            goto 70
         endif

   80    continue

         jd = -1
         goto 50
      endif

  100 continue
  110 continue

      if(nn.gt.lea) then
         write(lout,*)'ERROR IN DAINI, NN = ',nn,' EXCEEDS LEA'
         call dadeb(31,'ERR DAINI ',1)
      endif

!     ALLOCATING SCRATCH VARIABLES
!     ****************************

      iall(1) = 0
      call daall(iall(1),1,'$$UNPACK$$',nomax,nvmax)

      do i=0,nomax
        aa = '$$MUL   $$'
        write(aa(6:10),'(I5)') i
        iall(1) = 0
!       CALL DAALL(IALL,1,AA,I,NVMAX)
        call daall(iall(1),1,aa,nomax,nvmax)
      end do

      idall(1) = nmmax

!     DOUBLE CHECKING ARRAYS IE1,IE2,IA1,IA2
!     **************************************

      do i=1,nmmax
        jjj = ia1(ie1(i)) + ia2(ie2(i))
        if(jjj.ne.i) then
          write(lout,*) 'ERROR IN DAINI IN ARRAYS IE1,IE2,IA1,IA2 AT I = ',i
          call dadeb(31,'ERR DAINI ',1)
        endif
      end do

      if(iunit.eq.0) return

      write(lout,*)'ARRAY SETUP DONE, BEGIN PRINTING'

      iout = 32
      call f_open(unit=iout,file="daini.dat",formatted=.true.,mode="rw",status="new")
      write(iout,'(/A/A/)') ' ARRAYS I1 THROUGH I20, IE1,IE2,IEO',' **********************************'
      do i=1,nmmax
        call dancd(ie1(i),ie2(i),jj)
        write(iout,'(1X,I5,2X,4(5I2,1X),3I6)') i,(jj(jjjj),jjjj=1,lnv),ie1(i),ie2(i),ieo(i)
      end do

      write(iout,'(/A/A/)') ' ARRAYS IA1,IA2',' **************'

      do i=0,icmax
        write(iout,'(3I10)') i,ia1(i),ia2(i)
      end do

      return
end subroutine daini

subroutine daexter
      use floatPrecision
      use mod_lie_dab, only : allvec,lda
      implicit none
      integer i
!     *****************************
!
!-----------------------------------------------------------------------------1

!-----------------------------------------------------------------------------9

      do i=1, lda
        allvec(i)=.false.
      end do

      return
end subroutine daexter

subroutine dallsta(ldanow)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : allvec,lda
      implicit none
      integer i,ldanow
!     *****************************
!
!-----------------------------------------------------------------------------1

!-----------------------------------------------------------------------------9

      ldanow=0
      do 5 i=1, lda
        if(allvec(i)) ldanow=ldanow+1
5     continue

      write(lout,*) ' ALLOCATED ',ldanow

      return
end subroutine dallsta

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
subroutine daallno(ic,l,ccc)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : allvec,daname,nda,ndamaxi,nst,nomax,nvmax,nmmax,idano,idanv,idapo,    &
        idalm,idall,lda,lst
      implicit none
      integer i,ind,l,ndanum,no,nv
      real(kind=fPrec) x
!     ********************************
!
!     THIS SUBROUTINE ALLOCATES STORAGE FOR A DA VECTOR WITH
!     ORDER NOmax AND NUMBER OF VARIABLES NVmax
!
!-----------------------------------------------------------------------------1


      integer ic(*)
      logical incnda
      character(len=10) c,ccc

      no=nomax
      nv=nvmax
      ind = 1
      do 10 i=1,l
        if(ic(i).gt.0.and.ic(i).le.nda) then
!         DANAME(IC(I)) = C
!         IF(IDANO(IC(I)).EQ.NO.AND.IDANV(IC(I)).EQ.NV) THEN
        else
          if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
             write(lout,*)'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ', no,nv,' NOMAX, NVMAX = ',nomax,nvmax
            call dadeb(31,'ERR DAALL ',1)
          endif

20        if (allvec(ind)) then
           ind = ind + 1
           goto 20
          endif

          incnda = .false.
          if (ind .gt. nda) then
             incnda = .true.
             nda = nda + 1
             if(nda.gt.lda) then
         write(lout,*) 'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
             call dadeb(31,'ERR DAALL ',1)
             endif
          endif

          allvec(ind) = .true.

          ic(i) = ind

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

          if(nst.gt.lst) then
            x=-one
            write(lout,*)'ERROR IN DAALL, STACK EXHAUSTED '
            write(lout,*) ' NST,LST '
            write(lout,*)  nst,lst
            write(lout,*) ' NDA,NDANUM,NDA*NDANUM '
            write(lout,*)  nda,ndanum,nda*ndanum
!            X=DSQRT(X)
            call dadeb(31,'ERR DAALL ',1)
          endif

          if(nv.eq.0.or.nomax.eq.1) then
            call daclr(ic(i))
            idall(ic(i)) = idalm(ic(i))
          endif
        endif
  10  continue

      if(nda.gt.ndamaxi) ndamaxi=nda

      return
end subroutine daallno

subroutine daall(ic,l,ccc,no,nv)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : ndat,allvec,daname,nda,ndamaxi,nst,nomax,nvmax,nmmax,idano,idanv,     &
        idapo,idalm,idall,lda,lst
      implicit none
      integer i,ind,l,ndanum,no,nv
      real(kind=fPrec) x
!     ********************************
!
!     THIS SUBROUTINE ALLOCATES STORAGE FOR A DA VECTOR WITH
!     ORDER NO AND NUMBER OF VARIABLES NV
!
!-----------------------------------------------------------------------------1


      integer ic(*)
      logical incnda
      character(len=10) c,ccc

      ind = 1

      do 10 i=1,l
        if(ic(i).gt.0.and.ic(i).le.nda) then
!         DANAME(IC(I)) = C
!         IF(IDANO(IC(I)).EQ.NO.AND.IDANV(IC(I)).EQ.NV) THEN
        else
          if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
             write(lout,*)'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ', no,nv,' NOMAX, NVMAX = ',nomax,nvmax
            call dadeb(31,'ERR DAALL ',1)
          endif

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
        write(lout,*) 'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
             call dadeb(31,'ERR DAALL ',1)
             endif
          endif

          allvec(ind) = .true.

          ic(i) = ind

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

          if(nst.gt.lst) then
            x=-one
            write(lout,*)'ERROR IN DAALL, STACK EXHAUSTED '
            write(lout,*) ' NST,LST '
            write(lout,*)  nst,lst
            write(lout,*) ' NDA,NDANUM,NDA*NDANUM '
            write(lout,*)  nda,ndanum,nda*ndanum
!            X=DSQRT(X)
            call dadeb(31,'ERR DAALL ',1)
          endif

!          IF(NV.EQ.0) THEN
          if(nv.eq.0.or.nomax.eq.1) then
            call daclr(ic(i))
            idall(ic(i)) = idalm(ic(i))
          endif
        endif
  10  continue

      if(nda.gt.ndamaxi) ndamaxi=nda

      return
end subroutine daall

subroutine dadal(idal,l)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : allvec,daname,nda,nomax,nst,idapo,idall
      implicit none
      integer i,l
!     ************************
!
!     THIS SUBROUTINE DEALLOCATES THE VECTORS IDAL
!


      integer idal(*)

      do 10 i=l,1,-1
        if(idal(i).le.nomax+2.or.idal(i).gt.nda) then
          write(lout,*) 'ERROR IN ROUTINE DADAL, IDAL(I),NDA = ',idal(i),nda
          call dadeb(31,'ERR DADAL ',1)
        endif
        if(idal(i).eq.nda) then
          nst = idapo(nda) - 1
          nda = nda - 1
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
end subroutine dadal

subroutine davar(ina,ckon,i)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : eps,cc,nomax,nvmax,idall,i1,i2
      implicit none
      integer i,ibase,ic1,ic2,illa,ilma,ina,inoa,inva,ipoa
      real(kind=fPrec) ckon
!     ****************************
!
!     THIS SUBROUTINE DECLARES THE DA VECTOR
!     AS THE INDEPENDENT VARIABLE NUMBER I.
!
!-----------------------------------------------------------------------------1


      call dainf(ina,inoa,inva,ipoa,ilma,illa)


      if(i.gt.inva) then
         write(lout,*)'ERROR IN DAVAR, I = ',i,' EXCEEDS INVA = ',inva
         call dadeb(31,'ERR DAVAR ',1)
      endif

      if(nomax.eq.1) then
         if(i.gt.inva) then
            write(lout,*) 'ERROR IN DAVAR, I = ',i,' EXCEEDS INVA = ',inva
!           CALL DADEB(31,'ERR DAVAR3',1)
         endif
         call daclr(ina)
         cc(ipoa) = ckon
         cc(ipoa+i) = one
         return
      endif

      ibase = nomax+1

      if(i.gt.(nvmax+1)/2) then
        ic1 = 0
        ic2 = ibase**((i-(nvmax+1)/2)-1)                                 !hr10
      else
        ic1 = ibase**(i-1)
        ic2 = 0
      endif

      if(abs(ckon).gt.eps) then
         idall(ina) = 2
         cc(ipoa) = ckon
         i1(ipoa) = 0
         i2(ipoa) = 0

         cc(ipoa+1) = one
         i1(ipoa+1) = ic1
         i2(ipoa+1) = ic2
      else
         idall(ina) = 1
         cc(ipoa) = one
         i1(ipoa) = ic1
         i2(ipoa) = ic2
      endif

      return
end subroutine davar

subroutine dacon(ina,ckon)
      use floatPrecision
      use mod_lie_dab, only : eps,cc,nomax,idall,i1,i2
      implicit none
      integer illa,ilma,ina,inoa,inva,ipoa
      real(kind=fPrec) ckon
!     **************************
!
!     THIS SUBROUTINE SETS THE VECTOR C TO THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1


      call dainf(ina,inoa,inva,ipoa,ilma,illa)


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

      return
end subroutine dacon

subroutine danot(not)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : nomax,nocut
      implicit none
      integer not
!     *********************
!
!     THIS SUBROUTINE RESETS THE TRUNCATION ORDER NOCUT TO A NEW VALUE
!
!-----------------------------------------------------------------------------1


      if(not.gt.nomax) then
         write(lout,*)'ERROR, NOCUT = ',nocut,' EXCEEDS NOMAX = ',nomax
         call dadeb(31,'ERR DANOT ',1)
      endif

      nocut = not

      return
end subroutine danot

subroutine getdanot(not)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : nomax,nocut
      implicit none
      integer not
!     *********************
!
!     THIS SUBROUTINE RESETS THE TRUNCATION ORDER NOCUT TO A NEW VALUE
!
!-----------------------------------------------------------------------------1


      if(not.gt.nomax) then
         write(lout,*)'ERROR, NOCUT = ',nocut,' EXCEEDS NOMAX = ',nomax
         call dadeb(31,'ERR DANOT ',1)
      endif

      not=nocut

      return
end subroutine getdanot

subroutine daeps(deps)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : eps
      implicit none
      real(kind=fPrec) deps
!     **********************
!
!     THIS SUBROUTINE RESETS THE TRUNCATION ORDER NOCUT TO A NEW VALUE
!
!-----------------------------------------------------------------------------1


      if(deps.ge.zero) then
        eps = deps
      else
        deps=eps
      endif

      return
end subroutine daeps

subroutine dapek(ina,jj,cjj)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : cc,nomax,nvmax,i1,i2,ia1,ia2,lia,lnv
      implicit none
      integer i,ibase,ic,ic1,ic2,icu,icz,ii1,ikk,illa,ilma,ina,         &
        inoa,inva,ipek,ipoa,iu,iz,jj,jj1,mchk
      real(kind=fPrec) cjj
!     ****************************
!
!     THIS SUBROUTINE DETERMINES THE COEFFICIENT OF THE ARRAY
!     OF EXPONENTS JJ AND RETURNS IT IN CJJ
!
!-----------------------------------------------------------------------------1


      dimension jj(lnv)

      call dainf(ina,inoa,inva,ipoa,ilma,illa)

      if(illa.eq.0) then   ! etienne shit
        cjj = zero                                                        !hr10
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
            write(lout,*)                                               &
     &'ERROR IN DAPEK, INDEX OUTSIDE RANGE, JJ(1) = ',jj1
!           CALL DADEB(31,'ERR DAPEK1',1)
         endif
         ipek = ipoa + jj1 - 1
         cjj = cc(ipek)
         return
      endif

      ii1 = (nvmax+1)/2
      ibase = nomax+1

!     DETERMINE INDEX TO BE SEARCHED FOR
!     **********************************

      call dadcd(jj,ic1,ic2)

!ETIENNE
      if(ic1.gt.lia.or.ic2.gt.lia) then
       write(lout,*) 'DISASTER IN DAPEK, INA= ',ina
       write(lout,*) ic1,ic2
       write(lout,*) (jj(ikk),ikk=1,lnv)
      endif
!ETIENNE
      ic = ia1(ic1) + ia2(ic2)

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

      iu = ipoa
      iz = ipoa + illa - 1
      icu = ia1(i1(iu))+ia2(i2(iu))
      icz = ia1(i1(iz))+ia2(i2(iz))

      if(illa.eq.0) then
         cjj = zero                                                       !hr10
         return
      elseif(ic.eq.icu) then
         cjj = cc(iu)
         return
      elseif(ic.eq.icz) then
         cjj = cc(iz)
         return
      elseif(ic.lt.icu.or.ic.gt.icz) then
         cjj = zero                                                       !hr10
         return
      endif

!     SEARCHING PROPER MONOMIAL
!     *************************

 10   continue
      if(iz-iu.le.1) then
         cjj = 0
         return
      endif
      i = (iu+iz)/2

!     if(ia1(i1(i))+ia2(i2(i)) - ic) 20,30,40
      mchk=(ia1(i1(i))+ia2(i2(i))) - ic                                  !hr10
      if(mchk.lt.0) goto 20
      if(mchk.eq.0) goto 30
      if(mchk.gt.0) goto 40
 20   iu = i
      goto 10
 30   cjj = cc(i)

      return
 40   iz = i
      goto 10

end subroutine dapek

subroutine dapok(ina,jj,cjj)
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : cc,eps,nomax,idalm,idall,i1,i2,ia1,ia2,lnv
      implicit none
      integer i,ic,ic1,ic2,icu,icz,ii,illa,ilma,ina,inoa,inva,ipoa,ipok,iu,iz,jj,jj1,mchk
      real(kind=fPrec) cjj
!     ****************************
!
!     THIS SUBROUTINE SETS THE COEFFICIENT OF THE ARRAY
!     OF EXPONENTS JJ TO THE VALUE CJJ
!
!-----------------------------------------------------------------------------1


      dimension jj(lnv)

      call dainf(ina,inoa,inva,ipoa,ilma,illa)

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
            write(lout,*) 'ERROR IN DAPOK, INDEX OUTSIDE RANGE, JJ(1) = ',jj1
!           CALL DADEB(31,'ERR DAPOK1',1)
         endif
         ipok = ipoa + jj1 - 1
         cc(ipok) = cjj
         return
      endif

!     DETERMINE INDEX TO BE SEARCHED FOR
!     **********************************

      call dadcd(jj,ic1,ic2)

      ic = ia1(ic1) + ia2(ic2)
!
!     DETERMINE IF MONOMIAL TO BE POKED CONFORMS WITH INOA, INVA,NOCUT
!     ****************************************************************
!
      if(illa.ne.0) then ! etienne shit
      iu = ipoa
      iz = ipoa + illa - 1

!     DETERMINE IF MONOMIAL IS INSIDE FIRST AND LAST MONOMIALS OF A
!     *************************************************************

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


!     SEARCHING PLACE TO POKE INTO OR BEFORE WHICH TO POKE
!     ****************************************************

      iu = ipoa
      iz = ipoa + illa

 10   continue
      if(iz-iu.le.1) then
         i = iz
         goto 100
      endif
      i = (iu+iz)/2

!      if(ia1(i1(i))+ia2(i2(i)) - ic) 20,30,40
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

!     INSERTING THE MONOMIAL, MOVING THE REST
!     ***************************************

 100  continue

      if(abs(cjj).lt.eps) return

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
         write(lout,*)'ERROR IN DAPOK '
         call dadeb(31,'ERR DAPOK ',1)
      endif
!
      return
!
!     CASE OF CJJ = 0 WHICH MEANS MOVING THE REST
!     *********************************************
!
 200  continue
      if(abs(cjj).lt.eps) then
         do ii=i,ipoa+illa-2
           cc(ii) = cc(ii+1)
           i1(ii) = i1(ii+1)
           i2(ii) = i2(ii+1)
         end do
         idall(ina) = illa - 1
      endif
      return

end subroutine dapok

subroutine daclr(inc)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : cc
      implicit none
      integer i,illc,ilmc,inc,inoc,invc,ipoc
!     *********************
!
!     THIS SUBROUTINE SETS ALL THE STACK SPACE RESERVED FOR VARIABLE
!     C TO ZERO
!
!-----------------------------------------------------------------------------1


      call dainf(inc,inoc,invc,ipoc,ilmc,illc)

      do i=ipoc,ipoc+ilmc-1
        cc(i) = zero
      end do

      return
end subroutine daclr

subroutine dacop(ina,inb)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only :  cc,nomax,nocut,idalm,idall,i1,i2,ieo,ia1,ia2
      implicit none
      integer ia,ib,iif,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,ipoa,ipob
!     *************************
!
!     THIS SUBROUTINE COPIES THE DA VECTOR A TO THE DA VECTOR B
!
!-----------------------------------------------------------------------------1


      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)

!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)

      ib = ipob - 1

      iif = 0
      if(nomax.eq.1.or.inva.eq.0) iif = 1

      do 100 ia = ipoa,ipoa+illa-1

      if(iif.eq.0) then
        if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 100
      endif
      ib = ib + 1
      cc(ib) = cc(ia)
      i1(ib) = i1(ia)
      i2(ib) = i2(ia)

 100  continue

      idall(inb) = (ib - ipob) + 1                                       !hr10
      if(idall(inb).gt.idalm(inb)) then
         write(lout,*)'ERROR IN DACOP'
         call dadeb(31,'ERR DACOP ',1)
      endif

      return
end subroutine dacop

subroutine datrashn(idif,ina,inbb)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : nomax,i1,i2,lnv
      implicit none
      integer i,ia,idif,illa,ilma,ina,inbb,inoa,inva,ipoa
      integer inb(1)
      real(kind=fPrec) rr
!     *************************
!
!     THIS SUBROUTINE COPIES THE DA VECTOR A TO THE DA VECTOR B
!
!-----------------------------------------------------------------------------1

      integer jd(lnv)

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      inb(1)=0

      if(inbb.eq.ina) then
         call daall(inb(1),1,'$$DAADD $$',inoa,inva)
        else
         inb(1)=inbb
      endif

      call daclr(inb(1))

      do ia = ipoa,ipoa+illa-1

        if(nomax.ne.1) then
          call dancd(i1(ia),i2(ia),jd)
        else
          do i=1,lnv
            jd(i)=0
          end do

          if(ia.ne.ipoa) then
            jd(ia-ipoa+1)=1
          endif
        endif

        call dapek(ina,jd,rr)
        jd(idif)=0
        if(abs(rr).gt.zero) call dapok(inb(1),jd,rr)
      end do


      if(inbb.eq.ina) then
         call dacop(inb(1),inbb)
         call dadal(inb(1),1)
      endif

      return
end subroutine datrashn

subroutine daadd(ina,inb,inc)
      use floatPrecision
      use numerical_constants
      implicit none
      integer idaadd(1),illc,ilmc,ina,inb,inc,inoc,invc,ipoc
!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA ADDITION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C.

      if(ina.ne.inc.and.inb.ne.inc) then
         call dalin(ina,+one,inb,+one,inc)
      else
         idaadd(1) = 0
         call dainf(inc,inoc,invc,ipoc,ilmc,illc)
         call daall(idaadd(1),1,'$$DAADD $$',inoc,invc)
         call dalin(ina,+one,inb,+one,idaadd(1))
         call dacop(idaadd(1),inc)
         call dadal(idaadd(1),1)
      endif

      return
end subroutine daadd

subroutine dasub(ina,inb,inc)
      use floatPrecision
      use numerical_constants
      implicit none
      integer idasub(1),illc,ilmc,ina,inb,inc,inoc,invc,ipoc
!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA SUBTRACTION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C.

      if(ina.ne.inc.and.inb.ne.inc) then
         call dalin(ina,+one,inb,-one,inc)
      else
         idasub(1) = -1
         call dainf(inc,inoc,invc,ipoc,ilmc,illc)
         call daall(idasub(1),1,'$$DASUB $$',inoc,invc)
         call dalin(ina,+one,inb,-one,idasub(1))
         call dacop(idasub(1),inc)
         call dadal(idasub(1),1)
      endif

      return
end subroutine dasub

subroutine damulin(ina,inb,coe1,inc,ind,coe2,ine)
      use floatPrecision
      implicit none
      integer ina,inb,inc,incc(1),ind,ine,inoc,invc
      real(kind=fPrec) coe1,coe2
!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA MULTIPLICATION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C. AS TEMPORARY STORAGE, THE STACK SPACE
!     OF THE (NOMAX+2) SCRATCH VARIABLES ALLOCATED BY DAINI IS USED.
!
!-----------------------------------------------------------------------------1


      call daall(incc(1),1,'$$DAJUNK$$',inoc,invc)
      call damul(ina,inb,incc(1))
      call damul(inc,ind,ine)
      call dalin(incc(1),coe1,ine,coe2,ine )
      call dadal(incc(1),1)

      return
end subroutine damulin

! ANFANG UNTERPROGRAMM
subroutine daexx(ina,inb,inc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      implicit none
      integer illc,ilmc,ina,inb,inc,inoc,invc,ipoc
      integer inaa(1),inbb(1)

!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INE WITH THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
       write(lout,*) "daexx"
      if(ina.ne.inc.and.inb.ne.inc) then
         call daexxt(ina,inb,inc)
      else
         inaa(1) = 0
         inbb(1) = 0
         call dainf(inc,inoc,invc,ipoc,ilmc,illc)
         call daall(inaa(1),1,'$$DAADD $$',inoc,invc)
         call daall(inbb(1),1,'$$DAADD $$',inoc,invc)
         call dacop(ina,inaa(1))
         call dacop(inb,inbb(1))
         call daexxt(inaa(1),inbb(1),inc)
         call dadal(inaa(1),1)
         call dadal(inbb(1),1)
      endif

      return
end subroutine daexx

! ANFANG UNTERPROGRAMM

subroutine daexxt(ina,inb,inc)
      use floatPrecision
      implicit none
      integer idaexx(1),illa,illb,illc,ilma,ilmb,ilmc,ina,inb,inc,inoa, &
     &inob,inoc,inva,invb,invc,ipoa,ipob,ipoc
!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INA WITH INB
!
!-----------------------------------------------------------------------------1

!-----------------------------------------------------------------------------9
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)

      idaexx(1) = 0
      call daall(idaexx(1),1,'$$DAEXX $$',inoc,invc)
      call dafun('LOG   ',ina,inc)
      call damul(inb,inc,idaexx(1))
      call dafun('EXP   ',idaexx(1),inc)
      call dadal(idaexx(1),1)

      return
end subroutine daexxt

! ANFANG UNTERPROGRAMM

subroutine dacex(ina,ckon,inb)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      implicit none
      integer illc,ilmc,ina,inb,inc,inoc,invc,ipoc
      integer incc(1)
      real(kind=fPrec) ckon

!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INE WITH THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
        write(lout,*) "dacex"
      if(ina.eq.inb) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc(1)=0
        call daall(incc(1),1,'$$DAJUNK$$',inoc,invc)
        call dacext(ina,ckon,incc(1))
        call dacop(incc(1),inb)
        call dadal(incc(1),1)
      else
        call dacext(ina,ckon,inb)
      endif

      return
end subroutine dacex


subroutine dacext(ina,ckon,inb)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      implicit none
      integer idacex(1),illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,&
     &ipoa,ipob
      real(kind=fPrec) ckon
!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES THE CONSTANT CKON WITH INA
!
!-----------------------------------------------------------------------------1

!-----------------------------------------------------------------------------9
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)

      if(ckon.le.0) then
         write(lout,*)'ERROR IN DACEX, CKON NOT POSITIVE'
!        CALL DADEB(31,'ERR DACEX1',1)
      endif

      idacex(1) = 0
      call daall(idacex(1),1,'$$DACEX $$',inob,invb)
      ckon = log_mb(ckon)
      call dacmu(ina,ckon,idacex(1))
      call dafun('EXP   ',idacex(1),inb)
      call dadal(idacex(1),1)

      return
end subroutine dacext


subroutine daexc(ina,ckon,inb)
      use floatPrecision
      implicit none
      integer illc,ilmc,ina,inb,inc,inoc,invc,ipoc
      integer incc(1)
      real(kind=fPrec) ckon

!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INE WITH THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1

      if(ina.eq.inb) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc(1)=0
        call daall(incc(1),1,'$$DAJUNK$$',inoc,invc)
        call daexct(ina,ckon,incc(1))
        call dacop(incc(1),inb)
        call dadal(incc(1),1)
      else
        call daexct(ina,ckon,inb)
      endif

      return
end subroutine daexc


subroutine daexct(ina,ckon,inb)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : eps
      implicit none
      integer i,ic,idaexc(1),illa,illb,ilma,ilmb,ina,inb,inoa,inob,     &
     &inva,invb,ipoa,ipob
      real(kind=fPrec) ckon,xic

!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INE WITH THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)

      idaexc(1) = 0
      call daall(idaexc(1),1,'$$DAEXC $$',inob,invb)

      if(ckon.lt.zero) then
        call dafun('LOG   ',ina,inb)
        call dacmu(inb,ckon,idaexc(1))
        call dafun('EXP   ',idaexc(1),inb)
      else
        xic=abs(ckon-real(int(ckon),fPrec))                               !hr10
        if(xic.gt.eps) then
          call dafun('LOG   ',ina,inb)
          call dacmu(inb,ckon,idaexc(1))
          call dafun('EXP   ',idaexc(1),inb)
        else
          ic=int(ckon)
          call dacon(idaexc(1),one)
          do i=1,ic
            call damul(idaexc(1),ina,idaexc(1))
          enddo
          call dacop(idaexc(1),inb)
        endif
      endif
      call dadal(idaexc(1),1)

      return
end subroutine daexct


subroutine damul(ina,inb,inc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      implicit none
      integer illc,ilmc,ina,inb,inc,inoc,invc,ipoc,incc(1)
!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA MULTIPLICATION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C. AS TEMPORARY STORAGE, THE STACK SPACE
!     OF THE (NOMAX+2) SCRATCH VARIABLES ALLOCATED BY DAINI IS USED.
!
!-----------------------------------------------------------------------------1

!
!

      if(ina.eq.inc.or.inb.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc(1)=0
        call daall(incc(1),1,'$$DAJUNK$$',inoc,invc)
        call damult(ina,inb,incc(1))
        call dacop(incc(1),inc)
        call dadal(incc(1),1)
      else
        call damult(ina,inb,inc)
      endif

      return
end subroutine damul


subroutine damult(ina,inb,inc)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : cc,nomax,nocut,idapo,idall,i1,i2,ieo,ia1,ia2,lno
      implicit none
      integer i,i1ia,i2ia,ia,ib,ic,illa,illb,illc,ilma,ilmb,ilmc,ina,   &
     &inb,inc,inoa,inob,inoc,inva,invb,invc,ioffb,ipno,ipoa,ipob,ipoc,  &
     &ipos,minv,noff,noib,nom
      real(kind=fPrec) ccia,ccipoa,ccipob
      logical checkMultUnderflow
!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA MULTIPLICATION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C. AS TEMPORARY STORAGE, THE STACK SPACE
!     OF THE (NOMAX+2) SCRATCH VARIABLES ALLOCATED BY DAINI IS USED.
!
!-----------------------------------------------------------------------------1

!
      dimension ipno(0:lno),noff(0:lno)

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)

!      CALL DACHK(INA,INOA,INVA, INB,INOB,INVB, INC,INOC,INVC)

!     CASE OF FIRST ORDER ONLY
!     ************************

      if(nomax.eq.1) then
        minv = min(inva,invb,invc)
        ccipoa = cc(ipoa)
        ccipob = cc(ipob)
        cc(ipoc) = ccipoa*ccipob

        do i=1,minv
          cc(ipoc+i) = ccipoa*cc(ipob+i) + ccipob*cc(ipoa+i)
        end do

         do i=ipoc+minv+1,ipoc+invc
           cc(i) = zero
         end do

         return
      endif

!     GENERAL CASE
!     ************

      do i=0,nomax
        noff(i) = idapo(i+2)
        ipno(i) = 0
      end do

      call daclr(1)

!     RE-SORTING THE VECTOR B INTO PIECES THAT ARE OF ONLY ONE ORDER
!     *************************************************************

      do ib=ipob,ipob+illb-1
        noib = ieo(ia1(i1(ib))+ia2(i2(ib)))
        ipos = ipno(noib) + 1
        ipno(noib) = ipos
        inob = noff(noib) + ipos

        cc(inob) = cc(ib)
        i1(inob) = i1(ib)
        i2(inob) = i2(ib)
      end do

      do i=0,nomax
        idall(i+2) = ipno(i)
      end do

!     PERFORMING ACTUAL MULTIPLICATION
!     ********************************

      nom = min(nocut,inoc)

      do ia=ipoa,ipoa+illa-1
        i1ia = i1(ia)
        i2ia = i2(ia)
        ccia = cc(ia)

        do noib = 0,nom-ieo(ia1(i1(ia))+ia2(i2(ia)))
          ioffb = noff(noib)

          do ib = ioffb+1,ioffb+ipno(noib)
            ic = ia2(i2ia+i2(ib)) + ia1(i1ia + i1(ib))
            cc(ic) = cc(ic) + ccia*cc(ib)
          end do
        end do
      end do

      call dapac(inc)

      return
end subroutine damult


subroutine dadiv(ina,inb,inc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      implicit none
      integer idadiv(1),illc,ilmc,ina,inb,inc,inoc,invc,ipoc
!     *****************************

!     THIS SUBROUTINE PERFORMS A DA DIVISION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C.

      idadiv(1) = 0
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
      call daall(idadiv(1),1,'$$DADIV $$',inoc,invc)
      call dafun('INV ',inb,idadiv(1))
      call damul(ina,idadiv(1),inc)
      call dadal(idadiv(1),1)

      return
end subroutine dadiv


subroutine dasqr(ina,inc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      implicit none
      integer illc,ilmc,ina,inc,inoc,invc,ipoc,incc(1)
!     *************************
!
!     THIS SUBROUTINE SQUARES THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1


      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc(1)=0
        call daall(incc(1),1,'$$DAJUNK$$',inoc,invc)
        call dasqrt(ina,incc(1))
        call dacop(incc(1),inc)
        call dadal(incc(1),1)
      else
        call dasqrt(ina,inc)
      endif

      return
end subroutine dasqr


subroutine dasqrt(ina,inc)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : cc,nomax,nocut,idapo,idalm,idall,i1,i2,ieo,ia1,ia2,lno
      implicit none
      integer i,i1ia,i2ia,ia,ib,ib1,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,ioffa,ioffb, &
        ipno,ipoa,ipoc,ipos,minv,noff,noia,noib,nom
      real(kind=fPrec) ccia,ccipoa
!     *************************
!
!     THIS SUBROUTINE SQUARES THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1


      dimension ipno(0:lno),noff(0:lno)

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)

!      CALL DACHK(INA,INOA,INVA,'          ',-1,-1,INC,INOC,INVC)

      if(inva+invc.eq.0) then
         do i=0,illa-1
           cc(ipoc+i) = cc(ipoa+i)**2                                      !hr10
         end do

         idall(inc) = idall(ina)
         if(idall(inc).gt.idalm(inc)) then
            write(lout,*)'ERROR IN DASQR '
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
         cc(ipoc) = ccipoa**2                                            !hr10

         do i=1,minv
           cc(ipoc+i) = two*ccipoa*cc(ipoa+i)
         end do

         do i=ipoc+minv+1,ipoc+invc
           cc(i) = zero
         end do

         return
      endif

!     GENERAL CASE
!     ************

      do i=0,nomax
        noff(i) = idapo(i+2)
        ipno(i) = 0
      end do

      call daclr(1)

!     RESORTING THE VECTOR A INTO PIECES THAT ARE OF ONLY ONE ORDER
!     *************************************************************

      do ia=ipoa,ipoa+illa-1
        noia = ieo(ia1(i1(ia))+ia2(i2(ia)))
        ipos = ipno(noia) + 1
        ipno(noia) = ipos
        inoa = noff(noia) + ipos

        cc(inoa) = cc(ia)
        i1(inoa) = i1(ia)
        i2(inoa) = i2(ia)
      end do

      do i=0,nomax
        idall(i+2) = ipno(i)
      end do

!     PERFORMING ACTUAL MULTIPLICATION
!     ********************************

      nom = min(nocut,inoc)

      do noia = 0,nom/2
        ioffa = noff(noia)

        do ia=ioffa+1,ioffa+ipno(noia)
          i1ia = i1(ia)
          i2ia = i2(ia)
          ccia = cc(ia)

          ic = ia2(i2ia+i2ia) + ia1(i1ia+i1ia)
          cc(ic) = cc(ic) + ccia**2                                          !hr10
          ccia = ccia + ccia

          do noib = noia,nom-noia

            ioffb = noff(noib)
            if(noib.eq.noia) then
              ib1 = ia + 1
            else
              ib1 = ioffb + 1
            endif

            do ib = ib1,ioffb+ipno(noib)
              ic = ia2(i2ia+i2(ib)) + ia1(i1ia + i1(ib))
              cc(ic) = cc(ic) + ccia*cc(ib)
            end do
          end do
        end do
      end do

      call dapac(inc)

      return
end subroutine dasqrt


subroutine dacad(ina,ckon,inb)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : cc,nomax,lnv
      implicit none
      integer illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,ipoa,ipob
      real(kind=fPrec) ckon,const
!     ******************************
!
!     THIS SUBROUTINE ADDS THE CONSTANT CKON TO THE VECTOR A
!
!-----------------------------------------------------------------------------1

      integer jj(lnv)
      data jj / lnv*0 /

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)

      call dacop(ina,inb)

      if(nomax.eq.1) then
         cc(ipob) = cc(ipob) + ckon
         return
      endif

      call dapek(inb,jj,const)
      call dapok(inb,jj,const+ckon)

      return
end subroutine dacad


subroutine dacsu(ina,ckon,inb)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : cc,nomax,lnv
      implicit none
      integer illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,ipoa,ipob
      real(kind=fPrec) ckon,const
!     ******************************
!
!     THIS SUBROUTINE SUBTRACTS THE CONSTANT CKON FROM THE VECTOR A
!
!-----------------------------------------------------------------------------1

      integer jj(lnv)
      data jj / lnv*0 /

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
      call dacop(ina,inb)


      if(nomax.eq.1) then
         cc(ipob) = cc(ipob) - ckon
         return
      endif
!
      call dapek(inb,jj,const)
      call dapok(inb,jj,const-ckon)

      return
      end


subroutine dasuc(ina,ckon,inb)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      implicit none
      integer illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,ipoa,ipob
      real(kind=fPrec) ckon
!     ******************************
!
!     THIS SUBROUTINE SUBTRACTS THE VECTOR INA FROM THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1


      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)

      call dacsu(ina,ckon,inb)
      call dacmu(inb,-one,inb)

      return
      end


subroutine dacmu(ina,ckon,inc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      implicit none
      integer illc,ilmc,ina,inc,inoc,invc,ipoc,incc(1)
      real(kind=fPrec) ckon
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1


      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc(1)=0
        call daall(incc(1),1,'$$DAJUNK$$',inoc,invc)
        call dacmut(ina,ckon,incc(1))
        call dacop(incc(1),inc)
        call dadal(incc(1),1)
      else
        call dacmut(ina,ckon,inc)
      endif
      return
      end


subroutine dacmut(ina,ckon,inb)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : cc,eps,nomax,nocut,idalm,idall,i1,i2,ieo,ia1,ia2
      implicit none
      integer i,ia,ib,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,ipoa,ipob,minv
      real(kind=fPrec) ckon
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)

!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)

      if(nomax.eq.1) then
         minv = min(inva,invb)

         do i=0,minv
           cc(ipob+i) = cc(ipoa+i) * ckon
         end do

         do i=ipob+minv+1,ipob+invb
           cc(i) = zero
         end do

         return
      endif

      if(abs(ckon).lt.eps) then
         idall(inb) = 0
         return
      endif

      ib = ipob - 1

      do 100 ia=ipoa,ipoa+illa-1

      if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 100
      ib = ib + 1
      cc(ib) = cc(ia)*ckon
      i1(ib) = i1(ia)
      i2(ib) = i2(ia)

 100  continue

      idall(inb) = (ib-ipob)+1                                           !hr10
      if(idall(inb).gt.idalm(inb)) then
        write(lout,*)'ERROR IN DACMU '
        call dadeb(31,'ERR DACMU ',1)
      endif

      return
      end


subroutine dacdi(ina,ckon,inb)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      implicit none
      integer illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,ipoa,ipob
      real(kind=fPrec) ckon
!     ******************************
!
!     THIS SUBROUTINE DIVIDES THE VECTOR INA BY THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1


      if(ckon.eq.zero) then
         write(lout,*)'ERROR IN DACDI, CKON IS ZERO'
         call dadeb(31,'ERR DACDI ',1)
      endif

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)

      call dacmu(ina,one/ckon,inb)

      return
      end


subroutine dadic(ina,ckon,inc)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : eps
      implicit none
      integer idadic(1),illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,ipoa,ipoc
      real(kind=fPrec) ckon
!      parameter(zero=0.0_fPrec)
!     ******************************
!
!     THIS SUBROUTINE DIVIDES THE CONSTANT CKON BY THE VECTOR INA
!
!-----------------------------------------------------------------------------1


      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)

      if(abs(ckon).lt.eps) then
        call dacon(inc,zero)
        return
      endif

      idadic(1) = 0
      call daall(idadic(1),1,'$$DADIC $$',inoc,invc)

      if(ckon.eq.zero) then
        write(lout,*)'ERROR IN DACDI and DADIC, CKON IS ZERO'
        call dadeb(31,'ERR DACDI ',1)
      endif
      call dacdi(ina,ckon,idadic(1))
      call dafun('INV ',idadic(1),inc)
      call dadal(idadic(1),1)

      return
      end


subroutine dacma(ina,inb,bfac,inc)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      implicit none
      integer illc,ilmc,ina,inb,inc,inoc,invc,ipoc,idacma(1)
      real(kind=fPrec) bfac
!     **********************************
!
!     THIS SUBROUTINE PERFORMS THE OPERATIONS C = A + B*BFAC, WHERE A,B,C ARE
!     DA VECTORS AND BFAC IS A DOUBLE PRECISION. A AND C CAN BE IDENTICAL.
!     CAN LATER BE REPLACED BY SOMETHING LIKE DAADD WITH MINOR CHANGES.
!
!-----------------------------------------------------------------------------1


      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
      idacma(1) = 0
      call daall(idacma(1),1,'$$DACMA $$',inoc,invc)
      call dalin(ina,+one,inb,bfac,idacma(1))
      call dacop(idacma(1),inc)
      call dadal(idacma(1),1)

      return
      end


subroutine dalin(ina,afac,inb,bfac,inc)
      use floatPrecision
      implicit none
      integer illc,ilmc,ina,inb,inc,inoc,invc,ipoc,incc(1)
      real(kind=fPrec) afac,bfac
!     ***************************************
!
!     THIS SUBROUTINE COMPUTES THE LINEAR COMBINATION
!     C = AFAC*A + BFAC*B. IT IS ALSO USED TO ADD AND SUBTRACT.
!
!-----------------------------------------------------------------------------1


      if(ina.eq.inc.or.inb.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc(1)=0
        call daall(incc(1),1,'$$DAJUNK$$',inoc,invc)
        call dalint(ina,afac,inb,bfac,incc(1))
        call dacop(incc(1),inc)
        call dadal(incc(1),1)
      else
        call dalint(ina,afac,inb,bfac,inc)
      endif

      return
      end


subroutine dalint(ina,afac,inb,bfac,inc)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : cc,eps,nomax,nocut,idalm,idall,i1,i2,ieo,ia1,ia2
      implicit none
      integer i,ia,iamax,ib,ibmax,ic,icmax,illa,illb,illc,ilma,ilmb,ilmc,ina,inb,inc,inoa,inob,inoc,&
        inva,invb,invc,ipoa,ipob,ipoc,is,ismax,ismin,ja,jb,minv,mchk
      real(kind=fPrec) afac,bfac,ccc,copf
!     ***************************************
!
!     THIS SUBROUTINE COMPUTES THE LINEAR COMBINATION
!     C = AFAC*A + BFAC*B. IT IS ALSO USED TO ADD AND SUBTRACT.
!
!-----------------------------------------------------------------------------1


      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)

!      CALL DACHK(INA,INOA,INVA, INB,INOB,INVB, INC,INOC,INVC)

      if(nomax.eq.1) then
        minv = min(inva,invb,invc)
        do i=0,minv
          cc(ipoc+i) = cc(ipoa+i) * afac + cc(ipob+i) * bfac
        end do

        do i=ipoc+minv+1,ipoc+invc
          cc(i) = zero
        end do

        return
      endif

      ia = ipoa
      ib = ipob
      ic = ipoc - 1
      iamax = (ipoa+illa)-1                                              !hr10
      ibmax = (ipob+illb)-1                                              !hr10
      icmax = (ipoc+ilmc)-1                                              !hr10
      ja = ia1(i1(ia)) + ia2(i2(ia))
      jb = ia1(i1(ib)) + ia2(i2(ib))

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

!     COMPARING
!     *********

  10  continue
!      if(ja-jb) 30,20,40
      mchk=ja-jb
      if(mchk.lt.0) goto 30
      if(mchk.eq.0) goto 20
      if(mchk.gt.0) goto 40

!     ADDING TWO TERMS
!     ****************

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

!     STORING TERM A
!     **************

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

!     STORING TERM B
!     **************

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

!     COPYING THE REST
!     ****************

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

      idall(inc) = (ic - ipoc) + 1                                       !hr10

      if(idall(inc).gt.idalm(inc)) then
        write(lout,*)'ERROR IN DALIN, RESULT HAS TOO MANY TERMS '
        call dadeb(31,'ERR DALIN ',1)
      endif

      return
      end


subroutine dafun(cf,ina,inc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      implicit none
      integer illc,ilmc,ina,inc,inoc,invc,ipoc,incc(1)
!     ****************************
!
!     THIS SUBROUTINE COMPUTES THE FUNCTION CF OF THE DA VECTOR A
!     AND STORES THE RESULT IN C.
!     AT PRESENT, SOME FUNCTIONS CAN BE COMPUTED ONLY TO FIFTH ORDER.
!     THIS HAS TO BE FIXED IN THE FUTURE.
!
!-----------------------------------------------------------------------------1

      character(len=4) cf

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc(1)=0
        call daall(incc(1),1,'$$DAJUNK$$',inoc,invc)
        call dafunt(cf,ina,incc(1))
        call dacop(incc(1),inc)
        call dadal(incc(1),1)
      else
        call dafunt(cf,ina,inc)
      endif

      return
      end


subroutine dafunt(cf,ina,inc)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : epsmac,nocut,lno,lnv
      implicit none
      integer i,illa,illc,ilma,ilmc,ina,inc,ind,inoa,inoc,inon(1),inva, &
     &invc,ipoa,ipoc,ipow(1),iscr(1),jj,lfun,no

      real(kind=fPrec) a0,a1,a2,a3,a4,a5,ca,e1,e2,ea,era,p,ra,rpi4,sa,scr,t,xf
!     ****************************
!
!     THIS SUBROUTINE COMPUTES THE FUNCTION CF OF THE DA VECTOR A
!     AND STORES THE RESULT IN C.
!     AT PRESENT, SOME FUNCTIONS CAN BE COMPUTED ONLY TO FIFTH ORDER.
!     THIS HAS TO BE FIXED IN THE FUTURE.
!
!-----------------------------------------------------------------------------1


      character(len=4) cf,cfh
      character(len=26) abcs,abcc
      dimension xf(0:lno),jj(lnv)

      data jj /lnv*0/
      data abcs /'abcdefghijklmnopqrstuvwxyz'/
      data abcc /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      if(cf(1:1).eq.' ') then
        cfh(1:3) = cf(2:4)
        cfh(1:4) = ' '
        cf = cfh
      endif

      do i=1,4
        ind = index(abcs,cf(i:i))
        if(ind.ne.0) cf(i:i) = abcc(ind:ind)
      end do

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)

!     CASE OF NV = 0 WHICH MEANS COORDINATEWISE OPERATION
!     ***************************************************
!
!     CASE OF NV > 0 WHICH MEANS DIFFERENTIAL ALGEBRAIC OPERATION
!     ***********************************************************

      if(cf.eq.'SQR ') then
        call dasqr(ina,inc)
        return
      endif

!     ALLOCATE VARIABLES, PICK ZEROTH ORDER TERM
!     ******************************************

      ipow(1) = 0
      inon(1) = 0
      iscr(1) = 0

      call daall(ipow(1),1,'$$DAFUN1$$',inoc,invc)
      call daall(inon(1),1,'$$DAFUN2$$',inoc,invc)
      call daall(iscr(1),1,'$$DAFUN3$$',inoc,invc)

!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)

      call dapek(ina,jj,a0)

      no = min(nocut,inoa,inoc)

!     BRANCHING TO DIFFERENT FUNCTIONS
!     ********************************
!
      if(cf.eq.'INV ') then
!        1/(A0+P) = 1/A0*(1-(P/A0)+(P/A0)**2-...)
         if(a0.eq.zero) then                                              !hr10
            write(lout,1000) cf,ina,a0
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
         xf(0) = one/a0
         do i=1,no
           xf(i) = (-one*xf(i-1))/a0
         end do

      elseif(cf.eq.'SQRT') then
!        SQRT(A0+P) = SQRT(A0)*(1+1/2(P/A0)-1/8*(P/A0)**2+...)
         if(a0.le.zero) then                                              !hr10
            write(lout,1000) cf,ina,a0
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
         ra = sqrt(a0)                                                   !hr10
         xf(0) = ra
         do i=1,no
           xf(i) = (((-one*xf(i-1))/a0)/real(2*i,fPrec))*real(2*i-3,fPrec) !hr10
         end do

      elseif(cf.eq.'ISRT') then
!        1/SQRT(A0+P) = 1/SQRT(A0)*(1-1/2(P/A0)+3/8*(P/A0)**2-...)
         if(a0.le.zero) then                                              !hr10
            write(lout,1000) cf,ina,a0
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
         era = one/sqrt(a0)                                             !hr10
         xf(0) = era
         do i=1,no
           xf(i) = (((-one*xf(i-1))/a0)/real(2*i,fPrec))*real(2*i-1,fPrec) !hr10
         end do

      elseif(cf.eq.'EXP ') then
!        EXP(A0+P) = EXP(A0)*(1+P+P**2/2!+...)
         ea  = exp_mb(a0)
         xf(0) = ea
         do i=1,no
           xf(i) = xf(i-1)/real(i,fPrec)
         end do

      elseif(cf.eq.'LOG ') then
!        LOG(A0+P) = LOG(A0) + (P/A0) - 1/2*(P/A0)**2 + 1/3*(P/A0)**3 - ...)
         if(a0.le.zero) then                                              !hr10
            write(lout,1000) cf,ina,a0
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
         ea  = log_mb(a0)
         xf(0) = ea
         xf(1) = one/a0
         do i=2,no
           xf(i) = (((-one*xf(i-1))/a0)/real(i,fPrec))*real((i-1),fPrec)                !hr10
         end do

      elseif(cf.eq.'SIN ') then
!        SIN(A0+P) = SIN(A0)*(1-P**2/2!+P**4/4!) + COS(A0)*(P-P**3/3!+P**5/5!)
         sa  = sin_mb(a0)
         ca  = cos_mb(a0)
         xf(0) = sa
         xf(1) = ca
         do i=2,no
           xf(i) = (-one*xf(i-2))/real(i*(i-1),fPrec)                       !hr10
         end do

      elseif(cf.eq.'COS ') then
!        COS(A0+P) = COS(A0)*(1-P**2/2!+P**4/4!) - SIN(A0)*(P-P**3/3!+P**5/5!)
         sa  = sin_mb(a0)
         ca  = cos_mb(a0)
         xf(0) = ca
         xf(1) = -one*sa                                                 !hr10
         do i=2,no
           xf(i) = (-one*xf(i-2))/real(i*(i-1),fPrec)                      !hr10
         end do

      elseif(cf.eq.'SIRX') then
!        SIN(SQRT(P))/SQRT(P) = 1 - P/3! + P**2/5! - P**3/7! + ...
         if(a0.ne.zero) then                                              !hr10
            write(lout,1000) cf,ina,a0
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
         xf(0)=one
         do i=1,no
           xf(i) = (-one*xf(i-1))/real((2*i)*(2*i+1),fPrec)                !hr10
         end do

      elseif(cf.eq.'CORX') then
!        COS(SQRT(P)) = 1 - P/2! + P**2/4! - P**3/6! + ...
         if(a0.ne.zero) then                                              !hr10
            write(lout,1000) cf,ina,a0
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
         xf(0)=one
         do i=1,no
           xf(i) = (-one*xf(i-1))/real((2*i)*(2*i-1),fPrec)                !hr10
         end do

      elseif(cf.eq.'SIDX') then
!        SIN(P)/P = 1 - P**2/3! + P**4/5! - P**6/7! + ...
         if(a0.ne.zero) then                                              !hr10
            write(lout,1000) cf,ina,a0
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
         xf(0)=one
         xf(1)=zero
         do i=2,no
           xf(i) = (-one*xf(i-2))/real(i*(i+1),fPrec)                      !hr10
         end do

      elseif(cf.eq.'TAN ') then
         if(abs(cos_mb(a0)).lt.epsmac) then
            write(lout,1000) cf,ina,a0
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
         sa  = sin_mb(a0)
         ca  = cos_mb(a0)
         xf(0) = sa/ca
         xf(1) = (one/ca)/ca                                            !hr10
         xf(2) = ((((two*sa)/ca)/ca)/ca)/two
         xf(3) = (((((two*ca**2+six*sa**2)/ca)/ca)/ca)/ca)/six          !hr10
         xf(4) = ((((((16.0_fPrec*sa+eight*sa**3)/ca)/ca)/ca)/ca)/ca)/24.0_fPrec    !hr10
         xf(5) = (((((((((16.0_fPrec*ca**2+(24.0_fPrec*ca**2)*sa**2)+80.0_fPrec*sa**2)+40.0_fPrec*sa**4)/ca)/ca)/ca)/ca)/ca)/ca)/&
     &120.0_fPrec !hr10
         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
            call dadeb(31,'ERR DAFUN ',1)
#ifdef CR
      call abend('                                                  ')
#else
            stop
#endif
         endif
      elseif(cf.eq.'COT ') then
         if(abs(sin_mb(a0)).lt.epsmac) then
            write(lout,1000) cf,ina,a0
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
         sa  = sin_mb(a0)
         ca  = cos_mb(a0)
         xf(0) = ca/sa
         xf(1) = (-one/sa)/sa
         xf(2) = ((((two*ca)/sa)/sa)/sa)/two                           !hr10
         xf(3) = (((((-one*(two*sa**2+six*ca**2))/sa)/sa)/sa)/sa)/six    !hr10
         xf(4) = ((((((16.0_fPrec*ca+eight*ca**3)/sa)/sa)/sa)/sa)/sa)/24.0_fPrec !hr10
         xf(5) = (((((((-one*(((16.0_fPrec*sa**2+(24.0_fPrec*sa**2)*ca**2)+80.0_fPrec*ca**2)+ 40.0_fPrec*ca**4))/sa)/sa)/sa)/sa)/&
     &sa)/sa)/120.0_fPrec         !hr10

         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
            call dadeb(31,'ERR DAFUN ',1)
#ifdef CR
      call abend('                                                  ')
#else
            stop
#endif
         endif
      elseif(cf.eq.'ASIN') then
         if((one-abs(a0)).lt.zero) then
            write(lout,1000) cf,ina,a0
            call dadeb(31,'ERR DAFUN ',1)
            lfun = 0
            return
         endif
         xf(0) = asin_mb(a0)
!hr10 This code is not tested so leave **(-0.5d0) as it is.
!hr10 lf95 opt 1 gives a different result to opt 0 so should be changed to SQRT.
!        xf(1) = (1.d0-a0**2)**(-0.5d0)                                  !hr10
         xf(1) = sqrt(one-a0*a0)                                        !eric
!        xf(2) = (a0*xf(1)**3.d0)/2.d0                                   !hr10
         xf(2) = (a0*(xf(1)*xf(1)*xf(1)))/two                            !eric
!        xf(3) = ((1.d0+2.d0*a0**2)*xf(1)**5.d0)/6.d0                    !hr10
         xf(3) = ((one+two*(a0*a0))*(xf(1)*xf(1)*xf(1)*xf(1)*xf(1)))/six                   !eric
!        xf(4) = ((9.d0*a0+6.d0*a0**3)*xf(1)**7.d0)/24.d0                !hr10
         xf(4) = ((nine*a0+six*(a0*a0*a0))*(xf(1)*xf(1)*xf(1)*xf(1)*xf(1)*xf(1)*xf(1)))/24.0_fPrec !eric
!        xf(5) = ((9.d0+72.d0*a0**2+24.d0*a0**4)*xf(1)**9.d0)/120.d0     !hr10
         xf(5) = ((nine+72.0_fPrec*(a0*a0)+24.0_fPrec*(a0*a0*a0*a0))*(xf(1)*xf(1)*xf(1)*xf(1)*xf(1)*xf(1)*xf(1)*xf(1)*xf(1)))    &
     &/120.0_fPrec !eric
         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
#ifdef CR
      call abend('                                                  ')
#else
            stop
#endif
         endif
      elseif(cf.eq.'ACOS')then
         if((one-abs(a0)).lt.zero) then
            call dadeb(31,'ERR DAFUN ',1)
            write(lout,1000) cf,ina,a0
            lfun = 0
            return
         endif
         xf(0) =  acos_mb(a0)
!hr10 This code is not tested so leave **(-0.5d0) as it is.
!hr10 lf95 opt 1 gives a different result to opt 0 so should be changed to SQRT.
         scr =  (one-a0**2)**(-half)                                   !hr10
         xf(1) =  -one*scr
         xf(2) = ((-one*a0)*scr**three)/two                              !hr10
         xf(3) = ((-one*(one+two*a0**2))*scr**five)/six               !hr10
         xf(4) = ((-one*(nine*a0+six*a0**3))*scr**seven)/24.0_fPrec           !hr10
         xf(5) =((-one*(nine+72.0_fPrec*a0**2+24.0_fPrec*a0**4))*scr**nine)/120.0_fPrec !hr10
         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'ATAN') then
!        ATAN(A0+P) = ATAN(A0)+1/(1+A0**2)*P-A0/(1+A0**2)**2*P**2+....)
         xf(0) = atan_mb(a0)
         xf(1) = one/(one+a0*a0)
         xf(2) = -a0*(xf(1)*xf(1))
         xf(3) = (a0*a0-one/three)*xf(1)**3
         xf(4) = (a0-a0*a0*a0)*xf(1)**4
         xf(5) = (one/five+a0**4-two*a0*a0)*xf(1)**5
         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'ACOT') then
         xf(0) = two*atan_mb(one)-atan_mb(a0)
         scr = one/(one+a0*a0)
         xf(1) = -scr
         xf(2) = a0*(scr*scr)
         xf(3) = -(a0*a0-one/three)*scr**3
         xf(4) = -(a0-a0*a0*a0)*scr**4
         xf(5) = -(one/five+a0**4-two*a0*a0)*scr**5
         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'SINH') then
         sa  = sinh_mb(a0)
         ca  = cosh_mb(a0)
         xf(0) = sa
         xf(1) = ca
         xf(2) = sa/two
         xf(3) = ca/six
         xf(4) = sa/24.0_fPrec
         xf(5) = ca/120.0_fPrec
         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'COSH') then
         sa  = sinh_mb(a0)
         ca  = cosh_mb(a0)
         xf(0) = ca
         xf(1) = sa
         xf(2) = ca/two
         xf(3) = sa/six
         xf(4) = ca/24.0_fPrec
         xf(5) = sa/120.0_fPrec
         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'TANH') then
         sa  = sinh_mb(a0)
         ca  = cosh_mb(a0)
         xf(0) = sa/ca
         xf(1) = one/ca/ca
         xf(2) = -two*sa/ca/ca/ca/two
         xf(3) = (-two*ca*ca+six*sa*sa)/ca/ca/ca/ca/six
         xf(4) =(16.0_fPrec*sa-eight*sa*sa*sa)/ca/ca/ca/ca/ca/24.0_fPrec
         xf(5) = (16.0_fPrec*ca*ca-24.0_fPrec*ca*ca*sa*sa-80.0_fPrec*sa*sa+40.0_fPrec*sa*sa*sa*sa)/ca/ca/ca/ca/ca/ca/120.0_fPrec
         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'COTH') then
         if(abs(sinh_mb(a0)).lt.epsmac) then
            lfun = 0
            write(lout,1000) cf,ina,a0
            call dadeb(31,'ERR DAFUN ',1)
            return
         endif
         sa  = sinh_mb(a0)
         ca  = cosh_mb(a0)
         xf(0) = ca/sa
         xf(1) = -one/sa/sa
         xf(2) =  two*ca/sa/sa/sa/two
         xf(3) = (two*sa*sa-six*ca*ca)/sa/sa/sa/sa/six
         xf(4) = (16.0_fPrec*ca+eight*ca*ca*ca)/sa/sa/sa/sa/sa/24.0_fPrec
         xf(5) = (16.0_fPrec*sa*sa+24.0_fPrec*sa*sa*ca*ca-80.0_fPrec*ca*ca-40.0_fPrec*ca*ca*ca*ca)/sa/sa/sa/sa/sa/sa/120.0_fPrec
         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'ASNH') then
         xf(0) = log_mb(a0+sqrt(a0*a0+one))
         xf(1) = (one+a0*a0)**(-half)
         xf(2) = -a0*xf(1)**three/two
         xf(3) = (two*a0*a0-one)*xf(1)**five/six
         xf(4) = (nine*a0-six*a0*a0*a0)*xf(1)**seven/24.0_fPrec
         xf(5) = (nine-72.0_fPrec*a0*a0+24.0_fPrec*a0*a0*a0*a0)*xf(1)**nine/120.0_fPrec
         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'ACSH') then
         if((one-a0).ge.zero) then
            lfun = 0
            write(lout,1000) cf,ina,a0
            call dadeb(31,'ERR DAFUN ',1)
            return
         endif
         xf(0) = log_mb(a0+sqrt(a0*a0-one))
         xf(1) = (a0*a0-one)**(-half)
         xf(2) = -a0*xf(1)**three/two
         xf(3) = (two*a0*a0+one)*xf(1)**five/six
         xf(4) = (-nine*a0-six*a0*a0*a0)*xf(1)**seven/24.0_fPrec
         xf(5) = (nine+72.0_fPrec*a0*a0+24.0_fPrec*a0*a0*a0*a0)*xf(1)**nine/120.0_fPrec
         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'ATNH') then
         if((abs(a0)-one).ge.zero) then
            lfun = 0
            write(lout,1000) cf,ina,a0
            call dadeb(31,'ERR DAFUN ',1)
            return
         endif
         xf(0) =  half*log_mb((1+a0)/(1-a0))
         xf(1) =  one/(one-a0*a0)
         xf(2) =  a0*(xf(1)*xf(1))
         xf(3) = (a0*a0+one/three)*xf(1)**3
         xf(4) = (a0+a0*a0*a0)*xf(1)**4
         xf(5) = (one/five+a0**4+two*a0*a0)*xf(1)**5
         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
            call dadeb(31,'ERR DAFUN ',1)
         endif
      elseif(cf.eq.'ACTH') then
         if(one-abs(a0).ge.zero) then
            lfun = 0
            write(lout,1000) cf,ina,a0
            call dadeb(31,'ERR DAFUN ',1)
            return
         endif
         xf(0) =  half*log_mb((a0+1)/(a0-1))
         scr =  one/(-one+a0*a0)
         xf(1) = -scr
         xf(2) =  a0*(scr*scr)
         xf(3) = (-a0*a0-one/three)*scr**three
         xf(4) = (a0+a0*a0*a0)*scr**four
         xf(5) = (-one/five-a0**4-two*a0*a0)*scr**five
         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
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
         e1 = exp_mb(-a0*a0)
         a1 = .254829592_fPrec
         a2 = -.284496736_fPrec
         a3 = 1.421413741_fPrec
         a4 = -1.453152027_fPrec
         a5 = 1.061405429_fPrec
         p  = .3275911_fPrec
         rpi4 = sqrt(atan_mb(one))
         t  = one/(one+p*a0)
         e2 = one-t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))*e1
         xf(0)= e2
         xf(1) = e1/rpi4
         xf(2) = -a0*e1/rpi4
         xf(3) = (-two+four*a0*a0)/six*e1/rpi4
         xf(4) = (12.0_fPrec*a0-eight*a0*a0*a0)/24.0_fPrec*e1/rpi4
         xf(5) = (16.0_fPrec*a0*a0*a0*a0-48.0_fPrec*a0*a0+12.0_fPrec)/120.0_fPrec*e1/rpi4
         if(no.gt.5) then
            write(lout,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
            call dadeb(31,'ERR DAFUN ',1)
         endif
      else
         write(lout,*)'ERROR, UNSOPPORTED FUNCTION ',cf
      endif

      call dacon(inc,xf(0))
      call dacop(ina,inon(1))
      call dapok(inon(1),jj,zero)
      call dacon(ipow(1),one)

      do 800 i=1,min(no,nocut)

      call damul(inon(1),ipow(1),iscr(1))
      call dacop(iscr(1),ipow(1))
      call dacma(inc,ipow(1),xf(i),inc)

 800  continue

 1000 format('DAFUN> ERROR ',a4,' does not exist for vector ',i0,' const term = ',e12.5)

      call dadal(iscr(1),1)
      call dadal(inon(1),1)
      call dadal(ipow(1),1)

      return
end subroutine dafunt


subroutine daabs(ina,anorm)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : cc
      implicit none
      integer i,illa,ilma,ina,inoa,inva,ipoa
      real(kind=fPrec) anorm
!     ***************************
!
!     THIS SUBROUTINE COMPUTES THE NORM OF THE DA VECTOR A
!
!-----------------------------------------------------------------------------1


      call dainf(ina,inoa,inva,ipoa,ilma,illa)

      anorm = zero
      do i=ipoa,ipoa+illa-1
        anorm = anorm + abs(cc(i))
      end do

      return
      end


subroutine dacom(ina,inb,dnorm)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      implicit none
      integer idacom(1),illc,ilmc,ina,inb,inc,inoc,invc,ipoc
      real(kind=fPrec) dnorm
!     *******************************
!
!     THIS SUBROUTINE COMPARES TWO DA VECTORS BY RETURNING THE NORM
!     OF THE DIFFERENCE

      idacom(1) = 0
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
      call daall(idacom(1),1,'$$DACOM $$',inoc,invc)
      call dasub(ina,inb,idacom(1))
      call daabs(idacom(1),dnorm)
      call dadal(idacom(1),1)

      return
      end


subroutine dapos(ina,inb)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : cc,nocut,idalm,idall,i1,i2,ieo,ia1,ia2
      implicit none
      integer ia,ib,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,    &
     &ipoa,ipob
!     *************************
!
!     THIS SUBROUTINE MAKES THE SIGNS OF ALL THE COEFFICIENTS OF A POSITIVE
!
!-----------------------------------------------------------------------------1


      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)

!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)

      ib = ipob - 1

      do 100 ia = ipoa,ipoa+illa-1

      if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 100
      ib     = ib + 1
      cc(ib) = abs(cc(ia))
      i1(ib) = i1(ia)
      i2(ib) = i2(ia)

 100  continue

      idall(inb) = ib - ipob + 1
      if(idall(inb).gt.idalm(inb)) then
         write(lout,*)'ERROR IN DAPOS '
         call dadeb(31,'ERR DAPOS ',1)
      endif

      return
      end


subroutine dacct(ma,ia,mb,ib,mc,ic)
      use floatPrecision
      use mod_lie_dab, only : lnv
      implicit none
      integer i,ia,ib,ic,ij,illc,ilmc,inoc,invc,ipoc
!     ***********************************
!
!     THIS SUBROUTINE PERFORMS A CONCATENATION MA = MB o MC
!     WHERE MA, MB AND MC ARE MATRICES CONSISTING OF IA, IB AND IC
!     DA VECTORS, RESPECTIVELY.
!
!-----------------------------------------------------------------------------1


      integer mon(lnv),mb(*),mc(*),ma(*)

      if(ma(1).eq.mc(1).or.mb(1).eq.mc(1)) then
        call dainf(mc(1),inoc,invc,ipoc,ilmc,illc)

        do ij=1,ic
          mon(ij)=0
        end do

        call daall(mon,ic,'$$DAJUNK$$',inoc,invc)
        call dacctt(ma,ia,mb,ib,mon,ic)

        do i=1,ic
          call dacop(mon(i),mc(i))
        end do

        call dadal(mon,ic)
      else
        call dacctt(ma,ia,mb,ib,mc,ic)
      endif

      return
      end


subroutine dacctt(mb,ib,mc,ic,ma,ia)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : cc,eps,nomax,nvmax,idapo,idall,i1,i2,lno
      implicit none
      integer i,ia,ib,ic,iia,iib,iic,illa,illb,illc,ilma,ilmb,ilmc,inoa,&
     &inob,inoc,inva,invb,invc,ipoa,ipob,ipoc,iv,jl,jv
      real(kind=fPrec) ccf
!     ***********************************
!
!     THIS SUBROUTINE PERFORMS A CONCATENATION MA = MB o MC
!     WHERE MA, MB AND MC ARE MATRICES CONSISTING OF IA, IB AND IC
!     DA VECTORS, RESPECTIVELY.
!
!-----------------------------------------------------------------------------1

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
         write(lout,*)'ERROR IN DACCT, IA .NE. IB'
         call dadeb(31,'ERR DACCT1',1)
      elseif(ic.ne.invb) then
         write(lout,*)'ERROR IN DACCT, IC.NE.INVB'
         call dadeb(31,'ERR DACCT2',1)
      endif

!     ALLOCATING LOCAL VECTORS AND CALLING MTREE
!     ******************************************

      do i=1,ib
        icc(i) = 0
      end do

      do i=1,nomax+1
        mon(i) = 0
      end do

      call daall(icc,ib,'$$DACCT $$',nomax,nvmax)
      call daall(mon,nomax+1,'$$DAMON $$',inoc,invc)

      call mtree(mb,ib,icc,ib)

!     PERFORMING CONCATENATION
!     ************************

      do i=1,ia
        call dacon(ma(i),cc(idapo(icc(i))))
      end do

      call dacon(mon(1),one)

      do i=1,idall(icc(1))-1
        jl = i1(idapo(icc(1))+i)
        jv = i2(idapo(icc(1))+i)

        call damul(mon(jl),mc(jv),mon(jl+1))

        do iv=1,ia
          ccf = cc(idapo(icc(iv))+i)
          if(abs(ccf).gt.eps) call dacma(ma(iv),mon(jl+1),ccf,ma(iv))
        end do
      end do

      call dadal(mon,nomax+1)
      call dadal(icc,ib)

      return
      end


subroutine mtree(mb,ib,mc,ic)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : cc,epsmac,nomax,nvmax,nmmax,idapo,idalm,idall,i1,i2,ieo,ia1,ia2,lno,lnv
      implicit none
      integer i,ib,ib1,ibi,ic,ic1,ic2,icc,ichk(1),ii,iib,iic,illb,illc, &
     &ilmb,ilmc,inob,inoc,invb,invc,ipob,ipoc,j,jl,jnon,nterm,ntermf
      real(kind=fPrec) apek,bbijj,chkjj
!     *****************************
!
!     THIS SUBROUTINE IS USED FOR CONCATENATION AND TRACKING OF VECTORS
!     THROUGH A DA MAP. IT COMPUTES THE TREE THAT HAS TO BE TRANSVERSED
!     MB IS THE DA MATRIX WITH IA TERMS. THE OUTPUT MC IS A CA MATRIX WHICH
!     CONTAINS COEFFICIENTS AND CONTROL INTEGERS USED FOR THE TRAVERSAL.
!
!-----------------------------------------------------------------------------1


      integer jj(lnv),jv(0:lno),mb(*),mc(*)

!     CONSISTENCY CHECKS
!     ******************

      iib = mb(1)
      iic = mc(1)
      call dainf(iib,inob,invb,ipob,ilmb,illb)
      call dainf(iic,inoc,invc,ipoc,ilmc,illc)

      call damch(mb,ib)
      call damch(mc,ic)

      if(ib.ne.ic) then
         write(lout,*)'ERROR IN MTREE, IB .NE. IC'
         call dadeb(31,'ERR MTREE1',1)
      endif

!     ALLOCATING LOCAL VECTORS
!     ************************

      ichk(1) = 0
      call daall(ichk(1),1,'$$MTREE $$',nomax,nvmax)

!     FIND ALL THE ENTRIES TO BE LOOKED FOR
!     *************************************

      call daclr(1)

      cc(1) = one

      do 100 i=1,ib
      if(nomax.eq.1) then
      do 91 ib1 = 2,7
      cc(ib1) = one
   91 continue
      else
      do  90 ibi = idapo(mb(i)),idapo(mb(i))+idall(mb(i))-1
      icc = ia1(i1(ibi)) + ia2(i2(ibi))
      if(ieo(icc).gt.inob) goto 90
      cc(icc) = one
   90 continue
      endif
  100 continue

      do 150 ii=1,inob

!     SEARCHING FOR FATHER FOR EACH TERM

      do 140 i=1,nmmax
      if(cc(i).lt.half) goto 140

      jnon = 0
      call dancd(i1(i),i2(i),jj)
      do 130 j=1,invb
      if(jj(j).eq.0) goto 130
      jnon = j
      jj(j) = jj(j) - 1
      call dadcd(jj,ic1,ic2)
      apek = cc(ia1(ic1)+ia2(ic2))
      jj(j) = jj(j) + 1
      if(apek.ge.half) goto 140
  130 continue

      if(jnon.eq.0) goto 140

!     TERM IS AN ORPHAN, SO CREATE FOSTER FATHER

      jj(jnon) = jj(jnon) - 1
      call dadcd(jj,ic1,ic2)
      cc(ia1(ic1)+ia2(ic2)) = one

  140 continue
  150 continue

      call dapac(ichk(1))
!ETIENNE      CALL DAPRI(ICHK,32)
!
!     SETTING UP TREE STRUCTURE
!     *************************

      ntermf = idall(ichk(1))

!     ZEROTH ORDER TERMS
!     ******************

      do i=1,lnv
        jj(i) = 0
      end do

      do 170 i=1,ib
      call dapek(mb(i),jj,bbijj)
      i1(idapo(mc(i))) = 0
      i2(idapo(mc(i))) = 0
      cc(idapo(mc(i))) = bbijj
 170  continue

      call dapek(ichk(1),jj,chkjj)
      if(chkjj.gt.half) then
         call dapok(ichk(1),jj,-one)
      else
         write(lout,*) 'ERROR IN MTREE, ZEROTH ORDER TERM OF ICHK IS ZERO'
         call dadeb(31,'ERR MTREE2',1)
      endif

      nterm = 1

!     HIGHER ORDER TERMS
!     ******************

      do jl=1,inob
        jv(jl) = 0
      end do

      jl = 0
      chkjj = one

 200  continue
      if(jl.eq.0.and.chkjj.le.half) goto 250
      if(jl.lt.inob.and.chkjj.gt.half) then
         jl = jl + 1
         jj(1) = jj(1) + 1
         jv(jl) = 1
      elseif(jv(jl).eq.invb) then
         jj(jv(jl)) = jj(jv(jl)) - 1
         jv(jl) = 0
         jl = jl - 1
         chkjj = zero
         goto 200
      else
         jj(jv(jl)) = jj(jv(jl)) - 1
         jv(jl) = jv(jl) + 1
         jj(jv(jl)) = jj(jv(jl)) + 1
      endif

      call dapek(ichk(1),jj,chkjj)

      if(chkjj.le.half) goto 200

      nterm = nterm + 1
      if(nterm.gt.idalm(mc(1))) then
         write(lout,*)'ERROR IN MTREE, NTERM TOO LARGE'
         call dadeb(31,'ERR MTREE3',1)
      endif

      call dapok(ichk(1),jj,-one)

      do 210 i=1,ib
      call dapek(mb(i),jj,bbijj)
      i1(idapo(mc(i))+nterm-1) = jl
      i2(idapo(mc(i))+nterm-1) = jv(jl)
      cc(idapo(mc(i))+nterm-1) = bbijj
 210  continue

      goto 200

 250  continue

      do i=1,ib
        idall(mc(i)) = nterm
      end do

!     PERFORMING CROSS CHECKS
!     ***********************

      if(nterm.ne.ntermf.or.nterm.ne.idall(ichk(1))) then
         write(lout,*)'ERROR IN MTREE, NTERM, NTERMF, IDALL(ICHK) =  ',nterm,ntermf,idall(ichk(1))
         call dadeb(31,'ERR MTREE4',1)
      endif

      do 270 i=idapo(ichk(1)),idapo(ichk(1))+nterm-1
      if(abs(cc(i)+one).gt.epsmac) then
         write(lout,*)'ERROR IN MTREE, NOT ALL TERMS IN ICHK ARE -1'
         call dadeb(31,'ERR MTREE5',1)
      endif
 270  continue

      call dadal(ichk(1),1)

      return
      end


subroutine ppushpri(mc,ic,mf,jc,line)
      use floatPrecision
      use mod_lie_dab, only : cc,idapo,idall,i1,i2
      implicit none
      integer i,ic,iv,jc,jl,jv,mc,mf

      dimension mc(*)
      character(len=20) line
      if(mf.le.0) return
      write(mf,*) 0,0,jc+1,0,line

      do i=1,ic
        jc=1+jc
        write(mf,*) jc,jl,jv,cc(idapo(mc(i)))
      end do

!     xf(i) = cc(idapo(mc(i)))
!      xm(1) = 1.d0
      do i=1,idall(mc(1))-1
        jl = i1(idapo(mc(1))+i)
        jv = i2(idapo(mc(1))+i)
!       xx = xm(jl)*xi(jv)
!       xm(jl+1) = xx

        do iv=1,ic
          jc=1+jc
          write(mf,*) jc,jl,jv,cc(idapo(mc(iv))+i)
!         xf(iv) = xf(iv) + cc(idapo(mc(iv))+i) * xx
        end do
      end do

      return
      end


subroutine ppush(mc,ic,xi,xf)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : cc,idapo,idall,i1,i2,lno
      implicit none
      integer i,ic,iv,jl,jv,mc
      real(kind=fPrec) xf,xi,xm,xt,xx
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE MATRIX WHOSE TREE IS STORED IN CA VECTOR MC
!     TO THE COORDINATES IN XI AND STORES THE RESULT IN XF
!
!-----------------------------------------------------------------------------1


      dimension mc(*),xf(*),xi(*),xm(lno+1) ,xt(lno)

      do i=1,ic
        xt(i)=xi(i)
      end do

      do i=1,ic
        xf(i) = cc(idapo(mc(i)))
      end do

      xm(1) = one

      do i=1,idall(mc(1))-1
        jl = i1(idapo(mc(1))+i)
        jv = i2(idapo(mc(1))+i)
        xx = xm(jl)*xt(jv)
        xm(jl+1) = xx

        do iv=1,ic
          xf(iv) = xf(iv) + cc(idapo(mc(iv))+i) * xx
        end do
      end do

      return
      end


subroutine ppush1(mc,xi,xf)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : cc,nvmax,idapo,idall,i1,i2,lno
      implicit none
      integer i,jl,jv,mc
      real(kind=fPrec) xf,xi,xm,xt,xx
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE MATRIX WHOSE TREE IS STORED IN CA VECTOR MC
!     TO THE COORDINATES IN XI AND STORES THE RESULT IN XF
!
!-----------------------------------------------------------------------------1


      dimension xi(*),xm(lno+1) ,xt(lno)

      do i=1,nvmax
        xt(i)=xi(i)
      end do

      xf = cc(idapo(mc))

      xm(1) = one

      do i=1,idall(mc)-1
        jl = i1(idapo(mc)+i)
        jv = i2(idapo(mc)+i)
        xx = xm(jl)*xt(jv)
        xm(jl+1) = xx
        xf = xf + cc(idapo(mc)+i) * xx
      end do

      return
      end


subroutine dainv(ma,ia,mb,ib)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : lnv
      implicit none
      integer i,ia,ib,ij,illb,ilmb,inob,invb,ipob
      real(kind=fPrec) x
!     *****************************
!
!     THIS SUBROUTINE INVERTS THE MATRIX MA WITH IA DA VECTORS AND
!     STORES THE RESULT IN MI
!
!-----------------------------------------------------------------------------1


      integer jj(lnv),ml(lnv),ma(*),mb(*)
      dimension x(lnv)

      do i=1,lnv
        jj(i)=0
      end do

      if(ma(1).eq.mb(1)) then
        call dainf(mb(1),inob,invb,ipob,ilmb,illb)
        do i=1,ia
          call dapok(ma(i),jj,zero)
        end do

        do ij=1,ib
          ml(ij)=0
        end do

        call daall(ml,ib,'$$DAJUNK$$',inob,invb)
        call dainvt(ma,ia,ml,ib)

        do i=1,ib
          call dacop(ml(i),mb(i))
        end do

        call dadal(ml,ib)
      else
        do i=1,ia
          call dapek(ma(i),jj,x(i))
          call dapok(ma(i),jj,zero)
        end do

        call dainvt(ma,ia,mb,ib)
        do i=1,ia
          call dapok(ma(i),jj,x(i))
        end do
      endif

      return
      end


subroutine dainvt(ma,ia,mb,ib)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : eps,epsmac,nocut,lnv
      implicit none
      integer i,ia,ib,ie,ier,illa,illb,ilma,ilmb,inoa,inob,inva,invb,   &
     &ipoa,ipob,j,k,nocut0
      real(kind=fPrec) aa,ai,amjj,amsjj,prod
!     *****************************
!
!     THIS SUBROUTINE INVERTS THE MATRIX MA WITH IA DA VECTORS AND
!     STORES THE RESULT IN MI
!
!-----------------------------------------------------------------------------1


      integer jj(lnv),ms(lnv),ml(lnv),ma(*),mb(*)
      dimension aa(lnv,lnv),ai(lnv,lnv)

      call dainf(ma(1),inoa,inva,ipoa,ilma,illa)
      call dainf(mb(1),inob,invb,ipob,ilmb,illb)

!     CONSISTENCY CHECKS
!     ******************

      call damch(ma,ia)
      call damch(mb,ib)
!etienne
      do ie=1,ib
        call dacon(mb(ie),zero)
      end do
!etienne

      if(ia.ne.ib) then
         write(lout,*)'ERROR IN DAINV, IA .NE. IB'
         call dadeb(31,'ERR DAINV1',1)
      elseif(ia.ne.inva.or.ib.ne.invb) then
         write(lout,*)'ERROR IN DAINV, IA.NE.INVA.OR.IB.NE.INVB'
         call dadeb(31,'ERR DAINV2',1)
      endif

!     ALLOCATING LOCAL VECTORS
!     ************************

      do 10 i=1,ia
      ms(i) = 0
      ml(i) = 0
  10  continue

      call daall(ms,ia,'$$INV   $$',inoa,inva)
      call daall(ml,ia,'$$INVL  $$',inoa,inva)

!     EXTRACTING LINEAR MATRIX, GENERATING NONLINEAR PART OF A
!     ********************************************************

      do i=1,ib
        do j=1,ib
          do k=1,ib
            jj(k) = 0
          end do

          jj(j) = 1
          call dapek(ma(i),jj,amjj)
          if(abs(amjj).gt.eps) call dapok(ma(i),jj,zero)
          aa(i,j) = amjj
        end do
        call dacmu(ma(i),-one,ma(i))
      end do

!     INVERTING LINEAR MATRIX, CHECKING RESULT AND STORING IN ML
!     **********************************************************

      call matinv(aa,ai,ia,lnv,ier)

      if(ier.eq.132) then
         write(lout,*)'ERROR IN ROUTINE DAINV'
         call dadeb(31,'ERR DAINV3',1)
      endif

      ier = 0

      do i=1,ib
        do j=1,ib
        prod = zero

        do k=1,ib
          jj(k) = 0
          prod = prod + aa(i,k)*ai(k,j)
        end do

        if(i.eq.j) prod = prod - one

        if(abs(prod).gt.100*epsmac) then
          write(lout,*) 'ERROR IN DAINV, INVERSION DID NOT WORK,I,J,PROD = ',i,j,prod,epsmac,eps
          ier = 1
!ETIENNE
          return
!ETIENNE
        endif

        jj(j) = 1
        call dapok(mb(i),jj,ai(i,j))
        call dapok(ml(i),jj,ai(i,j))
        end do
      end do

      if(ier.eq.1) call dadeb(31,'ERR DAINV4',1)
!
!     ITERATIVELY COMPUTING DIFFERENT PARTS OF THE INVERSE
!     ****************************************************
!
!     MB (OF ORDER I) = A1^-1 o [ E - ANL (NONLINEAR) o MB (OF ORDER I) ]

      nocut0 = nocut

      do i=2,nocut
        nocut = i
        call dacct(ma,ia,mb,ib,ms,ia)
        do j=1,ib
          do k=1,ib
            jj(k) = 0
          end do
          jj(j) = 1
          call dapek(ms(j),jj,amsjj)
          call dapok(ms(j),jj,amsjj+one)
        end do
        call dacct(ml,ia,ms,ia,mb,ib)
      end do

      nocut = nocut0

!     FLIPPING BACK SIGN OF A, FILLING UP FIRST ORDER PART AGAIN
!     **********************************************************

      do i=1,ib
        call dacmu(ma(i),-one,ma(i))
        do j=1,ib
          do k=1,ib
            jj(k) = 0
          end do
        jj(j) = 1
        call dapok(ma(i),jj,aa(i,j))
        end do
      end do

      call dadal(ml,ia)
      call dadal(ms,ia)

      return
      end


subroutine matinv(a,ai,n,nmx,ier)
      use floatPrecision
      use numerical_constants
      implicit none
      integer i,ier,indx,j,n,nmax,nmx
      real(kind=fPrec) a,ai,aw,d
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

      do i=1,n
        do j=1,n
          aw(i,j) = a(i,j)
          ai(i,j) = zero
        end do

        ai(i,i) = one
      end do

      call ludcmp(aw,n,nmax,indx,d,ier)
      if (ier .eq. 132) return

      do j=1,n
        call lubksb(aw,n,nmax,indx,ai(1,j),nmx)
      end do

      return
      end


subroutine ludcmp(a,n,np,indx,d,ier)
      use floatPrecision
      use numerical_constants
      implicit none
      integer i,ier,imax,indx,j,k,n,nmax,np
      real(kind=fPrec) a,aamax,d,dum,sum,tiny,vv
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
      parameter (nmax = 400, tiny = 1.0e-20_fPrec)
      dimension a(np,np), indx(np), vv(nmax)
      ier=0.
      d=one
      imax = 0 ! -Wmaybe-uninitialized
      do 12 i=1,n
         aamax=zero
         do 11 j=1,n
            if(abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11       continue
         if(aamax.eq.zero) then
            ier=132
            return
         endif
         vv(i)=one/aamax
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
         aamax=zero
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
            if(a(j,j).eq.zero) a(j,j)=tiny
            dum=1./a(j,j)
            do 18 i=j+1,n
               a(i,j)=a(i,j)*dum
18          continue
         endif
19    continue
      if(a(n,n).eq.zero) a(n,n)=tiny
      return
      end


subroutine lubksb(a,n,np,indx,b,nmx)
      use floatPrecision
      use numerical_constants
      implicit none
      integer i,ii,indx,j,ll,n,nmx,np
      real(kind=fPrec) a,b,sum
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
         else if (sum.ne.zero) then
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


subroutine dapin(ma,ia,mb,ib,jx)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : lnv
      implicit none
      integer i,ia,ib,ij,illb,ilmb,inob,invb,ipob
      real(kind=fPrec) x
!     *****************************
!
!     THIS SUBROUTINE INVERTS THE MATRIX MA WITH IA DA VECTORS AND
!     STORES THE RESULT IN MI
!
!-----------------------------------------------------------------------------1


      integer jj(lnv),ml(lnv),ma(*),mb(*),jx(*)
      dimension x(lnv)

      do i=1,lnv
        jj(i)=0
      end do

      if(ma(1).eq.mb(1)) then
        call dainf(mb(1),inob,invb,ipob,ilmb,illb)

        do i=1,ia
          call dapok(ma(i),jj,zero)
        end do

        do ij=1,ib
          ml(ij)=0
        end do

        call daall(ml,ib,'$$DAJUNK$$',inob,invb)
        call dapint(ma,ia,ml,ib,jx)

        do i=1,ib
          call dacop(ml(i),mb(i))
        end do

        call dadal(ml,ib)
      else

        do i=1,ia
          call dapek(ma(i),jj,x(i))
          call dapok(ma(i),jj,zero)
        end do

        call dapint(ma,ia,mb,ib,jx)

        do i=1,ia
          call dapok(ma(i),jj,x(i))
        end do

      endif

      return
      end


subroutine dapint(ma,ia,mb,ib,jind)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : nvmax,lnv
      implicit none
      integer i,ia,ib,illa,ilma,inoa,inva,ipoa,k
!     **********************************
!
!     THIS SUBROUTINE PERFORMS A PARTIAL INVERSION OF THE ROWS MARKED WITH
!     NONZERO ENTRIES IN JJ OF THE MATRIX A. THE RESULT IS STORED IN B.
!
!-----------------------------------------------------------------------------1


      integer jj(lnv),jind(*),ma(*),mb(*),mn(lnv),mi(lnv),me(lnv)

      call dainf(ma(1),inoa,inva,ipoa,ilma,illa)

      do i=1,ia
        mn(i) = 0
        mi(i) = 0
        me(i) = 0
      end do

      call daall(mn,ia,'$$PIN1  $$',inoa,inva)
      call daall(mi,ia,'$$PIN2  $$',inoa,inva)
      call daall(me,ia,'$$PIN3  $$',inoa,inva)

      do i=1,ia
        do k=1,nvmax
          jj(k) = 0
        end do
        jj(i) = 1
        call dapok(me(i),jj,one)
      end do

      do i=1,ia
        call dacop(ma(i),mn(i))
        if(jind(i).eq.0) call dacop(me(i),mn(i))
      end do

      call dainv(mn,ia,mi,ia)

      do i=1,ia
        if(jind(i).eq.0) call dacop(ma(i),me(i))
      end do

      call dacct(me,ia,mi,ia,mb,ib)

      call dadal(me,ia)
      call dadal(mi,ia)
      call dadal(mn,ia)

      return
      end


subroutine dader(idif,ina,inc)
      use floatPrecision
      implicit none
      integer idif,illc,ilmc,ina,inc,incc(1),inoc,invc,ipoc
!     ******************************
!
!     THIS SUBROUTINE COMPUTES THE DERIVATIVE WITH RESPECT TO VARIABLE I
!     OF THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1



      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc(1)=0
        call daall(incc(1),1,'$$DAJUNK$$',inoc,invc)
        call dadert(idif,ina,incc(1))
        call dacop(incc(1),inc)
        call dadal(incc(1),1)
      else
        call dadert(idif,ina,inc)
      endif

      return
      end


subroutine dadert(idif,ina,inc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : cc,epsmac,nomax,nvmax,idalm,idall,i1,i2,lnv
      implicit none
      integer i,ibase,ic,ider1,ider1s,ider2,ider2s,idif,iee,ifac,illa,  &
     &illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,ipoa,ipoc,jj
      real(kind=fPrec) rr,x,xdivi
!     ******************************
!
!     THIS SUBROUTINE COMPUTES THE DERIVATIVE WITH RESPECT TO VARIABLE I
!     OF THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1

      integer jd(lnv)
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!
      if(nomax.eq.1) then
!         PRINT*,'ERROR, DADER CALLED WITH NOMAX = 1'
!        CALL DADEB(31,'ERR DADER1',1)
#ifdef CR
!     call abend('666                                               ')
#else
!        stop 666
#endif
      do i=1,lnv
        jd(i)=0
        enddo
        jd(idif)=1
        call dapek(ina,jd,rr)
        call dacon(inc,rr)
        return
      endif

!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)

      ibase = nomax + 1

      if(idif.gt.(nvmax+1)/2) then
         ider1  = 0
         ider1s = 0
         ider2  = idif-(nvmax+1)/2
         ider2s = 1
         do jj=1,ider2-1
           ider2s = ider2s*ibase
         end do
         xdivi  = ider2s*ibase
      else
         ider1  = idif
         ider1s = 1
         do jj=1,ider1-1
           ider1s = ider1s*ibase
         end do
         ider2  = 0
         ider2s = 0
         xdivi  = ider1s*ibase
      endif

      ibase = nomax+1

      ic = ipoc-1

      do 100 i=ipoa,ipoa+illa-1

      if(ider1.eq.0) then
         iee = i2(i)
      else
         iee = i1(i)
      endif

      x = iee/xdivi
      ifac = int(ibase*(x-int(x+epsmac)+epsmac))

      if(ifac.eq.0) goto 100

      ic = ic + 1
      cc(ic) = cc(i)*ifac
      i1(ic) = i1(i) - ider1s
      i2(ic) = i2(i) - ider2s

 100  continue

      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
         write(lout,*)'ERROR IN DADER '
         call dadeb(31,'ERR DADER2',1)
      endif

      return
      end


subroutine dapoi(ina,inb,inc,n)
      use floatPrecision
      use mod_lie_dab, only : nomax,nvmax
      implicit none
      integer i,ina,inb,inc,n
!     *******************************
!
!     THIS SUBROUTINE COMPUTES THE POISSON BRACKET OF THE VECTORS A AND
!     B AND STORES THE RESULT IN C. N IS THE DEGREE OF FREEDOM OF THE SYSTEM.
!
!-----------------------------------------------------------------------------1

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

subroutine dacfur(ina,fun,inc)
      use floatPrecision
      implicit none
      integer illc,ilmc,ina,inc,inoc,invc,ipoc
      integer incc(1)
      complex(kind=fPrec) fun
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1

!
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc(1)=0
        call daall(incc(1),1,'$$DAJUNK$$',inoc,invc)
        call dacfurt(ina,fun,incc(1))
        call dacop(incc(1),inc)
        call dadal(incc(1),1)
      else
        call dacfurt(ina,fun,inc)
      endif

      return
      end

subroutine dacfurt(ina,fun,inc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : cc,eps,nomax,idalm,idall,i1,i2,lnv
      implicit none
      integer i,ia,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,ipoa,ipoc,j
      real(kind=fPrec) cfac,rr
      complex(kind=fPrec) fun
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1

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
      cfac = real(fun(j))
      rr=cfac*rr
      call dapok(inc,j,rr)
      do i=1,lnv
        j(i)=1
        call dapek(ina,j,rr)
        cfac = real(fun(j))
        rr=cfac*rr
        call dapok(inc,j,rr)
        j(i)=0
      enddo
!         PRINT*,'ERROR, DACFU CALLED WITH NOMAX = 1'
!        CALL DADEB(31,'ERR DACFU ',1)
#ifdef CR
!     call abend('667                                               ')
#else
!        stop 667
#endif
      return
      endif
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      ic = ipoc - 1
!
      do ia=ipoa,ipoa+illa-1
!
      call dancd(i1(ia),i2(ia),j)
      cfac = real(fun(j))
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
         write(lout,*)'ERROR IN DACFU '
         call dadeb(31,'ERR DACFU ',1)
      endif
!
      return
      end
!

subroutine dacfu(ina,fun,inc)
      use floatPrecision
      implicit none
      integer illc,ilmc,ina,inc,inoc,invc,ipoc
      integer incc(1)
      real(kind=fPrec) fun
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE PRECISION FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1

!
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc(1)=0
        call daall(incc(1),1,'$$DAJUNK$$',inoc,invc)
        call dacfut(ina,fun,incc(1))
        call dacop(incc(1),inc)
        call dadal(incc(1),1)
      else
        call dacfut(ina,fun,inc)
      endif

      return
      end

subroutine dacfui(ina,fun,inc)
      use floatPrecision
      implicit none
      integer illc,ilmc,ina,inc,inoc,invc,ipoc
      complex(kind=fPrec) fun
      integer incc(1)
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1

!
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc(1)=0
        call daall(incc(1),1,'$$DAJUNK$$',inoc,invc)
        call dacfuit(ina,fun,incc(1))
        call dacop(incc(1),inc)
        call dadal(incc(1),1)
      else
        call dacfuit(ina,fun,inc)
      endif

      return
      end

subroutine dacfuit(ina,fun,inc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : cc,eps,nomax,idalm,idall,i1,i2,lnv
      implicit none
      integer i,ia,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,ipoa,ipoc,j
      real(kind=fPrec) cfac,rr
      complex(kind=fPrec) fun
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1

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
      cfac = aimag(fun(j))
      rr=cfac*rr
      call dapok(inc,j,rr)
      do i=1,lnv
        j(i)=1
        call dapek(ina,j,rr)
        cfac = aimag(fun(j))
        rr=cfac*rr
        call dapok(inc,j,rr)
        j(i)=0
      enddo
!         PRINT*,'ERROR, DACFU CALLED WITH NOMAX = 1'
!        CALL DADEB(31,'ERR DACFU ',1)
#ifdef CR
!     call abend('667                                               ')
#else
!        stop 667
#endif
      return
      endif
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)

      ic = ipoc - 1

      do ia=ipoa,ipoa+illa-1

      call dancd(i1(ia),i2(ia),j)
      cfac = aimag(fun(j))
!      IF(dABS(CFAC).LT.EPS) GOTO 100
!      IF(dABS(CFAC*CC(IA)).LT.EPS) GOTO 100
       if(abs(cfac*cc(ia)).lt.eps.or.abs(cc(ia)).lt.eps) goto 100

      ic = ic + 1
      cc(ic) = cc(ia)*cfac
      i1(ic) = i1(ia)
      i2(ic) = i2(ia)

 100  continue
      enddo

      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
         write(lout,*)'ERROR IN DACFU '
         call dadeb(31,'ERR DACFU ',1)
      endif

      return
      end


subroutine dacfut(ina,fun,inc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : cc,eps,nomax,idalm,idall,i1,i2,lnv
      implicit none
      integer i,ia,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,ipoa,ipoc,j
      real(kind=fPrec) cfac,fun,rr
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE PRECISION FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1


      dimension j(lnv)

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)

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
#ifdef CR
!     call abend('667                                               ')
#else
!        stop 667
#endif
      return
      endif
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)

      ic = ipoc - 1

      do 100 ia=ipoa,ipoa+illa-1

      call dancd(i1(ia),i2(ia),j)
      cfac = fun(j)
!      IF(dABS(CFAC).LT.EPS) GOTO 100
!      IF(dABS(CFAC*CC(IA)).LT.EPS) GOTO 100
       if(abs(cfac*cc(ia)).lt.eps.or.abs(cc(ia)).lt.eps) goto 100

      ic = ic + 1
      cc(ic) = cc(ia)*cfac
      i1(ic) = i1(ia)
      i2(ic) = i2(ia)

 100  continue

      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
         write(lout,*)'ERROR IN DACFU '
         call dadeb(31,'ERR DACFU ',1)
      endif

      return
      end


subroutine dapri(ina,iunit)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : daname,cc,eps,nda,nomax,nvmax,idano,idanv,idapo,idalm,idall,i1,i2,ieo,&
        ia1,ia2,lnv
      implicit none
      integer i,ii,iii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,iunit,j,k
!     ***************************
!       Frank
!     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
!
!-----------------------------------------------------------------------------1

!-----------------------------------------------------------------------------9
      dimension j(lnv)

      if(ina.lt.1.or.ina.gt.nda) then
         write(lout,*)'ERROR IN DAPRI, INA = ',ina
!        X = SQRT(-ONE)
!        PRINT*,X
#ifdef CR
      call abend('                                                  ')
#else
         stop
#endif
      endif

      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)

      write(iunit,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)') daname(ina),', NO =',inoa,', NV =',inva,', INA =',ina, &
     &'***********'//'**********************************'

      iout = 0
      ioa = 0

      if(inva.eq.0) then
         write(iunit,'(A)') '    I  VALUE  '
#ifdef CRLIBM
#ifndef LF95
                                                  call enable_xp()
#endif
#endif
         do i = ipoa,ipoa+illa-1
           write(iunit,'(I6,2X,G21.14)') i-ipoa, cc(i)
           !Eric
           write(111) cc(i)
         end do
#ifdef CRLIBM
#ifndef LF95
                                                  call disable_xp()
#endif
#endif
      elseif(nomax.eq.1) then
         if(illa.ne.0) write(iunit,'(A)') '    I  COEFFICIENT          ORDER   EXPONENTS'
         if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
#ifdef CRLIBM
#ifndef LF95
                                                  call enable_xp()
#endif
#endif
         do 90 i=1,illa
           do k=1,inva
             j(k)=0
           enddo
             iout=iout+1
           if(i.ne.1) then
             j(i-1)=1
             ioa=1
           endif
         write(iunit,'(I6,2X,G21.14,I5,4X,18(2I2,1X))') iout,cc(ipoa+i-1),ioa,(j(iii),iii=1,nvmax)
         write(111) cc(ipoa+i-1)
         write(iunit,'(G21.14)') cc(ipoa+i-1)
         write(111) cc(ipoa+i-1)
 90     continue
#ifdef CRLIBM
#ifndef LF95
                                                  call disable_xp()
#endif
#endif
      else
         if(illa.ne.0) write(iunit,'(A)') '    I  COEFFICIENT          ORDER   EXPONENTS'
         if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
      do ioa = 0,inoa
        do ii=ipoa,ipoa+illa-1
          if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
          call dancd(i1(ii),i2(ii),j)
!ETIENNE
          if(abs(cc(ii)).gt.eps) then
!ETIENNE
          iout = iout+1
#ifdef CRLIBM
#ifndef LF95
                                                  call enable_xp()
#endif
#endif
          write(iunit,'(I6,2X,G21.14,I5,4X,18(2I2,1X))') iout,cc(ii),ioa,(j(iii),iii=1,nvmax)
!Eric
          write(111) cc(ii)
!ETIENNE
          write(iunit,'(G21.14)') cc(ii)
          write(111) cc(ii)
          endif
!ETIENNE
!
#ifdef CRLIBM
#ifndef LF95
                                                  call disable_xp()
#endif
#endif
 100      continue
        end do
      end do

      endif

      write(iunit,'(A)') '                                      '

!Eric
      write(111) zero
      return
      end


subroutine dapri77(ina,iunit)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : daname,cc,eps,nda,nomax,nocut,idano,idanv,idapo,idalm,idall,i1,i2,ieo,&
        ia1,ia2,lnv
      implicit none
      integer i,ii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,iunit,j
      character(len=10) c10,k10
!     ***************************
!       Etienne
!     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
!
!-----------------------------------------------------------------------------1

      dimension j(lnv)

      if(iunit.eq.0) return

      if(ina.lt.1.or.ina.gt.nda) then
         write(lout,*)'ERROR IN DAPRI, INA = ',ina
#ifdef CR
      call abend('                                                  ')
#else
         stop
#endif
      endif

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
      write(iunit,'(/1X,A10,A6,I5,A6,I5,A7,I5/1X,A/)') daname(ina),', NO =',inoa,', NV =',inva,', INA =',ina, &
     &'***********'//'**********************************'

      if(illa.ne.0) write(iunit,'(A)') '    I  COEFFICIENT          ORDER   EXPONENTS'
      if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '

      c10='      NO ='
      k10='      NV ='

      write(iunit,'(A10,I6,A10,I6)') c10,inoa,k10,inva

      iout = 0

!      DO 100 IOA = 0,INOA
       do ioa = 0,nocut
         do ii=ipoa,ipoa+illa-1
           if(nomax.ne.1) then
             if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
           end if
          if(abs(cc(ii)).gt.eps) then
            if(nomax.ne.1) then
              call dancd(i1(ii),i2(ii),j)
              iout = iout+1
            else
              if(ii.eq.ipoa.and.ioa.eq.1) goto 100
              if(ii.gt.ipoa.and.ioa.eq.0) goto 100
              do i=1,lnv
                j(i)=0
              end do
              if(ii.ne.ipoa) j(ii-ipoa)=1
              iout = iout+1
            end if

#ifdef CRLIBM
#ifndef LF95
                                                  call enable_xp()
#endif
#endif
!      WRITE(IUNIT,*) IOA,CC(II),(J(I),I=1,INVA)
      if(abs(cc(ii)).gt.eps) then
      if(eps.gt.1.e-37_fPrec) then
       write(iunit,501) ioa,cc(ii),(j(i),i=1,inva)
      else
       write(iunit,503) ioa,cc(ii),(j(i),i=1,inva)
      endif
      endif
 501  format(' ', i3,1x,g23.16,1x,100(1x,i2))
 503  format(' ', i3,1x,g23.16,1x,100(1x,i2))
 502  format(' ', i5,1x,g23.16,1x,100(1x,i2))
#ifdef CRLIBM
#ifndef LF95
                                                  call disable_xp()
#endif
#endif
      end if
!ETIENNE
 100      continue
        end do
      end do

      do i=1,lnv
        j(i)=0
      end do

      if(iout.eq.0) iout=1
#ifdef CRLIBM
#ifndef LF95
                                                  call enable_xp()
#endif
#endif

      write(iunit,502) -iout,zero,(j(i),i=1,inva)

#ifdef CRLIBM
#ifndef LF95
                                                  call disable_xp()
#endif
#endif
      return
      end


subroutine dashift(ina,inc,ishift)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : daname,cc,eps,nda,nomax,nocut,idano,idanv,idapo,idalm,idall,i1,i2,ieo,&
        ia1,ia2,lnv
      implicit none
      integer i,ii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,j

!-----------------------------------------------------------------------------9
      dimension j(lnv)
!
      integer inb(1),ishift,ich,ik,jd(lnv),inc
!-----------------------------------------------------------------------------3

      inb(1)=0
      if(ina.lt.1.or.ina.gt.nda) then
         write(lout,*)'ERROR IN DAPRI, INA = ',ina
#ifdef CR
      call abend('                                                  ')
#else
         stop
#endif
      endif

      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
      call daall(inb(1),1,'$$DAJUNK$$',inoa,inva)

!      WRITE(IUNIT,*) INA, ' in dapri ', DANAME(INA)
!      WRITE(6,*) INA, ' in dapri ', DANAME(INA)
! 611  WRITE(6,*) ' MORE '
!        READ(5,*) MORE
!        IF(MORE.GT.0) THEN
!        WRITE(6,*) MORE,' ',DANAME(MORE)
!        GOTO 611
!        ENDIF
      iout = 0

!      DO 100 IOA = 0,INOA
       do ioa = 0,nocut
       do ii=ipoa,ipoa+illa-1
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

!      WRITE(IUNIT,*) IOA,CC(II),(J(I),I=1,INVA)
      if(abs(cc(ii)).gt.eps) then
      if(eps.gt.1.e-37_fPrec) then
#ifdef CRLIBM
!                                                 call enable_xp()
#endif
!       write(iunit,501) ioa,cc(ii),(j(i),i=1,inva)
#ifdef CRLIBM
!                                                 call disable_xp()
#endif
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
       call dapok(inb(1),jd,cc(ii))
      else
#ifdef CRLIBM
!                                                 call enable_xp()
#endif
!       write(iunit,503) ioa,cc(ii),(j(i),i=1,inva)
#ifdef CRLIBM
!                                                 call disable_xp()
#endif
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
        call dapok(inb(1),jd,cc(ii))
      endif
      endif
 501  format(' ', i3,1x,g23.16,1x,100(1x,i2))
 503  format(' ', i3,1x,g23.16,1x,100(1x,i2))
 502  format(' ', i5,1x,g23.16,1x,100(1x,i2))

      endif
!ETIENNE

 100      continue
        end do
      end do

      do i=1,lnv
        j(i)=0
      end do

      call dacop(inb(1),inc)
      call dadal(inb(1),1)

      return
      end


subroutine darea(ina,iunit)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : daname,cc,nda,nomax,idano,idanv,idapo,idalm,idall,ia1,ia2,lnv
      implicit none
      integer i,ic,iche,ii,ii1,ii2,iin,illa,ilma,ina,inoa,inva,io,io1,ipoa,iunit,iwarin,iwarno,iwarnv,j,nno
      real(kind=fPrec) c
!       Frank

!-----------------------------------------------------------------------------9

      character(len=10) c10
      dimension j(lnv)

      if(ina.lt.1.or.ina.gt.nda) then
         write(lout,*)'ERROR IN DAREA, INA = ',ina
!        X = SQRT(-ONE)
!        PRINT*,X
#ifdef CR
      call abend('                                                  ')
#else
         stop
#endif
      endif

      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)

      do i=1,lnv
        j(i) = 0
      end do

      call daclr(1)

      ic = 0

      iwarno = 0
      iwarnv = 0
      iwarin = 0

      read(iunit,'(A10)') c10
      read(iunit,'(18X,I4)') nno
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10

      iin = 0

  10  continue
      iin = iin + 1
#ifdef CRLIBM
      call enable_xp()
#endif
      read(iunit,'(I6,2X,G21.14,I5,4X,18(2I2,1X))') ii,c,io,(j(i),i=1,inva)
#ifdef CRLIBM
      call disable_xp()
#endif
!Eric
      read(111) c
!
      if(ii.eq.0) goto 20
!ETIENNE
#ifdef CRLIBM
      call enable_xp()
#endif
!Eric
      read(iunit,'(G21.14)') c
!Eric
      read(111) c
#ifdef CRLIBM
      call disable_xp()
#endif
!ETIENNE
      if(ii.ne.iin) then
         iwarin = 1
      endif
      io1 = 0
      do i=1,inva
        io1 = io1 + j(i)
      end do

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

      if(nomax.ne.1) then
        ic = ic + 1
        call dadcd(j,ii1,ii2)
        ic = ia1(ii1) + ia2(ii2)
        cc(ic) = c
        goto 10
      else
        iche=0
        do i=1,inva
          if(j(i).eq.1) iche=i
        enddo
        cc(ipoa+iche)=c
        goto 10
      endif

  20  continue

      if(nomax.ne.1) call dapac(ina)

      return
      end
!FF


subroutine darea77(ina,iunit)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : daname,cc,nda,nomax,idano,idanv,idapo,idalm,idall,ia1,ia2,lnv
      implicit none
      integer i,ic,iche,ii,ii1,ii2,iin,illa,ilma,ina,inoa,inva,ipoa,iunit,j,k,nojoh,nvjoh
      real(kind=fPrec) c
!     ***************************
!     Etienne
!     THIS SUBROUTINE READS THE DA VECTOR INA FROM UNIT IUNIT.
!
!-----------------------------------------------------------------------------1


      character(len=10) c10,k10
      dimension j(lnv)

      if(ina.lt.1.or.ina.gt.nda) then
         write(lout,*)'ERROR IN DAPRI, INA = ',ina
#ifdef CR
      call abend('                                                  ')
#else
         stop
#endif
      endif

      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)

      do i=1,lnv
        j(i) = 0
      end do

      call daclr(1)

      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10,I6,A10,I6)') c10,nojoh,k10,nvjoh

      iin = 0

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

  20  continue

      if(nomax.ne.1) call dapac(ina)

      return
      end


subroutine dadeb(iunit,c,istop)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : daname
      implicit none
      integer istop,iunit
!     *******************************
!
!     THIS SUBROUTINE SERVES AS A DEBUGGING TOOL. IT PRINTS ALL
!     NONZERO INFORMATION IN THE COMMON BLOCKS AND ALL DA  VECTORS.
!
!-----------------------------------------------------------------------------1


      character(len=10) c

!etienne
      write(lout,"(a)") "DABNEW> ERROR "//c(5:10)
#ifdef CR
      call abend('                                                  ')
#else
      call prror(-1)
#endif
end subroutine dadeb


subroutine danum(no,nv,numda)
      use floatPrecision
      implicit none
      integer i,mm,no,numda,nv
!     *****************************
!
!     THIS SUBROUTINE COMPUTES THE NUMBER OF MONOMIALS OF
!     ORDER NO AND NUMBER OF VARIABLES NV

      numda = 1
      mm = max(nv,no)

      do i=1,min(nv,no)
        numda = (numda*(mm+i))/i
      end do

      return
      end


subroutine dainf(inc,inoc,invc,ipoc,ilmc,illc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : nda,idano,idanv,idapo,idalm,idall
      implicit none
      integer illc,ilmc,inc,inoc,invc,ipoc
!     **********************************************
!
!     THIS SUBROUTINE SEARCHES THE NUMBER OF DA VECTOR C
!     AND RETURS THE INFORMATION IN COMMON DA
!
!-----------------------------------------------------------------------------1


      if(inc.ge.1.and.inc.le.nda) then
         inoc = idano(inc)
         invc = idanv(inc)
         ipoc = idapo(inc)
         ilmc = idalm(inc)
         illc = idall(inc)
         return
      endif

      write(lout,*) 'ERROR IN DAINF, DA VECTOR ',inc,' NOT FOUND '
      call dadeb(31,'ERR DAINF ',1)

      return
      end


subroutine dapac(inc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : cc,eps,nmmax,idalm,idall,i1,i2,ie1,ie2
      implicit none
      integer i,ic,illc,ilmc,inc,inoc,invc,ipoc
      real(kind=fPrec) ccc
!     ************************
!
!     THIS SUBROUTINE PACKS THE INFORMATION IN THE SCRATCH VECTOR 1
!     INTO THE VECTOR INC. IF LF = 1, THE FILTERING (CF DAMUF) IS
!     PERFORMED.
!     INVERSE IS DAUNP.
!
!-----------------------------------------------------------------------------1


      call dainf(inc,inoc,invc,ipoc,ilmc,illc)

      ic = ipoc - 1

         do 100 i=1,nmmax
         ccc = cc(i)
         if(abs(ccc).lt.eps) goto 100
         ic = ic + 1
         cc(ic) = ccc
         i1(ic) = ie1(i)
         i2(ic) = ie2(i)
 100     continue

      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
         write(lout,*)'ERROR IN DAPAC '
         call dadeb(31,'ERR DAPAC ',1)
      endif

      return
      end
!
!

subroutine dachk(ina,inoa,inva, inb,inob,invb, inc,inoc,invc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      implicit none
      integer ierr,ina,inb,inc,inoa,inob,inoc,inva,invb,invc,invsum,lsw
!     *************************************************************
!
!     THIS SUBROUTINE CHECKS IF THE VECTORS A, B AND C
!     HAVE COMPATIBLE ATTRIBUTES

      parameter(lsw=1)

      if(lsw.eq.1) return

      ierr = 0

!     CASE OF A UNARY OPERATION
!     *************************

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
            write(lout,*)'ERROR IN DACHK, ',ina,' AND ',inc,' ARE INCOMPATIBLE',inoa,inva,inoc,invc
            call dadeb(31,'ERR DACHK1',1)
         endif

!     CASE OF A BINARY OPERATION
!     **************************

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
            write(lout,*)'ERROR IN DACHK, ',ina,',',inb,' AND ',inc,' ARE INCOMPATIBLE'
            call dadeb(31,'ERR DACHK2',1)
         endif
      endif

      return
      end


subroutine damch(iaa,ia)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      implicit none
      integer i,ia,iaa,illa,ilma,ino1,inoi,inv1,invi,ipoa
!     ************************

!     THIS SUBROUTINE CHECKS IF THE IA VECTORS IN THE MATRIX IA HAVE
!     IDENTICAL ATTRIBUTES.

      dimension iaa(*)

      call dainf(iaa(1),ino1,inv1,ipoa,ilma,illa)

      do 10 i=2,ia
      call dainf(iaa(i),inoi,invi,ipoa,ilma,illa)
      if(ino1.ne.inoi.or.inv1.ne.invi) then
         write(lout,*)'ERROR IN DAMCH, VECTORS ',iaa(1),' AND ',iaa(i),' ARE INCOMPATIBLE '
#ifdef CR
      call abend('                                                  ')
#else
         stop
#endif
      endif
  10  continue

      return
      end


subroutine dadcd(jj,ic1,ic2)
      use floatPrecision
      use mod_lie_dab, only : nomax,nvmax,lnv
      implicit none
      integer i,ibase,ic1,ic2,isplit
!     ****************************
!
!     THIS SUBROUTINE CODES THE EXPONENTS IN JJ INTO THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1


      integer jj(lnv)
      ibase = nomax + 1
      isplit = (nvmax+1)/2
      ic1 = 0
      ic2 = 0

      do i=nvmax,isplit+1,-1
        ic2 = ic2*ibase + jj(i)
      end do

      do i=isplit,1,-1
        ic1 = ic1*ibase + jj(i)
      end do

      return
      end


subroutine dancd(ic1,ic2,jj)
      use floatPrecision
      use mod_lie_dab, only : epsmac,nomax,nvmax,lnv
      implicit none
      integer i,ibase,ic,ic1,ic2,isplit
      real(kind=fPrec) x
!     ****************************
!
!     THIS SUBROUTINE ENCODES THE EXPONENTS IN JJ FROM THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1


      integer jj(*)
      ibase = nomax + 1
      isplit = (nvmax+1)/2

      ic = ic1
      do i=1,isplit
        x  = real(ic,fPrec)/real(ibase,fPrec)
        ic = int(x+epsmac)
        jj(i) = nint(ibase*(x-ic))
      end do

      ic = ic2

      do i=isplit+1,nvmax
        x  = real(ic,fPrec)/real(ibase,fPrec)
        ic = int(x+epsmac)
        jj(i) = nint(ibase*(x-ic))
      end do

      do i=nvmax+1,lnv
        jj(i) = 0
      end do
      return
      end

!ETIENNE

subroutine datra(idif,ina,inc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : cc,epsmac,nomax,nvmax,idalm,idall,i1,i2
      implicit none
      integer i,ibase,ic,ider1,ider1s,ider2,ider2s,idif,iee,ifac,illa,  &
     &illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,ipoa,ipoc,jj
      real(kind=fPrec) x,xdivi
!     ******************************
!
!     THIS SUBROUTINE COMPUTES THE PSEUDO DERIVATIVE WITH RESPECT TO VARIABLE I
!     OF THE VECTOR A AND STORES THE RESULT IN C.
!
!     dx^n/dx= x^(n-1)
!-----------------------------------------------------------------------------1


      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)

!       CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)

      if(nomax.eq.1) then
        call dader(idif,ina,inc)
        return
      endif
      ibase = nomax + 1

      if(idif.gt.(nvmax+1)/2) then
         ider1  = 0
         ider1s = 0
         ider2  = idif-(nvmax+1)/2
         ider2s = 1
         do jj=1,ider2-1
           ider2s = ider2s*ibase
         end do
         xdivi  = ider2s*ibase
      else
         ider1  = idif
         ider1s = 1
         do jj=1,ider1-1
           ider1s = ider1s*ibase
         end do
         ider2  = 0
         ider2s = 0
         xdivi  = ider1s*ibase
      endif

      ibase = nomax+1

      ic = ipoc-1

      do 100 i=ipoa,ipoa+illa-1

      if(ider1.eq.0) then
        iee = i2(i)
      else
        iee = i1(i)
      endif

      x = iee/xdivi
      ifac = int(ibase*(x-int(x+epsmac)+epsmac))

      if(ifac.eq.0) goto 100

!etienne      IFAC = INT(IBASE*(X-INT(X)+1.D-8))

!etienne      IF(IFAC.EQ.0) GOTO 100

      ic = ic + 1
      cc(ic) = cc(i)
      i1(ic) = i1(i) - ider1s
      i2(ic) = i2(i) - ider2s

 100  continue

      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
         write(lout,*)'ERROR IN DADTRA'
         call dadeb(111,'ERR DADTRA',1)
      endif

      return
      end


subroutine etred(no1,nv1,ic1,ic2,no2,nv2,i11,i21)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : lnv
      implicit none
      integer i,i11,i21,ic,ic1,ic2,no1,no2,nv1,nv2
!     ****************************
!
!     THIS SUBROUTINE CODES THE EXPONENTS IN JJ INTO THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1

!
      integer jj(lnv)

      if(nv1.gt.lnv.or.nv2.gt.lnv) then
      write(lout,*) ' ERROR IN RECODING '
#ifdef CR
      call abend('123                                               ')
#else
      stop 123
#endif
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


subroutine hash(no1,nv1,jj,ic1,ic2)
      use floatPrecision
      implicit none
      integer i,ibase,ic1,ic2,isplit,no1,nv1
!     ****************************
!
!     THIS SUBROUTINE CODES THE EXPONENTS IN JJ INTO THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1

      integer jj(*)

      ibase = no1 + 1
      isplit = (nv1+1)/2
      ic1 = 0
      ic2 = 0

      do i=nv1,isplit+1,-1
        ic2 = ic2*ibase + jj(i)
      end do

      do i=isplit,1,-1
        ic1 = ic1*ibase + jj(i)
      end do

      return
      end


subroutine dehash(no1,nv1,ic1,ic2,jj)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : epsmac
      implicit none
      integer i,ibase,ic,ic1,ic2,isplit,no1,nv1
      real(kind=fPrec) x
!     ****************************
!
!     THIS SUBROUTINE ENCODES THE EXPONENTS IN JJ FROM THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1


      integer jj(*)

      epsmac=c1m7
      ibase = no1 + 1
      isplit = (nv1+1)/2

      ic = ic1
      do i=1,isplit
        x  = real(ic,fPrec)/real(ibase,fPrec)
        ic = int(x+epsmac)
        jj(i) = nint(ibase*(x-ic))
      end do

      ic = ic2
      do i=isplit+1,nv1
        x  = real(ic,fPrec)/real(ibase,fPrec)
        ic = int(x+epsmac)
        jj(i) = nint(ibase*(x-ic))
      end do

      return
      end


subroutine daswap(j1,j2,inb)
      use floatPrecision
      use mod_lie_dab, only : cc,nomax,nvmax,i1,i2,ia1,ia2,lnv
      implicit none
      integer ia,ic,ic1,ic2,illb,ilmb,inb,inob,invb,ipob,j1,j2,jj,k1,k2
!     *************************
!
!     SWAP A DA VECTOR
!
!-----------------------------------------------------------------------------1


      dimension jj(lnv)

      call dainf(inb,inob,invb,ipob,ilmb,illb)
      call daclr(1)

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

 100  continue

      call dapac(inb)
      return
      end


subroutine dagauss(ina,inb,nd2,anorm)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : cc,facint,i1,i2,lnv
      implicit none
      integer i,ia,ib,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,ipoa,ipob,ja,jb,nd2
      real(kind=fPrec) anorm,gau
!     ***************************
!
!     THIS SUBROUTINE COMPUTES THE NORM OF THE DA VECTOR A
!
!-----------------------------------------------------------------------------1


      dimension ja(lnv),jb(lnv)

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)

      anorm = zero

      do ia=ipoa,ipoa+illa-1
        do ib=ipob,ipob+illb-1
          call dancd(i1(ia),i2(ia),ja)
          call dancd(i1(ib),i2(ib),jb)
          gau=one
          do i=1,nd2
            gau= facint(ja(i)+jb(i))*gau
          end do
          anorm = anorm + cc(ia)*cc(ib)*gau
        end do
      end do

      return
      end


subroutine daran(ina,cm,xran)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : cc,nomax,nvmax,nmmax,idalm,idall
      implicit none
      integer i,illa,ilma,ina,inoa,inva,ipoa
      real(kind=fPrec) bran,cm,xran
!     ************************
!
!     THIS SUBROUTINE FILLS THE DA VECTOR A WITH RANDOM ENTRIES.
!     FOR CM > 0, THE VECTOR IS FILLED WITH REALS,
!     FOR CM < 0, THE VECTOR IS FILLED WITH SINGLE DIGIT INTEGERS
!     ABS(CM) IS THE FILLING FACTOR
!
!-----------------------------------------------------------------------------1


      call dainf(ina,inoa,inva,ipoa,ilma,illa)

      if(inva.eq.0.or.nomax.eq.1) then
         do 10 i=ipoa,ipoa+ilma-1
         if(cm.gt.zero) then
            cc(i) = bran(xran)
            if(cc(i).gt.cm) cc(i) = zero
         elseif(cm.lt.zero) then
            cc(i) = int(1+10*bran(xran))
            if(cc(i).gt.-c1e1*cm) cc(i) = zero
         endif
  10     continue
         idall(ina) = idalm(ina)
         return
      endif

      if(inoa.ne.nomax.or.inva.ne.nvmax) then
         write(lout,*)'ERROR IN DARAN, ONLY VECTORS WITH NO = NOMAX AND'//' NV = NVMAX ALLOWED'
         call dadeb(31,'ERR DARAN1',1)
      endif

      call daclr(1)

      do 100 i=1,nmmax
      if(cm.gt.zero) then
         cc(i) = bran(xran)
         if(cc(i).gt.cm) cc(i) = zero
      elseif(cm.lt.zero) then
         cc(i) = int(1+10*bran(xran))
         if(cc(i).gt.-c1e1*cm) cc(i) = zero
      else
         write(lout,*)'ERROR IN ROUTINE DARAN'
         call dadeb(31,'ERR DARAN2',1)
      endif
 100  continue

      call dapac(ina)

      return
      end
!

      real(kind=fPrec) function bran(xran)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      implicit none
      real(kind=fPrec) xran
!     ************************************
!
!     VERY SIMPLE RANDOM NUMBER GENERATOR
!
      xran = xran + c1e1
      if(xran.gt.c1e4) xran = xran - 9999.12345_fPrec
      bran = abs(sin_mb(xran))
      bran = c1e1*bran
      bran = bran - int(bran)
      return
      end


subroutine danorm2(ina,inc)
      use floatPrecision
      implicit none
      integer illc,ilmc,ina,inc,inoc,invc,ipoc
      integer incc(1)
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1

!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc(1)=0
        call daall(incc(1),1,'$$DAJUNK$$',inoc,invc)
        call danorm2t(ina,incc(1))
        call dacop(incc(1),inc)
        call dadal(incc(1),1)
      else
        call danorm2t(ina,inc)
      endif

      return
      end


subroutine danorm2t(ina,inb)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : cc,nocut,idalm,idall,i1,i2,ieo,ia1,ia2
      implicit none
      integer ia,ib,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,    &
     &ipoa,ipob
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1

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
         write(lout,*)'ERROR IN DANORM'
         call dadeb(31,'ERR DANOR1',1)
      endif
!
      return
      end


subroutine danormr(ina,inc)
      use floatPrecision
      implicit none
      integer illc,ilmc,ina,inc,inoc,invc,ipoc
      integer incc(1)
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1

!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc(1)=0
        call daall(incc(1),1,'$$DAJUNK$$',inoc,invc)
        call danormrt(ina,incc(1))
        call dacop(incc(1),inc)
        call dadal(incc(1),1)
      else
        call danormrt(ina,inc)
      endif

      return
      end


subroutine danormrt(ina,inb)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : cc,nocut,idalm,idall,i1,i2,ieo,ia1,ia2
      implicit none
      integer ia,ib,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,    &
     &ipoa,ipob
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1

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
      cc(ib) = sqrt(cc(ia))
      i1(ib) = i1(ia)
      i2(ib) = i2(ia)
!
 100  continue
!
      idall(inb) = ib-ipob+1
      if(idall(inb).gt.idalm(inb)) then
         write(lout,*)'ERROR IN DANORM '
         call dadeb(31,'ERR DANOR2',1)
      endif
!
      return
      end

subroutine dakey(c)
      use floatPrecision
      implicit none
      character c*(*)
!

      return
!
      end
! ANFANG UNTERPROGRAMM

subroutine dapri6(ina,result,ien,i56)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : daname,cc,nda,nomax,idano,idanv,idapo,idalm,idall,i1,i2,ieo,ia1,ia2,lnv
      implicit none
      integer i,i56,ien,ihp,ii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,j
      real(kind=fPrec) result
!     *************************************
!
!     THIS SUBROUTINE IS FOR REDUCED STORAGE DA VERSION JULY 91
!     RESULT CONTAINS THE (IEN-1) TH DERIVATIVE TO ENERGY AFTER EXECUTION
!     I56 SAYS WHETHER THE 5TH OR THE 6TH COORDINATE IS THE ENERGY
!     AND MUST HAVE THE VALUE 5 OR 6 ACCORDINGLY
!-----------------------------------------------------------------------------1

      dimension j(lnv)
      result=0.
      if(ina.lt.1.or.ina.gt.nda) then
        write(lout,*)'ERROR IN DAPRI6, INA = ',ina
#ifdef CR
      call abend('                                                  ')
#else
        stop
#endif
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
        do ioa=0,inoa
        do ii=ipoa,ipoa+illa-1
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
      end do
      end do
      end if
      return
      end
! ANFANG UNTERPROGRAMM


subroutine darea6(ina,zfeld,i56)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : daname,cc,nda,nomax,idano,idanv,idapo,idalm,idall,ia1,ia2,lnv
      implicit none
      integer i,i56,ic,ii1,ii2,iin,illa,ilma,ina,inoa,inva,io,io1,ip,ipoa,iwarin,iwarno,iwarnv,j
      real(kind=fPrec) zfeld
!     *************************************
!
!     THIS SUBROUTINE IS FOR REDUCED STORAGE DA VERSION JULY 91
!     I56 SAYS WHETHER THE 5TH OR THE 6TH COORDINATE IS THE ENERGY
!     AND MUST HAVE THE VALUE 5 OR 6 ACCORDINGLY
!
!-----------------------------------------------------------------------------1

      dimension zfeld(100)
!-----------------------------------------------------------------------------3
      dimension j(lnv)
      if(ina.lt.1.or.ina.gt.nda) then
        write(lout,*)'ERROR IN DAREA6, INA = ',ina
#ifdef CR
      call abend('                                                  ')
#else
        stop
#endif
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
        if(iwarnv.eq.0) write(lout,*)'WARNING IN DAREA6, FILE ', 'CONTAINS MORE VARIABLES THAN VECTOR'
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

      real(kind=fPrec) function dare(ina)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : cc,nomax,idano,idanv,idapo,idalm,idall,i1,i2,ieo,ia1,ia2,lnv
      implicit none
      integer ii,illa,ilma,ina,inoa,inva,ioa,ipoa,j,jj
!     ***********************************
!     NEW VERSION OF DARE, AUGUST 1992
!     SUPPOSED TO TREAT THE 0TH COMPONENT ACCURATELY
!
!     30.10 1997 E.Mcintosh & F.Schmidt
!
!-----------------------------------------------------------------------------1

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
      do ii=ipoa,ipoa+illa-1
        if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
        call dancd(i1(ii),i2(ii),j)
        do jj=1,inva
          if(j(jj).ne.0) goto 100
        end do
        dare = cc(ipoa)
        return
 100    continue
      end do

      dare = zero
      return
      end
! ANFANG UNTERPROGRAMM

subroutine daprimax(ina,iunit)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : daname,cc,eps,nda,nomax,nvmax,idano,idanv,idapo,idalm,idall,i1,i2,ieo,&
        ia1,ia2,lnv
      implicit none
      integer i,ii,iii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,iunit,j
!     ***************************
!
!     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
!
!-----------------------------------------------------------------------------1

!-----------------------------------------------------------------------------9

      dimension j(lnv)

      if(ina.lt.1.or.ina.gt.nda) then
         write(lout,*)'ERROR IN DAPRI, INA = ',ina
#ifdef CR
      call abend('                                                  ')
#else
         stop
#endif
      endif

      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)

      iout = 0
      ioa = 0

      if(inva.eq.0) then
#ifdef CRLIBM
#ifndef LF95
                                                  call enable_xp()
#endif
#endif
         do i = ipoa,ipoa+illa-1
           write(iunit,'(I6,2X,G21.14)') i-ipoa, cc(i)
         end do
#ifdef CRLIBM
#ifndef LF95
                                                  call disable_xp()
#endif
#endif
      elseif(nomax.eq.1) then
         do 90 i=1,illa
             iout=iout+1
           if(i.ne.1) then
             j(i-1)=1
             ioa=1
           endif
#ifdef CRLIBM
#ifndef LF95
                                                  call enable_xp()
#endif
#endif
         write(iunit,'(I6,2X,G21.14,I5,4X,18(2I2,1X))')                 &
     &iout,cc(ipoa+i-1),ioa,(j(iii),iii=1,nvmax)
         write(iunit,*) cc(ipoa+i-1)
#ifdef CRLIBM
#ifndef LF95
                                                  call disable_xp()
#endif
#endif
 90      continue
      else
        iout = 0
        do ioa = 0,inoa
        do ii=ipoa,ipoa+illa-1
          if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
          call dancd(i1(ii),i2(ii),j)
!ETIENNE
          if(abs(cc(ii)).gt.eps) then
!ETIENNE
          iout = iout+1
#ifdef CRLIBM
#ifndef LF95
                                                  call enable_xp()
#endif
#endif
          write(iunit,'(I6,2X,G21.14,I5,4X,18(2I2,1X))')                &
     &iout,cc(ii),ioa,(j(iii),iii=1,nvmax)
!ETIENNE
          write(iunit,*) cc(ii)
#ifdef CRLIBM
#ifndef LF95
                                                  call disable_xp()
#endif
#endif
          endif
!ETIENNE
 100       continue
        end do
      end do

      endif

      return
      end

!  unknown stuff

subroutine damono(ina,jd,cfac,istart,inc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : cc,idalm,idall,i1,i2
      implicit none
      integer ia,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,ipoa,ipoc,istart,jd
      real(kind=fPrec) cfac
!     *****************************
!
!     THIS SUBROUTINE RETURNS THE MONOMIALS ONE BY ONE
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1


      dimension jd(*)

      if(ina.eq.inc) then
       write(lout,*) ' USE DIFFERENT POWER SERIES IN DAMONO '
#ifdef CR
      call abend('999                                               ')
#else
       stop 999
#endif
      endif
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
      if(istart.eq.0) then
        istart=illa
        return
      endif

!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)

      ic = ipoc - 1

      ia=ipoa+istart-1

      call dancd(i1(ia),i2(ia),jd)

      ic = ic + 1
      cc(ic) = cc(ia)
       cfac=cc(ia)
      i1(ic) = i1(ia)
      i2(ic) = i2(ia)
 100  continue


      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
         write(lout,*)'ERROR IN DAMONO'
         call dadeb(31,'ERR DAMONO',1)
      endif

      return
      end


subroutine dacycle(ina,ipresent,value,j,illa)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : daname,cc,nda,idano,idanv,idapo,idalm,idall,i1,i2,lnv
      implicit none
      integer i,ii,illa,ilma,ina,inoa,inva,iout,ipoa,ipresent,j
      real(kind=fPrec) value
!     ***************************
!
!     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
!
!-----------------------------------------------------------------------------1

      dimension j(lnv)
!
!
      if(ina.lt.1.or.ina.gt.nda) then
        write(lout,*)'ERROR IN DAPRI, INA = ',ina
#ifdef CR
      call abend('                                                  ')
#else
        stop
#endif
      endif

      if(ina.eq.0) then
        value=zero
        illa=0

        do i=1,lnv
          j(i)=0
        end do

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

subroutine daorder(ina,iunit,jx,invo,nchop)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : daname,cc,nda,idano,idanv,idapo,idalm,idall,ia1,ia2,lnv
      implicit none
      integer i,ic,ii,ii1,ii2,iin,illa,ilma,ina,inoa,inva,invo,io,io1,ipoa,iunit,iwarin,iwarno,     &
        iwarnv,j,jh,jt,jx,nchop
      real(kind=fPrec) c
!     ***************************
!
!     THIS SUBROUTINE READS THE DA VECTOR INA FROM UNIT IUNIT.
!
!-----------------------------------------------------------------------------1


      character(len=10) c10
      dimension j(lnv),jx(lnv),jt(lnv)

      if(ina.lt.1.or.ina.gt.nda) then
         write(lout,*)'ERROR IN DAPRI, INA = ',ina
#ifdef CR
      call abend('                                                  ')
#else
         stop
#endif
      endif

      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)

      do i=1,lnv
        jt(i)=0
        j(i) = 0
      end do

      call daclr(1)

      ic = 0

      iwarno = 0
      iwarnv = 0
      iwarin = 0

      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10

      iin = 0

  10  continue
      iin = iin + 1
#ifdef CRLIBM
#ifndef LF95
                                                  call enable_xp()
#endif
#endif
      read(iunit,'(I6,2X,G21.14,I5,4X,18(2I2,1X))') ii,c,io,(jt(i),i=1,invo)
#ifdef CRLIBM
#ifndef LF95
                                                  call disable_xp()
      read(111) c
#endif
#endif

      if(ii.eq.0) goto 20
!etienne
!Eric
!     read(iunit,*) c
#ifdef CRLIBM
#ifndef LF95
                                                  call enable_xp()
#endif
#endif
      read(iunit,'(G21.14)') c
#ifdef CRLIBM
#ifndef LF95
                                                  call disable_xp()
#endif
#endif
      read(111) c

      do jh=1,invo
        j(jh)=jt(jx(jh))
      end do

      do jh=nchop+1,inva
        j(jh)=0
      end do
!etienne
      io1 = 0

      do i=1,inva
        io1 = io1 + j(i)
      end do

      ic = ic + 1
      call dadcd(j,ii1,ii2)
      ic = ia1(ii1) + ia2(ii2)
      cc(ic) = c
      goto 10

  20  continue

      call dapac(ina)

      return
      end

!ETIENNE

subroutine datrash(idif,ina,inc)
      use floatPrecision
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : cc,nomax,nvmax,idalm,idall,i1,i2,lnv
      implicit none
      integer i,ibase,ic,ider1,ider1s,ider2,ider2s,idif,ikil1,ikil2,illa,illc,ilma,ilmc,ina,inc,    &
        inoa,inoc,inva,invc,ipoa,ipoc,jj
      real(kind=fPrec) xdivi
!     ******************************
!
!     THIS SUBROUTINE COMPUTES THE DERIVATIVE WITH RESPECT TO VARIABLE I
!     OF THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1


      integer jx(lnv)

!      call daclr(1)
!      call dacop(ina,1)
!      call dapac(ina)

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)

!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)

      ibase = nomax + 1

      if(idif.gt.(nvmax+1)/2) then
         ider1  = 0
         ider1s = 0
         ider2  = idif-(nvmax+1)/2
         ider2s = 1
         do jj=1,ider2-1
           ider2s = ider2s*ibase
         end do
         xdivi  = ider2s*ibase
      else
         ider1  = idif
         ider1s = 1
         do jj=1,ider1-1
           ider1s = ider1s*ibase
         end do
         ider2  = 0
         ider2s = 0
         xdivi  = ider1s*ibase
      endif

      ibase = nomax+1

      ic = ipoc-1

      do 100 i=ipoa,ipoa+illa-1

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
         write(lout,*)'ERROR IN DATRASH '
         call dadeb(111,'ERR DATRAS',1)
      endif
!
      return
      end
