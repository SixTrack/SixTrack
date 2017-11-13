! ANFANG - HAUPTPROGRAMM -
  
      PROGRAM FOXY
!     ************
!
!     THIS IS THE LANGUAGE PROCESSOR FOR THE FOX FORTRAN EXTENSION
!
!     WUENSCHE
!     -ALLE OUTPUTS VON ZENTRALER ROUTINE, DIE RICHTIG FORMATIERT

      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT, ERROR_UNIT

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      !For command line argument parsing
      INTEGER :: cmdarg_i, cmdarg_length, cmdarg_status
      CHARACTER(len=100) :: cmdarg_arg
      CHARACTER(len=100) :: fname_in
      CHARACTER(len=100) :: fname_out

      WRITE(OUTPUT_UNIT,*)' '
      WRITE(OUTPUT_UNIT,*)' '
      WRITE(OUTPUT_UNIT,*)' '
      WRITE(OUTPUT_UNIT,*)' '
      WRITE(OUTPUT_UNIT,*)' '
      WRITE(OUTPUT_UNIT,*)'               **************************************'
      WRITE(OUTPUT_UNIT,*)'               *                                    *'
      WRITE(OUTPUT_UNIT,*)'               *               F O X                *'
      WRITE(OUTPUT_UNIT,*)'               *                                    *'
      WRITE(OUTPUT_UNIT,*)'               *    EXTENDED FORTRAN PRECOMPILER    *'
      WRITE(OUTPUT_UNIT,*)'               *                                    *'
      WRITE(OUTPUT_UNIT,*)'               *             VERSION  2             *'
      WRITE(OUTPUT_UNIT,*)'               *                                    *'
      WRITE(OUTPUT_UNIT,*)'               *            M. BERZ 1989            *'
      WRITE(OUTPUT_UNIT,*)'               *                                    *'
      WRITE(OUTPUT_UNIT,*)'               *     Minor Modifications for        *'
      WRITE(OUTPUT_UNIT,*)'               *             SIXTRACK               *'
      WRITE(OUTPUT_UNIT,*)'               *                                    *'
      WRITE(OUTPUT_UNIT,*)'               *         F. Schmidt 1997            *'
      WRITE(OUTPUT_UNIT,*)'               *                                    *'
      WRITE(OUTPUT_UNIT,*)'               **************************************'
      WRITE(OUTPUT_UNIT,*)' '
      WRITE(OUTPUT_UNIT,*)' '
      WRITE(OUTPUT_UNIT,*)' '
      WRITE(OUTPUT_UNIT,*)' '
      WRITE(OUTPUT_UNIT,*)' '

      ! Read command line arguments -- inspired by read90
      cmdarg_i = 0
      do
         call get_command_argument(cmdarg_i, cmdarg_arg,cmdarg_length,cmdarg_status)
         ! write (6,*) cmdarg_i, cmdarg_arg
         if (len_trim(cmdarg_arg)==0) then
          EXIT !Finished
         end if
         if (cmdarg_i.eq.0) then
            ! Skip first argument (command name)
            continue
         else if (cmdarg_i.eq.1) then
            if (cmdarg_status.ne.0) then
               WRITE(ERROR_UNIT,*) "Error: Input file name too long."
               stop 1
            end if
            fname_in = trim(cmdarg_arg)
         else if (cmdarg_i.eq.2) then
            if (cmdarg_status.ne.0) then
               WRITE(ERROR_UNIT,*) "Error: Output file name too long."
               stop 1
            end if
            fname_out = trim(cmdarg_arg)
         else
            WRITE(ERROR_UNIT,*) "Error: Expected either no arguments OR "
            WRITE(ERROR_UNIT,*) "two arguements (input and output file name)"
            stop 1
         end if
         cmdarg_i = cmdarg_i+1
      end do

      if (cmdarg_i.eq.1) then
         fname_in = "fort.1"
         fname_out = "fort.2"
      else if(cmdarg_i.ne.3) then
         WRITE(ERROR_UNIT,*) "Error: Expected either no arguments OR "
         WRITE(ERROR_UNIT,*) "two arguements (input and output file name)"
         stop 2
      end if
      if (fname_in.eq.fname_out) then
         WRITE(ERROR_UNIT,*) "Error: input and output filename should not be the same!"
         stop 3
      end if

      open(1,file=fname_in,form='formatted',status='old')
      open(2,file=fname_out,form='formatted',status='unknown') !status unknown -> allow overwrites

      CALL predata

      CALL PRECOM

      close(1)
      close(2)
      
      STOP
      END

! ANFANG UNTERPROGRAMM
  
      subroutine predata

! Eric: use local LL DATA and copy it to the COMMON blocks
!-----MEMORY MANAGEMENT ----------------------------------------------------- ! 1
      PARAMETER(LNAM=10000,LTEX=4000,LCC=10000)                               ! 2
      CHARACTER(LEN=8), DIMENSION(LNAM) :: CNAM
      CHARACTER(LEN=132), DIMENSION(LTEX) :: CTEX
      CHARACTER(LEN=132) CBLA                                                  ! 3
      INTEGER NPAR(LNAM,17)                                                   ! 4
      DOUBLE PRECISION CC(LCC)                                                ! 5
      COMMON / CMEM / CNAM,CTEX,CBLA                                          ! 6
      COMMON /  MEM / NPAR, CC, INAM, ITEX, ICC                               ! 7
!---------------------------------------------------------------------------- ! 8
!
!     CNAM: CONTAINS INAM NAMES OF VARIABLES, FUNCTIONS, SUBROUTINES
!     NPAR: PARAMETERS DESCRIBING THE VARIABLES
!           1: KIND (1:VAR, 2: FUN, 3:SUB)
!           2: TYPE (1:REAL, 2:INTEGER, 3:CHARACTER, 4:DA)
!           3: ADDRESS
!           4: INDICES (FOR ARRAY), NUMBER OF ARGUMENTS (FOR FUN,SUB)
!           11-17: BOUNDS (FOR ARRAYS)
!
!-----CODE ------------------------------------------------------------------! 1
      PARAMETER (LARI=10000)                                                 ! 2
      INTEGER NARI(LARI,11), IARI                                            ! 3
      COMMON / CODE / NARI, IARI                                             ! 4
!----------------------------------------------------------------------------! 5
!-----SYMBOL---------------------------------------------------------------- ! 1
      PARAMETER (LFUNC=100,LOPER=16)                                         ! 2
      INTEGER KFUN(LFUNC)
      INTEGER llKFUN(LFUNC)
      CHARACTER(LEN=6), DIMENSION(LOPER) :: OPER
      CHARACTER(LEN=6), DIMENSION(LFUNC) :: FUNC                             ! 3
      CHARACTER(LEN=6), DIMENSION(LOPER) :: llOPER
      CHARACTER(LEN=6), DIMENSION(LFUNC) :: llFUNC                           ! 3
      COMMON / SYMBOL1 / KFUN
      COMMON / SYMBOL / OPER, FUNC                                           ! 4
!----------------------------------------------------------------------------! 5
      integer i,j
      character(len=132) blanks
      character(len=8) blank8

      PARAMETER(ZERO=0.0D0)

!     OPER     : CONTAINS LOPER NAMES OF SUPPORTED BINARY OPERTORS
!     FUNC     : CONTAINS LFUNC NAMES OF SUPPORTED INTRINSIC FUNCTIONS

      DATA llOPER /  '+     ','-     ','*     ','/     ','^     ', &
                     'DIM   ','DIST  ','MIN   ','MAX   ','MOD   ', & !2
                     'EQ    ','NE    ','GT    ','LT    ','GE    ', & !3
                     'LE    ' /                                      !4

      DATA llFUNC /  'COS   ','SIN   ','TAN   ','ACOS  ','ASIN  ', &
                     'ATAN  ','COSD  ','SIND  ','TAND  ','ACOSD ', & !2
                     'ASIND ','ATAND ','COSH  ','SINH  ','TANH  ', & !3
                     'COSHD ','SINHD ','TANHD ','EXP   ','LOG   ', & !4
                     'LOG2  ','LOG10 ','NINT  ','ABS   ','SIGN  ', & !5
                     'SQRT  ','FAC   ','SQR   ','ISRT  ','REAL  ', & !6
                     'DBLE  ','FLOAT ','SNGL  ','DBLE  ','DCOS  ', & !7
                     'DSIN  ','DTAN  ','DACOS ','DASIN ','DATAN ', & !8
                     'DCOSD ','DSIND ','DTAND ','DACOSD','DASIND', & !9
                     'DATAND','DCOSH ','DSINH ','DTANH ','DEXP  ', & !10
                     'DLOG  ','DLOG2 ','DLOG10','DNINT ','DABS  ', & !11
                     'DSIGN ','DSQRT ','INT   ','IFIX  ','IDINT ', & !12
                     'AINT  ','DINT  ','ANINT ','DNINT ','IDNINT', & !13
                     'ABS   ','DABS  ','IABS  ','DSQRT ','DEXP  ', & !14
                     'ALOG  ','LOG10 ','ALOG10','DLOG10','DMIN1 ', & !15
                     'MAX   ','AMAX1 ','AMAX  ','MAX0  ','MAX1  ', & !16
                     'AMAX0 ','MIN   ','AMIN1 ','AMIN  ','MIN0  ', & !17
                     'MIN1  ','AMIN0 ','DIM   ','IDIM  ','MOD   ', & !18
                     'AMOD  ','SIGN  ','ISIGN ','DPROD ','ATAN2 ', & !19
                     'DATAN2','DMAX1 ','DDIM  ','DSIGN ','DMOD  '/   !20

      DATA llKFUN / 74*1,26*2 /
!   Eric left unchanged, but trying a DO loop later
!     DATA llCBLA / '
!    *                           ' /

      DATA llINAM / 0 /

      do i=1,132
        blanks(i:i)=' '     
      enddo
!      blank6='      '
      blank8='        '

      do i=1,lnam
        cnam(i)=blank8
      enddo

      do i=1,ltex
        ctex(i)=blanks
      enddo

      cbla=blanks

      do i=1,lnam
        do j=1,17
          npar(i,j)=0
        enddo
      enddo

      do i=1,lcc
        cc(i)=ZERO
      enddo

      itex=0
      icc=0

      do i=1,lari
        do j=1,11
          nari(i,j)=0
        enddo
      enddo
      iari=0

      do i=1,100
        kfun(i)=llkfun(i)
      enddo

      do i=1,lnam
        do j=1,17
          npar(i,j)=0
        enddo
      enddo

      do i=1,lfunc
        func(i)=llfunc(i)
      enddo

      do i=1,loper
        oper(i)=lloper(i)
      enddo

      inam=llINAM

      END

! ANFANG UNTERPROGRAMM
  
      SUBROUTINE PRECOM
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT, ERROR_UNIT
!     *****************
!
!     THIS IS THE FOX EXTENDED FORTRAN PRECOMILER.
!     IT READS A FILE AND TRANSFORMS IT TO REGULAR FORTRAN.
!
!-----MEMORY MANAGEMENT -----------------------------------------------------! 1
      PARAMETER(LNAM=10000,LTEX=4000,LCC=10000)                              ! 2
      CHARACTER CNAM(LNAM)*8,CTEX(LTEX)*132,CBLA*132                         ! 3
      INTEGER NPAR(LNAM,17)                                                  ! 4
      DOUBLE PRECISION CC(LCC)                                               ! 5
      COMMON / CMEM / CNAM,CTEX,CBLA                                         ! 6
      COMMON /  MEM / NPAR, CC, INAM, ITEX, ICC                              ! 7
!----------------------------------------------------------------------------! 8
!-----CODE ------------------------------------------------------------------! 1
      PARAMETER (LARI=10000)                                                 ! 2
      INTEGER NARI(LARI,11), IARI                                            ! 3
      COMMON / CODE / NARI, IARI                                             ! 4
      CHARACTER(LEN=480) NAMEDAL
      PARAMETER(MNAME=1000)
      COMMON / DADAL / NAMEDAL(MNAME),icount
!----------------------------------------------------------------------------! 5

      CHARACTER(LEN=10000) A
      CHARACTER(LEN=2) CID
      CHARACTER(LEN=8) BLANK,DNAM
      CHARACTER(LEN=6) PREC
      CHARACTER(LEN=64) IDAT
    
      INTEGER NA(50)

      DATA  BLANK / '        ' /

      IDN = 0
 100  CONTINUE

      CALL GETCOM('PRECOM',A,10000,LA,IEND)

      CALL POSFRA(A,1,LA,' ',NA,200,IA)
      IF(IEND.EQ.1) RETURN

      CID = A(NA(1):NA(1))//A(NA(2):NA(2))

      IF(INDEX(A(1:LA),'=').NE.0) THEN
         IARI = 0
         CALL SYNTAX(A,1,LA, 0,0,2,LERR)
         IF(LERR.EQ.0) CALL ARIFOR(1,IARI,2)
      ELSEIF(CID.EQ.'BD') THEN
         INAM = 0
         ITEX = 0
      ELSEIF(CID.EQ.'ED') THEN
         WRITE(2,'(A)')      '!FOX{'

         DO 30 I=1,INAM
         IF(NPAR(I,1).NE.1) GOTO 30
         IF(NPAR(I,2).NE.4) GOTO 30
         IT = NPAR(I,5) - 1
         ID = NPAR(I,4)
         IF(ID.EQ.0) THEN
            WRITE(2,'(A)') '      INTEGER '//CNAM(I)
         ELSE
            WRITE(2,'(2A,44A1,5(/''     *   '',60A1))')                 &
                   '      INTEGER ',CNAM(I),'(',                        &
            ( ( CTEX(IT+J+2)(K:K),K=1,ILAST(CTEX(IT+J+2),1,132) ),',',   &
              J=1,ID-1 ),                                               &
             ( CTEX(IT+ID+2)(K:K),K=1,ILAST(CTEX(IT+ID+2),1,132) ),')'
         ENDIF
  30     CONTINUE

         WRITE(2,'(6X,A)') 'INTEGER ISCRDA, ISCRRI,IDAO'
         WRITE(2,'(6X,A)') 'REAL(KIND=fPrec) RSCRRI'
         WRITE(2,'(6X,A)') 'COMMON/DASCR/ISCRDA(100),RSCRRI(100)'       &
                               //',ISCRRI(100),IDAO'

         DO 60 I=1,INAM
!         IF((NPAR(I,1).EQ.1).AND.(NPAR(I,2).EQ.4).AND.(NPAR(I,6).EQ.0))
!     *      WRITE(2,'(A)') '      SAVE '//CNAM(I)
  60     CONTINUE

!frs         READ(1,'(A6,i1,1X,A64)',END=60) PREC,INITIAL,IDAT
         READ(1,'(A6,i1,1X,A64)',END=600) PREC,INITIAL,IDAT
         IF(PREC(1:4).ne.'!FOX'.and.(INITIAL.ne.0.or.INITIAL.ne.1)) then
           WRITE(ERROR_UNIT,*) 'Sorry, you are not using the format introduced ',&
           'by F.Schmidt'//' which ',                                   &
           'requires an additional line right after the line: ',        &
           '!FOX  E D; for a properly allocation'//                     &
           ' and also deallocation of variables. '//                    &
           'Please insert a line: !FOX 0 or ',                          &
           'e.g.: !FOX 1 if(ierr.eq.0.or.ierr.eq.2) then '//            &
           'This is very handy when variables have to be reinitialised'
           STOP 4
         endif
         icount=0
         if(INITIAL.eq.0) then
!           WRITE(2,'(A)') '      SAVE LFOX0, LFOX1'
           WRITE(2,'(A)') '      DATA LFOX0   / 0 / '
           WRITE(2,'(A)') '      IF(LFOX0.EQ.0) THEN'
           WRITE(2,'(A)') '         LFOX0 = 1'
         else
           WRITE(2,'(A6,A64)') '      ',IDAT
         endif
         WRITE(2,'(A)') '         CALL DAKEY(''FOX V2.1'')'

         DO 70 I=1,INAM
         ID = NPAR(I,4)
         IT = NPAR(I,5) - 1
         IF(NPAR(I,6).EQ.3) THEN
            CNAM(I) = 'LFOX1   '
         ENDIF
         IF((NPAR(I,1).EQ.1).AND.(NPAR(I,2).EQ.4).AND.(NPAR(I,6).NE.1)) then
! Eric's attempt to initialise DAALL varibales (length 1)
!              WRITE(2,'(A)')
!    *         '         '//CNAM(I)//'= 0'
               WRITE(2,'(A,38A1,5(/''     *   '',60A1))')               &
               '         CALL DAALL('//CNAM(I)//',1',                   &
         ('*','(',(CTEX(IT+J+2)(K:K),K=1,ILAST(CTEX(IT+J+2),1,132)),')', &
                  J=1,ID),                                              &
             ',','''',(CNAM(I)(K:K),K=1,8),' ',' ','''',',',            &
             (CTEX(IT+1)(K:K),K=1,ILAST(CTEX(IT+1),1,132)),              &
             ',',(CTEX(IT+2)(K:K),K=1,ILAST(CTEX(IT+2),1,132)),')'
           icount=icount+1
           if(icount.gt.mname) then
             write(2,*) 'C   Number of variables to be allocated and ' 
             write(2,*) 'C   deallocated is larger than the parameter ', 'MNAME: ',MNAME
             write(2,*) 'C   Change in program dafor.f'  
             WRITE(OUTPUT_UNIT,*) 'C   Number of variables to be allocated and ' 
             WRITE(OUTPUT_UNIT,*) 'C   deallocated is larger than the parameter ', 'MNAME: ',MNAME
             WRITE(OUTPUT_UNIT,*) 'C   Change in program dafor.f'  
             stop 5
           endif
               write(NAMEDAL(icount),'(A,38A1,5(/''     *   '',60A1))') &
               '        CALL DADAL('//CNAM(I)//',1',                    &
         ('*','(',(CTEX(IT+J+2)(K:K),K=1,ILAST(CTEX(IT+J+2),1,132)),')', J=1,ID),')'
         endif
         IF(NPAR(I,6).EQ.3) THEN
            CNAM(I) = DNAM
         ENDIF
  70     CONTINUE

         WRITE(2,'(A)') '      ENDIF'
         WRITE(2,'(A)') '      IDAA = IDAO'
         IF(IDN.NE.0) WRITE(2,'(A)') '      '//DNAM//' = LFOX1'
         WRITE(2,'(A)') '!FOX}'
         IDN = 0
      ELSEIF(CID.EQ.'DV') THEN
         INAM = INAM + 1
         NPAR(INAM,1) = 1
         IF(A(NA(3):NA(3)+1).EQ.'RI') THEN
            NPAR(INAM,2) = 1
         ELSEIF(A(NA(3):NA(3)+1).EQ.'RE') THEN
            NPAR(INAM,2) = 1
         ELSEIF(A(NA(3):NA(3)+1).EQ.'IN') THEN
            NPAR(INAM,2) = 2
         ELSEIF(A(NA(3):NA(3)+1).EQ.'DA') THEN
            NPAR(INAM,2) = 4
         ELSE
            WRITE(2,'(1X,A)') '### ERROR, UNKNOWN TYPE'
            WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, UNKNOWN TYPE'
            INAM = INAM - 1
            GOTO 100
         ENDIF
         IF(A(NA(4):NA(4)+2).EQ.'INT') THEN
            NPAR(INAM,6) = 0
         ELSEIF(A(NA(4):NA(4)+2).EQ.'EXT') THEN
            NPAR(INAM,6) = 1
         ELSEIF(A(NA(4):NA(4)+2).EQ.'COM') THEN
            NPAR(INAM,6) = 2
         ELSEIF(A(NA(4):NA(4)+2).EQ.'FUN') THEN
            NPAR(INAM,6) = 3
            IDN = INAM
         ELSE
            WRITE(2,'(1X,A)') '### ERROR, EXT, INT OR COM NOT SPECIFIED'
            WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, EXT, INT OR COM NOT SPECIFIED'
            INAM = INAM - 1
            GOTO 100
         ENDIF
         IL = ILAST(A,NA(5),NA(6)-1)
         IF(IVCHK(A,NA(5),IL).NE.1) THEN
            WRITE(2,'(1X,A)') '### ERROR, NO VARIABLE FOUND'
            WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, NO VARIABLE FOUND'
            INAM = INAM - 1
            GOTO 100
         ENDIF
         CNAM(INAM) = BLANK
         CNAM(INAM)(1:IL-NA(5)+1) = A(NA(5):IL)
         IF(NPAR(INAM,6).EQ.3) THEN
            DNAM = BLANK
            DNAM = A(NA(5):IL)
         ENDIF
         IF((NPAR(INAM,2).EQ.1).OR.(NPAR(INAM,2).EQ.2)) THEN
            NPAR(INAM,4) = IA - 6
         ELSEIF(NPAR(INAM,2).EQ.4) THEN
            NPAR(INAM,4) = IA - 8
         ENDIF
         IF(NPAR(INAM,4).LT.0) THEN
            WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, MORE ENTRIES REQUIRED'
            WRITE(2,'(1X,A)') '### ERROR, MORE ENTRIES REQUIRED'
            INAM = INAM - 1
            GOTO 100
         ELSEIF(NPAR(INAM,4).GT.7) THEN
            WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, TOO MANY DIMENSIONS'
            WRITE(2,'(1X,A)') '### ERROR, TOO MANY DIMENSIONS'
            INAM = INAM - 1
            GOTO 100
         ELSEIF((NPAR(INAM,4).NE.0).AND.(NPAR(INAM,6).EQ.3)) THEN
            WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, NO ARRAYS IN FUNCTION NAME'
            WRITE(2,'(1X,A)') '### ERROR, NO ARRAYS IN FUNCTION NAME'
            INAM = INAM - 1
            GOTO 100
         ENDIF
         IF(ITEX+IA-5.GT.LTEX) THEN
            WRITE(ERROR_UNIT,'(1X,A)') '!!! ERROR IN PRECOM, CTEX EXHAUSTED'
            WRITE(2,'(1X,A)') '!!! ERROR IN PRECOM, CTEX EXHAUSTED'
            STOP 6
         ENDIF
         NPAR(INAM,5)  = ITEX + 1
         DO 80 I=6,IA - 1
         ITEX = ITEX + 1
         CTEX(ITEX) = CBLA
         IL = ILAST(A,NA(I),NA(I+1)-1)
         CTEX(ITEX)(1:IL-NA(I)+1)= A(NA(I):IL)
  80     CONTINUE
      ELSEIF(CID.EQ.'DF') THEN
         INAM = INAM + 1
         NPAR(INAM,1) = 2
         IF(A(NA(3):NA(3)+1).EQ.'RI') THEN
            NPAR(INAM,2) = 1
         ELSEIF(A(NA(3):NA(3)+1).EQ.'RE') THEN
            NPAR(INAM,2) = 1
         ELSEIF(A(NA(3):NA(3)+1).EQ.'IN') THEN
            NPAR(INAM,2) = 2
         ELSEIF(A(NA(3):NA(3)+1).EQ.'DA') THEN
            NPAR(INAM,2) = 4
         ELSE
            WRITE(2,'(1X,A)') '### ERROR, UNKNOWN TYPE'
            WRITE(OUTPUT_UNIT,'(1X,A)') '### ERROR, UNKNOWN TYPE'
            INAM = INAM - 1
            GOTO 100
         ENDIF
         IL = ILAST(A,NA(4),NA(5)-1)
         IF(IVCHK(A,NA(4),IL).NE.1) THEN
            WRITE(2,'(1X,A)') '### ERROR, NO FUNCTION FOUND'
            WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, NO FUNCTION FOUND'
            INAM = INAM - 1
            GOTO 100
         ENDIF
         CNAM(INAM) = BLANK
         CNAM(INAM)(1:IL-NA(4)+1) = A(NA(4):IL)
         NPAR(INAM,4) = INDEX('123456789',A(NA(5):NA(5)))
         IF(NPAR(INAM,4).EQ.0) THEN
            WRITE(2,'(1X,A)') '### ERROR IN NUMBER OF ARGUMENTS'
            WRITE(ERROR_UNIT,'(1X,A)') '### ERROR IN NUMBER OF ARGUMENTS'
            INAM = INAM - 1
            GOTO 100
         ENDIF
      ELSE
         LERR = 1
         WRITE(ERROR_UNIT,'(1X,2A)') '### ERROR, UNKNOWN COMMAND ',CID
         WRITE(2,'(1X,2A)') '### ERROR, UNKNOWN COMMAND ',CID
      ENDIF

      GOTO 100
 600  CONTINUE
      END

! ANFANG UNTERPROGRAMM
  
      SUBROUTINE GETCOM(CMOD,A,NA,IA,IEND)
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : ERROR_UNIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     ************************************
!
!     THIS SUBROUTINE GETS THE NEXT COMMAND FROM THE INPUT FILE AND STORES
!     IT IN THE CHARACTER A

      INTEGER MA(50)

!-----MEMORY MANAGEMENT -----------------------------------------------------! 1
      PARAMETER(LNAM=10000,LTEX=4000,LCC=10000)                              ! 2
      CHARACTER CNAM(LNAM)*8,CTEX(LTEX)*132,CBLA*132                         ! 3
      INTEGER NPAR(LNAM,17)                                                  ! 4
      DOUBLE PRECISION CC(LCC)                                               ! 5
      COMMON / CMEM / CNAM,CTEX,CBLA                                         ! 6
      COMMON /  MEM / NPAR, CC, INAM, ITEX, ICC                              ! 7
      CHARACTER NAMEDAL*480
      PARAMETER(MNAME=1000)
      COMMON / DADAL / NAMEDAL(MNAME),icount
!----------------------------------------------------------------------------! 8

      CHARACTER CMOD*6,PREC*6,REST*8,A*(*),ALIN*10000

      SAVE ALIN,IAMAX
      DATA IAMAX / 0 /

      IEND = 0
      IA   = 0

      DO I=1,NA
        A(I:I) = ' '
      END DO

  20  CONTINUE
      IF(IAMAX.NE.0) THEN
         IF(INDEX(ALIN(1:IAMAX),';').NE.0) GOTO 40
      ENDIF

      IF(IAMAX+132.GT.10000) THEN
         WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, COMMAND TOO LONG'
         WRITE(2,'(1X,A)') '### ERROR, COMMAND TOO LONG'
         STOP 7
      ELSEIF(CMOD.EQ.'PRECOM') THEN
         ALIN(IAMAX+1:IAMAX+1) = ' '
         READ(1,'(A6,A125,A8)',END=60) PREC,ALIN(IAMAX+2:IAMAX+126),REST
         if(alin(IAMAX+2:IAMAX+6).eq.'DADAL') then
           if(icount.gt.0) then
             do 70 i=icount,1,-1
               write(2,'(a132)') namedal(i)
 70          continue
           endif
         endif

!        CHECKING IF LINE CONTAINS DA VARIABLE
!        *************************************

         IF((PREC(1:1).NE.'*'.AND.PREC(1:1).NE.'C').OR.PREC(1:4).EQ.'!FOX'.OR.PREC(1:3).EQ.'!DA') THEN
            DO 30 I=1,INAM
            ICN = ILAST(CNAM(I),1,8)
            IF(NPAR(I,2).NE.4) GOTO 30
            IPOS = IAMAX+2
   25       CONTINUE
            INDI = INDEX(ALIN(IPOS:IAMAX+126),CNAM(I)(1:ICN))
            IF(INDI.EQ.0) GOTO 30
            IPOS = IPOS + INDI
            IF(IPOS.NE.IAMAX+3) THEN
               IF(INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ_1234567890', ALIN(IPOS-2:IPOS-2)).NE.0) GOTO 25
            ENDIF
            IF(IPOS+ICN-1.NE.IAMAX+126) THEN
               IF(INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ_1234567890', ALIN(IPOS+ICN-1:IPOS+ICN-1)).NE.0) GOTO 25
            ENDIF
            REST = '!FOX'
            GOTO 35
   30       CONTINUE
   35       CONTINUE
         ENDIF

!        CHECKING IF LINE IS END STATEMENT
!        *********************************

         CALL POSFRA(ALIN,IAMAX+1,IAMAX+126,' ',MA,200,IA)
         IF(IA.EQ.1) THEN
            IF(ALIN(MA(1):MA(1)+3).EQ.'END ') INAM = 0
         ENDIF

!        WRITING LINE TO OUTPUT
!        **********************

         ILL = ILAST(ALIN(IAMAX+1:IAMAX+126),2,126)
         IF(REST.NE.'        ') THEN
            WRITE(2,'(A6,A125,A)') PREC,ALIN(IAMAX+2:IAMAX+126), REST(1:ILAST(REST,1,8))
         ELSE
            WRITE(2,'(A6,A)') PREC,ALIN(IAMAX+2:IAMAX+ILL)
         ENDIF
         IF(PREC(1:4).EQ.'!FOX'.OR.PREC(1:3).EQ.'!DA') IAMAX = IAMAX + ILL

!     DIRECT COMPILE MODE
!     *******************

      ELSE
         READ(1,'(A132)',END=60) PREC,ALIN(IAMAX+1:IAMAX+132)
         IAMAX = IAMAX + ILAST(ALIN(IAMAX+1:IAMAX+132),1,132)
      ENDIF
      GOTO 20

!     EXTRACTING FOX COMMAND
!     **********************

  40  CONTINUE
      IA = INDEX(ALIN(1:IAMAX),';')
      IF(IA.GT.NA) THEN
         WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, COMMAND TOO LONG'
         WRITE(2,'(1X,A)') '### ERROR, COMMAND TOO LONG'
         STOP 8
      ELSEIF(IA.EQ.0) THEN
         GOTO 20
      ENDIF
      A(1:IA) = ALIN(1:IA)
      DO I=IA+1,IAMAX
        ALIN(I-IA:I-IA) = ALIN(I:I)
      END DO
      IAMAX = IAMAX - IA

!     REPLACING '**' BY '^ '
!     **********************

      DO 55 I=1,IA
      IF(A(I:I).EQ.'*') THEN
         DO 54 J=I+1,IA
         IF(A(J:J).NE.' ') THEN
            IF(A(J:J).EQ.'*') THEN
               A(I:I) = '^'
               A(J:J) = ' '
            ENDIF
            GOTO 55
         ENDIF
  54     CONTINUE
      ENDIF
  55  CONTINUE

      RETURN

  60  CONTINUE
      IA = 0
      IEND = 1
      RETURN
      END

! ANFANG UNTERPROGRAMM
  
      SUBROUTINE POSFRA(A,IA1,IA2,CDEL,NA,LA,IA)
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : ERROR_UNIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     ******************************************
!
!     THIS SUBROUTINE TAKES THE INFORMATION IN CHARACTER A BETWEEN IA1 AND
!     IA2 AND DETERMINES THE BEGINNING POSITIONS OF FRAGMENTS, SEPERATED BY
!     DELIMITER CDEL, STORING THEM IN NA

      CHARACTER A*(*),CDEL*1
      INTEGER NA(LA)

      IF(A(IA1:IA1).NE.CDEL) THEN
         IA = 1
         NA(IA) = IA1
      ELSE
         IA = 0
      ENDIF

      DO 100 I=IA1+1,IA2
      IF(A(I-1:I-1).EQ.CDEL.AND.A(I:I).NE.CDEL) THEN
         IA = IA + 1
         IF(IA.GT.LA) THEN
            WRITE(ERROR_UNIT,'(1X,A)') '!!! ERROR, TOO MANY FRAGS IN POSFRA'
            WRITE(2,'(1X,A)') '!!! ERROR, TOO MANY FRAGS IN POSFRA'
            STOP 9
         ENDIF
         NA(IA) = I
      ENDIF
 100  CONTINUE

      RETURN
      END

      INTEGER FUNCTION ILAST(A,IA1,IA2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     *********************************
!
!     THIS FUNCTION DETERMINES THE LAST NONBLANK POSITION IN CHARACTER A
!     BETWEEN POSITIONS IA1 AND IA2

      CHARACTER A*(*)

! Eric sets ILAST when line has no trailing blanks
      ILAST = IA2

      DO 10 ILAST = IA2,IA1,-1
      IF(A(ILAST:ILAST).NE.' ') RETURN
  10  CONTINUE

      ILAST = IA1
      RETURN
      END

! ANFANG UNTERPROGRAMM
  
      SUBROUTINE CAP(A,IA1,IA2)
!     *************************
!
!     THIS SUBROUTINE CAPITALIZES POSITIONS IA1 THROUGH IA2 IN CHARACTER A

      CHARACTER A*(*)
      SAVE ICALL,ICSA,ICSZ,IDIF
      DATA ICALL / 0 /

      IF(ICALL.EQ.0) THEN
         ICALL = 1
         ICSA = ICHAR('A')
         ICSZ = ICHAR('Z')
         ICCA = ICHAR('A')
         IDIF = ICCA - ICSA
      ENDIF

      DO 100 I=IA1,IA2
      ICI = ICHAR(A(I:I))
      IF(ICI.GE.ICSA.AND.ICI.LE.ICSZ) THEN
         ICI = ICI + IDIF
         A(I:I) = CHAR(ICI)
      ENDIF
 100  CONTINUE

      RETURN
      END

      INTEGER FUNCTION IVCHK(A,IA1,IA2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     *********************************
!
!     THIS FUNCTION RETURNS 1 IF THE SUBSTRING A(IA1:IA2) IS A SYNTACTICALLY
!     LEGAL NAME OF A VARIABLE, ZERO OTHERWISE.

      CHARACTER A*(*)

      IVCHK = 0

      IF(IA2.LT.IA1) RETURN
      IF(IA2-IA1.GT.7) RETURN
      IF(INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ',A(IA1:IA1)).EQ.0) RETURN

      ICH = 0

      DO 100 I=IA1+1,IA2
      IF(INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_', A(IA1:IA1)).EQ.0) THEN
        RETURN
      ELSE
        ICH = ICH + 1
      ENDIF
  100 CONTINUE

      IVCHK = 1
      RETURN
      END

! ANFANG UNTERPROGRAMM
  
      SUBROUTINE SYNTAX(B, IB1,IB2, IX, IN, IU, LERR)
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT, ERROR_UNIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     ***********************************************
!
!     THIS SUBROUTINE DECODES THE CODE STORED ON CHARACTER B FROM IB1 TO IB2
!     INTO ELEMENTARY OPERATIONS AND STORES THE RESULT IN NARI.
!
!     IF IX = 0, THE CODE IS A FULL ASSIGNMENT INCLUDING EQUAL SIGN;
!     IF IX <>0, THE CODE IS AN ARITHMETIC EXPRESSION WITHOUT EQUAL SIGN, AND
!                THE RESULT IS STORED IN IX.
!     IF IN = 1, INTERPRETABLE CODE IS REQUESTED,
!     IF IN = 0, COMPILED CODE IS REQUESTED.
!     IF IU > 0, SYNTAX ERRORS ARE PRINTED TO UNIT IU.
!
!     ON RETURN,
!     LERR = 0   MEANS THAT NO SYNTAX ERRORS WERE DETECTED IN B,
!     LERR = 1   MEANS THAT THERE WERE SYNTAX ERRORS AND NO CODE WAS GENERATED.
!
!     NARI     : CONTAINS REQUESTED IARI ELEMENTARY OPERATIONS
!                1       : ADDRESS WHERE RESULT IS TO BE STORED
!                2       : TYPE (1:ASSIGNMENT, 2:OPERATOR, 3:FUNCTION, 4:ARRAY)
!                3       : ID # OF TYPE
!                4       : ADDRESS OF SINGLE OPERAND (IF NEEDED)
!                5,6,... : ADDRESSES OF SEVERAL MORE OPERANDS (IF NEEDED)
!                ADDRESSES IA ARE CODED IN THE FOLLOWING WAY:
!                       IA < 0:      -IA ADDRESS OF SCRATCH VARIABLE
!                0    < IA < LNAM:   IA ADDRESS OF DECLARED QUANTITY
!                LNAM < IA:          IA - LNAM ADDRESS OF CONSTANT
!
!
!-----MEMORY MANAGEMENT -----------------------------------------------------! 1
      PARAMETER(LNAM=10000,LVAR=8,LTEX=4000,LCC=10000)                       ! 2
      CHARACTER CNAM(LNAM)*8,CTEX(LTEX)*132,CBLA*132                         ! 3
      INTEGER NPAR(LNAM,17)                                                  ! 4
      DOUBLE PRECISION CC(LCC)                                               ! 5
      COMMON / CMEM / CNAM,CTEX,CBLA                                         ! 6
      COMMON /  MEM / NPAR, CC, INAM, ITEX, ICC                              ! 7
!----------------------------------------------------------------------------! 8
!-----CODE ------------------------------------------------------------------! 1
      PARAMETER (LARI=10000)                                                 ! 2
      INTEGER NARI(LARI,11), IARI                                            ! 3
      COMMON / CODE / NARI, IARI                                             ! 4
!----------------------------------------------------------------------------! 5
!-----SYMBOL-----------------------------------------------------------------! 1
      PARAMETER (LFUNC=100,LOPER=16)                                         ! 2
      INTEGER KFUN(LFUNC)
      CHARACTER OPER(LOPER)*6, FUNC(LFUNC)*6                                 ! 3
      COMMON / SYMBOL1 / KFUN
      COMMON / SYMBOL / OPER, FUNC                                           ! 4
!----------------------------------------------------------------------------! 5
!
!frs      PARAMETER(LANA=100,LA=10000,LCHECK=0,LDEC=37)
      PARAMETER(LANA=100,LA=10000,LDEC=37)
      PARAMETER(LSCR=99)

      PARAMETER(ZERO=0D0,ONE=1.0D0,TEN=10.0D0)

      INTEGER NANA(LANA,5)
      CHARACTER B*(*),A*10000
      CHARACTER OPS*6,CANA*10!,CTYP*7!,CCNA*1
      CHARACTER BLANK*8,NUM*10,LET*26,SEARCH*8,AER*50
      DOUBLE PRECISION CDEC(0:LDEC)
      SAVE IFI, CDEC

      DATA OPS  / '+-*/^#' /
      DATA CANA / '=OFAVSC#,;' /
!      DATA CTYP / 'RICDGPB' /

      DATA NUM  / '1234567890' /
      DATA LET  / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
      DATA AER  / '                                                  ' /
      DATA IFI  / 0 /

      INTEGER IPRIO
! Eric
      IB = 0

      IPRIO = -4242

      IFUVA = 0
      IF(IFI.EQ.0) THEN
         IFI = 1
         CDEC(0) = ONE
         DO J=1,LDEC
           CDEC(J) = CDEC(J-1) * TEN
         END DO
      ENDIF

!-----------------------------------------------------------------------

      IF(IB2.GE.LA) THEN
         WRITE(ERROR_UNIT,'(1X,A)') '!!! ERROR IN SYNTAX, CHARACTER B TOO LARGE'
         WRITE(2,'(1X,A)') '!!! ERROR IN SYNTAX, CHARACTER B TOO LARGE'
         STOP 10
      ENDIF

!frs
!frs      IF(LCHECK.EQ.1) WRITE(OUTPUT_UNIT,'(1X,A)') B(IB1:IB2)
!frs
!
!     DELETING BLANKS IN CODE LINE AND CAPITALIZING
!     *********************************************

      IA = 0
      DO 10 I=IB1,IB2
      IF(B(I:I).NE.' ') THEN
         IA = IA + 1
         IF(IA.GT.LA) THEN
            AER = 'ERROR, INPUT LINE TOO LONG'
            GOTO 1000
         ENDIF
                 A(IA:IA) = B(I:I)
      ENDIF
  10  CONTINUE
      CALL CAP(A,1,IA)
      IA = IA + 1
      A(IA:IA) = ';'

!     SETTING START VALUES
!     ********************

      DO I=1,LVAR
        BLANK(I:I) = ' '
      END DO

!      IARIA = IARI+1
      ILEFT = IX
      ISCR  = 0
      LERR  = 0
      IOP   = 1
      IPAR  = 0
      IEQU  = 0
      I     = 0
!
      IANA  = 1
      NANA(IANA,1) = INDEX(CANA,';')
      NANA(IANA,2) = 0
      NANA(IANA,3) = 0
      NANA(IANA,4) = 0
      NANA(IANA,5) = IANA + 1

!     PROCESSING LEFT SIDE OF EQUAL SIGN
!     **********************************

      IF(ILEFT.NE.0) THEN
         IEQU = 1
         GOTO 100
      ENDIF

      IF(INDEX(LET,A(1:1)).EQ.0) THEN
         AER = 'ERROR, ILLEGAL VARIABLE NAME'
         GOTO 1000
      ENDIF
      IS = 0
      SEARCH = BLANK
  20  I = I + 1
      IF(INDEX(LET//NUM//'_',A(I:I)).EQ.0.OR.IS.GT.LVAR) GOTO 21
      IS = IS + 1
      SEARCH(IS:IS) = A(I:I)
      GOTO 20
  21  CONTINUE
      DO 22 ILEFT=1,INAM
      IF(SEARCH.EQ.CNAM(ILEFT)) GOTO 23
  22  CONTINUE
      AER = 'ERROR, UNKNOWN VARIABLE'
      GOTO 1000
  23  CONTINUE
      IF(A(I:I).EQ.'=') THEN
         IF(NPAR(ILEFT,4).NE.0) THEN
            AER = 'ERROR, VARIABLE IS DECLARED AS ARRAY'
            GOTO 1000
         ENDIF
      ELSEIF(A(I:I).EQ.'(') THEN
         IFUVA = 1
         IF(NPAR(ILEFT,4).EQ.0) THEN
            AER = 'ERROR, VARIABLE IS NOT DECLARED AS ARRAY'
            GOTO 1000
         ENDIF
      ELSE
         AER = 'ERROR, EXPECTED OPEN PARENTHESIS OR EQUAL SIGN'
         GOTO 1000
      ENDIF
      GOTO 100

!     EQUAL SIGN
!     **********

  50  CONTINUE
      IF(A(I:I).NE.'=') THEN
         AER = 'ERROR, EXPECTED EQUAL SIGN'
         GOTO 1000
      ELSEIF(IEQU.NE.0) THEN
         AER = 'ERROR, EQUAL SIGN NOT ALLOWED'
         GOTO 1000
      ELSEIF(IPAR.NE.0) THEN
         AER = 'ERROR, UNBALANCED PARENTHESIS AT EQUAL SIGN'
         GOTO 1000
      ENDIF
      IANA = IANA + 1
      NANA(IANA,1) = INDEX(CANA,'=')
      NANA(IANA,2) = ILEFT
      NANA(IANA,3) = 1
      NANA(IANA,4) = IANA - 1
      NANA(IANA,5) = IANA + 1
      IEQU = 1
      I = I + 1
      GOTO 100

!     SWITCH TO EQUAL, OPEN PARENTHESIS, SIGN, CONSTANT, FUNCTION OR VARIABLE
!     ***********************************************************************

 100  CONTINUE
      IF(A(I:I).EQ.'=') GOTO 50
      IF(A(I:I).EQ.'(') GOTO 110
      IF(A(I:I).EQ.'+'.OR.A(I:I).EQ.'-') GOTO 170
      IF(INDEX(NUM,A(I:I)).NE.0.OR.A(I:I).EQ.'.') GOTO 120
      IF(A(I:I).EQ.'[') GOTO 130
      GOTO 101

!     DECODING NAME OF FUNCTION OR VARIABLE
!     *************************************

 101  IF(INDEX(LET,A(I:I)).EQ.0) THEN
         AER = 'ERROR, ILLEGAL NAME OF VARIABLE OR FUNCTION'
         GOTO 1000
      ENDIF
      SEARCH = BLANK
      SEARCH(1:1) = A(I:I)
      IS = 1
 102  I = I + 1
      IS = IS + 1
      IF(INDEX(LET//NUM//'_',A(I:I)).EQ.0.OR.IS.GT.LVAR) GOTO 103
      SEARCH(IS:IS) = A(I:I)
      GOTO 102
 103  CONTINUE
      IS = IS - 1
      DO 104 II=1,INAM
      IF(SEARCH.EQ.CNAM(II)) THEN
         IF(NPAR(II,1).EQ.1) THEN
            GOTO 160
         ELSEIF(NPAR(II,1).EQ.2) THEN
            GOTO 150
         ENDIF
      ENDIF
 104  CONTINUE
      DO 105 II=1,LFUNC
      IF((SEARCH.EQ.FUNC(II)).AND.(IS.LE.6)) GOTO 150
 105  CONTINUE
      AER = 'ERROR, UNKNOWN VARIABLE OR FUNCTION'
      GOTO 1000

!     OPEN PARENTHESIS
!     ****************

 110  CONTINUE
      IF(A(I:I).NE.'(') THEN
         AER = 'ERROR, EXPECTED OPENING PARENTHESIS'
         GOTO 1000
      ENDIF
      I = I + 1
      IPAR = IPAR + 1
      IF(IPAR.GT.20) THEN
         AER = 'ERROR, PARENTHESES NESTED TOO DEEPLY'
         GOTO 1000
      ENDIF
      GOTO 100

!     CONSTANT
!     ********

 120  CONTINUE
      CVAL = ZERO
      IFRA = 0
      LPOI = 0
      IEXP = 0
      LEXP = 0
      IPME = 0
      IOLD = I
      I = I - 1
 121  I = I + 1
      IF(I-IOLD.GE.132) THEN
         AER = 'ERROR, TOO MANY DIGITS IN CONSTANT'
         GOTO 1000
      ELSEIF(INDEX(NUM,A(I:I)).NE.0) THEN
         IF(LEXP.NE.0) THEN
            IEXP = IEXP*10 + INDEX(NUM(1:9),A(I:I))
         ELSEIF(LPOI.EQ.0) THEN
            CVAL = CVAL*TEN + INDEX(NUM(1:9),A(I:I))
         ELSEIF(LPOI.EQ.1) THEN
            IFRA = IFRA + 1
            IF(IFRA.GT.LDEC) THEN
               AER = 'ERROR, TOO MANY DIGITS AFTER DECIMAL POINT'
               GOTO 1000
            ENDIF
            CVAL = CVAL + INDEX(NUM(1:9),A(I:I)) / CDEC(IFRA)
         ELSE
            AER = 'ERROR, LPOI NOT 0 OR 1'
            GOTO 1000
         ENDIF
      ELSEIF(A(I:I).EQ.'.') THEN
         IF(LPOI.NE.0.OR.LEXP.NE.0) THEN
            AER = 'ERROR IN SYNTAX OF CONSTANT'
            GOTO 1000
         ENDIF
         LPOI = 1
      ELSEIF(INDEX('DE',A(I:I)).NE.0) THEN
         IF(LEXP.NE.0) THEN
            AER = 'ERROR IN SYNTAX OF CONSTANT'
            GOTO 1000
         ENDIF
         LEXP = 1
      ELSEIF(INDEX('+-',A(I:I)).NE.0) THEN
         IF(LEXP.EQ.0.OR.(LEXP.NE.0.AND.INDEX('DE',A(I-1:I-1)).EQ.0)) THEN
            GOTO 122
         ELSEIF(A(I:I).EQ.'+') THEN
            IPME = +1
         ELSEIF(A(I:I).EQ.'-') THEN
            IPME = -1
         ENDIF
      ELSE
         GOTO 122
      ENDIF
      GOTO 121
 122  CONTINUE
      IF(LEXP.NE.0) THEN
         IF(IEXP.LT.0.OR.IEXP.GT.LDEC) THEN
            AER = 'ERROR, EXPONENT OF CONSTANT OUT OF RANGE'
            GOTO 1000
         ENDIF
         IF((IPME.EQ.0).OR.(IPME.EQ.+1)) THEN
            CVAL = CVAL * CDEC(IEXP)
         ELSEIF(IPME.EQ.-1) THEN
            CVAL = CVAL / CDEC(IEXP)
         ELSE
            AER = 'ERROR, WRONG VALUE OF IPME'
            GOTO 1000
         ENDIF
      ENDIF
      ITEX = ITEX + 1
      ICC = ICC + 1
      IF(ICC.GT.LCC) THEN
         AER = 'ERROR, CONSTANT MEMORY OVERFLOW'
         GOTO 1000
      ENDIF
      CC(ICC) = CVAL
      CTEX(ITEX)(1:I-IOLD) = A(IOLD:I-1)
      DO J=I-IOLD+1,132
        CTEX(ITEX)(J:J) = ' '
      END DO
      IANA = IANA + 1
      NANA(IANA,1) = INDEX(CANA,'C')
      IF(IN.EQ.1) NANA(IANA,2) = ICC
      IF(IN.EQ.0) NANA(IANA,2) = ITEX + LNAM
      NANA(IANA,3) = 0
      NANA(IANA,4) = IANA - 1
      NANA(IANA,5) = IANA + 1
      GOTO 200

!     ADDITIONAL OPTION
!     *****************

 130  CONTINUE
      AER = 'ERROR, OPERATION NOT SUPPORTED YET'
      GOTO 1000

!     CLOSING PARENTHESIS
!     *******************

 140  CONTINUE
      IF(A(I:I).NE.')') THEN
         AER = 'ERROR, EXPECTED CLOSING PARENTHESIS'
         GOTO 1000
      ENDIF
      IF(IPAR.EQ.IFUVA) IFUVA = 0
      IPAR = IPAR - 1
      IF(IPAR.LT.0) THEN
         AER = 'ERROR, PARENTHESIS CLOSED TOO EARLY'
         GOTO 1000
      ENDIF
      IF((IPAR.EQ.0).AND.(IEQU.EQ.0).AND.(A(I+1:I+1).NE.'=')) THEN
         AER = 'ERROR, EXPECTED EQUAL SIGN'
         GOTO 1000
      ENDIF
      I = I + 1
      GOTO 200

!     FUNCTION
!     ********

 150  CONTINUE
      IF(SEARCH.EQ.FUNC(II)) II = -II
      IF(A(I:I).NE.'(') THEN
         AER = 'ERROR, EXPECTED OPENING PARENTHESIS'
         GOTO 1000
      ENDIF
      IF(IFUVA.EQ.0) IFUVA = IPAR + 1
      IANA = IANA + 1
      NANA(IANA,1) = INDEX(CANA,'F')
      NANA(IANA,2) = II
      NANA(IANA,3) = 4*(IPAR+1)+4
      NANA(IANA,4) = IANA - 1
      NANA(IANA,5) = IANA + 1
      IOP = MAX(IOP,4*(IPAR+1)+4)
      GOTO 110

!     VARIABLE
!     ********

 160  CONTINUE
      IF(A(I:I).NE.'(') THEN
         IF(NPAR(II,4).NE.0) THEN
            AER = 'ERROR, VARIABLE IS DECLARED AS ARRAY'
            GOTO 1000
         ENDIF
         IANA = IANA + 1
         NANA(IANA,1) = INDEX(CANA,'V')
         NANA(IANA,2) = II
         NANA(IANA,3) = 0
         NANA(IANA,4) = IANA - 1
         NANA(IANA,5) = IANA + 1
         GOTO 200
      ELSE
         IF(NPAR(II,4).EQ.0) THEN
            AER = 'ERROR, VARIABLE IS NOT DECLARED AS ARRAY'
            GOTO 1000
         ENDIF
         IF(IFUVA.EQ.0) IFUVA = IPAR + 1
         IANA = IANA + 1
         NANA(IANA,1) = INDEX(CANA,'A')
         NANA(IANA,2) = II
         NANA(IANA,3) = 4*(IPAR+1)+4
         NANA(IANA,4) = IANA - 1
         NANA(IANA,5) = IANA + 1
         IOP = MAX(IOP,4*(IPAR+1)+4)
         GOTO 110
      ENDIF

!     PLUS OR MINUS SIGN IN FRONT OF VARIABLE,FUNCTION OR CONSTANT
!     ************************************************************

 170  CONTINUE
      IF(A(I:I).EQ.'+') THEN
         I = I + 1
         GOTO 100
      ELSEIF(A(I:I).EQ.'-') THEN
         I = I + 1
         ITEX = ITEX + 1
         CTEX(ITEX)(1:4) = '-ONE'
         DO J=5,132
           CTEX(ITEX)(J:J) = ' '
         END DO
         ICC = ICC + 1
         CC(ICC) = -ONE
         IANA = IANA + 1
         NANA(IANA,1) = INDEX(CANA,'C')
         IF(IN.EQ.1) NANA(IANA,2) = ICC
         IF(IN.EQ.0) NANA(IANA,2) = ITEX + LNAM
         NANA(IANA,3) = 0
         NANA(IANA,4) = IANA - 1
         NANA(IANA,5) = IANA + 1
         IANA = IANA + 1
         NANA(IANA,1) = INDEX(CANA,'O')
         NANA(IANA,2) = 3
         NANA(IANA,3) = (IPAR+1)*4+2
         NANA(IANA,4) = IANA - 1
         NANA(IANA,5) = IANA + 1
         IB = IB + 10
         IOP = MAX(IOP,4*(IPAR+1)+2)
         GOTO 100
      ELSE
         AER = 'ERROR WITH PLUS OR MINUS SIGN'
         GOTO 1000
      ENDIF

!     COMMA (IN FUNCTION OR ARRAY)
!     ****************************

 180  CONTINUE
      IF(A(I:I).NE.',') THEN
         AER = 'ERROR, COMMA EXPECTED'
         GOTO 1000
      ENDIF
      IF(IFUVA.EQ.0) THEN
         AER = 'ERROR, COMMA MISPLACED'
         GOTO 1000
      ENDIF
      I = I + 1
      IANA = IANA + 1
      NANA(IANA,1) = INDEX(CANA,',')
      NANA(IANA,2) = 0
      NANA(IANA,3) = 0
      NANA(IANA,4) = IANA - 1
      NANA(IANA,5) = IANA + 1
      GOTO 100

!     SWITCH TO OPERATOR, CLOSING PARENTHESIS, COMMA OR END
!     *****************************************************

 200  CONTINUE
      SEARCH = BLANK
      IF(INDEX(OPS,A(I:I)).NE.0) GOTO 210
      IF(A(I:I).EQ.'=') GOTO 50
      IF(A(I:I).EQ.')') GOTO 140
      IF(A(I:I).EQ.',') GOTO 180
      IF(A(I:I).EQ.';') GOTO 220
      AER = 'ERROR, OPERATOR, COMMA OR CLOSING PARENS EXPECTED'
      GOTO 1000

!     OPERATORS
!     *********

 210  CONTINUE
      SEARCH = BLANK
      IF(A(I:I).EQ.'#') THEN
         IPRIO = 4
         IS = 0
 211     I = I + 1
         IF(A(I:I).EQ.'#') GOTO 212
         IS = IS + 1
         IF(IS.GT.6) GOTO 219
         SEARCH(IS:IS) = A(I:I)
         GOTO 211
      ELSE
         SEARCH(1:1) = A(I:I)
         IF(A(I:I).EQ.'+'.OR.A(I:I).EQ.'-') IPRIO = 1
         IF(A(I:I).EQ.'*'.OR.A(I:I).EQ.'/') IPRIO = 2
         IF(A(I:I).EQ.'^') IPRIO = 3
      ENDIF
 212  CONTINUE
      DO 213 II=1,LOPER
      IF(SEARCH(1:6).EQ.OPER(II)) GOTO 214
 213  CONTINUE
      GOTO 219
 214  CONTINUE
      if (IPRIO.eq.-4242) then
         write (6,*) "ERROR in routine SYNTAX"
         write (6,*) "IPRIO not set."
         stop 17
      endif
      IANA = IANA + 1
      NANA(IANA,1) = INDEX(CANA,'O')
      NANA(IANA,2) = II
      NANA(IANA,3) = 4*(IPAR+1)+IPRIO
      NANA(IANA,4) = IANA - 1
      NANA(IANA,5) = IANA + 1
      IOP = MAX(IOP,4*(IPAR+1)+IPRIO)
      I = I + 1
      GOTO 100
 219  CONTINUE
      AER = 'ERROR, UNSUPPORTED OPERATOR'
      GOTO 1000

!     END OF SOURCE CODE
!     ******************

 220  CONTINUE
      IF(IPAR.NE.0) THEN
         AER = 'ERROR, UNBALANCED PARENTHESES'
         GOTO 1000
      ENDIF
      IF(IEQU.EQ.0) THEN
         AER = 'ERROR, ASSIGNMENT INCOMPLETE'
         GOTO 1000
      ENDIF
      IANA = IANA + 1
      NANA(IANA,1) = INDEX(CANA,';')
      NANA(IANA,2) = 0
      NANA(IANA,3) = 0
      NANA(IANA,4) = IANA - 1
      NANA(IANA,5) = 0

      DO 350 IO = IOP,1,-1
!     ********************

!     EXTRACTING NEXT ELEMENTARY OPERATION
!     ************************************

      IA = 1
 310  CONTINUE
      IF(NANA(IA,3).NE.IO) THEN
         IA = NANA(IA,5)
         IF(IA.EQ.IANA) GOTO 350
         GOTO 310
      ENDIF

!frs
!frs      IF(LCHECK.EQ.1) THEN
!frs         WRITE(OUTPUT_UNIT,'(A,I3,71A)') ' IO = ',IO , ' ',('-',J=1,70)
!frs         WRITE(OUTPUT_UNIT,'(100(6X,8(A1,I3,I3,A1)/))')
!frs     *   (CANA(NANA(J,1):NANA(J,1)),NANA(J,2),NANA(J,3),'|',J=1,IANA)
!frs      ENDIF
!frs

      IARI = IARI + 1

      DO J=1,11
        NARI(IARI,J) = 0
      END DO

!     FILLING IN FIRST THREE ENTRIES OF NARI
!     **************************************

      ITYP = NANA(IA,1)
      IF(CANA(ITYP:ITYP).EQ.'=') THEN
         NARI(IARI,1) = ILEFT
      ELSE
         ISCR = ISCR + 1
         IF(ISCR.GE.LSCR) THEN
            AER = 'ERROR, TOO MANY SCRATCH VARIABLES'
            GOTO 2000
         ENDIF
         NARI(IARI,1) = -ISCR
      ENDIF

      NARI(IARI,2) = ITYP
      NARI(IARI,3) = NANA(IA,2)

!     FINDING LEFT OPERANDS
!     *********************

      IF(CANA(ITYP:ITYP).EQ.'O') THEN
         NL = 1
      ELSEIF(CANA(ITYP:ITYP).EQ.'=') THEN
         NL = NPAR(ILEFT,4)
      ELSEIF(CANA(ITYP:ITYP).EQ.'A'.OR.CANA(ITYP:ITYP).EQ.'F') THEN
         NL = 0
      ELSE
         AER = 'ERROR, UNKNOWN TYPE'
         GOTO 2000
      ENDIF

      JL = NANA(IA,4)
      DO 330 IL=1,NL
      IF(IL.GE.2) THEN
         IF(CANA(NANA(JL,1):NANA(JL,1)).NE.',') THEN
            AER = 'ERROR, EXPECTED COMMA'
            GOTO 1000
         ENDIF
         NANA(JL,1) = INDEX(CANA,'#')
         JL = NANA(JL,4)
      ENDIF
      IF(INDEX('CVS',CANA(NANA(JL,1):NANA(JL,1))).EQ.0) THEN
         AER = 'ERROR, LEFT OPERAND NOT FOUND'
         GOTO 2000
      ENDIF
      IF(CANA(ITYP:ITYP).EQ.'=') THEN
         NARI(IARI,5+NPAR(ILEFT,4)-IL) = NANA(JL,2)
      ELSE
         NARI(IARI,4) = NANA(JL,2)
      ENDIF

      NANA(IA,4) = NANA(JL,4)
      NANA(NANA(JL,4),5) = IA

      NANA(JL,1) = INDEX(CANA,'#')
      NANA(JL,2) = 0
      NANA(JL,3) = 0
      NANA(JL,4) = 0
      NANA(JL,5) = 0

      JL = NANA(IA,4)

 330  CONTINUE

!     FINDING RIGHT OPERANDS
!     **********************

      IF(CANA(ITYP:ITYP).EQ.'O') THEN
         NR = 1
      ELSEIF(CANA(ITYP:ITYP).EQ.'A') THEN
         NR = NPAR(NANA(IA,2),4)
         IF(NR.EQ.0) THEN
            AER = 'ERROR, NR = 0 FOR ARRAY'
            GOTO 2000
         ENDIF
      ELSEIF(CANA(ITYP:ITYP).EQ.'F') THEN
         IF(NANA(IA,2).GT.0) THEN
            NR = NPAR(NANA(IA,2),4)
         ELSE
            NR = KFUN(-NANA(IA,2))
         ENDIF
      ELSEIF(CANA(ITYP:ITYP).EQ.'=') THEN
         NR = 1
      ELSE
         AER = 'ERROR, UNKNOWN TYPE'
         GOTO 2000
      ENDIF

      JR = NANA(IA,5)
      DO 340 IR = 1,NR
      IF(IR.GE.2) THEN
         IF(CANA(NANA(JR,1):NANA(JR,1)).NE.',') THEN
            AER = 'ERROR, EXPECTED COMMA'
            GOTO 1000
         ENDIF
         NANA(JR,1) = INDEX(CANA,'#')
         JR = NANA(JR,5)
      ENDIF
      IF(INDEX('CVS',CANA(NANA(JR,1):NANA(JR,1))).EQ.0) THEN
         AER = 'ERROR, RIGHT OPERAND NOT FOUND'
         GOTO 2000
      ENDIF
      IF(CANA(ITYP:ITYP).EQ.'=') THEN
         NARI(IARI,4) = NANA(JR,2)
      ELSE
         NARI(IARI,4+IR) = NANA(JR,2)
      ENDIF

      NANA(IA,5) = NANA(JR,5)
      NANA(NANA(JR,5),4) = IA

      NANA(JR,1) = INDEX(CANA,'#')
      NANA(JR,2) = 0
      NANA(JR,3) = 0
      NANA(JR,4) = 0
      NANA(JR,5) = 0

      JR = NANA(IA,5)

 340  CONTINUE

!     REPLACING OPERATION JUST PROCESSED BY SCRATCH ADDRESS
!     *****************************************************

      NANA(IA,1) = INDEX(CANA,'S')
      NANA(IA,2) = -ISCR
      NANA(IA,3) = 0

      GOTO 310

 350  CONTINUE

!frs
!frs      IF(LCHECK.EQ.1) THEN
!frs         WRITE(OUTPUT_UNIT,'(79A)') ('-',J=1,79)
!frs         WRITE(OUTPUT_UNIT,'(100(6X,8(A1,I3,I3,A1)/))')
!frs     *   (CANA(NANA(J,1):NANA(J,1)),NANA(J,2),NANA(J,3),'|',J=1,IANA)
!frs      ENDIF
!frs

!     ALL OPERATIONS EXTRACTED
!     ************************

! 400  CONTINUE

      IF(NANA(1,5).NE.NANA(IANA,4)) THEN
         WRITE(ERROR_UNIT,'(1X,2A)') '### ERROR, CODE NOT FULLY PROCESSED ', '(CHECK DIMENSIONS)'
         IF(IU.NE.0) WRITE(IU,'(1X,2A)') '### ERROR, CODE NOT FULLY ', 'PROCESSED (CHECK DIMENSIONS)'
         LERR = 1
      ENDIF

!frs
!frs      IF(LCHECK.GE.1) THEN
!frs         DO 410 I=IARIA,IARI
!frs  410    WRITE(OUTPUT_UNIT,'(5I5)') (NARI(I,J),J=1,5)
!frs      ENDIF
!frs

      RETURN

!     SYNTAX ERROR EXIT
!     *****************

 1000 CONTINUE
      WRITE(OUTPUT_UNIT,'(1X,2A)') ' ### ', AER
      WRITE(OUTPUT_UNIT,'(1X,2A)') ' ', A(1:IA-1)
      WRITE(OUTPUT_UNIT,'(1X,255A)') ' ', (' ',K=1,I-1),'*'
      IF(IU.NE.0) THEN
         WRITE(IU,'(1X,2A)') ' ### ', AER
         WRITE(IU,'(1X,2A)') ' ', A(1:IA-1)
         WRITE(IU,'(1X,255A)') ' ', (' ',K=1,I-1),'*'
      ENDIF
      LERR = 1
      RETURN

 2000 CONTINUE
!frs
!frs      IF(LCHECK.EQ.1) THEN
!frs         WRITE(OUTPUT_UNIT,'(1X,A)') '### INTERNAL ERROR WHILE EXTRACTING CODE:'
!frs         WRITE(OUTPUT_UNIT,'(1X,A)') AER
!frs         WRITE(OUTPUT_UNIT,'(100(6X,8(A1,I3,I3,A1)/))')
!frs     *     (CANA(NANA(J,1):NANA(J,1)),NANA(J,2),NANA(J,3),'|',J=1,IANA)
!frs         WRITE(OUTPUT_UNIT,'(100(6X,8A8/))') ('        ',K=1,IA-1),'*       '
!frs      ENDIF
!frs
      IF(IU.NE.0) THEN
         WRITE(IU,'(1X,A)') '### INTERNAL ERROR WHILE EXTRACTING CODE:'
         WRITE(IU,'(1X,A)') AER
         WRITE(IU,'(100(6X,8(A1,I3,I3,A1)/))') (CANA(NANA(J,1):NANA(J,1)),NANA(J,2),NANA(J,3),'|',J=1,IANA)
         WRITE(IU,'(100(6X,8A8/))') ('        ',K=1,IA-1),'*     '
      ENDIF
      RETURN
      END

! ANFANG UNTERPROGRAMM
  
      SUBROUTINE ARIFOR(IARI1,IARI2,IU)
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : ERROR_UNIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     **********************************
!
!     THIS SUBROUTINE OUTPUTS THE ARITMETIC INSTRUCTIONS CODED IN
!     THE ARRAY CODE BY THE SUBROUTINE SYNTAX AS REGULAR FORTRAN
!
!-----MEMORY MANAGEMENT -----------------------------------------------------! 1
      PARAMETER(LNAM=10000,LTEX=4000,LCC=10000)                              ! 2
      CHARACTER CNAM(LNAM)*8,CTEX(LTEX)*132,CBLA*132                         ! 3
      INTEGER NPAR(LNAM,17)                                                  ! 4
      DOUBLE PRECISION CC(LCC)                                               ! 5
      COMMON / CMEM / CNAM,CTEX,CBLA                                         ! 6
      COMMON /  MEM / NPAR, CC, INAM, ITEX, ICC                              ! 7
!----------------------------------------------------------------------------! 8
!-----CODE ------------------------------------------------------------------! 1
      PARAMETER (LARI=10000)                                                 ! 2
      INTEGER NARI(LARI,11), IARI                                            ! 3
      COMMON / CODE / NARI, IARI                                             ! 4
!----------------------------------------------------------------------------! 5
!-----SYMBOL---------------------------------------------------------------- ! 1
      PARAMETER (LFUNC=100,LOPER=16)                                         ! 2
      INTEGER KFUN(LFUNC)
      CHARACTER OPER(LOPER)*6, FUNC(LFUNC)*6                                 ! 3
      COMMON / SYMBOL1 / KFUN
      COMMON / SYMBOL / OPER, FUNC                                           ! 4
!----------------------------------------------------------------------------! 5
!
      PARAMETER (LSCR=99)
      CHARACTER A*800,AL*808,AR*800,AS*800,AC*800,BLANK*800!,ABC*26
      INTEGER ISCRTY(LSCR),INDEX(7)

      INTEGER ITS

      SAVE ICALL
      SAVE BLANK

      DATA ICALL / 0 /
!      DATA ABC / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
      IF(ICALL.EQ.0) THEN
         ICALL = 1
         DO I=1,800
           BLANK(I:I) = ' '
         END DO
      ENDIF

      ITS = -4242

      ISOFF = 0
      ISOUT = 1

      DO 900 I=IARI1,IARI2
!     ********************

      IS  = NARI(I,1)
      ITY = NARI(I,2)

!     DETERMINE AL AND AR AND TYPE OF STORAGE ADDRESS
!     -----------------------------------------------

      IF(ITY.EQ.1) THEN
         IR  = NARI(I,4)
         CALL VNAM(IR, INDEX, AR,800,LR,ITR,ISCRTY)
         ITS = ITR
      ELSEIF(ITY.EQ.2) THEN
         IL  = NARI(I,4)
         CALL VNAM(IL, INDEX, AL,800,LL,ITL,ISCRTY)
         IR  = NARI(I,5)
         CALL VNAM(IR, INDEX, AR,800,LR,ITR,ISCRTY)
         IF(ITR.EQ.1) THEN
            IF((ITL.EQ.1).OR.(ITL.EQ.2)) THEN
               ITS = 1
            ELSEIF(ITL.EQ.4) THEN
               ITS = 4
            ELSE
               WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, TYPES INCOMPATIBLE1'
               WRITE(2,'(1X,A)') '### ERROR, TYPES INCOMPATIBLE1'
            ENDIF
         ELSEIF(ITR.EQ.2) THEN
            IF(ITL.EQ.1) THEN
               ITS = 1
            ELSEIF(ITL.EQ.2) THEN
               ITS = 2
            ELSEIF(ITL.EQ.4) THEN
               ITS = 4
            ELSE
               WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, TYPES INCOMPATIBLE1'
               WRITE(2,'(1X,A)') '### ERROR, TYPES INCOMPATIBLE1'
            ENDIF
         ELSEIF(ITR.EQ.4) THEN
            IF((ITL.EQ.1).OR.(ITL.EQ.2)) THEN
               ITS = 4
            ELSEIF(ITL.EQ.4) THEN
               ITS = 4
            ELSE
               WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, TYPES INCOMPATIBLE2'
               WRITE(2,'(1X,A)') '### ERROR, TYPES INCOMPATIBLE2'
            ENDIF
         ELSE
            WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, TYPES INCOMPATIBLE3'
            WRITE(2,'(1X,A)') '### ERROR, TYPES INCOMPATIBLE3'
         ENDIF
      ELSEIF(ITY.EQ.3) THEN
         IF(NARI(I,3).GT.0) THEN
            IR  = NARI(I,3)
            DO J=1,NPAR(IR,4)
              INDEX(J) = NARI(I,4+J)
            END DO
            CALL VNAM(IR, INDEX, AR,800,LR,ITR,ISCRTY)
            ITS = ITR
         ELSE
            IF(KFUN(-NARI(I,3)).EQ.1) THEN
               IR  = NARI(I,5)
               CALL VNAM(IR, INDEX, AR,800,LR,ITR,ISCRTY)
            ELSE
               IL = NARI(I,5)
               CALL VNAM(IL, INDEX, AL,800,LL,ITL,ISCRTY)
               IR = NARI(I,6)
               CALL VNAM(IR, INDEX, AR,800,LR,ITR,ISCRTY)
               IF(ITL.NE.ITR) THEN
                  WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, TYPES INCOMPATIBLE9'
                  WRITE(2,'(1X,A)') '### ERROR, TYPES INCOMPATIBLE9'
               ENDIF
            ENDIF
            ITS = ITR
         ENDIF
      ELSEIF(ITY.EQ.4) THEN
         IR  = NARI(I,3)
         DO J=1,NPAR(IR,4)
           INDEX(J) = NARI(I,4+J)
         END DO
         CALL VNAM(IR, INDEX, AR,800,LR,ITR,ISCRTY)
         ITS = ITR
      ENDIF

!     CHECKING IF STORAGE TYPE IS CONSISTENT AND DETERMINING AS
!     ---------------------------------------------------------

      ISOUT = 0
      ICON = 0

      IF(IS.LT.0) THEN
         ISCRTY(-IS) = ITS
         ISOFF = -IS
         CALL VNAM(IS, INDEX, AS,800,LS,ITSS,ISCRTY)
      ELSEIF(IS.GT.0.AND.IS.LE.LNAM) THEN
         DO J=1,NPAR(IS,4)
           INDEX(J) = NARI(I,4+J)
         END DO
         CALL VNAM(IS, INDEX, AS,800,LS,ITSS,ISCRTY)
         IF(ITS.NE.NPAR(IS,2)) THEN
            IF((ITS.EQ.1).AND.(NPAR(IS,2).EQ.4)) THEN
               ICON = 1
               AC = BLANK
               AC = 'CALL DACON('//AS(1:LS)//','//'RSCRRI(100))'
!               LC = ILAST(AC,1,800)
               AS(1:11) = 'RSCRRI(100)'
               LS = 11
            ELSEIF((ITS.EQ.2).AND.(NPAR(IS,2).EQ.4)) THEN
               ICON = 1
               AC = BLANK
               AC = 'CALL DACON('//AS(1:LS)//',ONE*'//'ISCRRI(100))'
!               LC = ILAST(AC,1,800)
               AS(1:11) = 'ISCRRI(100)'
               LS = 11
            ELSE
               WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, TYPES INCOMPATIBLE4'
               WRITE(2,'(1X,A)') '### ERROR, TYPES INCOMPATIBLE4'
            ENDIF
         ENDIF
      ELSE
         WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, TYPES INCOMPATIBILE5'
         WRITE(2,'(1X,A)') '### ERROR, TYPES INCOMPATIBILE5'
      ENDIF

!     ASSIGNMENT
!     **********

      IF(ITY.EQ.1) THEN

         IF((ITR.EQ.1).OR.(ITR.EQ.2)) THEN
            IF((ITS.EQ.1).OR.(ITS.EQ.2)) THEN
               A = AS(1:LS)//' = '//AR(1:LR)
            ELSEIF(ITS.EQ.4) THEN
               A = 'CALL DACON('//AS(1:LS)//',ONE*'//AR(1:LR)//')'
            ELSE
               GOTO 1000
            ENDIF
         ELSEIF(ITR.EQ.4) THEN
            IF(ITS.EQ.4) THEN
               A = 'CALL DACOP('//AR(1:LR)//','//AS(1:LS)//')'
            ELSE
               GOTO 1000
            ENDIF
         ELSE
            GOTO 1000
         ENDIF

!     BINARY OPERATOR
!     ****************

      ELSEIF(ITY.EQ.2) THEN
         IOP = NARI(I,3)

         IF(IOP.EQ.1) THEN
!        -----------------
            IF((ITL.EQ.1).OR.(ITL.EQ.2)) THEN
               IF((ITR.EQ.1).OR.(ITR.EQ.2)) THEN
                  ITS = 1
                  A = AS(1:LS)//' = '//AL(1:LL)//' + '//AR(1:LR)
               ELSEIF(ITR.EQ.4) THEN
                  ITS = 4
                  A = 'CALL DACAD('//AR(1:LR)//',ONE*'//AL(1:LL)//','// AS(1:LS)//')'
               ELSE
                  GOTO 1000
               ENDIF
            ELSEIF(ITL.EQ.4) THEN
               IF((ITR.EQ.1).OR.(ITR.EQ.2)) THEN
                  ITS = 4
                  A = 'CALL DACAD('//AL(1:LL)//',ONE*'//AR(1:LR)//','// AS(1:LS)//')'
               ELSEIF(ITR.EQ.4) THEN
                  ITS = 4
                  A = 'CALL DAADD('//AL(1:LL)//','//AR(1:LR)//','// AS(1:LS)//')'
               ELSE
                  GOTO 1000
               ENDIF
            ELSE
               GOTO 1000
            ENDIF

         ELSEIF(IOP.EQ.2) THEN
!        ---------------------
            IF((ITL.EQ.1).OR.(ITL.EQ.2)) THEN
               IF((ITR.EQ.1).OR.(ITR.EQ.2)) THEN
                  A = AS(1:LS)//' = '//AL(1:LL)//' - '//AR(1:LR)
               ELSEIF(ITR.EQ.4) THEN
                  A = 'CALL DASUC('//AR(1:LR)//',ONE*'//AL(1:LL)//','// AS(1:LS)//')'
               ELSE
                  GOTO 1000
               ENDIF
            ELSEIF(ITL.EQ.4) THEN
               IF((ITR.EQ.1).OR.(ITR.EQ.2)) THEN
                  A = 'CALL DACSU('//AL(1:LL)//',ONE*'//AR(1:LR)//','// AS(1:LS)//')'
               ELSEIF(ITR.EQ.4) THEN
                  A = 'CALL DASUB('//AL(1:LL)//','//AR(1:LR)//','// AS(1:LS)//')'
               ELSE
                  GOTO 1000
               ENDIF
            ELSE
               GOTO 1000
            ENDIF

         ELSEIF(IOP.EQ.3) THEN
!        ---------------------
            IF((ITL.EQ.1).OR.(ITL.EQ.2)) THEN
               IF((ITR.EQ.1).OR.(ITR.EQ.2)) THEN
                  A = AS(1:LS)//' = '//AL(1:LL)//' * '//AR(1:LR)
               ELSEIF(ITR.EQ.4) THEN
                  A = 'CALL DACMU('//AR(1:LR)//',ONE*'//AL(1:LL)//','// AS(1:LS)//')'
               ELSE
                  GOTO 1000
               ENDIF
            ELSEIF(ITL.EQ.4) THEN
               IF((ITR.EQ.1).OR.(ITR.EQ.2)) THEN
                  A = 'CALL DACMU('//AL(1:LL)//',ONE*'//AR(1:LR)//','// AS(1:LS)//')'
               ELSEIF(ITR.EQ.4) THEN
                  A = 'CALL DAMUL('//AL(1:LL)//','//AR(1:LR)//','// AS(1:LS)//')'
               ELSE
                  GOTO 1000
               ENDIF
            ENDIF

         ELSEIF(IOP.EQ.4) THEN
!        ---------------------
            IF((ITL.EQ.1).OR.(ITL.EQ.2)) THEN
               IF((ITR.EQ.1).OR.(ITR.EQ.2)) THEN
                  A = AS(1:LS)//' = '//AL(1:LL)//' / '//AR(1:LR)
               ELSEIF(ITR.EQ.4) THEN
                  A = 'CALL DADIC('//AR(1:LR)//',ONE*'//AL(1:LL)//','// AS(1:LS)//')'
               ELSE
                  GOTO 1000
               ENDIF
            ELSEIF(ITL.EQ.4) THEN
               IF((ITR.EQ.1).OR.(ITR.EQ.2)) THEN
                  A = 'CALL DACDI('//AL(1:LL)//',ONE*'//AR(1:LR)//','// AS(1:LS)//')'
               ELSEIF(ITR.EQ.4) THEN
                  A = 'CALL DADIV('//AL(1:LL)//','//AR(1:LR)//','// AS(1:LS)//')'
               ELSE
                  GOTO 1000
               ENDIF
            ELSE
               GOTO 1000
            ENDIF
!
         ELSEIF(IOP.EQ.5) THEN
!        ---------------------
            IF((ITL.EQ.1).OR.(ITL.EQ.2)) THEN
               IF((ITR.EQ.1).OR.(ITR.EQ.2)) THEN
                  A = AS(1:LS)//' = '//AL(1:LL)//'** '//AR(1:LR)
               ELSEIF(ITR.EQ.4) THEN
                  A = 'CALL DACEX('//AR(1:LR)//',ONE*'//AL(1:LL)//',' //AS(1:LS)//')'
               ELSE
                  GOTO 1000
               ENDIF
            ELSEIF(ITL.EQ.4) THEN
               IF((ITR.EQ.1).OR.(ITR.EQ.2)) THEN
                  A = 'CALL DAEXC('//AL(1:LL)//',ONE*'//AR(1:LR)//',' //AS(1:LS)//')'
               ELSEIF(ITR.EQ.4) THEN
                  A = 'CALL DAEXX('//AL(1:LL)//','//AR(1:LR)//',' //AS(1:LS)//')'
               ELSE
                  GOTO 1000
               ENDIF
            ELSE
               GOTO 1000
            ENDIF
         ELSE
!        ----
            GOTO 1000
         ENDIF

!     FUNCTIONS
!     *********

      ELSEIF(ITY.EQ.3) THEN
         IOP = NARI(I,3)
         IF(IOP.GT.0) THEN
            WRITE(2,'(6X,A,I5)') 'IDAO = IDAA + ',ISOFF
            ISOUT = 1
         ENDIF

         IF((ITR.EQ.1).OR.(ITR.EQ.2)) THEN
!        ---------------------------------

            IF(IOP.GT.0) THEN
               A = AS(1:LS)//' = '//AR(1:LR)
            ELSE
               IF(KFUN(-IOP).EQ.1) THEN
                  A = AS(1:LS)//' = '//FUNC(-IOP)(1:4)// '('//AR(1:LR)//')'
               ELSE
                  A = AS(1:LS)//' = '//FUNC(-IOP)(1:4)// '('//AL(1:LL)//','//AR(1:LR)//')'
               ENDIF
            ENDIF

         ELSEIF(ITR.EQ.4) THEN
!        ---------------------

            IF(IOP.GT.0) THEN
               A = 'CALL DACOP('//AR(1:LR)//','//AS(1:LS)//')'
            ELSE
               IF(KFUN(-IOP).EQ.1) THEN
                  A = 'CALL DAFUN('//''''//FUNC(-IOP)(1:4)//''''//','// AR(1:LR)//','//AS(1:LS)//')'
               ELSE
                  A = 'CALL DAFUN2('//''''//FUNC(-IOP)(1:4)//''''//','// AL(1:LL)//','//AR(1:LR)//','//AS(1:LS)//')'
               ENDIF
            ENDIF

         ELSE
            GOTO 1000
         ENDIF

!     ARRAYS
!     ******

      if (ITS.eq.-4242) then
         WRITE(ERROR_UNIT,*) "ERROR in ARIFOR - ITS not set"
         stop 16
      endif
      
      ELSEIF(ITY.EQ.4) THEN
         IF((ITS.EQ.1).OR.(ITS.EQ.2)) THEN
            A = AS(1:LS)//' = '//AR(1:LR)
         ELSEIF(ITS.EQ.4) THEN
            A = 'CALL DACOP('//AR(1:LR)//','//AS(1:LS)//')'
         ELSE
            GOTO 1000
         ENDIF
      ENDIF

!     OUTPUTTING FORTRAN STATEMENT
!     ****************************

      ILIN = ILAST(A,1,800) / 125
      WRITE(2,'(A)') '      '//A(1:125), ('     !'//A(J*125+1:(J+1)*125),J=1,ILIN)
! Eric Moved the ILIN = ....
!     ILIN = ILAST(AC,1,800) / 66
      IF(ICON.NE.0) THEN
        ILIN = ILAST(AC,1,800) / 125
         WRITE(2,'(A)') '      '//AC(1:125), ('     !'//AC(J*125+1:(J+1)*125),J=1,ILIN)
         ICON = 0
      ELSEIF(ISOUT.NE.0) THEN
         WRITE(2,'(6X,A,I5)') 'IDAO = IDAA - ',ISOFF
         ISOUT = 0
      ENDIF


  900 CONTINUE

      RETURN

 1000 CONTINUE
      WRITE(IU,'(1X,A)') '!!! ERROR IN ARIFOR, CODE NOT SUPPORTED'
      WRITE(ERROR_UNIT,'(1X,A)') '!!! ERROR IN ARIFOR, CODE NOT SUPPORTED:'
      WRITE(ERROR_UNIT,'(1X,5I5)') (NARI(I,J),J=1,5)
  
      RETURN
      END

! ANFANG UNTERPROGRAMM
  
      SUBROUTINE VNAM(IA, IN, A,NA,LA, ITA, ISCRTY)
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : ERROR_UNIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     *********************************************
!
!     THIS SUBROUTINE DETERMINES THE PROPER CHARACTER NAME A BELONGING
!     TO AN ADDRESS IA. IN ARE THE UP TO SEVEN INDICES (FOR ARRAYS).
!     ITA IS THE DATA TYPE OF THE ADDRESS, ISCRTY CONTAINS THE DATA TYPES
!     OF THE SCRATCH VARIABLES.
!     THE SUBROUTINE IS CALLED FROM ARIFOR.
!
!-----MEMORY MANAGEMENT -----------------------------------------------------! 1
      PARAMETER(LNAM=10000,LTEX=4000,LCC=10000)                              ! 2
      CHARACTER CNAM(LNAM)*8,CTEX(LTEX)*132,CBLA*132                         ! 3
      INTEGER NPAR(LNAM,17)                                                  ! 4
      DOUBLE PRECISION CC(LCC)                                               ! 5
      COMMON / CMEM / CNAM,CTEX,CBLA                                         ! 6
      COMMON /  MEM / NPAR, CC, INAM, ITEX, ICC                              ! 7
!----------------------------------------------------------------------------! 8
!-----CODE ------------------------------------------------------------------! 1
      PARAMETER (LARI=10000)                                                 ! 2
      INTEGER NARI(LARI,11), IARI                                            ! 3
      COMMON / CODE / NARI, IARI                                             ! 4
!----------------------------------------------------------------------------! 5

      PARAMETER (LSCR=99)
      CHARACTER A*(*),CNUM(99)*2
      INTEGER ISCRTY(*),IN(*)

      INTEGER ND

      DATA CNUM / ' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9','10', &
                  '11','12','13','14','15','16','17','18','19','20', &
                  '21','22','23','24','25','26','27','28','29','30', &
                  '31','32','33','34','35','36','37','38','39','40', &
                  '41','42','43','44','45','46','47','48','49','50', &
                  '51','52','53','54','55','56','57','58','59','60', &
                  '61','62','63','64','65','66','67','68','69','70', &
                  '71','72','73','74','75','76','77','78','79','80', &
                  '81','82','83','84','85','86','87','88','89','90', &
                  '91','92','93','94','95','96','97','98','99'/

      ND = -4242
      
      IF(IA.LT.0) THEN
         IA = -IA
         ITA = ISCRTY(IA)
         ND = 0
         IF(IA.LT.1.OR.IA.GT.100.OR.IA.GT.LSCR) THEN
            WRITE(ERROR_UNIT,'(1X,A,I4)') '!!! ERROR IN VNAM, IA = ',IA
            WRITE(2,'(1X,A,I4)') '!!! ERROR IN VNAM, IA = ',IA
            STOP 11
         ENDIF
         IF(ITA.EQ.1) THEN
            LA = 16
            A(1:LA) = 'RSCRRI( '//CNUM(IA)//'+IDAA)'
         ELSEIF(ITA.EQ.2) THEN
            LA = 16
            A(1:LA) = 'ISCRRI( '//CNUM(IA)//'+IDAA)'
         ELSE
            LA = 16
            A(1:LA) = 'ISCRDA( '//CNUM(IA)//'+IDAA)'
         ENDIF

      ELSEIF((IA.GE.1).AND.(IA.LE.LNAM)) THEN
         ITA = NPAR(IA,2)
         LA = 11
         A(1:LA) = CNAM(IA)//'   '
         ND = NPAR(IA,4)
      ELSEIF(IA.GT.LNAM) THEN
         IA = IA - LNAM
         LA = MAX(11,ILAST(CTEX(IA),1,132))
         A = '('//CTEX(IA)(1:LA)//')'//CBLA
         LA = LA + 2
         IF(INDEX(A(1:LA),'.').EQ.0) THEN
            ITA = 2
         ELSE
            ITA = 1
         ENDIF
         ND = 0
      ELSE
         WRITE(ERROR_UNIT,'(1X,A)') '!!! ERROR, IA = 0 IN VNAM'
         WRITE(2,'(1X,A)') '!!! ERROR, IA = 0 IN VNAM'
      ENDIF
      if (ND .eq. -4242) then
         WRITE(ERROR_UNIT,*) "!!! ERROR, ND was not initialized in VNAM"
         stop 15
      end if
      

!     CASE OF ARRAYS OR FUNCTIONS
!     ***************************

      IF(ND.EQ.0) RETURN

!      IF(NPAR(IA,1).EQ.1) THEN
!         NT = NPAR(IA,5)
!         A(LA+1:LA+11) = '(NINT(ONE*'
!         LA = LA + 11
!      ELSEIF(NPAR(IA,1).EQ.2) THEN
         A(LA+1:LA+1) = '('
         LA = LA + 1
!      ENDIF

      DO 100 ID=1,ND

      IF(ID.NE.1) THEN
!         IF(NPAR(IA,1).EQ.1) THEN
!            LT = ILAST(CTEX(NT+ID),1,80)
!            A(LA+1:LA+5+LT) = ' + '//CTEX(NT+ID)(1:LT)//'*('
!            LA = LA + 5 + LT
!         ELSEIF(NPAR(IA,1).EQ.2) THEN
            A(LA+1:LA+1) = ','
            LA = LA + 1
!         ENDIF
      ENDIF

      II = IN(ID)
      IF(II.LT.0) THEN
         II = -II
         ITI = ISCRTY(II)
         IF(II.LT.1.OR.II.GT.40.OR.II.GT.LSCR) THEN
            WRITE(ERROR_UNIT,'(1X,A,I4)') '!!! ERROR IN VNAM, II = ',II
            WRITE(2,'(1X,A,I4)') '!!! ERROR IN VNAM, II = ',II
            STOP 12
         ENDIF
         IF(ITI.EQ.1) THEN
            A(LA+1:LA+16) = 'RSCRRI( '//CNUM(II)//'+IDAA)'
            LA = LA  + 16
         ELSEIF(ITI.EQ.2) THEN
            A(LA+1:LA+16) = 'ISCRRI( '//CNUM(II)//'+IDAA)'
            LA = LA  + 16
         ELSE
            A(LA+1:LA+16) = 'ISCRDA( '//CNUM(II)//'+IDAA)'
            LA = LA + 16
         ENDIF

      ELSEIF(II.LE.LNAM.AND.II.GE.1) THEN
         A(LA+1:LA+11) = CNAM(II)//'   '
         LA = LA + 11
         ITI = NPAR(II,2)
      ELSEIF(II.GT.LNAM) THEN
         II = II - LNAM
         IL = ILAST(CTEX(II),1,132)
         IF(INDEX(CTEX(II)(1:IL),'.').EQ.0) THEN
            ITI = 2
         ELSE
            ITI = 1
         ENDIF
         A(LA+1:LA+IL+2) = '('//CTEX(II)(1:IL)//')'
         IL = IL + 2
         LA = LA + IL
      ELSE
         WRITE(ERROR_UNIT,'(1X,A)') '!!! ERROR IN VNAM, II = 0'
         WRITE(2,'(1X,A)') '!!! ERROR IN VNAM, II = 0'
         STOP 13
      ENDIF
      IF((NPAR(IA,1).EQ.1).AND.(ITI.NE.1).AND.(ITI.NE.2)) THEN
         WRITE(ERROR_UNIT,'(1X,A)') '### ERROR, ARRAY INDEX NOT INTEGER OR REAL'
         WRITE(2,'(1X,A)') '### ERROR, ARRAY INDEX NOT INTEGER OR REAL'
         RETURN
      ENDIF

!      IF(ID.NE.1.AND.NPAR(IA,1).EQ.1) THEN
!         A(LA+1:LA+2) = '-1'
!         LA = LA + 2
!      ENDIF

 100  CONTINUE

!      IF(NPAR(IA,1).EQ.1) THEN
!         DO 110 ID=2,ND+2
!         A(LA+1:LA+1) = ')'
!         LA = LA + 1
! 110     CONTINUE
!      ELSEIF(NPAR(IA,1).EQ.2) THEN
         A(LA+1:LA+1) = ')'
         LA = LA + 1
!      ENDIF

      IF(LA.GT.NA) THEN
         WRITE(ERROR_UNIT,'(1X,A)') '!!! ERROR IN ROUTINE VNAM, IA > NA'
         WRITE(2,'(1X,A)') '!!! ERROR IN ROUTINE VNAM, IA > NA'
         STOP 14
      ENDIF

      RETURN
      END
