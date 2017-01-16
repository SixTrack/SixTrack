+cd crlibco
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
+cd crcoall
      integer lout
      common /crflags/lout

+dk plato_seq
*
* $Id: cfft.F,v 1.1.1.1 1996/02/15 17:48:48 mclareni Exp $
*
* $Log: cfft.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:48  mclareni
* Kernlib
* Conversion single->double precission by K.Sjobak, Dec. 2016
*
      SUBROUTINE CFFT(A,MSIGN)
      IMPLICIT NONE
      COMPLEX*16 A(*),U,W,T
      INTEGER MSIGN, M, N, NV2, NM1, J, I, K, L, LE, LE1, IP
      double precision C,S
      IF(MSIGN.EQ.0) RETURN
      M=IABS(MSIGN)
      N=2**M
      NV2=N/2
      NM1=N-1
      J=1
      DO 7 I=1,NM1
      IF(I.GE.J) GO TO 5
      T=A(J)
      A(J)=A(I)
      A(I)=T
 5    K=NV2
 6    IF(K.GE.J) GO TO 7
      J=J-K
      K=K/2
      GO TO 6
 7    J=J+K
      DO 8 I=1,N,2
      T=A(I+1)
      A(I+1)=A(I)-T
 8    A(I )=A(I)+T
      IF(M.EQ.1) RETURN
      C=0.
      S=ISIGN(1,MSIGN)
      LE=2
      DO 20 L=2,M
      W=CMPLX(C,S,kind(1.d0))
      U=W
      C=SQRT(C*.5+.5)
      S=DIMAG(W)/(C+C)
      LE1=LE
      LE=LE1+LE1
      DO 9 I=1,N,LE
      IP=I+LE1
      T=A(IP)
      A(IP)=A(I)-T
 9    A(I) =A(I)+T
      DO 20 J=2,LE1
      DO 10 I=J,N,LE
      IP=I+LE1
      T=A(IP)*U
      A(IP)=A(I)-T
 10   A(I) =A(I)+T
 20   U=U*W
      RETURN
      END

CDECK  ID>, TUNENEWT1.
C=============================================================
C COMPUTES THE TUNE USING A DISCRETE VERSION OF LASKAR METHOD.
C IT INCLUDES A NEWTON METHOD FOR THE SEARCH OF THE FREQUENCY.
C X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE 
C LENGTH IF THE ORBIT.
C
C AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
C             M. GIOVANNOZZI - CERN HAS INTRODUCED SOME 
C                                   MODIFICATIONS
C

      DOUBLE PRECISION FUNCTION TUNENEWT1(X,XP,MAXN)
      
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      DOUBLE PRECISION X,XP, DUEPI, STEP, SUM, FTMAX, TUNEFOU, DELTAT,
     &     TUNE1, TUNE
      INTEGER MAXN, MFT, NPOINT, MAXN2, MF, NPMIN, NPMAX, NFTMAX, NFT
      COMPLEX*16 Z
      COMPLEX*16 ZSING(MAXITER) ! Temp Z for CFFT, used to be SINGle precission
      
      DIMENSION X(MAXITER),XP(MAXITER)
      DIMENSION Z(MAXITER)
      
+ca crlibco
+if cr
+ca crcoall
+ei
C.............................................................
      IF (MAXN.GT.MAXITER) THEN
+if cr
         WRITE(lout,'(a)') '***ERROR(TUNENEWT1): TOO MANY ITERATIONS'
+ei
+if .not.cr
         WRITE(*,   '(a)') '***ERROR(TUNENEWT1): TOO MANY ITERATIONS'
+ei
         call prror(-1)
      ENDIF
C.............................................................
C    ESTIMATION OF TUNE WITH FFT 
C.............................................................
+if crlibm
      DUEPI=ATAN_RN(1D0)*8D0
      MFT=INT(LOG_RN(DBLE(MAXN))/LOG_RN(2D0))
+ei
+if .not.crlibm
      DUEPI=ATAN(1D0)*8D0
      MFT=INT(LOG(DBLE(MAXN))/LOG(2D0)) 
+ei
      NPOINT=2**MFT
      MAXN2=MAXN/2
      STEP=DUEPI/MAXN
C.............................................................
      SUM=0D0
      DO MF=1,MAXN
        Z(MF)=DCMPLX(X(MF),XP(MF)) !Returns COMPLEX*16 / COMPLEX(8)
        ZSING(MF)=Z(MF)
        SUM=SUM+XP(MF)
      ENDDO 
      CALL CFFT(ZSING,-MFT)
C.........COMPUTATION OF THE MAXIMUM OF THE FOURIER TRANSFORM
      NPMIN=1 
      IF (SUM.EQ.0D0) THEN
        NPMAX=NPOINT/2  !..REAL FFT ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
C............................................................
      FTMAX=0D0
      NFTMAX=0
C.............................REJECTS FREQUENCIES NEAR 0 OR 1
      DO NFT=NPMIN,NPMAX
        IF (ABS(ZSING(NFT)).GT.FTMAX) THEN
          FTMAX=ABS(ZSING(NFT))
          NFTMAX=NFT
        END IF
      ENDDO
      TUNEFOU=DBLE(NFTMAX-1)/DBLE(NPOINT)
      DELTAT=1D0/NPOINT
      TUNE1=TUNEFOU-DELTAT
      CALL ZFUN(TUNE,Z,MAXN,TUNE1,DELTAT)
      TUNENEWT1=1D0-TUNE
C............................................................  
      RETURN 
C............................................................  
      END

CDECK  ID>, TUNEFIT.
C============================================================
C COMPUTES THE TUNE USING A MODIFIED APA ALGORITHM. THE FIRST
C STEP CONSISTS OF TAKING THE AVERAGE OF TUNE COMPUTED WITH THE
C APA METHOD, THEN A BEST FIT IS PERFORMED.
C
C     AUTHORS: R. BARTOLINI A. BAZZANI - BOLOGNA UNIVERSITY
C

      DOUBLE PRECISION FUNCTION TUNEFIT(X,XP,MAX)
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER MAX,MAX1,MAX2,N
      DOUBLE PRECISION X,XP,DUEPI,C,SUMAPM,PA,R1,R2,CTHETA,STHETA,THETA,&
     &TUNEAPA,U,ATUNE,TUNEF,TUNE
      DIMENSION X(*),XP(*),TUNE(MAXITER),U(MAXITER) 

+ca crlibco
+if cr
+ca crcoall
+ei
C............................................................
      IF (MAX.LE.0) THEN
+if cr
         WRITE(lout,'(a)')
+ei
+if .not.cr
         WRITE(*,   '(a)')
+ei
     &        '***ERROR(TUNEFIT): THIRD PARAMETER OUT OF BOUNDS'
        call prror(-1)
      ENDIF
C............................................................
+if crlibm
      DUEPI=8*ATAN_RN(1D0)
+ei
+if .not.crlibm
      DUEPI=8*DATAN(1D0)
+ei
C.............................................................
      C=DBLE(MAX)/2D0-DBLE(INT(MAX/2D0))
      MAX1=MAX
      IF(C.GT.1D-1) MAX1=MAX-1
      MAX2=MAX1/2
C.............................................................
      SUMAPM=0D0
      PA=0D0
C.............................................................
      DO N=1,MAX-1
C.............................................................
        R1=DSQRT(X(N)**2+XP(N)**2)
        R2=DSQRT(X(N+1)**2+XP(N+1)**2)
        CTHETA=((X(N)*X(N+1)+XP(N)*XP(N+1))/R1)/R2
        STHETA=((-X(N)*XP(N+1)+X(N+1)*XP(N))/R1)/R2
+if crlibm
        IF (STHETA.GE.0) THEN
          THETA=ACOS_RN(CTHETA)
          ELSE
          THETA=DUEPI-ACOS_RN(CTHETA) 
        END IF
+ei
+if .not.crlibm
        IF (STHETA.GE.0) THEN
          THETA=ACOS(CTHETA)
          ELSE
          THETA=DUEPI-ACOS(CTHETA) 
        END IF
+ei
        PA=PA+THETA
        TUNEAPA=PA/DUEPI
        SUMAPM=SUMAPM+TUNEAPA
        IF (N.GE.MAX2) THEN
          TUNE(N-MAX2+1)=((SUMAPM/DBLE(N))/DBLE(N+1))*2.D0
          U(N-MAX2+1)=1.D0/N
        ENDIF
C.............................................................
      ENDDO
      CALL FIT(U,TUNE,MAX1-MAX2,ATUNE,TUNEF)
      TUNEFIT=TUNEF
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, FIT.
C============================================================
C COMPUTES THE STRAIGHT LINE WHICH FITS A SET OF N 2D 
C OBSERVATIONS.
C THE LINE IS REPRESENTED BY
C
C Y=AX+B
C
C     AUTHORS: R. BARTOLINI A. BAZZANI - BOLOGNA UNIVERSITY
C

      SUBROUTINE FIT(X,Y,N,A,B)
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION X2M,XM,YM,XYM,ERR,DELTA,A,B,X,Y
      DIMENSION X(*),Y(*)
+if cr
+ca crcoall
+ei
C............................................................
      IF (N.LE.0) THEN
+if cr
        WRITE(lout,'(a)')
+ei
+if .not.cr
        WRITE(*,   '(a)')
+ei
     &       '***ERROR(TUNEAPA): THIRD PARAMETER OUT OF BOUNDS'
        call prror(-1)
      ENDIF
C.............................................................
      X2M=0.D0
      XM=0.D0
      YM=0.D0
      XYM=0.D0
      ERR=0.D0
C.............................................................
      DO I=1,N
        X2M=X2M+X(I)*X(I)
        XYM=XYM+X(I)*Y(I)
        XM=XM+X(I)        
        YM=YM+Y(I)        
      END DO
C.............................................................
      X2M=X2M/N
      XYM=XYM/N
      XM=XM/N
      YM=YM/N
      DELTA=X2M-XM*XM
      A=(XYM-XM*YM)/DELTA
      B=(-XM*XYM+X2M*YM)/DELTA
C.............................................................
      RETURN
C.............................................................
      END

CDECK  ID>, TUNEABT.
C=============================================================
C COMPUTES THE TUNE USING FFT INTERPOLATED METHOD.
C X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE 
C LENGTH OF THE ORBIT.
C
C AUTHOR:     E. TODESCO - INFN AND CERN 
C

      DOUBLE PRECISION FUNCTION TUNEABT(X,XP,MAXN)              
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER MAXN,NPOINT,MF,NPMIN,NPMAX,NFTMAX,NFT
      DOUBLE PRECISION X,XP,PI,DUEPI,MFT,STEP,SUM,FTMAX,CF1,CF2,CF3,ASSK
      COMPLEX*16 Z,ZSING
      DIMENSION X(MAXITER),XP(MAXITER)
      DIMENSION Z(MAXITER),ZSING(MAXITER)

+ca crlibco

C..................................ESTIMATION OF TUNE WITH FFT 
+if crlibm
      PI=ATAN_RN(1D0)*4D0
      MFT=INT(LOG_RN(DBLE(MAXN))/LOG_RN(2D0)) 
+ei
+if .not.crlibm
      PI=ATAN(1D0)*4D0
      MFT=INT(LOG(DBLE(MAXN))/LOG(2D0)) 
+ei
      DUEPI=2*PI
      NPOINT=2**MFT
      STEP=DUEPI/NPOINT/2D+0
C.............................................................
      SUM=0D0            !..CHECKS FOR COMPLEX OR REAL DATA
      DO MF=1,NPOINT
        Z(MF)=DCMPLX(X(MF),XP(MF))
        ZSING(MF)=Z(MF)
        SUM=SUM+XP(MF)
      ENDDO
      CALL FFT_PLATO(ZSING,NPOINT,-1)
C.......................SEARCH FOR MAXIMUM OF FOURIER SPECTRUM
      NPMIN=1
      IF (SUM.EQ.0D0) THEN
        NPMAX=NPOINT/2         !..REAL FFT_PLATO ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
C.............................................................
      FTMAX=0D0
      NFTMAX=0
      DO NFT=NPMIN,NPMAX
        IF (ABS(ZSING(NFT)).GT.FTMAX) THEN
          FTMAX=ABS(ZSING(NFT))
          NFTMAX=NFT
        END IF
      ENDDO
C................................................INTERPOLATION
      CF1=ABS(ZSING(NFTMAX-1))
      CF2=ABS(ZSING(NFTMAX))
      CF3=ABS(ZSING(NFTMAX+1))
+if crlibm
      IF (CF3.GT.CF1) THEN      
        ASSK=DBLE(NFTMAX)+(NPOINT/PI)*
     .       ATAN2_RN(CF3*SIN_RN(PI/NPOINT),CF2+CF3*COS_RN(PI/NPOINT))
      ELSEIF (CF3.LE.CF1) THEN                   
        ASSK=DBLE(NFTMAX-1)+(NPOINT/PI)*
     .       ATAN2_RN(CF2*SIN_RN(PI/NPOINT),CF1+CF2*COS_RN(PI/NPOINT))
      ENDIF
+ei
+if .not.crlibm
      IF (CF3.GT.CF1) THEN      
        ASSK=DBLE(NFTMAX)+(NPOINT/PI)*
     .       ATAN2(CF3*SIN(PI/NPOINT),CF2+CF3*COS(PI/NPOINT))
      ELSEIF (CF3.LE.CF1) THEN                   
        ASSK=DBLE(NFTMAX-1)+(NPOINT/PI)*
     .       ATAN2(CF2*SIN(PI/NPOINT),CF1+CF2*COS(PI/NPOINT))
      ENDIF
+ei
      TUNEABT=1D+0-(ASSK-1D+0)/DBLE(NPOINT) !1D+0 = 1D0, i.e. double precision 1.0?
C............................................................  
      RETURN 
C............................................................  
      END
CDECK  ID>, TUNEABT2.
C=============================================================
C COMPUTES THE TUNE USING THE INTERPOLATED FFT METHOD
C WITH HANNING FILTER.
C X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE 
C LENGTH OF THE ORBIT.
C
C AUTHOR:     E. TODESCO - INFN AND CERN 
C

      DOUBLE PRECISION FUNCTION TUNEABT2(X,XP,MAXN)              
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER MAXN,NPOINT,MFT,MF,NPMIN,NPMAX,NFTMAX,NN,NFT
      DOUBLE PRECISION X,XP,PI,DUEPI,STEP,SUM,FTMAX,CF1,CF2,CF3,P1,P2,
     &CO,SI,SCRA1,SCRA2,SCRA3,SCRA4,ASSK
      COMPLEX*16 Z
      COMPLEX*16 ZSING(MAXITER)
      DIMENSION X(MAXITER),XP(MAXITER)
      DIMENSION Z(MAXITER)

+ca crlibco

C..................................ESTIMATION OF TUNE WITH FFT 
+if crlibm
      PI=ATAN_RN(1D0)*4D0
      MFT=INT(LOG_RN(DBLE(MAXN))/LOG_RN(2D0)) 
+ei
+if .not.crlibm
      PI=ATAN(1D0)*4D0
      MFT=INT(LOG(DBLE(MAXN))/LOG(2D0)) 
+ei
      DUEPI=2*PI
      NPOINT=2**MFT
      STEP=(DUEPI/NPOINT)/2D+0
C.............................................................
      SUM=0D0            !..CHECKS FOR COMPLEX OR REAL DATA
      DO MF=1,NPOINT
+if crlibm
        Z(MF)=DCMPLX(X(MF),XP(MF))*SIN_RN(STEP*MF)**2
+ei
+if .not.crlibm
        Z(MF)=DCMPLX(X(MF),XP(MF))*SIN(STEP*MF)**2
+ei
        ZSING(MF)=Z(MF)
        SUM=SUM+XP(MF)
      ENDDO 
      CALL FFT_PLATO(ZSING,NPOINT,-1)
C.......................SEARCH FOR MAXIMUM OF FOURIER SPECTRUM
      NPMIN=1
      IF (SUM.EQ.0D0) THEN
        NPMAX=NPOINT/2         !..REAL FFT_PLATO ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
C.............................................................
      FTMAX=0D0
      NFTMAX=0
      DO NFT=NPMIN,NPMAX
        IF (ABS(ZSING(NFT)).GT.FTMAX) THEN
          FTMAX=ABS(ZSING(NFT))
          NFTMAX=NFT
        END IF
      ENDDO
      CF1=ABS(ZSING(NFTMAX-1))
      CF2=ABS(ZSING(NFTMAX))
      CF3=ABS(ZSING(NFTMAX+1))
      IF (CF3.GT.CF1) THEN      
        P1=CF2
        P2=CF3
        NN=NFTMAX
      ELSEIF (CF3.LE.CF1) THEN                   
        P1=CF1
        P2=CF2
        NN=NFTMAX-1
      ENDIF
C..........................................INTERPOLATION
+if crlibm
      CO=COS_RN((2*PI)/DBLE(NPOINT))
      SI=SIN_RN((2*PI)/DBLE(NPOINT))
+ei
+if .not.crlibm
      CO=COS((2*PI)/DBLE(NPOINT))
      SI=SIN((2*PI)/DBLE(NPOINT))
+ei
      SCRA1=CO**2*(P1+P2)**2-((2*P1)*P2)*((2*CO**2-CO)-1)
      SCRA2=(P1+P2*CO)*(P1-P2)
      SCRA3=(P1**2+P2**2)+((2*P1)*P2)*CO
      SCRA4=(-SCRA2+P2*SQRT(SCRA1))/SCRA3
+if crlibm
      ASSK=DBLE(NN)+((NPOINT/2)/PI)*ASIN_RN(SI*SCRA4)
+ei
+if .not.crlibm
      ASSK=DBLE(NN)+((NPOINT/2)/PI)*ASIN(SI*SCRA4)
+ei
      TUNEABT2=1D+0-(ASSK-1D+0)/DBLE(NPOINT)
C............................................................  
      RETURN 
C............................................................  
      END

CDECK  ID>, FFT_PLATO.
C============================================================
C           COMPUTES THE FFT_PLATO   (DOUBLE PRECISION)
C           AUTHOR: NUMERICAL RECEIPES, PG. 395
C           NN IS THE NUMBER OF DATA: MUST BE A POWER OF 2
C           ISIGN=1: DIRECT FT
C           ISIGN=-1: INVERSE FT
C           DATA IS A COMPLEX ARRAY WITH THE SIGNAL IN INPUT
C                                   WITH THE FT IN OUTPUT
C           NOTE THAT REAL*8 DATA(*) THAT CONTAINS THE COMPLEX
C             NUMBER IN THE FOLLOWING ORDER:
C             REAL(DATA(1)),IMAG(DATA(1)),REAL(DATA(2)),IMAG(DATA(2)),...
C           GOOD LUCK, BABY
C

      SUBROUTINE FFT_PLATO(CDATA,NN,ISIGN)
      IMPLICIT NONE
      INTEGER I,J,N,NN,M,MMAX,ISTEP,ISIGN
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,TEMPR,TEMPI
      COMPLEX*16 CDATA(NN)
      REAL*8 DATA(2*NN)

+ca crlibco

      N=2*NN
C create real array DATA out of complex array CDATA
      DO I=1,N,2
        DATA(I)=DREAL(CDATA(I/2+1))
        DATA(I+1)=DIMAG(CDATA(I/2+1))
      ENDDO
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I) THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        END IF
        M=N/2
 1      IF((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GOTO 1
        END IF
        J=J+M
 11     CONTINUE
      MMAX=2
 2    IF(N.GT.MMAX) THEN
        ISTEP=2*MMAX
+if crlibm
        THETA=(8.D0*ATAN_RN(1.D0))/(ISIGN*MMAX)
        WPR=-2.D0*SIN_RN(0.5D0*THETA)**2
        WPI=SIN_RN(THETA)
+ei
+if .not.crlibm
        THETA=(8.D0*ATAN(1.D0))/(ISIGN*MMAX)
        WPR=-2.D0*SIN(0.5D0*THETA)**2
        WPI=SIN(THETA)
+ei
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=WR*DATA(J)-WI*DATA(J+1)
            TEMPI=WR*DATA(J+1)+WI*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
 12       CONTINUE
          WTEMP=WR
          WR=(WR*WPR-WI*WPI)+WR
          WI=(WI*WPR+WTEMP*WPI)+WI
 13     CONTINUE
        MMAX=ISTEP
        GOTO 2
      END IF
      DO I=1,N,2
        CDATA(I/2+1)=CMPLX(DATA(I),DATA(I+1),kind(1.d0))
      ENDDO
C.............................................................      
      END

CDECK  ID>, FFT_PLATO_REAL.
C============================================================
C           COMPUTES THE FFT_PLATO_REAL   (DOUBLE PRECISION)
C           AUTHOR: NUMERICAL RECEIPES, PG. 395
C           NN IS THE NUMBER OF DATA: MUST BE A POWER OF 2
C           ISIGN=1: DIRECT FT
C           ISIGN=-1: INVERSE FT
C           DATA IS A COMPLEX ARRAY WITH THE SIGNAL IN INPUT
C                                   WITH THE FT IN OUTPUT
C           NOTE THAT REAL*8 DATA(*) THAT CONTAINS THE COMPLEX
C             NUMBER IN THE FOLLOWING ORDER:
C             REAL(DATA(1)),IMAG(DATA(1)),REAL(DATA(2)),IMAG(DATA(2)),...
C           GOOD LUCK, BABY
C

      SUBROUTINE FFT_PLATO_REAL(DATA,NN,ISIGN)
      IMPLICIT NONE
      INTEGER I,J,N,NN,M,MMAX,ISTEP,ISIGN
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,TEMPR,TEMPI
      REAL*8 DATA(*)

+ca crlibco
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I) THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        END IF
        M=N/2
 1      IF((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GOTO 1
        END IF
        J=J+M
 11     CONTINUE
      MMAX=2
 2    IF(N.GT.MMAX) THEN
        ISTEP=2*MMAX
+if crlibm
        THETA=(8.D0*ATAN_RN(1.D0))/(ISIGN*MMAX)
        WPR=-2.D0*SIN_RN(0.5D0*THETA)**2
        WPI=SIN_RN(THETA)
+ei
+if .not.crlibm
        THETA=(8.D0*ATAN(1.D0))/(ISIGN*MMAX)
        WPR=-2.D0*SIN(0.5D0*THETA)**2
        WPI=SIN(THETA)
+ei
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=WR*DATA(J)-WI*DATA(J+1)
            TEMPI=WR*DATA(J+1)+WI*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
 12       CONTINUE
          WTEMP=WR
          WR=(WR*WPR-WI*WPI)+WR
          WI=(WI*WPR+WTEMP*WPI)+WI
 13     CONTINUE
        MMAX=ISTEP
        GOTO 2
      END IF
C.............................................................      
      END

CDECK  ID>, TUNENEWT.
C=============================================================
C COMPUTES THE TUNE USING A DISCRETE VERSION OF LASKAR METHOD.
C IT INCLUDES A NEWTON METHOD FOR THE SEARCH OF THE FREQUENCY.
C X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE 
C LENGTH IF THE ORBIT.
C
C AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
C             M. GIOVANNOZZI - CERN HAS INTRODUCED SOME 
C                                   MODIFICATIONS
C

      DOUBLE PRECISION FUNCTION TUNENEWT(X,XP,MAXN)
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER MAXN,MFT,NPOINT,MAXN2,MF,NPMAX,NPMIN,NFTMAX,NFT
      DOUBLE PRECISION X,XP,DUEPI,STEP,SUM,FTMAX,TUNEFOU,DELTAT,TUNE1,
     &TUNE
      COMPLEX*16 Z
      COMPLEX*16 ZSING(MAXITER)  ! Temp Z for CFFT, used to be SINGle precission
      DIMENSION X(MAXITER),XP(MAXITER)
      DIMENSION Z(MAXITER)

+ca crlibco
+if cr
+ca crcoall
+ei
C.............................................................
      IF (MAXN.GT.MAXITER) THEN
+if cr
        WRITE(lout,'(a)') '***ERROR(TUNENEWT): TOO MANY ITERATIONS'
+ei
+if .not.cr
        WRITE(*,   '(a)') '***ERROR(TUNENEWT): TOO MANY ITERATIONS'
+ei
        call prror(-1)
      ENDIF
C.............................................................
C    ESTIMATION OF TUNE WITH FFT 
C.............................................................
+if crlibm
      DUEPI=ATAN_RN(1D0)*8D0
      MFT=INT(LOG_RN(DBLE(MAXN))/LOG_RN(2D0)) 
+ei
+if .not.crlibm
      DUEPI=ATAN(1D0)*8D0
      MFT=INT(LOG(DBLE(MAXN))/LOG(2D0)) 
+ei
      NPOINT=2**MFT
      MAXN2=MAXN/2
      STEP=DUEPI/MAXN
C.............................................................
      SUM=0D0
      DO MF=1,MAXN
+if crlibm
        Z(MF)=DCMPLX(X(MF),XP(MF))*(1D0+COS_RN(STEP*(MF-MAXN2)))
+ei
+if .not.crlibm
        Z(MF)=DCMPLX(X(MF),XP(MF))*(1D0+COS(STEP*(MF-MAXN2)))
+ei
        ZSING(MF)=Z(MF)
        SUM=SUM+XP(MF)
      ENDDO 
      CALL CFFT(ZSING,-MFT)
C.........COMPUTATION OF THE MAXIMUM OF THE FOURIER TRANSFORM
      NPMIN=1
      IF (SUM.EQ.0D0) THEN
        NPMAX=NPOINT/2     !..REAL FFT ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
C............................................................
      FTMAX=0D0
      NFTMAX=0
C.............................REJECTS FREQUENCIES NEAR 0 OR 1
      DO NFT=NPMIN,NPMAX
        IF (ABS(ZSING(NFT)).GT.FTMAX) THEN
          FTMAX=ABS(ZSING(NFT))
          NFTMAX=NFT
        END IF
      ENDDO
      TUNEFOU=DBLE(NFTMAX-1)/DBLE(NPOINT)
      DELTAT=1D0/NPOINT
      TUNE1=TUNEFOU-DELTAT
      CALL ZFUN(TUNE,Z,MAXN,TUNE1,DELTAT)
      TUNENEWT=1D0-TUNE
C............................................................  
      RETURN 
C............................................................  
      END
C=============================================================
C AUXILIARY ROUTINE USED BY TUNENEWT.      
C
C AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
C

      SUBROUTINE ZFUN(TUNE,Z,MAXN,TUNEA1,DELTAT)
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER ND,MAXN,NUM,NTEST,NCONT,NFT,NC
      DOUBLE PRECISION DUEPI,ERR,DELTAT,TUNEA1,TUNEA2,DTUNEA1,DTUNEA2,
     &TUNE1,TUNE2,DTUNE1,DTUNE2,RATIO,TUNE3,DTUNE3,TUNETEST,TUNEVAL,
     &TUNEVMAX,TUNE
      COMPLEX*16 ZU,ZD,ZF,Z,ZTUNE1,ZTUNE2,ZTUNE3,ZFD
      DIMENSION Z(*),ZD(MAXITER),TUNETEST(10),TUNEVAL(10)

+ca crlibco

C............................................................  
+if crlibm
      DUEPI=ATAN_RN(1D0)*8D0
+ei
+if .not.crlibm
      DUEPI=ATAN(1D0)*8D0
+ei
      ERR=1D-10
      ZU=DCMPLX(0D0,1D0)
C............................................................
C.... we divide DELTAT in 5 parts
C............................................................
      DELTAT=DELTAT/5.D0
C............................................................  
      DO ND=1,MAXN
        ZD(ND)=(ZU*ND)*Z(ND)
      ENDDO
C............................................................  
+if crlibm
      ! EXP_RN expects a DOUBLE PRECISION, not COMPLEX -> rewrite expression for crlibm.
      ZTUNE1=cos_rn(DUEPI*TUNEA1) + ZU*sin_rn(DUEPI*TUNEA1)
+ei
+if .not.crlibm
      ZTUNE1=EXP((-ZU*DUEPI)*TUNEA1)
+ei
      CALL CALC(ZTUNE1,ZF,Z,MAXN)
      CALL CALC(ZTUNE1,ZFD,ZD,MAXN)
      DTUNEA1=DREAL(ZF)*DREAL(ZFD)+DIMAG(ZF)*DIMAG(ZFD)
      NUM=1
      DO NTEST=1, 10
        TUNEA2=TUNEA1+DELTAT
        ZTUNE2=EXP((-ZU*DUEPI)*TUNEA2)
        CALL CALC(ZTUNE2,ZF,Z,MAXN)
        CALL CALC(ZTUNE2,ZFD,ZD,MAXN)
        DTUNEA2=DREAL(ZF)*DREAL(ZFD)+DIMAG(ZF)*DIMAG(ZFD)
        IF ((DTUNEA1.LE.0D0).AND.(DTUNEA2.GE.0D0)) THEN
           TUNE1=TUNEA1
           TUNE2=TUNEA2
           DTUNE1=DTUNEA1
           DTUNE2=DTUNEA2
           DO NCONT=1,100
              RATIO=-DTUNE1/DTUNE2
              TUNE3=(TUNE1+RATIO*TUNE2)/(1.D0+RATIO)
              ZTUNE3=EXP((-ZU*DUEPI)*TUNE3)
              CALL CALC(ZTUNE3,ZF,Z,MAXN)
              CALL CALC(ZTUNE3,ZFD,ZD,MAXN)
              DTUNE3=DREAL(ZF)*DREAL(ZFD)+DIMAG(ZF)*DIMAG(ZFD)
              IF (DTUNE3.LE.0D0) THEN
                 IF(TUNE1.EQ.TUNE3) GOTO 100
                 TUNE1=TUNE3
                 DTUNE1=DTUNE3
              ELSE
                 IF(TUNE2.EQ.TUNE3) GOTO 100
                 TUNE2=TUNE3
                 DTUNE2=DTUNE3         
              ENDIF
              IF (ABS(TUNE2-TUNE1).LE.ERR) GOTO 100
           ENDDO
100        TUNETEST(NUM)=TUNE3
           TUNEVAL(NUM)=ABS(ZF)
           NUM=NUM+1
        ENDIF
        TUNEA1=TUNEA2
        DTUNEA1=DTUNEA2
      ENDDO
      TUNE=TUNETEST(1)
      TUNEVMAX=TUNEVAL(1)
      DO NC=2, NUM-1
         IF(TUNEVMAX.LE.TUNEVAL(NC)) THEN
            TUNEVMAX=TUNEVAL(NC)
            TUNE=TUNETEST(NC)
         ENDIF
      ENDDO
C............................................................  
      RETURN
C............................................................  
      END
C=============================================================
C AUXILIARY ROUTINE USED BY TUNENEWT.      
C
C AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
C

      SUBROUTINE CALC(ZV,ZPP,ZP,MAXD)
      IMPLICIT NONE
      INTEGER MAXD,NP
      COMPLEX*16 ZV,ZPP,ZP
      DIMENSION ZP(*)
      ZPP=ZP(MAXD)
C............................................................  
      DO NP=MAXD-1,1, -1
        ZPP=ZPP*ZV+ZP(NP)
      ENDDO
C............................................................  
      RETURN
C............................................................  
      END

CDECK  ID>, TUNEAPA.
C=============================================================
C COMPUTES THE TUNE AS THE AVERAGE PHASE ADVANCE ON A TWO
C DIMENSIONAL PLANE, GIVEN N ITERATES OF A MAP. TUNEAPA IS
C IN [0,1].
C X IS AN ARRAY CONTAINING THE FIRST COMPONENT OF THE ITERATES
C P IS AN  ARRAY CONTAINING THE SECOND COMPONENT OF THE
C   ITERATES
C N IS THE TOTAL NUMBER OF ITERATES USED IN COMPUTATION
C
C AUTHOR:    E. TODESCO - UNIVERSITY OF BOLOGNA
C
 
      DOUBLE PRECISION FUNCTION TUNEAPA(X,P,N)
C............................................................
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION X,P,PI,ADV,ADV1,ADV2,ADVS,S1,S2,S3,
     &ADVSIG,ADVMIN,ADVMAX
      DIMENSION X(*),P(*)

+ca crlibco
+if cr
+ca crcoall
+ei
C............................................................
      COMMON/TUNEPAR/ADVSIG,ADVMIN,ADVMAX
C............................................................
      IF (N.LE.0) THEN
+if cr
        WRITE(lout,'(a)')
+ei
+if .not.cr
        WRITE(*,   '(a)')
+ei
     &       '***ERROR(TUNEAPA): THIRD PARAMETER OUT OF BOUNDS'
        call prror(-1)
      ENDIF
C............................................................
+if crlibm
      PI=4*ATAN_RN(1.D0)
+ei
+if .not.crlibm
      PI=4*ATAN(1.D0)
+ei
      ADV=0.D0
      ADV1=10D0
      ADV2=-10D0
      ADVS=0D0
C.....................EVALUATION OF THE AVERAGE PHASE ADVANCE
      DO I=1,N-1
        S1=X(I+1)*X(I)+P(I+1)*P(I)
        S2=SQRT((X(I)*X(I)+P(I)*P(I))*(X(I+1)*X(I+1)+P(I+1)*P(I+1)))
+if crlibm
        S3=ACOS_RN(S1/S2)
+ei
+if .not.crlibm
        S3=ACOS(S1/S2)
+ei
        IF(-X(I+1)*P(I)+P(I+1)*X(I).LT.0) S3=-S3
        ADV=ADV+S3
C............................................................
        ADVS=ADVS+S3*S3
C............................................................
        ADV1=DMIN1(ADV1,S3)
        ADV2=DMAX1(ADV2,S3)
C............................................................
      ENDDO
C...............................................COMPUTES SIGMA
      ADVSIG=DSQRT((((ADVS/(N-1))/4)/PI)/PI-
     .            (((ADV/(N-1))/2)/PI)*(((ADV/(N-1))/2)/PI))
C......................................NORMALIZATION TO [0,1]
      ADV1=-ADV1
      IF(ADV1.LT.0) ADV1=1+ADV1
      ADV2=-ADV2
      IF(ADV2.LT.0) ADV2=1+ADV2
C...........................................FINDS MIN AND MAX
      ADVMIN=(MIN(ADV1,ADV2)/2)/PI
      ADVMAX=(MAX(ADV1,ADV2)/2)/PI
C......................................NORMALIZATION TO [0,1]
      ADV=-ADV
      IF(ADV.LT.0) THEN
        TUNEAPA=1+((ADV/(N-1))/2)/PI
      ELSE
        TUNEAPA=((ADV/(N-1))/2)/PI
      ENDIF
C.............................................................
      END
CDECK  ID>, TUNEFFT.
C=============================================================
C COMPUTES THE TUNE AS THE FFT ON A TWO DIMENSIONAL PLANE,
C GIVEN N ITERATES OF A MAP. THE FFT IS PERFORMED OVER THE
C MAXIMUM MFT WHICH SATIFIES 2**MFT.LE.N
C THE MAXIMUM NUMBER OF ITERATES IS FIXED IN A PARAMETER.
C TUNEFFT IS IN [0,1].
C X IS AN ARRAY CONTAINING THE FIRST COMPONENT OF THE ITERATES
C P IS AN ARRAY CONTAINING THE SECOND COMPONENT OF THE ITERATES
C N IS THE TOTAL NUMBER OF ITERATES USED IN COMPUTATION
C
C AUTHOR:    E. TODESCO - UNIVERSITY OF BOLOGNA
C
 
      DOUBLE PRECISION FUNCTION TUNEFFT(X,P,N)
C............................................................
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER N,M,NPOINT,I,NPMIN,NPMAX
      DOUBLE PRECISION SUM,AMAX
      DOUBLE PRECISION X(*),P(*)
      COMPLEX*16  Z(MAXITER) ! Temp Z for CFFT, used to be SINGle precission
+if cr
+ca crcoall
+ei
C..................................................CHECK OF N
      IF(N.GT.MAXITER) THEN
+if cr
        write(lout,'(a)')
+ei
+if .not.cr
        write(*,   '(a)')
+ei
     &       '***ERROR(TUNEFFT): TOO MANY ITERATES'
        call prror(-1)
      ENDIF
C............................................................
      IF (N.LE.0) THEN
+if cr
        write(lout,'(a)')
+ei
+if .not.cr
        write(*,   '(a)')
+ei
     &       '***ERROR(TUNEFFT): THIRD PARAMETER OUT OF BOUNDS'
        call prror(-1)
      ENDIF
C............................................COMPUTATION OF M
      DO M=1,50
        IF(2**M.GT.N) GOTO 10
      ENDDO
 10   M=M-1
C............................................................
      NPOINT=2**M
      SUM=0D0
C..................................SWITCH TO COMPLEX NOTATION
      DO I=1,NPOINT
        Z(I)=X(I)+(0.,1.)*P(I)
        SUM=SUM+P(I)
      ENDDO
C...........................................FOURIER TRANSFORM
      CALL CFFT(Z,-M)
C.........COMPUTATION OF THE MAXIMUM OF THE FOURIER TRANSFORM
      AMAX=0
      NPMIN=1
      IF (SUM.EQ.0D0) THEN
        NPMAX=NPOINT/2     !..REAL FFT ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
C.............................REJECTS FREQUENCIES NEAR 0 OR 1
      DO I=NPMIN,NPMAX
        IF(ABS(Z(I)).GT.AMAX) THEN
          TUNEFFT=FLOAT(I-1)/NPOINT
          AMAX=ABS(Z(I))
        ENDIF
      ENDDO
C......................................NORMALIZATION TO [0,1]
      TUNEFFT=-TUNEFFT
      IF(TUNEFFT.LE.0) THEN
        TUNEFFT=1+TUNEFFT
      ENDIF
C............................................................
      END
CDECK  ID>, TUNEFFTI.
C=============================================================
C COMPUTES THE TUNE AS THE FFT ON A TWO DIMENSIONAL PLANE,
C GIVEN N ITERATES OF A MAP. THE FFT IS PERFORMED OVER THE
C MAXIMUM MFT WHICH SATIFIES 2**MFT.LE.N. THEN, THE FFT IS
C INTERPOLATED FITTING THE THREE POINTS AROUND THE MAXIMUM
C USING A GAUSSIAN. THE TUNE IS COMPUTED AS THE MAXIMUM OF
C THE GAUSSIAN. THE MAXIMUM NUMBER OF ITERATES IS FIXED IN A
C PARAMETER. THE FUNCTION IS IN [0,1].
C X IS AN ARRAY CONTAINING THE FIRST COMPONENT OF THE ITERATES
C P IS AN ARRAY CONTAINING THE SECOND COMPONENT OF THE ITERATES
C N IS THE TOTAL NUMBER OF ITERATES USED IN COMPUTATION
C
C AUTHOR:     E. TODESCO - UNIVERSITY OF BOLOGNA
C
 
      DOUBLE PRECISION FUNCTION TUNEFFTI(X,P,N)
C............................................................
      IMPLICIT NONE
      INTEGER*4 MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER*4 N,M,NPOINT,NPMAX,NPMIN,ITUNE,I
      DOUBLE PRECISION SUM,AMAX,X1,X2,X3,Y1,Y2,Y3,X12,X13,Y12,Y13,X212,
     &X213,A,B
      DOUBLE PRECISION X(*),P(*)
      COMPLEX*16 Z(MAXITER)  ! Temp Z for CFFT, used to be SINGle precission

+ca crlibco
+if cr
+ca crcoall
+ei
C..................................................CHECK OF N
      IF(N.GT.MAXITER) THEN
+if cr
        write(lout,'(a)')
+ei
+if .not.cr
        write(*,   '(a)')
+ei
     &       '***ERROR(TUNEFFTI): TOO MANY ITERATES'
        call prror(-1)
      ENDIF
C............................................................
      IF (N.LE.0) THEN
+if cr
        write(lout,'(a)')
+ei
+if .not.cr
        write(*,   '(a)')
+ei
     &       '***ERROR(TUNEFFTI): THIRD PARAMETER OUT OF BOUNDS'
        call prror(-1)
      ENDIF
C............................................COMPUTATION OF M
      DO M=1,50
        IF(2**M.GT.N) GOTO 10
      ENDDO
 10   M=M-1
C............................................................
      NPOINT=2**M
      SUM=0D0
C..................................SWITCH TO COMPLEX NOTATION
      DO I=1,NPOINT
        Z(I)=X(I)+(0.,1.)*P(I)
        SUM=SUM+P(I)
      ENDDO
C...........................................FOURIER TRANSFORM
      CALL CFFT(Z,-M)
C.........COMPUTATION OF THE MAXIMUM OF THE FOURIER TRANSFORM
      AMAX=0
      NPMIN=1
      IF (SUM.EQ.0D0) THEN
        NPMAX=NPOINT/2     !..REAL FFT ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
C.............................REJECTS FREQUENCIES NEAR 0 OR 1
      DO I=NPMIN,NPMAX
        IF(ABS(Z(I)).GT.AMAX) THEN
          ITUNE=I
          AMAX=ABS(Z(I))
        ENDIF
      ENDDO
C..............................EVALUATION OF THE NEARBY PEAKS
      X1=ITUNE-1
      X2=ITUNE
      X3=ITUNE+1
      Y1=ABS(Z(ITUNE-1))
      Y2=ABS(Z(ITUNE))
      Y3=ABS(Z(ITUNE+1))
C...............................INTERPOLATION WITH A GAUSSIAN
      X12=X1-X2
      X13=X1-X3
+if crlibm
      Y12=LOG_RN(DBLE(Y1/Y2))
      Y13=LOG_RN(DBLE(Y1/Y3))
+ei
+if .not.crlibm
      Y12=LOG(Y1/Y2)
      Y13=LOG(Y1/Y3)
+ei
      X212=X1*X1-X2*X2
      X213=X1*X1-X3*X3
C........COMPUTATION OF THE POSITION OF THE INTERPOLATED PEAK
      A=X212*Y13-X213*Y12
      B=X12*Y13-X13*Y12
      TUNEFFTI=((A/2)/B-1)/NPOINT
C......................................NORMALIZATION TO [0,1]
      TUNEFFTI=-TUNEFFTI
      IF(TUNEFFTI.LE.0) THEN
        TUNEFFTI=1+TUNEFFTI
      ENDIF
C............................................................
      END

CDECK  ID>, TUNELASK.
C=============================================================
C COMPUTE THE TUNE OF A 2D MAP BY MEANS OF LASKAR METHOD.
C A FIRST INDICATION OF THE POSITION OF THE TUNE IS OBTAINED
C BY MEANS OF A FFT. REFINEMENT IS OBTAINED THROUGH A NEWTON
C PROCEDURE. THE MAXIMUM NUMBER OF ITERATIONS IS FIXED IN A
C PARAMETER. THE FUNCTION IS IN [0,1].
C X IS AN ARRAY CONTAINING THE ITERATION OF THE X COORDINATE OF
C   THE MAP
C PX IS AN ARRAY CONTAINING THE ITERATION OF THE PX COORDINATE
C   OF THE MAP
C MAX IS THE NUMBER OF ITERATIONS OF THE MAP
C
C AUTHOR:     R. BARTOLINI - BOLOGNA UNIVERSITY. 
C             M. GIOVANNOZZI - CERN HAS INTRODUCED SIGNIFICANT 
C             MODFICATIONS
C

      DOUBLE PRECISION FUNCTION TUNELASK(X,PX,MAX)
C............................................................
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER MAX,NPOINT,MF,NFTMAX,NPMIN,NPMAX,NFT,JITER,JMAX,MAX1,JIT,
     &I,K,JOM,J,N,MBODE,MFT
      DOUBLE PRECISION DUEPI,SUM,FTMAX,TUNEFOU,FSIGMA,OMEMIN,STEP,
     &OMEMAX,FOMEGA,OME,ABSFOM,TMPR,TMPI
      DOUBLE PRECISION X(*),PX(*)
      COMPLEX*16 ZSING(MAXITER) ! Temp Z for CFFT, used to be SINGle precission
      COMPLEX*16 Z(MAXITER),FOME,ZC,SD,SP

+ca crlibco
+if cr
+ca crcoall
+ei
+if crlibm
      DUEPI=8*ATAN_RN(1D+0)
+ei
+if .not.crlibm
      DUEPI=8*ATAN(1D+0)
+ei
C...............................CHECK OF THE ITERATION NUMBER
      IF(MAX.GT.MAXITER) THEN
+if cr
        write(lout,'(a)')
+ei
+if .not.cr
        write(*,   '(a)')
+ei
     &       '***ERROR(TUNELASK): TOO MANY ITERATIONS'
        call prror(-1)
      ENDIF
C............................................................
      IF (MAX.LE.0) THEN
+if cr
        write(lout,'(a)')
+ei
+if .not.cr
        write(*,   '(a)')
+ei
     &       '***ERROR(TUNELASK): THIRD PARAMETER OUT OF BOUNDS'
        call prror(-1)
      ENDIF
C.................................ESTIMATION OF TUNE WITH FFT
      SUM=0D0
+if crlibm
      MFT=INT(LOG_RN(DBLE(MAX))/LOG_RN(2D+0))
+ei
+if .not.crlibm
      MFT=INT(LOG(DBLE(MAX))/LOG(2D+0))
+ei
      NPOINT=2**MFT
      DO MF=1,NPOINT
        ZSING(MF)=X(MF)-(0.,1.)*PX(MF)
        SUM=SUM+PX(MF)
      ENDDO
      CALL CFFT(ZSING,-MFT)
C......................SEARCH FOR MAXIMUM OF FOURIER SPECTRUM
      FTMAX=0D+0
      NFTMAX=0D+0
C............................................................
      NPMIN=1
      IF (SUM.EQ.0D0) THEN
        NPMAX=NPOINT/2 !..REAL FFT ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
C.............................REJECTS FREQUENCIES NEAR 0 OR 1
      DO NFT=NPMIN,NPMAX
        IF (ABS(ZSING(NFT)).GT.FTMAX) THEN
          FTMAX=ABS(ZSING(NFT))
          NFTMAX=NFT
        ENDIF
      ENDDO
      TUNEFOU=FLOAT(NFTMAX-1)/NPOINT
C..............DEFINITION OF THE FIRST INTERVAL FOR BISECTION
C..............PROCEDURE AND INIZIALIZATION OF BISECTION
C..............PROCEDURE: THE INTERVAL (OMEMIN,OMEMAX) IS
C..............DIVIDED IN JMAX-1 SUBINTERVAL AND THIS
C..............SUBDIVISION IS ITERATED JITER TIMES
      FSIGMA=1D+0
      OMEMIN=TUNEFOU*DUEPI-(FSIGMA*DUEPI)/NPOINT
      OMEMAX=TUNEFOU*DUEPI+(FSIGMA*DUEPI)/NPOINT
C.....JITER=8 PROVIDES A PRECISION WHICH IS 1.5E-4 * 1/2**MFT
      JITER=8
C................JMAX=7 IS THE VALUE WHICH MINIMIZES CPU TIME
      JMAX=7
      MAX1=MAX-1
      STEP=(DUEPI*.5D0)/MAX1
C.........................................BISECTION PROCEDURE
      DO JIT=1,JITER
        FOMEGA=0D+0
        DO J=1,JMAX
          OME=OMEMIN+((OMEMAX-OMEMIN)/(JMAX-1D+0))*(J-1D+0)
          DO N=1,MAX
+if crlibm
            ZC=(X(N)-(0D+0,1D+0)*PX(N))
     .        *(1D0+COS_RN(STEP*(2*N-MAX1)))
            TMPR=DREAL((-(0D+0,1D+0)*OME)*N)
            TMPI=DIMAG((-(0D+0,1D+0)*OME)*N)
            Z(N)=ZC*(EXP_RN(TMPR)*DCMPLX(COS_RN(TMPI),SIN_RN(TMPI))) !exp_rn is only defined for real numbers -> decompose in real and imaginary part
+ei
+if .not.crlibm
            ZC=(X(N)-(0D+0,1D+0)*PX(N))
     .        *(1D0+COS(STEP*(2*N-MAX1)))
            Z(N)=ZC*EXP(-(0D+0,1D+0)*OME*N)
+ei
          ENDDO
C..COMPUTATION OF SCALAR PRODUCT WITH ITERATED BODE ALGORITHM
          FOME=(0D0,0D0)
          MBODE=(MAX-5)/4
          DO I=0,MBODE
            K=4*I
            FOME=FOME+((7D0*(Z(1+K)+Z(5+K))+32D0*(Z(2+K)+Z(4+K)))
     .               +12D0*Z(3+K))
          ENDDO
          FOME=((.5D0*FOME)/45D0)/DBLE(MBODE+1)
C..........SEARCH FOR MAXIMUM OF SCALAR PRODUCT AND DEFINITION
C..........OF THE NEW INTERVAL (OMEMIN,OMEMAX) WHERE TO
C........................................RESTART THE PROCEDURE
          ABSFOM=ABS(FOME)
          IF (ABSFOM.GT.FOMEGA) THEN
            FOMEGA=ABSFOM
            TUNELASK=-OME/DUEPI
            JOM=J
          ENDIF
        ENDDO
        OMEMIN=OMEMIN+((OMEMAX-OMEMIN)/(JMAX-1D+0))*(JOM-2D+0)
        OMEMAX=OMEMIN+((OMEMAX-OMEMIN)/(JMAX-1D+0))*JOM
      ENDDO
C......................................NORMALIZATION TO [0,1]
      TUNELASK=-TUNELASK
      IF(TUNELASK.LE.0) THEN
        TUNELASK=1+TUNELASK
      ENDIF
C............................................................
      END
