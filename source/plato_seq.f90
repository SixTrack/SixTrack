module platoFMA
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants, only : zero, half, one, two, four, five, seven, eight, c1e1, c1m1, c1m10

      implicit none

      contains
!*
!* $Id: cfft.F,v 1.1.1.1 1996/02/15 17:48:48 mclareni Exp $
!*
!* $Log: cfft.F,v $
!* Revision 1.1.1.1  1996/02/15 17:48:48  mclareni
!* Kernlib
!* Conversion single->double precission by K.Sjobak, Dec. 2016
!*
      SUBROUTINE CFFT(A,MSIGN)
      IMPLICIT NONE
      COMPLEX(kind=fPrec) A(*),U,W,T
      INTEGER MSIGN, M, N, NV2, NM1, J, I, K, L, LE, LE1, IP
      real(kind=fPrec) C,S
      IF(MSIGN.EQ.0) RETURN
      M=IABS(MSIGN)
      N=2**M
      NV2=N/2
      NM1=N-1
      J=1
      DO I=1,NM1
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
      END DO

      DO I=1,N,2
        T=A(I+1)
        A(I+1)=A(I)-T
        A(I )=A(I)+T
      END DO

      IF(M.EQ.1) RETURN
      C=zero
      S=ISIGN(1,MSIGN)
      LE=2

      DO L=2,M
        W=CMPLX(C,S,fPrec)
        U=W
        C=SQRT(C*half+half)
        S=AIMAG(W)/(C+C)
        LE1=LE
        LE=LE1+LE1

        DO I=1,N,LE
          IP=I+LE1
          T=A(IP)
          A(IP)=A(I)-T
          A(I) =A(I)+T
        END DO

        DO J=2,LE1
          DO I=J,N,LE
            IP=I+LE1
            T=A(IP)*U
            A(IP)=A(I)-T
            A(I) =A(I)+T
          END DO
          U=U*W
        END DO
      END DO

      RETURN
      END SUBROUTINE

!CDECK  ID>, TUNENEWT1.
!C=============================================================
!C COMPUTES THE TUNE USING A DISCRETE VERSION OF LASKAR METHOD.
!C IT INCLUDES A NEWTON METHOD FOR THE SEARCH OF THE FREQUENCY.
!C X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE
!C LENGTH IF THE ORBIT.
!C
!C AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
!C             M. GIOVANNOZZI - CERN HAS INTRODUCED SOME
!C                                   MODIFICATIONS
!C

      REAL(KIND=fPrec) FUNCTION TUNENEWT1(X,XP,MAXN)
      use crcoall
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      REAL(KIND=fPrec) X,XP, DUEPI, STEP, SUM, FTMAX, TUNEFOU, DELTAT, TUNE1, TUNE
      INTEGER MAXN, MFT, NPOINT, MAXN2, MF, NPMIN, NPMAX, NFTMAX, NFT
      COMPLEX(kind=fPrec) Z
      COMPLEX(kind=fPrec), DIMENSION(MAXITER) :: ZSING ! Temp Z for CFFT, used to be SINGLE precision

      DIMENSION X(MAXITER),XP(MAXITER)
      DIMENSION Z(MAXITER)

!C.............................................................
      IF (MAXN.GT.MAXITER) THEN
         WRITE(lout,'(a)') '***ERROR(TUNENEWT1): TOO MANY ITERATIONS'
         call prror(-1)
      ENDIF
!C.............................................................
!C    ESTIMATION OF TUNE WITH FFT
!C.............................................................
      DUEPI=ATAN_MB(one)*eight
      MFT=INT(LOG_MB(REAL(MAXN,fPrec))/LOG_MB(two))
      NPOINT=2**MFT
      MAXN2=MAXN/2
      STEP=DUEPI/MAXN
!.............................................................
      SUM=zero
      DO MF=1,MAXN
        Z(MF)=CMPLX(X(MF),XP(MF),fPrec) !Returns COMPLEX*16 / COMPLEX(8)
        ZSING(MF)=Z(MF)
        SUM=SUM+XP(MF)
      ENDDO
      CALL CFFT(ZSING,-MFT)
!.........COMPUTATION OF THE MAXIMUM OF THE FOURIER TRANSFORM
      NPMIN=1
      IF (SUM.EQ.zero) THEN
        NPMAX=NPOINT/2  !..REAL FFT ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
!............................................................
      FTMAX=zero
      NFTMAX=0
!C.............................REJECTS FREQUENCIES NEAR 0 OR 1
      DO NFT=NPMIN,NPMAX
        IF (ABS(ZSING(NFT)).GT.FTMAX) THEN
          FTMAX=ABS(ZSING(NFT))
          NFTMAX=NFT
        END IF
      ENDDO
      TUNEFOU=REAL(NFTMAX-1,fPrec)/REAL(NPOINT,fPrec)
      DELTAT=one/NPOINT
      TUNE1=TUNEFOU-DELTAT
      CALL ZFUN(TUNE,Z,MAXN,TUNE1,DELTAT)
      TUNENEWT1=one-TUNE
!C............................................................
      RETURN
!C............................................................
      END FUNCTION

!CDECK  ID>, TUNEFIT.
!C============================================================
!C COMPUTES THE TUNE USING A MODIFIED APA ALGORITHM. THE FIRST
!C STEP CONSISTS OF TAKING THE AVERAGE OF TUNE COMPUTED WITH THE
!C APA METHOD, THEN A BEST FIT IS PERFORMED.
!C
!C     AUTHORS: R. BARTOLINI A. BAZZANI - BOLOGNA UNIVERSITY
!C

      REAL(KIND=fPrec) FUNCTION TUNEFIT(X,XP,MAX)
      use crcoall
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER MAX,MAX1,MAX2,N
      REAL(KIND=fPrec) X,XP,DUEPI,C,SUMAPM,PA,R1,R2,CTHETA,STHETA,THETA,TUNEAPA,U,ATUNE,TUNEF,TUNE
      DIMENSION X(*),XP(*),TUNE(MAXITER),U(MAXITER)

!C............................................................
      IF (MAX.LE.0) THEN
         WRITE(lout,'(a)') '***ERROR(TUNEFIT): THIRD PARAMETER OUT OF BOUNDS'
        call prror(-1)
      ENDIF
!C............................................................
      DUEPI=eight*ATAN_MB(one)
!C.............................................................
      C=REAL(MAX,fPrec)/TWO-REAL(INT(MAX/TWO),fPrec)
      MAX1=MAX
      IF(C.GT.C1M1) MAX1=MAX-1
      MAX2=MAX1/2
!C.............................................................
      SUMAPM=zero
      PA=zero
!C.............................................................
      DO N=1,MAX-1
!C.............................................................
        R1=SQRT(X(N)**2+XP(N)**2)
        R2=SQRT(X(N+1)**2+XP(N+1)**2)
        CTHETA=((X(N)*X(N+1)+XP(N)*XP(N+1))/R1)/R2
        STHETA=((-X(N)*XP(N+1)+X(N+1)*XP(N))/R1)/R2
        IF (STHETA.GE.0) THEN
          THETA=ACOS_MB(CTHETA)
          ELSE
          THETA=DUEPI-ACOS_MB(CTHETA)
        END IF
        PA=PA+THETA
        TUNEAPA=PA/DUEPI
        SUMAPM=SUMAPM+TUNEAPA
        IF (N.GE.MAX2) THEN
          TUNE(N-MAX2+1)=((SUMAPM/REAL(N,fPrec))/REAL(N+1,fPrec))*TWO
          U(N-MAX2+1)=ONE/N
        ENDIF
!C.............................................................
      ENDDO
      CALL FIT(U,TUNE,MAX1-MAX2,ATUNE,TUNEF)
      TUNEFIT=TUNEF
!C.............................................................
      RETURN
!C.............................................................
      END FUNCTION
!CDECK  ID>, FIT.
!C============================================================
!C COMPUTES THE STRAIGHT LINE WHICH FITS A SET OF N 2D
!C OBSERVATIONS.
!C THE LINE IS REPRESENTED BY
!C
!C Y=AX+B
!C
!C     AUTHORS: R. BARTOLINI A. BAZZANI - BOLOGNA UNIVERSITY
!C

      SUBROUTINE FIT(X,Y,N,A,B)
      use crcoall
      IMPLICIT NONE
      INTEGER N,I
      REAL(KIND=fPrec) X2M,XM,YM,XYM,ERR,DELTA,A,B,X,Y
      DIMENSION X(*),Y(*)
!C............................................................
      IF (N.LE.0) THEN
        WRITE(lout,'(a)') '***ERROR(TUNEAPA): THIRD PARAMETER OUT OF BOUNDS'
        call prror(-1)
      ENDIF
!C.............................................................
      X2M=zero
      XM=zero
      YM=zero
      XYM=zero
      ERR=zero
!C.............................................................
      DO I=1,N
        X2M=X2M+X(I)*X(I)
        XYM=XYM+X(I)*Y(I)
        XM=XM+X(I)
        YM=YM+Y(I)
      END DO
!C.............................................................
      X2M=X2M/N
      XYM=XYM/N
      XM=XM/N
      YM=YM/N
      DELTA=X2M-XM*XM
      A=(XYM-XM*YM)/DELTA
      B=(-XM*XYM+X2M*YM)/DELTA
!C.............................................................
      RETURN
!C.............................................................
      END SUBROUTINE

!CDECK  ID>, TUNEABT.
!C=============================================================
!C COMPUTES THE TUNE USING FFT INTERPOLATED METHOD.
!C X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE
!C LENGTH OF THE ORBIT.
!C
!C AUTHOR:     E. TODESCO - INFN AND CERN
!C

      REAL(KIND=fPrec) FUNCTION TUNEABT(X,XP,MAXN)
      use numerical_constants, only : pi
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER MAXN,NPOINT,MF,NPMIN,NPMAX,NFTMAX,NFT
      REAL(KIND=fPrec) X,XP,DUEPI,MFT,STEP,SUM,FTMAX,CF1,CF2,CF3,ASSK
      COMPLEX(kind=fPrec) Z,ZSING
      DIMENSION X(MAXITER),XP(MAXITER)
      DIMENSION Z(MAXITER),ZSING(MAXITER)

!C..................................ESTIMATION OF TUNE WITH FFT
      MFT=INT(LOG_MB(REAL(MAXN,fPrec))/LOG_MB(two))
      DUEPI=two*PI
      NPOINT=2**MFT
      STEP=DUEPI/NPOINT/two
!C.............................................................
      SUM=zero            !..CHECKS FOR COMPLEX OR REAL DATA
      DO MF=1,NPOINT
        Z(MF)=CMPLX(X(MF),XP(MF),fPrec)
        ZSING(MF)=Z(MF)
        SUM=SUM+XP(MF)
      ENDDO
      CALL FFT_PLATO(ZSING,NPOINT,-1)
!C.......................SEARCH FOR MAXIMUM OF FOURIER SPECTRUM
      NPMIN=1
      IF (SUM.EQ.zero) THEN
        NPMAX=NPOINT/2         !..REAL FFT_PLATO ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
!C.............................................................
      FTMAX=zero
      NFTMAX=0
      DO NFT=NPMIN,NPMAX
        IF (ABS(ZSING(NFT)).GT.FTMAX) THEN
          FTMAX=ABS(ZSING(NFT))
          NFTMAX=NFT
        END IF
      ENDDO
!C................................................INTERPOLATION
      CF1=ABS(ZSING(NFTMAX-1))
      CF2=ABS(ZSING(NFTMAX))
      CF3=ABS(ZSING(NFTMAX+1))
      IF (CF3.GT.CF1) THEN
        ASSK=REAL(NFTMAX,fPrec)+(NPOINT/PI)*ATAN2_MB(CF3*SIN_MB(PI/NPOINT),CF2+CF3*COS_MB(PI/NPOINT))
      ELSE
        ASSK=REAL(NFTMAX-1,fPrec)+(NPOINT/PI)*ATAN2_MB(CF2*SIN_MB(PI/NPOINT),CF1+CF2*COS_MB(PI/NPOINT))
      ENDIF
      TUNEABT=one-(ASSK-one)/REAL(NPOINT,fPrec) !1D+0 = 1D0, i.e. real(kind=fPrec) 1.0?
!C............................................................
      RETURN
!C............................................................
      END FUNCTION
!CDECK  ID>, TUNEABT2.
!C=============================================================
!C COMPUTES THE TUNE USING THE INTERPOLATED FFT METHOD
!C WITH HANNING FILTER.
!C X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE
!C LENGTH OF THE ORBIT.
!C
!C AUTHOR:     E. TODESCO - INFN AND CERN
!C

      REAL(KIND=fPrec) FUNCTION TUNEABT2(X,XP,MAXN)
      use numerical_constants, only : pi
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER MAXN,NPOINT,MFT,MF,NPMIN,NPMAX,NFTMAX,NN,NFT
      REAL(KIND=fPrec) X,XP,DUEPI,STEP,SUM,FTMAX,CF1,CF2,CF3,P1,P2,CO,SI,SCRA1,SCRA2,SCRA3,SCRA4,ASSK
      COMPLEX(kind=fPrec) Z
      COMPLEX(kind=fPrec) ZSING(MAXITER)
      DIMENSION X(MAXITER),XP(MAXITER)
      DIMENSION Z(MAXITER)

!C..................................ESTIMATION OF TUNE WITH FFT
      MFT=INT(LOG_MB(REAL(MAXN,fPrec))/LOG_MB(two))
      DUEPI=2*PI
      NPOINT=2**MFT
      STEP=(DUEPI/NPOINT)/two
!C.............................................................
      SUM=zero            !..CHECKS FOR COMPLEX OR REAL DATA
      DO MF=1,NPOINT
        Z(MF)=CMPLX(X(MF),XP(MF),fPrec)*SIN_MB(STEP*MF)**2
        ZSING(MF)=Z(MF)
        SUM=SUM+XP(MF)
      ENDDO
      CALL FFT_PLATO(ZSING,NPOINT,-1)
!C.......................SEARCH FOR MAXIMUM OF FOURIER SPECTRUM
      NPMIN=1
      IF (SUM.EQ.zero) THEN
        NPMAX=NPOINT/2         !..REAL FFT_PLATO ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
!C.............................................................
      FTMAX=zero
      NFTMAX=0
      DO NFT=NPMIN,NPMAX
        IF (ABS(ZSING(NFT)).GT.FTMAX) THEN
          FTMAX=ABS(ZSING(NFT))
          NFTMAX=NFT
        END IF
      ENDDO

      ! VKBO Bugfix for Debug build type
      ! Make sure nftmax is greater than 1, or we may get a segfault in the next lines
      if(nftmax <= 1) then
        tuneabt2 = one
        return
      end if

      CF1=ABS(ZSING(NFTMAX-1))
      CF2=ABS(ZSING(NFTMAX))
      CF3=ABS(ZSING(NFTMAX+1))
      IF (CF3.GT.CF1) THEN
        P1=CF2
        P2=CF3
        NN=NFTMAX
      ELSE
        P1=CF1
        P2=CF2
        NN=NFTMAX-1
      ENDIF
!C..........................................INTERPOLATION
      CO=COS_MB((2*PI)/REAL(NPOINT,fPrec))
      SI=SIN_MB((2*PI)/REAL(NPOINT,fPrec))
      SCRA1=CO**2*(P1+P2)**2-((2*P1)*P2)*((2*CO**2-CO)-1)
      SCRA2=(P1+P2*CO)*(P1-P2)
      SCRA3=(P1**2+P2**2)+((2*P1)*P2)*CO
      SCRA4=(-SCRA2+P2*SQRT(SCRA1))/SCRA3
      ASSK=REAL(NN,fPrec)+((NPOINT/2)/PI)*ASIN_MB(SI*SCRA4)
      TUNEABT2=one-(ASSK-one)/REAL(NPOINT,fPrec)
!C............................................................
      RETURN
!C............................................................
      END FUNCTION

!CDECK  ID>, FFT_PLATO.
!C============================================================
!C           COMPUTES THE FFT_PLATO   (REAL(KIND=fPrec))
!C           AUTHOR: NUMERICAL RECEIPES, PG. 395
!C           NN IS THE NUMBER OF DATA: MUST BE A POWER OF 2
!C           ISIGN=1: DIRECT FT
!C           ISIGN=-1: INVERSE FT
!C           DATA IS A COMPLEX ARRAY WITH THE SIGNAL IN INPUT
!C                                   WITH THE FT IN OUTPUT
!C           NOTE THAT REAL*8 DATA(*) THAT CONTAINS THE COMPLEX
!C             NUMBER IN THE FOLLOWING ORDER:
!C             REAL(DATA(1)),IMAG(DATA(1)),REAL(DATA(2)),IMAG(DATA(2)),...
!C           GOOD LUCK, BABY
!C

      SUBROUTINE FFT_PLATO(CDATA,NN,ISIGN)
      IMPLICIT NONE
      INTEGER I,J,N,NN,M,MMAX,ISTEP,ISIGN
      REAL(kind=fPrec) WR,WI,WPR,WPI,WTEMP,THETA,TEMPR,TEMPI
      COMPLEX(kind=fPrec) CDATA(NN)
      REAL(kind=fPrec) DATA(2*NN)

      N=2*NN
!C create real array DATA out of complex array CDATA
      DO I=1,N,2
        DATA(I)=REAL(CDATA(I/2+1))
        DATA(I+1)=AIMAG(CDATA(I/2+1))
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
        THETA=(EIGHT*ATAN_MB(ONE))/(ISIGN*MMAX)
        WPR=-TWO*SIN_MB(HALF*THETA)**2
        WPI=SIN_MB(THETA)
        WR=ONE
        WI=ZERO
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
        CDATA(I/2+1)=CMPLX(DATA(I),DATA(I+1),fPrec)
      ENDDO
!C.............................................................
      END SUBROUTINE

!CDECK  ID>, FFT_PLATO_REAL.
!C============================================================
!C           COMPUTES THE FFT_PLATO_REAL   (REAL(KIND=fPrec))
!C           AUTHOR: NUMERICAL RECEIPES, PG. 395
!C           NN IS THE NUMBER OF DATA: MUST BE A POWER OF 2
!C           ISIGN=1: DIRECT FT
!C           ISIGN=-1: INVERSE FT
!C           DATA IS A COMPLEX ARRAY WITH THE SIGNAL IN INPUT
!C                                   WITH THE FT IN OUTPUT
!C           NOTE THAT REAL*8 DATA(*) THAT CONTAINS THE COMPLEX
!C             NUMBER IN THE FOLLOWING ORDER:
!C             REAL(DATA(1)),IMAG(DATA(1)),REAL(DATA(2)),IMAG(DATA(2)),...
!C           GOOD LUCK, BABY
!C

      SUBROUTINE FFT_PLATO_REAL(DATA,NN,ISIGN)
      IMPLICIT NONE
      INTEGER I,J,N,NN,M,MMAX,ISTEP,ISIGN
      REAL(kind=fPrec) WR,WI,WPR,WPI,WTEMP,THETA,TEMPR,TEMPI
      REAL(kind=fPrec) DATA(*)

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
        THETA=(EIGHT*ATAN_MB(ONE))/(ISIGN*MMAX)
        WPR=-TWO*SIN_MB(HALF*THETA)**2
        WPI=SIN_MB(THETA)
        WR=ONE
        WI=ZERO
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
!C.............................................................
      END SUBROUTINE

!CDECK  ID>, TUNENEWT.
!C=============================================================
!C COMPUTES THE TUNE USING A DISCRETE VERSION OF LASKAR METHOD.
!C IT INCLUDES A NEWTON METHOD FOR THE SEARCH OF THE FREQUENCY.
!C X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE
!C LENGTH IF THE ORBIT.
!C
!C AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
!C             M. GIOVANNOZZI - CERN HAS INTRODUCED SOME
!C                                   MODIFICATIONS
!C

      REAL(KIND=fPrec) FUNCTION TUNENEWT(X,XP,MAXN)
      use crcoall
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER MAXN,MFT,NPOINT,MAXN2,MF,NPMAX,NPMIN,NFTMAX,NFT
      REAL(KIND=fPrec) X,XP,DUEPI,STEP,SUM,FTMAX,TUNEFOU,DELTAT,TUNE1,TUNE
      COMPLEX(kind=fPrec) Z
      COMPLEX(kind=fPrec), DIMENSION(MAXITER) :: ZSING  ! Temp Z for CFFT, used to be SINGle precission
      DIMENSION X(MAXITER),XP(MAXITER)
      DIMENSION Z(MAXITER)

!C.............................................................
      IF (MAXN.GT.MAXITER) THEN
        WRITE(lout,'(a)') '***ERROR(TUNENEWT): TOO MANY ITERATIONS'
        call prror(-1)
      ENDIF
!C.............................................................
!C    ESTIMATION OF TUNE WITH FFT
!C.............................................................
      DUEPI=ATAN_MB(ONE)*EIGHT
      MFT=INT(LOG_MB(REAL(MAXN,fPrec))/LOG_MB(TWO))
      NPOINT=2**MFT
      MAXN2=MAXN/2
      STEP=DUEPI/MAXN
!C.............................................................
      SUM=zero
      DO MF=1,MAXN
        Z(MF)=CMPLX(X(MF),XP(MF),fPrec)*(one+COS_MB(STEP*(MF-MAXN2)))
        ZSING(MF)=Z(MF)
        SUM=SUM+XP(MF)
      ENDDO
      CALL CFFT(ZSING,-MFT)
!C.........COMPUTATION OF THE MAXIMUM OF THE FOURIER TRANSFORM
      NPMIN=1
      IF (SUM.EQ.ZERO) THEN
        NPMAX=NPOINT/2     !..REAL FFT ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
!C............................................................
      FTMAX=ZERO
      NFTMAX=0
!C.............................REJECTS FREQUENCIES NEAR 0 OR 1
      DO NFT=NPMIN,NPMAX
        IF (ABS(ZSING(NFT)).GT.FTMAX) THEN
          FTMAX=ABS(ZSING(NFT))
          NFTMAX=NFT
        END IF
      ENDDO
      TUNEFOU=REAL(NFTMAX-1,fPrec)/REAL(NPOINT)
      DELTAT=ONE/NPOINT
      TUNE1=TUNEFOU-DELTAT
      CALL ZFUN(TUNE,Z,MAXN,TUNE1,DELTAT)
      TUNENEWT=ONE-TUNE
!C............................................................
      RETURN
!C............................................................
      END FUNCTION
!C=============================================================
!C AUXILIARY ROUTINE USED BY TUNENEWT.
!C
!C AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
!C

      SUBROUTINE ZFUN(TUNE,Z,MAXN,TUNEA1,DELTAT)
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER ND,MAXN,NUM,NTEST,NCONT,NFT,NC
      REAL(KIND=fPrec) DUEPI,ERR,DELTAT,TUNEA1,TUNEA2,DTUNEA1,DTUNEA2, &
     &TUNE1,TUNE2,DTUNE1,DTUNE2,RATIO,TUNE3,DTUNE3,TUNETEST,TUNEVAL,TUNEVMAX,TUNE
      COMPLEX(kind=fPrec) ZU,ZD,ZF,Z,ZTUNE1,ZTUNE2,ZTUNE3,ZFD
      DIMENSION Z(*),ZD(MAXITER),TUNETEST(10),TUNEVAL(10)

!C............................................................
      DUEPI=ATAN_MB(one)*eight
      ERR=c1m10
      ZU=CMPLX(zero,one,fPrec)
!C............................................................
!C.... we divide DELTAT in 5 parts
!C............................................................
      DELTAT=DELTAT/FIVE
!C............................................................
      DO ND=1,MAXN
        ZD(ND)=(ZU*ND)*Z(ND)
      ENDDO
!C............................................................
#ifdef CRLIBM
      ! EXP_MB expects a REAL(KIND=fPrec), not COMPLEX -> rewrite expression for crlibm.
      ZTUNE1=cos_mb(DUEPI*TUNEA1) + ZU*sin_mb(DUEPI*TUNEA1)
#else
      ZTUNE1=EXP((-ZU*DUEPI)*TUNEA1)
#endif
      CALL CALC(ZTUNE1,ZF,Z,MAXN)
      CALL CALC(ZTUNE1,ZFD,ZD,MAXN)
      DTUNEA1=REAL(ZF)*REAL(ZFD)+AIMAG(ZF)*AIMAG(ZFD)
      NUM=1
      DO NTEST=1, 10
        TUNEA2=TUNEA1+DELTAT
        ZTUNE2=EXP((-ZU*DUEPI)*TUNEA2)
        CALL CALC(ZTUNE2,ZF,Z,MAXN)
        CALL CALC(ZTUNE2,ZFD,ZD,MAXN)
        DTUNEA2=REAL(ZF)*REAL(ZFD)+AIMAG(ZF)*AIMAG(ZFD)
        IF ((DTUNEA1.LE.ZERO).AND.(DTUNEA2.GE.ZERO)) THEN
           TUNE1=TUNEA1
           TUNE2=TUNEA2
           DTUNE1=DTUNEA1
           DTUNE2=DTUNEA2
           DO NCONT=1,100
              RATIO=-DTUNE1/DTUNE2
              TUNE3=(TUNE1+RATIO*TUNE2)/(ONE+RATIO)
              ZTUNE3=EXP((-ZU*DUEPI)*TUNE3)
              CALL CALC(ZTUNE3,ZF,Z,MAXN)
              CALL CALC(ZTUNE3,ZFD,ZD,MAXN)
              DTUNE3=REAL(ZF)*REAL(ZFD)+AIMAG(ZF)*AIMAG(ZFD)
              IF (DTUNE3.LE.ZERO) THEN
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
!C............................................................
      RETURN
!C............................................................
      END SUBROUTINE
!C=============================================================
!C AUXILIARY ROUTINE USED BY TUNENEWT.
!C
!C AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
!C

      SUBROUTINE CALC(ZV,ZPP,ZP,MAXD)
      IMPLICIT NONE
      INTEGER MAXD,NP
      COMPLEX(kind=fPrec) ZV,ZPP,ZP
      DIMENSION ZP(*)
      ZPP=ZP(MAXD)
!C............................................................
      DO NP=MAXD-1,1, -1
        ZPP=ZPP*ZV+ZP(NP)
      ENDDO
!C............................................................
      RETURN
!C............................................................
      END SUBROUTINE

!CDECK  ID>, TUNEAPA.
!C=============================================================
!C COMPUTES THE TUNE AS THE AVERAGE PHASE ADVANCE ON A TWO
!C DIMENSIONAL PLANE, GIVEN N ITERATES OF A MAP. TUNEAPA IS
!C IN [0,1].
!C X IS AN ARRAY CONTAINING THE FIRST COMPONENT OF THE ITERATES
!C P IS AN  ARRAY CONTAINING THE SECOND COMPONENT OF THE
!C   ITERATES
!C N IS THE TOTAL NUMBER OF ITERATES USED IN COMPUTATION
!C
!C AUTHOR:    E. TODESCO - UNIVERSITY OF BOLOGNA
!C

      REAL(KIND=fPrec) FUNCTION TUNEAPA(X,P,N)
!C............................................................
      use crcoall
      use numerical_constants, only : pi
      IMPLICIT NONE
      INTEGER N,I
      REAL(KIND=fPrec) X,P,ADV,ADV1,ADV2,ADVS,S1,S2,S3,ADVSIG,ADVMIN,ADVMAX
      DIMENSION X(*),P(*)
!C............................................................
      COMMON/TUNEPAR/ADVSIG,ADVMIN,ADVMAX
!C............................................................
      IF (N.LE.0) THEN
        WRITE(lout,'(a)') '***ERROR(TUNEAPA): THIRD PARAMETER OUT OF BOUNDS'
        call prror(-1)
      ENDIF
!C............................................................
      ADV=ZERO
      ADV1=C1E1
      ADV2=-C1E1
      ADVS=ZERO
!C.....................EVALUATION OF THE AVERAGE PHASE ADVANCE
      DO I=1,N-1
        S1=X(I+1)*X(I)+P(I+1)*P(I)
        S2=SQRT((X(I)*X(I)+P(I)*P(I))*(X(I+1)*X(I+1)+P(I+1)*P(I+1)))
        S3=ACOS_MB(S1/S2)
        IF(-X(I+1)*P(I)+P(I+1)*X(I).LT.0) S3=-S3
        ADV=ADV+S3
!C............................................................
        ADVS=ADVS+S3*S3
!C............................................................
        ADV1=MIN(ADV1,S3)
        ADV2=MAX(ADV2,S3)
!C............................................................
      ENDDO
!C...............................................COMPUTES SIGMA
      ADVSIG=SQRT((((ADVS/(N-1))/4)/PI)/PI-(((ADV/(N-1))/2)/PI)*(((ADV/(N-1))/2)/PI))
!C......................................NORMALIZATION TO [0,1]
      ADV1=-ADV1
      IF(ADV1.LT.0) ADV1=1+ADV1
      ADV2=-ADV2
      IF(ADV2.LT.0) ADV2=1+ADV2
!C...........................................FINDS MIN AND MAX
      ADVMIN=(MIN(ADV1,ADV2)/2)/PI
      ADVMAX=(MAX(ADV1,ADV2)/2)/PI
!C......................................NORMALIZATION TO [0,1]
      ADV=-ADV
      IF(ADV.LT.0) THEN
        TUNEAPA=1+((ADV/(N-1))/2)/PI
      ELSE
        TUNEAPA=((ADV/(N-1))/2)/PI
      ENDIF
!C.............................................................
      END FUNCTION
!CDECK  ID>, TUNEFFT.
!C=============================================================
!C COMPUTES THE TUNE AS THE FFT ON A TWO DIMENSIONAL PLANE,
!C GIVEN N ITERATES OF A MAP. THE FFT IS PERFORMED OVER THE
!C MAXIMUM MFT WHICH SATIFIES 2**MFT.LE.N
!C THE MAXIMUM NUMBER OF ITERATES IS FIXED IN A PARAMETER.
!C TUNEFFT IS IN [0,1].
!C X IS AN ARRAY CONTAINING THE FIRST COMPONENT OF THE ITERATES
!C P IS AN ARRAY CONTAINING THE SECOND COMPONENT OF THE ITERATES
!C N IS THE TOTAL NUMBER OF ITERATES USED IN COMPUTATION
!C
!C AUTHOR:    E. TODESCO - UNIVERSITY OF BOLOGNA
!C

      REAL(KIND=fPrec) FUNCTION TUNEFFT(X,P,N)
!C............................................................
      use crcoall
      use numerical_constants, only : zero
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER N,M,NPOINT,I,NPMIN,NPMAX
      REAL(KIND=fPrec) SUM,AMAX
      REAL(KIND=fPrec) X(*),P(*)
      COMPLEX(kind=fPrec)  Z(MAXITER) ! Temp Z for CFFT, used to be SINGle precission

      TUNEFFT = ZERO

!C..................................................CHECK OF N
      IF(N.GT.MAXITER) THEN
        write(lerr,'(a)') '***ERROR(TUNEFFT): TOO MANY ITERATES'
        call prror(-1)
      ENDIF
!C............................................................
      IF (N.LE.0) THEN
        write(lerr,'(a)') '***ERROR(TUNEFFT): THIRD PARAMETER OUT OF BOUNDS'
        call prror(-1)
      ENDIF
!C............................................COMPUTATION OF M
      DO M=1,50
        IF(2**M.GT.N) GOTO 10
      ENDDO
 10   M=M-1
!C............................................................
      NPOINT=2**M
      SUM=ZERO
!C..................................SWITCH TO COMPLEX NOTATION
      DO I=1,NPOINT
        Z(I)=X(I)+(ZERO,ONE)*P(I)
        SUM=SUM+P(I)
      ENDDO
!C...........................................FOURIER TRANSFORM
      CALL CFFT(Z,-M)
!C.........COMPUTATION OF THE MAXIMUM OF THE FOURIER TRANSFORM
      AMAX=ZERO
      NPMIN=1
      IF (SUM.EQ.ZERO) THEN
        NPMAX=NPOINT/2     !..REAL FFT ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
!C.............................REJECTS FREQUENCIES NEAR 0 OR 1
      DO I=NPMIN,NPMAX
        IF(ABS(Z(I)).GT.AMAX) THEN
          TUNEFFT=FLOAT(I-1)/NPOINT
          AMAX=ABS(Z(I))
        ENDIF
      ENDDO
!C......................................NORMALIZATION TO [0,1]
      TUNEFFT=-TUNEFFT
      IF(TUNEFFT.LE.0) THEN
        TUNEFFT=1+TUNEFFT
      ENDIF
!C............................................................
      END FUNCTION
!CDECK  ID>, TUNEFFTI.
!C=============================================================
!C COMPUTES THE TUNE AS THE FFT ON A TWO DIMENSIONAL PLANE,
!C GIVEN N ITERATES OF A MAP. THE FFT IS PERFORMED OVER THE
!C MAXIMUM MFT WHICH SATIFIES 2**MFT.LE.N. THEN, THE FFT IS
!C INTERPOLATED FITTING THE THREE POINTS AROUND THE MAXIMUM
!C USING A GAUSSIAN. THE TUNE IS COMPUTED AS THE MAXIMUM OF
!C THE GAUSSIAN. THE MAXIMUM NUMBER OF ITERATES IS FIXED IN A
!C PARAMETER. THE FUNCTION IS IN [0,1].
!C X IS AN ARRAY CONTAINING THE FIRST COMPONENT OF THE ITERATES
!C P IS AN ARRAY CONTAINING THE SECOND COMPONENT OF THE ITERATES
!C N IS THE TOTAL NUMBER OF ITERATES USED IN COMPUTATION
!C
!C AUTHOR:     E. TODESCO - UNIVERSITY OF BOLOGNA
!C

      REAL(KIND=fPrec) FUNCTION TUNEFFTI(X,P,N)
!C............................................................
      use crcoall
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER N,M,NPOINT,NPMAX,NPMIN,ITUNE,I
      REAL(KIND=fPrec) SUM,AMAX,X1,X2,X3,Y1,Y2,Y3,X12,X13,Y12,Y13,X212,X213,A,B
      REAL(KIND=fPrec) X(*),P(*)
      COMPLEX(kind=fPrec) Z(MAXITER)  ! Temp Z for CFFT, used to be SINGle precission

      ITUNE = 0

!C..................................................CHECK OF N
      IF(N.GT.MAXITER) THEN
        write(lerr,'(a)') '***ERROR(TUNEFFTI): TOO MANY ITERATES'
        call prror(-1)
      ENDIF
!C............................................................
      IF (N.LE.0) THEN
        write(lerr,'(a)') '***ERROR(TUNEFFTI): THIRD PARAMETER OUT OF BOUNDS'
        call prror(-1)
      ENDIF
!C............................................COMPUTATION OF M
      DO M=1,50
        IF(2**M.GT.N) GOTO 10
      ENDDO
 10   M=M-1
!C............................................................
      NPOINT=2**M
      SUM=ZERO
!C..................................SWITCH TO COMPLEX NOTATION
      DO I=1,NPOINT
        Z(I)=X(I)+(ZERO,ONE)*P(I)
        SUM=SUM+P(I)
      ENDDO
!C...........................................FOURIER TRANSFORM
      CALL CFFT(Z,-M)
!C.........COMPUTATION OF THE MAXIMUM OF THE FOURIER TRANSFORM
      AMAX=0
      NPMIN=1
      IF (SUM.EQ.ZERO) THEN
        NPMAX=NPOINT/2     !..REAL FFT ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
!C.............................REJECTS FREQUENCIES NEAR 0 OR 1
      DO I=NPMIN,NPMAX
        IF(ABS(Z(I)).GT.AMAX) THEN
          ITUNE=I
          AMAX=ABS(Z(I))
        ENDIF
      ENDDO

      ! VKBO Bugfix for Debug build type
      ! Make sure itune is greater than 1, or we may get a segfault in the next lines
      if(itune <= 1) then
        tuneffti = one
        return
      end if
!C..............................EVALUATION OF THE NEARBY PEAKS
      X1=ITUNE-1
      X2=ITUNE
      X3=ITUNE+1
      Y1=ABS(Z(ITUNE-1))
      Y2=ABS(Z(ITUNE))
      Y3=ABS(Z(ITUNE+1))
!C...............................INTERPOLATION WITH A GAUSSIAN
      X12=X1-X2
      X13=X1-X3
      Y12=LOG_MB(REAL(Y1/Y2,fPrec))
      Y13=LOG_MB(REAL(Y1/Y3,fPrec))
      X212=X1*X1-X2*X2
      X213=X1*X1-X3*X3
!C........COMPUTATION OF THE POSITION OF THE INTERPOLATED PEAK
      A=X212*Y13-X213*Y12
      B=X12*Y13-X13*Y12
      TUNEFFTI=((A/2)/B-1)/NPOINT
!C......................................NORMALIZATION TO [0,1]
      TUNEFFTI=-TUNEFFTI
      IF(TUNEFFTI.LE.0) THEN
        TUNEFFTI=1+TUNEFFTI
      ENDIF
!C............................................................
      END FUNCTION

!CDECK  ID>, TUNELASK.
!C=============================================================
!C COMPUTE THE TUNE OF A 2D MAP BY MEANS OF LASKAR METHOD.
!C A FIRST INDICATION OF THE POSITION OF THE TUNE IS OBTAINED
!C BY MEANS OF A FFT. REFINEMENT IS OBTAINED THROUGH A NEWTON
!C PROCEDURE. THE MAXIMUM NUMBER OF ITERATIONS IS FIXED IN A
!C PARAMETER. THE FUNCTION IS IN [0,1].
!C X IS AN ARRAY CONTAINING THE ITERATION OF THE X COORDINATE OF
!C   THE MAP
!C PX IS AN ARRAY CONTAINING THE ITERATION OF THE PX COORDINATE
!C   OF THE MAP
!C MAX IS THE NUMBER OF ITERATIONS OF THE MAP
!C
!C AUTHOR:     R. BARTOLINI - BOLOGNA UNIVERSITY.
!C             M. GIOVANNOZZI - CERN HAS INTRODUCED SIGNIFICANT
!C             MODFICATIONS
!C

      REAL(KIND=fPrec) FUNCTION TUNELASK(X,PX,MAX)
!C............................................................
      use crcoall
      use numerical_constants, only : zero
      IMPLICIT NONE
      INTEGER MAXITER
      PARAMETER(MAXITER=100000)
      INTEGER MAX,NPOINT,MF,NFTMAX,NPMIN,NPMAX,NFT,JITER,JMAX,MAX1,JIT,I,K,JOM,J,N,MBODE,MFT
      REAL(KIND=fPrec) DUEPI,SUM,FTMAX,TUNEFOU,FSIGMA,OMEMIN,STEP,OMEMAX,FOMEGA,OME,ABSFOM,TMPR,TMPI
      REAL(KIND=fPrec) X(*),PX(*)
      COMPLEX(kind=fPrec) ZSING(MAXITER) ! Temp Z for CFFT, used to be SINGle precission
      COMPLEX(kind=fPrec) Z(MAXITER),FOME,ZC,SD,SP

      DUEPI=EIGHT*ATAN_MB(ONE)
      TUNELASK = zero
      JOM = 0
!C...............................CHECK OF THE ITERATION NUMBER
      IF(MAX.GT.MAXITER) THEN
        write(lerr,'(a)') '***ERROR(TUNELASK): TOO MANY ITERATIONS'
        call prror(-1)
      ENDIF
!C............................................................
      IF (MAX.LE.0) THEN
        write(lerr,'(a)') '***ERROR(TUNELASK): THIRD PARAMETER OUT OF BOUNDS'
        call prror(-1)
      ENDIF
!C.................................ESTIMATION OF TUNE WITH FFT
      SUM=ZERO
      MFT=INT(LOG_MB(REAL(MAX,fPrec))/LOG_MB(TWO))
      NPOINT=2**MFT
      DO MF=1,NPOINT
        ZSING(MF)=X(MF)-(ZERO,ONE)*PX(MF)
        SUM=SUM+PX(MF)
      ENDDO
      CALL CFFT(ZSING,-MFT)
!C......................SEARCH FOR MAXIMUM OF FOURIER SPECTRUM
      FTMAX=ZERO
      NFTMAX=ZERO
!C............................................................
      NPMIN=1
      IF (SUM.EQ.ZERO) THEN
        NPMAX=NPOINT/2 !..REAL FFT ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
!C.............................REJECTS FREQUENCIES NEAR 0 OR 1
      DO NFT=NPMIN,NPMAX
        IF (ABS(ZSING(NFT)).GT.FTMAX) THEN
          FTMAX=ABS(ZSING(NFT))
          NFTMAX=NFT
        ENDIF
      ENDDO
      TUNEFOU=FLOAT(NFTMAX-1)/NPOINT
!C..............DEFINITION OF THE FIRST INTERVAL FOR BISECTION
!C..............PROCEDURE AND INIZIALIZATION OF BISECTION
!C..............PROCEDURE: THE INTERVAL (OMEMIN,OMEMAX) IS
!C..............DIVIDED IN JMAX-1 SUBINTERVAL AND THIS
!C..............SUBDIVISION IS ITERATED JITER TIMES
      FSIGMA=ONE
      OMEMIN=TUNEFOU*DUEPI-(FSIGMA*DUEPI)/NPOINT
      OMEMAX=TUNEFOU*DUEPI+(FSIGMA*DUEPI)/NPOINT
!C.....JITER=8 PROVIDES A PRECISION WHICH IS 1.5E-4 * 1/2**MFT
      JITER=8
!C................JMAX=7 IS THE VALUE WHICH MINIMIZES CPU TIME
      JMAX=7
      MAX1=MAX-1
      STEP=(DUEPI*HALF)/MAX1
!C.........................................BISECTION PROCEDURE
      DO JIT=1,JITER
        FOMEGA=ZERO
        DO J=1,JMAX
          OME=OMEMIN+((OMEMAX-OMEMIN)/(JMAX-ONE))*(J-ONE)
          DO N=1,MAX
#ifdef CRLIBM
            ZC=(X(N)-(ZERO,ONE)*PX(N))*(ONE+COS_MB(STEP*(2*N-MAX1)))
            TMPR=REAL((-(ZERO,ONE)*OME)*N)
            TMPI=AIMAG((-(ZERO,ONE)*OME)*N)
            !exp_mb is only defined for real numbers -> decompose in real and imaginary part
            Z(N)=ZC*(EXP_MB(TMPR)*CMPLX(COS_MB(TMPI),SIN_MB(TMPI),fPrec))
#else
            ZC=(X(N)-(ZERO,ONE)*PX(N))*(ONE+COS(STEP*(2*N-MAX1)))
            Z(N)=ZC*EXP(-(ZERO,ONE)*OME*N)
#endif
          ENDDO
!C..COMPUTATION OF SCALAR PRODUCT WITH ITERATED BODE ALGORITHM
          FOME=(ZERO,ZERO)
          MBODE=(MAX-5)/4
          DO I=0,MBODE
            K=4*I
            FOME=FOME+((SEVEN*(Z(1+K)+Z(5+K))+32.0_fPrec*(Z(2+K)+Z(4+K)))+12.0_fPrec*Z(3+K))
          ENDDO
          FOME=((HALF*FOME)/45.0_fPrec)/REAL(MBODE+1,fPrec)
!C..........SEARCH FOR MAXIMUM OF SCALAR PRODUCT AND DEFINITION
!C..........OF THE NEW INTERVAL (OMEMIN,OMEMAX) WHERE TO
!C........................................RESTART THE PROCEDURE
          ABSFOM=ABS(FOME)
          IF (ABSFOM.GT.FOMEGA) THEN
            FOMEGA=ABSFOM
            TUNELASK=-OME/DUEPI
            JOM=J
          ENDIF
        ENDDO
        OMEMIN=OMEMIN+((OMEMAX-OMEMIN)/(JMAX-ONE))*(JOM-TWO)
        OMEMAX=OMEMIN+((OMEMAX-OMEMIN)/(JMAX-ONE))*JOM
      ENDDO
!C......................................NORMALIZATION TO [0,1]
      TUNELASK=-TUNELASK
      IF(TUNELASK.LE.0) THEN
        TUNELASK=1+TUNELASK
      ENDIF
!C............................................................
      END FUNCTION

end module platoFMA
