CDECK  ID>, EDTENG.
C======================================================================
C SYMPLECTIC ROTATION PARAMETRIZATION FOR A SYMPLECTIC MATRIX R
C FROM EDWARDS AND TENG. PHI = COUPLING ANGLE.
C PARA(1:7) = (TUNE1,BETA1,APLHA1,TUNE2,BETA2,ALPHA2,PHI IN DEGREE)
C PARA(11:15) = (D11,D12,D21,D22,PHI IN RADIANS)
C
 
      SUBROUTINE EDTENG(R,PARA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(TWOPI=6.28318530717958623D0,PI=0.5D0)
      DIMENSION R(4,4),PARA(20)
C
      TRMN=0.5D0*(R(1,1)-R(3,3)+R(2,2)-R(4,4))
      IF (ABS(TRMN).LT.1D-7) THEN
        WRITE(6,*) '***ERROR(EDTENG): EQUAL TUNES'
        PHI=0D0
        D11=0D0
        D12=0D0
        D21=0D0
        D22=0D0
        GOTO 2000
      ENDIF
C
      DETM=R(3,1)*R(4,2)-R(3,2)*R(4,1)
      TR=R(1,3)*R(3,1)+R(1,4)*R(4,1)+R(2,3)*R(3,2)+R(2,4)*R(4,2)
      CCMU=SQRT(TRMN*TRMN+2D0*DETM+TR)
      IF (ABS(DETM).LT.1D-10.OR.ABS(TR).LT.1D-10) CCMU=TRMN
C
      QQ=TRMN/CCMU
C     WRITE(6,*) 'trmn =',TRMN,CCMU,QQ
      IF (ABS(QQ).GT.1D0) QQ=SIGN(1D0,QQ)
      PHI=0.5D0*ACOS(QQ)      !...COUPLING ANGLE
C     IF (ABS(CCMU-TRMN).LT.1D-7.OR.ABS(CCMU+TRMN).LT.1D-7) PHI=0D0
      DENOM=CCMU*SIN(2D0*PHI)
      IF (ABS(DENOM).GT.1D-10) THEN
        D11=-(R(3,1)+R(2,4))/DENOM
        D12=-(R(3,2)-R(1,4))/DENOM
        D21=-(R(4,1)-R(2,3))/DENOM
        D22=-(R(4,2)+R(1,3))/DENOM
      ELSE
        D11=0D0
        D12=0D0
        D21=0D0
        D22=0D0
      ENDIF
C
 2000 CONTINUE
      A11=R(1,1)-(D22*R(3,1)-D12*R(4,1))*TAN(PHI)
      A12=R(1,2)-(D22*R(3,2)-D12*R(4,2))*TAN(PHI)
      A21=R(2,1)-(D11*R(4,1)-D21*R(3,1))*TAN(PHI)
      A22=R(2,2)-(D11*R(4,2)-D21*R(3,2))*TAN(PHI)
C
      B11=R(3,3)+(D11*R(1,3)+D12*R(2,3))*TAN(PHI)
      B12=R(3,4)+(D11*R(1,4)+D12*R(2,4))*TAN(PHI)
      B21=R(4,3)+(D21*R(1,3)+D22*R(2,3))*TAN(PHI)
      B22=R(4,4)+(D21*R(1,4)+D22*R(2,4))*TAN(PHI)
C
C      WRITE(6,'(A,4F16.6)') ' & ',R(1,1),R(1,2),R(2,1),R(2,2)
C      WRITE(6,'(A,4F16.6)') '   ',A11,A12,A21,A22
      AUX=0.5D0*(A11+A22)
      IF (ABS(AUX).LE.1D0) THEN
        PHASE1=ACOS(0.5D0*(A11+A22))
        BETA1=A12/SIN(PHASE1)
        IF (BETA1.LT.0D0) THEN
          PHASE1=-PHASE1
          BETA1=A12/SIN(PHASE1)
        ENDIF
        ALPHA1=0.5D0*(A11-A22)/SIN(PHASE1)
        TUNE1=PHASE1/TWOPI
        IF (TUNE1.LT.0D0) TUNE1=1D0+TUNE1
      ELSE
C       WRITE(6,*) '***MESSAGE(EDTENG): MODE 1 UNSTABLE'
        PHASE1=0D0
        BETA1=0D0
        ALPHA1=0D0
        TUNE1=0D0
      ENDIF
      IF (BETA1.LE.0D0) WRITE(6,'(A,G15.6)')
     .  ' ***ERROR(EDTENG): NEGATIVE BETA1 = ',BETA1
C
      AUX=0.5D0*(B11+B22)
      IF (ABS(AUX).LE.1D0) THEN
        PHASE2=ACOS(0.5D0*(B11+B22))
        BETA2=B12/SIN(PHASE2)
        IF (BETA2.LT.0D0) THEN
          PHASE2=-PHASE2
          BETA2=B12/SIN(PHASE2)
        ENDIF
        ALPHA2=0.5D0*(B11-B22)/SIN(PHASE2)
        TUNE2=PHASE2/TWOPI
        IF (TUNE2.LT.0D0) TUNE2=1D0+TUNE2
      ELSE
C       WRITE(6,*) '***MESSAGE(EDTENG): MODE 2 UNSTABLE'
        PHASE2=0D0
        BETA2=0D0
        ALPHA2=0D0
        TUNE2=0D0
      ENDIF
      IF (BETA2.LE.0D0) WRITE(6,'(A,G15.6)')
     .  ' ***ERROR(EDTENG): NEGATIVE BETA2 = ',BETA2
C
      PARA(1)=TUNE1
      PARA(2)=BETA1
      PARA(3)=ALPHA1
      PARA(4)=TUNE2
      PARA(5)=BETA2
      PARA(6)=ALPHA2
      PARA(7)=PHI*57.295779513082D0       !..TO DEGREE
C
      PARA(11)=D11
      PARA(12)=D12
      PARA(13)=D21
      PARA(14)=D22
      PARA(15)=PHI
C
      RETURN
      END
CDECK  ID>, GAUSSJ.
C=================================================================
C GAUSS-JORDAN SOLVER FOR LINEAR EQUATIONS FROM NUMERICAL RECIPES.
C ON EXIT A CONTAINS THE INVERSE AND B THE SOLUTION VECTOR.
C
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=1000)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      IF (NP.GT.NMAX) STOP '***ERROR(GAUSSJ): NP.GT.NMAX'
      IF (N.GT.NMAX) STOP '***ERROR(GAUSSJ): N.GT.NMAX'
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSEIF (IPIV(K).GT.1) THEN
                WRITE(6,*) '***ERROR(GAUSSJ): Singular Matrix'
                DO 26 I1=1,N
                  DO 26 I2=1,N
                    A(I1,I2)=0.0
26                CONTINUE
                A(1,1)=1.234567890
                RETURN
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) THEN
          WRITE(6,*) '***ERROR(GAUSSJ): Singular Matrix'
          DO 25 I1=1,N
            DO 25 I2=1,N
              A(I1,I2)=0.0
25        CONTINUE
          A(1,1)=1.234567890
          RETURN
        ENDIF
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
CDECK  ID>, INPS.
C=====================================================================
C GIVEN A FULL-TURN-MAP R THIS ROUTINE CALCULATES THE SYMPLECTIC
C TRANSFER MATRIX S THAT MAPS A STATE VECTOR X,X',Y,Y' INTO NORMALIZED
C PHASE SPACE.
C USES SUBROUTINE EDTENG
C
 
      SUBROUTINE INPS(R,S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TWOPI=6.28318530717958623D0)
      LOGICAL DEBUG
      DIMENSION R(4,4),S(4,4),T(4,4),Q(4,4),P20(20)
      DIMENSION SINV(4,4)
      DATA DEBUG/.FALSE./
C
      IF (DEBUG)
     .   WRITE(6,*) '***MESSAGE(INPS): INTO NORMALIZED PHASE SPACE'
      CALL EDTENG(R,P20)
      IF (DEBUG) THEN
        WRITE(6,'(3F18.8,/,3F18.8)') (P20(I),I=1,6)
        WRITE(6,'(F12.5,/,4G16.6)') P20(7),(P20(I),I=11,14)
      ENDIF
      IF (P20(2).LE.0D0.OR.P20(5).LE.0D0) THEN
        WRITE(6,*) '***ERROR(INPS): NEGATIVE BETAS RETURNED FROM EDTENG'
        WRITE(6,'(5G14.6)') (P20(I),I=1,20)
        CALL UNITM(S)
        RETURN
      ENDIF
C..................................FIRST ASSEMBLE THE NORMALIZING PART
      CALL UNITM(Q)
      SQBX=SQRT(P20(2))
      AX=P20(3)
      SQBY=SQRT(P20(5))
      AY=P20(6)
      Q(1,1)=1D0/SQBX
      Q(2,1)=AX/SQBX
      Q(2,2)=SQBX
      Q(3,3)=1D0/SQBY
      Q(4,3)=AY/SQBY
      Q(4,4)=SQBY
      IF (DEBUG) THEN
        WRITE(6,*) ' NORMALIZING MATRIX = '
        WRITE(6,'(4F18.8)') ((Q(I,J),J=1,4),I=1,4)
      ENDIF
C..................................SECOND ASSEMBLE THE COUPLING MATRIX
      CALL UNITM(T)
      IF (ABS(P20(15)).GE.1D-20) THEN
        CPHI=COS(P20(15))
        SPHI=SIN(P20(15))
        T(1,1)=CPHI
        T(2,2)=CPHI
        T(3,3)=CPHI
        T(4,4)=CPHI
        T(3,1)=P20(11)*SPHI
        T(3,2)=P20(12)*SPHI
        T(4,1)=P20(13)*SPHI
        T(4,2)=P20(14)*SPHI
        T(1,3)=-P20(14)*SPHI
        T(1,4)=P20(12)*SPHI
        T(2,3)=P20(13)*SPHI
        T(2,4)=-P20(11)*SPHI
        IF (DEBUG) THEN
          WRITE(6,*) ' DECOUPLING MATRIX T = '
          WRITE(6,'(4F18.8)') ((T(I,J),J=1,4),I=1,4)
        ENDIF
        CALL VMMUL(S,Q,T)
      ELSE
        CALL VMCOPY(S,Q)
      ENDIF
      IF (DEBUG) THEN
        WRITE(6,*) ' FULL MATRIX S ='
        WRITE(6,'(4F18.8)') ((S(I,J),J=1,4),I=1,4)
      ENDIF
C..................................CHECK THAT THE WHOLE THING IS RIGHT
      IF (DEBUG) THEN
        WRITE(6,*) ' ORIGINAL MATRIX R ='
        WRITE(6,'(4F18.8)') ((R(I,J),J=1,4),I=1,4)
        CALL VISYMP(SINV,S)
        CALL UNITM(Q)
        Q(1,1)=COS(TWOPI*P20(1))
        Q(1,2)=SIN(TWOPI*P20(1))
        Q(2,1)=-Q(1,2)
        Q(2,2)=Q(1,1)
        Q(3,3)=COS(TWOPI*P20(4))
        Q(3,4)=SIN(TWOPI*P20(4))
        Q(4,3)=-Q(3,4)
        Q(4,4)=Q(3,3)
        CALL VMMUL(T,Q,S)
        CALL VMMUL(Q,SINV,T)
        WRITE(6,*) ' RECONSTRUCTED MATRIX R ='
        WRITE(6,'(4F18.8)') ((Q(I,J),J=1,4),I=1,4)
      ENDIF
C
      RETURN
      END
CDECK  ID>, LUMPINPS.
C====================================================================
C ADDS THE TRANSFORMATION OUT OF NORMALIZED PHASE SPACE TO THE FIRST
C ELEMENT TO IN VZLUMP AND INTO NORMALIZED PHASE SPACE AT THE END OF
C THE LAST TRANSFER MATRIX IN VZLUMP ONLY
C
 
      SUBROUTINE LUMPINPS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NELE=80000)
      PARAMETER(NLUMP=NELE/10)
      DIMENSION R(4,4),S(4,4),X(4)
      CHARACTER*10 ALBL
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/ACCMAT/RA(4,4,0:NELE),IMAX2,IUPD,IUSE
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
C
      CALL INPS(RA(1,1,IMAX),S)  !..S MAPS INTO NORMALIZED PHASE SPACE
C..............................................TREAT THE LAST MATRIX
      CALL VMMUL(R,S,RL(1,1,NLMAX+1))
      CALL VMCOPY(RL(1,1,NLMAX+1),R)
      CALL PROPV4(X,S,RK(1,NLMAX+1))    !..TREAT THE KICK
      DO I=1,4
        RK(I,NLMAX+1)=X(I)
      ENDDO
C.............................................TREAT THE FIRST MATRIX
C NO KICKS NEEDED HERE, BECAUSE THEY ARE AFTER THE EMERGENCE OUT OF
C NORMALIZED PHASE SPACE
C
      CALL VISYMP(R,S)
      CALL VMMUL(S,RL(1,1,1),R)
      CALL VMCOPY(RL(1,1,1),S)
C
      RETURN
      END
CDECK  ID>, LUMPIT.
C===========================================================================
C THIS ROUTINE LUMPS THE LINEAR MAP BETWEEN NON-LINEAR ELEMENTS WHICH ARE
C CHARACTERIZED BY CODES LARGER THAN LIMIT. THE LUMPED TRANSFER MATRICES
C AND KICKS ARE PUT INTO THE COMMON BLOCK VZLUMP.RL IS THE 4 BY 4 TRANSFER
C MATRIX OF A GIVEN LUMP AND RK IS THE ACCUMULATED KICK OF THAT LUMP.
C INL CONTAINS THE STARTING AND ENDING POSITIONS OF THE NON-LINEAR BLOCKS.
C AND NLMAX IS THE NUMBER OF NON-LINEAR BLOCKS IN THE BEAM LINE. NOTE THAT
C THERE ARE NLMAX+1 LUMPS WITH TRANSFER MATRICES.
C IF NLMAX2.EQ.1 ON INPUT THERE WILL BE SOME DEBUGGING OUTPUT
C ELEMENTS WITH (IA(I).GE.LIMIT) ARE TREATED AS LUMPING SEPARATORS
C ELEMENTS WITH (IA(I).EQ.IKICK) ARE TREATED AS CENTROID DISPLACEMENT
C CARDS IN THE TRANSPORT SENSE, WHERE (IKICK.EQ.7).
C
 
      SUBROUTINE LUMPIT(NLMAX2,LIMIT,IKICK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NELE=80000)
      PARAMETER(NLUMP=NELE/10)
      DIMENSION R(4,4)
      CHARACTER*10 ALBL
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/ACCMAT/RA(4,4,0:NELE),IMAX2,IUPD,IUSE
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
C
      WRITE(6,*) '    LIMIT = ',LIMIT,'    IKICK = ',IKICK
      IDEBUG=NLMAX2
      NLMAX=0                          !..COUNTER OF NON-LINEAR BLOCKS
      IC=0
 1000 CONTINUE                         !!!..REPEAT
        IC=IC+1
        IF (IC.GT.IMAX) GOTO 1999      !..LAST ELEMENT REACHED, EXIT
        IF (IA(IC).GE.LIMIT) THEN       !..SPECIAL ELEMENT FOUND
          NLMAX=NLMAX+1
          IF (NLMAX.GT.NLUMP) STOP '***ERROR(LUMPIT): TOO MANY LUMPS!'
          INL(1,NLMAX)=IC
          JC=IC
 1100     CONTINUE
            JC=JC+1
            IF (JC.GE.IMAX) THEN
              INL(2,NLMAX)=IMAX
              GOTO 1999                !..LAST ELEMENT REACHED, EXIT
            ENDIF
            IF (IA(JC).LT.LIMIT) THEN  !..NORMAL ELEMENT FOUND AGAIN
              INL(2,NLMAX)=JC-1
              IC=JC
              GOTO 1000
            ENDIF
          IF (JC.LE.IMAX) GOTO 1100
        ENDIF
      IF (IC.LE.IMAX) GOTO 1000        !!!..UNTIL
 1999 CONTINUE
      INL(1,NLMAX+1)=IMAX+1        !..FIXES END-OF-BEAMLINE BOOK KEEPING
C
      WRITE(6,'(A,I6,A)') ' ***MESSAGE(LUMPIT): ',NLMAX,' BLOCKS FOUND'
      IF (IDEBUG.EQ.1) WRITE(6,'(I4,2I7)')
     .                      (I,INL(1,I),INL(2,I),I=1,NLMAX)
C..........................................NOW GET THE TRANSFER MATRICES
      IF (NLMAX.EQ.0) THEN                 !..NO NON-LINEAR ELEMENTS
        CALL VMCOPY(RL(1,1,1),RA(1,1,IMAX))
      ELSE
        CALL VLATCAL(RL(1,1,1),1,INL(1,1),1)            !..FIRST LUMP
        DO 1 IL=2,NLMAX+1
          CALL VLATCAL(RL(1,1,IL),MIN(IMAX,INL(2,IL-1)+1),INL(1,IL)-1,1)
    1   CONTINUE
      ENDIF
C................................................GET THE INHOMOGENEITIES
      DO 2 IL=1,NLMAX+1
        J1=1                         !..SEARCH FROM END OF PREVIOUS BLOCK
        IF (IL.GE.2) J1=MIN(IMAX,INL(2,IL-1)+1 )
        J2=MAX(1,INL(1,IL)-1)        !..TO THE START OF THE CURRENT BLOCK
        IF (IDEBUG.EQ.1) WRITE(6,'(A,I4,A,I7,A,I7)')
     .                        ' LUMP ',IL,' FROM ',J1,' TO ',J2
        DO II=1,4
          RK(II,IL)=0D0
        ENDDO
        DO 3 J=J1,J2
          IF (IA(J).EQ.IKICK) THEN   !..7-CARD FOUND
            CALL VLATCAL(R,J,J2,1)
            DO II=1,4
              DO JJ=1,4
                RK(II,IL)=RK(II,IL)+R(II,JJ)*A(J,JJ)
              ENDDO
            ENDDO
          ENDIF
    3   CONTINUE
        IF (IDEBUG.EQ.1) WRITE(6,'(4G15.6,2X,G15.6)')
     .                        ((RL(I,J,IL),J=1,4),RK(I,IL),I=1,4)
    2 CONTINUE
C
      NLMAX2=NLMAX
      END
CDECK  ID>, MANYTURN.
C============================================================================
C MAPS THE VECTOR X1=(X,X',Y,Y') THROUGH THE LUMPED LATTICE AS GIVEN IN
C COMMON BLOCK VZLUMP FOR NTURN TURNS. NLFUN IS AN EXTERNALLY DECLARED
C SUBROUTINE WITH THE FORMAT "SUBROUTINE NLFUN(XIN,XOUT,IL,NDIM)" WITH
C XIN(4) AND XOUT(4). IL THE THE LUMP NUMBER AT WHOSE END THE NON-LINEAR
C BLOCK BESIDES. NLFUN IS SUPPOSED TO HAVE ACCESS TO INL IN VZLUMP AND THE
C COMMON BLOCK LATICE AND ALSO KNOW HOW TO MAP THE PARTICLE THROUGH THE
C SECTION DEFINED BY INL(1,IL) AND INL(2,IL).
C IF THE PARTICLE WITH INITIAL COORDINATES X1(4) SURVIVES NTURN TURNS
C MANYTURN RETURNS NTURN. IF THE PARTICLE IS LOST IT RETURNS THE NEGATIVE
C OF THE TURN NUMBER WHEN IT WAS LOST.
C
      FUNCTION MANYTURN(X1,NTURN,NDIM,NLFUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NELE=80000)
      PARAMETER(NLUMP=NELE/10)
      DIMENSION X1(4),X2(4)
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
      EXTERNAL NLFUN
C......................................................LOOP OVER TURNS
      DO 1 ITURN=1,NTURN
C......................................................LOOP OVER LUMPS
        DO 2 IL=1,NLMAX  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C.............................................LINEAR TRANSPORT OF LUMP
          DO I=1,NDIM
            X2(I)=RK(I,IL)
            DO J=1,NDIM
              X2(I)=X2(I)+RL(I,J,IL)*X1(J)
            ENDDO
          ENDDO
C........................................................NON-LINEAR MAP
          CALL NLFUN(X2,X1,IL,NDIM)
C................................NOW X1 CONTAINS THE ACTUAL COORDINATES
C................................................TEST FOR PARTICLE LOSS
          IF (ABS(X1(1))+ABS(X1(2))+ABS(X1(3))+ABS(X1(4)).GT.1D10) THEN
            MANYTURN=-ITURN
C           WRITE(6,*) '***MESSAGE(MANYTURN): LOST AT ',ITURN,'/',IL
            GOTO 9999   !..EXIT
          ENDIF
    2   CONTINUE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C................................................TREAT THE LAST SECTION
        DO I=1,NDIM
          X2(I)=RK(I,NLMAX+1)
          DO J=1,NDIM
            X2(I)=X2(I)+RL(I,J,NLMAX+1)*X1(J)
          ENDDO
        ENDDO
        DO I=1,NDIM
          X1(I)=X2(I)
        ENDDO
    1 CONTINUE
      MANYTURN=NTURN
C
 9999 CONTINUE
      RETURN
      END
CDECK  ID>, PROPV4.
C====================================================================
C PROPAGATES THE (X,X',Y,Y') STATE VECTOR XIN BY THE TRANSFER MATRIX R
C
 
      SUBROUTINE PROPV4(XOUT,R,XIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NP=4)
      DIMENSION XIN(NP),XOUT(NP),R(NP,NP)
C
      DO 1 I=1,NP
        XOUT(I)=0.D0
        DO 1 J=1,NP
          XOUT(I)=XOUT(I)+R(I,J)*XIN(J)
    1 CONTINUE
C
      RETURN
      END
CDECK  ID>, UNITM.
C============================================================
C RETURNS THE 4 BY 4 UNIT MATRIX AS R.
C
 
      SUBROUTINE UNITM(R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(4,4)
C
      DO 1 I=1,4
        DO 1 J=1,4
          R(I,J)=0D0
          IF (I.EQ.J) R(I,J)=1D0
    1 CONTINUE
C
      RETURN
      END
CDECK  ID>, VFRESH2.
C===================================================================
C UPDATES THE LATTICE
 
      SUBROUTINE VFRESH2(R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NELE=80000)
      DIMENSION R(4,4)
      CHARACTER*10 ALBL
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
C
      CALL VUPUSE(1,0)                 !..SET UPDATE FLAG
      CALL VLATCAL(R,1,IMAX,1)         !..GENERATE TRANSFER MATRICES
      CALL VUPUSE(0,1)                 !..SET USE FLAG
      CALL VLATCAL(R,1,IMAX,1)         !..RECALCULATE THE TRANSFER MAP
C
      RETURN
      END
CDECK  ID>, VISYMP.
C================================================================
C EFFICIENT WAY TO CALCULATE THE INVERSE OF A 4 BY 4 SYMPLECTIC
C MATRIX. FROM R. TALMAN, SPRINGER LECTURE NOTES 343.
C
 
      SUBROUTINE VISYMP(RINV,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(4,4),RINV(4,4)
C...........................................TOP LEFT SUBMATRIX
      RINV(1,1)=R(2,2)
      RINV(1,2)=-R(1,2)
      RINV(2,1)=-R(2,1)
      RINV(2,2)=R(1,1)
C..........................................TOP RIGHT SUBMATRIX
      RINV(1,3)=R(4,2)
      RINV(1,4)=-R(3,2)
      RINV(2,3)=-R(4,1)
      RINV(2,4)=R(3,1)
C........................................BOTTOM LEFT SUBMATRIX
      RINV(3,1)=R(2,4)
      RINV(3,2)=-R(1,4)
      RINV(4,1)=-R(2,3)
      RINV(4,2)=R(1,3)
C.......................................BOTTOM RIGHT SUBMATRIX
      RINV(3,3)=R(4,4)
      RINV(3,4)=-R(3,4)
      RINV(4,3)=-R(4,3)
      RINV(4,4)=R(3,3)
C
      RETURN
      END
CDECK  ID>, VLATCAL.
C==================================================================
C RETURNS THE TRANSFER MATRIX R BETWEEN FROM THE START OF ELEMENT I1
C TO THE END OF ELEMENT I2 IF I2 IS GREATER THAN I1 AND THE ZERO
C MATRIX OTHERWISE.
C   IDIR =  1 : THE NORMAL FORWARD MATRIX IS RETURNED
C   IDIR =  0 : SAME AS 1, BUT NO RA(.,.,0) CALCULATION
C   IDIR = -1 : THE MATRIX FOR PARTICLES OF THE OPPOSITE CHARGE,
C               TRAVERSING THE LATTICE IN THE OPPOSITE DIRECTION
C               IS RETURNED
C   IDIR = -2 : THE INVERSE OF THE FORWARD MATRIX IS RETURNED
C
 
      SUBROUTINE VLATCAL(R,I1,I2,IDIR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NELE=80000)
C
      DIMENSION R(4,4),TMP(4,4),Q(4,4)
      CHARACTER*10 ALBL
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/ACCMAT/RA(4,4,0:NELE),IMAX2,IUPD,IUSE
C
      IF (I1.GT.I2) THEN
        CALL ZEROM(R)
        RETURN
      ENDIF
C
      IF (IUSE.EQ.1) THEN
        IF (I1.EQ.1) THEN
          DO J=1,4
            DO K=1,4
              R(J,K)=RA(J,K,I2)
            ENDDO
          ENDDO
        ELSE
          DO J=1,4
            DO K=1,4
              R(J,K)=RA(J,K,I1-1)
              Q(J,K)=RA(J,K,I2)
            ENDDO
          ENDDO
          CALL VISYMP(TMP,R)
          CALL VMMUL(R,Q,TMP)
        ENDIF
        GOTO 7000    !..SKIP THE DRUDGY CALCULATION
      ELSEIF (IUPD.EQ.1) THEN
        IF (I1.NE.1.OR.I2.NE.IMAX2.OR.IDIR.LT.0.OR.IUSE.NE.0)
     .    STOP '***ERROR(VLATCAL): ILLEGAL ATTEMPT TO UPDATE'
        DO J=1,4
          RA(J,J,1)=1D0
          DO K=J+1,4
            RA(J,K,1)=0D0
            RA(K,J,1)=0D0
          ENDDO
        ENDDO
      ENDIF
C
      CALL UNITM(TMP)
      CALL UNITM(R)
C
      DO 1 I=I1,I2
        IDO=1
        IF (IA(I).EQ.2) THEN                      ! POLE FACE ROTATION
          IF (IA(I-1).EQ.4) THEN
            III=I-1
          ELSEIF (IA(I+1).EQ.4) THEN
            III=I+1
          ELSEIF (IA(I-2).EQ.4) THEN
            III=I-2
          ELSEIF (IA(I+2).EQ.4) THEN
            III=I+2
          ELSE
            WRITE(6,'(A5,I5)') ' I = ',I
            STOP '***ERROR(VLATCAL): ICODE=2 NEEDS ICODE=4 AS NEIGHBOR'
          ENDIF
          CALL UNITM(Q)
          ONERHO=0.0299792458D0*A(III,2)/A(III,4)
          IF (A(I,3).GE.0D0) THEN
            TPHI=TAN(A(I,2)/57.295779513082D0)
          ELSE
            TPHI=TAN(0.5D0*ONERHO*A(III,1))
          ENDIF
          Q(2,1)=TPHI*ONERHO
          Q(4,3)=-Q(2,1)
        ELSEIF (IA(I).EQ.3) THEN                  ! DRIFT
          IF (ABS(A(I,1)).GT.1.D-8) THEN
            CALL UNITM(Q)
            Q(1,2)=A(I,1)
            Q(3,4)=A(I,1)
          ELSE
            IDO=0
          ENDIF
C       ELSEIF (IA(I).EQ.4) THEN              ! SECTOR BEND
C         CALL UNITM(Q)
C         ONERHO=0.0299792458D0*A(I,2)/A(I,4)
C         PHI=ONERHO*A(I,1)
C         IF (ABS(PHI).GT.1.D-8) THEN
C           Q(1,1)=COS(PHI)
C           Q(1,2)=SIN(PHI)/ONERHO
C           Q(2,1)=-ONERHO*SIN(PHI)
C           Q(2,2)=COS(PHI)
C           Q(3,4)=A(I,1)
C         ELSE
C           Q(1,2)=A(I,1)
C           Q(3,4)=A(I,1)
C         ENDIF
        ELSEIF (IA(I).EQ.4) THEN              ! GRADIENT SECTOR BEND
          CALL UNITM(Q)
          ONERHO=0.0299792458D0*A(I,2)/A(I,4)
          YK2=-0.0299792458D0*A(I,3)/A(I,4)
          XK2=ONERHO*ONERHO-YK2
          IF (XK2.GT.1D-10) THEN
            XK=SQRT(XK2)
            PHIX=XK*A(I,1)
            Q(1,1)=COS(PHIX)
            Q(1,2)=SIN(PHIX)/XK
            Q(2,1)=-XK*SIN(PHIX)
            Q(2,2)=COS(PHIX)
            IF (YK2.GT.1D-10) THEN
              YK=SQRT(YK2)
              PHIY=YK*A(I,1)
              Q(3,3)=COS(PHIY)
              Q(3,4)=SIN(PHIY)/YK
              Q(4,3)=-YK*SIN(PHIY)
              Q(4,4)=COS(PHIY)
            ELSEIF (ABS(YK2).LE.1D-10) THEN
              Q(3,4)=A(I,1)
            ELSEIF (YK2.LT.-1D-10) THEN
              YK=SQRT(-YK2)
              PHIY=YK*A(I,1)
              Q(3,3)=COSH(PHIY)
              Q(3,4)=SINH(PHIY)/YK
              Q(4,3)=YK*SINH(PHIY)
              Q(4,4)=COSH(PHIY)
            ENDIF
          ELSEIF (ABS(XK2).LE.1D-10) THEN
            Q(1,2)=A(I,1)
            IF (YK2.GT.1D-10) THEN
              YK=SQRT(YK2)
              PHIY=YK*A(I,1)
              Q(3,3)=COS(PHIY)
              Q(3,4)=SIN(PHIY)/YK
              Q(4,3)=-YK*SIN(PHIY)
              Q(4,4)=COS(PHIY)
            ELSEIF (ABS(YK2).LE.1D-10) THEN
              Q(3,4)=A(I,1)
            ELSEIF (YK2.LT.-1D-10) THEN
              YK=SQRT(-YK2)
              PHIY=YK*A(I,1)
              Q(3,3)=COSH(PHIY)
              Q(3,4)=SINH(PHIY)/YK
              Q(4,3)=YK*SINH(PHIY)
              Q(4,4)=COSH(PHIY)
            ENDIF
          ELSEIF (XK2.LT.-1D-10) THEN
            XK=SQRT(-XK2)
            PHIX=XK*A(I,1)
            Q(1,1)=COSH(PHIX)
            Q(1,2)=SINH(PHIX)/XK
            Q(2,1)=XK*SINH(PHIX)
            Q(2,2)=COSH(PHIX)
            IF (YK2.GT.1D-10) THEN
              YK=SQRT(YK2)
              PHIY=YK*A(I,1)
              Q(3,3)=COS(PHIY)
              Q(3,4)=SIN(PHIY)/YK
              Q(4,3)=-YK*SIN(PHIY)
              Q(4,4)=COS(PHIY)
            ELSEIF (ABS(YK2).LE.1D-10) THEN
              Q(3,4)=A(I,1)
            ELSEIF (YK2.LT.-1D-10) THEN
              YK=SQRT(-YK2)
              PHIY=YK*A(I,1)
              Q(3,3)=COSH(PHIY)
              Q(3,4)=SINH(PHIY)/YK
              Q(4,3)=YK*SINH(PHIY)
              Q(4,4)=COSH(PHIY)
            ENDIF
          ENDIF
        ELSEIF (IA(I).EQ.44) THEN             ! RECTANGULAR BEND
          CALL UNITM(Q)
          ONERHO=0.0299792458D0*A(I,2)/A(I,4)
          PHI=ONERHO*A(I,1)
          XXX=ONERHO*TAN(0.5D0*PHI)
          IF (ABS(PHI).GT.1.D-8) THEN
            Q(1,1)=COS(PHI)+XXX*SIN(PHI)/ONERHO
            Q(1,2)=SIN(PHI)/ONERHO
            Q(2,1)=-ONERHO*SIN(PHI)+2D0*XXX*COS(PHI)
     .             +XXX*XXX*SIN(PHI)/ONERHO
            Q(2,2)=Q(1,1)
            Q(3,3)=1D0-XXX*A(I,1)
            Q(3,4)=A(I,1)
            Q(4,3)=-2D0*XXX+XXX*XXX*A(I,1)
            Q(4,4)=Q(3,3)
          ELSE
            Q(1,2)=A(I,1)
            Q(3,4)=A(I,1)
          ENDIF
        ELSEIF (IA(I).EQ.5) THEN              ! QUADRUPOLE
          CALL ZEROM(Q)
          XK=0.0299792458D0*A(I,2)/(A(I,3)*A(I,4))
          IF (XK.GT.1.D-8) THEN
            XKR=SQRT(XK)
            PHI=XKR*A(I,1)
            Q(1,1)=COS(PHI)
            Q(1,2)=SIN(PHI)/XKR
            Q(2,1)=-XKR*SIN(PHI)
            Q(2,2)=COS(PHI)
            Q(3,3)=COSH(PHI)
            Q(3,4)=SINH(PHI)/XKR
            Q(4,3)=XKR*SINH(PHI)
            Q(4,4)=COSH(PHI)
          ELSEIF (XK.LT.-1.D-8) THEN
            XKR=SQRT(-XK)
            PHI=XKR*A(I,1)
            Q(1,1)=COSH(PHI)
            Q(1,2)=SINH(PHI)/XKR
            Q(2,1)=XKR*SINH(PHI)
            Q(2,2)=COSH(PHI)
            Q(3,3)=COS(PHI)
            Q(3,4)=SIN(PHI)/XKR
            Q(4,3)=-XKR*SIN(PHI)
            Q(4,4)=COS(PHI)
          ELSE
            CALL UNITM(Q)
            Q(1,2)=A(I,1)
            Q(3,4)=A(I,1)
          ENDIF
        ELSEIF (IA(I).EQ.7) THEN              ! CORRECTORS
          IDO=0
        ELSEIF (IA(I).EQ.18) THEN             ! SEXTUPOLE
          IF (ABS(A(I,1)).GT.1.D-8) THEN
            CALL UNITM(Q)
            Q(1,2)=A(I,1)
            Q(3,4)=A(I,1)
          ELSE
            IDO=0
          ENDIF
        ELSEIF (IA(I).EQ.19) THEN             ! SOLENOID
          CALL ZEROM(Q)
          XK=0.5D0*0.0299792458D0*A(I,2)/A(I,4)
          IF (ABS(XK).LT.1D-10) THEN
            CALL UNITM(Q)
            Q(1,2)=A(I,1)
            Q(3,4)=A(I,1)
          ELSE
            ALPHA=XK*A(I,1)
            C=COS(ALPHA)
            S=SIN(ALPHA)
            Q(1,1)=C*C
            Q(1,2)=S*C/XK
            Q(1,3)=S*C
            Q(1,4)=S*S/XK
            Q(2,1)=-S*C*XK
            Q(2,2)=C*C
            Q(2,3)=-S*S*XK
            Q(2,4)=S*C
            Q(3,1)=-S*C
            Q(3,2)=-S*S/XK
            Q(3,3)=C*C
            Q(3,4)=S*C/XK
            Q(4,1)=S*S*XK
            Q(4,2)=-S*C
            Q(4,3)=-S*C*XK
            Q(4,4)=C*C
          ENDIF
        ELSEIF (IA(I).EQ.20) THEN             ! ROTATION
          CALL ZEROM(Q)
          ALPHA=A(I,2)/57.295779513082D0
          C=COS(ALPHA)
          S=SIN(ALPHA)
          Q(1,1)=C
          Q(1,3)=S
          Q(2,2)=C
          Q(2,4)=S
          Q(3,1)=-S
          Q(3,3)=C
          Q(4,2)=-S
          Q(4,4)=C
        ELSEIF (IA(I).EQ.31) THEN                ! RF CAVITY
          CALL UNITM(Q)
          Q(1,2)=A(I,1)
          Q(3,4)=A(I,1)
        ELSEIF (IA(I).EQ.40) THEN                ! THIN QUAD
          IF (ABS(A(I,2)).LT.1D-6) THEN
            IDO=0                                !..A(I,2)=FOCAL LENGTH
          ELSE                                   !..POSITIVE->FOCUSING
            CALL UNITM(Q)
            Q(2,1)=-1D0/A(I,2)
            Q(4,3)=1D0/A(I,2)
          ENDIF
        ELSEIF (IA(I).EQ.41) THEN                ! THIN X-QUAD ONLY
          IF (ABS(A(I,2)).LT.1D-6) THEN
            IDO=0
          ELSE
            CALL UNITM(Q)                        !..A(I,2)=-FOCAL LENGTH
            Q(2,1)=1D0/A(I,2)                    !..POSITIVE=DEFOCUSING
          ENDIF
        ELSEIF (IA(I).EQ.42) THEN                ! THIN Y-QUAD ONLY
          IF (ABS(A(I,2)).LT.1D-6) THEN
            IDO=0
          ELSE
            CALL UNITM(Q)
            Q(4,3)=1D0/A(I,2)
          ENDIF
        ELSEIF (IA(I).EQ.55) THEN              !..ROLLED QUADRUPOLE
          CALL ZEROM(Q)                    ! A(I,1)=LENGTH [m]
          ALPHA=-A(I,3)/57.295779513082D0  ! A(I,2)=GRADIENT [kG/m]
          C=COS(ALPHA)                     ! A(I,3)=ROLL ANGLE [DEGREE]
          S=SIN(ALPHA)                     ! A(I,4)=ENERGY [GeV]
          C2=C*C
          S2=S*S
          CS=C*S
          XK=0.0299792458D0*A(I,2)/A(I,4)
          IF (XK.GT.1.D-8) THEN
            IF (ABS(A(I,1)).LT.1D-8) THEN    !..THIN QUAD
              A11=1D0                        !..A(I,2)=INT. GRADIENT
              A12=0D0
              A21=-XK
              A22=1D0
              A33=1D0
              A34=0D0
              A43=XK
              A44=1D0
            ELSE                             !..THICK QUAD
              XKR=SQRT(XK)
              PHI=XKR*A(I,1)
              A11=COS(PHI)
              A12=SIN(PHI)/XKR
              A21=-XKR*SIN(PHI)
              A22=A11
              A33=COSH(PHI)
              A34=SINH(PHI)/XKR
              A43=XKR*SINH(PHI)
              A44=A33
            ENDIF
            Q(1,1)=C2*A11+S2*A33
            Q(1,2)=C2*A12+S2*A34
            Q(1,3)=CS*(A33-A11)
            Q(1,4)=CS*(A34-A12)
            Q(2,1)=C2*A21+S2*A43
            Q(2,2)=C2*A22+S2*A44
            Q(2,3)=CS*(A43-A21)
            Q(2,4)=CS*(A44-A22)
            Q(3,1)=Q(1,3)
            Q(3,2)=Q(1,4)
            Q(3,3)=S2*A11+C2*A33
            Q(3,4)=S2*A12+C2*A34
            Q(4,1)=Q(2,3)
            Q(4,2)=Q(2,4)
            Q(4,3)=S2*A21+C2*A43
            Q(4,4)=S2*A22+C2*A44
          ELSEIF (XK.LT.-1.D-8) THEN
            IF (ABS(A(I,1)).LT.1D-8) THEN     !..THIN QUAD
              A11=1D0                        !..A(I,2)=INT. GRADIENT
              A12=0D0
              A21=-XK
              A22=1D0
              A33=1D0
              A34=0D0
              A43=XK
              A44=1D0
            ELSE                             !..THICK QUAD
              XKR=SQRT(-XK)
              PHI=XKR*A(I,1)
              A11=COSH(PHI)
              A12=SINH(PHI)/XKR
              A21=XKR*SINH(PHI)
              A22=COSH(PHI)
              A33=COS(PHI)
              A34=SIN(PHI)/XKR
              A43=-XKR*SIN(PHI)
              A44=COS(PHI)
            ENDIF
            Q(1,1)=C2*A11+S2*A33
            Q(1,2)=C2*A12+S2*A34
            Q(1,3)=CS*(A33-A11)
            Q(1,4)=CS*(A34-A12)
            Q(2,1)=C2*A21+S2*A43
            Q(2,2)=C2*A22+S2*A44
            Q(2,3)=CS*(A43-A21)
            Q(2,4)=CS*(A44-A22)
            Q(3,1)=Q(1,3)
            Q(3,2)=Q(1,4)
            Q(3,3)=S2*A11+C2*A33
            Q(3,4)=S2*A12+C2*A34
            Q(4,1)=Q(2,3)
            Q(4,2)=Q(2,4)
            Q(4,3)=S2*A21+C2*A43
            Q(4,4)=S2*A22+C2*A44
          ELSE
            CALL UNITM(Q)
            Q(1,2)=A(I,1)
            Q(3,4)=A(I,1)
          ENDIF
        ELSEIF (IA(I).EQ.90) THEN     !..JUST TO BE COMPATIBLE WITH 6D
          CALL UNITM(Q)
          Q(1,1)=-1D0
          Q(2,2)=-1D0
        ELSEIF (IA(I).EQ.101) THEN              !..X BPM
          IDO=0
        ELSEIF (IA(I).EQ.102) THEN              !..Y BPM
          IDO=0
        ELSEIF (IA(I).EQ.103) THEN              !..X/Y BPM
          IDO=0
        ELSEIF (IA(I).EQ.107) THEN              !..CENTROID SHIFT
          IDO=0
        ELSEIF (IA(I).EQ.155) THEN              !..GENERAL THIN LENS QUAD
          CALL ZEROM(Q)
          ALPHA=-A(I,3)/57.295779513082D0    ! A(I,2)=FOCAL LENGTH IN METER
          C=COS(ALPHA)                       ! A(I,3)=ROLL ANGLE IN DEGREE
          S=SIN(ALPHA)
          C2=C*C
          S2=S*S
          CS=C*S
          IF (ABS(A(I,2)).LT.1D-10) THEN
            IDO=0
          ELSE
            XK=1D0/A(I,2)
            A11=1D0
            A12=0D0
            A21=-XK
            A22=1D0
            A33=1D0
            A34=0D0
            A43=XK
            A44=1D0
            Q(1,1)=C2*A11+S2*A33
            Q(1,2)=C2*A12+S2*A34
            Q(1,3)=CS*(A33-A11)
            Q(1,4)=CS*(A34-A12)
            Q(2,1)=C2*A21+S2*A43
            Q(2,2)=C2*A22+S2*A44
            Q(2,3)=CS*(A43-A21)
            Q(2,4)=CS*(A44-A22)
            Q(3,1)=Q(1,3)
            Q(3,2)=Q(1,4)
            Q(3,3)=S2*A11+C2*A33
            Q(3,4)=S2*A12+C2*A34
            Q(4,1)=Q(2,3)
            Q(4,2)=Q(2,4)
            Q(4,3)=S2*A21+C2*A43
            Q(4,4)=S2*A22+C2*A44
          ENDIF
        ELSEIF (IA(I).EQ.200) THEN    !..INTO NORMALIZED PHASE SPACE
          CALL ZEROM(Q)
          SQBX=SQRT(A(I,1))    !..SQRT(BETAX)
          AX=A(I,2)            !..ALPHAX
          SQBY=SQRT(A(I,3))    !..SQRT(BETAY)
          AY=A(I,4)            !..ALPHAY
          Q(1,1)=1D0/SQBX
          Q(2,1)=AX/SQBX
          Q(2,2)=SQBX
          Q(3,3)=1D0/SQBY
          Q(4,3)=AY/SQBY
          Q(4,4)=SQBY
        ELSEIF (IA(I).EQ.201) THEN                   !..PHASE ADVANCE
          CALL ZEROM(Q)
          PHIX=A(I,2)/57.295779513082D0     !..IN DEGREE
          PHIY=A(I,3)/57.295779513082D0     !..IN DEGREE
          Q(1,1)=COS(PHIX)
          Q(1,2)=SIN(PHIX)
          Q(2,1)=-SIN(PHIX)
          Q(2,2)=COS(PHIX)
          Q(3,3)=COS(PHIY)
          Q(3,4)=SIN(PHIY)
          Q(4,3)=-SIN(PHIY)
          Q(4,4)=COS(PHIY)
        ELSEIF (IA(I).EQ.202) THEN    !..OUT OF NORMALIZED PHASE SPACE
          CALL ZEROM(Q)
          SQBX=SQRT(A(I,1))    !..SQRT(BETAX)
          AX=A(I,2)            !..ALPHAX
          SQBY=SQRT(A(I,3))    !..SQRT(BETAY)
          AY=A(I,4)            !..ALPHAY
          Q(1,1)=SQBX
          Q(2,1)=-AX/SQBX
          Q(2,2)=1D0/SQBX
          Q(3,3)=SQBY
          Q(4,3)=-AY/SQBY
          Q(4,4)=1D0/SQBY
        ELSEIF (IA(I).EQ.203) THEN !..COUPLING INTO NORMALIZED PHASE SPACE
          CALL UNITM(Q)
          IF (ABS(A(I,1)).GE.1D-10) THEN
            CPHI=COS(A(I,1)/57.295779513082D0)
            SPHI=SIN(A(I,1)/57.295779513082D0)
            D11=A(I,2)
            D12=A(I,3)
            D21=A(I,4)
            D22=(1D0+D12*D21)/D11
            Q(1,1)=CPHI
            Q(2,2)=CPHI
            Q(3,3)=CPHI
            Q(4,4)=CPHI
            Q(3,1)=D11*SPHI
            Q(3,2)=D12*SPHI
            Q(4,1)=D21*SPHI
            Q(4,2)=D22*SPHI
            Q(1,3)=-D22*SPHI
            Q(1,4)=D12*SPHI
            Q(2,3)=D21*SPHI
            Q(2,4)=-D11*SPHI
          ENDIF
        ELSEIF (IA(I).EQ.204) THEN !..COUPLING OUT OF NORMALIZED PHASE SPACE
          CALL UNITM(Q)
          IF (ABS(A(I,1)).GE.1D-10) THEN
            CPHI=COS(A(I,1)/57.295779513082D0)
            SPHI=SIN(A(I,1)/57.295779513082D0)
            D11=A(I,2)
            D12=A(I,3)
            D21=A(I,4)
            D22=(1D0+D12*D21)/D11
            Q(1,1)=CPHI
            Q(2,2)=CPHI
            Q(3,3)=CPHI
            Q(4,4)=CPHI
            Q(3,1)=-D11*SPHI
            Q(3,2)=-D12*SPHI
            Q(4,1)=-D21*SPHI
            Q(4,2)=-D22*SPHI
            Q(1,3)=D22*SPHI
            Q(1,4)=-D12*SPHI
            Q(2,3)=-D21*SPHI
            Q(2,4)=D11*SPHI
          ENDIF
        ELSEIF (IA(I).EQ.1000) THEN    !..THIN QUADRUPOLE FROM MULTIPOLE
           CALL UNITM(Q)
           Q(2,1)=-A(I,1)
           Q(2,3)=-A(I,2)
           Q(4,1)=-A(I,2)
           Q(4,3)=A(I,1)
        ELSEIF (IA(I).LE.0.OR.IA(I).GT.1000) THEN
          IDO=0
        ELSE
          IDO=0
          WRITE(6,'(A,2I7,1X,A)')
     .      ' ***ERROR(VLATCAL): CODE NOT SUPPORTED AT ',I,IA(I),ALBL(I)
        ENDIF
C
        IF (IDO.NE.0) THEN
          DO 2 J=1,4
            DO 2 K=1,4
              R(J,K)=0.
              DO 2 L=1,4
                R(J,K)=R(J,K)+Q(J,L)*TMP(L,K)
    2     CONTINUE
          DO 3 J=1,4
            DO 3 K=1,4
              TMP(J,K)=R(J,K)
    3     CONTINUE
        ENDIF
        IF (IUPD.EQ.1) THEN        !..UPDATE THE ACCUMULATED TRANSFER
          DO J=1,4                 !..MATRICES
            DO K=1,4
              RA(J,K,I)=R(J,K)
            ENDDO
          ENDDO
        ENDIF
    1 CONTINUE
C..................................UPDATE THE INV(1-R) TRANSFER MATRIX
      IF (IUPD.EQ.1.AND.IDIR.EQ.1) THEN
        CALL VMCOPY(RA(1,1,0),RA(1,1,IMAX))
        CALL VMDIN(RA(1,1,0))
        IF (ABS(RA(1,1,0)-1.234567890).LT.1D-10.AND.
     .      ABS(RA(1,2,0)).LT.1D-10) THEN
          WRITE(6,*) '***ERROR(VLATCAL): UNABLE TO CALCULATE INV(1-R)'
        ENDIF
      ENDIF
C
C IF IDIR=-1 THE REVERSED MATRIX IS REQUESTED. IT IS DEFINED BY
C THE INVERSE AND THE COLUMNS GET A FACTOR OF (-1)**(I+J)
C IF IDIR=-2 SIMPLY CALCULATE THE INVERSE OF THE TRANSFER MATRIX
C
 7000 CONTINUE
      IF (IDIR.EQ.-1) THEN
        CALL VISYMP(TMP,R)
        DO 10 I=1,4
          DO 10 J=1,4
            R(I,J)=TMP(I,J)*(-1)**(I+J)
   10   CONTINUE
      ELSEIF (IDIR.EQ.-2) THEN
        CALL VISYMP(TMP,R)
        DO 11 I=1,4
          DO 11 J=1,4
            R(I,J)=TMP(I,J)
   11   CONTINUE
      ENDIF
C
      RETURN
      END
CDECK  ID>, VMCOPY.
C======================================================================
C COPIES THE 4 BY 4 MATRIX R INTO RCOPY
C
 
      SUBROUTINE VMCOPY(RCOPY,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  R(4,4),RCOPY(4,4)
C
      DO 1 I=1,4
        DO 1 J=1,4
          RCOPY(I,J)=R(I,J)
    1 CONTINUE
C
      RETURN
      END
CDECK  ID>, VMDIN.
C======================================================================
C CALCULATES INV(I-R) FOR 4 BY 4 MATRICES
C
      SUBROUTINE VMDIN(R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  R(4,4),TMP(4)
C
      DO 1 I=1,4
        TMP(I)=0D0
        DO 1 J=1,4
          R(I,J)=-R(I,J)
          IF (I.EQ.J) R(I,J)=1D0+R(I,J)
    1 CONTINUE
      CALL GAUSSJ(R,4,4,TMP,1,1)
C
      RETURN
      END
CDECK  ID>, VMMUL.
C======================================================================
C MATRIX MULTIPLICATION OF 4 BY 4 MATRICES, R3=R2*R1
C
 
      SUBROUTINE VMMUL(R3,R2,R1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  R1(4,4),R2(4,4),R3(4,4)
C
      DO 1 I=1,4
        DO 1 J=1,4
          R3(I,J)=0D0
          DO 1 K=1,4
            R3(I,J)=R3(I,J)+R2(I,K)*R1(K,J)
    1 CONTINUE
C
      RETURN
      END
CDECK  ID>, VUPUSE.
C=======================================================================
C SET THE UPDATE AND USE FLAGS FOR THE ACCUMULATED TRANSFER MATRIX
C CALCULATIONS.
C
 
      SUBROUTINE VUPUSE(IUPD1,IUSE1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NELE=80000)
      COMMON/ACCMAT/RA(4,4,0:NELE),IMAX2,IUPD,IUSE
      IUPD=IUPD1
      IUSE=IUSE1
      IF (IUSE.NE.0) WRITE(6,*)
     .  '***MESSAGE(VUPUSE): LATTICE SPEEDER ACTIVE'
      IF (IUSE.EQ.0) WRITE(6,*)
     .  '***MESSAGE(VUPUSE): LATTICE SPEEDER TURNED OFF'
      IF (IUPD.NE.0) WRITE(6,*)
     .  '***MESSAGE(VUPUSE): LATTICE SPEEDER UPDATE FLAG TURNED ON'
      IF (IUPD.EQ.0) WRITE(6,*)
     .  '***MESSAGE(VUPUSE): LATTICE SPEEDER UPDATE FLAG TURNED OFF'
      RETURN
      END
CDECK  ID>, WRILAT.
C===========================================================================
C WRITES OUT A LATTICE LINE I TO UNIT NU IN THE PROPER FORMAT
C
 
      SUBROUTINE WRILAT(NU,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NELE=80000)
      CHARACTER*10 ALBL
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
C
      IF (IABS(IA(I)).GE.100.AND.IABS(IA(I)).LE.1000.OR.
     .    IABS(IA(I)).GE.2000) THEN
        WRITE(NU,'(I5,4(1X,G12.6),4X,A1,A10,A1)')
     .    IA(I),(A(I,J),J=1,4),'''',ALBL(I),''''
      ELSE
        WRITE(NU,'(I5,F12.7,E18.10,G13.6,F12.3,1X,A1,A10,A1)')
     .    IA(I),(A(I,J),J=1,4),'''',ALBL(I),''''
      ENDIF
      RETURN
      END
CDECK  ID>, ZEROM.
C============================================================
C RETURNS THE 4 BY 4 ZERO MATRIX AS R
C
 
      SUBROUTINE ZEROM(R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(4,4)
C
      DO 1 I=1,4
        DO 1 J=1,4
          R(I,J)=0D0
    1 CONTINUE
C
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
      PARAMETER(MAXITER=100000)
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      COMPLEX*8 ZSING(MAXITER)
      DIMENSION X(MAXITER),XP(MAXITER)
      DIMENSION Z(MAXITER)
C.............................................................
      IF (MAXN.GT.MAXITER) THEN
        WRITE(6,*) '***ERROR(TUNENEWT1): TOO MANY ITERATIONS'
        STOP
      ENDIF
C.............................................................
C    ESTIMATION OF TUNE WITH FFT 
C.............................................................
      DUEPI=DATAN(1D0)*8D0
      MFT=INT(LOG(FLOAT(MAXN))/LOG(2D0)) 
      NPOINT=2**MFT
      MAXN2=MAXN/2
      STEP=DUEPI/MAXN
C.............................................................
      SUM=0D0
      DO MF=1,MAXN
        Z(MF)=DCMPLX(X(MF),XP(MF))
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
      TUNEFOU=DFLOAT(NFTMAX-1)/DFLOAT(NPOINT)
      DELTAT=1D0/NPOINT
      TUNE1=TUNEFOU-DELTAT
      CALL ZFUN(TUNE,Z,MAXN,TUNE1,DELTAT)
      TUNENEWT1=1D0-TUNE
C............................................................  
      RETURN 
C............................................................  
      END
CDECK  ID>, TUNSFFT.
C=============================================================
C COMPUTES THE TUNE USING FFT.
C IT WRITES THE FOURIER SPECTRUM ON UNIT IUNIT.
C X, XP ARE THE COORDINATES OF THE ORBIT AND N IS THE 
C LENGTH OF THE ORBIT.
C

      DOUBLE PRECISION FUNCTION TUNSFFT(X,XP,N,IUNIT)
C............................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(100000),XP(100000)
      COMPLEX  Z(100000)
C..................................................CHECK OF N
      IF(N.GT.100000) THEN
        WRITE(6,*) '***ERROR(TUNSFFT): TOO MANY ITERATES'
         STOP
      ENDIF
C............................................................
      IF (N.LE.0) THEN
        WRITE(6,*) '***ERROR(TUNSFFT): THIRD PARAMETER OUT OF BOUNDS'
        STOP
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
        Z(I)=X(I)+(0.,1.)*XP(I)
        SUM=SUM+XP(I)
      ENDDO
C...........................................FOURIER TRANSFORM
      CALL CFFT(Z,-M)
C.........COMPUTATION OF THE MAXIMUM OF THE FOURIER TRANSFORM
      AMAX=0
      NPMIN=DMAX1(1D0,NPOINT*.05D0) !..ZERO EXCLUDED
      IF (SUM.EQ.0D0) THEN
        NPMAX=NPOINT/2*0.95     !..REAL FFT ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT*.95        !..ONE EXCLUDED
      ENDIF
C.............................REJECTS FREQUENCIES NEAR 0 OR 1
      DO I=1,NPOINT
C............................................................
        IF (I.GE.NPMIN.AND.I.LE.NPMAX) THEN
C............................................................
          IF (ABS(Z(I)).GT.AMAX) THEN
C............................................................
            TUNSFFT=FLOAT(I-1)/NPOINT
            AMAX=ABS(Z(I))
C............................................................
          ENDIF
C............................................................
        ENDIF
C............................................................
      ENDDO
      DO I=1,NPOINT
C............................................................
        IF (IUNIT.GT.0) THEN
          FREQ=-DFLOAT(I-1)/DFLOAT(NPOINT)
          IF (FREQ.LE.0) THEN
            FREQ=1D0+FREQ
          ENDIF
          WRITE(IUNIT,*) FREQ,CABS(Z(I))/AMAX
        ENDIF
C............................................................
      ENDDO
C......................................NORMALIZATION TO [0,1]
      TUNSFFT=-TUNSFFT
      IF (TUNSFFT.LE.0) THEN
        TUNSFFT=1+TUNSFFT
      ENDIF
C............................................................
      END
CDECK  ID>, SMOOTH.
C=============================================================
C SMOOTHS A SIGNAL BY USING THE NUMERICAL RECIPES ROUTINE
C SAVGOL.
C X, XP ARE THE TWO TIME SERIES REPRESENTING THE TWO 
C COORDINATES.
C ITMAX IS THE MAXIMUM NUMBER OF TERMS IN THE TIME SERIES.
C NL, NR IS THE NUMBER OF POINTS TO THE LEFT AND TO THE RIGHT
C USED IN THE FILTER.
C ORDER IS THE ORDER OF THE FILTER
C
C AUTHOR: M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C

      SUBROUTINE SMOOTH(X,XP,ITMAX,NL,NR,NORD) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(100000),XP(100000),TEMPX(100000),TEMPXP(100000)
      DIMENSION C(100)
C.............................................................
      NP=100
      LD=0
      CALL SAVGOL(C,NP,NL,NR,LD,NORD)
C.............................................................
      DO IDAT=1+NL,ITMAX-NR
C.............................................................
        SUMX=0D0
        SUMXP=0D0
C.............................................................
        DO K=-NL,NR
          KK=MOD(NP-K,NP)+1
          SUMX=SUMX+C(KK)*X(IDAT+K)
          SUMXP=SUMXP+C(KK)*XP(IDAT+K)
        ENDDO
C.............................................................
        TEMPX(IDAT-NL)=SUMX
        TEMPXP(IDAT-NL)=SUMXP
C.............................................................
      ENDDO
C.............................................................
      ITMAX=ITMAX-NL-NR
C.............................................................
      DO IDAT=1,ITMAX
        X(IDAT)=TEMPX(IDAT)
        XP(IDAT)=TEMPXP(IDAT)
      ENDDO
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, SAVGOL.
C=================================================================

      SUBROUTINE SAVGOL(C,NP,NL,NR,LD,M)
      INTEGER LD,M,NL,NP,NR,MMAX
      DOUBLE PRECISION C(NP)
      PARAMETER (MMAX=6)
CU    USES LUBKSB,LUDCMP
      INTEGER IMJ,IPJ,J,K,KK,MM,INDX(MMAX+1)
      DOUBLE PRECISION D,FAC,SUM,A(MMAX+1,MMAX+1),B(MMAX+1)
      IF(NP.LT.NL+NR+
     *1.OR.NL.LT.0.OR.NR.LT.0.OR.LD.GT.M.OR.M.GT.MMAX.OR.NL+NR.LT.M) 
     *PAUSE 'BAD ARGS IN SAVGOL'
      DO 14 IPJ=0,2*M
        SUM=0.D0
        IF(IPJ.EQ.0)SUM=1.D0
        DO 11 K=1,NR
          SUM=SUM+ DFLOAT(K)**IPJ
 11     CONTINUE
        DO 12 K=1,NL
          SUM=SUM+ DFLOAT(-K)**IPJ
 12     CONTINUE
        MM=MIN(IPJ,2*M-IPJ)
        DO 13 IMJ=-MM,MM,2
          A(1+(IPJ+IMJ)/2,1+(IPJ-IMJ)/2)=SUM
 13     CONTINUE
 14   CONTINUE
      CALL LUDCMP(A,M+1,MMAX+1,INDX,D)
      DO 15 J=1,M+1
        B(J)=0.D0
 15   CONTINUE
      B(LD+1)=1.D0
      CALL LUBKSB(A,M+1,MMAX+1,INDX,B)
      DO 16 KK=1,NP
        C(KK)=0.D0
 16   CONTINUE
      DO 18 K=-NL,NR
        SUM=B(1)
        FAC=1.D0
        DO 17 MM=1,M
          FAC=FAC*K
          SUM=SUM+B(MM+1)*FAC
 17     CONTINUE
        KK=MOD(NP-K,NP)+1
        C(KK)=SUM
 18   CONTINUE
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software B2..
CDECK  ID>, LUBKSB.
C=================================================================

      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      INTEGER N,NP,INDX(N)
      DOUBLE PRECISION A(NP,NP),B(N)
      INTEGER I,II,J,LL
      DOUBLE PRECISION SUM
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
 11       CONTINUE
        ELSE IF (SUM.NE.0.D0) THEN
          II=I
        ENDIF
        B(I)=SUM
 12   CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        DO 13 J=I+1,N
          SUM=SUM-A(I,J)*B(J)
 13     CONTINUE
        B(I)=SUM/A(I,I)
 14   CONTINUE
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software B2..
CDECK  ID>, LUDCMP.
C=================================================================

      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      INTEGER N,NP,INDX(N),NMAX
      DOUBLE PRECISION D,A(NP,NP),TINY
      PARAMETER (NMAX=500,TINY=1.0D-20)
      INTEGER I,IMAX,J,K
      DOUBLE PRECISION AAMAX,DUM,SUM,VV(NMAX)
      D=1.D0
      DO 12 I=1,N
        AAMAX=0.D0
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
 11     CONTINUE
        IF (AAMAX.EQ.0.D0) PAUSE 'SINGULAR MATRIX IN LUDCMP'
        VV(I)=1.D0/AAMAX
 12   CONTINUE
      DO 19 J=1,N
        DO 14 I=1,J-1
          SUM=A(I,J)
          DO 13 K=1,I-1
            SUM=SUM-A(I,K)*A(K,J)
 13       CONTINUE
          A(I,J)=SUM
 14     CONTINUE
        AAMAX=0.D0
        DO 16 I=J,N
          SUM=A(I,J)
          DO 15 K=1,J-1
            SUM=SUM-A(I,K)*A(K,J)
 15       CONTINUE
          A(I,J)=SUM
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
 16     CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
 17       CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(A(J,J).EQ.0.D0)A(J,J)=TINY
        IF(J.NE.N)THEN
          DUM=1.D0/A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
 18       CONTINUE
        ENDIF
 19   CONTINUE
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software B2..
CDECK  ID>, SIGNAL1.
C=============================================================
C GENERATES A SIMPLE SIGNAL WITH TWO MAIN FREQUENCIES PLUS 
C SOME HARMONICS.
C X, XP ARE THE TWO TIME SERIES REPRESENTING THE TWO 
C COORDINATES.
C MAX IS THE MAXIMUM NUMBER OF TERMS IN THE TIME SERIES.
C OME1, OME2 ARE THE TWO FREQUENCIES.
C

      SUBROUTINE SIGNAL1(X,XP,MAX,OME1,OME2,NHARM) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(100000),XP(100000)
C.............................................................
      DUEPI=8D0*DATAN(1D0)
      DO J=1,MAX
        X(J)=0D0
        XP(J)=0D0
      ENDDO   
C.............................................................
      DO J=1,MAX
C.............................................................
        FF1=.4D0
        FF2=.1D0
C.............................................................
        IF (OME1.EQ.0D0) FF1=0D0
        IF (OME2.EQ.0D0) FF2=0D0
C.............................................................
        DO NORD=1,NHARM
C.............................................................
          ALFA1=NORD*OME1*J
          IALFA1=ALFA1
          ALFA1=ALFA1-IALFA1
          ALFA1=ALFA1*DUEPI 
          ALFA2=NORD*OME2*J
          IALFA2=ALFA2
          ALFA2=ALFA2-IALFA2
          ALFA2=ALFA2*DUEPI 
          X(J)=X(J)+FF1*DCOS(ALFA1)+FF2*DCOS(ALFA2)
          XP(J)=XP(J)-FF1*DSIN(ALFA1)-FF2*DSIN(ALFA2)
          FF1=FF1*FF1            
          FF2=FF2*FF2
C.............................................................
        ENDDO
C.............................................................
      ENDDO
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, SIGNAL2.
C=============================================================
C GENERATES A SIMPLE SIGNAL WITH ONE MAIN FREQUENCY PLUS 
C MODULATION.
C X, XP ARE THE TWO TIME SERIES REPRESENTING THE TWO 
C COORDINATES.
C MAX IS THE MAXIMUM NUMBER OF TERMS IN THE TIME SERIES.
C OME1 IS THE MAIN FREQUENCY.
C OME2 IS THE MODULATION FREQUENCY AND AMP IS THE MODULATION
C AMPLITUDE.
C

           SUBROUTINE SIGNAL2(X,XP,MAX,OME1,OME2,AMP) 
           IMPLICIT REAL*8(A-H,O-Z)
           DIMENSION X(100000),XP(100000)
C.............................................................
           DUEPI=8*DATAN(1D+0)
C.............................................................
           DO J=1,MAX
             X(J)=0D0
             XP(J)=0D0
           ENDDO   
C.............................................................
           DO J=1,MAX
             ALFA2=OME2*J
             IALFA2=ALFA2
             ALFA2=ALFA2-IALFA2
             ALFA2=ALFA2*DUEPI 
             ALFA1=OME1*J+AMP*DSIN(ALFA2)
             IALFA1=ALFA1
             ALFA1=ALFA1-IALFA1 
             ALFA1=ALFA1*DUEPI
             X(J)=DCOS(ALFA1)
             XP(J)=DSIN(ALFA1)
           ENDDO
C.............................................................
           RETURN
C.............................................................
           END
CDECK  ID>, SIGNAL3.
C=============================================================
C GENERATES A SIMPLE SIGNAL WITH ONE MAIN FREQUENCY PLUS 
C DUMPING.
C X, XP ARE THE TWO TIME SERIES REPRESENTING THE TWO 
C COORDINATES.
C MAX IS THE MAXIMUM NUMBER OF TERMS IN THE TIME SERIES.
C OME1 IS THE INITIAL FREQUENCY AND OME2 IS THE FINAL 
C FREQUENCY. DUMP IS THE DUMPING TIME.
C

           SUBROUTINE SIGNAL3(X,XP,MAX,OME1,OME2,DUMP) 
           IMPLICIT REAL*8(A-H,O-Z)
           DIMENSION X(100000),XP(100000)
C.............................................................
           DUEPI=8*DATAN(1D+0)
C.............................................................
           DO J=1,MAX
             X(J)=0D0
             XP(J)=0D0
           ENDDO   
C.............................................................
           DO J=1,MAX
             ALFA1=(OME1-OME2)*DEXP(DFLOAT(J)/DUMP)-OME1
             IALFA1=ALFA1
             ALFA1=ALFA1-IALFA1
             ALFA1=ALFA1*DUEPI 
             X(J)=DCOS(ALFA1)
             XP(J)=DSIN(ALFA1)
           ENDDO
C.............................................................
           RETURN
C.............................................................
           END
CDECK  ID>, SIXPOST2.
C=============================================================
C PERFORMS A POST PROCESSING OF SIXTRACK TRACKING RESULTS. 
C THE FILE WRITTEN ON UNIT 6 IS ANALYZED AND THE PARTICLE 
C STABILITY IS EXTRACTED. THEN THE PROGRAM READS THE BINARY 
C FILES CONTAINING THE TRACKING DATA. 
C THE SIXTRACK OUTPUT FILE IS READ FROM UNIT 3. THE BINARY 
C FILES ARE READ FROM UNITS 90, 89, ...(SEE SIXTRACK MANUAL 
C FOR MORE DETAILS).
C NPART IS THE NUMBER OF PARTICLES FOR THE RUN. 
C NTURN IS THE MAXIMUM NUMBER OF TURNS FOR THE RUN. IF NTURN
C IS BIGGER THAN 8192, THEN THE LYAPUNOV ARE COMPUTED FOR THIS
C TIME INTERVAL AND NOT ON THE FOUR WINDOWS.
C INORM IS A FLAG TO NORMALIZE COORDINATES:
C INORM=0 ===> PHYSICAL COORDINATES
C INORM=1 ===> COURANT-SNYDER COORDINATES
C IPOST IS A FLAG TO SELECT THE POSTPROCESSING:
C IPOST=0 ===> NO POSTPROCESSING
C IPOST=1 ===> PARTICLE STABILITY 
C IPOST=2 ===> TUNE OVER DIFFERENT TIME WINDOWS
C IPOST=3 ===> LYAPUNOV OVER DIFFERENT TIME WINDOWS
C ITUNE IS A FLAG FOR TUNE COMPUTATIONS (USED ONLY IF IPOST=2):
C ITUNE=1 ===> TUNELASK IS USED
C ITUNE=2 ===> TUNENEWT IS USED
C ITUNE=3 ===> TUNEABT2 IS USED
C LPART IS A FLAG TO SELECT THE PARTICLES TO BE ANALYSED:
C LPART=0 ===> THE TWINS (PARTICLES) ARE USED
C LPART=1 ===> ONLY ONE PARTICLE OUT OF THE TWO IS USED
C IPARIT IS A FLAG TO SELECT THE PARTICLES TO BE ANALYSED:
C IPARIT=0 ===> ODD PARTICLES ARE USED
C IPARIT=1 ===> EVEN PARTICLES ARE USED
C IBATCH IS A FLAG TO SELECT THE TYPE OF EXECUTION:
C IBATCH=0 ===> INTERACTIVE
C IBATCH=1 ===> BATCH
C
C     AUTHOR: M. GIOVANNOZZI - CERN
C

      SUBROUTINE SIXPOST2(NPART,NTURN,INORM,IPOST,ITUNE,LPART,
     .                    IPARIT,IBATCH)
C.............................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.............................................................
      IMPLICIT INTEGER (I-N)
C.............................................................
      PARAMETER (NDIM=200000)
C.............................................................
      DIMENSION QWC(3),CLO(2),CLOP(2),TA(4,4),BET0(2),ALF0(2)
      DIMENSION X1(NDIM),XP1(NDIM),Y1(NDIM),YP1(NDIM),Z1(NDIM),
     .          ZP1(NDIM)
      DIMENSION X2(NDIM),XP2(NDIM),Y2(NDIM),YP2(NDIM),Z2(NDIM),
     .          ZP2(NDIM)
      DIMENSION XDUMMY(NDIM),X3(NDIM),Y3(NDIM),Z3(NDIM)
      DIMENSION XP3(NDIM),YP3(NDIM),ZP3(NDIM)
      DIMENSION DQ(4,64),BLYAP(4,64)
      DIMENSION NT(4)
C.............................................................
      CHARACTER*80 SIXTIT,COMMENT,TOPTIT(4),TITLE(12)
      CHARACTER*8 CDATE,CTIME,PROGRAM
C.............................................................
      COMMON/STAB/LOSS(64),ANGLE,RHO(64)
C.............................................................
      DATA XDUMMY/NDIM*0D0/QZ1,QZ2/2*0D0/
C.............................................................
      IF (NPART.LE.0) THEN
        WRITE(6,*) '***ERROR(SIXPOST2): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (NTURN.LE.0) THEN
        WRITE(6,*) '***ERROR(SIXPOST2): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (INORM.NE.0.AND.INORM.NE.1) THEN
        WRITE(6,*) '***ERROR(SIXPOST2): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IPOST.LT.0.OR.IPOST.GT.3) THEN
        WRITE(6,*) '***ERROR(SIXPOST2): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (ITUNE.LE.0.OR.ITUNE.GE.4) THEN
        WRITE(6,*) '***ERROR(SIXPOST2): FIFTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (LPART.NE.0.AND.LPART.NE.1) THEN
        WRITE(6,*) '***ERROR(SIXPOST2): SIXTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IPARIT.NE.0.AND.IPARIT.NE.1) THEN
        WRITE(6,*) '***ERROR(SIXPOST2): EIGHTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IBATCH.NE.0.AND.IBATCH.NE.1) THEN
        WRITE(6,*) '***ERROR(SIXPOST2): NINETH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IPOST.EQ.0) THEN
        RETURN
      ELSEIF (IPOST.EQ.1) THEN
        IPRINT=0
        CALL SIXPOST1(NPART,NTURN,IPRINT,LPART,IPARIT,IBATCH) !..PART. STAB.
        GOTO 30
      ENDIF
C..........................................LOOP OVER THE FILES
      ICOUNT=0   !..COUNTER FOR PARTICLES
      LCOUNT=0   !..COUNTER FOR PARTICLES
C.............................................................
      DO JFILE=1,NPART/2 
        READ(90-JFILE+1) SIXTIT,COMMENT,CDATE,CTIME,
     .  PROGRAM,IFIPA,ILAPA,ITOPA,ICODE,NUML,QWC(1),QWC(2),QWC(3),
     .  CLO(1),CLOP(1),CLO(2),CLOP(2),DUMMY,DUMMY,
     .  DI0X,DIP0X,DI0Z,DIP0Z,DUMMY,DUMMY,
     .  TA(1,1),TA(1,2),TA(1,3),TA(1,4),DUMMY,DUMMY,
     .  TA(2,1),TA(2,2),TA(2,3),TA(2,4),DUMMY,DUMMY,
     .  TA(3,1),TA(3,2),TA(3,3),TA(3,4),DUMMY,DUMMY,
     .  TA(4,1),TA(4,2),TA(4,3),TA(4,4),DUMMY,DUMMY,
     .  DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     .  DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     .  DMMAC,DNMS,DIZU0,DNUMLR,
     .  DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     .  DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     .  DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     .  DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     .  DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     .  DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     .  DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     .  DUMMY,DUMMY,DUMMY,DUMMY
C................................COMPUTES THE LINEAR FUNCTIONS
        IF (INORM.EQ.1) THEN
          BET0(1)=TA(1,1)*TA(1,1)+TA(1,2)*TA(1,2)
          BET0X2=TA(1,3)*TA(1,3)+TA(1,4)*TA(1,4)
          ALF0(1)=-(TA(1,1)*TA(2,1)+TA(1,2)*TA(2,2))
          ALF0X2=-(TA(1,3)*TA(2,3)+TA(1,4)*TA(2,4))
          BET0(2)=TA(3,3)*TA(3,3)+TA(3,4)*TA(3,4)
          BET0Z2=TA(3,1)*TA(3,1)+TA(3,2)*TA(3,2)
          ALF0(2)=-(TA(3,3)*TA(4,3)+TA(3,4)*TA(4,4))
          ALF0Z2=-(TA(3,1)*TA(4,1)+TA(3,2)*TA(4,2))
        ENDIF
C........................................READS THE COORDINATES
        DO ITER=1,1000000
          READ(90-JFILE+1,END=10) I1,I2,A,B,C,D,E,F,G,H,
     .                          I3,A1,B1,C1,D1,E1,F1,G1,H1
C...............................................FIRST PARTICLE
          X1D=B-CLO(1)               !..HOR. CLOSED ORBIT
          XXP1=C-CLOP(1)
C.............................................................
          IF (ICODE.GE.4) THEN
            X1D=X1D-DI0X*G           !..HOR. DISPERSION
            XXP1=XXP1-DIP0X*G
          ENDIF
C.............................................................
          IF (INORM.EQ.1) XXP1=BET0(1)*XXP1+ALF0(1)*X1D
C.............................................................
          Z1D=D-CLO(2)               !..VER. CLOSED ORBIT
          ZZP1=E-CLOP(2)
C.............................................................
          IF (ICODE.GE.4) THEN
            Z1D=Z1D-DI0Z*G           !..VER. DISPERSION 
            ZZP1=ZZP1-DIP0Z*G
          ENDIF
C.............................................................
          IF (INORM.EQ.1) ZZP1=BET0(2)*ZZP1+ALF0(2)*Z1D
C.............................................................
          X1(I1+1)=X1D
          XP1(I1+1)=XXP1
          Y1(I1+1)=Z1D
          YP1(I1+1)=ZZP1
C.............................................................
          IF (ICODE.GE.4) THEN
            Z1(I1+1)=F
            ZP1(I1+1)=G
          ELSE
            Z1(I1+1)=0D0
            ZP1(I1+1)=0D0
          ENDIF  
C..............................................SECOND PARTICLE
          X1D=B1-CLO(1)              !..HOR. CLOSED ORBIT
          XXP1=C1-CLOP(1)
C.............................................................
          IF (ICODE.GE.4) THEN
            X1D=X1D-DI0X*G           !..HOR. DISPERSION
            XXP1=XXP1-DIP0X*G
           ENDIF
C.............................................................
          IF (INORM.EQ.1) XXP1=BET0(1)*XXP1+ALF0(1)*X1D
C.............................................................
          Z1D=D1-CLO(2)              !..VER. CLOSED ORBIT
          ZZP1=E1-CLOP(2)
C.............................................................
          IF (ICODE.GE.4) THEN
            Z1D=Z1D-DI0Z*G           !..VER. DISPERSION
            ZZP1=ZZP1-DIP0Z*G
          ENDIF
C.............................................................
          IF (INORM.EQ.1) ZZP1=BET0(2)*ZZP1+ALF0(2)*Z1D
C.............................................................
          X2(I1+1)=X1D
          XP2(I1+1)=XXP1
          Y2(I1+1)=Z1D
          YP2(I1+1)=ZZP1
C.............................................................
          IF (ICODE.GE.4) THEN
            Z2(I1+1)=F1
            ZP2(I1+1)=G1
          ELSE
            Z2(I1+1)=0D0
            ZP2(I1+1)=0D0
          ENDIF  
C.............................................................
        ENDDO
C.............................................................
 10     IF (IPOST.EQ.2) THEN !..TUNE ANALYSIS
C.............................................................
          IPRINT=0
          CALL SIXPOST1(NPART,NTURN,IPRINT,LPART,IPARIT,IBATCH) !..PART. STAB.
C.............................................................
          NT(1)=64
          NT(2)=256
          NT(3)=1024
          NT(4)=4096
          IWMAX=4
C.............................................................
          ICOUNT=ICOUNT+1
C.............................................................
          IF (NTURN.GT.4096) THEN
C.............................................................
            NT(1)=NTURN/2
            NT(2)=NTURN/2
            NT(3)=NTURN/2
            NT(4)=NTURN/2
            IWMAX=1
C.............................................................
          ENDIF
C.............................................................
          DO IWIN=1,IWMAX    !..FIRST PARTICLE
            IF (LOSS(ICOUNT).GE.2*NT(IWIN)) THEN      
              IF (INORM.EQ.0) THEN
                IF (ITUNE.EQ.1) THEN      !..TUNE COMPUTATIONS
C.............................................................
                  QX1=TUNELASK(X1,XDUMMY,NT(IWIN)) !..WINDOW1
                  QY1=TUNELASK(Y1,XDUMMY,NT(IWIN))
                  IF (ICODE.GE.4) QZ1=TUNELASK(Z1,XDUMMY,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.2) THEN
C.............................................................
                  QX1=TUNENEWT(X1,XDUMMY,NT(IWIN)) !..WINDOW1
                  QY1=TUNENEWT(Y1,XDUMMY,NT(IWIN))
                  IF (ICODE.GE.4) QZ1=TUNENEWT(Z1,XDUMMY,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.3) THEN
C.............................................................
                  QX1=TUNEABT2(X1,XDUMMY,NT(IWIN)) !..WINDOW1
                  QY1=TUNEABT2(Y1,XDUMMY,NT(IWIN))
                  IF (ICODE.GE.4) QZ1=TUNEABT2(Z1,XDUMMY,NT(IWIN))
C.............................................................
                ENDIF
C.............................................................
              ELSE
C.............................................................
                IF (ITUNE.EQ.1) THEN      !..TUNE COMPUTATIONS
C.............................................................
                  QX1=TUNELASK(X1,XP1,NT(IWIN)) !..WINDOW1
                  QY1=TUNELASK(Y1,YP1,NT(IWIN))
                  IF (ICODE.GE.4) QZ1=TUNELASK(Z1,ZP1,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.2) THEN
C.............................................................
                  QX1=TUNENEWT(X1,XP1,NT(IWIN)) !..WINDOW1
                  QY1=TUNENEWT(Y1,YP1,NT(IWIN))
                  IF (ICODE.GE.4) QZ1=TUNENEWT(Z1,ZP1,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.3) THEN
C.............................................................
                  QX1=TUNEABT2(X1,XP1,NT(IWIN)) !..WINDOW1
                  QY1=TUNEABT2(Y1,YP1,NT(IWIN))
                  IF (ICODE.GE.4) QZ1=TUNEABT2(Z1,ZP1,NT(IWIN))
C.............................................................
                ENDIF
C.............................................................
              ENDIF 
C.............................................................
              CALL MCOPYV(X1,NT(IWIN)+1,2*NT(IWIN),X3)     !..COPIES
              CALL MCOPYV(Y1,NT(IWIN)+1,2*NT(IWIN),Y3)
              IF (ICODE.GE.4) 
     .          CALL MCOPYV(Z1,NT(IWIN)+1,2*NT(IWIN),Z3)
C.............................................................
              IF (INORM.EQ.0) THEN
C.............................................................
                IF (ITUNE.EQ.1) THEN      !..TUNE COMPUTATIONS
C.............................................................
                  QX2=TUNELASK(X3,XDUMMY,NT(IWIN)) !..WINDOW2
                  QY2=TUNELASK(Y3,XDUMMY,NT(IWIN))
                  IF (ICODE.GE.4) QZ2=TUNELASK(Z3,XDUMMY,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.2) THEN
C.............................................................
                  QX2=TUNENEWT(X3,XDUMMY,NT(IWIN)) !..WINDOW2
                  QY2=TUNENEWT(Y3,XDUMMY,NT(IWIN))
                  IF (ICODE.GE.4) QZ2=TUNENEWT(Z3,XDUMMY,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.3) THEN
C.............................................................
                  QX2=TUNEABT2(X3,XDUMMY,NT(IWIN)) !..WINDOW2
                  QY2=TUNEABT2(Y3,XDUMMY,NT(IWIN))
                  IF (ICODE.GE.4) QZ2=TUNEABT2(Z3,XDUMMY,NT(IWIN))
C.............................................................
                ENDIF
C.............................................................
              ELSE
C.............................................................
                CALL MCOPYV(XP1,NT(IWIN)+1,2*NT(IWIN),XP3)  !..COPIES
                CALL MCOPYV(YP1,NT(IWIN)+1,2*NT(IWIN),YP3)
                IF (ICODE.GE.4) 
     .            CALL MCOPYV(ZP1,NT(IWIN)+1,2*NT(IWIN),ZP3)
C.............................................................
                IF (ITUNE.EQ.1) THEN      !..TUNE COMPUTATIONS
C.............................................................
                  QX2=TUNELASK(X3,XP3,NT(IWIN)) !..WINDOW2
                  QY2=TUNELASK(Y3,YP3,NT(IWIN))
                  IF (ICODE.GE.4) QZ2=TUNELASK(Z3,YP3,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.2) THEN
C.............................................................
                  QX2=TUNENEWT(X3,XP3,NT(IWIN)) !..WINDOW2
                  QY2=TUNENEWT(Y3,YP3,NT(IWIN))
                  IF (ICODE.GE.4) QZ2=TUNENEWT(Z3,ZP3,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.3) THEN
C.............................................................
                  QX2=TUNEABT2(X3,XP3,NT(IWIN)) !..WINDOW2
                  QY2=TUNEABT2(Y3,YP3,NT(IWIN))
                  IF (ICODE.GE.4) QZ2=TUNEABT2(Z3,ZP3,NT(IWIN))
C.............................................................
                ENDIF
C.............................................................
              ENDIF
C....................................COMPUTES DISTANCE IN TUNE
              DQX=QX2-QX1
              DQY=QY2-QY1
              DQZ=QZ2-QZ1
              DQ(IWIN,ICOUNT)=DSQRT(DQX*DQX+DQY*DQY+DQZ*DQZ)
C.............................................................
            ELSE
C.............................................................
              DQ(IWIN,ICOUNT)=1D0
C.............................................................
            ENDIF
C.............................................................
          ENDDO
C.............................................................
          ICOUNT=ICOUNT+1
C.............................................................
          DO IWIN=1,IWMAX    !..SECOND PARTICLE
            IF (LOSS(ICOUNT).GE.2*NT(IWIN)) THEN      
              IF (INORM.EQ.0) THEN
                IF (ITUNE.EQ.1) THEN      !..TUNE COMPUTATIONS
C.............................................................
                  QX1=TUNELASK(X2,XDUMMY,NT(IWIN)) !..WINDOW1
                  QY1=TUNELASK(Y2,XDUMMY,NT(IWIN))
                  IF (ICODE.GE.4) QZ1=TUNELASK(Z2,XDUMMY,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.2) THEN
C.............................................................
                  QX1=TUNENEWT(X2,XDUMMY,NT(IWIN)) !..WINDOW1
                  QY1=TUNENEWT(Y2,XDUMMY,NT(IWIN))
                  IF (ICODE.GE.4) QZ1=TUNENEWT(Z2,XDUMMY,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.3) THEN
C.............................................................
                  QX1=TUNEABT2(X2,XDUMMY,NT(IWIN)) !..WINDOW1
                  QY1=TUNEABT2(Y2,XDUMMY,NT(IWIN))
                  IF (ICODE.GE.4) QZ1=TUNEABT2(Z2,XDUMMY,NT(IWIN))
C.............................................................
                ENDIF
C.............................................................
              ELSE
C.............................................................
                IF (ITUNE.EQ.1) THEN      !..TUNE COMPUTATIONS
C.............................................................
                  QX1=TUNELASK(X2,XP2,NT(IWIN)) !..WINDOW1
                  QY1=TUNELASK(Y2,YP2,NT(IWIN))
                  IF (ICODE.GE.4) QZ1=TUNELASK(Z2,ZP2,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.2) THEN
C.............................................................
                  QX1=TUNENEWT(X2,XP2,NT(IWIN)) !..WINDOW1
                  QY1=TUNENEWT(Y2,YP2,NT(IWIN))
                  IF (ICODE.GE.4) QZ1=TUNENEWT(Z2,ZP2,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.3) THEN
C.............................................................
                  QX1=TUNEABT2(X2,XP2,NT(IWIN)) !..WINDOW1
                  QY1=TUNEABT2(Y2,YP2,NT(IWIN))
                  IF (ICODE.GE.4) QZ1=TUNEABT2(Z2,ZP2,NT(IWIN))
C.............................................................
                ENDIF
C.............................................................
              ENDIF 
C.............................................................
              CALL MCOPYV(X2,NT(IWIN)+1,2*NT(IWIN),X3)     !..COPIES
              CALL MCOPYV(Y2,NT(IWIN)+1,2*NT(IWIN),Y3)
              IF (ICODE.GE.4) 
     .          CALL MCOPYV(Z2,NT(IWIN)+1,2*NT(IWIN),Z3)
C.............................................................
              IF (INORM.EQ.0) THEN
C.............................................................
                IF (ITUNE.EQ.1) THEN      !..TUNE COMPUTATIONS
C.............................................................
                  QX2=TUNELASK(X3,XDUMMY,NT(IWIN)) !..WINDOW2
                  QY2=TUNELASK(Y3,XDUMMY,NT(IWIN))
                  IF (ICODE.GE.4) QZ2=TUNELASK(Z3,XDUMMY,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.2) THEN
C.............................................................
                  QX2=TUNENEWT(X3,XDUMMY,NT(IWIN)) !..WINDOW2
                  QY2=TUNENEWT(Y3,XDUMMY,NT(IWIN))
                  IF (ICODE.GE.4) QZ2=TUNENEWT(Z3,XDUMMY,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.3) THEN
C.............................................................
                  QX2=TUNEABT2(X3,XDUMMY,NT(IWIN)) !..WINDOW2
                  QY2=TUNEABT2(Y3,XDUMMY,NT(IWIN))
                  IF (ICODE.GE.4) QZ2=TUNEABT2(Z3,XDUMMY,NT(IWIN))
C.............................................................
                ENDIF
C.............................................................
              ELSE
C.............................................................
                CALL MCOPYV(XP2,NT(IWIN)+1,2*NT(IWIN),XP3)  !..COPIES
                CALL MCOPYV(YP2,NT(IWIN)+1,2*NT(IWIN),YP3)
                IF (ICODE.GE.4) 
     .            CALL MCOPYV(ZP2,NT(IWIN)+1,2*NT(IWIN),ZP3)
C.............................................................
                IF (ITUNE.EQ.1) THEN      !..TUNE COMPUTATIONS
C.............................................................
                  QX2=TUNELASK(X3,XP3,NT(IWIN)) !..WINDOW2
                  QY2=TUNELASK(Y3,YP3,NT(IWIN))
                  IF (ICODE.GE.4) QZ2=TUNELASK(Z3,YP3,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.2) THEN
C.............................................................
                  QX2=TUNENEWT(X3,XP3,NT(IWIN)) !..WINDOW2
                  QY2=TUNENEWT(Y3,YP3,NT(IWIN))
                  IF (ICODE.GE.4) QZ2=TUNENEWT(Z3,ZP3,NT(IWIN))
C.............................................................
                ELSEIF (ITUNE.EQ.3) THEN
C.............................................................
                  QX2=TUNEABT2(X3,XP3,NT(IWIN)) !..WINDOW2
                  QY2=TUNEABT2(Y3,YP3,NT(IWIN))
                  IF (ICODE.GE.4) QZ2=TUNEABT2(Z3,ZP3,NT(IWIN))
C.............................................................
                ENDIF
C.............................................................
              ENDIF
C....................................COMPUTES DISTANCE IN TUNE
              DQX=QX2-QX1
              DQY=QY2-QY1
              DQZ=QZ2-QZ1
              DQ(IWIN,ICOUNT)=DSQRT(DQX*DQX+DQY*DQY+DQZ*DQZ)
C.............................................................
            ELSE
C.............................................................
              DQ(IWIN,ICOUNT)=1D0
C.............................................................
            ENDIF
C.............................................................
          ENDDO
C.............................................................
        ELSEIF (IPOST.EQ.3) THEN !..LYAPUNOV ANALYSIS
C.............................................................
          IPRINT=0
          CALL SIXPOST1(NPART,NTURN,IPRINT,LPART,IPARIT,IBATCH) !..PART. STAB.
C.............................................................
          NT(1)=64
          NT(2)=256
          NT(3)=1024
          NT(4)=4096
          IWMAX=4
C.............................................................
          LCOUNT=LCOUNT+2
C.............................................................
          IF (NTURN.GT.4096) THEN
C.............................................................
            NT(1)=NTURN
            NT(2)=NTURN
            NT(3)=NTURN
            NT(4)=NTURN
            IWMAX=1
C.............................................................
          ENDIF
C.............................................................
          IWIN=0
C.............................................................
          DELTA1=DSQRT((X1(1)-X2(1))*(X1(1)-X2(1))+
     .                 (XP1(1)-XP2(1))*(XP1(1)-XP2(1))+
     .                 (Y1(1)-Y2(1))*(Y1(1)-Y2(1))+
     .                 (YP1(1)-YP2(1))*(YP1(1)-YP2(1))+
     .                 (Z1(1)-Z2(1))*(Z1(1)-Z2(1))+
     .                 (ZP1(1)-ZP2(1))*(ZP1(1)-ZP2(1))) 
C.............................................................
          CHECK=1D-3*DSQRT(X1(1)*X1(1)+XP1(1)*XP1(1)+
     .                     Y1(1)*Y1(1)+YP1(1)*YP1(1)+
     .                     Z1(1)*Z1(1)+ZP1(1)*ZP1(1))
C.............................................................
          DO ITER=1,NT(4)
C.............................................................
            DELTA=DSQRT((X1(ITER)-X2(ITER))*(X1(ITER)-X2(ITER))+
     .                  (XP1(ITER)-XP2(ITER))*(XP1(ITER)-XP2(ITER))+
     .                  (Y1(ITER)-Y2(ITER))*(Y1(ITER)-Y2(ITER))+
     .                  (YP1(ITER)-YP2(ITER))*(YP1(ITER)-YP2(ITER))+
     .                  (Z1(ITER)-Z2(ITER))*(Z1(ITER)-Z2(ITER))+
     .                  (ZP1(ITER)-ZP2(ITER))*(ZP1(ITER)-ZP2(ITER))) 
C.............................................................
            IF (DELTA.GT.CHECK) THEN !..THRESHOLD CHECK
C.............................................................
              ALYAP=1D0/DFLOAT(ITER-1)*LOG(DELTA/DELTA1)
C.............................................................
              IF (NT(1).GT.ITER) BLYAP(1,JFILE)=ALYAP
C.............................................................
              IF (NT(2).GT.ITER) BLYAP(2,JFILE)=ALYAP
C.............................................................
              IF (NT(3).GT.ITER) BLYAP(3,JFILE)=ALYAP
C.............................................................
              IF (NT(4).GT.ITER) BLYAP(4,JFILE)=ALYAP
C.............................................................
            ENDIF
C.............................................................
            IF (LOSS(LCOUNT).LE.ITER.OR.LOSS(LCOUNT-1).LE.ITER) THEN
C.............................................................
              ALYAP=1D0
C.............................................................
              IF (NT(1).GT.ITER) BLYAP(1,JFILE)=ALYAP
C.............................................................
              IF (NT(2).GT.ITER) BLYAP(2,JFILE)=ALYAP
C.............................................................
              IF (NT(3).GT.ITER) BLYAP(3,JFILE)=ALYAP
C.............................................................
              IF (NT(4).GT.ITER) BLYAP(4,JFILE)=ALYAP
C.............................................................
              GOTO 20
            ENDIF
C.............................................................
            ALYAP=1D0/DFLOAT(ITER-1)*LOG(DELTA/DELTA1)
C.............................................................
            IF (NT(1).EQ.ITER) BLYAP(1,JFILE)=ALYAP
C.............................................................
            IF (NT(2).EQ.ITER) BLYAP(2,JFILE)=ALYAP
C.............................................................
            IF (NT(3).EQ.ITER) BLYAP(3,JFILE)=ALYAP
C.............................................................
            IF (NT(4).EQ.ITER) BLYAP(4,JFILE)=ALYAP
C.............................................................
          ENDDO
C.............................................................
 20       CONTINUE
C.............................................................
        ENDIF
C.............................................................
      ENDDO
C.............................................................
 30   IF (IPOST.EQ.2) THEN
        DO IDAT=1+IPARIT,NPART,1+LPART
          WRITE(25,40) ANGLE,RHO(IDAT),LOSS(IDAT)
          WRITE(26,50) ANGLE,RHO(IDAT),(DQ(IWIN,IDAT),IWIN=1,IWMAX)
        ENDDO
      ENDIF
C.............................................................
      IF (IPOST.EQ.3) THEN
        DO IDAT=1,NPART/2
          JDAT=2*IDAT-1
          WRITE(25,40) ANGLE,RHO(JDAT),LOSS(JDAT)
          WRITE(26,50) ANGLE,RHO(JDAT),(BLYAP(IWIN,IDAT),IWIN=1,IWMAX)
        ENDDO
      ENDIF
C.............................................................
 40   FORMAT(1X,2F8.4,I12)
 50   FORMAT(1X,2F8.4,4F13.9)
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, SIXPOST1.
C=============================================================
C PERFORMS A POST PROCESSING OF THE SIXTRACK OUTPUT FILE IN 
C CASE OF LONG TERM TRACKING.
C THE FILE WRITTEN ON UNIT 6 IS ANALYZED AND THE PARTICLE 
C STABILITY IS EXTRACTED. THE SIXTRACK OUTPUT FILE IS READ 
C FROM UNIT 3. THE OUTCOME OF THE ANALYSIS IS WRITTEN ON 
C UNIT 25.
C NPART IS THE NUMBER OF PARTICLES FOR THE RUN. 
C NTURN IS THE MAXIMUM NUMBER OF TURNS FOR THE RUN.
C IPRINT IS A FLAG FOR PRINTOUT:
C IPRINT=0 ===> NO PRINTOUT
C IPRINT=1 ===> THE ARRAY LOSS IS PRINTED ON UNIT 25
C THE STABILITY TIME OF DIFFERENT PARTICLES IS ALSO  STORED IN THE
C ARRAY LOSS IN THE COMMON BLOCK /STAB/
C LPART IS A FLAG FOR THE ANALYSIS:
C LPART=0 ===> BOTH PARTICLES ARE ANALYZED
C LPART=1 ===> ONLY ONE PARTICLE IS ANALYZED
C IPARIT IS A FLAG TO SELECT WHETHER ODD OR EVEN PARTICLES ARE ANALYSED
C IPARIT=0 ===> ODD PARTICLES ARE ANALIZED
C IPARIT=1 ===> EVEN PARTICLES ARE ANALIZED
C IBATCH IS A FLAG FOR THE BATCH EXECUTION
C IBATCH=0 ===> INTERACTIVE EXECUTION. THE PROGRAM WILL ASK TO
C               OPEN A FILE AND IT WILL ANALYZE MANY FILES 
C               (EVENTUALLY).
C IBATCH=1 ===> BATCH EXECUTION
C
C     AUTHOR: M. GIOVANNOZZI - CERN
C

      SUBROUTINE SIXPOST1(NPART,NTURN,IPRINT,LPART,IPARIT,IBATCH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION X(64),PX(64),Y(64),PY(64)
      CHARACTER*80 ALINE,FILEIN
      CHARACTER*17 DUMMY0
      CHARACTER*27 DUMMY1
      CHARACTER*20 DUMMY2
      CHARACTER*103 DUMMY3
      CHARACTER*1 ANS
C.............................................................
      COMMON/STAB/LOSS(64),ANGLE,RHO(64)
C.............................................................
      IF (NPART.LE.0) THEN
        WRITE(6,*) '***ERROR(SIXPOST1): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (NTURN.LE.0) THEN
        WRITE(6,*) '***ERROR(SIXPOST1): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IPRINT.LT.0) THEN
        WRITE(6,*) '***ERROR(SIXPOST1): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (LPART.LT.0.AND.LPART.GT.1) THEN
        WRITE(6,*) '***ERROR(SIXPOST1): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IPARIT.LT.0.AND.IPARIT.GT.1) THEN
        WRITE(6,*) '***ERROR(SIXPOST1): FIFTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IBATCH.LT.0.AND.IBATCH.GT.1) THEN
        WRITE(6,*) '***ERROR(SIXPOST1): SIXTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      PI2=8D0*DATAN(1D0)
      IPARIT=LPART*IPARIT
C.............................................................
 10   IF (IBATCH.EQ.0) THEN
        WRITE(6,*) 'ENTER THE FILE NAME OF THE SIXTRACK OUTPUT ',
     .             '(< 80 CHAR)'
        READ(*,'(A)') FILEIN
        OPEN(3,FILE=FILEIN,FORM='FORMATTED',STATUS='UNKNOWN')
      ENDIF
C.............................................................
      DO IDAT=1,NPART  !..INITIALIZATION
        LOSS(IDAT)=NTURN+1
      ENDDO
C.............................................................
      KPART=0
      IBETAH=0
      IBETAV=0
      IALFAH=0
      IALFAV=0
C.............................................................
      DO ILINE=1,100000  !..READS OUTPUT FILE
        READ(3,'(A)',END=20) ALINE
        IF (INDEX(ALINE,'ABNORMALLY').NE.0) THEN !..THE PARTICLE IS LOST
C.............................................................
          READ(3,'(1X,A17,I3)',END=20) DUMMY0,IPART
          READ(3,'(1X,A27,1X,I8)',END=20) DUMMY1,LOSS(IPART)
C.............................................................
        ELSEIF (INDEX(ALINE,'HORIZONTAL BETA ').NE.0.AND.IBETAH.EQ.0) 
     .                                          THEN !..VALUE OF BETA
C.............................................................
          BACKSPACE(3)
          READ(3,'(A103,F16.10)',END=20) DUMMY3,BETAH
          WRITE(6,'(A103,F16.10)') DUMMY3,BETAH
          IBETAH=1
C.............................................................
        ELSEIF (INDEX(ALINE,'VERTICAL BETA ').NE.0.AND.IBETAV.EQ.0) 
     .                                        THEN !..VALUE OF BETA
C.............................................................
          BACKSPACE(3)
          READ(3,'(A103,F16.10)',END=20) DUMMY3,BETAV
          WRITE(6,'(A103,F16.10)') DUMMY3,BETAV
          IBETAV=1
C.............................................................
        ELSEIF (INDEX(ALINE,'HORIZONTAL ALFA ').NE.0.AND.IALFAH.EQ.0) 
     .                                          THEN !..VALUE OF ALFA
C.............................................................
          BACKSPACE(3)
          READ(3,'(A103,F16.10)',END=20) DUMMY3,ALFAH
          WRITE(6,'(A103,F16.10)') DUMMY3,ALFAH
          IALFAH=1
C.............................................................
        ELSEIF (INDEX(ALINE,'VERTICAL ALFA ').NE.0.AND.IALFAV.EQ.0) 
     .                                       THEN !..VALUE OF ALFA
C.............................................................
          BACKSPACE(3)
          READ(3,'(A103,F16.10)',END=20) DUMMY3,ALFAV
          WRITE(6,'(A103,F16.10)') DUMMY3,ALFAV
          IALFAV=1
C.............................................................
        ELSEIF(INDEX(ALINE,'NO CL.ORBIT').NE.0) THEN
C.............................................................
          KPART=KPART+1                !..FIRST PARTICLE
C.............................................................
          READ(3,'(1X,A20,F47.33)',END=20) DUMMY2,X(KPART)
          READ(3,'(1X,A20,F47.33)',END=20) DUMMY2,PX(KPART)
          READ(3,'(1X,A20,F47.33)',END=20) DUMMY2,Y(KPART)
          READ(3,'(1X,A20,F47.33)',END=20) DUMMY2,PY(KPART)
          READ(3,'(1X,A20,F47.33)',END=20) DUMMY2,SIGMA
          READ(3,'(1X,A20,F47.33)',END=20) DUMMY2,DELTAP
C.............................................................
          KPART=KPART+1                !..SECOND PARTICLE
C.............................................................
          READ(3,'(1X,A20,F47.33)',END=20) DUMMY2,X(KPART)
          READ(3,'(1X,A20,F47.33)',END=20) DUMMY2,PX(KPART)
          READ(3,'(1X,A20,F47.33)',END=20) DUMMY2,Y(KPART)
          READ(3,'(1X,A20,F47.33)',END=20) DUMMY2,PY(KPART)
C.............................................................
        ENDIF
      ENDDO   
C.............................................................
      ANGLE=0D0
 20   DO IDAT=1,NPART !..COMPUTES THE ANGLE IN NORM PHASE SPACE
C.............................................................
        IF (X(IDAT).NE.0D0.AND.BETAV.NE.0D0) THEN
C.............................................................
          ANGLE=DATAN(Y(IDAT)/X(IDAT)*DSQRT(BETAH/BETAV))/PI2
C.............................................................
        ENDIF
C.............................................................
        RHO(IDAT)=DSQRT(X(IDAT)*X(IDAT)/BETAH+Y(IDAT)*Y(IDAT)/BETAV)
C.............................................................
      ENDDO
C.............................................................
      IF (ANGLE.EQ.0D0) THEN
        WRITE(6,*) '***ERROR(SIXPOST1): X COORDINATE OR ',
     .             'HOR. BETA FUNCTION ARE ZERO!'
        RETURN
      ENDIF
C.............................................................
      IF (IPRINT.EQ.1) THEN
C.............................................................
        DO IDAT=1+IPARIT,NPART,1+LPART  !..WRITES THE RESULTS 
C.............................................................
          WRITE(25,30) ANGLE,RHO(IDAT),LOSS(IDAT)
C.............................................................
          WRITE(26,*) 'PX NORM= ',
     .     ALFAH/DSQRT(BETAH)*X(IDAT)+DSQRT(BETAH)*PX(IDAT)
          WRITE(26,*) 'PY NORM= ',
     .     ALFAV/DSQRT(BETAV)*Y(IDAT)+DSQRT(BETAV)*PY(IDAT)
C.............................................................
        ENDDO
      ENDIF
C.............................................................
      IF (IBATCH.EQ.0) THEN
        WRITE(6,*) 'DO YOU WANT TO CONTINUE (Y/N)?'
        READ(*,'(A)') ANS
C.............................................................
        IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') GOTO 10
C.............................................................
      ENDIF
C.............................................................
      CLOSE(3)
      CLOSE(25)
      CLOSE(26)
C.............................................................
      RETURN
C.............................................................
 30   FORMAT(1X,2F8.4,I12)
C.............................................................
      END
CDECK  ID>, SIXGEN.
C=============================================================
C GENERATES THE INITIAL COORDINATES FOR SIMULATIONS BASED ON
C SIXTRACK, SUCH THAT THE GRID IN THE NORMALIZED COORDINATES 
C IS UNIFORM AND ALSO THE SECOND PARTICLE BELONGS TO THE SAME
C LINE IN NORMALIZED COORDINATES.
C TWISS IS AN ARRAY CONTAINING THE TWISS PARAMETERS:
C TWISS(1)=BETAH
C TWISS(2)=ALFAH
C TWISS(3)=BETAV
C TWISS(4)=ALFAH
C RHOMIN, RHOMAX AND NRHO ARE THE MIN. RADIUS, MAX RADIUS AND 
C NUMBER OF STEPS IN THE NORMALIZED GRID.
C ANGMIN, ANGMAX AND NANG ARE THE MIN. ANGLE, MAX ANGLE AND 
C NUMBER OF STEPS IN THE NORMALIZED GRID. THE ANGLES ARE 
C EXPRESSED IN UNITS OF 2PI.
C X IS AN ARRAY CONTAINING THE MIN AND MAX AMPLITUDE TO BE
C SPECIFIED IN THE TRACKING BLOCK OF SIXTRACK.
C DELTA IS AN ARRAY CONTAINING THE DISPLACEMENTS OF THE SECOND
C PARTICLE TO BE SPECIFIED IN THE INITIAL COORDINATES BLOCK
C OF SIXTRACK.
C
C     AUTHOR: M. GIOVANNOZZI - CERN
C

      SUBROUTINE SIXGEN(TWISS,RHOMIN,RHOMAX,NRHO,ANGMIN,ANGMAX,NANG,
     .                  X,DELTA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION TWISS(4),X(2),DELTA(4)
C.............................................................
      IF (RHOMIN.GT.RHOMAX) THEN
        WRITE(6,*) '***WARNING(SIXGEN): BOUNDS INVERTED FOR INTERVAL ',
     .             'IN RHO'
        TEMP=RHOMIN
        RHOMIN=RHOMAX
        RHOMAX=TEMP
      ENDIF
C.............................................................
      IF (NRHO.LE.0) THEN
        WRITE(6,*) '***ERROR(SIXGEN): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (ANGMIN.GT.ANGMAX) THEN
        WRITE(6,*) '***WARNING(SIXGEN): BOUNDS INVERTED FOR INTERVAL ',
     .             'IN ANGLE'
        TEMP=ANGMIN
        ANGMIN=ANGMAX
        ANGMAX=TEMP
      ENDIF
C.............................................................
      IF (NANG.LE.0) THEN
        WRITE(6,*) '***ERROR(SIXGEN): SEVENTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      PI2=8D0*DATAN(1D0)
C.............................................................
      NRHO=NRHO*.5
      DANG=(ANGMAX-ANGMIN)/DFLOAT(NANG+1)
C.............................................................
      DO IANG=1,NANG                          !..LOOP ON ALL THE ANGLES
C.............................................................
        ERATIO=DTAN(PI2*(ANGMIN+DFLOAT(IANG)*DANG))*
     .         DTAN(PI2*(ANGMIN+DFLOAT(IANG)*DANG)) !..EMITTANCE RATIO
C.............................................................
        XMAX=DSQRT(TWISS(1)/(1D0+ERATIO))*RHOMAX
        XMIN=DSQRT(TWISS(1)/(1D0+ERATIO))*RHOMIN
        DAMP=(XMAX-XMIN)/DFLOAT(NRHO-1)
C.............................................................
        X(1)=XMIN                                         !..MIN AMP
        X(2)=XMAX                                         !..MAX AMP
C.............................................................
        DELTA(1)=DAMP*.5                                  !..DX
        DELTA(2)=-TWISS(2)/TWISS(1)*DELTA(1)              !..DPX
        DELTA(3)=DSQRT(TWISS(3)/TWISS(1)*ERATIO)*DELTA(1) !..DY
        DELTA(4)=-TWISS(4)/TWISS(3)*DELTA(3)              !..DPY
C.............................................................
        WRITE(26,*) 'ANGLE # ',IANG
        WRITE(26,*) 'MIN. AMP.= ',X(1)
        WRITE(26,*) 'MAX. AMP.= ',X(2)
        WRITE(26,*) 'NUMBER OF AMP. VAR.= ',NRHO
        WRITE(26,*) 'EMITTANCE RATIO= ',ERATIO
        WRITE(26,*) 'DX = ',DELTA(1)
        WRITE(26,*) 'DPX= ',DELTA(2)
        WRITE(26,*) 'DY = ',DELTA(3)
        WRITE(26,*) 'DPY= ',DELTA(4)
C.............................................................
      ENDDO
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
      PARAMETER(MAXITER=100000)
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      COMPLEX*16 ZSING(MAXITER)
      DIMENSION X(MAXITER),XP(MAXITER)
      DIMENSION Z(MAXITER)
C..................................ESTIMATION OF TUNE WITH FFT 
      PI=DATAN(1D0)*4D0
      DUEPI=2*PI
      MFT=INT(LOG(FLOAT(MAXN))/LOG(2D0)) 
      NPOINT=2**MFT
      STEP=DUEPI/NPOINT/2D+0
C.............................................................
      SUM=0D0            !..CHECKS FOR COMPLEX OR REAL DATA
      DO MF=1,NPOINT
        Z(MF)=DCMPLX(X(MF),XP(MF))
        ZSING(MF)=Z(MF)
        SUM=SUM+XP(MF)
      ENDDO 
      CALL FFT(ZSING,NPOINT,-1)
C.......................SEARCH FOR MAXIMUM OF FOURIER SPECTRUM
      NPMIN=1
      IF (SUM.EQ.0D0) THEN
        NPMAX=NPOINT/2         !..REAL FFT ONLY HALF COEFFICIENTS
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
      IF (CF3.GT.CF1) THEN      
        ASSK=DFLOAT(NFTMAX)+NPOINT/PI*
     .       ATAN2(CF3*DSIN(PI/NPOINT),CF2+CF3*DCOS(PI/NPOINT))
      ELSEIF (CF3.LE.CF1) THEN                   
        ASSK=DFLOAT(NFTMAX-1)+NPOINT/PI*
     .       ATAN2(CF2*DSIN(PI/NPOINT),CF1+CF2*DCOS(PI/NPOINT))
      ENDIF
      TUNEABT=1D+0-(ASSK-1D+0)/DFLOAT(NPOINT)
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
      PARAMETER(MAXITER=100000)
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      COMPLEX*16 ZSING(MAXITER)
      DIMENSION X(MAXITER),XP(MAXITER)
      DIMENSION Z(MAXITER)
C..................................ESTIMATION OF TUNE WITH FFT 
      PI=DATAN(1D0)*4D0
      DUEPI=2*PI
      MFT=INT(LOG(FLOAT(MAXN))/LOG(2D0)) 
      NPOINT=2**MFT
      STEP=DUEPI/NPOINT/2D+0
C.............................................................
      SUM=0D0            !..CHECKS FOR COMPLEX OR REAL DATA
      DO MF=1,NPOINT
        Z(MF)=DCMPLX(X(MF),XP(MF))*DSIN(STEP*MF)**2
        ZSING(MF)=Z(MF)
        SUM=SUM+XP(MF)
      ENDDO 
      CALL FFT(ZSING,NPOINT,-1)
C.......................SEARCH FOR MAXIMUM OF FOURIER SPECTRUM
      NPMIN=1
      IF (SUM.EQ.0D0) THEN
        NPMAX=NPOINT/2         !..REAL FFT ONLY HALF COEFFICIENTS
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
      CO=DCOS(2*PI/DFLOAT(NPOINT))
      SI=DSIN(2*PI/DFLOAT(NPOINT))
      SCRA1=CO**2*(P1+P2)**2-2*P1*P2*(2*CO**2-CO-1)       
      SCRA2=(P1+P2*CO)*(P1-P2)
      SCRA3=P1**2+P2**2+2*P1*P2*CO
      SCRA4=(-SCRA2+P2*SQRT(SCRA1))/SCRA3
      ASSK=DFLOAT(NN)+NPOINT/2/PI*ASIN(SI*SCRA4)
      TUNEABT2=1D+0-(ASSK-1D+0)/DFLOAT(NPOINT)
C............................................................  
      RETURN 
C............................................................  
      END
CDECK  ID>, FFT.
C============================================================
C           COMPUTES THE FFT   (DOUBLE PRECISION)
C           AUTHOR: NUMERICAL RECEIPES, PG. 395
C           NN IS THE NUMBER OF DATA: MUST BE A POWER OF 2
C           ISIGN=1: DIRECT FT
C           ISIGN=-1: INVERSE FT
C           DATA IS A COMPLEX ARRAY WITH THE SIGNAL IN INPUT
C                                   WITH THE FT IN OUTPUT
C           GOOD LUCK, BABY
C

      SUBROUTINE FFT(DATA,NN,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,TEMPI,TEMPR
      REAL*8 DATA(*)

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
        THETA=8.D0*DATAN(1.D0)/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
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
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
 13     CONTINUE
        MMAX=ISTEP
        GOTO 2
      END IF
C.............................................................      
      END
CDECK  ID>, MLYAP.
C============================================================
C COMPUTES THE LYAPUNOV EXPONENT FOR AN INITIAL CONDITION USING
C TWO NEARBY PARTICLES. 
C IMAP IS A PARAMETER SPECIFING THE MAP TYPE:
C IMAP = 1 ===> HENON MAP
C IMAP = 2 ===> ACCELERATOR MODEL
C NTURN IS THE TURN NUMBER USED TO COMPUTE THE STABILITY 
C (INPUT). IN OUTPUT IT REPRESENTS THE LOSS TIME IN CASE THE 
C PARTICLE IS UNSTABLE.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C X IS AN ARRAY CONTAINING THE INITIAL CONDITION.
C DELTA IS THE DISTANCE BETWEEN THE TWO PARTICLES USED TO COMPUTE
C THE LYAPUNOV EXPONENT.
C BLYAP IS AN ARRAY CONTAINING THE ESTIMATED LYAPUNOV EXPONENT 
C AS A FUNCTION OF THE ITERATION.
C IOUT IS A FLAG FOR OUTPUT:
C IOUT = 0 ===> NO OUTPUT
C IOUT = 1 ===> THE ARRAY BLYAP IS WRITTEN ON A FILE IN UNIT 31
C ISTEP IS USED TO SELECT THE OUTPUT FREQUENCY. IT IS IGNORED IN
C CASE IOUT = 0.
C
C AUTHOR : M. GIOVANNOZZI & E. TODESCO - CERN 
C

      SUBROUTINE MLYAP(IMAP,NTURN,IDIME,X,DELTA,BLYAP,IOUT,ISTEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      PARAMETER (NMAX=200000)
      DIMENSION X(4),X1(4),X2(4),BLYAP(NMAX)
C.............................................................
      EXTERNAL NLFUN
C.............................................................
      IF (IMAP.NE.1.AND.IMAP.NE.2) THEN
        WRITE(6,*) '***ERROR(MLYAP): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (NTURN.LT.0.OR.NTURN.GT.NMAX) THEN
        WRITE(6,*) '***ERROR(MLYAP): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MLYAP): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IOUT.NE.0.AND.IOUT.NE.1) THEN
        WRITE(6,*) '***ERROR(MLYAP): SEVENTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (ISTEP.LT.0) THEN
        WRITE(6,*) '***ERROR(MLYAP): EIGHTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      DO IDIM=1,4                     !..DEFINES INITIAL CONDITIONS
        X1(IDIM)=X(IDIM)
C.............................................................
        X2(IDIM)=X(IDIM)+DELTA/DSQRT(DFLOAT(IDIME))
      ENDDO
C.............................................................
      EXPAC=0D0                       !..INITIALIZATION
      LTURN=NTURN
      EPSLY=DELTA
C.............................................................
      DO ITURN=2,LTURN                !..TRACKING
        IF (IMAP.EQ.1) THEN
          KTURN1=MANY_HENON(X1,1)    !..FIRST PARTICLE
          KTURN2=MANY_HENON(X2,1)    !..COMPANION
        ELSE
          KTURN1=MANYTURN(X1,1,IDIME,NLFUN)  !..FIRST PARTICLE
          KTURN2=MANYTURN(X2,1,IDIME,NLFUN)  !..COMPANION
        ENDIF
C.............................................................
        IF (KTURN1.LT.0.OR.KTURN2.LT.0) THEN  !..CHECK PARTICLE LOSS
C.............................................................
          NTURN=ITURN
          GOTO 20
C.............................................................
        ELSE
C............................................COMPUTES DISTANCE
          DX=DSQRT((X1(1)-X2(1))*(X1(1)-X2(1))
     .      +(X1(2)-X2(2))*(X1(2)-X2(2))
     .      +(X1(3)-X2(3))*(X1(3)-X2(3))
     .      +(X1(4)-X2(4))*(X1(4)-X2(4)))
C.............................................................
          IRESCA=0
C.............................................................
          IF (DX.GT.1D-2) THEN   !..RESCALES THE DISTANCE
C.............................................................
            EXPAC=EXPAC+LOG(DX/EPSLY)
            X2(1)=X1(1)+(X2(1)-X1(1))/DX*EPSLY
            X2(2)=X1(2)+(X2(2)-X1(2))/DX*EPSLY
            X2(3)=X1(3)+(X2(3)-X1(3))/DX*EPSLY
            X2(4)=X1(4)+(X2(4)-X1(4))/DX*EPSLY
C.............................................................
            DX=EPSLY
            IRESCA=1       
C.............................................................
          END IF
C.............................................................
          BLYAP(ITURN)=1/DFLOAT(ITURN-1)*(LOG(DX/EPSLY)+EXPAC)
C.............................................................
        ENDIF
C.............................................................
      ENDDO
C.............................................................
      NTURN=NTURN+1
C.............................................................
 20   IF (IOUT.EQ.1) THEN
C.............................................................
        DO IDAT=1,NTURN-1,ISTEP
C.............................................................
          WRITE(31,*) IDAT,BLYAP(IDAT)
C.............................................................
        ENDDO
C.............................................................
      ENDIF
C.............................................................
      RETURN
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
      PARAMETER(MAXITER=100000)
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      COMPLEX*8 ZSING(MAXITER)
      DIMENSION X(MAXITER),XP(MAXITER)
      DIMENSION Z(MAXITER)
C.............................................................
      IF (MAXN.GT.MAXITER) THEN
        WRITE(6,*) '***ERROR(TUNENEWT): TOO MANY ITERATIONS'
        STOP
      ENDIF
C.............................................................
C    ESTIMATION OF TUNE WITH FFT 
C.............................................................
      DUEPI=DATAN(1D0)*8D0
      MFT=INT(LOG(FLOAT(MAXN))/LOG(2D0)) 
      NPOINT=2**MFT
      MAXN2=MAXN/2
      STEP=DUEPI/MAXN
C.............................................................
      SUM=0D0
      DO MF=1,MAXN
        Z(MF)=DCMPLX(X(MF),XP(MF))*(1D0+DCOS(STEP*(MF-MAXN2)))
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
      TUNEFOU=DFLOAT(NFTMAX-1)/DFLOAT(NPOINT)
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
      PARAMETER(MAXITER=100000)
      IMPLICIT DOUBLE PRECISION (A-H,O-U)
      IMPLICIT COMPLEX*16 (V-Z)
      DIMENSION Z(*),ZD(MAXITER),TUNETEST(10),TUNEVAL(10)
C............................................................  
      DUEPI=DATAN(1D0)*8D0
      ERR=1D-10
      ZU=DCMPLX(0D0,1D0)
C............................................................
C.... we divide DELTAT in 5 parts
C............................................................
      DELTAT=DELTAT/5.D0
C............................................................  
      DO ND=1,MAXN
        ZD(ND)=ZU*ND*Z(ND)
      ENDDO
C............................................................  
      ZTUNE1=CDEXP(-ZU*DUEPI*TUNEA1)
      CALL CALC(ZTUNE1,ZF,Z,MAXN)
      CALL CALC(ZTUNE1,ZFD,ZD,MAXN)
      DTUNEA1=DREAL(ZF)*DREAL(ZFD)+DIMAG(ZF)*DIMAG(ZFD)
      NUM=1
      DO NTEST=1, 10
        TUNEA2=TUNEA1+DELTAT
        ZTUNE2=CDEXP(-ZU*DUEPI*TUNEA2)
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
              ZTUNE3=CDEXP(-ZU*DUEPI*TUNE3)
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
              IF (DABS(TUNE2-TUNE1).LE.ERR) GOTO 100
           ENDDO
100        TUNETEST(NUM)=TUNE3
           TUNEVAL(NUM)=CDABS(ZF)
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
      IMPLICIT DOUBLE PRECISION (A-H,O-T)
      IMPLICIT COMPLEX*16 (U-Z)
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
CDECK  ID>, RULESORT.
C=========================================================================
C RULESORT : SORTING 48 DIPOLI 8 CELLE 2 GRUPPI 2PI ADVANCE PHASE FOR GROUP
C            SECONDO REGOLE DI WALTER MODIFIED FOR LHC II VERSION
C            V VETTORE DISTRIBUZIONE RANDOM  ERROR MEDIA = 0
C            POPOLAZIONE BIN SIMMETRICA RISPETTO MEDIA
C            VIENE FORNITA POPOLAZIONE E SEQUENZA PRELIEVO PER I  PRIMI 12
C            DIPOLI INDICE DISPARI DEL PRIMO GRUPPO.
C RULE 1  :  12 DIPOLI INDICE PARI PRIMO GRUPPO HA SEQUENZA SIMMETRICA
C            DI PRELIEVO NEI BIN ,BIN 1 = LOW QUALITY FIELD ERROR (<0)
C            SEQUENZA BIN PRELIEVO DIPOLI PARI = NBIN-SEQ+1
C RULE 2  :  II GRUPPO : SEQ GRUPPO I INVERTITA NEL PRELIEVO DAI BIN
C
C SUBROUTINE RULESORT
C
C AUTHOR : R.GRASSI 5/6/94-9/6/94 TESTED & OK FROM WALTER
C
C NOTE : OCCHIO AL PASSAGGIO PUNTATORE ARRAY IN INPUT
C        POPOLAZIONE 48 DIPOLI DERIVA DA QUELLA A 24 ==> PUO' ESSERE MIGLIORATA
C        NSEQ1 = SEQUENZA SIMILE ALLA G IN TABLE 2 SPS/AMS/Note/88-12 L.WANG
C        NGRUP1 = SEQUENZA PER 48 DIPOLI DERIVA DA SEQ1 +RULE1 + RULE2.
C

      SUBROUTINE RULESORT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NELE=80000,NLUMP=NELE/10)
      DIMENSION DUMMY(48),NSEQ1(12),NGRUP1(48)
      PARAMETER (NBIN=10)    ! N. BIN CAMPIONATURA V-ARRAY
      INTEGER NFORB(NBIN)    ! POPOLAZIONE BIN
      INTEGER INBIN(NBIN)    ! INDICE INIZIO BIN
      COMMON/RANERR/ERR(10,NELE,0:1),IND(NELE,10),IEL
C...................................................
C     SORTING SU 48 DIPOLI
C...................................................
      DATA NFORB/2,4,4,6,8,8,6,4,4,2/      ! POPOLAZIONE BIN SU 48 DIPOLI
      DATA NSEQ1/6,6,7,8,9,10,9,8,7,6,7,6/ ! SEQUENZA PRELIEVO 12 DIPOLI EVEN
C..........................................
C     COSTRUZIONE SEQUENZA NGRUP1
C..........................................
      DO I=1,12                            ! SEQUENZA NGRUP1 ( 48 D ) DA NSEQ1
         NGRUP1(2*I-1)=NSEQ1(I)
         NGRUP1(2*I)=NBIN-NSEQ1(I)+1
         NGRUP1(24+2*I-1)=NBIN-NSEQ1(I)+1  ! II GRUPPO SEGNI INVERTITI IN ERR
         NGRUP1(24+2*I)=NSEQ1(I)
      ENDDO
      INBIN(1)=1                           ! INIZIALIZZO INDICE INIZIO BIN
      DO I=2,NBIN
         INBIN(I)=INBIN(I-1)+NFORB(I-1)
      ENDDO
      DO I=1,48                            ! COPIA ERR IN DUMMY
         DUMMY(I)=ERR(2,I,IND(IEL+1,2))
      ENDDO
      CALL VSORTD(DUMMY,1,48) ! SORTING 48 DIPOLES FOR INCREASING VALUES
      DO I=1,NBIN             ! SORTING BINS ACCORDINGLY TO MAGNITUDE OF THE ERR
         CALL VABSORT(DUMMY,INBIN(I),INBIN(I)-1+NFORB(I))
      ENDDO
C..........................................
C     ASSEGNAZIONE SEQUENZA NGRUP1
C..........................................
      DO I=1,48 ! ASSEGNAZIONE DIPOLI GRUPPO I E II
         ERR(2,I,IND(IEL+1,2))=DUMMY(INBIN(NGRUP1(I)))
         INBIN(NGRUP1(I))=INBIN(NGRUP1(I))+1   ! PUNTA AL SUCCESSIVO EL. NEL BIN
      ENDDO
      RETURN
      END
CDECK  ID>, VSORTD.
C=======================================================================
C VSORTD : ORDINA PER INCREASING VALUES GLI ELEMENTI DELL'ARRAY
C          NEL RANGE SPECIFICATO DA I1 E 12
C SUBROUTINE VSORTD(V,I1,I2)
C
C V     (DOUBLE) ARRAY DIMENSIONE MAX = 10000
C I1,I2 (INTEGER) RANGE PER IL SORTING 1 <= I1 < I2
C
      SUBROUTINE VSORTD(V,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=10000)
      DIMENSION V(N)
      DO I=I1,I2-1
         DO II=I+1,I2
            IF (V(I).LE.V(II)) GOTO 10
            A=V(I)
            V(I)=V(II)
            V(II)=A
 10      ENDDO
      ENDDO
      RETURN
      END
CDECK  ID>, VABSORT.
C=======================================================================
C VABSORT : ORDINA PER ABSOLUTE VALUES GLI ELEMENTI DELL'ARRAY
C          NEL RANGE SPECIFICATO DA I1 E 12
C SUBROUTINE VABSORT(V,I1,I2)
C
C V     (DOUBLE) ARRAY DIMENSIONE MAX = 10000
C I1,I2 (INTEGER) RANGE PER IL SORTING 1 <= I1 < I2
C
      SUBROUTINE VABSORT(V,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=10000)
      DIMENSION V(N)
      DO I=I1,I2-1
         DO II=I+1,I2
            IF (DABS(V(I)).LE.DABS(V(II))) GOTO 10
            A=V(I)
            V(I)=V(II)
            V(II)=A
 10      ENDDO
      ENDDO
      RETURN
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
      PARAMETER(MAXITER=100000)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),XP(*),TUNE(MAXITER),U(MAXITER) 
C............................................................
      IF (MAX.LE.0) THEN
        WRITE(6,*) '***ERROR(TUNEFIT): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C............................................................
      DUEPI=8*DATAN(1D0)
C.............................................................
      C=DFLOAT(MAX)/2D0-DFLOAT(INT(MAX/2D0))
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
        CTHETA=(X(N)*X(N+1)+XP(N)*XP(N+1))/R1/R2
        STHETA=(-X(N)*XP(N+1)+X(N+1)*XP(N))/R1/R2
        IF (STHETA.GE.0) THEN
          THETA=DACOS(CTHETA)
          ELSE
          THETA=DUEPI-DACOS(CTHETA) 
        END IF
        PA=PA+THETA
        TUNEAPA=PA/DUEPI
        SUMAPM=SUMAPM+TUNEAPA
        IF (N.GE.MAX2) THEN
          TUNE(N-MAX2+1)=SUMAPM/DFLOAT(N)/DFLOAT(N+1)*2.D0
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),Y(*)
C............................................................
      IF (N.LE.0) THEN
        WRITE(6,*) '***ERROR(TUNEAPA): THIRD PARAMETER OUT OF BOUNDS'
        STOP
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
CDECK  ID>, DINTEG.
C=============================================================
C PERFORMES THE INTEGRATION ON THETA1 AND THETA2: THIS IS 
C SUBSTITUTED BY THE INTEGRATION OVER THE ORBIT. 
C
C     AUTHOR : E. TODESCO - BOLOGNA UNIVERSITY
C

      DOUBLE PRECISION FUNCTION DINTEG(X,PX,Y,PY,NN)
C.............................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(*),PX(*),Y(*),PY(*),R(500,500)
      INTEGER ICONT(500,500)
      PI=DATAN(1.D0)*4
C.............................................................
      N=NN
      NSTEP=SQRT(FLOAT(N))    !  FIRST GUESS FOR THE NUMBER OF INTERVALS
 10   CONTINUE
      DO J=1,NSTEP
        DO JJ=1,NSTEP
          ICONT(J,JJ)=0       !  INIZIALITATION
          R(J,JJ)=0
        END DO
      END DO
C.............................................................
      RAV=0
      DO I=1,N
        THETA1=DATAN2(PX(I),X(I))+PI !  EVALUATES THE ANGLES OF THE ITERATES
        THETA2=DATAN2(PY(I),Y(I))+PI
        INTER1=THETA1/2/PI*NSTEP+1
        INTER2=THETA2/2/PI*NSTEP+1
        ICONT(INTER1,INTER2)=ICONT(INTER1,INTER2)+1
        R(INTER1,INTER2)=R(INTER1,INTER2)+
     .              SQRT(X(I)*X(I)+PX(I)*PX(I)+Y(I)*Y(I)+PY(I)*PY(I))
      END DO
      DO J1=1,NSTEP
        DO J2=1,NSTEP
          IF(ICONT(J1,J2).EQ.0) THEN
            NSTEP=NSTEP/1.1
            GOTO 10
          END IF
          RAV=RAV+R(J1,J2)/ICONT(J1,J2)
        END DO
      END DO
      DINTEG=RAV/NSTEP/NSTEP
      N=NSTEP
      END
CDECK  ID>, MCOPYV.
C============================================================
C COPIES A SUBSET OF AN ARRAY INTO ANOTHER ONE:
C B(K)=A(I+K)     1< K < J-I
C
C AUTHOR: M. GIOVANNOZZI - CERN
C
       SUBROUTINE MCOPYV(A,I,J,B)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION A(*),B(*)
C............................................................
       IF (I.LE.0) THEN
         WRITE(6,*) '***ERROR(MCOPYV): SECOND PARAMETER OUT OF BOUNDS'
         STOP
       ENDIF
       IF (I.GT.J) THEN
         WRITE(6,*) '***ERROR(MCOPYV): THIRD PARAMETER OUT OF BOUNDS'
         STOP
       ENDIF
C............................................................
       DO K=I,J
         B(K-I+1)=A(K)
       ENDDO
C............................................................
       RETURN
C............................................................
       END
CDECK  ID>, MLINOPT.
C==================================================================
C COMPUTES THE LINEAR FUNCTIONS OF AN ACCELERATOR STRUCTURE. IT
C IS BASED ON THE EDWARD-TENG FORMALISM, SO THAT THE PRESENCE OF
C LINEAR COUPLING DOES NOT AFFECT THE CALCULATIONS.
C IBLOC IS A FLAG TO SELECT WHERE TO COMPUTE THE LINEAR FUNCTIONS:
C IF IBLOC = 0 THE LINEAR FUNCTIONS ARE COMPUTED AT EVERY ELEMENT
C            LOCATION.
C IF IBLOC = 1 THE LINEAR FUNCTIONS ARE COMPUTED AT EVERY LUMPED
C            BLOC.
C IDUMP IS A FLAG USED TO WRITE THE FUNCTIONS ON A FILE.
C IF IDUMP = 0 THEN PRINTOUT IS NOT PERFORMED.
C IF IDUMP = 1 THEN THE LINEAR FUNCTIONS ARE WRITTEN ON UNIT NUMBER
C IUNIT.
C INDEPENDENTLY FROM THE VALUE OF THE VARIABLE IDUMP, THE OUTCOME
C OF THE CALCULATIONS IS WRITTEN IN THE ARRAY OPTICS OF THE COMMON
C BLOC LINFUN. THE ARRAY OPTICS HAS THE FOLLOWING STRUCTURE:
C OPTICS(7,NELE), WHERE:
C OPTICS(1,J) = BETA_H  AT J-TH LOCATION
C OPTICS(2,J) = ALPHA_H AT J-TH LOCATION
C OPTICS(3,J) = MU_H    AT J-TH LOCATION
C OPTICS(4,J) = BETA_V  AT J-TH LOCATION
C OPTICS(5,J) = ALPHA_V AT J-TH LOCATION
C OPTICS(6,J) = MU_V    AT J-TH LOCATION
C OPTICS(7,J) = PHI     AT J-TH LOCATION
C PHI IS THE COUPLING WIDTH IN DEGREES.
C
C AUTHOR: M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
  
      SUBROUTINE MLINOPT(IBLOC,IDUMP,IUNIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      CHARACTER*10 NAME,ALBL
      DIMENSION R(4,4),B(4,4),BINV(4,4),TEM1(4,4),TEM2(4,4),PARA(20)
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
      COMMON/LINOPT/OPTICS(7,NELE),DIST(NELE),NAME(NELE),IFILL,NREC
C.............................................................
      IF (IBLOC.NE.0.AND.IBLOC.NE.1) THEN
        WRITE(6,*) '***ERROR(MLINOPT): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IDUMP.NE.0.AND.IDUMP.NE.1) THEN
        WRITE(6,*) '***ERROR(MLINOPT): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IUNIT.LE.0) THEN
        WRITE(6,*) '***ERROR(MLINOPT): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IDIME=4
      IFILL=19
      TWOPI=8*DATAN(1D0)
C.....INIZIALIZATION OF THE ONE TURN MAP PARAMETERS
      CALL VLATCAL(R,1,IMAX,1)  !..COMPUTES TRANSFER MAP
      CALL EDTENG(R,PARA)         !..COMPUTES OPTICAL PARAMETERS
      AMUH_F=PARA(1)              !..FRACTIONAL TUNE X
      AMUV_F=PARA(4)              !..FRACTIONAL TUNE Y
C.............................................................
      IF (IBLOC.EQ.0) THEN
C...............................................INITIALIZATION
        NREC=IMAX+1
C.............................................................
        DIST(1)=0D0
        NAME(1)=ALBL(1)
        OPTICS(1,1)=PARA(2)
        OPTICS(2,1)=PARA(3)
        OPTICS(3,1)=0D0
        OPTICS(4,1)=PARA(5)
        OPTICS(5,1)=PARA(6)
        OPTICS(6,1)=0D0
        OPTICS(7,1)=PARA(7)
C.............................................................
        DO IEL=2,IMAX
C.......OPTICAL FUNCTIONS FROM THE ONE TURN MAP 
C.......CALCULATED AT DIFFERENT SECTIONS
          CALL VLATCAL(B,1,IEL,1)
          CALL VISYMP(BINV,B)
          CALL MMUL(TEM1,R,BINV,IDIME)
          CALL MMUL(TEM2,B,TEM1,IDIME)
          CALL EDTENG(TEM2,PARA)
          DIST(IEL)=A(IEL,5)
          NAME(IEL)=ALBL(IEL)
          OPTICS(1,IEL)=PARA(2)
          OPTICS(2,IEL)=PARA(3)
          OPTICS(4,IEL)=PARA(5)
          OPTICS(5,IEL)=PARA(6)
          OPTICS(7,IEL)=PARA(7)
C.......CLEANING
          CALL ZEROM(TEM1)
          CALL ZEROM(TEM2)
          CALL ZEROM(B)
          CALL ZEROM(BINV)
C.......PHASE ADVANCE FROM THE TRANSFER MAP 
C.......COMPUTED BETWEEN TWO SUCCESSIVE ELEMENTS
          CALL VLATCAL(B,IEL,IEL,1)
C.......HORIZONTAL
          BX0=OPTICS(1,IEL-1)
          BX=OPTICS(1,IEL)
          AX0=OPTICS(2,IEL-1)
          AX=OPTICS(2,IEL)
          CX1=DSQRT(BX0*BX)
          CX2=BX0+BX
          CX3=BX*AX0-BX0*AX
          AUX1=B(1,2)/CX1
          AUX2=(B(1,1)+B(2,2)-AUX1*CX3/CX1)/CX2*CX1
          PHASE1=DATAN2(AUX1,AUX2)
          AMUH=PHASE1/TWOPI
          OPTICS(3,IEL)=OPTICS(3,IEL-1)+AMUH
C.......VERTICAL
          BY0=OPTICS(4,IEL-1)
          BY=OPTICS(4,IEL)
          AY0=OPTICS(5,IEL-1)
          AY=OPTICS(5,IEL)
          CY1=DSQRT(BY0*BY)
          CY2=BY0+BY
          CY3=BY*AY0-BY0*AY
          AUY1=B(3,4)/CY1
          AUY2=(B(3,3)+B(4,4)-AUY1*CY3/CY1)/CY2*CY1
          PHASE2=DATAN2(AUY1,AUY2)
          AMUV=PHASE2/TWOPI
          OPTICS(6,IEL)=OPTICS(6,IEL-1)+AMUV
C.......CLEANING
          CALL ZEROM(TEM1)
          CALL ZEROM(TEM2)
          CALL ZEROM(B)
          CALL ZEROM(BINV)
C.............................................................
100     ENDDO
C.............................................................
        DIST(IMAX+1)=A(1,5)
        NAME(IMAX+1)=ALBL(1)
        OPTICS(1,IMAX+1)=OPTICS(1,1)
        OPTICS(2,IMAX+1)=OPTICS(2,1)
        OPTICS(3,IMAX+1)=AMUH_F
        OPTICS(4,IMAX+1)=OPTICS(4,1)
        OPTICS(5,IMAX+1)=OPTICS(5,1)
        OPTICS(6,IMAX+1)=AMUV_F
        OPTICS(7,IMAX+1)=OPTICS(7,1)
C.............................................................
      ELSE
C...............................................INITIALIZATION
        NREC=NLMAX+2
C.............................................................
        DIST(1)=0D+0
        NAME(1)='START'
C.............................................................
        OPTICS(1,1)=PARA(2)
        OPTICS(2,1)=PARA(3)
        OPTICS(3,1)=0D0
        OPTICS(4,1)=PARA(5)
        OPTICS(5,1)=PARA(6)
        OPTICS(6,1)=0D0
        OPTICS(7,1)=PARA(7)
C.............................................................
        NLMAX2=0
        LIMIT=2000
        IKICK=107
        CALL LUMPIT(NLMAX2,LIMIT,IKICK)   !..LUMPS THE LATTICE
C.............................................................
C N.B.: INL CONTAINS THE STARTING AND ENDING POSITIONS OF THE 
C       NON-LINEAR BLOCKS AND RL IS THE MATRIX OF A GIVEN LUMP
C       THE NON-LINEAR BLOCKS MUST BE OF ZERO LENGTH.
C.............................................................
c        DIST(2)=0D+0
c        DO J=1,INL(1,1)
c          DIST(2)=DIST(2)+A(J,1)
c        ENDDO
c        DO ILUM=2,NLMAX !..COMPUTES THE LENGTH OF THE BLOCKS
c          DIST(ILUM+1)=DIST(ILUM)
c          DO JLUM=INL(1,ILUM-1),INL(1,ILUM)
c            DIST(ILUM+1)=DIST(ILUM+1)+A(JLUM,1)
c          ENDDO
c        ENDDO
        DO ILUM=1,NLMAX !..COMPUTES THE LENGTH OF THE BLOCKS
          DIST(ILUM+1)=A(INL(1,ILUM),5)
        ENDDO
C.............................................................
        DO IEL=2,NLMAX+1
C.......OPTICAL FUNCTIONS FROM THE ONE TURN MAP
C.......CALCULATED AT DIFFERENT SECTIONS
          DO IROW=1,IDIME
            DO ICOL=1,IDIME
              B(IROW,ICOL)=RL(IROW,ICOL,IEL-1)
            ENDDO
          ENDDO
          CALL VISYMP(BINV,B)
          CALL MMUL(TEM1,R,BINV,IDIME)
          CALL MMUL(TEM2,B,TEM1,IDIME)
          CALL EDTENG(TEM2,PARA)
C.............................................................
          OPTICS(1,IEL)=PARA(2)
          OPTICS(2,IEL)=PARA(3)
          OPTICS(4,IEL)=PARA(5)
          OPTICS(5,IEL)=PARA(6)
          OPTICS(7,IEL)=PARA(7)
C.......UPDATING R TO THE NEXT BLOCK AND CLEANING
          DO IROW=1,IDIME
            DO ICOL=1,IDIME
              R(IROW,ICOL)=TEM2(IROW,ICOL)
            ENDDO
          ENDDO
          CALL ZEROM(TEM1)
          CALL ZEROM(TEM2)
C.......RATIONAL PART OF THE PHASE ADVANCE FROM THE TRANSFER MAP 
C.......COMPUTED BETWEEN TWO SUCCESSIVE BLOCKS
C.......HORIZONTAL
          BX0=OPTICS(1,IEL-1)
          BX=OPTICS(1,IEL)
          AX0=OPTICS(2,IEL-1)
          AX=OPTICS(2,IEL)
          CX1=DSQRT(BX0*BX)
          CX2=BX0+BX
          CX3=BX*AX0-BX0*AX
          AUX1=B(1,2)/CX1
          AUX2=(B(1,1)+B(2,2)-AUX1*CX3/CX1)/CX2*CX1
          PHASE1=DATAN2(AUX1,AUX2)
          AMUH=PHASE1/TWOPI
          OPTICS(3,IEL)=OPTICS(3,IEL-1)+AMUH
C.......VERTICAL
          BY0=OPTICS(4,IEL-1)
          BY=OPTICS(4,IEL)
          AY0=OPTICS(5,IEL-1)
          AY=OPTICS(5,IEL)
          CY1=DSQRT(BY0*BY)
          CY2=BY0+BY
          CY3=BY*AY0-BY0*AY
          AUY1=B(3,4)/CY1
          AUY2=(B(3,3)+B(4,4)-AUY1*CY3/CY1)/CY2*CY1
          PHASE2=DATAN2(AUY1,AUY2)
          AMUV=PHASE2/TWOPI
          OPTICS(6,IEL)=OPTICS(6,IEL-1)+AMUV
C.............................................................
          NAME(IEL)='BLOCK'
          WRITE(NAME(IEL)(6:10),'(I4)') IEL-1
C.......CLEANING
          CALL ZEROM(B)
          CALL ZEROM(BINV)
C.............................................................
        ENDDO
C.............................................................
C        DIST(NLMAX+2)=0D+0
C        DO KF=1,IMAX
C          DIST(NLMAX+2)=DIST(NLMAX+2)+A(KF,1)
C        ENDDO
        DIST(NLMAX+2)=A(IMAX,5)
        NAME(NLMAX+2)='END'
C.............................................................
        OPTICS(1,NLMAX+2)=OPTICS(1,1)
        OPTICS(2,NLMAX+2)=OPTICS(2,1)
        OPTICS(3,NLMAX+2)=AMUH_F
        OPTICS(4,NLMAX+2)=OPTICS(4,1)
        OPTICS(5,NLMAX+2)=OPTICS(5,1)
        OPTICS(6,NLMAX+2)=AMUV_F
        OPTICS(7,NLMAX+2)=OPTICS(7,1)
C.............................................................
      ENDIF
C.............................................................
      IF (IDUMP.EQ.1) CALL MWRIOPT(IUNIT) !..WRITES OUT RESULTS
C,............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MERR.
C=================================================================
C GENERATES A SEQUENCE OF RANDOM ERRORS TO BE ASSIGNED TO
C SPECIFIED ELEMENTS. THE ERRORS CAN BE ASSIGNED TO ZERO LENGTH
C ELEMENTS ONLY, STARTING FROM THE QUADRUPOLAR COMPONENT.
C MULTMI IS THE MINIMUM ORDER OF THE MULTIPOLAR COMPONENT TO BE
C GENERATED.
C MULTMA IS THE MAXIMUM ORDER OF THE MULTIPOLAR COMPONENT TO BE
C GENERATED.
C RADIUS IS THE RADIUS USED FOR THE ERROR MEASUREMENT. IT IS
C EXPRESSED IN mm.
C PHI IS THE BENDING ANGLE OF THE DIPOLE USED AS A REFERENCE FOR
C THE ERROR EXPRESSED IN mrad.
C BMEAN IS AN ARRAY CONTAINING THE SYSTEMATIC PART OF THE NORMAL
C ERRORS.
C BSIG IS AN ARRAY CONTAINING THE RANDOM PART OF THE NORMAL
C ERRORS.
C AMEAN IS AN ARRAY CONTAINING THE SYSTEMATIC PART OF THE SKEW
C ERRORS.
C ASIG IS AN ARRAY CONTAINING THE RANDOM PART OF THE SKEW ERRORS.
C ISEED IS A SEED USED FOR THE RANDOM NUMBER GENERATOR.
C IASSIGN IS A FLAG TO SPECIFY WHETHER THE GENERATED ERRORS HAVE
C TO BE ASSIGNED TO A PHYSICAL ELEMENT OR TO A COMMON BLOCK FOR
C SUBSEQUENT ANALYSIS (REORDERING).
C IASSIGN = 0 THE ARRAY ERR IN THE COMMON BLOCK RANERR IS FILLED.
C IASSIGN = 1 THE ERRORS ARE ASSIGNED TO THE PROPER MULTIPOLAR
C ELEMENT.
C ELEMENT IS A STRING CONTAINING THE NAME OF THE PHYSICAL
C MAGNET WHO IS SUPPOSED TO BE AFFECTED BY THE GENERATED ERRORS.
C
C AUTHOR: M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MERR(MULTMI,MULTMA,RADIUS,PHIDIP,BMEAN,BSIG,
     .                AMEAN,ASIG,ISEED,IASSIGN,ELEMENT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      CHARACTER*10 ELEMENT
      DIMENSION BMEAN(10),BSIG(10),AMEAN(10),ASIG(10)
      COMMON/RANERR/ERR(10,NELE,0:1),IND(NELE,10),IEL
C.............................................................
      IF (MULTMI.LT.0.OR.MULTMI.GT.MULTMA) THEN
        WRITE(6,*) '***ERROR(MERR): FIRST PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (MULTMA.LT.0) THEN
        WRITE(6,*) '***ERROR(MERR): SECOND PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (RADIUS.LT.0D0) THEN
        WRITE(6,*) '***ERROR(MERR): SECOND PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (PHIDIP.LT.0D0) THEN
        WRITE(6,*) '***ERROR(MERR): THIRD PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (ISEED.LT.0) THEN
        WRITE(6,*) '***ERROR(MERR): NINETH PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (IASSIGN.NE.0.AND.IASSIGN.NE.1) THEN
        WRITE(6,*) '***ERROR(MERR): TENTH PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
C.............................................................
      KSEED=ISEED-1!..INITIALIZATION FOR RANDOM NUMBER GENERATION
C.............................................................
      DO MULT=MULTMI,MULTMA
C.............................................................
        KSEED=KSEED+1
        CALL MRANERR(MULT,RADIUS,PHIDIP,BMEAN(MULT),BSIG(MULT),
     .               0,KSEED,IASSIGN,ELEMENT)
C.............................................................
        KSEED=KSEED+1
        CALL MRANERR(MULT,RADIUS,PHIDIP,AMEAN(MULT),ASIG(MULT),
     .               1,KSEED,IASSIGN,ELEMENT)
C.............................................................
      ENDDO
C.............................................................
      RETURN
C.............................................................
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(*),P(*)
C............................................................
      COMMON/TUNEPAR/ADVSIG,ADVMIN,ADVMAX
C............................................................
      IF (N.LE.0) THEN
        WRITE(6,*) '***ERROR(TUNEAPA): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C............................................................
      PI=4*DATAN(1.D0)
      ADV=0.D0
      ADV1=10D0
      ADV2=-10D0
      ADVS=0D0
C.....................EVALUATION OF THE AVERAGE PHASE ADVANCE
      DO I=1,N-1
        S1=X(I+1)*X(I)+P(I+1)*P(I)
        S2=SQRT((X(I)*X(I)+P(I)*P(I))*(X(I+1)*X(I+1)+P(I+1)*P(I+1)))
        S3=DACOS(S1/S2)
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
      ADVSIG=DSQRT(ADVS/(N-1)/4/PI/PI-
     .            (ADV/(N-1)/2/PI)*(ADV/(N-1)/2/PI))
C......................................NORMALIZATION TO [0,1]
      ADV1=-ADV1
      IF(ADV1.LT.0) ADV1=1+ADV1
      ADV2=-ADV2
      IF(ADV2.LT.0) ADV2=1+ADV2
C...........................................FINDS MIN AND MAX
      ADVMIN=MIN(ADV1,ADV2)/2/PI
      ADVMAX=MAX(ADV1,ADV2)/2/PI
C......................................NORMALIZATION TO [0,1]
      ADV=-ADV
      IF(ADV.LT.0) THEN
        TUNEAPA=1+ADV/(N-1)/2/PI
      ELSE
        TUNEAPA=ADV/(N-1)/2/PI
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
      PARAMETER(MAXITER=100000)
      DOUBLE PRECISION X(*),P(*)
      COMPLEX  Z(MAXITER)
C..................................................CHECK OF N
      IF(N.GT.MAXITER) THEN
        WRITE(6,*) '***ERROR(TUNEFFT): TOO MANY ITERATES'
         STOP
      ENDIF
C............................................................
      IF (N.LE.0) THEN
        WRITE(6,*) '***ERROR(TUNEFFT): THIRD PARAMETER OUT OF BOUNDS'
        STOP
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
      PARAMETER(MAXITER=100000)
      DOUBLE PRECISION X(*),P(*)
      COMPLEX  Z(MAXITER)
C..................................................CHECK OF N
      IF(N.GT.MAXITER) THEN
        WRITE(6,*) '***ERROR(TUNEFFTI): TOO MANY ITERATES'
         STOP
      ENDIF
C............................................................
      IF (N.LE.0) THEN
        WRITE(6,*) '***ERROR(TUNEFFTI): THIRD PARAMETER OUT OF BOUNDS'
        STOP
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
      Y12=LOG(Y1/Y2)
      Y13=LOG(Y1/Y3)
      X212=X1*X1-X2*X2
      X213=X1*X1-X3*X3
C........COMPUTATION OF THE POSITION OF THE INTERPOLATED PEAK
      A=X212*Y13-X213*Y12
      B=X12*Y13-X13*Y12
      TUNEFFTI=(A/2/B-1)/NPOINT
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
      PARAMETER(MAXITER=100000)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION X(*),PX(*)
      COMPLEX*8 ZSING(MAXITER)
      COMPLEX*16 Z(MAXITER),FOME,ZC,SD,SP
      DUEPI=8*DATAN(1D+0)
C...............................CHECK OF THE ITERATION NUMBER
      IF(MAX.GT.MAXITER) THEN
        WRITE(6,*) '***ERROR(TUNELASK): TOO MANY ITERATIONS'
        STOP
      ENDIF
C............................................................
      IF (MAX.LE.0) THEN
        WRITE(6,*) '***ERROR(TUNELASK): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.................................ESTIMATION OF TUNE WITH FFT
      SUM=0D0
      MFT=INT(LOG(FLOAT(MAX))/LOG(2D+0))
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
      OMEMIN=TUNEFOU*DUEPI-FSIGMA*DUEPI/NPOINT
      OMEMAX=TUNEFOU*DUEPI+FSIGMA*DUEPI/NPOINT
C.....JITER=8 PROVIDES A PRECISION WHICH IS 1.5E-4 * 1/2**MFT
      JITER=8
C................JMAX=7 IS THE VALUE WHICH MINIMIZES CPU TIME
      JMAX=7
      MAX1=MAX-1
      STEP=DUEPI*.5D0/MAX1
C.........................................BISECTION PROCEDURE
      DO JIT=1,JITER
        FOMEGA=0D+0
        DO J=1,JMAX
          OME=OMEMIN+(OMEMAX-OMEMIN)/(JMAX-1D+0)*(J-1D+0)
          DO N=1,MAX
            ZC=(X(N)-(0D+0,1D+0)*PX(N))
     .        *(1D0+DCOS(STEP*(2*N-MAX1)))
            Z(N)=ZC*CDEXP(-(0D+0,1D+0)*OME*N)
          ENDDO
C..COMPUTATION OF SCALAR PRODUCT WITH ITERATED BODE ALGORITHM
          FOME=(0D0,0D0)
          MBODE=(MAX-5)/4
          DO I=0,MBODE
            K=4*I
            FOME=FOME+(7D0*(Z(1+K)+Z(5+K))+32D0*(Z(2+K)+Z(4+K))
     .               +12D0*Z(3+K))
          ENDDO
          FOME=.5D0*FOME/45D0/DFLOAT(MBODE+1)
C..........SEARCH FOR MAXIMUM OF SCALAR PRODUCT AND DEFINITION
C..........OF THE NEW INTERVAL (OMEMIN,OMEMAX) WHERE TO
C........................................RESTART THE PROCEDURE
          ABSFOM=CDABS(FOME)
          IF (ABSFOM.GT.FOMEGA) THEN
            FOMEGA=ABSFOM
            TUNELASK=-OME/DUEPI
            JOM=J
          ENDIF
        ENDDO
        OMEMIN=OMEMIN+(OMEMAX-OMEMIN)/(JMAX-1D+0)*(JOM-2D+0)
        OMEMAX=OMEMIN+(OMEMAX-OMEMIN)/(JMAX-1D+0)*JOM
      ENDDO
C......................................NORMALIZATION TO [0,1]
      TUNELASK=-TUNELASK
      IF(TUNELASK.LE.0) THEN
        TUNELASK=1+TUNELASK
      ENDIF
C............................................................
      END
CDECK  ID>, TUNEDIAG.
C=============================================================
C CREATES TWO FILES IN ORDER TO DRAW THE RESONANT LINES IN THE
C TUNE DIAGRAM. ALL THE RESONANT LINES FROM ORDER IORDMIN TO
C ORDER IORDMAX AND DRAWN. XMIN, XMAX, YMIN AND YMAX ARE THE
C COORDINATES OF THE SQUARE IN THE TUNE WHERE ONE WANTS TO
C PLOT THE DIAGRAM. TWO FILES ARE PRODUCED: IN UNIT 20 THE
C LIST OF THE ORDERED POINTS WHICH PRODUCE THE DIAGRAM IS
C GIVEN IN UNIT 21 A PAW EXECUTABLE FILE IS GENERATED
C
C AUTHOR:    E. TODESCO - INFN OF BOLOGNA
C
 
      SUBROUTINE TUNEDIAG(IORDMIN,IORDMAX,XMIN,XMAX,YMIN,YMAX)
C.............................................................
      WRITE(21,*) 'FORTRAN/FILE 33 TUDIA.PS'
      WRITE(21,*) 'IGSET LTYP 2'
      WRITE(21,*) 'SET CHHE 0.22'
      WRITE(21,*) 'GR/ME 33 -111'
C.............................................................
C     PLEASE PUT -111 IF YOU WANT PS
C     PLEASE PUT -113 IF YOU WANT EPS
C.............................................................
      WRITE(21,*) 'V/RE A,B FOR020.DAT'
      WRITE(21,97) XMIN,XMAX,YMIN,YMAX
      EPS=(XMAX-XMIN)/80
C.............................................................
 97   FORMAT(1X,'NULL ',4F10.5)
 98   FORMAT(1X,'IGSET TANG ',I5)
 99   FORMAT(1X,'ITX',2F10.5,' ''(',I1,',',I1,')'' ')
 100  FORMAT(1X,'ITX',2F10.5,' ''(',I1,',',I2,')'' ')
C.............................................................
      WRITE(20,*) XMIN,YMIN
      WRITE(20,*) XMIN,YMAX
      WRITE(20,*) XMAX,YMAX
      WRITE(20,*) XMAX,YMIN
      WRITE(20,*) XMIN,YMIN
      ICONT=5
C.............................................................
      WRITE(21,98) 90
      DO IQ=IORDMIN,IORDMAX
        DO IP=1,IQ-1
          IF(FLOAT(IP)/IQ.GT.XMIN.AND.FLOAT(IP)/IQ.LT.XMAX) THEN
            CALL SIMPL1(IP,IQ,IPP,IQQ)
            IF(IQQ.LT.IORDMIN.OR.IQQ.EQ.IQ) THEN
              WRITE(20,*) FLOAT(IP)/IQ,YMIN
              WRITE(20,*) FLOAT(IP)/IQ,YMAX
              WRITE(20,*) XMIN,YMAX
              WRITE(20,*) XMIN,YMIN
              WRITE(21,99) FLOAT(IP)/IQ-EPS,YMIN+EPS,IQ,0
              ICONT=ICONT+4
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C.............................................................
      WRITE(21,98) 0
      DO IQ=IORDMIN,IORDMAX
        DO IP=1,IQ-1
          IF(FLOAT(IP)/IQ.GT.YMIN.AND.FLOAT(IP)/IQ.LT.YMAX) THEN
            CALL SIMPL1(IP,IQ,IPP,IQQ)
            IF(IQQ.LT.IORDMIN.OR.IQQ.EQ.IQ) THEN
              WRITE(20,*) XMIN,FLOAT(IP)/IQ
              WRITE(20,*) XMAX,FLOAT(IP)/IQ
              WRITE(20,*) XMAX,YMIN
              WRITE(20,*) XMIN,YMIN
              WRITE(21,99) XMIN+EPS,FLOAT(IP)/IQ+EPS,0,IQ
              ICONT=ICONT+4
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C.............................................................
      DO I=IORDMIN,IORDMAX
        DO IP=1,I-1
          WRITE(20,*) XMIN,YMAX
          ICONT=ICONT+1
          IQ=I-IP
          DO L=1,IP+IQ-1
            CALL SIMPL2(IP,IQ,L,IPP,IQQ,LL)
            IF(IPP+IQQ.LT.IORDMIN.OR.IQQ.EQ.IQ) THEN
              IANG=-ATAN2(FLOAT(IP),FLOAT(IQ))/3.141592654*180
              X=(L-IQ*YMAX)/IP
              Y=(L-IP*XMIN)/IQ
              IF(X.GE.XMIN.AND.X.LE.XMAX) THEN
                WRITE(21,98) IANG
                WRITE(21,99) X+2*EPS*SIN(IANG/180.*3.1415),
     .                       YMAX-2*EPS*COS(IANG/180.*3.1415),IP,IQ
                WRITE(20,*) X,YMAX
              ELSEIF (Y.GE.YMIN.AND.Y.LE.YMAX) THEN
                WRITE(21,98) IANG
                WRITE(21,99) XMIN-EPS*SIN(IANG/180.*3.1415),
     .                       Y+EPS*COS(IANG/180.*3.1415),IP,IQ
                WRITE(20,*) XMIN,Y
              ENDIF
              X=(L-IQ*YMIN)/IP
              Y=(L-IP*XMAX)/IQ
              IF(X.GE.XMIN.AND.X.LE.XMAX) THEN
                WRITE(20,*) X,YMIN
                WRITE(20,*) XMIN,YMIN
                WRITE(20,*) XMIN,YMAX
                ICONT=ICONT+4
              ELSEIF (Y.GE.YMIN.AND.Y.LE.YMAX) THEN
                WRITE(20,*) XMAX,Y
                WRITE(20,*) XMAX,YMAX
                WRITE(20,*) XMIN,YMAX
                ICONT=ICONT+4
              ENDIF
            ENDIF
          ENDDO
          IQ=IP-I
          WRITE(20,*) XMIN,YMIN
          ICONT=ICONT+1
          DO L=IQ+1,IP-1
            CALL SIMPL2(IP,IQ,L,IPP,IQQ,LL)
            IF(IPP-IQQ.LT.IORDMIN.OR.IQQ.EQ.IQ) THEN
              IANG=ATAN2(FLOAT(IP),-FLOAT(IQ))/3.141592654*180
              X=(L-IQ*YMIN)/IP
              Y=(L-IP*XMIN)/IQ
              IF(X.GE.XMIN.AND.X.LE.XMAX) THEN
                WRITE(21,98) IANG
                WRITE(21,100) X-EPS*SIN(IANG/180.*3.1415),
     .                       YMIN+EPS*COS(IANG/180.*3.1415),IP,IQ
                WRITE(20,*) X,YMIN
              ELSEIF (Y.GE.YMIN.AND.Y.LE.YMAX) THEN
                WRITE(21,98) IANG
                WRITE(21,100) XMIN+2*EPS*SIN(IANG/180.*3.1415),
     .                       Y-2*EPS*COS(IANG/180.*3.1415),IP,IQ
                WRITE(20,*) XMIN,Y
              ENDIF
              X=(L-IQ*YMAX)/IP
              Y=(L-IP*XMAX)/IQ
              IF(X.GE.XMIN.AND.X.LE.XMAX) THEN
                WRITE(20,*) X,YMAX
                WRITE(20,*) XMIN,YMAX
                WRITE(20,*) XMIN,YMIN
                ICONT=ICONT+4
              ELSEIF (Y.GE.YMIN.AND.Y.LE.YMAX) THEN
                WRITE(20,*) XMAX,Y
                WRITE(20,*) XMAX,YMIN
                WRITE(20,*) XMIN,YMIN
                ICONT=ICONT+4
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C.............................................................
      WRITE(*,*) ICONT
      WRITE(21,*) 'GRAPH',ICONT,' A  B  L'
      WRITE(21,*) 'GR/ME 33 -111'
      WRITE(21,*) 'V/DE A,B'
      WRITE(21,*) 'FORTRAN/CLOSE 33'
C.............................................................
      END
CDECK  ID>, SIMPL1.
C=============================================================
C GIVEN TWO INTEGER NUMBERS IP AND IQ, WRITES IN IPP AND IQQ
C  1)  IP AND IQ, IF HAVE NO COMMON DIVISOR
C  2)  IP/I AND IQ/I, WHERE I IS THE SMALLEST COMMON DIVISOR
C      (GREATER THAN ONE)
C
C AUTHOR:    E. TODESCO - INFN OF BOLOGNA
C
 
      SUBROUTINE SIMPL1(IP,IQ,IPP,IQQ)
C.............................................................
      IPP=IP
      IQQ=IQ
      DO I=2,IP
        ISCR1=IP/I
        ISCR2=IQ/I
        IF(ISCR1*I.EQ.IP.AND.ISCR2*I.EQ.IQ) THEN
          IPP=ISCR1
          IQQ=ISCR2
          GOTO 20
        ENDIF
      ENDDO
 20   CONTINUE
C.............................................................
      END
CDECK  ID>, SIMPL2.
C=============================================================
C GIVEN THREE INTEGER NUMBERS IP,IQ AND L, WRITES IN IPP,
C IQQ AND LL
C  1)  IP, IQ AND LL, IF THEY HAVE NO COMMON DIVISOR
C  2)  IP/I, IQ/I AND LL/I, WHERE I IS THE SMALLEST COMMON
C      DIVISOR (GREATER  THAN ONE)
C
C AUTHOR:    E. TODESCO - INFN OF BOLOGNA
C
 
      SUBROUTINE SIMPL2(IP,IQ,L,IPP,IQQ,LL)
C.............................................................
      IPP=IP
      IQQ=IQ
      LL=L
      DO I=2,IP
        ISCR1=IP/I
        ISCR2=IQ/I
        ISCR3=L/I
        IF(ISCR1*I.EQ.IP.AND.ISCR2*I.EQ.IQ.AND.ISCR3*I.EQ.L) THEN
          IPP=ISCR1
          IQQ=ISCR2
          LL=ISCR3
          GOTO 20
        ENDIF
      ENDDO
 20   CONTINUE
C.............................................................
      END
CDECK  ID>, MAPDYNNORM.
C=============================================================
C IT COMPUTES THE DYNAMIC APERTURE USING A COMBINED METHOD
C BASED ON TRACKING AND NORMAL FORMS. TRACKING IS COMPUTED
C OVER THE GRID
C              X(1)= R COS(ALPHA)
C              X(2)= 0
C              X(3)= R SIN(ALPHA)
C              X(4)= 0
C WHERE 0<R<RMAX AND 0<ALPHA<PI/2. FOR EACH ALPHA, THE
C LAST STABLE R IS EVALUATED, AND NONLINEAR INVARIANTS
C ARE COMPUTED USING NORMAL FORMS. THE DYNAMIC APERTURE
C IS THEN COMPUTED AS AN INTEGRAL OVER ALPHA.
C
C AUTHOR:    E. TODESCO - INFN & CERN (M. GIOVANNOZZI - CERN - HAS
C INTRODUCED SOME CHANGES)
C
 
      SUBROUTINE MAPDYNNORM(IMAP,IORD,NTURN,INORM,IOUT,RMAX,
     .                      RSTEP,NANG,AMPLNORM,ERRNORM,INORMBEST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      PARAMETER (NMAX=50000,NPOINT=50000)
      DIMENSION X(4),RHO0(1000)
      COMPLEX*16 Z1,Z2,ZETA1,ZETA2,ZZ1,ZZ2,PSI,AUTORES
      COMMON/INVERSA/PSI(NDIM,4)
      COMMON/RISONANZA/IRES1(2),IRES2(2),AUTORES(4),IRESON,ICASE
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C.............................................................
      EXTERNAL NLFUN
C.............................................................
      IF (NTURN.GT.NMAX) THEN
        WRITE(6,*) '***ERROR(MAPDYNNORM): TOO MANY ITERATIONS'
        WRITE(6,*) '***MESSAGE(MAPDYNNORM): THE NUMBER OF ',
     .             'ITERATIONS MUST BE < ',NMAX
        STOP
      ENDIF
C.............................................................
      IF (INORM.LT.2.OR.INORM.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MAPDYNNORM): NORMAL FORM ORDER ',
     .             'OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IOUT.NE.0.AND.IOUT.NE.1) THEN
        WRITE(6,*) '***ERROR(MAPDYNNORM): IOUT VARIABLE OUT ',
     .             'OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (RMAX.LT.0) THEN
        WRITE(6,*) '***ERROR(MAPDYNNORM): THE MAXIMUM RHO IS ',
     .             'NEGATIVE'
        STOP
      ENDIF
C.............................................................
      IF (RSTEP.LT.0) THEN
        WRITE(6,*) '***ERROR(MAPDYNNORM): THE STEP IN RHO IS ',
     .             'NEGATIVE'
        STOP
      ENDIF
C.............................................................
      IF (NANG.LT.1) THEN
        WRITE(6,*) '***ERROR(MAPDYNNORM): THE NUMBER OF ',
     .   'ANGLES IS <1 '
        STOP
      ENDIF
C.............................................................
      IF (NANG.GT.1000) THEN
        WRITE(6,*) '***ERROR(MAPDYNNORM): TOO MANY ANGLES'
        WRITE(6,*) 'THE NUMBER OF ANGLES MUST BE <= 1000'
        STOP
      ENDIF
C.............................................................
      IF (IMAP.NE.1.AND.IMAP.NE.2) THEN
        WRITE(6,*) '***ERROR(MAPDYNNORM): THE LAST PARAMETER IS ',
     .             'OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IMAP.EQ.1) CALL INIZ_PUNT(NORMAX)
C.............................................................
      PI=4D0*DATAN(1D0)
      TWOPI=8D0*DATAN(1D0)
      IDIME=4
C.................................NORMAL FORM INPUT PARAMETERS
      IRESON=0
      ICASE=0
      IRES1(1)=0
      IRES1(2)=0
      IRES2(1)=0
      IRES2(2)=0
      AUTORES(1)=DCMPLX(1D0,0D0)
      AUTORES(2)=DCMPLX(1D0,0D0)
      AUTORES(3)=DCMPLX(1D0,0D0)
      AUTORES(4)=DCMPLX(1D0,0D0)
C......................................NORMAL FORM COMPUTATION
      INTERACT=0
      IINV=1
      IOUTG=3
      IAZAN=1
      CALL GENERALE_4D(RMAP,IORD,INORM,INTERACT,IINV,IOUTG,IAZAN)
C............................................TRACKING ANALYSIS
      DYNAP=0
      DO I3=1,NANG                !..LOOP ON ALPHA
        A3=FLOAT(I3)/NANG*PI/2
        CA3=COS(A3)
        SA3=SIN(A3)
        DO RHO=RSTEP,RMAX,RSTEP   !..LOOP ON RHO
          X(1)=RHO*CA3
          X(2)=0.
          X(3)=RHO*SA3
          X(4)=0.
C.............................................................
          IF (IMAP.EQ.1) THEN
            ITURN=MANY_HENON(X,NTURN) !..TRACKING HENON
          ELSE
            ITURN=MANYTURN(X,NTURN,IDIME,NLFUN) !..GENERAL TRACKING
          ENDIF
C.............................................................
          IF (ITURN.NE.NTURN) THEN
            RHO0(I3)=RHO-RSTEP
            GOTO 12
          ENDIF
        ENDDO
        WRITE(6,*) '***ERROR(MAPDYNNORM): RMAX STABLE'
 12     CONTINUE
      ENDDO
C........................................ NORMAL FORM ANALYSIS
      ERRNORM=1E5
      AMPLNORM=0
      DO II=2,INORM
        DYNAP=0
        DO I3=1,NANG
          A3=FLOAT(I3)/NANG*PI/2
          CA3=COS(A3)
          SA3=SIN(A3)
          Z1=RHO0(I3)*CA3
          Z2=RHO0(I3)*SA3
          CALL EVAL_PSI(II,Z1,Z2,ZETA1,ZETA2)
          RHONORM=SQRT(ABS(ZETA1)**2+ABS(ZETA2)**2)
          DYNAP=DYNAP+(RHONORM**4)*SIN(2*A3)
        ENDDO
        DYNAP=DYNAP*PI*PI*PI/4/NANG
        AMPL=SQRT(SQRT(2*DYNAP/PI/PI))
C.............................................................
        Z1=RHO0(NANG/2)/SQRT(2.)
        Z2=Z1
        CALL EVAL_PSI(II,Z1,Z2,ZETA1,ZETA2)
        CALL EVAL_PHI(II,ZETA1,ZETA2,ZZ1,ZZ2)
        ERR=SQRT(ABS(Z1-ZZ1)**2+ABS(Z2-ZZ2)**2)/RHO0(NANG)
C.............................................................
        IF(ERRNORM.GT.ERR) THEN
          AMPLNORM=AMPL
          ERRNORM=ERR
          INORMBEST=II
        ENDIF
C.............................................................
        IF(IOUT.EQ.1) THEN
          WRITE(6,99) II,AMPL,ERR
        ENDIF
 99     FORMAT(1X,'ORDER',I3,'     AVERAGE RADIUS',F15.5,
     .       '      ERROR',F10.3)
C.............................................................
      ENDDO
C.............................................................
      END
CDECK  ID>, EVAL_PSI.
C=============================================================
C IT COMPUTES ZETA=PSI(Z) USING TRUNCATION OF ORDER IORD.
C PSI IS THE INVERSE CONJUGATING FUNCTION OF THE NORMAL FORM.
C
C AUTHOR:    E. TODESCO - INFN & CERN
C
 
      SUBROUTINE EVAL_PSI(IORD,Z1,Z2,ZETA1,ZETA2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      COMPLEX*16 PSI,Z1,ZC1,Z2,ZC2,ZETA1,ZETA2
      COMPLEX*16 Z1POW(0:NORMAX),ZC1POW(0:NORMAX)
      COMPLEX*16 Z2POW(0:NORMAX),ZC2POW(0:NORMAX)
      COMMON/INVERSA/PSI(NDIM,4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C.............................................................
      ZC1=DCONJG(Z1)
      ZC2=DCONJG(Z2)
      Z1POW(0)=1.D0
      ZC1POW(0)=1.D0
      Z2POW(0)=1.D0
      ZC2POW(0)=1.D0
C.............................................................
      DO N=1,IORD
        Z1POW(N)=Z1POW(N-1)*Z1
        ZC1POW(N)=ZC1POW(N-1)*ZC1
        Z2POW(N)=Z2POW(N-1)*Z2
        ZC2POW(N)=ZC2POW(N-1)*ZC2
      ENDDO
C.............................................................
      ZETA1=0.D0
C.............................................................
      DO N=1,IORD
C.............................................................
        DO K1=0,N
C.............................................................
          DO K2=0,N-K1
C.............................................................
            DO K3=0,N-K1-K2
C.............................................................
              K4=N-K1-K2-K3
              ZETA1=ZETA1+PSI(IPUNT(K1,K2,K3,K4),1)
     .               *Z1POW(K1)*ZC1POW(K2)*Z2POW(K3)*ZC2POW(K4)
            ENDDO
C.............................................................
          ENDDO
C.............................................................
        ENDDO
C.............................................................
      ENDDO
C.............................................................
      ZETA2=0.D0
C.............................................................
      DO N=1,IORD
C.............................................................
        DO K1=0,N
C.............................................................
          DO K2=0,N-K1
C.............................................................
            DO K3=0,N-K1-K2
C.............................................................
              K4=N-K1-K2-K3
              ZETA2=ZETA2+PSI(IPUNT(K1,K2,K3,K4),3)
     .               *Z1POW(K1)*ZC1POW(K2)*Z2POW(K3)*ZC2POW(K4)
            ENDDO
C.............................................................
          ENDDO
C.............................................................
        ENDDO
C.............................................................
      ENDDO
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, EVAL_PHI.
C=============================================================
C IT COMPUTES Z=PHI(ZETA) USING TRUNCATION OF ORDER IORD.
C PHI IS THE DIRECT CONJUGATING FUNCTION OF THE NORMAL FORM.
C
C AUTHOR:    E. TODESCO - INFN & CERN
C
 
      SUBROUTINE EVAL_PHI(IORD,ZETA1,ZETA2,Z1,Z2)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      COMPLEX*16 PHI1,Z1,Z2,ZETA1,ZETAC1,ZETA2,ZETAC2
      COMPLEX*16 ZETA1POW(0:NORMAX),ZETAC1POW(0:NORMAX)
      COMPLEX*16 ZETA2POW(0:NORMAX),ZETAC2POW(0:NORMAX)
      COMMON/DIRETTA/PHI1(NDIM,4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C.............................................................
      ZETAC1=DCONJG(ZETA1)
      ZETAC2=DCONJG(ZETA2)
      ZETA1POW(0)=1.D0
      ZETAC1POW(0)=1.D0
      ZETA2POW(0)=1.D0
      ZETAC2POW(0)=1.D0
C.............................................................
      DO N=1,IORD
        ZETA1POW(N)=ZETA1POW(N-1)*ZETA1
        ZETAC1POW(N)=ZETAC1POW(N-1)*ZETAC1
        ZETA2POW(N)=ZETA2POW(N-1)*ZETA2
        ZETAC2POW(N)=ZETAC2POW(N-1)*ZETAC2
      ENDDO
C.............................................................
      Z1=0.D0
C.............................................................
      DO N=1,IORD
C.............................................................
        DO K1=0,N
C.............................................................
          DO K2=0,N-K1
C.............................................................
            DO K3=0,N-K1-K2
C.............................................................
              K4=N-K1-K2-K3
              Z1=Z1+PHI1(IPUNT(K1,K2,K3,K4),1)
     .          *ZETA1POW(K1)*ZETAC1POW(K2)*ZETA2POW(K3)*ZETAC2POW(K4)
            ENDDO
C.............................................................
          ENDDO
C.............................................................
        ENDDO
C.............................................................
      ENDDO
C.............................................................
      Z2=0.D0
C.............................................................
      DO N=1,IORD
C.............................................................
        DO K1=0,N
C.............................................................
          DO K2=0,N-K1
C.............................................................
            DO K3=0,N-K1-K2
C.............................................................
              K4=N-K1-K2-K3
              Z2=Z2+PHI1(IPUNT(K1,K2,K3,K4),3)
     .          *ZETA1POW(K1)*ZETAC1POW(K2)*ZETA2POW(K3)*ZETAC2POW(K4)
            ENDDO
C.............................................................
          ENDDO
C.............................................................
        ENDDO
C.............................................................
      ENDDO
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MANY_HENON.
C=============================================================
C TRACKS THE GENERALIZED 4D HENON MAP. THE EXPLICIT FORM OF THE MAP IS:
C AMU IS THE ARRAY WITH THE MULTIPOLAR COEFFICIENTS.
C THE SEXTUPOLAR COEFFICIENT IS FIXED TO 1. 
C THE OCTUPOLE IS IN AMU(3), THE DECAPOLE IN AMU(4), ETC. ETC.
C MULTIPOLES UP TO NORMAX ARE INCLUDED
C THE EXPLICIT FORM OF THE MAP IS:
C
C  [X1 ]   [R(2*PI*OMEGA1)][  X                         ]
C  [PX1] = [              ][ X^2-Y^2 + SUM_{J=2} AMU(J)*(X^J-....)]
C  [Y1 ]   [R(2*PI*OMEGA2)][  Y                         ]
C  [PY1]   [              ][  2*X*Y  + SUM_{J=2} AMU(J)*(X^J-....)]
C
C
C AUTHOR: E. TODESCO - INFN.
C
 
      INTEGER FUNCTION MANY_HENON(X,NTURN)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(4)
      COMMON/HENON/C1,S1,C2,S2,AMU(NORMAX),COUPANGLE,
     .             CC11,CC12,CC13,CC14,CC21,CC22,CC23,CC24,
     .             CC31,CC32,CC33,CC34,CC41,CC42,CC43,CC44,
     .             MINHEN,MAXHEN
C.............................................................
      IF (ABS(C1)+ABS(C2)+ABS(S1)+ABS(S2).EQ.0D0) THEN
        WRITE(6,*) '***ERROR(MANY_HENON): THE COMMON BLOCK HENON '
     .             ,'IS EMPTY'
      ENDIF
C.............................................................
      DO I=1,NTURN
        IF(MAXHEN.EQ.2.AND.AMU(2).EQ.1.D0) THEN
          CALL NLHEN(X(1),X(3),SCR1,SCR2)
        ELSE
          CALL NLHENGEN(X(1),X(3),SCR1,SCR2)
        END IF
C.............................................ACTUAL ITERATION
        IF (COUPANGLE.EQ.0D0) THEN
C..........................................NO COUPLING PRESENT
          XS=C1*X(1)+S1*(X(2)+SCR1)     
          PXS=-S1*X(1)+C1*(X(2)+SCR1)
          YS=C2*X(3)+S2*(X(4)+SCR2)
          PYS=-S2*X(3)+C2*(X(4)+SCR2)
C.............................................................
        ELSE
C.............................................COUPLING PRESENT
          XS=CC11*X(1)+CC12*(X(2)+SCR1)+CC13*X(3)+CC14*(X(4)+SCR2)     
          PXS=CC21*X(1)+CC22*(X(2)+SCR1)+CC23*X(3)+CC24*(X(4)+SCR2)
          YS=CC31*X(1)+CC32*(X(2)+SCR1)+CC33*X(3)+CC34*(X(4)+SCR2)
          PYS=CC41*X(1)+CC42*(X(2)+SCR1)+CC43*X(3)+CC44*(X(4)+SCR2)
C.............................................................
        ENDIF
C.............................................................
        X(1)=XS                       !..BACK TO VECTOR X
        X(2)=PXS
        X(3)=YS
        X(4)=PYS
C.............................................................
        IF(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)+X(4)*X(4).GT.1E6) THEN
          MANY_HENON=-I
          GOTO 20
        ENDIF
C.............................................................
      ENDDO
C.............................................................
      MANY_HENON=NTURN
 20   CONTINUE
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, IBIN.
C=============================================================
C  COMPUTES THE BINOMIAL COEFFICIENT IBIN(K,J)= K! /J! / (K-J)!
C
C  AUTHOR: E. TODESCO - INFN
C
C

      INTEGER FUNCTION IBIN(K,J)
C.............................................................
      IBIN=1
      DO L=K-J+1,K
        IBIN=IBIN*L
      END DO
      DO L=1,J
        IBIN=IBIN/L
      END DO
C.............................................................
      END 
CDECK  ID>, NLHEN.
C=============================================================
C GENERATES THE NONLINEARITY OF THE HENON MAP
C
C AUTHOR: E. TODESCO - INFN.
C

      SUBROUTINE NLHEN(XX,YY,SCR1,SCR2)
      IMPLICIT REAL*8(A-H,O-Z)

      SCR1=XX*XX-YY*YY
      SCR2=-2*XX*YY

      END

CDECK  ID>, NLHENGEN.
C=============================================================
C GENERATES THE NONLINEARITY OF A GENERALIZED HENON MAP 
C WHOSE MULTIPOLAR COEFFICIENTS ARE STORED IN /HENON/
C
C AUTHOR: E. TODESCO - INFN.
C

      SUBROUTINE NLHENGEN(XX,YY,SCR1,SCR2)
      PARAMETER(NORMAX=10)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 XPOW(0:NORMAX),YPOW(0:NORMAX)
      COMMON/HENON/C1,S1,C2,S2,AMU(NORMAX),COUPANGLE,
     .             CC11,CC12,CC13,CC14,CC21,CC22,CC23,CC24,
     .             CC31,CC32,CC33,CC34,CC41,CC42,CC43,CC44,
     .             MINHEN,MAXHEN

      XPOW(0)=1
      YPOW(0)=1
      DO JPOW=1,MAXHEN
        XPOW(JPOW)=XPOW(JPOW-1)*XX
        YPOW(JPOW)=YPOW(JPOW-1)*YY
      END DO
      SCR1=0
      SCR2=0
      DO NN=MINHEN,MAXHEN
        IF(AMU(NN).NE.0) THEN
          DO J=NN,0,-4
            SCR1=SCR1+XPOW(J)*YPOW(NN-J)*IBIN(NN,J)*AMU(NN)
            IF(J-1.GE.0) SCR2=SCR2-
     .              XPOW(J-1)*YPOW(NN-J+1)*IBIN(NN,J-1)*AMU(NN)
            IF(J-2.GE.0) SCR1=SCR1-
     .              XPOW(J-2)*YPOW(NN-J+2)*IBIN(NN,J-2)*AMU(NN)
            IF(J-3.GE.0) SCR2=SCR2+
     .              XPOW(J-3)*YPOW(NN-J+3)*IBIN(NN,J-3)*AMU(NN)
          END DO
        END IF
      END DO

      END
CDECK  ID>, CREAHENON.
C=============================================================
C GENERATES A GENERALIZED 4D HENON MAP (MULTIPOLAR KICK) AND STORES IT IN THE 
C ARRAY RMAP. OMEGA1 AND OMEGA2 ARE THE LINEAR FREQUENCIES.
C AMU IS THE ARRAY WITH THE MULTIPOLAR COEFFICIENTS.
C THE SEXTUPOLAR COEFFICIENT IS FIXED TO 1. 
C THE OCTUPOLE IS IN AMU(3), THE DECAPOLE IN AMU(4), ETC. ETC.
C MULTIPOLES UP TO THE NORMAX ARE INCLUDED
C THE EXPLICIT FORM OF THE MAP IS:
C
C  [X1 ]   [R(2*PI*OMEGA1)][  X                         ]
C  [PX1] = [              ][ X^2-Y^2 + SUM_{J=2} AMU(J)*(X^J-....)]
C  [Y1 ]   [R(2*PI*OMEGA2)][  Y                         ]
C  [PY1]   [              ][  2*X*Y  + SUM_{J=2} AMU(J)*(X^J-....)]
C
C IT FILLS THE COMMON BLOCK HENON FOR SUBSEQUENT CALCULATIONS
C
C AUTHOR: E. TODESCO CERN AND INFN (M. GIOVANNOZZI - CERN - HAS
C INTRODUCED SOME CHANGES).
C
 
       SUBROUTINE CREAHENON(OMEGA1,OMEGA2,AMULT,COUPAN)
       PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION AMULT(NORMAX)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/HENON/C1,S1,C2,S2,AMU(NORMAX),COUPANGLE,
     .             CC11,CC12,CC13,CC14,CC21,CC22,CC23,CC24,
     .             CC31,CC32,CC33,CC34,CC41,CC42,CC43,CC44,
     .             MINHEN,MAXHEN
      COMMON/FACTOR/FAC(0:11)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
C.............................................................
      IDIME=4
      CALL INIZ_PUNT(NORMAX)   !..INIZIALIZES POINTER
      CALL MZERO(RMAP,IDIME)   !..SETS MAP TO ZERO
C.............................................................
      AMU(1)=0
      DO J=2,NORMAX
        AMU(J)=AMULT(J)
      END DO
      DO J=NORMAX,2,-1
        IF(AMU(J).NE.0) THEN
          MAXHEN=J
          GOTO 17
        END IF
      END DO
 17   CONTINUE
      DO J=2,NORMAX
        IF(AMU(J).NE.0) THEN
          MINHEN=J
          GOTO 18
        END IF
      END DO
 18   CONTINUE
C.............................................................
      PI2=8D0*DATAN(1D0)
      C1=DCOS(PI2*OMEGA1)
      S1=DSIN(PI2*OMEGA1)
      C2=DCOS(PI2*OMEGA2)
      S2=DSIN(PI2*OMEGA2)
C.............................................................
      COUPANGLE=PI2*COUPAN
      CC=DCOS(COUPANGLE)
      CS=DSIN(COUPANGLE)
C.......................................DEFINES COUPLED MATRIX
      CC11=CC*C1                        !..FIRST ROW
      CC12=-CC*S1
      CC13=-CS*C2
      CC14=CS*S2
C.............................................................
      CC21=CC*S1                        !..SECOND ROW
      CC22=CC*C1
      CC23=-CS*S2
      CC24=-CS*C2
C.............................................................
      CC31=CS*C1                        !..THIRD ROW
      CC32=-CS*S1
      CC33=CC*C2
      CC34=-CC*S2
C.............................................................
      CC41=CS*S1                        !..FOURTH ROW                    
      CC42=CS*C1
      CC43=CC*S2
      CC44=CC*C2
C...............................................LINEAR PART...
      RMAP(5,1)=C1
      RMAP(4,1)=S1
      RMAP(5,2)=-S1
      RMAP(4,2)=C1
      RMAP(3,3)=C2
      RMAP(2,3)=S2
      RMAP(3,4)=-S2
      RMAP(2,4)=C2
C.............................................................
      DO N=MINHEN,MAXHEN 
C...........................................NONLINEAR PART,KICK IN X,PX.........
        DO J=N,0,-2
          RMAP(IPUNT(J,0,N-J,0),1)=S1*IBIN(N,J)*(-1)**((N-J)/2)*AMU(N)
          RMAP(IPUNT(J,0,N-J,0),2)=C1*IBIN(N,J)*(-1)**((N-J)/2)*AMU(N)
        END DO
C...........................................NONLINEAR PART,KICK IN Y,PY.........
        DO J=N-1,0,-2
          RMAP(IPUNT(J,0,N-J,0),3)=S2*IBIN(N,J)*(-1)**((N-J+1)/2)*AMU(N)
          RMAP(IPUNT(J,0,N-J,0),4)=C2*IBIN(N,J)*(-1)**((N-J+1)/2)*AMU(N)
        END DO
      END DO 
C.............................................................
C..........................................COMPUTES FACTORIALS
C..............................USUALLY FACTORIALS ARE COMPUTED
C...............................IN MLATINI, BUT WITH THE HENON
C.................................MAP THIS ROUTINE IS NOT USED
      FAC(0)=1D0
      DO IFAC=1,NORMAX1
        FAC(IFAC)=DFLOAT(IFAC)*FAC(IFAC-1)
      ENDDO
C.............................................................
      END
CDECK  ID>, MRANERR.
C=============================================================
C GENERATES A SEQUENCE OF RANDOM ERRORS TO BE ASSIGNED TO
C SPECIFIED ELEMENTS. THE ERRORS CAN BE ASSIGNED TO ZERO LENGTH
C ELEMENTS ONLY, STARTING FROM THE QUADRUPOLAR COMPONENT.
C NMULT IS THE ORDER OF THE MULTIPOLAR COMPONENT TO BE GENERATED
C RADIUS IS THE RADIUS USED FOR THE ERROR MEASUREMENT. IT IS
C EXPRESSED IN mm.
C PHI IS THE BENDING ANGLE OF THE DIPOLE USED AS A REFERENCE FOR
C THE ERROR EXPRESSED IN mrad.
C RMEAN REPRESENT THE SO CALLED SYSTEMATIC PART OF THE ERROR.
C SIGMA IS THE RANDOM PART OF THE ERROR. IF SIGMA IS ZERO, THEN THE
C SYSTEMATIC ERROR IS ALWAYS ASSIGNED TO THE ELEMENTS AND THE COMMON
C BLOCK RANERR IS NOT FILLED.
C N.B: THE GAUSSIAN DISTRIBUTION USED TO GENERATE THE SEQUENCE OF
C ERRORS IS TRUNCATED AT 3 SIGMA.
C ISKEW IS A FLAG TO SPECIFY WHETHER THE ERROR IS TO BE CONSIDERED
C AS A NORMAL MULTIPOLE OR A SKEW ONE.
C IF ISKEW = 0 THEN THE ELEMENT IS A NORMAL MAGNET.
C IF ISKEW = 1 THEN THE ELEMENT IS A SKEW MAGNET.
C ISEED IS A SEED USED FOR THE RANDOM NUMBER GENERATOR.
C IASSIGN IS A FLAG TO SPECIFY WHETHER THE GENERATED ERRORS HAVE
C TO BE ASSIGNED TO A PHYSICAL ELEMENT OR TO A COMMON BLOCK FOR
C SUBSEQUENT ANALYSIS (REORDERING).
C IASSIGN = 0 THE ARRAY ERR IN THE COMMON BLOCK RANERR IS FILLED.
C IASSIGN = 1 THE ERRORS ARE ASSIGNED TO THE PROPER MULTIPOLAR
C ELEMENT.
C ELEMENT IS A STRING CONTAINING THE NAME OF THE PHYSICAL
C MAGNET WHO IS SUPPOSED TO BE AFFECTED BY THE GENERATED ERRORS.
C
C AUTHOR: M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MRANERR(NMULT,RADIUS,PHIDIP,RMEAN,SIGMA,ISKEW,ISEED,
     .                   IASSIGN,ELEMENT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      CHARACTER*10 ALBL,ELEMENT,STRING
      DIMENSION R(4,4),TEMP(NELE)
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/ACCMAT/RA(4,4,0:NELE),IMAX2,IUPD,IUSE
      COMMON/CSTRANS/CSM(4,4),CSMIN(4,4)
      COMMON/RANERR/ERR(10,NELE,0:1),IND(NELE,10),IEL
C.............................................................
      IF (NMULT.LT.0.OR.NMULT.GT.10) THEN
        WRITE(6,*) '***ERROR(MRANERR): FIRST PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (RADIUS.LT.0D0) THEN
        WRITE(6,*) '***ERROR(MRANERR): SECOND PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (SIGMA.LT.0D0) THEN
        WRITE(6,*) '***ERROR(MRANERR): FIFTH PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (ISKEW.NE.0.AND.ISKEW.NE.1) THEN
        WRITE(6,*) '***ERROR(MRANERR): SIXTH PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (ISEED.LT.0) THEN
        WRITE(6,*) '***ERROR(MRANERR): SEVENTH PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (IASSIGN.NE.0.AND.IASSIGN.NE.1) THEN
        WRITE(6,*) '***ERROR(MRANERR): EIGHTH PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
C.............................................................
      IEL=0         !..INITIALIZES THE NUMBER OF RANDOM ELEMENTS
C.............................................................
      DO J=1,9   !..CHECKS FOR TRAILING BLANK CHARACTERS IN ELEMENT
        IF (ELEMENT(J:J).EQ.' '.AND.ELEMENT(J+1:J+1).NE.' ') THEN
          L1=J+1
          GOTO 10
        ENDIF
      ENDDO
      L1=1
  10  DO J=L1,9
        IF (ELEMENT(J:J).NE.' '.AND.ELEMENT(J+1:J+1).EQ.' ') THEN
          L2=J
          GOTO 20
        ENDIF
      ENDDO
      L2=10
  20  CONTINUE
C.............................................................
      CALL CLTOU(ELEMENT(L1:L2))  !..CONVERTS TO UPPERCASE
C.............................................................
      STRING=ELEMENT(L1:L2)
C.............................................................
      DO JEL=1,IMAX
C.............................................................
        IF (INDEX(ALBL(JEL),STRING).NE.0) THEN
C.............................................................
          IF (IABS(IA(JEL)).EQ.1000) THEN !..THIN QUADRUPOLE
C.............................................................
            IF (NMULT.EQ.1) THEN    !..CHECKS THE CORRISPONDENCE
C.............................................................
              IEL=IEL+1
              IND(IEL,NMULT)=JEL
C.............................................................
            ENDIF
C.............................................................
          ELSEIF (IABS(IA(JEL)).EQ.2000) THEN  !..MULTIPOLE
C.............................................................
            IF (INT(A(JEL,3)).EQ.NMULT.AND.NMULT.GT.1) THEN
C.............................................................
              IEL=IEL+1
              IND(IEL,NMULT)=JEL
C.............................................................
            ENDIF
C.............................................................
          ELSEIF (IABS(IA(JEL)).EQ.(2000+NMULT)) THEN
C.............................................................
            IEL=IEL+1
            IND(IEL,NMULT)=JEL
C.............................................................
          ELSE
C.............................................................
            WRITE(6,*) '***ERROR(MRANERR): SPECIFIED ELEMENT IS NOT OF',
     .                 ' THE PROPER TYPE.'
C.............................................................
          ENDIF
C.............................................................
        ENDIF
C.............................................................
      ENDDO
C.............................................................
      IF (IEL.EQ.0) THEN  !..THE ELEMENT IS NOT PRESENT
C.............................................................
        WRITE(6,*) '***ERROR(MRANERR): SPECIFIED ELEMENT DOES NOT',
     .             ' EXIST.'
        RETURN
C.............................................................
      ELSE
C.............................................................
        IND(IEL+1,NMULT)=ISKEW
C.............................................................
      ENDIF
C.............................................................
      FACT=PHIDIP/(RADIUS**NMULT)
      GRADMU=FACT*RMEAN
      GRADSIG=FACT*SIGMA
C.............................................................
      IF (GRADSIG.GT.1D-20) THEN
C.............................................................
        CALL MGENRAN(ISEED,GRADMU,GRADSIG,IEL,3,TEMP)
C.............................................................
        DO IERR=1,IEL        !..COPIES TEMP IN ERR
C.............................................................
          ERR(NMULT,IERR,ISKEW)=TEMP(IERR)
C.............................................................
        ENDDO
C.............................................................
      ELSE
C.............................................................
        DO JERR=1,IEL
C.............................................................
          ERR(NMULT,JERR,ISKEW)=GRADMU
C.............................................................
        ENDDO
C.............................................................
      ENDIF
C.............................................................
  30  IF (IASSIGN.EQ.1) THEN
C.............................................................
        DO JEL=1,IEL
C.............................................................
          A(IND(JEL,NMULT),1+ISKEW)=ERR(NMULT,JEL,ISKEW)
          IA(IND(JEL,NMULT))=IABS(IA(IND(JEL,NMULT)))
C.............................................................
        ENDDO
C.............................................................
C..THE ERRORS REFER TO LINEAR MULTIPOLAR COMPONENTS. THIS MEANS
C..ALL THE CALCULATIONS FOR THE LINEAR OPTICS HAVE TO BE REDONE.
C.............................................................
        IMAX2=IMAX
        CALL VFRESH2(R)        !..ACCUMULATES THE MATRICES
C.............................................................
        NLMAX2=0
        LIMIT=2000
        IKICK=107
        CALL LUMPIT(NLMAX2,LIMIT,IKICK)   !..LUMPS THE LATTICE
C.............................................................
        IDIME=4
        TOL=1D-6
        NITER=50
        CALL MCLORB(TOL,NITER,IDIME,IDEB) !..COMPUTES THE CLOSED ORBIT
C.............................................................
        IF (NMULT.EQ.1) THEN
          CALL INPS(RA(1,1,IMAX),CSM)   !..COMPUTES THE C-S MAT.
          CALL VISYMP(CSMIN,CSM)   !..COMPUTES THE INV. C-S MAT.
        ENDIF
C.............................................................
        CALL LUMPINPS                     !..C-S COORDINATES
C.............................................................
      ENDIF
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MGENRAN.
C=================================================================
C GENERATES A SEQUENCE OF GAUSSIAN RANDOM NUMBERS.
C ISEED REPRESENTS THE SEED FOR THE GENERATION OF THE NUMBERS.
C RMEAN IS THE MEAN VALUE OF THE DISTRIBUTION.
C SIGMA IS THE STANDARD DEVIATION OF THE DISTRIBUTION.
C IEL REPRESENTS HOW MANY NUMBERS HAVE TO BE GENERATED.
C NSIGMA REPRESENTS THE CUT IN UNITS OF SIGMA:
C IF |X-RMEAN| > SIGMA*NSIGMA THEN X IS REJECTED.
C GAUSS IS THE ARRAY CONTAINING THE IEL RANDOM NUMBERS.
C THE REAL RANDOM NUMBER GENERATOR IS REPRESENTED BY THE ROUTINES:
C RLUXGO (INITIALIZATION OF THE RANDOM NUMBER GENERATOR).
C RNORMX (RANDOM NUMBER GENERATOR).
C FOR MORE DETAILS SEE CNL 213.
C
C AUTHOR: M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MGENRAN(ISEED,RMEAN,SIGMA,IEL,NSIGMA,GAUSS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      DIMENSION GAUSS(NELE)
      REAL RVEC(NELE)
C.............................................................
      EXTERNAL RANLUX
C.............................................................
      IF (ISEED.LT.0) THEN
        WRITE(6,*) '***ERROR(MGENRAN): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (SIGMA.LT.0D0) THEN
        WRITE(6,*) '***ERROR(MGENRAN): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IEL.LT.0.OR.IEL.GT.NELE) THEN
        WRITE(6,*) '***ERROR(MGENRAN): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (NSIGMA.LT.0) THEN
        WRITE(6,*) '***ERROR(MGENRAN): FIFTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (2*IEL.LT.NELE) THEN
        ILEN=2*IEL            !..ACTUAL NUMBER OF RANDOM NUMBERS
      ELSE
        ILEN=NELE
      ENDIF
C.............................................................
      LUX=4
      CALL RLUXGO(LUX,ISEED,0,0)   !..INITIALIZES THE RANDOM SEQUENCE
 10   CALL RNORMX(RVEC,ILEN,RANLUX) !..ACTUAL RANDOM NUMBER GENERATOR
C.............................................................
      IGAUSS=0
      DO JLEN=1,ILEN
        IF (IGAUSS.EQ.IEL) THEN
          GOTO 20
        ELSE
          TEMP=SIGMA*RVEC(JLEN)+RMEAN
          IF (DABS(TEMP-RMEAN).LT.SIGMA*NSIGMA) THEN
            IGAUSS=IGAUSS+1
            GAUSS(IGAUSS)=TEMP
          ENDIF
        ENDIF
C.............................................................
      ENDDO
C.............................................................
      GOTO 10
C.............................................................
 20   RETURN
C.............................................................
      END
CDECK  ID>, MSORT.
C==================================================================
C SORTS A SEQUENCE OF RANDOM ERROR ERRORS STORED IN THE COMMON
C BLOCK RANERR.
C ISEED IS USED TO GENERATE THE RANDOM PERMUTATION OF THE ELEMENTS
C OF THE ARRAY ERR.
C MULTMI AND MULTMA ARE THE MIN AND MAX MULTIPOLAR COMPONENTS
C IN THE COMMON BLOCK IPERM THE ARRAY KPERM CONTAINS THE PERMUTATION
C GENERATED.
C
C AUTHOR: M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MSORT(ISEED,MULTMI,MULTMA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000)
      DIMENSION TEMP0(NELE,10),TEMP1(NELE,10)
      REAL TEMP2(1)
      COMMON/RANERR/ERR(10,NELE,0:1),IND(NELE,10),IEL
      COMMON/IPERM/KPERM(NELE)
C.............................................................
      KEL=0
C.............................................................
      ILEN=1
      LUX=4
      CALL RLUXGO(LUX,ISEED,0,0)  !..INITIALIZES RANDOM NUMBER SEQUENCE
C.............................................................
      DO I=1,100000
        IF (KEL.EQ.IEL) GOTO 20
 10     CALL RANLUX(TEMP2,ILEN)    !..RANDOM NUMBER GENERATOR
        II=1+INT(IEL*TEMP2(1))
        DO KK=1,KEL
          IF (KPERM(KK).EQ.II) GOTO 10
        ENDDO
        KEL=KEL+1
        KPERM(KEL)=II
C.............................................................
        DO MULT=MULTMI,MULTMA
          TEMP0(KEL,MULT)=ERR(MULT,II,0)
          TEMP1(KEL,MULT)=ERR(MULT,II,1)
        ENDDO
C.............................................................
      ENDDO
C.............................................................
 20   DO MULT=MULTMI,MULTMA
C.............................................................
        DO JEL=1,IEL
C.............................................................
          ERR(MULT,JEL,0)=TEMP0(JEL,MULT)
          ERR(MULT,JEL,1)=TEMP1(JEL,MULT)
C.............................................................
        ENDDO
C.............................................................
      ENDDO
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MWRIOPT.
C==================================================================
C WRITES THE LINEAR FUNCTIONS COMPUTED BY THE ROUTINE MLINOPT.
C IUNIT REPRESENTS THE FORTRAN UNIT USED TO WRITE THE LINEAR
C PARAMETERS.
C
C AUTHOR: M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
      SUBROUTINE MWRIOPT(IUNIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000)
      CHARACTER*10 NAME
      COMMON/LINOPT/OPTICS(7,NELE),DIST(NELE),NAME(NELE),IFILL,NREC
C.............................................................
      IF (IUNIT.LT.0) THEN
        WRITE(6,*) '***ERROR(MWRIOPT): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IFILL.NE.19) THEN
        WRITE(6,*) '***ERROR(MWRIOPT): THE COMMON BLOCK LATTICE IS ',
     .             'EMPTY'
        STOP
      ENDIF
C.............................................................
      WRITE(IUNIT,*)
      WRITE(IUNIT,10)
      WRITE(IUNIT,*)
C.............................................................
      DO IEL=1,NREC
        WRITE(IUNIT,20) NAME(IEL),DIST(IEL),(OPTICS(IFUN,IEL),IFUN=1,7)
      ENDDO
C.............................................................
 10   FORMAT(1H ,'  NAME      L-TOTAL  BETAH   ALFAH  ',
     .           '   PHIH  BETAV   ALFAV     PHIV  COUPLING'/
     .       1H ,'              (M)     (M)           ',
     .           '   (RAD)  (M)              (RAD)   (RAD) ')
 20   FORMAT(1X,A10,1X,F8.2,7(1X,F7.2))
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MDUMPLAT.
C==================================================================
C WRITES THE CURRENT LATTICE TO IUNIT 16.
C IDUMP IS A PARAMETER TO CHOOSE THE ACTION.
C IF IDUMP = 0 THEN NO ACTION IS PERFORMED.
C IF IDUMP = 1 THEN THE LATTICE ARRAYS ARE WRITTEN ON UNIT 16.
C COMMENT IS A STRING USED AS COMMENT AT THE TOP OF THE FILE.
C
C AUTHOR: M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MDUMPLAT(IDUMP,COMMENT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000)
      CHARACTER*(*) COMMENT
      CHARACTER*10 ALBL
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
C.............................................................
      IF (IDUMP.NE.0.AND.IDUMP.NE.1) THEN
        WRITE(6,*) '***ERROR(MDUMPLAT): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IDUMP.EQ.1) THEN
        WRITE(16,*) '! ',COMMENT
        DO 1 I=1,IMAX
          CALL WRILAT(16,I)
 1        CONTINUE
        CLOSE(16)
        WRITE(6,*) '***MESSAGE(DUMPLAT): WRITING LATTICE FILE '
C.............................................................
      ENDIF
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MASSIGN.
C=================================================================
C ASSIGNS A SEQUENCE OF RANDOM ERRORS STORED IN THE ARRAY ERR
C (COMMON BLOCK RAN) THE LATTICE ARRAYS. IT IS MEANT TO VERIFY THE
C EFFECTS OF SORTING MAGNETS.
C THE ARRAY ERR(10,NELE,0:1) CONTAINS THE ERRORS: THE FIRST INDEX
C REFERS TO THE MULTIPOLAR COMPONENT (1 QUADRUPOLE...), THE SECOND
C REFERS TO THE ELEMENT NUMBER, WHILE THE THIRD INDICATES WHETHER
C THE ELEMENT IS A NORMAL (0) OR A SKEW ONE (1) (SEE THE ROUTINE
C MRANERR FOR MORE DETAILS).
C IN CASE THE ERRORS REPRESENT QUADRUPOLAR COMPONENTS, THE
C CALCULATIONS OF THE LINEAR OPTICS ARE PERFORMED.
C THE ARRAY IND INDICATES THE CORRISPONDENCE BETWEEN ERRORS AND
C ELEMENTS IN THE COMMON LATICE. 
C IEL IS THE NUMBER OF ERRORS GENERATED.
C MULTMIN REPRESENTS THE LOWEST MULTIPOLAR COMPONENT TO BE ASSIGNED
C MULTMAX REPRESENTS THE HIGHEST MULTIPOLAR COMPONENT TO BE ASSIGNED
C
C AUTHOR: M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MASSIGN(MULTMIN,MULTMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      CHARACTER*10 ALBL
      DIMENSION R(4,4)
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/ACCMAT/RA(4,4,0:NELE),IMAX2,IUPD,IUSE
      COMMON/CSTRANS/CSM(4,4),CSMIN(4,4)
      COMMON/RANERR/ERR(10,NELE,0:1),IND(NELE,10),IEL
C.............................................................
      IF (MULTMIN.LT.0.OR.MULTMIN.GT.MULTMAX) THEN
        WRITE(6,*) '***ERROR(MASSIGN): FIRST PARAMETER OUT OF ',
     .             'BOUNDS'
        STOP
      ENDIF
      IF (MULTMAX.GT.10) THEN
        WRITE(6,*) '***ERROR(MASSIGN): SECOND PARAMETER OUT OF',
     .             'BOUNDS'
        STOP
      ENDIF
C.............................................................
      DO MULT=MULTMIN,MULTMAX
        DO JEL=1,IEL
C.............................................................
          A(IND(JEL,MULT),1)=ERR(MULT,JEL,0)
          A(IND(JEL,MULT),2)=ERR(MULT,JEL,1)
          IA(IND(JEL,MULT))=IABS(IA(IND(JEL,MULT)))
C.............................................................
        ENDDO
      ENDDO
C.............................................................
C..THE ERRORS REFER TO LINEAR MULTIPOLAR COMPONENTS. THIS MEANS
C..ALL THE CALCULATIONS FOR THE LINEAR OPTICS HAVE TO BE REDONE.
C.............................................................
      IMAX2=IMAX
      CALL VFRESH2(R)        !..ACCUMULATES THE MATRICES
C.............................................................
      NLMAX2=0
      LIMIT=2000
      IKICK=107
      CALL LUMPIT(NLMAX2,LIMIT,IKICK)   !..LUMPS THE LATTICE
C.............................................................
      IDIME=4
      TOL=1D-6
      NITER=50
      CALL MCLORB(TOL,NITER,IDIME,IDEB) !..COMPUTES THE CLOSED ORBIT
C.............................................................
      IF (MULTMIN.LE.1) THEN
        CALL INPS(RA(1,1,IMAX),CSM)   !..COMPUTES THE C-S MAT.
        CALL VISYMP(CSMIN,CSM)   !..COMPUTES THE INV. C-S MAT.
      ENDIF
C.............................................................
      CALL LUMPINPS                     !..C-S COORDINATES
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MWRITIME.
C===================================================================
C WRITES THE USED CPU-TIME TO UNIT NU.
C
 
      SUBROUTINE MWRITIME(NU,COMMENT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) COMMENT
      WRITE(NU,'(A,F12.3,2A)')
     .  ' ***MESSAGE(MWRITIME): ',TIMEM(0),' s   ',COMMENT
      RETURN
      END
CDECK  ID>, TIMEM.
C=================================================================
C SYSTEM INDEPENDENT FUNCTION TO COMPUTE THE CPU-TIME ELAPSED
C SINCE THE PREVIOUS CALL. THIS VERSION RUNS ON ALL THE CERN
C PLATFORMS. IT USES THE FUNCTION TIMED OF THE KERNLIB LIBRARY.
C I IS A FLAG:
C I=0  TIMEM REPRESENTS THE ELAPSED TIME SINCE THE PREVIOUS CALL,
C      OTHERWISE IT RETURNS -1.
C
C AUTHOR: M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      DOUBLE PRECISION FUNCTION TIMEM(I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL TMP
      COMMON/ACCTIME/IFLAG
C.............................................................
      IF (I.NE.0) THEN
        TIMEM=-1D0
C.............................................................
      ELSE
C.............................................................
        IF (IFLAG.NE.1) THEN  !..FIRST CALL
C.............................................................
          IFLAG=1
          CALL TIMED(TMP)         !..CALL TO TIMING ROUTINE
          TIMEM=0D0
C.............................................................
        ELSE
C.............................................................
          CALL TIMED(TMP)         !..CALL TO TIMING ROUTINE
          TIMEM=TMP
C.............................................................
        ENDIF
C.............................................................
      ENDIF
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MQFACT.
C=================================================================
C COMPUTES THE QUALITY FACTOR FOR AN ACCELERATOR STRUCTURE. THERE
C ARE FOUR DIFFERENT CHOICES OF QUALITY FACTORS, AND THEY CAN BE
C SELECTED BY USING THE FLAGS IMAP,ITUNE,IRESO1,IRESO2.
C IF IMAP = 1 THE NORM OF THE NONLINEAR PART OF THE MAP IS
C COMPUTED. THE NORM USED IS THE SUM OF MODULI OF THE COEFFICIENTS
C OF THE MAP UP TO ORDER IORD.
C IF ITUNE = 1 THE AVERAGE TUNE OVER THE SUM OF THE INVARIANTS IS
C IS COMPUTED UP TO ORDER (INORM-1)/2 (SEE YELLOW REPORT 94-02 FOR
C MORE DETAILS).
C IF IRESO1 = 1 THE SUM OF THE DRIVING TERMS OF THE RESONANCES IN
C THE INTERVAL [MINRES,MAXRES] AND WITHIN A DISTANCE EPSMAX FROM
C THE LINEAR TUNES IS COMPUTED.
C IF IRESO2 = 1 A MORE REFINED VERSION OF THE PREVIOUS QUALITY
C FACTOR IS USED (NOT IMPLEMENTED YET).
C RADIUS REPRESENTS THE SQUARE ROOT OF THE EMITTANCE. EQUAL
C EMITTANCES IN THE TWO PLANES ARE ASSUMED.
C IORD IS THE ORDER OF THE POLYNOMIAL TRANSFER MAP.
C INORM IS THE ORDER OF THE NORMAL FORM.
C MINRES IS THE MINIMUM ORDER OF THE RESONANCES CONSIDERED IN THE
C QUALITY FACTOR CALCULATIONS (USED ONLY IF IRESO1 OR IRESO2 ARE
C EQUAL TO ONE).
C MAXRES IS THE MAXIMUM ORDER OF THE RESONANCES CONSIDERED IN THE
C QUALITY FACTOR CALCULATIONS (USED ONLY IF IRESO1 OR IRESO2 ARE
C EQUAL TO ONE).
C EPSMAX REPRESENTS THE CUT ON THE DISTANCE FROM THE LINEAR TUNES
C OF THE RESONANCES. ONLY RESONANCES WITHIN A DISC OF RADIUS EPSMAX
C ARE CONSIDERED (USED ONLY IF IRESO1 OR IRESO2 ARE EQUAL TO ONE).
C IOUT IS A FLAG FOR THE OUTPUT.
C IOUT = 1 THE OUTPUT IS ENABLED ON UNIT 6.
C IOUT = 0 THE OUTPUT IS SUPPRESSED.
C IN BOTH CASES THE RESULTS OF THE CALCULATIONS (NAMELY THE QUALITY
C FACTORS) ARE STORED IN THE COMMON BLOC /QFACT/ IN THE VARIABLES:
C QMAP,QTUNE,QRES1,QRES2.
C IOUTNF IS A FLAG TO BE USED BY THE ROUTINE GENERALE_4D. IT SELECTS
C THE OUTPUT TYPE (IF ANY). SEE THE ROUTINE GENERALE_4D FOR MORE
C DETAILS.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MQFACT(IMAP,ITUNE,IRESO1,IRESO2,RADIUS,IORD,INORM,
     .                  MINRES,MAXRES,EPSMAX,IOUT,IOUTNF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      COMMON/FACTOR/FAC(0:11)
      COMMON/QFACT/QMAP(NORMAX1),QTUNE(NORMAX1),
     .             QRES1(0:NORMAX1,-NORMAX1:NORMAX1,NORMAX1),
     .             IFLRES(0:NORMAX1,-NORMAX1:NORMAX1),
     .             QRES2(0:NORMAX1,-NORMAX1:NORMAX1)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
C.............................................................
      IF (IMAP.NE.0.AND.IMAP.NE.1) THEN
         WRITE(6,*) '***ERROR(MQFACT): FIRST ARGUMENT OUT OF BOUNDS'
         STOP
      ENDIF
      IF (ITUNE.NE.0.AND.ITUNE.NE.1) THEN
         WRITE(6,*) '***ERROR(MQFACT): SECOND ARGUMENT OUT OF BOUNDS'
         STOP
      ENDIF
      IF (IRESO1.NE.0.AND.IRESO1.NE.1) THEN
         WRITE(6,*) '***ERROR(MQFACT): THIRD ARGUMENT OUT OF BOUNDS'
         STOP
      ENDIF
      IF (IRESO2.NE.0.AND.IRESO2.NE.1) THEN
         WRITE(6,*) '***ERROR(MQFACT): FOURTH ARGUMENT OUT OF BOUNDS'
         STOP
      ENDIF
      IF (IORD.GT.NORMAX) THEN
         WRITE(6,*) '***ERROR(MQFACT): IORD TOO BIG. MAX ORDER IS ',
     .                      NORMAX
          STOP
      ENDIF
      IF (INORM.GT.NORMAX) THEN
         WRITE(6,*) '***ERROR(MQFACT): INORM TOO BIG. MAX ORDER IS ',
     .                      NORMAX
          STOP
      ENDIF
      IF (IOUT.NE.0.AND.IOUT.NE.1) THEN
         WRITE(6,*) '***ERROR(MQFACT): ELEVENTH ARGUMENT OUT OF',
     .                     ' BOUNDS'
         STOP
      ENDIF
      IF (IOUTNF.LT.0.OR.IOUTNF.GT.3) THEN
         WRITE(6,*) '***ERROR(MQFACT): TWELVTH ARGUMENT OUT OF',
     .                     ' BOUNDS'
         STOP
      ENDIF
C.............................................................
      PI2=8D0*DATAN(1D0)
      PI=PI2*.5D0
C.............................................................
      IF (IMAP.EQ.1) THEN     !..COMPUTES THE MAP NORM
C.............................................................
        IF (IOUT.EQ.1) THEN
C.............................................................
          WRITE(6,*) '***NORM OF THE MAP:                MAP ORDER = '
     .               ,IORD
          WRITE(6,*) '                                   AMPLITUDE = '
     .               ,RADIUS
          WRITE(6,*) '***QUALITY FACTOR                            = '
     .               ,Q1(RADIUS,IORD)
          WRITE(6,*) '***MESSAGE(QFACT): QUALITY FACTOR SUCCESSFULLY',
     .               ' COMPUTED'
C.............................................................
        ENDIF
C.............................................................
      ENDIF
C.............................................................
      IF (ITUNE.EQ.1) THEN     !..COMPUTES THE TUNE-SHIFT
C.............................................................
        CALL PREP_NF_NONRES    !..PREPARES NORMAL FORMS PARAMETERS
C.............................................................
        IF (INORM.EQ.0) THEN
          WRITE(6,*) '***MESSAGE(MQFACT): ORDER OF THE NORMAL FORM ZERO'
          RETURN
        ENDIF
        INTERACT=0
        IINV=0
        IAZAN=1
        CALL GENERALE_4D(RMAP,IORD,INORM,INTERACT,IINV,IOUTNF,
     .                         IAZAN) !..COMPUTES NORMAL FORM
C.............................................................
        IF (IOUT.EQ.1) THEN
C.............................................................
          WRITE(6,*) '***NON RESONANT TUNE-SHIFT:        MAP ORDER = '
     .               ,IORD
          WRITE(6,*) '                           NORMAL FORM ORDER = '
     .               ,INORM
          WRITE(6,*) '                                   AMPLITUDE = '
     .               ,RADIUS
          WRITE(6,*) '***QUALITY FACTOR                            = '
     .               ,Q2I(RADIUS)
          WRITE(6,*) '***MESSAGE(QFACT): QUALITY FACTOR SUCCESSFULLY',
     .               ' COMPUTED'
C.............................................................
        ENDIF
C.............................................................
      ENDIF
C.............................................................
      IF (IRESO1.EQ.1) THEN     !..COMPUTES THE RESONANCES
C.............................................................
        O1=DACOS(RMAP(5,1))/PI2
        O2=DACOS(RMAP(3,3))/PI2
C.............................................................
        DO K1=0,NORMAX1
          DO K2=-NORMAX1,NORMAX1
            IFLRES(K1,K2)=0
            DO JJORD=1,NORMAX1
              QRES1(K1,K2,JJORD)=0D0
            ENDDO
          ENDDO
        ENDDO
C.............................................................
        DO JRES=MINRES,MAXRES
          DO K2=-JRES+1,JRES
            K1=JRES-ABS(K2)
C.............................................................
            EPS=DIRES(O1,O2,K1,K2)
            IF (EPS.LT.EPSMAX) THEN 
              IFLRES(K1,K2)=1
              CALL PREP_NF_SINGRES(O1,O2,K1,K2) !..PREPARES N.F.
C.............................................................
              IF (INORM.EQ.0) THEN
                WRITE(6,*) '***MESSAGE(MQFACT): ORDER OF THE NORMAL '
     .                   ,'FORM ZERO'
                RETURN
              ENDIF
C.............................................................
              INTERACT=0
              IINV=0
              IAZAN=1
              CALL GENERALE_4D(RMAP,IORD,INORM,INTERACT,IINV,IOUTNF,
     .                         IAZAN) !..COMPUTES NORMAL FORM
C.............................................................
              IF (IOUT.EQ.1) THEN
                WRITE(6,*) K1,K2
                WRITE(6,*) 'DISTANCE TO THE RESONANCE', EPS
                WRITE(6,*) 'CONTRIBUTION FROM THE RESONANCE ',
     .                      Q3(RADIUS)
              ENDIF
C.............................................................
              DO JJORD=1,INORM+1
                QRES1(0,0,JJORD)=QRES1(0,0,JJORD)+QRES1(K1,K2,JJORD)
              ENDDO
C.............................................................
            ENDIF
C.............................................................
          ENDDO
C.............................................................
        ENDDO
C.............................................................
        IF (IOUT.EQ.1) THEN
C.............................................................
          WRITE(6,*) '***MESSAGE(QFACT): QUALITY FACTOR '
     .                     ,'SUCCESSFULLY COMPUTED'
          WRITE(6,*) '***RESONANCE WIDTHS:              '
     .                     ,' MAP ORDER = ',IORD
          WRITE(6,*) '                           NORMAL '
     .                     ,'FORM ORDER = ',INORM
          WRITE(6,*) '                                  '
     .                     ,' AMPLITUDE = ',RADIUS
          WRITE(6,*) 'INTERVAL IN RESONANCE ORDER       '
     .                     ,'           = ','[',MINRES,',',MAXRES,']'
          WRITE(6,*) 'MAXIMUM DISTANCE FROM LINEAR TUNES'
     .                     ,'           = ',EPSMAX
        ENDIF
C.............................................................
      ENDIF
C.............................................................
      IF (IRESO2.EQ.1) THEN        !..COMPUTES HYPERVOLUME
C.............................................................
        O1=DACOS(RMAP(5,1))/PI2
        O2=DACOS(RMAP(3,3))/PI2
C.............................................................
        DO K1=0,NORMAX1
          DO K2=-NORMAX1,NORMAX1
            QRES2(K1,K2)=0D0
          ENDDO
        ENDDO
C.............................................................
        DO JRES=MINRES,MAXRES
          DO K2=-JRES+1,JRES
            K1=JRES-ABS(K2)
C.............................................................
            EPS=DIRES(O1,O2,K1,K2)
            IF (EPS.LT.EPSMAX) THEN 
              CALL PREP_NF_SINGRES(O1,O2,K1,K2) !..PREPARES N.F.
C.............................................................
              IF (INORM.EQ.0) THEN
                WRITE(6,*) '***MESSAGE(MQFACT): ORDER OF THE NORMAL '
     .                   ,'FORM ZERO'
                RETURN
              ENDIF
C.............................................................
              INTERACT=0
              IINV=0
              IAZAN=1
              CALL GENERALE_4D(RMAP,IORD,INORM,INTERACT,IINV,IOUTNF,
     .                         IAZAN) !..COMPUTES NORMAL FORM
              CALL FIXLIN(RADIUS,IOUT,HYPVOL)
C.............................................................
              IF (IOUT.EQ.1) THEN
                WRITE(6,*) K1,K2
                WRITE(6,*) 'DISTANCE TO THE RESONANCE', EPS
                WRITE(6,*) 'HYPERVOLUME OF THE RESONANCE ',
     .                      HYPVOL
              ENDIF
C.............................................................
              QRES2(K1,K2)=HYPVOL
              QRES2(0,0)=QRES2(0,0)+HYPVOL
C.............................................................
            ENDIF
C.............................................................
          ENDDO
C.............................................................
        ENDDO
C.............................................................
        IF (IOUT.EQ.1) THEN
C.............................................................
          WRITE(6,*) '***MESSAGE(QFACT): QUALITY FACTOR '
     .                     ,'SUCCESSFULLY COMPUTED'
          WRITE(6,*) '***RESONANCE WIDTHS:              '
     .                     ,' MAP ORDER = ',IORD
          WRITE(6,*) '                           NORMAL '
     .                     ,'FORM ORDER = ',INORM
          WRITE(6,*) '                                  '
     .                     ,' AMPLITUDE = ',RADIUS
          WRITE(6,*) 'INTERVAL IN RESONANCE ORDER       '
     .                     ,'           = ','[',MINRES,',',MAXRES,']'
          WRITE(6,*) 'MAXIMUM DISTANCE FROM LINEAR TUNES'
     .                     ,'           = ',EPSMAX
        ENDIF
C.............................................................
      ENDIF
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MINIT.
C=================================================================
C INITIALIZES ALL THE COMMON BLOCKS USED BY THE OTHER ROUTINES IN
C THE LIBRARY MLAT4.
C IORD IS THE ORDER FOR POLYNOMIAL MANIPULATIONS.
C IDIME IS THE PHASE SPACE DIMENSION.
C FNAME IS THE FILE NAME CONTAINING THE ACCELERATOR STRUCTURE.
C IN CASE FNAME IS A BLANK THEN A MAP IS RED FROM UNIT 10.
C IDEB IS A DEBUGGING FLAG:
C IDEB = 1 DEBUGGING OUTPUT ENABLED.
C IDEB = 0 DEBUGGING OUTPUT SUPPRESSED.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MINIT(IORD,IDIME,FNAME,TITLE,IDEB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      CHARACTER*10 ALBL,ALBL1
      CHARACTER*80 TITLE,FNAME
      DIMENSION B(28),R(4,4)
      COMMON/FACTOR/FAC(0:11)
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/DISP/B,X,Y,Z,SUML,PHI,THETA,PSI,NEL,ALBL1
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
      COMMON/ACCMAT/RA(4,4,0:NELE),IMAX2,IUPD,IUSE
      COMMON/CSTRANS/CSM(4,4),CSMIN(4,4)
C.............................................................
      IF (IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MINIT): ORDER TOO BIG. MAX ORDER IS ',
     .                    NORMAX
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MINIT): SECOND PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (IDEB.NE.0.AND.IDEB.NE.1) THEN
        WRITE(6,*) '***ERROR(MINIT): FIFTH PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
C.............................................................
      CALL INIZ_PUNT(NORMAX) !..INITIALIZES THE POINTER
C.............................................................
      IF (FNAME.EQ.' ') THEN
C.............................................................
        IUNIT=10
        CALL MREAD(RMAP,IORD,IDIME,IUNIT)
C.............................................................
        IF (IDIME.EQ.2) THEN              !..COMPUTES C-S TRANSFORMATION
C.............................................................
          PI2=8D0*DATAN(1D0)
          PI2I=1D0/PI2
          TUNE=PI2I*DACOS(5D-1*(RMAP(IPUNT(1,0,0,0),1)+
     .                          RMAP(IPUNT(0,1,0,0),2))) !..TUNE
          CM=DCOS(PI2*TUNE)
          SM=DSIGN(DSIN(PI2*TUNE),RMAP(IPUNT(1,0,0,0),2))
C...........................................C-S TRANSFORMATION
          CSM(1,1)=DSQRT(RMAP(IPUNT(1,0,0,0),2)/SM)
          CSM(1,2)=0D0
          CSM(2,1)=(CM-RMAP(IPUNT(1,0,0,0),1))/CSM(1,1)/SM
          CSM(2,2)=1D0/CSM(1,1)
C...................................INVERSE C-S TRANSFORMATION
          CSMIN(1,1)=CSM(2,2)
          CSMIN(1,2)=-CSM(1,2)
          CSMIN(2,1)=-CSM(2,1)
          CSMIN(2,2)=CSM(1,1)
C.............................................................
        ELSE
C.............................................................
          ICOUNT=0
          DO I1=0,1
            DO I2=0,1-I1
              DO I3=0,1-I1-I2
                I4=1-I1-I2-I3
                ICOUNT=ICOUNT+1
                DO ICOMP=1,4
                  R(ICOMP,ICOUNT)=RMAP(IPUNT(I4,I3,I2,I1),ICOMP)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
C.............................................................
          CALL INPS(R,CSM)         !..COMPUTES THE C-S MAT.
          CALL VISYMP(CSMIN,CSM)   !..COMPUTES THE INV. C-S MAT.
C.............................................................
        ENDIF
C.............................................................
      ELSE
C.............................................................
        CALL MLATIN(FNAME,TITLE,IDIME)    !..READS LATTICE FILE
C.............................................................
        IMAX2=IMAX
        CALL VFRESH2(R)                   !..ACCUMULATES THE MATRICES
C.............................................................
        NLMAX2=IDEB
        LIMIT=2000
        IKICK=107
        CALL LUMPIT(NLMAX2,LIMIT,IKICK)   !..LUMPS THE LATTICE
C.............................................................
        CALL MRESCALE                     !..RESCALES THE MAP (mm,mrad)
C.............................................................
        TOL=1D-6
        NITER=50
        CALL MCLORB(TOL,NITER,IDIME,IDEB) !..COMPUTES THE CLOSED ORBIT
C.............................................................
        IF (IDIME.EQ.2) THEN              !..COMPUTES C-S TRANSFORMATION
C.............................................................
          PI2=8D0*DATAN(1D0)
          PI2I=1D0/PI2
          TUNE=PI2I*DACOS(5D-1*(RA(1,1,IMAX)+RA(2,2,IMAX))) !..TUNE
          CM=DCOS(PI2*TUNE)
          SM=DSIGN(DSIN(PI2*TUNE),RA(1,2,IMAX))
C...........................................C-S TRANSFORMATION
          CSM(1,1)=DSQRT(RA(1,2,IMAX)/SM)
          CSM(1,2)=0D0
          CSM(2,1)=(CM-RA(1,1,IMAX))/CSM(1,1)/SM
          CSM(2,2)=1D0/CSM(1,1)
C...................................INVERSE C-S TRANSFORMATION
          CSMIN(1,1)=CSM(2,2)
          CSMIN(1,2)=-CSM(1,2)
          CSMIN(2,1)=-CSM(2,1)
          CSMIN(2,2)=CSM(1,1)
C.............................................................
        ELSE
C.............................................................
          CALL INPS(RA(1,1,IMAX),CSM)   !..COMPUTES THE C-S MAT.
          CALL VISYMP(CSMIN,CSM)   !..COMPUTES THE INV. C-S MAT.
C.............................................................
        ENDIF
C.............................................................
        CALL LUMPINPS                     !..C-S COORDINATES
C.............................................................
      ENDIF
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MINITI.
C=================================================================
C INITIALIZES ALL THE COMMON BLOCKS USED BY THE OTHER ROUTINES IN
C THE LIBRARY MLAT4. SIMILAR TO MINIT.
C IORD IS THE ORDER FOR POLYNOMIAL MANIPULATIONS.
C IDIME IS THE PHASE SPACE DIMENSION.
C IOPT IS A FLAG TO SPECIFY THE ACTION TO BE DONE.
C IOPT = 1 THEN A STRUCTURE FILE IS READ FROM UNIT 2 USING THE
C          ROUTINE MLATINI.
C IOPT = 0 THEN A TRANSFER MAP IS RED FROM UNIT 10.
C IDEB IS A DEBUGGING FLAG:
C IDEB = 1 DEBUGGING OUTPUT ENABLED.
C IDEB = 0 DEBUGGING OUTPUT SUPPRESSED.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MINITI(IORD,IDIME,IOPT,TITLE,IDEB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      CHARACTER*10 ALBL,ALBL1
      CHARACTER*80 TITLE
      DIMENSION B(28),R(4,4)
      COMMON/FACTOR/FAC(0:11)
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/DISP/B,X,Y,Z,SUML,PHI,THETA,PSI,NEL,ALBL1
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
      COMMON/ACCMAT/RA(4,4,0:NELE),IMAX2,IUPD,IUSE
      COMMON/CSTRANS/CSM(4,4),CSMIN(4,4)
C.............................................................
      IF (IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MINITI): ORDER TOO BIG. MAX ORDER IS ',
     .                    NORMAX
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MINITI): SECOND PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (IOPT.NE.0.AND.IOPT.NE.1) THEN
        WRITE(6,*) '***ERROR(MINITI): THIRD PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (IDEB.NE.0.AND.IDEB.NE.1) THEN
        WRITE(6,*) '***ERROR(MINITI): FIFTH PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
C.............................................................
      CALL INIZ_PUNT(NORMAX) !..INITIALIZES THE POINTER
C.............................................................
      IF (IOPT.EQ.0) THEN
C.............................................................
        IUNIT=10
        CALL MREAD(RMAP,IORD,IDIME,IUNIT)
C.............................................................
        IF (IDIME.EQ.2) THEN              !..COMPUTES C-S TRANSFORMATION
C.............................................................
          PI2=8D0*DATAN(1D0)
          PI2I=1D0/PI2
          TUNE=PI2I*DACOS(5D-1*(RMAP(IPUNT(1,0,0,0),1)+
     .                          RMAP(IPUNT(0,1,0,0),2))) !..TUNE
          CM=DCOS(PI2*TUNE)
          SM=DSIGN(DSIN(PI2*TUNE),RMAP(IPUNT(1,0,0,0),2))
C...........................................C-S TRANSFORMATION
          CSM(1,1)=DSQRT(RMAP(IPUNT(1,0,0,0),2)/SM)
          CSM(1,2)=0D0
          CSM(2,1)=(CM-RMAP(IPUNT(1,0,0,0),1))/CSM(1,1)/SM
          CSM(2,2)=1D0/CSM(1,1)
C...................................INVERSE C-S TRANSFORMATION
          CSMIN(1,1)=CSM(2,2)
          CSMIN(1,2)=-CSM(1,2)
          CSMIN(2,1)=-CSM(2,1)
          CSMIN(2,2)=CSM(1,1)
C.............................................................
        ELSE
C.............................................................
          ICOUNT=0
          DO I1=0,1
            DO I2=0,1-I1
              DO I3=0,1-I1-I2
                I4=1-I1-I2-I3
                ICOUNT=ICOUNT+1
                DO ICOMP=1,4
                  R(ICOMP,ICOUNT)=RMAP(IPUNT(I4,I3,I2,I1),ICOMP)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
C.............................................................
          CALL INPS(R,CSM)         !..COMPUTES THE C-S MAT.
          CALL VISYMP(CSMIN,CSM)   !..COMPUTES THE INV. C-S MAT.
C.............................................................
        ENDIF
C.............................................................
      ELSE
C.............................................................
        CALL MLATINI(TITLE,IDIME)    !..READS LATTICE FILE
C.............................................................
        IMAX2=IMAX
        CALL VFRESH2(R)                   !..ACCUMULATES THE MATRICES
C.............................................................
        NLMAX2=IDEB
        LIMIT=2000
        IKICK=107
        CALL LUMPIT(NLMAX2,LIMIT,IKICK)   !..LUMPS THE LATTICE
C.............................................................
        CALL MRESCALE                     !..RESCALES THE MAP (mm,mrad)
C.............................................................
        TOL=1D-6
        NITER=50
        CALL MCLORB(TOL,NITER,IDIME,IDEB) !..COMPUTES THE CLOSED ORBIT
C.............................................................
        IF (IDIME.EQ.2) THEN              !..COMPUTES C-S TRANSFORMATION
C.............................................................
          PI2=8D0*DATAN(1D0)
          PI2I=1D0/PI2
          TUNE=PI2I*DACOS(5D-1*(RA(1,1,IMAX)+RA(2,2,IMAX))) !..TUNE
          CM=DCOS(PI2*TUNE)
          SM=DSIGN(DSIN(PI2*TUNE),RA(1,2,IMAX))
C...........................................C-S TRANSFORMATION
          CSM(1,1)=DSQRT(RA(1,2,IMAX)/SM)
          CSM(1,2)=0D0
          CSM(2,1)=(CM-RA(1,1,IMAX))/CSM(1,1)/SM
          CSM(2,2)=1D0/CSM(1,1)
C...................................INVERSE C-S TRANSFORMATION
          CSMIN(1,1)=CSM(2,2)
          CSMIN(1,2)=-CSM(1,2)
          CSMIN(2,1)=-CSM(2,1)
          CSMIN(2,2)=CSM(1,1)
C.............................................................
        ELSE
C.............................................................
          CALL INPS(RA(1,1,IMAX),CSM)   !..COMPUTES THE C-S MAT.
          CALL VISYMP(CSMIN,CSM)   !..COMPUTES THE INV. C-S MAT.
C.............................................................
        ENDIF
C.............................................................
        CALL LUMPINPS                     !..C-S COORDINATES
C.............................................................
      ENDIF
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MCLORB.
C=================================================================
C COMPUTES THE CLOSED ORBIT FOR AN ACCELERATOR STRUCTURE.
C TOL IS THE TOLERANCE IN THE ITERATIVE SEARCH FOR CLOSED ORBIT.
C NITER IS THE ITERATION NUMBER FOR CLOSED ORBIT CALCULATIONS.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C IDEB IS A FLAG FOR DEBUGGING PURPOSES:
C IDEB = 1 DEBUGGING OUTPUT ENABLED
C IDEB = 0 DEBUGGING OUTPUT SUPPRESSED
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MCLORB(TOL,NITER,IDIME,IDEB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      CHARACTER*10 ALBL
      DIMENSION R(4,4),XC(4),XNEW(4)
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
      COMMON/ACCMAT/RA(4,4,0:NELE),IMAX2,IUPD,IUSE
C.............................................................
      EXTERNAL NLFUN
C.............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MCLORB): THIRD PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (IDEB.NE.0.AND.IDEB.NE.1) THEN
        WRITE(6,*) '***ERROR(MCLORB): FOURTH PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
C.............................................................
      IF (IDIME.EQ.2) RETURN    !..NO CLOSED ORBIT FOR 2D
C.............................................................
      IF (IDIME.EQ.4) THEN
C.............................................................
        IMAX2=IMAX
        CALL VFRESH2(R)               !..COMPUTES (1-R)^-1
C.............................................................
        DO ICOMP=1,IDIME             !..INITIALIZATION
          XC(ICOMP)=0D0
          XNEW(ICOMP)=0D0
        ENDDO
        NTURN=1
C.............................................................
        DO ITER=1,NITER              !..STARTS ITERATIVE PROCESS
          ITURN=MANYTURN(XC,NTURN,IDIME,NLFUN)  !..TRACKING
C.............................................................
          DO IROW=1,IDIME            !..NEW CLOSED ORBIT
            DO ICOL=1,IDIME
              XNEW(IROW)=XNEW(IROW)+RA(IROW,ICOL,0)*XC(ICOL)
            ENDDO
          ENDDO
C.............................................................
          ANORM=0D0                  !..COMPUTES THE NORM OF THE SOLUTION
          DO I=1,IDIME
            ANORM=ANORM+DABS(XNEW(I))
          ENDDO
C.............................................................
          IF (ANORM.LT.TOL) THEN     !..SOLUTION FOUND
            WRITE(6,*) '***MESSAGE(MCLORB): CLOSED ORBIT FOUND IN ',
     .                 ITER,' ITERATIONS'
            GOTO 10
          ENDIF
C.............................................................
          DO I=1,IDIME
            XC(I)=XC(I)+XNEW(I)
          ENDDO
C.............................................................
        ENDDO
C.............................................................
        WRITE(6,*) '***MESSAGE(MCLORB): CONVERGENCE NOT ACHIEVED'
        RETURN
C.............................................................
 10     IF (IA(1).EQ.107) THEN  !..MODIFIES THE LATTICE DUE TO CLOSED ORBIT
C.............................................................
          A(1,1)=A(1,1)+XC(1) !..ADDS THE CLOSED ORBIT TO THE FIRST KICK
          A(1,3)=A(1,3)+XC(3)
          ALBL(1)='NEW KICK'
C.............................................................
        ELSE
C.............................................................
          IF (IMAX+1.GT.NELE) THEN
            WRITE(6,*) '***MESSAGE(MCLORB): TOO MANY ELEMENTS'
            RETURN
          ELSE
C.............................................................
            DO IEL=IMAX,1,-1    !..SHIFTS ALL THE ELEMENTS BY ONE
              IA(IEL+1)=IA(IEL)
              ALBL(IEL+1)=ALBL(IEL)
              DO JEL=1,4
                A(IEL+1,JEL)=A(IEL,JEL)
              ENDDO
            ENDDO
C.............................................................
            IA(1)=107           !..DEFINES THE CLOSED ORBIT KICK
            DO JEL=1,4
              A(1,JEL)=XC(JEL)
            ENDDO
            ALBL(1)='CLORB KICK'
C.............................................................
            IMAX=IMAX+1
            IMAX2=IMAX
C.............................................................
          ENDIF
C.............................................................
        ENDIF
C.............................................................
      ENDIF
C.............................................................
      CALL VFRESH2(R)          !..ACCUMULATE NEW MATRICES
C.............................................................
      NLMAX2=IDEB
      IKICK=107
      LIMIT=2000
      CALL LUMPIT(NLMAX2,LIMIT,IKICK)  !..LUMPS THE NEW STRUCTURE
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MAPDYN1.
C=================================================================
C COMPUTES THE BORDER OF THE STABILITY DOMAIN (DYNAMIC APERTURE).
C IMAP IS A PARAMETER SPECIFING THE MAP TYPE:
C IMAP = 1 ===> HENON MAP
C IMAP = 2 ===> ACCELERATOR MODEL
C NTURN IS THE TURN NUMBER USED TO COMPUTE THE STABILITY.
C TOL IS THE TOLERANCE OF THE BISECTION METHOD USED TO DETERMINE
C THE BORDER OF THE STABILITY DOMAIN .
C IDIME IS THE PHASE SPACE DIMENSION.
C (ANG1,ANG2) REPRESENTS THE INTERVAL FOR ANGLES (IN UNITS OF 2PI).
C NANG INITIAL CONDITIONS ARE CONSIDERED.
C BORDER IS THE ARRAY CONTAINING THE DYNAMIC APERTURE FOR THE
C DIFFERENT ANGLES.
C RADIUS REPRESENTS THE MINIMUM OF BORDER OVER THE ANGLES.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MAPDYN1(IMAP,NTURN,TOL,IDIME,ANG1,ANG2,NANG,BORDER,
     .                   RADIUS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      PARAMETER (NMAX=5000,RHOMAX=2D1)
      DIMENSION X(4),BORDER(NMAX,2),RHO1(100),RHO2(100)
C.............................................................
      EXTERNAL NLFUN
C.............................................................
      IF (IMAP.NE.1.AND.IMAP.NE.2) THEN
        WRITE(6,*) '***ERROR(MAPDYN1): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (NANG.GT.NMAX) THEN
        WRITE(6,*) '***ERROR(MAPDYN1): NUMBER OF ANGLES TOO BIG.'
        WRITE(6,*) '***MESSAGE(MAPDYN1): FIFTH PARAMETER MUST BE < ',
     .              NMAX
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MAPDYN1): FOURTH PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
C.............................................................
      IPOINT=0                          !..INITIALIZATION
      IFLAG=0
      PI2=8D0*DATAN(1D0)
      ASTEP=(ANG2-ANG1)/NANG
C.............................................................
      IF (ASTEP.LT.0D0) THEN
        ASTEP=-ASTEP
        WRITE(6,*) '***MESSAGE(MAPDYN1): WRONG ORDER IN INTERVAL'
        WRITE(6,*) '***MESSAGE(MAPDYN1): ORDER REVERSED'
      ENDIF
C.............................................................
      RHOMAX1=RHOMAX
C.............................................................
 10   X(1)=RHOMAX1                       !..CHECKS THE INITIAL GUESS
      X(2)=0D0
      X(3)=0D0
      X(4)=0D0
C.............................................................
      IF (IMAP.EQ.1) THEN
        ITURN=MANY_HENON(X,NTURN)
      ELSE
        ITURN=MANYTURN(X,NTURN,IDIME,NLFUN)      !..TRACKING
      ENDIF
C.............................................................
      IF (ITURN.EQ.NTURN) THEN         !..CHECKS THE STABILITY
C.............................................................
        RHOMAX1=2D0*RHOMAX1
        IFLAG=1
        WRITE(6,*) '***MESSAGE(MAPDYN1): INITIAL GUESS TOO SMALL'
        GOTO 10
C.............................................................
      ENDIF
C.............................................................
      IF (IFLAG.EQ.1) WRITE(6,*) '***MESSAGE(MAPDYN1): FINAL GUESS ',
     .                RHOMAX1
C.............................................................
      RADIUS=RHOMAX             !..INITIALIZES RADIUS
C.............................................................
      DO IANG=1,NANG                 !..LOOP ON THE ANGLES
        ANG=PI2*(ANG1+ASTEP*DFLOAT(IANG))
C.............................................................
        RHO1(1)=0D0
        RHO2(1)=RHOMAX1
C.............................................................
        DO ITER=1,99                 !..BISECTION LOOP
C.............................................................
          IF (DABS(RHO2(ITER)-RHO1(ITER)).LT.TOL)
     .        THEN   !..CHECKS TOLERANCE
C.............................................................
            BORDER(IANG,1)=RHO1(ITER)*DCOS(ANG)
            BORDER(IANG,2)=RHO1(ITER)*DSIN(ANG)
C.............................................................
            RADIUS=DMIN1(RADIUS,RHO1(ITER))
C.............................................................
            GOTO 20
C.............................................................
          ELSE
C.............................................................
            RHONEW=(RHO2(ITER)-RHO1(ITER))*5D-1+RHO1(ITER)
C.............................................................
            IF (IDIME.EQ.2) THEN
C.............................................................
              X(1)=RHONEW*DCOS(ANG)
              X(2)=RHONEW*DSIN(ANG)
              X(3)=0D0
              X(4)=0D0
C.............................................................
            ELSE
C.............................................................
              X(1)=RHONEW*DCOS(ANG)
              X(2)=0D0
              X(3)=RHONEW*DSIN(ANG)
              X(4)=0D0
C.............................................................
            ENDIF
C.............................................................
            IF (IMAP.EQ.1) THEN
              ITURN=MANY_HENON(X,NTURN)
            ELSE
              ITURN=MANYTURN(X,NTURN,IDIME,NLFUN)
            ENDIF
C.............................................................
            IF (ITURN.EQ.NTURN) THEN
C.............................................................
              RHO1(ITER+1)=RHONEW
              RHO2(ITER+1)=RHO2(ITER)
C.............................................................
            ELSE
C.............................................................
              RHO1(ITER+1)=RHO1(ITER)
              RHO2(ITER+1)=RHONEW
C.............................................................
            ENDIF
C.............................................................
          ENDIF
C.............................................................
        ENDDO
C.............................................................
 20     CONTINUE
C.............................................................
      ENDDO
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MLATIN.
C=================================================================
C READS A STRUCTURE FILE CALLED FNAME.
C TITLE IS A STRING RED FROM THE MAD FILE.
C IDIME IS THE DIMENSION OF THE PHASE SPACE:
C
C IDIME = 2 ALL THE SOURCES OF COUPLING ARE SUPPRESSED AND ONLY
C           THE HORIZONTAL MOTION IS CONSIDERED.
C
C IDIME = 4 FULL 4D MOTION WITH COUPLING.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MLATIN(FNAME,TITLE,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      CHARACTER*4 KEYWORD
      CHARACTER*8 PROGVRSN,DATAVRSN,DATE,TIME,JOBNAME
      CHARACTER*10 ALBL,ALBL1
      CHARACTER*16 NAME
      CHARACTER*80 TITLE,FNAME
      LOGICAL SYMM
      DIMENSION B(28)
      INTRINSIC INDEX
      INTEGER SUPER
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/DISP/B,X,Y,Z,SUML,PHI,THETA,PSI,NEL,ALBL1
      COMMON/FACTOR/FAC(0:11)
C.............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MLATIN): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.......................................OPEN FILE AND READ HEADER
      OPEN(UNIT=2,FILE=FNAME,STATUS='OLD')
      READ(2,'(5A8,I8,L8,I8/A80)')
     .  PROGVRSN,DATAVRSN,DATE,TIME,JOBNAME,SUPER,SYMM,NPOS,
     .  TITLE
      WRITE(6,*) PROGVRSN,'...',DATE,'...',TIME
      WRITE(6,*) DATAVRSN
      WRITE(6,*) TITLE
      WRITE(6,*) NPOS,' RECORDS EXPECTED'
C
      IF ((INDEX(DATAVRSN,'STRUCT')+INDEX(DATAVRSN,'STRUCT')).EQ.0)
     .  STOP '***ERROR(MLATIN): WRONG FILE FORMAT'
C..............................................DETERMINE THE ENERGY
        ENERGY=450D0  ! GEV PROTONS
        IF (((INDEX(TITLE,'LHC')+INDEX(TITLE,'lhc')).NE.0.)
     .    .AND.((INDEX(TITLE,'COL')+INDEX(TITLE,'col')).NE.0)) THEN
          ENERGY=8000D0       ! GEV
        ELSEIF (((INDEX(TITLE,'LHC')+INDEX(TITLE,'lhc')).NE.0)
     .    .AND.((INDEX(TITLE,'INJ')+INDEX(TITLE,'inj')).NE.0)) THEN
          ENERGY=450D0        ! GEV
        ELSEIF (((INDEX(TITLE,'LEP')+INDEX(TITLE,'lep')).NE.0)
     .    .AND.((INDEX(TITLE,'COL')+INDEX(TITLE,'col')).NE.0)) THEN
          ENERGY=45.6D0       ! GEV
        ELSEIF (((INDEX(TITLE,'LEP')+INDEX(TITLE,'lep')).NE.0)
     .    .AND.((INDEX(TITLE,'INJ')+INDEX(TITLE,'inj')).NE.0)) THEN
          ENERGY=20D0         ! GEV
        ELSE
 7000     ENERGY=1D0
        ENDIF
C.......................................................OUTPUT FILE
      NEL=0                             !..INITIALIZES COUNTER
C.............................................INITIALIZES ARRAYS
      DO JJ=1,4
        DO II=1,NELE
          A(II,JJ)=0D0
          IA(II)=0
        ENDDO
      ENDDO
C..........................................FACTORIALS
      FAC(0)=1D0
      DO JFACT=1,11
        FAC(JFACT)=FLOAT(JFACT)*FAC(JFACT-1)
      ENDDO
C.................LOOP OVER NPOS ELEMENTS, EACH RECORD 4 LINES LONG
      DO I=1,NPOS
C..............................................READ A SINGLE RECORD
        READ(2,'(A4,A16,F12.6,3E16.9)',END=1000)
     .    KEYWORD,NAME,(B(II),II=1,4)
C......................................................PROCESS DATA
        ALBL1=NAME(1:10)
        III=INDEX(ALBL1,' ')-1
        IF (INDEX(KEYWORD,'DRIF').EQ.1) THEN             !..DRIFT
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
        ELSEIF (INDEX(KEYWORD,'MARK').EQ.1) THEN         !..MARKER
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IF(INDEX(ALBL1,'CHARGEKICK').EQ.1) THEN        !..SPACE CHARGE KICK
            IA(NEL)=2300
          ELSE
            IA(NEL)=0
          ENDIF
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
        ELSEIF (INDEX(KEYWORD,'RBEN').EQ.1) THEN         !..RECT BEND
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          IF (ABS(X).NE.0D0.OR.ABS(Y).NE.0D0.OR.ABS(Z).NE.0D0.OR.
     .    ABS(THETA).NE.0D0.OR.ABS(PHI).NE.0D0)
     .    WRITE(6,*) '***ERROR(MLATIN): DISPLACEMENTS AT ',ALBL1,
     .    ' IGNORED'
C         IF (ABS(Z).GT.1D-10) THEN              !..LONGITUDINAL DISP.
C           NEL=NEL+1                           !..NOT IMPLEMENTED YET
C         IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
C           IA(NEL)=3
C           A(NEL,1)=-Z
C           A(NEL,4)=ENERGY
C           ALBL(NEL)=ALBL1
C         ENDIF
          IF (ABS(B(5)+PSI).GT.1D-10) THEN       !..TILT ANGLE
            NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,1)=(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='RBEND TILT '
          ENDIF
          IF (ABS(B(8)).GT.1D-10)                               !..H1
     .      WRITE(6,*) '***ERROR(MLATIN): RBEND  H1  AT ',ALBL1
          IF (ABS(B(9)).GT.1D-10)                               !..H2
     .      WRITE(6,*) '***ERROR(MLATIN): RBEND  H2  AT ',ALBL1
          IF (ABS(B(6)).GT.1D-10) THEN                          !..E1
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=2
            A(NEL,2)=B(6)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='WEDGE      '
          ENDIF
          IF (ABS(B(3)).GT.1D-10) THEN
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=2                           !..WEDGE OF THE RECT
            A(NEL,1)=0D0
            A(NEL,2)=0D0
            A(NEL,3)=-1D0
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=4
            A(NEL,1)=B(1)
            A(NEL,2)=ENERGY*B(2)/(0.0299792458D0*B(1))
            A(NEL,3)=B(3)*ENERGY/0.0299792458D0
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=2                           !..WEDGE OF THE RECT
            A(NEL,1)=0D0
            A(NEL,2)=0D0
            A(NEL,3)=-1D0
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
          ELSE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=44
            A(NEL,1)=B(1)
            A(NEL,2)=ENERGY*B(2)/(0.0299792458D0*B(1))
            A(NEL,3)=0D0
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
          ENDIF
          IF (ABS(B(7)).GT.1D-10) THEN                      !..E2
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=2
            A(NEL,2)=B(7)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='WEDGE      '
          ENDIF
          IF (ABS(B(5)+PSI).GT.1D-10) THEN              !..TILT ANGLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=-(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='RBEND TILT '
          ENDIF
C         IF (ABS(Z).GT.1D-10) THEN              !..LONGITUDINAL DISP.
C           NEL=NEL+1                           !..NOT IMPLMENTED YET
C           IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
C           IA(NEL)=3
C           A(NEL,1)=Z
C           A(NEL,4)=ENERGY
C           ALBL(NEL)=ALBL1
C         ENDIF
        ELSEIF (INDEX(KEYWORD,'SBEN').EQ.1) THEN        !..SECTOR BEND
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          IF (ABS(X).NE.0D0.OR.ABS(Y).NE.0D0.OR.ABS(Z).NE.0D0.OR.
     .    ABS(THETA).NE.0D0.OR.ABS(PHI).NE.0D0)
     .    WRITE(6,*) '***ERROR(MLATIN): DISPLACEMENTS AT ',ALBL1,
     .    ' IGNORED'
C         IF (ABS(Z).GT.1D-10) THEN              !..LONGITUDINAL DISP.
C           NEL=NEL+1                           !..NOT IMPLEMENTED YET
C           IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
C           IA(NEL)=3
C           A(NEL,1)=-Z
C           A(NEL,4)=ENERGY
C           ALBL(NEL)=ALBL1
C         ENDIF
          IF (ABS(B(5)+PSI).GT.1D-10) THEN              !..TILT ANGLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='SBEND TILT '
          ENDIF
          IF (ABS(B(8)).GT.1D-10)                                !..H1
     .      WRITE(6,*) '***ERROR(MLATIN): SBEND  H1  AT ',ALBL1
          IF (ABS(B(9)).GT.1D-10)                                !..H2
     .      WRITE(6,*) '***ERROR(MLATIN): SBEND  H2  AT ',ALBL1
          IF (ABS(B(6)).GT.1D-10) THEN                           !..E1
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=2
            A(NEL,2)=B(6)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='WEDGE      '
          ENDIF
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=4
          A(NEL,1)=B(1)
          A(NEL,2)=ENERGY*B(2)/(0.0299792458D0*B(1))
          A(NEL,3)=B(3)*ENERGY/0.0299792458D0
          A(NEL,4)=ENERGY
          IF (ABS(B(7)).GT.1D-10) THEN                           !..E2
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=2
            A(NEL,2)=B(7)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='WEDGE      '
          ENDIF
          IF (ABS(B(5)).GT.1D-10) THEN                  !..TILT ANGLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=-(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='SBEND TILT '
          ENDIF
C         IF (ABS(Z).GT.1D-10) THEN              !..LONGITUDINAL DISP.
C           NEL=NEL+1                           !..NOT IMPLEMENTED YET
C           IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
C           IA(NEL)=3
C           A(NEL,1)=Z
C           A(NEL,4)=ENERGY
C           ALBL(NEL)=ALBL1
C         ENDIF
        ELSEIF (INDEX(KEYWORD,'QUAD').EQ.1) THEN       !..QUADRUPOLE
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          CALL MISAL(1,ENERGY)                !..MISALIG.
          IF (ABS(B(5)+PSI).GT.1D-10) THEN       !..ROLLED QUAD
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=55
            A(NEL,3)=(B(5)+PSI)*57.295779513082D0
          ELSE                               !..UPRIGHT QUAD
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=5
            A(NEL,3)=1D0
          ENDIF
          A(NEL,1)=B(1)
          A(NEL,2)=ENERGY*B(3)/0.0299792458D0
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
          CALL MISAL(-1,ENERGY)               !..MISALIG. (INVERSE)
        ELSEIF (INDEX(KEYWORD,'SEXT').EQ.1) THEN      !..SEXTUPOLE
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          CALL MISAL(1,ENERGY)                !..MISALIG.
          IF (ABS(B(5)+PSI).GT.1D-10) THEN              !..TILT ANGLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='SEXT TILT '
          ENDIF
          IF (ABS(B(1)).GE.1D-10) THEN   !..THICK SEXTUPOLE (DR+SEX+DR)
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=5D-1*B(1)
            A(NEL,4)=ENERGY
            ALBL(NEL)='SEXT DRIFT'
          ENDIF
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=2002
          IF (B(4).EQ.0.) IA(NEL)=-2002      !..FLAG FOR LINEAR BLOCS
          A(NEL,1)=B(4)*B(1)/FAC(2)          !..K_2
          A(NEL,2)=0D0
          A(NEL,3)=2D0
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
          IF (ABS(B(1)).GE.1D-10) THEN   !..THICK SEXTUPOLE (DR+SEX+DR)
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=5D-1*B(1)
            A(NEL,4)=ENERGY
            ALBL(NEL)='SEXT DRIFT'
          ENDIF
          IF (ABS(B(5)+PSI).GT.1D-10) THEN      !..TILT SEXTUPOLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=-(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='SEXT TILT '
          ENDIF
          CALL MISAL(-1,ENERGY)               !..MISALIG. (INVERSE)
        ELSEIF (INDEX(KEYWORD,'OCTU').EQ.1) THEN     !..OCTUPOLE
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          CALL MISAL(1,ENERGY)                !..MISALIG.
          IF (ABS(B(5)+PSI).GT.1D-10) THEN      !..TILT OCTUPOLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='OCTUP TILT'
          ENDIF
          IF (ABS(B(1)).GE.1D-10) THEN   !..THICK OCTUPOLE (DR+OCT+DR)
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=5D-1*B(1)
            ALBL(NEL)=ALBL1
          ENDIF
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=2003
          IF (B(6).EQ.0.) IA(NEL)=-2003      !..FLAG FOR LINEAR BLOCS
          A(NEL,1)=B(6)*B(1)/FAC(3)          !..K_3
          A(NEL,2)=0D0
          A(NEL,3)=3D0
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
          IF (ABS(B(1)).GE.1D-10) THEN   !..THICK OCTUPOLE (DR+OCT+DR)
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=5D-1*B(1)
            ALBL(NEL)='OCTU DRIFT'
          ENDIF
          IF (ABS(B(5)+PSI).GT.1D-10) THEN      !..TILT OCTUPOLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=-(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='OCTUP TILT'
          ENDIF
          CALL MISAL(-1,ENERGY)               !..MISALIG. (INVERSE)
        ELSEIF (INDEX(KEYWORD,'MULT').EQ.1) THEN      !..MULTIPOLE
          BACKSPACE(2)
          READ(2,'(A4,A16,F12.6,2E16.9,I16)')
     .    KEYWORD,NAME,(B(II),II=1,3),NMULT
          READ(2,'(5E16.9)') (B(II),II=5,5+NMULT) !..K_N
          READ(2,'(5E16.9)') (B(II),II=6+NMULT,6+2*NMULT) !..J_N
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          CALL MISAL(1,ENERGY)                !..MISALIG.
          IF (ABS(PSI).GT.1D-10) THEN      !..TILT MULTIPOLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=PSI*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='MULTI TILT'
          ENDIF
          IF (ABS(B(2)+B(5)+B(3)+B(6+NMULT)).NE.0D0) THEN
            NEL=NEL+1                                !..DIPOLE KICK
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=107
            A(NEL,1)=B(2)+B(5)
            A(NEL,2)=0D0
            A(NEL,3)=B(3)+B(6+NMULT)
            A(NEL,4)=0D0
            ALBL(NEL)=ALBL1
          ENDIF
          NEL=NEL+1                                !..QUADRUPOLE KICK
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=-1000
          IF (ABS(B(6)+B(NMULT+7)).NE.0D0) IA(NEL)=1000
          A(NEL,1)=B(6)
          A(NEL,2)=B(NMULT+7)
          A(NEL,3)=1
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
          DO JMULT=0,1                   !..MULTIPOLE COMPONENTS
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=-2000
            A(NEL,1)=0D0
            A(NEL,2)=0D0
            A(NEL,3)=JMULT
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
          ENDDO
          DO JMULT=2,NMULT                   !..MULTIPOLE COMPONENTS
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=2000
            A(NEL,1)=B(JMULT+5)/FAC(JMULT)
            A(NEL,2)=B(NMULT+JMULT+6)/FAC(JMULT)
            A(NEL,3)=JMULT
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
            IF ((ABS(A(NEL,1))+ABS(A(NEL,2))).EQ.0D0) IA(NEL)=-2000
          ENDDO
          IF (ABS(PSI).GT.1D-10) THEN      !..TILT MULTIPOLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=-PSI*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='MULTI TILT'
          ENDIF
          CALL MISAL(-1,ENERGY)               !..MISALIG. (INVERSE)
        ELSEIF (INDEX(KEYWORD,'MONI').EQ.1.) THEN   ! MONITOR
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=5D-1*B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)='DRIFT_MON'
C............................................................
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=103                !..SPECIAL CARD FOR V.Z.
          ALBL(NEL)=ALBL1
C............................................................
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=5D-1*B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)='DRIFT_MON'
        ELSEIF (INDEX(KEYWORD,'HMON').EQ.1) THEN          !..BPM
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=5D-1*B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)='DRIFT_MON'
C............................................................
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=101                !..SPECIAL CARD FOR V.Z.
          ALBL(NEL)=ALBL1
C............................................................
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=5D-1*B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)='DRIFT_MON'
        ELSEIF (INDEX(KEYWORD,'VMON').EQ.1) THEN          !..BPM
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=5D-1*B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)='DRIFT_MON'
C............................................................
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=102                !..SPECIAL CARD FOR V.Z.
          ALBL(NEL)=ALBL1
C............................................................
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=5D-1*B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)='DRIFT_MON'
        ELSEIF (INDEX(KEYWORD,'KICK').EQ.1.OR.      !..ORBIT CORRECTOR
     .          INDEX(KEYWORD,'HKIC').EQ.1.OR.
     .          INDEX(KEYWORD,'VKIC').EQ.1) THEN
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          IF (ABS(B(1)).GT.1D-10) THEN
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=5D-1*B(1)
            A(NEL,4)=ENERGY
            ALBL(NEL)='DRIFT_COR'
          ENDIF
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=107
          ALBL(NEL)=ALBL1
          IF (ABS(B(1)).GT.1D-10) THEN
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=5D-1*B(1)
            A(NEL,4)=ENERGY
            ALBL(NEL)='DRIFT_COR'
          ENDIF
        ELSEIF (INDEX(KEYWORD,'RFCA').EQ.1) THEN       !..RF-CAVITY
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
C         CALL MISAL(1,ENERGY)                !..MISALIG.
C         IF (ABS(PSI).GT.1D-10) THEN     !..IGNORED FOR THE MOMENT
C           NEL=NEL+1
C           IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
C           IA(NEL)=20
C           A(NEL,2)=PSI*57.295779513082D0
C           A(NEL,4)=ENERGY
C           ALBL(NEL)='CAVIT TILT'
C         ENDIF
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=31
          A(NEL,1)=B(1)                  !..LENGTH
          A(NEL,2)=B(7)                  !..VOLTAGE IN MV
          A(NEL,3)=B(8)*360D0            !..LAG IN DEGREE
          A(NEL,4)=B(6)                  !..FREQUENCY IN MHZ
          ALBL(NEL)=ALBL1
C         IF (ABS(PSI).GT.1D-10) THEN     !..IGNORED FOR THE MOMENT
C           NEL=NEL+1
C           IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
C           IA(NEL)=20
C           A(NEL,2)=PSI*57.295779513082D0
C           A(NEL,4)=ENERGY
C           ALBL(NEL)='CAVIT TILT'
C         ENDIF
C         CALL MISAL(-1,ENERGY)               !..MISALIG. (INVERSE)
C         WRITE(2,'(I5,F12.7,G12.5,G13.6,G18.10,1X,A1,A10,A1)')
C    .      IA,(A(II),II=1,4),'''',ALBL,''''
        ELSEIF (INDEX(KEYWORD,'RCOL').EQ.1) THEN    !..COLLIMATOR
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=2100  !  RECTANGULAR COLLIMATOR
          A(NEL,1)=B(2)
          A(NEL,2)=B(3)
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
          IF (B(2).LT.1D-20.OR.B(3).LT.1D-20) IA(NEL)=-2100
          IF (B(1).GT.0D0) THEN   !  ADD A DRIFT AND A COLLIMATOR
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=B(1)
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=2100  !  RECTANGULAR COLLIMATOR
            A(NEL,1)=B(2)
            A(NEL,2)=B(3)
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
            IF (B(2).LT.1D-20.OR.B(3).LT.1D-20) IA(NEL)=-2100
          ENDIF
        ELSEIF (INDEX(KEYWORD,'ECOL').EQ.1) THEN    !..COLLIMATOR
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=2200  !  ELLIPTICAL COLLIMATOR
          A(NEL,1)=B(2)*B(2)
          A(NEL,2)=B(3)*B(3)
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
          IF (B(2).LT.1D-20.OR.B(3).LT.1D-20) IA(NEL)=-2200
          IF (B(1).GT.0D0) THEN   !  ADD A DRIFT AND A COLLIMATOR
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=B(1)
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=2200  !  ELLIPTICAL COLLIMATOR
            A(NEL,1)=B(2)*B(2)
            A(NEL,2)=B(3)*B(3)
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
            IF (B(2).LT.1D-20.OR.B(3).LT.1D-20) IA(NEL)=-2200
          ENDIF
        ELSEIF (INDEX(KEYWORD,'ELSE').EQ.1) THEN  !..E-STATIC SEPARATOR
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
        ELSEIF (INDEX(KEYWORD,'    ').EQ.1) THEN
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
        ELSEIF (INDEX(KEYWORD,'SOLE').EQ.1) THEN     !..SOLENOID
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          CALL MISAL(1,ENERGY)                !..MISALIG.
          IF (ABS(PSI).GT.1D-10) THEN
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=PSI*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='SOLE TILT '
          ENDIF
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=19
          A(NEL,1)=B(1)
          A(NEL,2)=ENERGY*B(6)/0.0299792458D0
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
          IF (ABS(PSI).GT.1D-10) THEN
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=-PSI*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='SOLE TILT '
          ENDIF
          CALL MISAL(-1,ENERGY)               !..MISALIG. (INVERSE)
        ELSEIF (INDEX(KEYWORD,'SROT').EQ.1) THEN     !..COORDINATE
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          IA(NEL)=20                                      !  ROTATION
          A(NEL,1)=B(1)
          A(NEL,2)=B(6)*57.295779513082D0
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
        ELSE
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATIN): TOO MANY ELEMENTS'
          WRITE(6,'(3A,I6)')
     .      ' ***MSG(MLATIN): UNSUPPORTED CODE = ',KEYWORD,' AT ',I
          WRITE(6,'(3A,I6)')
     .      ' ***MSG(MLATIN): UNSUPPORTED CODE = ',KEYWORD,' AT ',I
          IA(NEL)=-9999
        ENDIF
        A(NEL,5)=SUML
      ENDDO
 1000 IMAX=NEL                           !..NUMBER OF ELEMENTS
C...............................................................
C......... SUPPRESS ALL SOURCES OF COUPLING.
C......... SOLENOID            ---> DRIFT
C......... ROTATION            ---> MARKER
C......... ROLLED QUAD         ---> UPRIGHT QUAD
C......... HOR MONITOR         ---> MARKER
C......... VER MONITOR         ---> MARKER
C......... H/V MONITOR         ---> MARKER
C......... CENTROID KICK       ---> MARKER
C......... A_N (MULTIPOLE)     ---> ZERO
C.........
      IF (IDIME.EQ.2) THEN
C....................................LOOP OVER THE LATTICE ARRAY
        DO JLAT=1,IMAX
          IF (IA(JLAT).EQ.19) THEN     !..SOLENOID ---> DRIFT
            IA(JLAT)=3
            A(JLAT,2)=0D0
            A(JLAT,3)=0D0
          ELSEIF (IA(JLAT).EQ.20) THEN !..ROTATION ---> MARKER
            IA(JLAT)=0
            A(JLAT,1)=0D0
            A(JLAT,2)=0D0
            A(JLAT,3)=0D0
            A(JLAT,4)=ENERGY
          ELSEIF (IA(JLAT).EQ.55) THEN !..ROLLED QUAD ---> UPRIGHT QUAD
            IA(JLAT)=5
            A(JLAT,3)=1D0
          ELSEIF (IA(JLAT).EQ.101) THEN !..101 CARD ---> MARKER
            IA(JLAT)=0
            A(JLAT,1)=0D0
            A(JLAT,2)=0D0
            A(JLAT,3)=0D0
            A(JLAT,4)=0D0
          ELSEIF (IA(JLAT).EQ.102) THEN !..102 CARD ---> MARKER
            IA(JLAT)=0
            A(JLAT,1)=0D0
            A(JLAT,2)=0D0
            A(JLAT,3)=0D0
            A(JLAT,4)=0D0
          ELSEIF (IA(JLAT).EQ.103) THEN !..103 CARD ---> MARKER
            IA(JLAT)=0
            A(JLAT,1)=0D0
            A(JLAT,2)=0D0
            A(JLAT,3)=0D0
            A(JLAT,4)=0D0
          ELSEIF (IA(JLAT).EQ.107) THEN !..107 CARD ---> MARKER
            IA(JLAT)=0
            A(JLAT,2)=0D0
            A(JLAT,3)=0D0
            A(JLAT,4)=0D0
          ELSEIF (IA(JLAT).EQ.2000) THEN !..MULTIPOLE COMPONENT A_N=0
            JORD=INT(A(JLAT,3))
            IF(JORD.EQ.0) A(JLAT,1)=0D0
            A(JLAT,2)=0D0
          ENDIF
C...............................................................
        ENDDO
C...............................................................
      ENDIF
C...............................................................
      CLOSE(2)
C.........................................................FINISHED
      RETURN
      END
CDECK  ID>, MLATINI.
C=================================================================
C READS A STRUCTURE FILE FROM UNIT 2. IT DOES NOT NEED THE FILE
C NAME. THIS IMPROVES THE PORTABILITY.
C TITLE IS A STRING RED FROM THE MAD FILE.
C IDIME IS THE DIMENSION OF THE PHASE SPACE:
C
C IDIME = 2 ALL THE SOURCES OF COUPLING ARE SUPPRESSED AND ONLY
C           THE HORIZONTAL MOTION IS CONSIDERED.
C
C IDIME = 4 FULL 4D MOTION WITH COUPLING.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MLATINI(TITLE,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      CHARACTER*4 KEYWORD
      CHARACTER*8 PROGVRSN,DATAVRSN,DATE,TIME,JOBNAME
      CHARACTER*10 ALBL,ALBL1
      CHARACTER*16 NAME
      CHARACTER*80 TITLE
      LOGICAL SYMM
      DIMENSION B(28)
      INTRINSIC INDEX
      INTEGER SUPER
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/DISP/B,X,Y,Z,SUML,PHI,THETA,PSI,NEL,ALBL1
      COMMON/FACTOR/FAC(0:11)
C.............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MLATINI): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      READ(2,'(5A8,I8,L8,I8/A80)')   !..READ HEADER
     .  PROGVRSN,DATAVRSN,DATE,TIME,JOBNAME,SUPER,SYMM,NPOS,
     .  TITLE
      WRITE(6,*) PROGVRSN,'...',DATE,'...',TIME
      WRITE(6,*) DATAVRSN
      WRITE(6,*) TITLE
      WRITE(6,*) NPOS,' RECORDS EXPECTED'
C
      IF ((INDEX(DATAVRSN,'STRUCT')+INDEX(DATAVRSN,'STRUCT')).EQ.0)
     .  STOP '***ERROR(MLATINI): WRONG FILE FORMAT'
C..............................................DETERMINE THE ENERGY
        ENERGY=450D0  ! GEV PROTONS
        IF (((INDEX(TITLE,'LHC')+INDEX(TITLE,'lhc')).NE.0.)
     .    .AND.((INDEX(TITLE,'COL')+INDEX(TITLE,'col')).NE.0)) THEN
          ENERGY=8000D0       ! GEV
        ELSEIF (((INDEX(TITLE,'LHC')+INDEX(TITLE,'lhc')).NE.0)
     .    .AND.((INDEX(TITLE,'INJ')+INDEX(TITLE,'inj')).NE.0)) THEN
          ENERGY=450D0        ! GEV
        ELSEIF (((INDEX(TITLE,'LEP')+INDEX(TITLE,'lep')).NE.0)
     .    .AND.((INDEX(TITLE,'COL')+INDEX(TITLE,'col')).NE.0)) THEN
          ENERGY=45.6D0       ! GEV
        ELSEIF (((INDEX(TITLE,'LEP')+INDEX(TITLE,'lep')).NE.0)
     .    .AND.((INDEX(TITLE,'INJ')+INDEX(TITLE,'inj')).NE.0)) THEN
          ENERGY=20D0         ! GEV
        ELSE
 7000     ENERGY=1D0
        ENDIF
C.......................................................OUTPUT FILE
      NEL=0                             !..INITIALIZES COUNTER
C.............................................INITIALIZES ARRAYS
      DO JJ=1,4
        DO II=1,NELE
          A(II,JJ)=0D0
          IA(II)=0
        ENDDO
      ENDDO
C..........................................FACTORIALS
      FAC(0)=1D0
      DO JFACT=1,11
        FAC(JFACT)=FLOAT(JFACT)*FAC(JFACT-1)
      ENDDO
C.................LOOP OVER NPOS ELEMENTS, EACH RECORD 4 LINES LONG
      DO I=1,NPOS
C..............................................READ A SINGLE RECORD
        READ(2,'(A4,A16,F12.6,3E16.9)',END=1000)
     .    KEYWORD,NAME,(B(II),II=1,4)
C......................................................PROCESS DATA
        ALBL1=NAME(1:10)
        III=INDEX(ALBL1,' ')-1
        IF (INDEX(KEYWORD,'DRIF').EQ.1) THEN             !..DRIFT
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
        ELSEIF (INDEX(KEYWORD,'MARK').EQ.1) THEN         !..MARKER
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IF(INDEX(ALBL1,'CHARGEKICK').EQ.1) THEN        !..SPACE CHARGE KICK
            IA(NEL)=2300
          ELSE
            IA(NEL)=0
          ENDIF
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
        ELSEIF (INDEX(KEYWORD,'RBEN').EQ.1) THEN         !..RECT BEND
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          IF (ABS(X).NE.0D0.OR.ABS(Y).NE.0D0.OR.ABS(Z).NE.0D0.OR.
     .    ABS(THETA).NE.0D0.OR.ABS(PHI).NE.0D0)
     .    WRITE(6,*) '***ERROR(MLATINI): DISPLACEMENTS AT ',ALBL1,
     .    ' IGNORED'
C         IF (ABS(Z).GT.1D-10) THEN              !..LONGITUDINAL DISP.
C           NEL=NEL+1                           !..NOT IMPLEMENTED YET
C         IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
C           IA(NEL)=3
C           A(NEL,1)=-Z
C           A(NEL,4)=ENERGY
C           ALBL(NEL)=ALBL1
C         ENDIF
          IF (ABS(B(5)+PSI).GT.1D-10) THEN       !..TILT ANGLE
            NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,1)=(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='RBEND TILT '
          ENDIF
          IF (ABS(B(8)).GT.1D-10)                               !..H1
     .      WRITE(6,*) '***ERROR(MLATINI): RBEND  H1  AT ',ALBL1
          IF (ABS(B(9)).GT.1D-10)                               !..H2
     .      WRITE(6,*) '***ERROR(MLATINI): RBEND  H2  AT ',ALBL1
          IF (ABS(B(6)).GT.1D-10) THEN                          !..E1
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=2
            A(NEL,2)=B(6)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='WEDGE      '
          ENDIF
          IF (ABS(B(3)).GT.1D-10) THEN
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=2                           !..WEDGE OF THE RECT
            A(NEL,1)=0D0
            A(NEL,2)=0D0
            A(NEL,3)=-1D0
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=4
            A(NEL,1)=B(1)
            A(NEL,2)=ENERGY*B(2)/(0.0299792458D0*B(1))
            A(NEL,3)=B(3)*ENERGY/0.0299792458D0
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=2                           !..WEDGE OF THE RECT
            A(NEL,1)=0D0
            A(NEL,2)=0D0
            A(NEL,3)=-1D0
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
          ELSE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=44
            A(NEL,1)=B(1)
            A(NEL,2)=ENERGY*B(2)/(0.0299792458D0*B(1))
            A(NEL,3)=0D0
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
          ENDIF
          IF (ABS(B(7)).GT.1D-10) THEN                      !..E2
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=2
            A(NEL,2)=B(7)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='WEDGE      '
          ENDIF
          IF (ABS(B(5)+PSI).GT.1D-10) THEN              !..TILT ANGLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=-(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='RBEND TILT '
          ENDIF
C         IF (ABS(Z).GT.1D-10) THEN              !..LONGITUDINAL DISP.
C           NEL=NEL+1                           !..NOT IMPLMENTED YET
C           IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
C           IA(NEL)=3
C           A(NEL,1)=Z
C           A(NEL,4)=ENERGY
C           ALBL(NEL)=ALBL1
C         ENDIF
        ELSEIF (INDEX(KEYWORD,'SBEN').EQ.1) THEN        !..SECTOR BEND
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          IF (ABS(X).NE.0D0.OR.ABS(Y).NE.0D0.OR.ABS(Z).NE.0D0.OR.
     .    ABS(THETA).NE.0D0.OR.ABS(PHI).NE.0D0)
     .    WRITE(6,*) '***ERROR(MLATINI): DISPLACEMENTS AT ',ALBL1,
     .    ' IGNORED'
C         IF (ABS(Z).GT.1D-10) THEN              !..LONGITUDINAL DISP.
C           NEL=NEL+1                           !..NOT IMPLEMENTED YET
C           IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
C           IA(NEL)=3
C           A(NEL,1)=-Z
C           A(NEL,4)=ENERGY
C           ALBL(NEL)=ALBL1
C         ENDIF
          IF (ABS(B(5)+PSI).GT.1D-10) THEN              !..TILT ANGLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='SBEND TILT '
          ENDIF
          IF (ABS(B(8)).GT.1D-10)                                !..H1
     .      WRITE(6,*) '***ERROR(MLATINI): SBEND  H1  AT ',ALBL1
          IF (ABS(B(9)).GT.1D-10)                                !..H2
     .      WRITE(6,*) '***ERROR(MLATINI): SBEND  H2  AT ',ALBL1
          IF (ABS(B(6)).GT.1D-10) THEN                           !..E1
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=2
            A(NEL,2)=B(6)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='WEDGE      '
          ENDIF
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=4
          A(NEL,1)=B(1)
          A(NEL,2)=ENERGY*B(2)/(0.0299792458D0*B(1))
          A(NEL,3)=B(3)*ENERGY/0.0299792458D0
          A(NEL,4)=ENERGY
          IF (ABS(B(7)).GT.1D-10) THEN                           !..E2
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=2
            A(NEL,2)=B(7)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='WEDGE      '
          ENDIF
          IF (ABS(B(5)).GT.1D-10) THEN                  !..TILT ANGLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=-(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='SBEND TILT '
          ENDIF
C         IF (ABS(Z).GT.1D-10) THEN              !..LONGITUDINAL DISP.
C           NEL=NEL+1                           !..NOT IMPLEMENTED YET
C           IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
C           IA(NEL)=3
C           A(NEL,1)=Z
C           A(NEL,4)=ENERGY
C           ALBL(NEL)=ALBL1
C         ENDIF
        ELSEIF (INDEX(KEYWORD,'QUAD').EQ.1) THEN       !..QUADRUPOLE
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          CALL MISAL(1,ENERGY)                !..MISALIG.
          IF (ABS(B(5)+PSI).GT.1D-10) THEN       !..ROLLED QUAD
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=55
            A(NEL,3)=(B(5)+PSI)*57.295779513082D0
          ELSE                               !..UPRIGHT QUAD
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=5
            A(NEL,3)=1D0
          ENDIF
          A(NEL,1)=B(1)
          A(NEL,2)=ENERGY*B(3)/0.0299792458D0
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
          CALL MISAL(-1,ENERGY)               !..MISALIG. (INVERSE)
        ELSEIF (INDEX(KEYWORD,'SEXT').EQ.1) THEN      !..SEXTUPOLE
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          CALL MISAL(1,ENERGY)                !..MISALIG.
          IF (ABS(B(5)+PSI).GT.1D-10) THEN              !..TILT ANGLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='SEXT TILT '
          ENDIF
          IF (ABS(B(1)).GE.1D-10) THEN   !..THICK SEXTUPOLE (DR+SEX+DR)
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=5D-1*B(1)
            A(NEL,4)=ENERGY
            ALBL(NEL)='SEXT DRIFT'
          ENDIF
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=2002
          IF (B(4).EQ.0.) IA(NEL)=-2002      !..FLAG FOR LINEAR BLOCS
          A(NEL,1)=B(4)*B(1)/FAC(2)          !..K_2
          A(NEL,2)=0D0
          A(NEL,3)=2D0
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
          IF (ABS(B(1)).GE.1D-10) THEN   !..THICK SEXTUPOLE (DR+SEX+DR)
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=5D-1*B(1)
            A(NEL,4)=ENERGY
            ALBL(NEL)='SEXT DRIFT'
          ENDIF
          IF (ABS(B(5)+PSI).GT.1D-10) THEN      !..TILT SEXTUPOLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=-(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='SEXT TILT '
          ENDIF
          CALL MISAL(-1,ENERGY)               !..MISALIG. (INVERSE)
        ELSEIF (INDEX(KEYWORD,'OCTU').EQ.1) THEN     !..OCTUPOLE
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          CALL MISAL(1,ENERGY)                !..MISALIG.
          IF (ABS(B(5)+PSI).GT.1D-10) THEN      !..TILT OCTUPOLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='OCTUP TILT'
          ENDIF
          IF (ABS(B(1)).GE.1D-10) THEN   !..THICK OCTUPOLE (DR+OCT+DR)
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=5D-1*B(1)
            ALBL(NEL)=ALBL1
          ENDIF
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=2003
          IF (B(6).EQ.0.) IA(NEL)=-2003      !..FLAG FOR LINEAR BLOCS
          A(NEL,1)=B(6)*B(1)/FAC(3)          !..K_3
          A(NEL,2)=0D0
          A(NEL,3)=3D0
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
          IF (ABS(B(1)).GE.1D-10) THEN   !..THICK OCTUPOLE (DR+OCT+DR)
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=5D-1*B(1)
            ALBL(NEL)='OCTU DRIFT'
          ENDIF
          IF (ABS(B(5)+PSI).GT.1D-10) THEN      !..TILT OCTUPOLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=-(B(5)+PSI)*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='OCTUP TILT'
          ENDIF
          CALL MISAL(-1,ENERGY)               !..MISALIG. (INVERSE)
        ELSEIF (INDEX(KEYWORD,'MULT').EQ.1) THEN      !..MULTIPOLE
          BACKSPACE(2)
          READ(2,'(A4,A16,F12.6,2E16.9,I16)')
     .    KEYWORD,NAME,(B(II),II=1,3),NMULT
          READ(2,'(5E16.9)') (B(II),II=5,5+NMULT) !..K_N
          READ(2,'(5E16.9)') (B(II),II=6+NMULT,6+2*NMULT) !..J_N
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          CALL MISAL(1,ENERGY)                !..MISALIG.
          IF (ABS(PSI).GT.1D-10) THEN      !..TILT MULTIPOLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=PSI*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='MULTI TILT'
          ENDIF
          IF (ABS(B(2)+B(5)+B(3)+B(6+NMULT)).NE.0D0) THEN
            NEL=NEL+1                                !..DIPOLE KICK
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=107
            A(NEL,1)=B(2)+B(5)
            A(NEL,2)=0D0
            A(NEL,3)=B(3)+B(6+NMULT)
            A(NEL,4)=0D0
            ALBL(NEL)=ALBL1
          ENDIF
          NEL=NEL+1                                !..QUADRUPOLE KICK
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=-1000
          IF (ABS(B(6)+B(NMULT+7)).NE.0D0) IA(NEL)=1000
          A(NEL,1)=B(6)
          A(NEL,2)=B(NMULT+7)
          A(NEL,3)=1
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
          DO JMULT=0,1                   !..MULTIPOLE COMPONENTS
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=-2000
            A(NEL,1)=0D0
            A(NEL,2)=0D0
            A(NEL,3)=JMULT
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
          ENDDO
          DO JMULT=2,NMULT                   !..MULTIPOLE COMPONENTS
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=2000
            A(NEL,1)=B(JMULT+5)/FAC(JMULT)
            A(NEL,2)=B(NMULT+JMULT+6)/FAC(JMULT)
            A(NEL,3)=JMULT
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
            IF ((ABS(A(NEL,1))+ABS(A(NEL,2))).EQ.0D0) IA(NEL)=-2000
          ENDDO
          IF (ABS(PSI).GT.1D-10) THEN      !..TILT MULTIPOLE
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=-PSI*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='MULTI TILT'
          ENDIF
          CALL MISAL(-1,ENERGY)               !..MISALIG. (INVERSE)
        ELSEIF (INDEX(KEYWORD,'MONI').EQ.1.) THEN   ! MONITOR
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=5D-1*B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)='DRIFT_MON'
C............................................................
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=103                !..SPECIAL CARD FOR V.Z.
          ALBL(NEL)=ALBL1
C............................................................
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=5D-1*B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)='DRIFT_MON'
        ELSEIF (INDEX(KEYWORD,'HMON').EQ.1) THEN          !..BPM
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=5D-1*B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)='DRIFT_MON'
C............................................................
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=101                !..SPECIAL CARD FOR V.Z.
          ALBL(NEL)=ALBL1
C............................................................
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=5D-1*B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)='DRIFT_MON'
        ELSEIF (INDEX(KEYWORD,'VMON').EQ.1) THEN          !..BPM
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=5D-1*B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)='DRIFT_MON'
C............................................................
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=102                !..SPECIAL CARD FOR V.Z.
          ALBL(NEL)=ALBL1
C............................................................
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=5D-1*B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)='DRIFT_MON'
        ELSEIF (INDEX(KEYWORD,'KICK').EQ.1.OR.      !..ORBIT CORRECTOR
     .          INDEX(KEYWORD,'HKIC').EQ.1.OR.
     .          INDEX(KEYWORD,'VKIC').EQ.1) THEN
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          IF (ABS(B(1)).GT.1D-10) THEN
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=5D-1*B(1)
            A(NEL,4)=ENERGY
            ALBL(NEL)='DRIFT_COR'
          ENDIF
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=107
          ALBL(NEL)=ALBL1
          IF (ABS(B(1)).GT.1D-10) THEN
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=5D-1*B(1)
            A(NEL,4)=ENERGY
            ALBL(NEL)='DRIFT_COR'
          ENDIF
        ELSEIF (INDEX(KEYWORD,'RFCA').EQ.1) THEN       !..RF-CAVITY
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
C         CALL MISAL(1,ENERGY)                !..MISALIG.
C         IF (ABS(PSI).GT.1D-10) THEN     !..IGNORED FOR THE MOMENT
C           NEL=NEL+1
C           IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
C           IA(NEL)=20
C           A(NEL,2)=PSI*57.295779513082D0
C           A(NEL,4)=ENERGY
C           ALBL(NEL)='CAVIT TILT'
C         ENDIF
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=31
          A(NEL,1)=B(1)                  !..LENGTH
          A(NEL,2)=B(7)                  !..VOLTAGE IN MV
          A(NEL,3)=B(8)*360D0            !..LAG IN DEGREE
          A(NEL,4)=B(6)                  !..FREQUENCY IN MHZ
          ALBL(NEL)=ALBL1
C         IF (ABS(PSI).GT.1D-10) THEN     !..IGNORED FOR THE MOMENT
C           NEL=NEL+1
C           IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
C           IA(NEL)=20
C           A(NEL,2)=PSI*57.295779513082D0
C           A(NEL,4)=ENERGY
C           ALBL(NEL)='CAVIT TILT'
C         ENDIF
C         CALL MISAL(-1,ENERGY)               !..MISALIG. (INVERSE)
C         WRITE(2,'(I5,F12.7,G12.5,G13.6,G18.10,1X,A1,A10,A1)')
C    .      IA,(A(II),II=1,4),'''',ALBL,''''
 
        ELSEIF (INDEX(KEYWORD,'RCOL').EQ.1) THEN    !..COLLIMATOR
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=2100  !  RECTANGULAR COLLIMATOR
          A(NEL,1)=B(2)
          A(NEL,2)=B(3)
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
          IF (B(2).LT.1D-20.OR.B(3).LT.1D-20) IA(NEL)=-2100
          IF (B(1).GT.0D0) THEN   !  ADD A DRIFT AND A COLLIMATOR
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=B(1)
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=2100  !  RECTANGULAR COLLIMATOR
            A(NEL,1)=B(2)
            A(NEL,2)=B(3)
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
            IF (B(2).LT.1D-20.OR.B(3).LT.1D-20) IA(NEL)=-2100
          ENDIF
        ELSEIF (INDEX(KEYWORD,'ECOL').EQ.1) THEN    !..COLLIMATOR
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=2200  !  ELLIPTICAL COLLIMATOR
          A(NEL,1)=B(2)*B(2)
          A(NEL,2)=B(3)*B(3)
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
          IF (B(2).LT.1D-20.OR.B(3).LT.1D-20) IA(NEL)=-2200
          IF (B(1).GT.0D0) THEN   !  ADD A DRIFT AND A COLLIMATOR
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=3
            A(NEL,1)=B(1)
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=2200  !  ELLIPTICAL COLLIMATOR
            A(NEL,1)=B(2)*B(2)
            A(NEL,2)=B(3)*B(3)
            A(NEL,4)=ENERGY
            ALBL(NEL)=ALBL1
            IF (B(2).LT.1D-20.OR.B(3).LT.1D-20) IA(NEL)=-2200
          ENDIF
        ELSEIF (INDEX(KEYWORD,'ELSE').EQ.1) THEN  !..E-STATIC SEPARATOR
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
        ELSEIF (INDEX(KEYWORD,'    ').EQ.1) THEN
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=B(1)
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
        ELSEIF (INDEX(KEYWORD,'SOLE').EQ.1) THEN     !..SOLENOID
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          CALL MISAL(1,ENERGY)                !..MISALIG.
          IF (ABS(PSI).GT.1D-10) THEN
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=PSI*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='SOLE TILT '
          ENDIF
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=19
          A(NEL,1)=B(1)
          A(NEL,2)=ENERGY*B(6)/0.0299792458D0
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
          IF (ABS(PSI).GT.1D-10) THEN
            NEL=NEL+1
            IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
            IA(NEL)=20
            A(NEL,2)=-PSI*57.295779513082D0
            A(NEL,4)=ENERGY
            ALBL(NEL)='SOLE TILT '
          ENDIF
          CALL MISAL(-1,ENERGY)               !..MISALIG. (INVERSE)
        ELSEIF (INDEX(KEYWORD,'SROT').EQ.1) THEN     !..COORDINATE
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          IA(NEL)=20                                      !  ROTATION
          A(NEL,1)=B(1)
          A(NEL,2)=B(6)*57.295779513082D0
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
        ELSE
          READ(2,'(5E16.9)') (B(II),II=5,9)
          READ(2,'(4E16.9/3E16.9)') X,Y,Z,SUML,PHI,THETA,PSI
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MLATINI): TOO MANY ELEMENTS'
          WRITE(6,'(3A,I6)')
     .      ' ***MSG(MLATINI): UNSUPPORTED CODE = ',KEYWORD,' AT ',I
          WRITE(6,'(3A,I6)')
     .      ' ***MSG(MLATINI): UNSUPPORTED CODE = ',KEYWORD,' AT ',I
          IA(NEL)=-9999
        ENDIF
      ENDDO
 1000 IMAX=NEL                           !..NUMBER OF ELEMENTS
C...............................................................
C......... SUPPRESS ALL SOURCES OF COUPLING.
C......... SOLENOID            ---> DRIFT
C......... ROTATION            ---> MARKER
C......... ROLLED QUAD         ---> UPRIGHT QUAD
C......... HOR MONITOR         ---> MARKER
C......... VER MONITOR         ---> MARKER
C......... H/V MONITOR         ---> MARKER
C......... CENTROID KICK       ---> MARKER
C......... A_N (MULTIPOLE)     ---> ZERO
C.........
      IF (IDIME.EQ.2) THEN
C....................................LOOP OVER THE LATTICE ARRAY
        DO JLAT=1,IMAX
          IF (IA(JLAT).EQ.19) THEN     !..SOLENOID ---> DRIFT
            IA(JLAT)=3
            A(JLAT,2)=0D0
            A(JLAT,3)=0D0
          ELSEIF (IA(JLAT).EQ.20) THEN !..ROTATION ---> MARKER
            IA(JLAT)=0
            A(JLAT,1)=0D0
            A(JLAT,2)=0D0
            A(JLAT,3)=0D0
            A(JLAT,4)=ENERGY
          ELSEIF (IA(JLAT).EQ.55) THEN !..ROLLED QUAD ---> UPRIGHT QUAD
            IA(JLAT)=5
            A(JLAT,3)=1D0
          ELSEIF (IA(JLAT).EQ.101) THEN !..101 CARD ---> MARKER
            IA(JLAT)=0
            A(JLAT,1)=0D0
            A(JLAT,2)=0D0
            A(JLAT,3)=0D0
            A(JLAT,4)=0D0
          ELSEIF (IA(JLAT).EQ.102) THEN !..102 CARD ---> MARKER
            IA(JLAT)=0
            A(JLAT,1)=0D0
            A(JLAT,2)=0D0
            A(JLAT,3)=0D0
            A(JLAT,4)=0D0
          ELSEIF (IA(JLAT).EQ.103) THEN !..103 CARD ---> MARKER
            IA(JLAT)=0
            A(JLAT,1)=0D0
            A(JLAT,2)=0D0
            A(JLAT,3)=0D0
            A(JLAT,4)=0D0
          ELSEIF (IA(JLAT).EQ.107) THEN !..107 CARD ---> MARKER
            IA(JLAT)=0
            A(JLAT,2)=0D0
            A(JLAT,3)=0D0
            A(JLAT,4)=0D0
          ELSEIF (IA(JLAT).EQ.2000) THEN !..MULTIPOLE COMPONENT A_N=0
            JORD=INT(A(JLAT,3))
            IF(JORD.EQ.0) A(JLAT,1)=0D0
            A(JLAT,2)=0D0
          ENDIF
C...............................................................
        ENDDO
C...............................................................
      ENDIF
C...............................................................
      CLOSE(2)
C.........................................................FINISHED
      RETURN
      END
CDECK  ID>, MISAL.
C=================================================================
C GENERATES THE CARDS FOR MISALIGMENT (BEAM CENTROID KICK).
C INV INDICATES WHETHER THE DIRECT OR INVERSE TRANSFORMATION IS
C NEEDED.
C ENERGY IS THE ENERGY OF THE MACHINE
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MISAL(INV,ENERGY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      DIMENSION B(28)
      CHARACTER*10 ALBL,ALBL1
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/DISP/B,X,Y,Z,SUML,PHI,THETA,PSI,NEL,ALBL1
C...............................................................
      IF (INV.NE.1.AND.INV.NE.-1) THEN
        WRITE(6,*) '***ERROR(MISAL): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (ENERGY.LT.0D0) THEN
        WRITE(6,*) '***ERROR(MISAL): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      IF (INV.EQ.1) THEN
        IF (ABS(Z).GT.1D-10) THEN                   !..DRIFT
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MISAL): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=-Z
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
        ENDIF
        IF ((ABS(X)+ABS(Y)+ABS(THETA)+ABS(PHI)).GT.1D-10) THEN
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MISAL): TOO MANY ELEMENTS'
          IA(NEL)=107                                   !..107 CARD
          A(NEL,1)=X+5D-1*B(1)*THETA
          A(NEL,2)=THETA
          A(NEL,3)=Y+5D-1*B(1)*PHI
          A(NEL,4)=PHI
          ALBL(NEL)=ALBL1
        ENDIF
      ELSEIF (INV.EQ.-1) THEN
        IF ((ABS(X)+ABS(Y)+ABS(THETA)+ABS(PHI)).GT.1D-10) THEN
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MISAL): TOO MANY ELEMENTS'
          IA(NEL)=107                                   !..107 CARD
          A(NEL,1)=-X-5D-1*B(1)*THETA
          A(NEL,2)=-THETA
          A(NEL,3)=-Y-5D-1*B(1)*PHI
          A(NEL,4)=-PHI
          ALBL(NEL)=ALBL1
        ENDIF
        IF (ABS(Z).GT.1D-10) THEN                   !..DRIFT
          NEL=NEL+1
          IF (NEL.GT.NELE) STOP '***ERROR(MISAL): TOO MANY ELEMENTS'
          IA(NEL)=3
          A(NEL,1)=Z
          A(NEL,4)=ENERGY
          ALBL(NEL)=ALBL1
        ENDIF
C................................................................
      ENDIF
C...............................................................
      RETURN
      END
CDECK  ID>, MRESCALE.
C=================================================================
C RESCALES THE VARIABLES TO MM, MRAD.
C THE NON-LINEAR ELEMENT STRENGTHS ARE RESCALED ACCORDING TO
C K_N(NEW)=10^(3N-3) * K_N(OLD).
C THE INHOMOGENEITIES ARE RESCALED AS WELL
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MRESCALE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      CHARACTER*10 ALBL
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
C...................................................
C....................................LOOP OVER THE LUMPS
      DO JLUMP=1,NLMAX
        DO JDIME=1,4              !..RESCALES INHOMOGENEITIES
          RK(JDIME,JLUMP)=1D3*RK(JDIME,JLUMP)
        ENDDO
        DO JNL=INL(1,JLUMP),INL(2,JLUMP)
          IF (IA(JNL).EQ.2000) THEN       !..MULTIPOLE ELEMENT
            IMULT=INT(A(JNL,3))
            FACT=1D1**(3*(1-IMULT))
            A(JNL,1)=A(JNL,1)*FACT
            A(JNL,2)=A(JNL,2)*FACT
          ELSEIF (IA(JNL).EQ.2002) THEN   !..SEXTUPOLE
            A(JNL,1)=A(JNL,1)*1D-3
          ELSEIF (IA(JNL).EQ.2003) THEN   !..OCTUPOLE
            A(JNL,1)=A(JNL,1)*1D-6
          ENDIF
C...............................................................
        ENDDO
C...............................................................
      ENDDO
C...............................................................
      DO JDIME=1,4              !..RESCALES INHOMOGENEITIES
        RK(JDIME,NLMAX+1)=1D3*RK(JDIME,NLMAX+1) !..LAST BLOCK
      ENDDO
C...............................................................
      RETURN
      END
CDECK  ID>, NLFUN.
C=================================================================
C COMPUTES A NON-LINEAR KICK IN THE CASE OF PLAIN TRACKING.
C X2 IS THE INITIAL POINT
C X1 IS THE FINAL POINT
C IL IS THE LINEAR BLOCK NUMBER
C IDIME IS THE DIMENSION OF THE PHASE SPACE
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE NLFUN(X2,X1,IL,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      DIMENSION X1(4),X2(4)
      CHARACTER*10 ALBL
      COMPLEX*16 Z1,Z2
      COMMON/CSTRANS/CSM(4,4),CSMIN(4,4)
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
C...............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(NLFUN): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      IF (IDIME.EQ.2) THEN
        X2(3)=0D0
        X2(4)=0D0
      ENDIF
C...............................................................
      MINORD=INT(A(INL(1,IL),3))
      MAXORD=INT(A(INL(2,IL),3))
C...............................................................
      Z1=DCMPLX(X2(1),X2(3))        !..DEFINES COMPLEX VARIABLES
      Z2=DCMPLX(A(INL(2,IL),1),-A(INL(2,IL),2))  !..AGAIN COMPLEX
C...............................................................
      IF (IA(INL(1,IL)).EQ.2002) THEN
        Z2=Z2*Z1*Z1                             !..SEXTUPOLE KICK
      ELSEIF (IA(INL(1,IL)).EQ.2003) THEN
        Z2=Z2*Z1*Z1*Z1                          !..OCTUPOLE KICK
      ELSEIF (IA(INL(1,IL)).EQ.2000) THEN
        DO JCOMP=INL(2,IL)-1,INL(1,IL),-1       !..MULTIPOLE
          Z2=Z2*Z1+DCMPLX(A(JCOMP,1),-A(JCOMP,2))
        ENDDO
C...............................................................
        IF (MINORD.NE.0) THEN
C...............................................................
          JEL=INL(1,IL)-1
 10       IF (IA(JEL).EQ.-2000) THEN
            Z2=Z2*Z1+DCMPLX(A(JEL,1),-A(JEL,2))
          ELSE
            GOTO 20
          ENDIF
          JEL=JEL-1
          GOTO 10
C...............................................................
        ENDIF
C...............................................................
      ELSEIF (IA(INL(1,IL)).EQ.2100) THEN  !  RECTANGULAR COLLIMATOR
        IF (IDIME.EQ.2) THEN
          IF (ABS(X2(1)).GT.A(INL(1,IL),1)*CSM(1,1)) THEN
            WRITE(6,*) '***MESSAGE(NLFUN): LOST AT COLLIMATOR ',
     .                         ALBL(INL(1,IL))
            RETURN
          ENDIF
        ELSE
          IF (ABS(X2(1)).GT.A(INL(1,IL),1)*CSM(1,1).OR.
     .        ABS(X2(3)).GT.A(INL(1,IL),2)*CSM(3,3)) THEN
            WRITE(6,*) '***MESSAGE(NLFUN): LOST AT COLLIMATOR ',
     .                         ALBL(INL(1,IL))
            RETURN
          ENDIF
        ENDIF
      ELSEIF (IA(INL(1,IL)).EQ.2200) THEN  !  ELLIPTICAL COLLIMATOR
        IF (IDIME.EQ.2) THEN
          IF (X2(1)*X2(1)/(A(INL(1,IL),1)*CSM(1,1)*CSM(1,1)).GT.1D0)
     .       THEN
            WRITE(6,*) '***MESSAGE(NLFUN): LOST AT COLLIMATOR ',
     .                         ALBL(INL(1,IL))
            RETURN
          ENDIF
        ELSE
          IF (X2(1)*X2(1)/(A(INL(1,IL),1)*CSM(1,1)*CSM(1,1))+
     .        X2(3)*X2(3)/(A(INL(1,IL),2)*CSM(3,3)*CSM(3,3)).GT.1D0)
     .       THEN
            WRITE(6,*) '***MESSAGE(NLFUN): LOST AT COLLIMATOR ',
     .                         ALBL(INL(1,IL))
            RETURN
          ENDIF
        ENDIF
C...............................................................
      ELSEIF (IA(INL(1,IL)).EQ.2300) THEN  !  SPACE CHARGE KICK
        WRITE(*,*) 'FOUND A SPACE CHARGE KICK!'
      ENDIF
C...............................................................
 20   X1(1)=X2(1)         !..NEW POINT
      X1(2)=X2(2)-DREAL(Z2)
      X1(3)=X2(3)
      X1(4)=X2(4)+DIMAG(Z2)
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, NLFUNI.
C=================================================================
C COMPUTES THE INVERSE OF A NON-LINEAR KICK IN THE CASE OF PLAIN
C TRACKING.
C X2 IS THE INITIAL POINT
C X1 IS THE FINAL POINT
C IL IS THE LINEAR BLOCK NUMBER
C IDIME IS THE DIMENSION OF THE PHASE SPACE
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE NLFUNI(X2,X1,IL,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      DIMENSION X1(4),X2(4)
      CHARACTER*10 ALBL
      COMPLEX*16 Z1,Z2
      COMMON/CSTRANS/CSM(4,4),CSMIN(4,4)
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
C...............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(NLFUNI): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      IF (IDIME.EQ.2) THEN
        X2(3)=0D0
        X2(4)=0D0
      ENDIF
C...............................................................
      MINORD=INT(A(INL(1,IL),3))
      MAXORD=INT(A(INL(2,IL),3))
C...............................................................
      Z1=DCMPLX(X2(1),X2(3))        !..DEFINES COMPLEX VARIABLES
      Z2=DCMPLX(-A(INL(2,IL),1),A(INL(2,IL),2)) !..AGAIN COMPLEX
C...............................................................
      IF (IA(INL(1,IL)).EQ.2002) THEN
        Z2=Z2*Z1*Z1                             !..SEXTUPOLE KICK
      ELSEIF (IA(INL(1,IL)).EQ.2003) THEN
        Z2=Z2*Z1*Z1*Z1                          !..OCTUPOLE KICK
      ELSEIF (IA(INL(1,IL)).EQ.2000) THEN
        DO JCOMP=INL(2,IL)-1,INL(1,IL),-1       !..MULTIPOLE
          Z2=Z2*Z1+DCMPLX(-A(JCOMP,1),A(JCOMP,2))
        ENDDO
C...............................................................
        IF (MINORD.NE.0) THEN
C...............................................................
          JEL=INL(1,IL)-1
 10       IF (IA(JEL).EQ.-2000) THEN
            Z2=Z2*Z1+DCMPLX(-A(JEL,1),A(JEL,2))
          ELSE
            GOTO 20
          ENDIF
          JEL=JEL-1
          GOTO 10
C...............................................................
        ENDIF
C...............................................................
      ELSEIF (IA(INL(1,IL)).EQ.2100) THEN  !  RECTANGULAR COLLIMATOR
        IF (IDIME.EQ.2) THEN
          IF (ABS(X2(1)).GT.A(INL(1,IL),1)*CSM(1,1)) THEN
            WRITE(6,*) '***MESSAGE(NLFUNI): LOST AT COLLIMATOR ',
     .                         ALBL(INL(1,IL))
            RETURN
          ENDIF
        ELSE
          IF (ABS(X2(1)).GT.A(INL(1,IL),1)*CSM(1,1).OR.
     .        ABS(X2(3)).GT.A(INL(1,IL),2)*CSM(3,3)) THEN
            WRITE(6,*) '***MESSAGE(NLFUNI): LOST AT COLLIMATOR ',
     .                         ALBL(INL(1,IL))
            RETURN
          ENDIF
        ENDIF
      ELSEIF (IA(INL(1,IL)).EQ.2200) THEN  !  ELLIPTICAL COLLIMATOR
        IF (IDIME.EQ.2) THEN
          IF (X2(1)*X2(1)/(A(INL(1,IL),1)*CSM(1,1)*CSM(1,1)).GT.1D0)
     .       THEN
            WRITE(6,*) '***MESSAGE(NLFUNI): LOST AT COLLIMATOR ',
     .                         ALBL(INL(1,IL))
            RETURN
          ENDIF
        ELSE
          IF (X2(1)*X2(1)/(A(INL(1,IL),1)*CSM(1,1)*CSM(1,1))+
     .        X2(3)*X2(3)/(A(INL(1,IL),2)*CSM(3,3)*CSM(3,3)).GT.1D0)
     .       THEN
            WRITE(6,*) '***MESSAGE(NLFUNI): LOST AT COLLIMATOR ',
     .                         ALBL(INL(1,IL))
            RETURN
          ENDIF
        ENDIF
C...............................................................
      ELSEIF (IA(INL(1,IL)).EQ.2300) THEN  !  SPACE CHARGE KICK
        WRITE(*,*) 'FOUND A SPACE CHARGE KICK!'
        WRITE(*,*) 'IT IS NOT POSSIBLE TO INVERT SUCH A KICK'
        WRITE(*,*) 'PROGRAM STOPS'
        STOP
      ENDIF
C...............................................................
 20   X1(1)=X2(1)         !..NEW POINT
      X1(2)=X2(2)-DREAL(Z2)
      X1(3)=X2(3)
      X1(4)=X2(4)+DIMAG(Z2)
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MCONVMAT.
C=================================================================
C CONVERTS A MATRIX INTO A LINEAR POLYNOMIAL TO COMPUTE THE TRANSFER
C MAP USING OUR COMPOSITION PROGRAM.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C IL IS THE LINEAR BLOCK NUMBER.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MCONVMAT(IDIME,IL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      CHARACTER*10 ALBL
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
C...............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MCONVMAT): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      JORD=1               !..INITIALIZATION (LINEAR STUFF)
C...............................................................
      DO ICOMP=1,IDIME
        RMAT(IPUNT(0,0,0,0),ICOMP)=RK(ICOMP,IL)
C...............................................................
        ICOUNT=0           !..COUNTER FOR VARIABLES
        DO I1=0,JORD
C...............................................................
          DO I2=0,JORD-I1
C...............................................................
            DO I3=0,JORD-I1-I2
C...............................................................
              I4=JORD-I1-I2-I3
              ICOUNT=ICOUNT+1
              RMAT(IPUNT(I4,I3,I2,I1),ICOMP)=RL(ICOMP,ICOUNT,IL)
C...............................................................
            ENDDO
C...............................................................
          ENDDO
C...............................................................
        ENDDO
C...............................................................
      ENDDO
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MZEROC.
C=================================================================
C SETS A COMPLEX POLYNOMIAL TO ZERO.
C A IS THE POLYNOMIAL TO BE SET TO ZERO.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MZEROC(A,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      COMPLEX*16 A(NDIM,4)
C...............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MZEROC): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      DO ICOMP=1,IDIME
C...............................................................
        DO JEL=1,NDIM
          A(JEL,ICOMP)=DCMPLX(0D0,0D0)
C...............................................................
        ENDDO
C...............................................................
      ENDDO
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MZERO.
C=================================================================
C SETS A REAL POLYNOMIAL TO ZERO.
C A IS THE POLYNOMIAL TO BE SET TO ZERO.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MZERO(A,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION A(NDIM,4)
C...............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MZERO): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      DO ICOMP=1,IDIME
C...............................................................
        DO JEL=1,NDIM
          A(JEL,ICOMP)=0D0
C...............................................................
        ENDDO
C...............................................................
      ENDDO
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MONEC.
C=================================================================
C SETS A COMPLEX POLYNOMIAL TO IDENTITY.
C A IS THE POLYNOMIAL TO BE SET TO ONE.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MONEC(A,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      COMPLEX*16 A(NDIM,4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C...............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MONEC): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      CALL MZEROC(A,IDIME)
C...............................................................
      JORD=1
      ICOUNT=0
C...............................................................
      DO I1=0,JORD
C...............................................................
        DO I2=0,JORD-I1
C...............................................................
          DO I3=0,JORD-I1-I2
C...............................................................
            I4=JORD-I1-I2-I3
            ICOUNT=ICOUNT+1
            A(IPUNT(I4,I3,I2,I1),ICOUNT)=DCMPLX(1D0,0D0)
C...............................................................
          ENDDO
C...............................................................
        ENDDO
C...............................................................
      ENDDO
C...............................................................
      IF (IDIME.EQ.2) THEN
        A(IPUNT(0,0,1,0),3)=DCMPLX(0D0,0D0)
        A(IPUNT(0,0,0,1),4)=DCMPLX(0D0,0D0)
      ENDIF
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MONE.
C=================================================================
C SETS A REAL POLYNOMIAL TO IDENTITY.
C A IS THE POLYNOMIAL TO BE SET TO ONE.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MONE(A,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION A(NDIM,4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C...............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MONE): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      CALL MZERO(A,IDIME)
C...............................................................
      JORD=1
      ICOUNT=0
C...............................................................
      DO I1=0,JORD
C...............................................................
        DO I2=0,JORD-I1
C...............................................................
          DO I3=0,JORD-I1-I2
C...............................................................
            I4=JORD-I1-I2-I3
            ICOUNT=ICOUNT+1
            A(IPUNT(I4,I3,I2,I1),ICOUNT)=1D0
C...............................................................
          ENDDO
C...............................................................
        ENDDO
C...............................................................
      ENDDO
C...............................................................
      IF (IDIME.EQ.2) THEN
        A(IPUNT(0,0,1,0),3)=0D0
        A(IPUNT(0,0,0,1),4)=0D0
      ENDIF
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MSWAPC.
C=================================================================
C PERFORMS THE FOLLOWING OPERATION ON TWO COMPLEX POLYNOMIALS:
C A=B & B=0
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MSWAPC(A,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      COMPLEX*16 A(NDIM,4),B(NDIM,4)
C...............................................................
      DO ICOMP=1,4
C...............................................................
        DO JDIM=1,NDIM
C...............................................................
          A(JDIM,ICOMP)=B(JDIM,ICOMP)
          B(JDIM,ICOMP)=DCMPLX(0D0,0D0)
C...............................................................
        ENDDO
C...............................................................
      ENDDO
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MSWAP.
C=================================================================
C PERFORMS THE FOLLOWING OPERATION ON TWO COMPLEX POLYNOMIALS:
C A=B & B=0
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MSWAP(A,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION A(NDIM,4),B(NDIM,4)
C...............................................................
      DO ICOMP=1,4
C...............................................................
        DO JDIM=1,NDIM
C...............................................................
          A(JDIM,ICOMP)=B(JDIM,ICOMP)
          B(JDIM,ICOMP)=0D0
C...............................................................
        ENDDO
C...............................................................
      ENDDO
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MCONVKICK2.
C=================================================================
C GENERATES POLYNOMIALS FOR A KICK MAP TO COMPUTE THE TRANSFER MAP
C USING OUR COMPOSITION PROGRAM.
C TWO DIMENSIONAL VERSION.
C IL IS THE LINEAR BLOCK NUMBER.
C IORD IS THE MAXIMUM ORDER ALLOWED.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MCONVKICK2(IL,IORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      CHARACTER*10 ALBL
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
C...............................................................
      IDIME=2               !..DIMENSION OF PHASE SPACE
C...............................................................
      CALL MONE(RKICK,IDIME)    !..SETS KICK TO IDENTITY
C...............................................................
      DO JMULT=INL(1,IL),INL(2,IL)  !..MULTIPOLE COMPONENTS
C...............................................................
        JORD=INT(A(JMULT,3))
        IF (JORD.GT.IORD) GOTO 10
        RKICK(IPUNT(JORD,0,0,0),2)=-A(JMULT,1)
C...............................................................
      ENDDO
C...............................................................
 10   RETURN
C...............................................................
      END
CDECK  ID>, MCONVKICK4.
C=================================================================
C GENERATES POLYNOMIALS FOR A KICK MAP TO COMPUTE THE TRANSFER MAP
C USING OUR COMPOSITION PROGRAM.
C FOUR DIMENSIONAL VERSION.
C IL IS THE LINEAR BLOCK NUMBER.
C IORD IS THE MAXIMUM ORDER ALLOWED.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MCONVKICK4(IL,IORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      CHARACTER*10 ALBL
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
      COMMON/FACTOR/FAC(0:11)
C...............................................................
      IDIME=4               !..DIMENSION OF PHASE SPACE
C...............................................................
      CALL MZERO(RKICK,IDIME)   !..SETS KICK TO ZERO
C      CALL MONE(RKICK,IDIME)    !..SETS KICK TO IDENTITY
C...............................................................
      DO JMULT=INL(1,IL),INL(2,IL)  !..MULTIPOLE COMPONENTS
C...............................................................
        JORD=INT(A(JMULT,3))
        IF (JORD.GT.IORD) GOTO 10
        POT=-1D0
C...............................................................
        DO K=0,JORD
          COEF=FAC(JORD)/(FAC(K)*FAC(JORD-K))
          IF ((K/2)*2.EQ.K) THEN
            POT=POT*(-1D0)
            RKICK(IPUNT(JORD-K,0,K,0),2)= -(COEF*A(JMULT,1)*POT)
C...............................................................
            RKICK(IPUNT(JORD-K,0,K,0),4)= -(COEF*A(JMULT,2)*POT)
C...............................................................
          ELSE
            RKICK(IPUNT(JORD-K,0,K,0),2)= -(COEF*A(JMULT,2)*POT)
C...............................................................
            RKICK(IPUNT(JORD-K,0,K,0),4)= (COEF*A(JMULT,1)*POT)
C...............................................................
          ENDIF
C...............................................................
        ENDDO
C...............................................................
      ENDDO
C...............................................................
 10   RETURN
C...............................................................
      END
CDECK  ID>, MPRINTD.
C=================================================================
C PRINTS THE JACOBIAN OF A TRANSFER MAP USING BERZ FORMAT.
C AMAP IS THE ARRAY CONTAINING THE JACOBIAN COEFFICIENTS.
C IORD IS THE MAXIMUM ORDER ALLOWED.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C IUNIT IS THE FORTRAN UNIT USED TO WRITE THE RESULTS.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MPRINTD(AMAP,IORD,IDIME,IUNIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION AMAP(NDIM,16)
      INTEGER IN(16)
      CHARACTER*10 NAME(16)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MPRINTD): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MPRINTD): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IUNIT.LT.0.OR.IUNIT.GT.100) THEN
        WRITE(6,*) '***ERROR(MPRINTD): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      EPS=1D-30
C...............................................................
      IN(1)=112
      IN(2)=120
      IN(3)=113
      IN(4)=121
      IN(5)=114
      IN(6)=122
      IN(7)=115
      IN(8)=123
      IN(9)=116
      IN(10)=124
      IN(11)=117
      IN(12)=125
      IN(13)=118
      IN(14)=126
      IN(15)=119
      IN(16)=127
C...............................................................
      IDIME2=IDIME*IDIME
C...............................................................
      IF (IDIME.EQ.2) THEN         !..2D TRANSFER MAP
C...............................................................
        NAME(1)='DX1/DX'
        NAME(2)='DX1/DPX'
        NAME(3)='DPX1/DX'
        NAME(4)='DPX1/DPX'
C...............................................................
        DO ICOMP=1,IDIME2           !..LOOP OVER THE COMPONENTS
          WRITE(IUNIT,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)')
     .           NAME(ICOMP),', NO =',IORD,', NV =',IDIME,', INA =',
     .           IN(ICOMP),
     .              '***********'//'**********************************'
          WRITE(IUNIT,'(A)')
     .    '    I  COEFFICIENT          ORDER   EXPONENTS'
C...............................................................
          ICOEF=0
          DO JORD=0,IORD
            DO I1=0,JORD
              I2=JORD-I1
              IF (ABS(AMAP(IPUNT(I2,I1,0,0),ICOMP)).GT.EPS) THEN
                ICOEF=ICOEF+1
                WRITE(IUNIT,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')
     .          ICOEF,AMAP(IPUNT(I2,I1,0,0),ICOMP),JORD,I2,I1
                WRITE(IUNIT,*) AMAP(IPUNT(I2,I1,0,0),ICOMP)
              ENDIF
            ENDDO
          ENDDO
          WRITE(IUNIT,'(A)') '                                      '
        ENDDO
      ELSEIF (IDIME.EQ.4) THEN     !..4D TRANSFER MAP
C...............................................................
        NAME(1)='DX1/DX'
        NAME(2)='DX1/DPX'
        NAME(3)='DX1/DY'
        NAME(4)='DX1/DPY'
        NAME(5)='DPX1/DX'
        NAME(6)='DPX1/DPX'
        NAME(7)='DPX1/DY'
        NAME(8)='DPX1/DPY'
        NAME(9)='DY1/DX'
        NAME(10)='DY1/DPX'
        NAME(11)='DY1/DY'
        NAME(12)='DY1/DPY'
        NAME(13)='DPY1/DX'
        NAME(14)='DPY1/DPX'
        NAME(15)='DPY1/DY'
        NAME(16)='DPY1/DPY'
C...............................................................
        DO ICOMP=1,IDIME2           !..LOOP OVER THE COMPONENTS
          WRITE(IUNIT,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)')
     .           NAME(ICOMP),', NO =',IORD,', NV =',IDIME,', INA =',
     .           IN(ICOMP),
     .              '***********'//'**********************************'
          WRITE(IUNIT,'(A)')
     .    '    I  COEFFICIENT          ORDER   EXPONENTS'
C...............................................................
          ICOEF=0
          DO JORD=0,IORD
            DO I1=0,JORD
              DO I2=0,JORD-I1
                DO I3=0,JORD-I1-I2
                  I4=JORD-I1-I2-I3
                  IF (ABS(AMAP(IPUNT(I4,I3,I2,I1),ICOMP)).GT.EPS) THEN
                    ICOEF=ICOEF+1
                    WRITE(IUNIT,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')
     .              ICOEF,AMAP(IPUNT(I4,I3,I2,I1),ICOMP),JORD,I4,
     .              I3,I2,I1
                    WRITE(IUNIT,*) AMAP(IPUNT(I4,I3,I2,I1),ICOMP)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          WRITE(IUNIT,'(A)') '                                      '
        ENDDO
      ELSE
        WRITE(6,*) '***ERROR(MPRINT): WRONG DIMENSION OF PHASE SPACE'
      ENDIF
C...............................................................
      CLOSE(IUNIT)
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MPRINTAZAN.
C=================================================================
C  PRINTS OUT ON UNIT IUNIT THE HAMILTONIAN IN ACTION ANGLE 
C  VARIABLES. THE HAMILTONIAN MUST HAVE BEEN PREVIOUSLY COMPUTED 
C  (SUBROUTINE HAZANFILL) AND STORED IN THE ARRAY A*, ALPHA* IN 
C  THE COMMON BLOCK /HAMILT/
C
C  AUTHOR: E. TODESCO - INFN
C
C
 
      SUBROUTINE MPRINTAZAN(IUNIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      COMPLEX*16 AUTORES
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
      COMMON/RISONANZA/IRES1(2),IRES2(2),AUTORES(4),IRESON,ICASE
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C...............................................................
      IF (IUNIT.LT.0.OR.IUNIT.GT.100) THEN
        WRITE(6,*) '***ERROR(MPRINTAZAN): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
 664  FORMAT(1H ,2D24.16,4X,'(',2I2,')')
 665  FORMAT(1H ,2D24.16,4X,'(',3I2,')')
 666  FORMAT(1H ,2D24.16,4X,'(',4I2,')')
      WRITE(IUNIT,*) NMAXH
      IF(ICASE.EQ.0)  THEN
C...............................................NONRESONANT CASE 
        DO N=2,NMAXH,2
          DO K=0,N/2
             WRITE(IUNIT,664) A0(K,N/2-K),ALPHA0(K,N/2-K),K,N/2-K
          END DO
        END DO
      ELSE IF(ICASE.EQ.1)  THEN
C...............................................SINGLE RESONANCE
        DO N=1,NMAXH
          DO L=0,N/KSUM
            DO K1=0,(N-L*KSUM)/2
              K2=(N-L*KSUM-2*K1)/2
              IF(2*K2.EQ.N-L*KSUM-2*K1) THEN
                WRITE(IUNIT,665) A1(K1,K2,L),ALPHA1(K1,K2,L),K1,K2,L
              END IF
            END DO
          END DO
        END DO
      ELSE IF(ICASE.EQ.2)  THEN
C...............................................DOUBLE RESONANCE
        DO N=2,NMAXH
          DO L1=0,N
            DO L2=0,N
              I1=ABS(L1*IQ1+L2*IP2)
              I2=ABS(L1*IP1+L2*IQ2)
              ISCR=(N-I1-I2)/2
              IF(N.GE.I1+I2.AND.ISCR*2.EQ.N-I1-I2) THEN
                DO K1=0,(N-I1-I2)/2
                  K2=(N-I1-I2-2*K1)/2  
                  WRITE(IUNIT,666) A2(K1,K2,L1,L2),
     .                 ALPHA2(K1,K2,L1,L2),K1,K2,L1,L2
                END DO
              END IF
            END DO
          END DO
          DO L1=1,N
            DO L2=1,N
              I1=ABS(L1*IQ1-L2*IP2)
              I2=ABS(L1*IP1-L2*IQ2)
              ISCR=(N-I1-I2)/2
              IF(N.GE.I1+I2.AND.ISCR*2.EQ.N-I1-I2) THEN
                DO K1=0,(N-I1-I2)/2
                  K2=(N-I1-I2-2*K1)/2  
                  WRITE(IUNIT,666) A2(K1,K2,L1,-L2),
     .                 ALPHA2(K1,K2,L1,-L2),K1,K2,L1,-L2
                END DO
              END IF
            END DO
          END DO
        END DO
      END IF
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MHAZANFILL.
C=================================================================
C  TRANSFORMS THE HAMILTONIAN FROM COMPLEX COORDINATES TO ACTION 
C  ANGLE VARIABLES FIRSTLY, THE COMMON BLOCK /ORDER/ IS FILLED 
C  WITH THE INFORMATION RELATIVE TO THE TYPE OF RESONANCE. THEN 
C  THE CONVERSION IS CARRIED OUT AND THE COMMON BLOCK /HAMILT/ IS 
C  FILLED WITH THE HAMILTONIAN COEFFICIENTS
C  ICASE=0: NONRESONANT CASE
C  ICASE=1: SINGLE RESONANCE
C  ICASE=2: DOUBLE RESONANCE
C
C  AUTHOR: E. TODESCO - INFN
C
 
      SUBROUTINE MHAZANFILL(H,IORD,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      COMPLEX*16 H(NDIMPLUS),COEFF,AUTORES
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
      COMMON/RISONANZA/IRES1(2),IRES2(2),AUTORES(4),IRESON,ICASE
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX1) THEN
        WRITE(6,*) '***ERROR(MHAZANFILL): SECOND PARAMETER OUT OF ',
     .             'BOUNDS'
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MHAZANFILL): THIRD PARAMETER OUT OF ',
     .             'BOUNDS'
        STOP
      ENDIF
C...............................................................
      PI=4.D0*DATAN(1.D0)
      IF (IDIME.NE.4) RETURN
C......................INITIALIZATION OF THE COMMON BLOCS HAZAN*
      DO I1=0,NORMAX1
        DO I2=0,NORMAX1
          A0(I1,I2)=0D0
          ALPHA0(I1,I2)=0D0
          DO I3=0,NORMAX1
            A1(I1,I2,I3)=0D0
            ALPHA1(I1,I2,I3)=0D0
C...............................................................
            DO I4=-NORMAX1,NORMAX1
C...............................................................
              A2(I1,I2,I3,I4)=0D0
              ALPHA2(I1,I2,I3,I4)=0D0
C...............................................................
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C.......................................FILLS THE COMMON /ORDER/
      NMAX=IORD-1
      NMAXH=IORD
      K11=IRES1(1)
      K12=IRES2(1)
      K21=IRES1(2)
      K22=IRES2(2)
C...............................................................
      IF(ICASE.EQ.1) THEN
        IF(K11.LT.0) THEN
          K11=-K11
          K12=-K12
        END IF
        IF(K11.EQ.0.AND.K12.LT.0) K12=-K12
        IF(K11.EQ.0.AND.K12.EQ.0) 
     .         WRITE(*,*) '***ERROR IN THE RESONANT VECTORS'
        IAK12=ABS(K12)
        KSUM=K11+IAK12
      ELSE IF(ICASE.EQ.2) THEN
        IF(K11.EQ.0.OR.K22.EQ.0) THEN
          ISCR1=K11
          ISCR2=K12
          K11=K21
          K12=K22
          K21=ISCR1
          K22=ISCR2
        END IF
        IF (K11*K12.LT.0) THEN
          K12=-ABS(K12)
        ELSE
          K12=ABS(K12)
        END IF
        K11=ABS(K11)
        IF (K21*K22.LT.0) THEN
          K21=-ABS(K21)
        ELSE
          K21=ABS(K21)
        END IF
        K22=ABS(K22)
      END IF
      OR1=DATAN2(DIMAG(AUTORES(1)),DREAL(AUTORES(1)))/2/PI
      IF(OR1.LT.0) OR1=1+OR1
      OR2=DATAN2(DIMAG(AUTORES(3)),DREAL(AUTORES(3)))/2/PI
      IF(OR2.LT.0) OR2=1+OR2
C...............................................................
      IF(ICASE.EQ.0)  THEN
C...............................................NONRESONANT CASE
        DO N=2,NMAXH,2
          DO K=0,N/2
             a0(K,N/2-K)=CDABS(H(IPUNT(K,K,N/2-K,N/2-K))*(0.D0,-1.D0))
             ALPHA0(K,N/2-K)= DATAN2(DIMAG(H(IPUNT(K,K,N/2-K,N/2-K))*
     *  (0.D0,-1.D0)),DREAL(H(IPUNT(K,K,N/2-K,N/2-K))*(0.D0,-1.D0)))
          END DO
        END DO
      ELSE IF(ICASE.EQ.1)  THEN
C...............................................SINGLE RESONANCE
        DO N=1,NMAXH
          DO L=0,N/KSUM
            IF(L.NE.0) KA=2
            IF(L.EQ.0) KA=1
            DO K1=0,(N-L*KSUM)/2
              K2=(N-L*KSUM-2*K1)/2
              IF(2*K2.EQ.N-L*KSUM-2*K1) THEN
                A1(K1,K2,L)=0.D0
                ALPHA1(K1,K2,L)=0.D0
C..............SINGLE RESONANCE ON X,PX: RESONANT VECTOR (K11,0)
                IF(K11.GT.0.AND.K12.EQ.0) THEN
                  IF(H(IPUNT(K1+L*K11,K1,K2,K2)).NE.0.) THEN
                    A1(K1,K2,L)=CDABS(KA*(0.D0,-1.D0)*
     *                   H(IPUNT(K1+L*K11,K1,K2,K2)))
                    ALPHA1(K1,K2,L)=DATAN2(DIMAG((0.D0,-1.D0)*
     *                  H(IPUNT(K1+L*K11,K1,K2,K2))),DREAL((0.D0,
     *                   -1.D0)*H(IPUNT(K1+L*K11,K1,K2,K2))))
                  END IF
C...............SINGLE RESONANCE ON Y,PY: RESONANT VECTOR (0,K12)
                ELSE IF(K11.EQ.0.AND.K12.GT.0) THEN
                  IF(H(IPUNT(K1,K1,K2+L*K12,K2)).NE.0.) THEN
                    A1(K1,K2,L)=CDABS(KA*(0.D0,-1.D0)*
     *                   H(IPUNT(K1,K1,K2+L*K12,K2)))
                    ALPHA1(K1,K2,L)=DATAN2(DIMAG((0.D0,-1.D0)*
     *                  H(IPUNT(K1,K1,K2+L*K12,K2))),DREAL((0.D0,
     *                   -1.D0)*H(IPUNT(K1,K1,K2+L*K12,K2))))
                  END IF
C........SINGLE COUPLED SUM RESONANCE: RESONANT VECTOR (K11,K12)
                ELSE IF(K11.GT.0.AND.K12.GT.0) THEN
                  IF(H(IPUNT(K1+L*K11,K1,K2+L*K12,K2)).NE.0.) THEN
                    A1(K1,K2,L)=CDABS(KA*(0.D0,-1.D0)*
     *                   H(IPUNT(K1+L*K11,K1,K2+L*K12,K2)))
                    ALPHA1(K1,K2,L)=DATAN2(DIMAG((0.D0,-1.D0)*
     *                  H(IPUNT(K1+L*K11,K1,K2+L*K12,K2))),DREAL((0.D0,
     *                   -1.D0)*H(IPUNT(K1+L*K11,K1,K2+L*K12,K2))))
                  END IF
C.......SINGLE COUPLED DIFF. RESONANCE: RESONANT VECTOR (K11,K12)
                ELSE IF(K11.GT.0.AND.K12.LT.0) THEN
                  IF(H(IPUNT(K1+L*K11,K1,K2,K2-L*K12)).NE.0.) THEN
                    A1(K1,K2,L)=CDABS(KA*(0.D0,-1.D0)*
     *                   H(IPUNT(K1+L*K11,K1,K2,K2-L*K12)))
                    ALPHA1(K1,K2,L)=DATAN2(DIMAG((0.D0,-1.D0)*
     *                  H(IPUNT(K1+L*K11,K1,K2,K2-L*K12))),DREAL((0.D0,
     *                   -1.D0)*H(IPUNT(K1+L*K11,K1,K2,K2-L*K12))))
                  END IF
                ELSE
                  WRITE(*,*) 
     .            '***ERROR(HAZANFILL) IN THE RESONANT VECTORS'
                END IF
              END IF
            END DO
          END DO
        END DO
      ELSE 
C...............................................DOUBLE RESONANCE 
        IQ1=K11
        IP1=K12
        IAP1=ABS(IP1)
        IP2=K21
        IAP2=ABS(IP2)
        IQ2=K22
C.................................................DETUNING TERMS
        DO N=2,NMAXH,2
          DO K=0,N/2
             COEFF=H(IPUNT(K,K,N/2-K,N/2-K))*(0.D0,-1.D0)
             A2(K,N/2-K,0,0)=CDABS(COEFF)
             IF(COEFF.NE.0) THEN
               ALPHA2(K,N/2-K,0,0)= DATAN2(DIMAG(COEFF),DREAL(COEFF))
             ELSE
               ALPHA2(K,N/2-K,0,0)= 0
             END IF
          END DO
        END DO
C............................................................H_S
        I1=IAP1+IQ1
        IF(IP1.NE.IAP1) THEN
          IS1=-1
        ELSE
          IS1=1
        END IF
        DO N=2,NMAXH
          DO L1=1,N/I1
            ISCR=(N-L1*I1)/2
            IF(ISCR*2.EQ.N-L1*I1) THEN
              DO K1=0,(N-L1*I1)/2
                K2=(N-2*K1-L1*I1)/2
                IR1=L1*IAP1
                COEFF=-2*(0.,1.)*H(IPUNT(K1+L1*IQ1,K1,
     .                K2+IR1*(1+IS1)/2,K2+IR1*(1-IS1)/2))
                A2(K1,K2,L1,0)=CDABS(COEFF)
                IF(COEFF.NE.0) THEN
                  ALPHA2(K1,K2,L1,0)=DATAN2(DIMAG(COEFF),DREAL(COEFF))
                ELSE
                  ALPHA2(K1,K2,L1,0)=0
                END IF
              END DO
            END IF
          END DO
        END DO
C...............................................................
        I2=IAP2+IQ2
        IF(IP2.NE.IAP2) THEN
          IS2=-1
        ELSE
          IS2=1
        END IF
        DO N=2,NMAXH
          DO L2=1,N/I2
            ISCR=(N-L2*I2)/2
            IF(ISCR*2.EQ.N-L2*I2) THEN
              DO K1=0,(N-L2*I2)/2
                K2=(N-2*K1-L2*I2)/2
                IR2=L2*IAP2
                COEFF=-2*(0.,1.)*H(IPUNT(
     .            K1+IR2*(1+IS2)/2,K1+IR2*(1-IS2)/2,K2+L2*IQ2,K2))
                A2(K1,K2,0,L2)=CDABS(COEFF)
                IF(COEFF.NE.0) THEN
                  ALPHA2(K1,K2,0,L2)=DATAN2(DIMAG(COEFF),DREAL(COEFF))
                ELSE
                  ALPHA2(K1,K2,0,L2)=0.
                END IF
              END DO
            END IF
          END DO
        END DO
C............................................................H_D
        DO N=2,NMAXH
          DO L1=1,N
            DO L2=-N,N
              IF(L2.NE.0) THEN
                I1=ABS(L1*IQ1+L2*IP2)
                I2=ABS(L1*IP1+L2*IQ2)
                IF(I1+I2.LE.N) THEN
                  ISCR=(N-I1-I2)/2           
                  IF(ISCR*2.EQ.N-I1-I2) THEN
                    IF(L1*IQ1+L2*IP2.GE.0) IS1=1
                    IF(L1*IQ1+L2*IP2.LT.0) IS1=-1
                    IR1=ABS(L1*IQ1+L2*IP2)
                    IF(L1*IP1+L2*IQ2.GE.0) IS2=1
                    IF(L1*IP1+L2*IQ2.LT.0) IS2=-1
                    IR2=ABS(L1*IP2+L2*IQ2)
                    DO K1=0,(N-I1-I2)/2
                      K2=(N-I1-I2-2*K1)/2
                      COEFF=-2.*(0.,1.)*
     .                   H(IPUNT(K1+IR1*(1+IS1)/2,K1+IR1*(1-IS1)/2,
     .                   K2+IR2*(1+IS2)/2,K2+IR2*(1-IS2)/2))
                      A2(K1,K2,L1,L2)=CDABS(COEFF)
                      IF(COEFF.NE.0) THEN
                        ALPHA2(K1,K2,L1,L2)=DATAN2(DIMAG(COEFF),
     .                         DREAL(COEFF))
                      ELSE
                        ALPHA2(K1,K2,L1,L2)=0.
                      END IF
                    END DO
                  END IF
                END IF
              END IF
            END DO
          END DO
        END DO
C...............................................................
      END IF
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MPRINTHC.
C=================================================================
C PRINTS A COMPLEX HAMILTONIAN USING BERZ FORMAT.
C AMAP IS THE COMPLEX ARRAY CONTAINING THE HAMILTONIAN COEFFICIENTS.
C IORD IS THE MAXIMUM ORDER ALLOWED.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C IUNIT IS THE FORTRAN UNIT USED TO WRITE THE RESULTS.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MPRINTHC(AMAP,IORD,IDIME,IUNIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      COMPLEX*16 AMAP(NDIMPLUS)
      INTEGER IN(4)
      CHARACTER*10 NAME(4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MPRINTHC): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MPRINTHC): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IUNIT.LT.0.OR.IUNIT.GT.100) THEN
        WRITE(6,*) '***ERROR(MPRINTHC): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      EPS=1D-30
C...............................................................
      NAME(1)='H'
C...............................................................
      IN(1)=112
      ICOMP=1
C...............................................................
      IF (IDIME.EQ.2) THEN         !..2D HAMILTONIAN
        WRITE(IUNIT,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)')
     .       NAME(ICOMP),', NO =',IORD,', NV =',IDIME,', INA =',
     .       IN(ICOMP),
     .       '***********'//'**********************************'
        WRITE(IUNIT,'(A)')
     .     '    I            '//
     .     '          COEFFICIENT             ORDER  EXPONENTS'
C...............................................................
        ICOEF=0
        DO JORD=0,IORD
          DO I1=0,JORD
            I2=JORD-I1
            IF (CDABS(AMAP(IPUNT(I2,I1,0,0))).GT.EPS) THEN
              ICOEF=ICOEF+1
              WRITE(IUNIT,10)
     .        ICOEF,AMAP(IPUNT(I2,I1,0,0)),JORD,I2,I1
              WRITE(IUNIT,*) AMAP(IPUNT(I2,I1,0,0))
            ENDIF
          ENDDO
        ENDDO
        WRITE(IUNIT,'(A)') '                                      '
      ELSEIF (IDIME.EQ.4) THEN     !..4D HAMILTONIAN
        WRITE(IUNIT,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)')
     .       NAME(ICOMP),', NO =',IORD,', NV =',IDIME,', INA =',
     .       IN(ICOMP),
     .       '***********'//'**********************************'
        WRITE(IUNIT,'(A)')
     .     '    I            '//
     .     '          COEFFICIENT             ORDER  EXPONENTS'
C...............................................................
        ICOEF=0
        DO JORD=0,IORD
          DO I1=0,JORD
            DO I2=0,JORD-I1
              DO I3=0,JORD-I1-I2
                I4=JORD-I1-I2-I3
                IF (CDABS(AMAP(IPUNT(I4,I3,I2,I1))).GT.EPS) THEN
                  ICOEF=ICOEF+1
                  WRITE(IUNIT,10)
     .                 ICOEF,AMAP(IPUNT(I4,I3,I2,I1)),JORD,I4,
     .                 I3,I2,I1
                  WRITE(IUNIT,*) AMAP(IPUNT(I4,I3,I2,I1))
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        WRITE(IUNIT,'(A)') '                                      '
      ELSE
        WRITE(6,*) '***ERROR(MPRINT): WRONG DIMENSION OF PHASE SPACE'
      ENDIF
C...............................................................
      CLOSE(IUNIT)
C...............................................................
 10   FORMAT(I6,2X,G20.14,1X,G20.14,I5,4X,18(2I2,1X))
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MPRINTC.
C=================================================================
C PRINTS A COMPLEX TRANSFER MAP USING BERZ FORMAT.
C AMAP IS THE COMPLEX ARRAY CONTAINING THE MAP COEFFICIENTS.
C IORD IS THE MAXIMUM ORDER ALLOWED.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C IUNIT IS THE FORTRAN UNIT USED TO WRITE THE RESULTS.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MPRINTC(AMAP,IORD,IDIME,IUNIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      COMPLEX*16 AMAP(NDIM,4)
      INTEGER IN(4)
      CHARACTER*10 NAME(4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MPRINTC): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MPRINTC): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IUNIT.LT.0.OR.IUNIT.GT.100) THEN
        WRITE(6,*) '***ERROR(MPRINTC): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      EPS=1D-30
C...............................................................
      NAME(1)='X'
      NAME(2)='PX'
      NAME(3)='Y'
      NAME(4)='PY'
C...............................................................
      WRITE(NAME(1)(6:10),'(I5)') 1
      WRITE(NAME(2)(6:10),'(I5)') 1
      WRITE(NAME(3)(6:10),'(I5)') 2
      WRITE(NAME(4)(6:10),'(I5)') 2
C...............................................................
      IN(1)=112
      IN(2)=114
      IN(3)=113
      IN(4)=115
C...............................................................
      IF (IDIME.EQ.2) THEN         !..2D TRANSFER MAP
        DO ICOMP=1,IDIME           !..LOOP OVER THE COMPONENTS
          WRITE(IUNIT,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)')
     .           NAME(ICOMP),', NO =',IORD,', NV =',IDIME,', INA =',
     .           IN(ICOMP),
     .              '***********'//'**********************************'
          WRITE(IUNIT,'(A)')
     .     '    I            '//
     .     '          COEFFICIENT             ORDER  EXPONENTS'
C...............................................................
          ICOEF=0
          DO JORD=0,IORD
            DO I1=0,JORD
              I2=JORD-I1
              IF (CDABS(AMAP(IPUNT(I2,I1,0,0),ICOMP)).GT.EPS) THEN
                ICOEF=ICOEF+1
                WRITE(IUNIT,10)
     .          ICOEF,AMAP(IPUNT(I2,I1,0,0),ICOMP),JORD,I2,I1
                WRITE(IUNIT,*) AMAP(IPUNT(I2,I1,0,0),ICOMP)
              ENDIF
            ENDDO
          ENDDO
          WRITE(IUNIT,'(A)') '                                      '
        ENDDO
      ELSEIF (IDIME.EQ.4) THEN     !..4D TRANSFER MAP
        DO ICOMP=1,IDIME           !..LOOP OVER THE COMPONENTS
          WRITE(IUNIT,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)')
     .           NAME(ICOMP),', NO =',IORD,', NV =',IDIME,', INA =',
     .           IN(ICOMP),
     .              '***********'//'**********************************'
          WRITE(IUNIT,'(A)')
     .     '    I            '//
     .     '          COEFFICIENT             ORDER  EXPONENTS'
C...............................................................
          ICOEF=0
          DO JORD=0,IORD
            DO I1=0,JORD
              DO I2=0,JORD-I1
                DO I3=0,JORD-I1-I2
                  I4=JORD-I1-I2-I3
                  IF (CDABS(AMAP(IPUNT(I4,I3,I2,I1),ICOMP)).GT.EPS) THEN
                    ICOEF=ICOEF+1
                    WRITE(IUNIT,10)
     .              ICOEF,AMAP(IPUNT(I4,I3,I2,I1),ICOMP),JORD,I4,
     .              I3,I2,I1
                    WRITE(IUNIT,*) AMAP(IPUNT(I4,I3,I2,I1),ICOMP)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          WRITE(IUNIT,'(A)') '                                      '
        ENDDO
      ELSE
        WRITE(6,*) '***ERROR(MPRINT): WRONG DIMENSION OF PHASE SPACE'
      ENDIF
C...............................................................
      CLOSE(IUNIT)
C...............................................................
 10   FORMAT(I6,2X,G20.14,1X,G20.14,I5,4X,18(2I2,1X))
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MPRINTH.
C=================================================================
C PRINTS A HAMILTONIAN USING BERZ FORMAT.
C AMAP IS THE ARRAY CONTAINING THE HAMILTONIAN COEFFICIENTS.
C IORD IS THE MAXIMUM ORDER ALLOWED.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C IUNIT IS THE FORTRAN UNIT USED TO WRITE THE RESULTS.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MPRINTH(AMAP,IORD,IDIME,IUNIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DOUBLE PRECISION AMAP(NDIMPLUS)
      INTEGER IN(4)
      CHARACTER*10 NAME(4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MPRINTH): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MPRINTH): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IUNIT.LT.0.OR.IUNIT.GT.100) THEN
        WRITE(6,*) '***ERROR(MPRINTH): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      EPS=1D-30
C...............................................................
      NAME(1)='H'
C...............................................................
      IN(1)=112
      ICOMP=1
C...............................................................
      IF (IDIME.EQ.2) THEN         !..2D HAMILTONIAN
        WRITE(IUNIT,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)')
     .       NAME(ICOMP),', NO =',IORD,', NV =',IDIME,', INA =',
     .       IN(ICOMP),
     .       '***********'//'**********************************'
        WRITE(IUNIT,'(A)')
     .    '    I  COEFFICIENT          ORDER   EXPONENTS'
C...............................................................
        ICOEF=0
        DO JORD=0,IORD
          DO I1=0,JORD
            I2=JORD-I1
            IF (ABS(AMAP(IPUNT(I2,I1,0,0))).GT.EPS) THEN
              ICOEF=ICOEF+1
              WRITE(IUNIT,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')
     .        ICOEF,AMAP(IPUNT(I2,I1,0,0)),JORD,I2,I1
              WRITE(IUNIT,*) AMAP(IPUNT(I2,I1,0,0))
            ENDIF
          ENDDO
        ENDDO
        WRITE(IUNIT,'(A)') '                                      '
      ELSEIF (IDIME.EQ.4) THEN     !..4D HAMILTONIAN
        WRITE(IUNIT,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)')
     .       NAME(ICOMP),', NO =',IORD,', NV =',IDIME,', INA =',
     .       IN(ICOMP),
     .       '***********'//'**********************************'
        WRITE(IUNIT,'(A)')
     .    '    I  COEFFICIENT          ORDER   EXPONENTS'
C...............................................................
        ICOEF=0
        DO JORD=0,IORD
          DO I1=0,JORD
            DO I2=0,JORD-I1
              DO I3=0,JORD-I1-I2
                I4=JORD-I1-I2-I3
                IF (ABS(AMAP(IPUNT(I4,I3,I2,I1))).GT.EPS) THEN
                  ICOEF=ICOEF+1
                  WRITE(IUNIT,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')
     .                 ICOEF,AMAP(IPUNT(I4,I3,I2,I1)),JORD,I4,
     .                 I3,I2,I1
                  WRITE(IUNIT,*) AMAP(IPUNT(I4,I3,I2,I1))
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        WRITE(IUNIT,'(A)') '                                      '
      ELSE
        WRITE(6,*) '***ERROR(MPRINT): WRONG DIMENSION OF PHASE SPACE'
      ENDIF
C...............................................................
      CLOSE(IUNIT)
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MPRINT.
C=================================================================
C PRINTS A REAL TRANSFER MAP USING BERZ FORMAT.
C AMAP IS THE REAL ARRAY CONTAINING THE MAP COEFFICIENTS.
C IORD IS THE MAXIMUM ORDER ALLOWED.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C IUNIT IS THE FORTRAN UNIT USED TO WRITE THE RESULTS.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MPRINT(AMAP,IORD,IDIME,IUNIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION AMAP(NDIM,4)
      INTEGER IN(4)
      CHARACTER*10 NAME(4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MPRINT): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MPRINT): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IUNIT.LT.0.OR.IUNIT.GT.100) THEN
        WRITE(6,*) '***ERROR(MPRINT): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      EPS=1D-30
C...............................................................
      NAME(1)='X'
      NAME(2)='PX'
      NAME(3)='Y'
      NAME(4)='PY'
C...............................................................
      WRITE(NAME(1)(6:10),'(I5)') 1
      WRITE(NAME(2)(6:10),'(I5)') 1
      WRITE(NAME(3)(6:10),'(I5)') 2
      WRITE(NAME(4)(6:10),'(I5)') 2
C...............................................................
      IN(1)=112
      IN(2)=114
      IN(3)=113
      IN(4)=115
C...............................................................
      IF (IDIME.EQ.2) THEN         !..2D TRANSFER MAP
        DO ICOMP=1,IDIME           !..LOOP OVER THE COMPONENTS
          WRITE(IUNIT,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)')
     .           NAME(ICOMP),', NO =',IORD,', NV =',IDIME,', INA =',
     .           IN(ICOMP),
     .              '***********'//'**********************************'
          WRITE(IUNIT,'(A)')
     .    '    I  COEFFICIENT          ORDER   EXPONENTS'
C...............................................................
          ICOEF=0
          DO JORD=0,IORD
            DO I1=0,JORD
              I2=JORD-I1
              IF (ABS(AMAP(IPUNT(I2,I1,0,0),ICOMP)).GT.EPS) THEN
                ICOEF=ICOEF+1
                WRITE(IUNIT,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')
     .          ICOEF,AMAP(IPUNT(I2,I1,0,0),ICOMP),JORD,I2,I1
                WRITE(IUNIT,*) AMAP(IPUNT(I2,I1,0,0),ICOMP)
              ENDIF
            ENDDO
          ENDDO
          WRITE(IUNIT,'(A)') '                                      '
        ENDDO
      ELSEIF (IDIME.EQ.4) THEN     !..4D TRANSFER MAP
        DO ICOMP=1,IDIME           !..LOOP OVER THE COMPONENTS
          WRITE(IUNIT,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)')
     .           NAME(ICOMP),', NO =',IORD,', NV =',IDIME,', INA =',
     .           IN(ICOMP),
     .              '***********'//'**********************************'
          WRITE(IUNIT,'(A)')
     .    '    I  COEFFICIENT          ORDER   EXPONENTS'
C...............................................................
          ICOEF=0
          DO JORD=0,IORD
            DO I1=0,JORD
              DO I2=0,JORD-I1
                DO I3=0,JORD-I1-I2
                  I4=JORD-I1-I2-I3
                  IF (ABS(AMAP(IPUNT(I4,I3,I2,I1),ICOMP)).GT.EPS) THEN
                    ICOEF=ICOEF+1
                    WRITE(IUNIT,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')
     .              ICOEF,AMAP(IPUNT(I4,I3,I2,I1),ICOMP),JORD,I4,
     .              I3,I2,I1
                    WRITE(IUNIT,*) AMAP(IPUNT(I4,I3,I2,I1),ICOMP)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          WRITE(IUNIT,'(A)') '                                      '
        ENDDO
      ELSE
        WRITE(6,*) '***ERROR(MPRINT): WRONG DIMENSION OF PHASE SPACE'
      ENDIF
C...............................................................
      CLOSE(IUNIT)
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MREADC.
C=================================================================
C READS A COMPLEX TRANSFER MAP USING BERZ FORMAT.
C AMAP IS THE COMPLEX ARRAY CONTAINING THE MAP COEFFICIENTS.
C IORD IS THE MAXIMUM ORDER ALLOWED (INPUT). IORD IS
C THE ACTUAL MAP ORDER IN CASE OF DISCREPANCIES (OUTPUT).
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C IUNIT IS THE FORTRAN UNIT USED TO READ THE RESULTS.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MREADC(AMAP,IORD,IDIME,IUNIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      COMPLEX*16 AMAP(NDIM,4),C
      DIMENSION J(10)
      CHARACTER*10 C10
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MREADC): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MREADC): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IUNIT.LT.0.OR.IUNIT.GT.100) THEN
        WRITE(6,*) '***ERROR(MREADC): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      DO IVEC=1,10
        J(IVEC)=0D0
      ENDDO
C...............................................................
      IIORD=-1
C...............................................................
      DO ICOMP=1,IDIME
        READ(IUNIT,'(A10)') C10
        READ(IUNIT,'(18X,I4,7X,I4)') JORD,JDIM
        READ(IUNIT,'(A10)') C10
        READ(IUNIT,'(A10)') C10
        READ(IUNIT,'(A10)') C10
        IF (IDIME.NE.JDIM) THEN
          WRITE(6,*) '***ERROR(MREAD): WRONG PHASE SPACE DIMENSION'
          RETURN
        ENDIF
        IIN = 0
C
  10    CONTINUE
        IIN = IIN + 1
        READ(IUNIT,'(I6,2X,G20.14,1X,G20.14,I5,4X,18(2I2,1X))')
     .             II,C,IO,(J(I),I=1,IDIME)
C
        IF (II.EQ.0) GOTO 20
        READ(IUNIT,*) C
        IF (II.NE.IIN) WRITE(6,*) '***ERROR(MREAD) FILE ',
     .              'NUMBERING OUT OF ORDER '
C
        AMAP(IPUNT(J(1),J(2),J(3),J(4)),ICOMP)=C
C
        IIORD=MAX(IO,IIORD)
C
        GOTO 10
C
  20    CONTINUE
C
        AMAP(IPUNT(0,0,0,0),ICOMP)=0D0 !..THE CONSTANT TERM IS 0
C
      ENDDO
C...............................................................
      IORD=IIORD
C...............................................................
      CLOSE(IUNIT)
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MREAD.
C=================================================================
C READS A REAL TRANSFER MAP USING BERZ FORMAT.
C AMAP IS THE REAL ARRAY CONTAINING THE MAP COEFFICIENTS.
C IORD IS THE MAXIMUM ORDER ALLOWED (INPUT). IORD IS
C THE ACTUAL MAP ORDER IN CASE OF DISCREPANCIES (OUTPUT).
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C IUNIT IS THE FORTRAN UNIT USED TO READ THE RESULTS.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MREAD(AMAP,IORD,IDIME,IUNIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION AMAP(NDIM,4)
      DIMENSION J(10)
      CHARACTER*10 C10
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MREAD): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MREAD): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IUNIT.LT.0.OR.IUNIT.GT.100) THEN
        WRITE(6,*) '***ERROR(MREAD): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      DO IVEC=1,10
        J(IVEC)=0D0
      ENDDO
C...............................................................
      IIORD=-1
C...............................................................
      DO ICOMP=1,IDIME
        READ(IUNIT,'(A10)') C10
        READ(IUNIT,'(18X,I4,7X,I4)') JORD,JDIM
        READ(IUNIT,'(A10)') C10
        READ(IUNIT,'(A10)') C10
        READ(IUNIT,'(A10)') C10
        IF (IDIME.NE.JDIM) THEN
          WRITE(6,*) '***ERROR(MREAD): WRONG PHASE SPACE DIMENSION'
          RETURN
        ENDIF
        IIN = 0
C
  10    CONTINUE
        IIN = IIN + 1
        READ(IUNIT,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')
     .             II,C,IO,(J(I),I=1,IDIME)
C
        IF (II.EQ.0) GOTO 20
        READ(IUNIT,*) C
        IF (II.NE.IIN) WRITE(6,*) '***ERROR(MREAD) FILE ',
     .              'NUMBERING OUT OF ORDER '
C
        AMAP(IPUNT(J(1),J(2),J(3),J(4)),ICOMP)=C
C
        IIORD=MAX(IO,IIORD)
C
        GOTO 10
C
  20    CONTINUE
C
        AMAP(IPUNT(0,0,0,0),ICOMP)=0D0 !..THE CONSTANT TERM IS 0
C
      ENDDO
C...............................................................
      IORD=IIORD
C...............................................................
      CLOSE(IUNIT)
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MTRANSF.
C=================================================================
C COMPUTES THE TRANSFER MAP BY COMPOSING ALL THE SINGLE ELEMENT
C MAPS.
C IORD IS THE MAXIMUM ORDER ALLOWED.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MTRANSF(IORD,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      CHARACTER*10 ALBL
      DIMENSION TEMP1(NDIM,4)
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MTRANSF): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MTRANSF): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
C      CALL MZERO(RMAP,IDIME)     !..SETS THE MAP TO ZERO
      CALL MONE(RMAP,IDIME)      !..SETS THE MAP TO UNITY
C...............................................................
      IF (IDIME.EQ.2) THEN
C
        DO IL=1,NLMAX               !..LOOP OVER THE LINEAR BLOCKS
C
          CALL MCONVMAT(IDIME,IL)   !..SETUP OF THE POLYNOMIAL
          CALL COMPOSE(RMAT,RMAP,TEMP,1,IORD,IORD) !..COMPOSITION
C          CALL MSWAP(RMAP,TEMP)     !..SWAPS RMAP<->TEMP
C
          CALL MCONVKICK2(IL,IORD)     !..2D KICK
C
          KORD=MIN(INT(A(INL(2,IL),3)),IORD)
C
          CALL COMPOSE(RKICK,TEMP,TEMP1,KORD,IORD,IORD) !..COMPOSITION
C
          CALL MSWAP(RMAP,TEMP1)    !..COPIES RMAP<->TEMP1
C
        ENDDO
C.....................................!..LAST LINEAR SECTION
        CALL MCONVMAT(IDIME,NLMAX+1)    !..SETUP OF THE POLYNOMIAL
        CALL COMPOSE(RMAT,RMAP,TEMP,1,IORD,IORD) !..COMPOSITION
        CALL MSWAP(RMAP,TEMP)           !..SWAPS RMAP<->TEMP
C
      ELSEIF (IDIME.EQ.4) THEN
C.....................................!..LAST LINEAR SECTION
        DO IL=1,NLMAX               !..LOOP OVER THE LINEAR BLOCKS
C
          CALL MCONVMAT(IDIME,IL)   !..SETUP OF THE POLYNOMIAL
          CALL COMPOSEL(RMAT,RMAP,TEMP,1,IORD,IORD) !..COMPOSITION
C          CALL MSWAP(RMAP,TEMP)     !..SWAPS RMAP<->TEMP
C
          CALL MCONVKICK4(IL,IORD)     !..4D KICK
C
          KORD=MIN(INT(A(INL(2,IL),3)),IORD)
C
          CALL COMPOSEK(RKICK,TEMP,TEMP1,KORD,IORD,IORD) !..COMPOSITION
C
          CALL MSWAP(RMAP,TEMP1)    !..COPIES RMAP<->TEMP1
        ENDDO
C...............................................................
C.....................................!..LAST LINEAR SECTION
        CALL MCONVMAT(IDIME,NLMAX+1)    !..SETUP OF THE POLYNOMIAL
        CALL COMPOSEL(RMAT,RMAP,TEMP,1,IORD,IORD) !..COMPOSITION
        CALL MSWAP(RMAP,TEMP)           !..SWAPS RMAP<->TEMP
C...............................................................
C COMPOSITION BASED ON LIE OPERATORS. INTERESTING FOR HIGH ORDER
C KICKS
C        CALL MCONVMAT(IDIME,NLMAX+1)    !..SETUP OF THE POLYNOMIAL
C        CALL MSWAP(RMAP,RMAT)           !..SWAPS RMAP<->TEMP
C        DO IL=NLMAX,1,-1
C
C          CALL MCONVKICK4(IL,IORD)     !..4D KICK
C
C          KORD=MIN(INT(A(INL(2,IL),3)),IORD)    !..ORDER OF THE KICK
C..............................OPTIMIZED VERSION USING LIE TRANSF
C          CALL LIE_TRASF(RKICK,RMAP,TEMP,KORD,IORD,IORD)
C
C
C          CALL MCONVMAT(IDIME,IL)   !..SETUP OF THE POLYNOMIAL
C          CALL COMPOSEL(TEMP,RMAT,TEMP1,IORD,1,IORD) !..COMPOSITION
C          CALL MSWAP(RMAP,TEMP1)     !..COPIES RMAP<->TEMP1
C        ENDDO
C...............................................................
      ENDIF
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MDER.
C=================================================================
C COMPUTES THE DERIVATIVE OF THE POLYNOMIAL MAP CONTAINED IN THE
C ARRAY A. THE DERIVATIVE IS STORED IN THE ARRAY RDER.
C IORD IS THE MAXIMUM ORDER ALLOWED.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MDER(A,IORD,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION A(NDIM,4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MDER): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MDER): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      CALL MZERO(RDER,IDIME)
C...............................................................
      DO JCOMP=1,IDIME         !..LOOP OVER THE COMPONENTS
        DO JORD=1,IORD
          DO I1=0,JORD
            DO I2=0,JORD-I1
              DO I3=0,JORD-I1-I2
                I4=JORD-I1-I2-I3
                IF (I4.GT.0) RDER(IPUNT(I4-1,I3,I2,I1),1+IDIME*
     .            (JCOMP-1))=FLOAT(I4)*A(IPUNT(I4,I3,I2,I1),JCOMP)
                IF (I3.GT.0) RDER(IPUNT(I4,I3-1,I2,I1),2+IDIME*
     .            (JCOMP-1))=FLOAT(I3)*A(IPUNT(I4,I3,I2,I1),JCOMP)
                IF (I2.GT.0) RDER(IPUNT(I4,I3,I2-1,I1),3+IDIME*
     .            (JCOMP-1))=FLOAT(I2)*A(IPUNT(I4,I3,I2,I1),JCOMP)
                IF (I1.GT.0) RDER(IPUNT(I4,I3,I2,I1-1),4+IDIME*
     .            (JCOMP-1))=FLOAT(I1)*A(IPUNT(I4,I3,I2,I1),JCOMP)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MRTOC.
C=================================================================
C CONVERTS A REAL ARRAY REPRESENTING A REAL POLYNOMIAL IN REAL
C COORDINATES INTO A COMPLEX POLYNOMIAL EXPRESSED IN COMPLEX
C VARIABLES.
C A IS THE REAL ARRAY (INPUT). B IS THE COMPLEX ARRAY (OUTPUT).
C THE CONVENTION IS:
C REAL VARIABLES     ===>  (X,PX,Y,PY)
C COMPLEX VARIABLES  ===>  (Z1,Z1*,Z2,Z2*)
C WITH:
C Z1=X-IPX
C Z2=Y-IPY
C AND THE COMPLEX CONJUGATES.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MRTOC(A,B,IORD,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION A(NDIM,4)
      COMPLEX*16 B(NDIM,4),C(NDIM,4),R(NDIM,4),UNITY
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MRTOC): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MRTOC): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      NTERM=IPUNT(IORD,0,0,0)   !..NUMBER OF TERMS
C...............................................................
      CALL MZEROC(B,IDIME)      !..SETS THE ARRAY B TO ZERO
C...............................................................
      DO ICOMP=1,IDIME
        DO IEL=1,NTERM
          B(IEL,ICOMP)=DCMPLX(A(IEL,ICOMP),0D0)
        ENDDO
      ENDDO
C...............................................................
      CALL MZEROC(R,IDIME)    !..SETS THE ARRAY R TO ZERO
C...............................................................
      UNITY=DCMPLX(0D0,1D0)
      R(5,1)=.5D0             !..LINEAR TRANSFORMATION IN CPLX VAR.
      R(4,1)=.5D0
      R(5,2)=UNITY*.5D0
      R(4,2)=-UNITY*.5D0
      R(3,3)=.5D0
      R(2,3)=.5D0
      R(3,4)=UNITY*.5D0
      R(2,4)=-UNITY*.5D0
C...............................................................
      NR=1
      CALL COMPOSEC(B,R,C,IORD,NR,IORD)  !..COMPOSITION
C...............................................................
      DO ICONT=1, NTERM                  !..FINAL CHANGE OF VARIABLES
         B(ICONT,1)=C(ICONT,1)-UNITY*C(ICONT,2)
         B(ICONT,2)=C(ICONT,1)+UNITY*C(ICONT,2)
         B(ICONT,3)=C(ICONT,3)-UNITY*C(ICONT,4)
         B(ICONT,4)=C(ICONT,3)+UNITY*C(ICONT,4)
      ENDDO
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MCTOR.
C=================================================================
C CONVERTS A COMPLEX ARRAY REPRESENTING A COMPLEX POLYNOMIAL IN
C COMPLEX COORDINATES INTO A REAL POLYNOMIAL EXPRESSED IN REAL
C VARIABLES.
C A IS THE COMPLEX ARRAY (INPUT). B IS THE REAL ARRAY (OUTPUT).
C THE CONVENTION IS:
C REAL VARIABLES     ===>  (X,PX,Y,PY)
C COMPLEX VARIABLES  ===>  (Z1,Z1*,Z2,Z2*)
C WITH:
C Z1=X-IPX
C Z2=Y-IPY
C AND THE COMPLEX CONJUGATES.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MCTOR(A,B,IORD,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION B(NDIM,4)
      COMPLEX*16 A(NDIM,4),C(NDIM,4),R(NDIM,4),UNITY
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MCTOR): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MCTOR): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      NTERM=IPUNT(IORD,0,0,0)   !..NUMBER OF TERMS
C...............................................................
      CALL MZERO(B,IDIME)     !..SETS THE ARRAY B TO ZERO
      CALL MZEROC(C,IDIME)    !..SETS THE ARRAY C TO ZERO
      CALL MZEROC(R,IDIME)    !..SETS THE ARRAY R TO ZERO
C...............................................................
      UNITY=DCMPLX(0D0,1D0)
      R(5,1)=1D0             !..LINEAR TRANSFORMATION IN REAL VAR.
      R(4,1)=-UNITY
      R(5,2)=1D0
      R(4,2)=UNITY
      R(3,3)=1D0
      R(2,3)=-UNITY
      R(3,4)=1D0
      R(2,4)=UNITY
C...............................................................
      NR=1
      CALL COMPOSEC(A,R,C,IORD,NR,IORD)  !..COMPOSITION
C...............................................................
      DO ICONT=1,NTERM                  !..FINAL CHANGE OF VARIABLES
         A(ICONT,1)=.5D0*C(ICONT,1)+.5D0*C(ICONT,2)
         A(ICONT,2)=.5D0*UNITY*C(ICONT,1)-.5D0*UNITY*C(ICONT,2)
         A(ICONT,3)=.5D0*C(ICONT,3)+.5D0*C(ICONT,4)
         A(ICONT,4)=.5D0*UNITY*C(ICONT,3)-.5D0*UNITY*C(ICONT,4)
      ENDDO
C...............................................................
      DO ICOMP=1,IDIME
        DO IEL=1,NTERM
          B(IEL,ICOMP)=DREAL(A(IEL,ICOMP))
        ENDDO
      ENDDO
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MPHASHIFT.
C=================================================================
C INTRODUCES A PHASE SHIFT IN THE ACCELERATOR STRUCTURE.
C TUNEH IS THE HORIZONTAL TUNE.
C TUNEV IS THE VERTICAL TUNE.
C IDIME IS THE PHASE-SPACE DIMENSION.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MPHASHIFT(TUNEH,TUNEV,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION AMAT(4,4),BMAT(4,4),CMAT(4,4),RMAT(4,4),PARAM(20)
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
      COMMON/ACCMAT/RA(4,4,0:NELE),IMAX2,IUPD,IUSE
C...............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MPHASHIFT): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      IF (TUNEH.EQ.0D0) RETURN
      IF (TUNEV.EQ.0D0.AND.IDIME.EQ.4) RETURN
C...........................................TUNE CALCULATIONS
      CALL EDTENG(RA(1,1,IMAX2),PARAM)
      T1=PARAM(1)
      T2=PARAM(4)
C...............................................................
      IF (T1.NE.TUNEH.OR.T2.NE.TUNEV) THEN
C...............................................................
        D1=TUNEH-T1              !..HOR DIFFERENCE
        D2=TUNEV-T2              !..VER DIFFERENCE
        IF (DABS(D1).GT.5D0) D1=1D0-D1
        IF (DABS(D2).GT.5D0) D2=1D0-D1
C...............................................................
        PI2=8D0*DATAN(1D0)
C...............................................................
        CM1=DCOS(PI2*D1)
        SM1=DSIN(PI2*D1)
C...............................................................
        CM2=DCOS(PI2*D2)
        SM2=DSIN(PI2*D2)
C...........................................ROTATION MATRIX
        CALL ZEROM(CMAT)
C...............................................................
        CMAT(1,1)=CM1        !..FIRST BLOCK
        CMAT(1,2)=SM1
        CMAT(2,1)=-SM1
        CMAT(2,2)=CM1
        CMAT(3,3)=CM2        !..SECOND BLOCK
        CMAT(3,4)=SM2
        CMAT(4,3)=-SM2
        CMAT(4,4)=CM2
C...............................................................
        CALL ZEROM(AMAT)
        CALL ZEROM(BMAT)
        CALL MCOPY(NLMAX+1,AMAT)
        CALL VMMUL(BMAT,CMAT,AMAT)   !..LAST BLOCK OK
        CALL VMCOPY(RL(1,1,NLMAX+1),BMAT) !..NEW BLOCK
C......................................................TEST
        CALL ZEROM(AMAT)
        CALL ZEROM(BMAT)
        CALL ZEROM(CMAT)
        CALL UNITM(RMAT)
        DO IBLOCK=1,NLMAX+1
          CALL VMMUL(BMAT,RL(1,1,IBLOCK),RMAT)  !..MATRIX MULTIPLICATION
          CALL VMCOPY(RMAT,BMAT)
          CALL ZEROM(BMAT)
        ENDDO
C...........................................TUNE CALCULATIONS
      ENDIF
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MTRACK2.
C=================================================================
C TRACKS A POLYNOMIAL MAP.
C X1 ---> INPUT  ---> CONTAINS THE INITIAL COORDINATES.
C X1 ---> OUTPUT ---> CONTAINS THE FINAL COORDINATES.
C NTURN IS THE MAXIMUM TURN NUMBER.
C IORD IS THE MAXIMUM ORDER ALLOWED.
C N IS A FLAG:
C N >  0 ITERATION OF THE MAP RMAP.
C N =  0 FIRST ROW OF JACOBIAN.
C N = -1 SECOND ROW OF JACOBIAN.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      INTEGER FUNCTION MTRACK2(X1,NTURN,IORD,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION XPOT(4,0:20),XTEMP(4),X1(4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MTRACK2): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (N.LT.-1) THEN
        WRITE(6,*) '***ERROR(MTRACK2): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      IF (N.GT.0) THEN            !..MAP TRACKING
C...............................................................
        DO ICOMP=1,2                  !..INITIALIZATION
          XTEMP(ICOMP)=0D0
          DO IPOT=1,IORD
            XPOT(ICOMP,IPOT)=0D0
          ENDDO
        ENDDO
C...............................................................
        DO ITURN=1,NTURN
          XPOT(1,0)=1D0
          XPOT(2,0)=1D0
          DO JORD=1,IORD
            XPOT(1,JORD)=XPOT(1,JORD-1)*X1(1)
            XPOT(2,JORD)=XPOT(2,JORD-1)*X1(2)
            DO I1=0,JORD
              I2=JORD-I1
              XTEMP(1)=XTEMP(1)+RMAP(IPUNT(I2,I1,0,0),1)*
     .               XPOT(1,I2)*XPOT(2,I1)
              XTEMP(2)=XTEMP(2)+RMAP(IPUNT(I2,I1,0,0),2)*
     .               XPOT(1,I2)*XPOT(2,I1)
            ENDDO
          ENDDO
C.............................CHECKS FOR PARTICLE LOSS
          IF ((ABS(XTEMP(1))+ABS(XTEMP(2))).GT.1D10) THEN
            MTRACK2=-ITURN
C            WRITE(6,*) '***MSG(MTRACK2): LOST AT ',ITURN,'/',IL
            RETURN
          ENDIF
C...............................................................
          DO ICOMP=1,2                  !..REINITIALIZATION
            X1(ICOMP)=XTEMP(ICOMP)
            XTEMP(ICOMP)=0D0
            DO IPOT=1,IORD
              XPOT(ICOMP,IPOT)=0D0
            ENDDO
          ENDDO
C...............................................................
        ENDDO
        MTRACK2=NTURN
C...............................................................
      ELSEIF (N.EQ.0) THEN       !..FIRST JACOBIAN ROW
C...............................................................
        DO ICOMP=1,2                  !..INITIALIZATION
          XTEMP(ICOMP)=0D0
          DO IPOT=1,IORD-1
            XPOT(ICOMP,IPOT)=0D0
          ENDDO
        ENDDO
C...............................................................
        DO ITURN=1,NTURN
          XPOT(1,0)=1D0
          XPOT(2,0)=1D0
          DO JORD=1,IORD-1
            XPOT(1,JORD)=XPOT(1,JORD-1)*X1(1)
            XPOT(2,JORD)=XPOT(2,JORD-1)*X1(2)
            DO I1=0,JORD
              I2=JORD-I1
              XTEMP(1)=XTEMP(1)+RDER(IPUNT(I2,I1,0,0),1)*
     .               XPOT(1,I2)*XPOT(2,I1)
              XTEMP(2)=XTEMP(2)+RDER(IPUNT(I2,I1,0,0),2)*
     .               XPOT(1,I2)*XPOT(2,I1)
            ENDDO
          ENDDO
C.............................CHECKS FOR PARTICLE LOSS
          IF ((ABS(XTEMP(1))+ABS(XTEMP(2))).GT.1D10) THEN
            MTRACK2=-ITURN
C            WRITE(6,*) '***MSG(MTRACK2): LOST AT ',ITURN,'/',IL
            RETURN
          ENDIF
C...............................................................
          DO ICOMP=1,2                  !..REINITIALIZATION
            X1(ICOMP)=XTEMP(ICOMP)
            XTEMP(ICOMP)=0D0
            DO IPOT=1,IORD-1
              XPOT(ICOMP,IPOT)=0D0
            ENDDO
          ENDDO
C...............................................................
        ENDDO
        MTRACK2=NTURN
C...............................................................
      ELSEIF (N.EQ.-1) THEN      !..SECOND JACOBIAN ROW
C...............................................................
        DO ICOMP=1,2                  !..INITIALIZATION
          XTEMP(ICOMP)=0D0
          DO IPOT=1,IORD-1
            XPOT(ICOMP,IPOT)=0D0
          ENDDO
        ENDDO
C...............................................................
        DO ITURN=1,NTURN
          XPOT(1,0)=1D0
          XPOT(2,0)=1D0
          DO JORD=1,IORD-1
            XPOT(1,JORD)=XPOT(1,JORD-1)*X1(1)
            XPOT(2,JORD)=XPOT(2,JORD-1)*X1(2)
            DO I1=0,JORD
              I2=JORD-I1
              XTEMP(1)=XTEMP(1)+RDER(IPUNT(I2,I1,0,0),3)*
     .               XPOT(1,I2)*XPOT(2,I1)
              XTEMP(2)=XTEMP(2)+RDER(IPUNT(I2,I1,0,0),4)*
     .               XPOT(1,I2)*XPOT(2,I1)
            ENDDO
          ENDDO
C.............................CHECKS FOR PARTICLE LOSS
          IF ((ABS(XTEMP(1))+ABS(XTEMP(2))).GT.1D10) THEN
            MTRACK2=-ITURN
C            WRITE(6,*) '***MSG(MTRACK2): LOST AT ',ITURN,'/',IL
            RETURN
          ENDIF
C...............................................................
          DO ICOMP=1,2                  !..REINITIALIZATION
            X1(ICOMP)=XTEMP(ICOMP)
            XTEMP(ICOMP)=0D0
            DO IPOT=1,IORD-1
              XPOT(ICOMP,IPOT)=0D0
            ENDDO
          ENDDO
C...............................................................
        ENDDO
        MTRACK2=NTURN
C...............................................................
      ENDIF
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MTRACK4.
C=================================================================
C TRACKS THE TRUNCATED 4D TRANSFER MAP.
C X1 ---> INPUT  ---> CONTAINS THE INITIAL COORDINATES.
C X1 ---> OUTPUT ---> CONTAINS THE FINAL COORDINATES.
C NTURN IS THE MAXIMUM TURN NUMBER.
C IORD IS THE MAXIMIM ORDER ALLOWED.
C N IS A FLAG:
C N >  0 ITERATION OF THE MAP RMAP.
C N =  0 FIRST ROW OF JACOBIAN.
C N = -1 SECOND ROW OF JACOBIAN.
C N = -2 THIRD ROW OF JACOBIAN.
C N = -3 FOURTH ROW OF JACOBIAN.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      INTEGER FUNCTION MTRACK4(X1,NTURN,IORD,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION XPOT(4,0:20),XTEMP(4),X1(4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MTRACK4): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (N.LT.-3) THEN
        WRITE(6,*) '***ERROR(MTRACK4): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      IF (N.GT.0) THEN                !..MAP TRACKING
C...............................................................
        DO ICOMP=1,4                  !..INITIALIZATION
           XTEMP(ICOMP)=0D0
           DO IPOT=1,IORD
             XPOT(ICOMP,IPOT)=0D0
           ENDDO
        ENDDO
C...............................................................
        DO ITURN=1,NTURN
          XPOT(1,0)=1D0
          XPOT(2,0)=1D0
          XPOT(3,0)=1D0
          XPOT(4,0)=1D0
          DO JORD=1,IORD
            XPOT(1,JORD)=XPOT(1,JORD-1)*X1(1)
            XPOT(2,JORD)=XPOT(2,JORD-1)*X1(2)
            XPOT(3,JORD)=XPOT(3,JORD-1)*X1(3)
            XPOT(4,JORD)=XPOT(4,JORD-1)*X1(4)
            DO I1=0,JORD
              DO I2=0,JORD-I1
                DO I3=0,JORD-I1-I2
                  I4=JORD-I1-I2-I3
                  XTEMP(1)=XTEMP(1)+RMAP(IPUNT(I4,I3,I2,I1),1)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(2)=XTEMP(2)+RMAP(IPUNT(I4,I3,I2,I1),2)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(3)=XTEMP(3)+RMAP(IPUNT(I4,I3,I2,I1),3)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(4)=XTEMP(4)+RMAP(IPUNT(I4,I3,I2,I1),4)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
C.............................CHECKS FOR PARTICLE LOSS
          IF (ABS(XTEMP(1))+ABS(XTEMP(2))+ABS(XTEMP(3))+ABS(XTEMP(4))
     .       .GT.1D10) THEN
            MTRACK4=-ITURN
C            WRITE(6,*) '***MSG(MTRACK4): LOST AT ',ITURN,'/',IL
            RETURN
          ENDIF
C...............................................................
          DO ICOMP=1,4                  !..REINITIALIZATION
            X1(ICOMP)=XTEMP(ICOMP)
            XTEMP(ICOMP)=0D0
            DO IPOT=1,IORD
              XPOT(ICOMP,IPOT)=0D0
            ENDDO
          ENDDO
C...............................................................
        ENDDO
        MTRACK4=NTURN
C...............................................................
      ELSEIF (N.EQ.-1) THEN            !..FIRST JACOBIAN ROW
C...............................................................
        DO ICOMP=1,4                  !..INITIALIZATION
           XTEMP(ICOMP)=0D0
           DO IPOT=1,IORD
             XPOT(ICOMP,IPOT)=0D0
           ENDDO
        ENDDO
C...............................................................
        DO ITURN=1,NTURN
          XPOT(1,0)=1D0
          XPOT(2,0)=1D0
          XPOT(3,0)=1D0
          XPOT(4,0)=1D0
          DO JORD=1,IORD-1
            XPOT(1,JORD)=XPOT(1,JORD-1)*X1(1)
            XPOT(2,JORD)=XPOT(2,JORD-1)*X1(2)
            XPOT(3,JORD)=XPOT(3,JORD-1)*X1(3)
            XPOT(4,JORD)=XPOT(4,JORD-1)*X1(4)
            DO I1=0,JORD
              DO I2=0,JORD-I1
                DO I3=0,JORD-I1-I2
                  I4=JORD-I1-I2-I3
                  XTEMP(1)=XTEMP(1)+RDER(IPUNT(I4,I3,I2,I1),1)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(2)=XTEMP(2)+RDER(IPUNT(I4,I3,I2,I1),2)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(3)=XTEMP(3)+RDER(IPUNT(I4,I3,I2,I1),3)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(4)=XTEMP(4)+RDER(IPUNT(I4,I3,I2,I1),4)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
C.............................CHECKS FOR PARTICLE LOSS
          IF (ABS(XTEMP(1))+ABS(XTEMP(2))+ABS(XTEMP(3))+ABS(XTEMP(4))
     .       .GT.1D10) THEN
            MTRACK4=-ITURN
C            WRITE(6,*) '***MSG(MTRACK4): LOST AT ',ITURN,'/',IL
            RETURN
          ENDIF
C...............................................................
          DO ICOMP=1,4                  !..REINITIALIZATION
            X1(ICOMP)=XTEMP(ICOMP)
            XTEMP(ICOMP)=0D0
            DO IPOT=1,IORD
              XPOT(ICOMP,IPOT)=0D0
            ENDDO
          ENDDO
C...............................................................
        ENDDO
        MTRACK4=NTURN
C...............................................................
      ELSEIF (N.EQ.-2) THEN            !..SECOND JACOBIAN ROW
C...............................................................
        DO ICOMP=1,4                  !..INITIALIZATION
           XTEMP(ICOMP)=0D0
           DO IPOT=1,IORD
             XPOT(ICOMP,IPOT)=0D0
           ENDDO
        ENDDO
C...............................................................
        DO ITURN=1,NTURN
          XPOT(1,0)=1D0
          XPOT(2,0)=1D0
          XPOT(3,0)=1D0
          XPOT(4,0)=1D0
          DO JORD=1,IORD-1
            XPOT(1,JORD)=XPOT(1,JORD-1)*X1(1)
            XPOT(2,JORD)=XPOT(2,JORD-1)*X1(2)
            XPOT(3,JORD)=XPOT(3,JORD-1)*X1(3)
            XPOT(4,JORD)=XPOT(4,JORD-1)*X1(4)
            DO I1=0,JORD
              DO I2=0,JORD-I1
                DO I3=0,JORD-I1-I2
                  I4=JORD-I1-I2-I3
                  XTEMP(1)=XTEMP(1)+RDER(IPUNT(I4,I3,I2,I1),5)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(2)=XTEMP(2)+RDER(IPUNT(I4,I3,I2,I1),6)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(3)=XTEMP(3)+RDER(IPUNT(I4,I3,I2,I1),7)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(4)=XTEMP(4)+RDER(IPUNT(I4,I3,I2,I1),8)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
C.............................CHECKS FOR PARTICLE LOSS
          IF (ABS(XTEMP(1))+ABS(XTEMP(2))+ABS(XTEMP(3))+ABS(XTEMP(4))
     .       .GT.1D10) THEN
            MTRACK4=-ITURN
C            WRITE(6,*) '***MSG(MTRACK4): LOST AT ',ITURN,'/',IL
            RETURN
          ENDIF
C...............................................................
          DO ICOMP=1,4                  !..REINITIALIZATION
            X1(ICOMP)=XTEMP(ICOMP)
            XTEMP(ICOMP)=0D0
            DO IPOT=1,IORD
              XPOT(ICOMP,IPOT)=0D0
            ENDDO
          ENDDO
C...............................................................
        ENDDO
        MTRACK4=NTURN
C...............................................................
      ELSEIF (N.EQ.-3) THEN            !..THIRD JACOBIAN ROW
C...............................................................
        DO ICOMP=1,4                  !..INITIALIZATION
           XTEMP(ICOMP)=0D0
           DO IPOT=1,IORD
             XPOT(ICOMP,IPOT)=0D0
           ENDDO
        ENDDO
C...............................................................
        DO ITURN=1,NTURN
          XPOT(1,0)=1D0
          XPOT(2,0)=1D0
          XPOT(3,0)=1D0
          XPOT(4,0)=1D0
          DO JORD=1,IORD-1
            XPOT(1,JORD)=XPOT(1,JORD-1)*X1(1)
            XPOT(2,JORD)=XPOT(2,JORD-1)*X1(2)
            XPOT(3,JORD)=XPOT(3,JORD-1)*X1(3)
            XPOT(4,JORD)=XPOT(4,JORD-1)*X1(4)
            DO I1=0,JORD
              DO I2=0,JORD-I1
                DO I3=0,JORD-I1-I2
                  I4=JORD-I1-I2-I3
                  XTEMP(1)=XTEMP(1)+RDER(IPUNT(I4,I3,I2,I1),9)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(2)=XTEMP(2)+RDER(IPUNT(I4,I3,I2,I1),10)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(3)=XTEMP(3)+RDER(IPUNT(I4,I3,I2,I1),11)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(4)=XTEMP(4)+RDER(IPUNT(I4,I3,I2,I1),12)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
C.............................CHECKS FOR PARTICLE LOSS
          IF (ABS(XTEMP(1))+ABS(XTEMP(2))+ABS(XTEMP(3))+ABS(XTEMP(4))
     .       .GT.1D10) THEN
            MTRACK4=-ITURN
C            WRITE(6,*) '***MSG(MTRACK4): LOST AT ',ITURN,'/',IL
            RETURN
          ENDIF
C...............................................................
          DO ICOMP=1,4                  !..REINITIALIZATION
            X1(ICOMP)=XTEMP(ICOMP)
            XTEMP(ICOMP)=0D0
            DO IPOT=1,IORD
              XPOT(ICOMP,IPOT)=0D0
            ENDDO
          ENDDO
C...............................................................
        ENDDO
        MTRACK4=NTURN
C...............................................................
      ELSEIF (N.EQ.-4) THEN            !..FOURTH JACOBIAN ROW
C...............................................................
        DO ICOMP=1,4                  !..INITIALIZATION
           XTEMP(ICOMP)=0D0
           DO IPOT=1,IORD
             XPOT(ICOMP,IPOT)=0D0
           ENDDO
        ENDDO
C...............................................................
        DO ITURN=1,NTURN
          XPOT(1,0)=1D0
          XPOT(2,0)=1D0
          XPOT(3,0)=1D0
          XPOT(4,0)=1D0
          DO JORD=1,IORD-1
            XPOT(1,JORD)=XPOT(1,JORD-1)*X1(1)
            XPOT(2,JORD)=XPOT(2,JORD-1)*X1(2)
            XPOT(3,JORD)=XPOT(3,JORD-1)*X1(3)
            XPOT(4,JORD)=XPOT(4,JORD-1)*X1(4)
            DO I1=0,JORD
              DO I2=0,JORD-I1
                DO I3=0,JORD-I1-I2
                  I4=JORD-I1-I2-I3
                  XTEMP(1)=XTEMP(1)+RDER(IPUNT(I4,I3,I2,I1),13)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(2)=XTEMP(2)+RDER(IPUNT(I4,I3,I2,I1),14)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(3)=XTEMP(3)+RDER(IPUNT(I4,I3,I2,I1),15)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                  XTEMP(4)=XTEMP(4)+RDER(IPUNT(I4,I3,I2,I1),16)*
     .                   XPOT(1,I4)*XPOT(2,I3)*XPOT(3,I2)*XPOT(4,I1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
C.............................CHECKS FOR PARTICLE LOSS
          IF (ABS(XTEMP(1))+ABS(XTEMP(2))+ABS(XTEMP(3))+ABS(XTEMP(4))
     .       .GT.1D10) THEN
            MTRACK4=-ITURN
C            WRITE(6,*) '***MSG(MTRACK4): LOST AT ',ITURN,'/',IL
            RETURN
          ENDIF
C...............................................................
          DO ICOMP=1,4                  !..REINITIALIZATION
            X1(ICOMP)=XTEMP(ICOMP)
            XTEMP(ICOMP)=0D0
            DO IPOT=1,IORD
              XPOT(ICOMP,IPOT)=0D0
            ENDDO
          ENDDO
C...............................................................
        ENDDO
        MTRACK4=NTURN
C...............................................................
      ENDIF
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MPHYTONOR.
C=================================================================
C CONVERTS A TRUNCATED TRANSFER MAP EXPRESSED IN PHYSICAL
C COORDINATES INTO COURANT-SNYDER COORDINATES.
C IN INPUT A IS AN ARRAY CONTAINING THE TRANSFER MAP.
C IN OUTPUT IT CONTAINS THE MAP IN NORMALIZED COORDINATES.
C IDIME IS THE PHASE-SPACE DIMENSION.
C IORD IS THE TRUNCATION ORDER OF THE MAP.
C THE COURANT-SNYDER TRANSFORMATION IS STORED IN THE COMMON
C CSTRANSF.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MPHYTONOR(A,IORD,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION CSMAT(NDIM,4),CSMATIN(NDIM,4)
      DIMENSION A(NDIM,4),TEMP(NDIM,4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/CSTRANS/CSM(4,4),CSMIN(4,4)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MPHYTONOR): SECOND PARAMETER OUT OF ',
     .             'BOUNDS'
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MPHYTONOR): THIRD PARAMETER OUT OF ',
     .             'BOUNDS'
        STOP
      ENDIF
C...............................................................
      CALL MZERO(CSMAT,IDIME) !..INITIALIZES THE ARRAYS CSMAT, CSMATIN
      CALL MZERO(CSMATIN,IDIME)
C...............................................................
      IF (IDIME.EQ.2) THEN
C...........................................C-S TRANSFORMATION
        CSMAT(IPUNT(1,0,0,0),1)=CSM(1,1)
        CSMAT(IPUNT(0,1,0,0),1)=CSM(1,2)
        CSMAT(IPUNT(1,0,0,0),2)=CSM(2,1)
        CSMAT(IPUNT(0,1,0,0),2)=CSM(2,2)
C...................................INVERSE C-S TRANSFORMATION
        CSMATIN(IPUNT(1,0,0,0),1)=CSMATIN(1,1)
        CSMATIN(IPUNT(0,1,0,0),1)=CSMATIN(1,2)
        CSMATIN(IPUNT(1,0,0,0),2)=CSMATIN(2,1)
        CSMATIN(IPUNT(0,1,0,0),2)=CSMATIN(2,2)
C........................CSMAT,CSMATIN READY
        CALL MZERO(TEMP,IDIME+2)   !..INITIALIZES TEMP
        CALL COMPOSE(A,CSMAT,TEMP,IORD,1,IORD) !..COMPOSITION
        CALL MZERO(A,IDIME+2)      !..REINITIALIZES A
        CALL COMPOSE(CSMATIN,TEMP,A,1,IORD,IORD) !..COMPOSITION
C...............................................................
      ELSEIF (IDIME.EQ.4) THEN
C...............................................................
        DO ICOMP=1,IDIME
          IND=0
          DO I1=0,1
            DO I2=0,1-I1
              DO I3=0,1-I1-I2
                I4=1-I1-I2-I3
                IND=IND+1
                CSMAT(IPUNT(I4,I3,I2,I1),ICOMP)=CSM(ICOMP,IND)
                CSMATIN(IPUNT(I4,I3,I2,I1),ICOMP)=CSMIN(ICOMP,IND)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C........................CSMAT,CSMATIN READY
        CALL MZERO(TEMP,IDIME)     !..INITIALIZES TEMP
        CALL COMPOSE(A,CSMATIN,TEMP,IORD,1,IORD) !..COMPOSITION
        CALL MZERO(A,IDIME)        !..REINITIALIZES A
        CALL COMPOSE(CSMAT,TEMP,A,1,IORD,IORD) !..COMPOSITION
C...............................................................
      ENDIF
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MNORTOPHY.
C=================================================================
C CONVERTS A TRUNCATED TRANSFER MAP EXPRESSED IN COURANT-SNYDER
C COORDINATES INTO PHYSICAL COORDINATES.
C IN INPUT A IS AN ARRAY CONTAINING THE TRANSFER MAP.
C IN OUTPUT IT CONTAINS THE MAP IN PHYSICAL COORDINATES.
C IDIME IS THE PHASE-SPACE DIMENSION.
C IORD IS THE TRUNCATION ORDER OF THE MAP.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MNORTOPHY(A,IORD,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION CSMAT(NDIM,4),CSMATIN(NDIM,4)
      DIMENSION A(NDIM,4),TEMP(NDIM,4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/CSTRANS/CSM(4,4),CSMIN(4,4)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MNORTOPHY): SECOND PARAMETER OUT OF ',
     .             'BOUNDS'
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MNORTOPHY): THIRD PARAMETER OUT OF ',
     .             'BOUNDS'
        STOP
      ENDIF
C...............................................................
      CALL MZERO(CSMAT,IDIME)  !..INITIALIZES THE ARRAYS CSMAT, CSMATIN
      CALL MZERO(CSMATIN,IDIME)
C...............................................................
      IF (IDIME.EQ.2) THEN
C...........................................C-S TRANSFORMATION
        CSMAT(IPUNT(1,0,0,0),1)=CSM(1,1)
        CSMAT(IPUNT(0,1,0,0),1)=CSM(1,2)
        CSMAT(IPUNT(1,0,0,0),2)=CSM(2,1)
        CSMAT(IPUNT(0,1,0,0),2)=CSM(2,2)
C...................................INVERSE C-S TRANSFORMATION
        CSMATIN(IPUNT(1,0,0,0),1)=CSMATIN(1,1)
        CSMATIN(IPUNT(0,1,0,0),1)=CSMATIN(1,2)
        CSMATIN(IPUNT(1,0,0,0),2)=CSMATIN(2,1)
        CSMATIN(IPUNT(0,1,0,0),2)=CSMATIN(2,2)
C........................CSMAT,CSMATIN READY
        CALL MZERO(TEMP,IDIME+2)   !..INITIALIZES TEMP
        CALL COMPOSE(A,CSMATIN,TEMP,IORD,1,IORD) !..COMPOSITION
        CALL MZERO(A,IDIME+2)      !..REINITIALIZES A
        CALL COMPOSE(CSMAT,TEMP,A,1,IORD,IORD) !..COMPOSITION
C...............................................................
      ELSEIF (IDIME.EQ.4) THEN
C...............................................................
        DO ICOMP=1,IDIME
          IND=0
          DO I1=0,1
            DO I2=0,1-I1
              DO I3=0,1-I1-I2
                I4=1-I1-I2-I3
                IND=IND+1
                CSMAT(IPUNT(I4,I3,I2,I1),ICOMP)=CSM(ICOMP,IND)
                CSMATIN(IPUNT(I4,I3,I2,I1),ICOMP)=CSMIN(ICOMP,IND)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C........................CSMAT,CSMATIN READY
        CALL MZERO(TEMP,IDIME)     !..INITIALIZES TEMP
        CALL COMPOSE(A,CSMAT,TEMP,IORD,1,IORD) !..COMPOSITION
        CALL MZERO(A,IDIME)        !..REINITIALIZES A
        CALL COMPOSE(CSMATIN,TEMP,A,1,IORD,IORD) !..COMPOSITION
C...............................................................
      ENDIF
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, MINVLAT.
C=================================================================
C SUBROUTINE TO INVERT THE LINEAR BLOCS OF AN ACCELERATOR STRUCTURE.
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MINVLAT(IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION TEMP1(4,4),TEMP2(4,4)
      CHARACTER*10 ALBL
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
      COMMON/INVLAT/RLINV(4,4,NLUMP)
      COMMON/JACOBIAN/RJAC(4,4)
C...............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MINVLAT): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IDIME.EQ.2) THEN    !..2D INVERSION
C.............................................................
        DO IBLOC=1,NLMAX+1
          DET=RL(1,1,IBLOC)*RL(2,2,IBLOC)-
     .        RL(1,2,IBLOC)*RL(2,1,IBLOC)
          DETINV=1D0/DET
          RLINV(1,1,IBLOC)=DETINV*RL(2,2,IBLOC)
          RLINV(1,2,IBLOC)=-DETINV*RL(1,2,IBLOC)
          RLINV(2,1,IBLOC)=-DETINV*RL(2,1,IBLOC)
          RLINV(2,2,IBLOC)=DETINV*RL(1,1,IBLOC)
        ENDDO
C.............................................................
      ELSEIF (IDIME.EQ.4) THEN
        DO IBLOC=1,NLMAX+1
C.............................................................
          CALL MCOPY(IBLOC,TEMP2)
          CALL ZEROM(TEMP1)
C.............................................................
          CALL VISYMP(TEMP1,TEMP2)       !..INVERSION
C.............................................................
          DO IROW=1,4
            DO ICOL=1,4
              RLINV(IROW,ICOL,IBLOC)=TEMP1(IROW,ICOL)
            ENDDO
          ENDDO
C.............................................................
        ENDDO
C.............................................................
      ENDIF
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MANYTURNI.
C=============================================================
C MAPS THE VECTOR X1=(X,X',Y,Y') THROUGH THE INVERSE LUMPED
C LATTICE AS GIVEN IN COMMON BLOCK MINVLAT FOR NTURN TURNS.
C NLFUN IS AN EXTERNALLY DECLARED SUBROUTINE WITH THE FORMAT
C "SUBROUTINE NLFUNI(XIN,XOUT,IL,NDIM)" WITH XIN(4) AND XOUT(4).
C IL IS THE LUMP NUMBER AT WHOSE END THE NON-LINEAR BLOCK
C BESIDES. NLFUNI IS SUPPOSED TO HAVE ACCESS TO INL IN VZLUMP
C AND THE COMMON BLOCK LATICE AND ALSO KNOW HOW TO MAP THE
C PARTICLE THROUGH THE SECTION DEFINED BY INL(1,IL) AND
C INL(2,IL). IF THE PARTICLE WITH INITIAL COORDINATES X1(4)
C SURVIVES NTURN TURNS MANYTURNI RETURNS NTURN. IF THE PARTICLE
C IS LOST IT RETURNS THE NEGATIVE C OF THE TURN NUMBER WHEN IT
C WAS LOST.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      INTEGER FUNCTION MANYTURNI(X1,NTURN,NDIM,NLFUNI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000)
      PARAMETER (NLUMP=NELE/10)
      DIMENSION X1(4),X2(4)
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
      COMMON/INVLAT/RLINV(4,4,NLUMP)
      EXTERNAL NLFUNI
C...............................................................
      IF (NDIM.NE.2.AND.NDIM.NE.4) THEN
        WRITE(6,*) '***ERROR(MANYTURN): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C..................................................LOOP OVER TURNS
      DO 1 ITURN=1,NTURN
C.................................................................
        DO I=1,NDIM
          X2(I)=RK(I,NLMAX+1)
          DO J=1,NDIM
            X2(I)=X2(I)+RLINV(I,J,NLMAX+1)*X1(J)
          ENDDO
        ENDDO
C.........................................BACKWARD LOOP OVER LUMPS
        DO 2 IL=NLMAX,1,-1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C....................................................NON-LINEAR MAP
          CALL NLFUNI(X2,X1,IL,NDIM)
C............................NOW X1 CONTAINS THE ACTUAL COORDINATES
C............................................TEST FOR PARTICLE LOSS
          IF (ABS(X1(1))+ABS(X1(2))+ABS(X1(3))+ABS(X1(4)).GT.1D6) THEN
            MANYTURNI=-ITURN
C           WRITE(6,*) '***MSG(MANYTURNI): LOST AT ',ITURN,'/',IL
            GOTO 9999   !..EXIT
          ENDIF
C.........................................LINEAR TRANSPORT OF LUMP
          DO I=1,NDIM
            X2(I)=RK(I,IL)
            DO J=1,NDIM
              X2(I)=X2(I)+RLINV(I,J,IL)*X1(J)
            ENDDO
          ENDDO
    2   CONTINUE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C..................................................................
        DO I=1,NDIM
          X1(I)=X2(I)
        ENDDO
    1 CONTINUE
      MANYTURNI=NTURN
C
 9999 CONTINUE
      RETURN
      END
CDECK  ID>, MJAC.
C=================================================================
C SUBROUTINE TO COMPUTE THE JACOBIAN OF AN ACCELERATOR STRUCTURE.
C XC CONTAINS IN INPUT THE COORDINATES OF THE POINT USED TO COMPUTE
C THE JACOBIAN, IT HAS TO BE DECLARED XC(4) IN THE CALLING PROGRAM.
C IN OUTPUT IT CONTAINS THE ROWS OF THE JACOBIAN
C ACCORDING TO THE VALUES OF N.
C IDIME IS THE PHASE SPACE DIMENSION.
C N IS A FLAG USED TO OUTPUT THE JACOBIAN ROWS:
C N =  0 FIRST ROW
C N = -1 SECOND ROW
C N = -2 THIRD ROW (IDIME=4)
C AND SO ON.
C IFLAG IS USED TO SPECIFIY THE MAP:
C IFLAG =  1 THE EXACT DIRECT MAP IS USED.
C IFLAG = -1 THE EXACT INVERSE MAP IS USED (NOT ACTIVATED).
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MJAC(XC,IDIME,N,IFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION XC(4),X1(4),X2(4),TEMP1(4,4),TEMP2(4,4),TEMP3(4,4)
      CHARACTER*10 ALBL
      COMPLEX*16 UNITY,VAR,GRAD,SUM,POT
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
      COMMON/INVLAT/RLINV(4,4,NLUMP)
      COMMON/JACOBIAN/RJAC(4,4)
C...............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MJAC): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (N.GT.0.OR.N.LT.-(IDIME-1)) THEN
        WRITE(6,*) '***ERROR(MJAC): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (IFLAG.NE.1.AND.IFLAG.NE.-1) THEN
        WRITE(6,*) '***ERROR(MJAC): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IFLAG.EQ.1) THEN  !..DIRECT MAP COMPUTATIONS
C.............................................................
        IF (N.EQ.0) THEN     !..JACOBIAN CALCULATIONS AND 1 ROW
C.............................................................
          UNITY=DCMPLX(0D0,1D0)
C...............................................................
          CALL ZEROM(TEMP1)    !..SET MATRIX TO ZERO
          CALL ZEROM(TEMP2)    !..SET MATRIX TO ZERO
          CALL ZEROM(TEMP3)    !..SET MATRIX TO ZERO
          CALL UNITM(RJAC)     !..SET MATRIX TO ONE
C.............................................................
          DO IL=1,NLMAX
C...............................................................
            MINORD=INT(A(INL(1,IL),3))
            MAXORD=INT(A(INL(2,IL),3))
C.............................................................
            DO I=1,IDIME          !..TRACKING THROUGH LIN. BLOC.
              X1(I)=RK(I,IL)
              DO J=1,IDIME
                X1(I)=X1(I)+RL(I,J,IL)*XC(J)
              ENDDO
            ENDDO
            CALL MCOPY(IL,TEMP1)  !..LIN. PART. JACOBIAN
            CALL MMUL(TEMP2,TEMP1,RJAC,IDIME)
            CALL VMCOPY(RJAC,TEMP2)
C.............................................................
            CALL UNITM(TEMP1)     !..NON-LIN. PART. JACOBIAN
            IF (IA(INL(1,IL)).EQ.2002) THEN       !..SEXTUPOLE
C.............................................................
              GRAD=DCMPLX(A(INL(2,IL),1),-A(INL(2,IL),2))
              VAR=DCMPLX(X1(1),X1(3))
C.............................................................
              TEMP1(2,1)=-DREAL(2D0*GRAD*VAR)
              TEMP1(2,3)=-DREAL(2D0*GRAD*VAR*UNITY)
              TEMP1(4,1)=DIMAG(2D0*GRAD*VAR)
              TEMP1(4,3)=DIMAG(2D0*GRAD*VAR*UNITY)
C.............................................................
            ELSEIF (IA(INL(1,IL)).EQ.2003) THEN   !..OCTUPOLE
C.............................................................
              GRAD=DCMPLX(A(INL(2,IL),1),-A(INL(2,IL),2))
              VAR=DCMPLX(X1(1),X1(3))
C.............................................................
              TEMP1(2,1)=-DREAL(3D0*GRAD*VAR*VAR)
              TEMP1(2,3)=-DREAL(3D0*GRAD*VAR*VAR*UNITY)
              TEMP1(4,1)=DIMAG(3D0*GRAD*VAR*VAR)
              TEMP1(4,3)=DIMAG(3D0*GRAD*VAR*VAR*UNITY)
C.............................................................
            ELSEIF (IA(INL(1,IL)).EQ.2000) THEN   !..MULTIPOLE
C.............................................................
              SUM=DCMPLX(A(INL(2,IL),1),-A(INL(2,IL),2))*
     .            INT(A(INL(2,IL),3))
              VAR=DCMPLX(X1(1),X1(3))
              LOWER1=MAX(MINORD,1)
              LOWER2=INL(1,IL)
              IF (MINORD.EQ.0) LOWER2=LOWER2+1
C.............................................................
              DO JCOMP=INL(2,IL)-1,LOWER2,-1
                IND=INT(A(JCOMP,3))
                SUM=SUM*VAR+IND*DCMPLX(A(JCOMP,1),-A(JCOMP,2))
              ENDDO
C...............................................................
              POT=DCMPLX(1D0,0D0)
              DO JPOT=1,LOWER1-1
                POT=POT*VAR
              ENDDO
              SUM=SUM*POT
C.............................................................
              TEMP1(2,1)=-DREAL(SUM)
              TEMP1(2,3)=-DREAL(SUM*UNITY)
              TEMP1(4,1)=DIMAG(SUM)
              TEMP1(4,3)=DIMAG(SUM*UNITY)
C.............................................................
            ENDIF
C.............................................................
            CALL MMUL(TEMP2,TEMP1,RJAC,IDIME)
            CALL VMCOPY(RJAC,TEMP2)
C.............................................................
            CALL NLFUN(X1,X2,IL,IDIME)  !..TRACKING THROUGH N.L. EL.
            DO I=1,IDIME
              XC(I)=X2(I)          !..CLOSE ITERATION
            ENDDO
C.............................................................
          ENDDO
C.............................................................
          CALL MCOPY(NLMAX+1,TEMP1)    !..TREATS FINAL SECTION
          CALL MMUL(TEMP2,TEMP1,RJAC,IDIME)
          CALL VMCOPY(RJAC,TEMP2)
C.............................................................
          DO JEL=1,IDIME        !..OUTPUT
            XC(JEL)=RJAC(1,JEL)
          ENDDO
C.............................................................
        ELSEIF (N.EQ.-1) THEN !..2 ROW
C.............................................................
          DO JEL=1,IDIME        !..OUTPUT
            XC(JEL)=RJAC(2,JEL)
          ENDDO
C.............................................................
        ELSEIF (N.EQ.-2.AND.IDIME.EQ.4) THEN !..3 ROW
C.............................................................
          DO JEL=1,IDIME        !..OUTPUT
            XC(JEL)=RJAC(3,JEL)
          ENDDO
C.............................................................
        ELSEIF (N.EQ.-3.AND.IDIME.EQ.4) THEN !..4 ROW
C.............................................................
          DO JEL=1,IDIME        !..OUTPUT
            XC(JEL)=RJAC(4,JEL)
          ENDDO
C.............................................................
        ENDIF
C.............................................................
      ELSE                  !..INVERSE MAP COMPUTATIONS
C.............................................................
C        IF(N.EQ.0) THEN     !..JACOBIAN CALCULATIONS AND 1 ROW
C.............................................................
C.............................................................
C        ELSEIF (N.EQ.-1) THEN !..2 ROW
C.............................................................
C          DO JEL=1,IDIME        !..OUTPUT
C            XC(JEL)=RJAC(2,JEL)
C          ENDDO
C.............................................................
C        ELSEIF (N.EQ.-2.AND.IDIME.EQ.4) THEN !..3 ROW
C.............................................................
C          DO JEL=1,IDIME        !..OUTPUT
C            XC(JEL)=RJAC(3,JEL)
C          ENDDO
C.............................................................
C        ELSEIF (N.EQ.-3.AND.IDIME.EQ.4) THEN !..4 ROW
C.............................................................
C          DO JEL=1,IDIME        !..OUTPUT
C            XC(JEL)=RJAC(4,JEL)
C          ENDDO
C.............................................................
C        ENDIF
C.............................................................
      ENDIF
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MMUL.
C======================================================================
C MATRIX MULTIPLICATION OF IDIME BY IDIME  MATRICES, R3=R2*R1
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MMUL(R3,R2,R1,IDIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  R1(4,4),R2(4,4),R3(4,4)
C...............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MMUL): FOURTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C...............................................................
      DO 1 I=1,IDIME
        DO 1 J=1,IDIME
          R3(I,J)=0D0
          DO 1 K=1,IDIME
            R3(I,J)=R3(I,J)+R2(I,K)*R1(K,J)
    1 CONTINUE
C
      RETURN
      END
CDECK  ID>, MCOPY.
C=================================================================
C SUBROUTINE TO COPY THE MATRIX RL(4,4,IL) INTO R.
C IL IS THE INDEX IN THE VZLUMP COMMON.
C R IS THE NEW MATRIX.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MCOPY(IL,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
       DIMENSION R(4,4)
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
C.............................................................
      DO JCOL=1,4
        DO JROW=1,4
          R(JROW,JCOL)=RL(JROW,JCOL,IL)
        ENDDO
      ENDDO
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MGIOMAD.
C=================================================================
C SUBROUTINE TO READ A STRUCTURE FILE PRODUCED BY MAD INSIDE GIOTTO.
C FNAME IS THE FILE NAME WITH THE ACCELERATOR STRUCTURE.
C IORD IS THE ORDER OF THE MAP TRUNCATED MAP.
C ARRAY IS THE ARRAY CONTAINING THE COEFFICIENTS OF THE MAP.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MGIOMAD(FNAME,IORD,ARRAY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      CHARACTER*(*) FNAME
      CHARACTER*10 ALBL
      CHARACTER*80 TITLE
      DIMENSION ARRAY(*),R(4,4)
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
      COMMON/INVLAT/RLINV(4,4,NLUMP)
      COMMON/ACCMAT/RA(4,4,0:NELE),IMAX2,IUPD,IUSE
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MGIOMAD): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IDIME=2                !..DIMENSION OF THE PHASE SPACE
C.............................................................
      IUNIT=1
      OPEN(IUNIT,FILE=FNAME,STATUS='UNKNOWN',FORM='FORMATTED')
C.............................................................
      CALL MLATIN(FNAME,TITLE,IDIME)    !..READS LATTICE FILE
C.............................................................
      IMAX2=IMAX
      CALL VFRESH2(R)         !..ACCUMULATES THE MATRICES
C.............................................................
      NLMAX2=0
      LIMIT=2000
      IKICK=107
      CALL LUMPIT(NLMAX2,LIMIT,IKICK)   !..LUMPS THE LATTICE
      CALL LUMPINPS                     !..C-S COORDINATES
      CALL MINVLAT(IDIME)               !..INVERTS THE LATTICE
C.............................................................
      CALL MRESCALE                     !..RESCALES THE MAP
C.............................................................
      IF (IORD.NE.0) THEN
C.............................................................
        IF (IORD.GT.NORMAX) THEN
          IORD=NORMAX
          WRITE(6,*) '***ORDER GREATER THAN MAXIMUM ALLOWED'
          WRITE(6,*) '***NEW ORDER = ',NORMAX
        ENDIF
C.............................................................
        CALL INIZ_PUNT(NORMAX)      !..INITIALIZES THE POINTER
        CALL MTRANSF(IORD,IDIME)      !..COMPUTES TRANSFER MAP
        CALL MDER(RMAP,IDIME,IORD)      !..DERIVATIVE
        CALL MCONVARR(ARRAY,IORD)       !..CONVERTS THE ARRAYS
      ENDIF
C.............................................................
      CLOSE(IUNIT)
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MITER.
C=================================================================
C THIS FUNCTION PERFORMS THE TRACKING FOR AN ACCELERATOR STRUCTURE.
C XC IN INPUT IS THE STARTING VALUE OF THE COORDINATE AND IT HAS TO
C BE DECLARED XC(4) IN THE CALLING PROGRAM. IN OUTPUT IT IS THE
C FINAL VALUE.
C IMAP IS AN INTEGER SPECIFYING THE MAP TYPE:
C IMAP=0  ===>  EXACT MAP
C IMAP=1  ===>  TRUNCATED MAP
C N REPRESENTS THE NUMBER OF ITERATIONS. THERE ARE THE FOLLOWING
C POSSIBILITIES:
C N >   0 THEN MITER COMPUTES (MAP**N)(XC)
C N =   0 THEN MITER COMPUTES THE FIRST ROW OF THE JACOBIAN
C N =  -1 THEN MITER COMPUTES THE SECOND ROW OF THE JACOBIAN
C N <= -2 THEN MITER COMPUTES (MAP**(N+1))(XC)
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      INTEGER FUNCTION MITER(XC,IMAP,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION XC(4)
      CHARACTER*10 ALBL
      COMMON/LATICE/IA(NELE),A(NELE,5),ALBL(NELE),IMAX
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
      COMMON/VZLUMP/RL(4,4,NLUMP),RK(4,NLUMP),INL(2,NLUMP),NLMAX
      COMMON/INVLAT/RLINV(4,4,NLUMP)
      EXTERNAL NLFUN,NLFUNI
C...............................................................
      IF (IMAP.LT.0) THEN
        WRITE(6,*) '***ERROR(MITER): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      NTURN=1
      IDIME=2
      IFLAG=1
C.............................................................
      IF (N.GT.0) THEN                         !..DIRECT MAP
C.............................................................
        IF (IMAP.EQ.0) THEN
          ITURN=MANYTURN(XC,NTURN,IDIME,NLFUN) !..EXACT MAP
        ELSE
          ITURN=MTRACK2(XC,NTURN,IMAP,N) !..TRUNCATED MAP
        ENDIF
C.............................................................
      ELSEIF (N.EQ.0) THEN                    !..FIRST JAC
C.............................................................
        IF (IMAP.EQ.0) THEN
          CALL MJAC(XC,IDIME,N,IFLAG)       !..EXACT MAP
        ELSE
          ITURN=MTRACK2(XC,NTURN,IMAP,N)    !..TRUNCATED MAP
        ENDIF
C.............................................................
      ELSEIF (N.EQ.-1) THEN                   !..SECOND JAC
C.............................................................
        IF (IMAP.EQ.0) THEN
          CALL MJAC(XC,IDIME,N,IFLAG)             !..EXACT MAP
        ELSE
          ITURN=MTRACK2(XC,NTURN,IMAP,N)    !..TRUNCATED MAP
        ENDIF
C.............................................................
      ELSEIF (N.LE.-2) THEN    !..INVERSE FOR EXACT MAP ONLY
C.............................................................
        IF(IMAP.EQ.0) ITURN=MANYTURNI(XC,NTURN,IDIME,NLFUNI)
C.............................................................
      ENDIF
C.............................................................
      MITER=1
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MGIOSIX.
C=================================================================
C SUBROUTINE TO READ A MAP PRODUCED BY SIXTRACK INSIDE GIOTTO.
C FNAME IS THE FILE NAME WITH THE COEFFICIENTS.
C IORD IS THE ORDER OF THE MAP.
C ARRAY IS THE ARRAY CONTAINING THE COEFFICIENTS OF THE MAP.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MGIOSIX(FNAME,IORD,ARRAY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      CHARACTER*(*) FNAME
      DIMENSION ARRAY(*)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MGIOSOX): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IORD=2                 !..STARTING ORDER OF THE MAP
      IDIME=2                !..DIMENSION OF THE PHASE SPACE
C.............................................................
      CALL INIZ_PUNT(NORMAX)   !..INITIALIZES THE POINTER
C.............................................................
      IUNIT=1
      OPEN(IUNIT,FILE=FNAME,STATUS='UNKNOWN',FORM='FORMATTED')
C.............................................................
      CALL MREAD(RMAP,IORD,IDIME,IUNIT)  !..READS THE MAP
      CALL MPHYTONOR(RMAP,IORD,IDIME)    !..C-S COORDINATES
      CALL MDER(RMAP,IDIME,IORD)         !..DERIVATIVE
      CALL MCONVARR(ARRAY,IORD)          !..CONVERTS THE ARRAYS
C.............................................................
      CLOSE(IUNIT)
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, MCONVARR.
C=================================================================
C SUBROUTINE TO CONVERT AN ARRAY A(NDIM,4) IN A ONE DIMENSIONAL
C ARRAY.
C ARRAY IS THE ONE DIMENSIONAL ARRAY.
C IORD IS THE MAP ORDER.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MCONVARR(ARRAY,IORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION ARRAY(*)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
C...............................................................
      IF (IORD.LT.0.OR.IORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(MCONVARR): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      ICOUNT=0               !..INITIALIZES COUNTER
C.............................................................
      DO JORD=1,IORD         !..LOOP ON THE MAP ORDER
        DO I1=0,JORD
          I2=JORD-I1
          ICOUNT=ICOUNT+1
          ARRAY(ICOUNT)=RMAP(IPUNT(I2,I1,0,0),1) !..COMPONENT 1
          WRITE(20,*) ARRAY(ICOUNT)
          ICOUNT=ICOUNT+1
          ARRAY(ICOUNT)=RMAP(IPUNT(I2,I1,0,0),2) !..COMPONENT 2
          WRITE(20,*) ARRAY(ICOUNT)
        ENDDO
      ENDDO
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, COMPOSEK.
C=================================================================
C COMPUTES THE COMPOSITION OF TWO POLYNOMIALS.
C C=Q(B)
C
C AUTHOR A. BAZZANI - BOLOGNA UNIVERSITY
C
 
      SUBROUTINE COMPOSEK(Q, B, C, NQ, NB, MAXD)
C...............................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION Q(NDIM,4),B(NDIM,4),C(NDIM,4)
C...............................................................
C.....NQ AND NB REPRESENT THE DEGREES OF THE POLYNOMIALS Q E B
      INTEGER NQ,NB
C...............................................................
      CALL COMPONIK(Q(1,1),Q(1,2),Q(1,3),Q(1,4),B(1,1),B(1,2),
     .             B(1,3),B(1,4),C(1,1),C(1,2),C(1,3),C(1,4),
     .             NQ,NB,MAXD)
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, COMPONIK.
C=================================================================
C COMPUTES THE COMPOSITION OF TWO POLYNOMIALS.
C C=Q(B)
C
C AUTHOR A. BAZZANI - BOLOGNA UNIVERSITY
C
 
      SUBROUTINE COMPONIK(Q1, Q2, Q3, Q4, B1, B2, B3, B4,
     .    C1, C2, C3, C4, NQ, NB, MAXD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION Q1(NDIM),Q2(NDIM),Q3(NDIM),Q4(NDIM)
      DIMENSION B1(NDIM),B2(NDIM),B3(NDIM),B4(NDIM)
      DIMENSION C1(NDIM),C2(NDIM),C3(NDIM),C4(NDIM)
C...............................................................
C     NQ ED NB CONTENGONO RISPETTIVAMENTE I GRADI DEI POLINOMI Q E B
C     SI CONSIDERANO I B TUTTI DELLO STESSO GRADO PER NON INTERFERIRE
C     CON LA VETTORIZZAZIONE
C
      INTEGER NQ,NB
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/SCRATCH11/A1(NDIM),A2(NDIM),A3(NDIM),A4(NDIM)
      COMMON/SCRATCH22/D1(NDIM)
C...............................................................
      NCOPPIE(0)=0
      MAXQ=NQ
      MAXP=NB
C     LE SEGUENTI DEFINIZIONI SONO NECESSARIE PER EVITARE WARNING DI
C     COMPILAZIONE, IN QUANTO GLI ARRAYS Q1, Q3 NON SONO USATI IN
C     QUESTA VERSIONE DELLA ROUTINE.
      PIPPO=Q1(1)
      PLUTO=Q3(1)
C
C     CALCOLA GLI INDIRIZZI PER IL PRODOTTO P*Q DI GRADO MAXP E MAXQ
C     TENENDO I TERMINI FINO AL GRADO MAXD
C
      CALL INIZPRODO4(MAXD,MAXP,MAXD)
C
C     DA` IL NUMERO TOTALE DI TERMINI IN UN POLINOMIO DI GRADO MAXD
C
      NMAX=(MAXD+4)*(MAXD+3)*(MAXD+2)*(MAXD+1)/24
C
C AZZERAMENTO INIZIALE COMPLETO DELL'ARRAY RISULTATO:
C
      DO  1  I1=1, NMAX
      A1 (I1)=0.D0
      A2 (I1)=0.D0
      A3 (I1)=0.D0
      A4 (I1)=0.D0
      C1 (I1)=0.D0
      C2 (I1)=0.D0
      C3 (I1)=0.D0
      C4 (I1)=0.D0
  1   D1 (I1)=0.D0
C
C L'ARRAY "A1" SCRATCH E` INIZIALIZZATO ALLA COSTANTE 1.
C
      D1(1) = 1.D0
      A1(1) = 1.D0
      A2(1) = 1.D0
      A3(1) = 1.D0
      A4(1) = 1.D0
C
C INGRESSO NEL QUADRUPLO LOOP CHE NUMERA I CONTRIBUTI ALLA COMPOSIZIONE:
C
      DO  2  J1=0, MAXQ
C
C INGRESSO NEL SECONDO PASSO DEL QUADRUPLO LOOP:
C
        J2=0
C      DO  3  J2=0, 0
C
C INGRESSO NEL TERZO PASSO DEL QUADRUPLO LOOP:
C
      DO  4  J3=0, MAXQ-J1-J2
C
C INGRESSO NEL QUARTO PASSO DEL QUADRUPLO LOOP:
C
         J4=0
C      DO  5  J4=0, 0
C
C NELL'ARRAY A4 E' CONTENUTO IL CONTRIBUTO DEL MONOMIO DI INDICE
C "J1,J2,J3,J4".
C ORA TUTTI GLI ELEMENTI DELL'ARRAY "A4" SONO MOLTIPLICATI PER IL COEFFICIENTE
C DI INDICI "J1,J2,J3,J4" DELL'ARRAY COMPONENDO "Q1,Q2,Q3,Q4"
C E ACCUMULATI NELL'ARRAY RISULTATO
C
C      VAR=Q1(IPUNT(J1,J2,J3,J4))
C      IF (DABS(VAR).GE.1.E-14) THEN
C         DO   I1=1, NMAX
C              C1(I1)=C1(I1)+
C     .               A4(I1)*VAR
C         ENDDO
C      ENDIF
      VAR=Q2(IPUNT(J1,J2,J3,J4))
      IF (DABS(VAR).GE.1.E-14) THEN
         DO   I1=1, NMAX
              C2(I1)=C2(I1)+
     .               A4(I1)*VAR
         ENDDO
      ENDIF
C      VAR=Q3(IPUNT(J1,J2,J3,J4))
C      IF (DABS(VAR).GE.1.E-14) THEN
C         DO   I1=1, NMAX
C              C3(I1)=C3(I1)+
C     .               A4(I1)*VAR
C         ENDDO
C      ENDIF
      VAR=Q4(IPUNT(J1,J2,J3,J4))
      IF (DABS(VAR).GE.1.E-14) THEN
         DO  I1=1, NMAX
              C4(I1)=C4(I1)+
     .               A4(I1)*VAR
         ENDDO
      ENDIF
C
C L'ARRAY SCRATCH "A4",PROVENIENTE DALL'ESTERNO DEL LOOP,E` MOLTIPLICATO
C PER L'ARRAY "B4", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
C
C      IF (J4.LT.MAXQ-J3-J2-J1) THEN
C           CALL PRODO(A4, B4, D1, MAXD)
C
C DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "D1", E`
C RIMESSO IN "A4" CHE RISULTA COSI` AGGIORNATO AL VALORE
C B1**J1*B2**J2*B3**J3*B4**(J4+1)
C
C           DO  59  I1=1, NMAX
C 59        A4(I1) = D1(I1)
C      ENDIF
C 5    CONTINUE
C
C L'ARRAY SCRATCH "A3", PROVENIENTE DALL'ESTERNO DEL LOOP, E` MOLTIPLICATO
C PER L'ARRAY "B3", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
C
      IF (J3.LT.MAXQ-J1-J2) THEN
           CALL PRODO(A3, B3, A4, MAXD)
C
C DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "A4", E`
C RIMESSO IN "A3" CHE RISULTA COSI` AGGIORNATO AL VALORE
C B1**J1*B2**J2*B3**(J3+1)
C
          DO  29  I1=1, NMAX
 29       A3(I1) = A4(I1)
      ENDIF
C
 4    CONTINUE
C
C L'ARRAY SCRATCH "A2" ,PROVENIENTE DALL'ESTERNO DEL LOOP,E` MOLTIPLICATO
C PER L'ARRAY "B2", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
C
C      IF (J2.LT.MAXQ-J1) THEN
C           CALL PRODO(A2, B2, A4, MAXD)
C
C DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "D1", E`
C RIMESSO IN "A2" CHE RISULTA COSI` AGGIORNATO AL VALORE B1**J1*B2**(J2+1)
C
C L'ARRAY "A3" SCRATCH E` INIZIALIZZATO AL VALORE CORRENTE DELL'ARRAY "A2":
C
C           DO  26  I1=1, NMAX
C           A3(I1) = A4(I1)
C 26        A2(I1) = A3(I1)
C      ENDIF
C
C 3    CONTINUE
C
C L'ARRAY SCRATCH "A1" ,PROVENIENTE DALL'ESTERNO DEL LOOP,E` MOLTIPLICATO
C PER L'ARRAY "B1", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
C
      IF (J1.LT.MAXQ) THEN
           CALL PRODO(A1, B1, A4, MAXD)
C
C DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "D1", E`
C RIMESSO IN "A1" CHE RISULTA COSI` AGGIORNATO AL VALORE B1**J1
C
C L'ARRAY "A2" SCRATCH E` INIZIALIZZATO AL VALORE CORRENTE DELL'ARRAY "A1":
C
           DO  23  I1=1, NMAX
           A3(I1) = A4(I1)
           A2(I1) = A3(I1)
 23        A1(I1) = A2(I1)
      ENDIF
C
  2   CONTINUE
C...............................................................
      DO 24 I1=1,NMAX
        C1(I1)=C1(I1)+B1(I1)
        C2(I1)=C2(I1)+B2(I1)
        C3(I1)=C3(I1)+B3(I1)
        C4(I1)=C4(I1)+B4(I1)
 24   CONTINUE
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, COMPOSE.
C=================================================================
C COMPUTES THE COMPOSITION OF TWO POLYNOMIALS.
C C=Q(B)
C
C AUTHOR A. BAZZANI - BOLOGNA UNIVERSITY
C
 
      SUBROUTINE COMPOSE(Q, B, C, NQ, NB, MAXD)
C...............................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION Q(NDIM,4),B(NDIM,4),C(NDIM,4)
C...............................................................
C.....NQ AND NB REPRESENT THE DEGREES OF THE POLYNOMIALS Q E B
      INTEGER NQ,NB
C...............................................................
      CALL COMPONI(Q(1,1),Q(1,2),Q(1,3),Q(1,4),B(1,1),B(1,2),
     .             B(1,3),B(1,4),C(1,1),C(1,2),C(1,3),C(1,4),
     .             NQ,NB,MAXD)
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, COMPONI.
C=================================================================
C COMPUTES THE COMPOSITION OF TWO POLYNOMIALS.
C C=Q(B)
C
C AUTHOR A. BAZZANI - BOLOGNA UNIVERSITY
C
 
      SUBROUTINE COMPONI(Q1, Q2, Q3, Q4, B1, B2, B3, B4,
     .    C1, C2, C3, C4, NQ, NB, MAXD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION Q1(NDIM),Q2(NDIM),Q3(NDIM),Q4(NDIM)
      DIMENSION B1(NDIM),B2(NDIM),B3(NDIM),B4(NDIM)
      DIMENSION C1(NDIM),C2(NDIM),C3(NDIM),C4(NDIM)
C...............................................................
C     NQ ED NB CONTENGONO RISPETTIVAMENTE I GRADI DEI POLINOMI Q E B
C     SI CONSIDERANO I B TUTTI DELLO STESSO GRADO PER NON INTERFERIRE
C     CON LA VETTORIZZAZIONE
C
      INTEGER NQ,NB
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/SCRATCH11/A1(NDIM),A2(NDIM),A3(NDIM),A4(NDIM)
      COMMON/SCRATCH22/D1(NDIM)
C...............................................................
      NCOPPIE(0)=0
      MAXQ=NQ
      MAXP=NB
C
C     CALCOLA GLI INDIRIZZI PER IL PRODOTTO P*Q DI GRADO MAXP E MAXQ
C     TENENDO I TERMINI FINO AL GRADO MAXD
C
      CALL INIZPRODO4(MAXD,MAXP,MAXD)
C
C     DA` IL NUMERO TOTALE DI TERMINI IN UN POLINOMIO DI GRADO MAXD
C
      NMAX=(MAXD+4)*(MAXD+3)*(MAXD+2)*(MAXD+1)/24
C
C AZZERAMENTO INIZIALE COMPLETO DELL'ARRAY RISULTATO:
C
      DO  1  I1=1, NMAX
      A1 (I1)=0.D0
      A2 (I1)=0.D0
      A3 (I1)=0.D0
      A4 (I1)=0.D0
      C1 (I1)=0.D0
      C2 (I1)=0.D0
      C3 (I1)=0.D0
      C4 (I1)=0.D0
  1   D1 (I1)=0.D0
C
C L'ARRAY "A1" SCRATCH E` INIZIALIZZATO ALLA COSTANTE 1.
C
      D1(1) = 1.D0
      A1(1) = 1.D0
      A2(1) = 1.D0
      A3(1) = 1.D0
      A4(1) = 1.D0
C
C INGRESSO NEL QUADRUPLO LOOP CHE NUMERA I CONTRIBUTI ALLA COMPOSIZIONE:
C
      DO  2  J1=0, MAXQ
C
C INGRESSO NEL SECONDO PASSO DEL QUADRUPLO LOOP:
C
      DO  3  J2=0, MAXQ-J1
C
C INGRESSO NEL TERZO PASSO DEL QUADRUPLO LOOP:
C
      DO  4  J3=0, MAXQ-J1-J2
C
C INGRESSO NEL QUARTO PASSO DEL QUADRUPLO LOOP:
C
      DO  5  J4=0, MAXQ-J1-J2-J3
C
C NELL'ARRAY A4 E' CONTENUTO IL CONTRIBUTO DEL MONOMIO DI INDICE
C "J1,J2,J3,J4".
C ORA TUTTI GLI ELEMENTI DELL'ARRAY "A4" SONO MOLTIPLICATI PER IL COEFFICIENTE
C DI INDICI "J1,J2,J3,J4" DELL'ARRAY COMPONENDO "Q1,Q2,Q3,Q4"
C E ACCUMULATI NELL'ARRAY RISULTATO
C
      VAR=Q1(IPUNT(J1,J2,J3,J4))
      IF (DABS(VAR).GE.1.E-14) THEN
         DO   I1=1, NMAX
              C1(I1)=C1(I1)+
     .               A4(I1)*VAR
         ENDDO
      ENDIF
      VAR=Q2(IPUNT(J1,J2,J3,J4))
      IF (DABS(VAR).GE.1.E-14) THEN
         DO   I1=1, NMAX
              C2(I1)=C2(I1)+
     .               A4(I1)*VAR
         ENDDO
      ENDIF
      VAR=Q3(IPUNT(J1,J2,J3,J4))
      IF (DABS(VAR).GE.1.E-14) THEN
         DO   I1=1, NMAX
              C3(I1)=C3(I1)+
     .               A4(I1)*VAR
         ENDDO
      ENDIF
      VAR=Q4(IPUNT(J1,J2,J3,J4))
      IF (DABS(VAR).GE.1.E-14) THEN
         DO  I1=1, NMAX
              C4(I1)=C4(I1)+
     .               A4(I1)*VAR
         ENDDO
      ENDIF
C
C L'ARRAY SCRATCH "A4",PROVENIENTE DALL'ESTERNO DEL LOOP,E` MOLTIPLICATO
C PER L'ARRAY "B4", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
C
      IF (J4.LT.MAXQ-J3-J2-J1) THEN
           CALL PRODO(A4, B4, D1, MAXD)
C
C DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "D1", E`
C RIMESSO IN "A4" CHE RISULTA COSI` AGGIORNATO AL VALORE
C B1**J1*B2**J2*B3**J3*B4**(J4+1)
C
           DO  59  I1=1, NMAX
 59        A4(I1) = D1(I1)
      ENDIF
 5    CONTINUE
C
C L'ARRAY SCRATCH "A3", PROVENIENTE DALL'ESTERNO DEL LOOP, E` MOLTIPLICATO
C PER L'ARRAY "B3", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
C
      IF (J3.LT.MAXQ-J1-J2) THEN
           CALL PRODO(A3, B3, A4, MAXD)
C
C DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "A4", E`
C RIMESSO IN "A3" CHE RISULTA COSI` AGGIORNATO AL VALORE
C B1**J1*B2**J2*B3**(J3+1)
C
          DO  29  I1=1, NMAX
 29       A3(I1) = A4(I1)
      ENDIF
C
 4    CONTINUE
C
C L'ARRAY SCRATCH "A2" ,PROVENIENTE DALL'ESTERNO DEL LOOP,E` MOLTIPLICATO
C PER L'ARRAY "B2", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
C
      IF (J2.LT.MAXQ-J1) THEN
           CALL PRODO(A2, B2, A4, MAXD)
C
C DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "D1", E`
C RIMESSO IN "A2" CHE RISULTA COSI` AGGIORNATO AL VALORE B1**J1*B2**(J2+1)
C
C L'ARRAY "A3" SCRATCH E` INIZIALIZZATO AL VALORE CORRENTE DELL'ARRAY "A2":
C
           DO  26  I1=1, NMAX
           A3(I1) = A4(I1)
 26        A2(I1) = A3(I1)
      ENDIF
C
 3    CONTINUE
C
C L'ARRAY SCRATCH "A1" ,PROVENIENTE DALL'ESTERNO DEL LOOP,E` MOLTIPLICATO
C PER L'ARRAY "B1", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
C
      IF (J1.LT.MAXQ) THEN
           CALL PRODO(A1, B1, A4, MAXD)
C
C DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "D1", E`
C RIMESSO IN "A1" CHE RISULTA COSI` AGGIORNATO AL VALORE B1**J1
C
C L'ARRAY "A2" SCRATCH E` INIZIALIZZATO AL VALORE CORRENTE DELL'ARRAY "A1":
C
           DO  23  I1=1, NMAX
           A3(I1) = A4(I1)
           A2(I1) = A3(I1)
 23        A1(I1) = A2(I1)
      ENDIF
C
  2   CONTINUE
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, PRODO.
C=================================================================
C COMPUTES THE COMPOSITION OF TWO POLYNOMIALS.
C C=Q(B)
C
C AUTHOR A. BAZZANI - BOLOGNA UNIVERSITY
C
 
      SUBROUTINE PRODO(Q, P, R, MAXD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION P(NDIM),Q(NDIM),R(NDIM)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      MAXNT=IPUNT(MAXD,0,0,0)
      DO 1 J=1, MAXNT
1      R(J)=0.D0
      DO 2 J=1, MAXNT
      DO 2 NJ=NCOPPIE(J-1)+1, NCOPPIE(J)
2           R(J) = R(J) + P(INDP(NJ))*Q(INDQ(NJ))
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, INIZPRODO4.
C=================================================================
C COMPUTES THE COMPOSITION OF TWO POLYNOMIALS.
C C=Q(B)
C
C AUTHOR A. BAZZANI - BOLOGNA UNIVERSITY
C
 
      SUBROUTINE INIZPRODO4(MAXQ,MAXP,MAXD)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      NCOPPIE(0)=0
      LCONT=0
      MAXP=MIN0(MAXD,MAXP)
      MAXQ=MIN0(MAXD,MAXQ)
      IP=0
         DO 2 J=0, MAXD
         DO 2 K1=0, J
         DO 2 K2=0, J-K1
         DO 2 K3=0, J-K1-K2
            K4=J-K1-K2-K3
            IP=IP+1
            DO 3 J1=MAX0(0,K1-MAXQ), MIN0(MAXP,K1)
                 I1=K1-J1
            DO 3 J2=MAX0(0,K2-MAXQ+I1), MIN0(MAXP-J1,K2)
                 I2=K2-J2
            DO 3 J3=MAX0(0,K3-MAXQ+I1+I2), MIN0(MAXP-J1-J2,K3)
                 I3=K3-J3
            DO 3 J4=MAX0(0,K4-MAXQ+I1+I2+I3), MIN0(MAXP-J1-J2-J3,K4)
                 I4=K4-J4
                LCONT=LCONT+1
                INDP(LCONT)=IPUNT(J1,J2,J3,J4)
                INDQ(LCONT)=IPUNT(I1,I2,I3,I4)
3            CONTINUE
      NCOPPIE(IP)=LCONT
2     CONTINUE
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, INIZ_PUNT.
C=================================================================
C COMPUTES THE COMPOSITION OF TWO POLYNOMIALS.
C C=Q(B)
C
C AUTHOR A. BAZZANI - BOLOGNA UNIVERSITY
C
 
      SUBROUTINE INIZ_PUNT(MAXD)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
C...............................................................
      IF (MAXD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(INIZ_PUNT): FIRST PARAMETER OUT OF ',
     .             'BOUNDS.'
        STOP
      ENDIF
C...............................................................
      ICONT=0
      DO 1 NORD=0, MAXD+1
         DO 1 N1=0, NORD
            DO 1 N2=0, NORD-N1
               DO 1 N3=0, NORD-N1-N2
                  N4=NORD-N1-N2-N3
                  ICONT=ICONT+1
1                 IPUNT(N1,N2,N3,N4)=ICONT
      DO NVAR=1, 4
         DO ICONT=1, IPUNT(MAXD,0,0,0)
            IDERIV(ICONT,NVAR)=0
         ENDDO
      ENDDO
      ICONT=0
      DO 2 NORD=0, MAXD+1
         DO 2 N1=0, NORD
            DO 2 N2=0, NORD-N1
               DO 2 N3=0, NORD-N1-N2
                  N4=NORD-N1-N2-N3
                  ICONT=ICONT+1
              IF (N1.GE.1) IDERIV(ICONT,1)=IPUNT(N1-1,N2,N3,N4)
              IF (N2.GE.1) IDERIV(ICONT,2)=IPUNT(N1,N2-1,N3,N4)
              IF (N3.GE.1) IDERIV(ICONT,3)=IPUNT(N1,N2,N3-1,N4)
2             IF (N4.GE.1) IDERIV(ICONT,4)=IPUNT(N1,N2,N3,N4-1)
C...............................................................
        RETURN
C...............................................................
        END
CDECK  ID>, COMPOSEL.
C=================================================================
C COMPUTES THE COMPOSITION OF TWO POLYNOMIALS.
C C=Q(B)
C
C AUTHOR A. BAZZANI - BOLOGNA UNIVERSITY
C
 
      SUBROUTINE COMPOSEL(Q, B, C, NQ, NB, MAXD)
C...............................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION Q(NDIM,4),B(NDIM,4),C(NDIM,4)
C...............................................................
C.....NQ AND NB REPRESENT THE DEGREES OF THE POLYNOMIALS Q E B
      INTEGER NQ,NB
C...............................................................
      CALL COMPONIL(Q(1,1),Q(1,2),Q(1,3),Q(1,4),B(1,1),B(1,2),
     .             B(1,3),B(1,4),C(1,1),C(1,2),C(1,3),C(1,4),
     .             NQ,NB,MAXD)
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, COMPONIL.
C=================================================================
C COMPUTES THE COMPOSITION OF TWO POLYNOMIALS.
C C=Q(B)
C
C AUTHOR A. BAZZANI - BOLOGNA UNIVERSITY
C
 
      SUBROUTINE COMPONIL(Q1, Q2, Q3, Q4, B1, B2, B3, B4,
     .    C1, C2, C3, C4, NQ, NB, MAXD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION Q1(NDIM),Q2(NDIM),Q3(NDIM),Q4(NDIM)
      DIMENSION B1(NDIM),B2(NDIM),B3(NDIM),B4(NDIM)
      DIMENSION C1(NDIM),C2(NDIM),C3(NDIM),C4(NDIM)
C...............................................................
C     NQ ED NB CONTENGONO RISPETTIVAMENTE I GRADI DEI POLINOMI Q E B
C     SI CONSIDERANO I B TUTTI DELLO STESSO GRADO PER NON INTERFERIRE
C     CON LA VETTORIZZAZIONE
C
      INTEGER NQ,NB
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/SCRATCH11/A1(NDIM),A2(NDIM),A3(NDIM),A4(NDIM)
      COMMON/SCRATCH22/D1(NDIM)
C...............................................................
      NCOPPIE(0)=0
      MAXQ=NQ
      MAXP=NB
C
C     CALCOLA GLI INDIRIZZI PER IL PRODOTTO P*Q DI GRADO MAXP E MAXQ
C     TENENDO I TERMINI FINO AL GRADO MAXD
C
C      CALL INIZPRODO4(MAXD,MAXP,MAXD)
C
C     DA` IL NUMERO TOTALE DI TERMINI IN UN POLINOMIO DI GRADO MAXD
C
      NMAX=(MAXD+4)*(MAXD+3)*(MAXD+2)*(MAXD+1)/24
C
C AZZERAMENTO INIZIALE COMPLETO DELL'ARRAY RISULTATO:
C
      DO  1  I1=1, NMAX
      A1 (I1)=0.D0
      A2 (I1)=0.D0
      A3 (I1)=0.D0
      A4 (I1)=0.D0
      C1 (I1)=0.D0
      C2 (I1)=0.D0
      C3 (I1)=0.D0
      C4 (I1)=0.D0
  1   D1 (I1)=0.D0
C
C L'ARRAY "A1" SCRATCH E` INIZIALIZZATO ALLA COSTANTE 1.
C
      D1(1) = 1.D0
      A1(1) = 1.D0
      A2(1) = 1.D0
      A3(1) = 1.D0
      A4(1) = 1.D0
C
C INGRESSO NEL QUADRUPLO LOOP CHE NUMERA I CONTRIBUTI ALLA COMPOSIZIONE:
C
      DO  2  J1=0, MAXQ
C
C INGRESSO NEL SECONDO PASSO DEL QUADRUPLO LOOP:
C
      DO  3  J2=0, MAXQ-J1
C
C INGRESSO NEL TERZO PASSO DEL QUADRUPLO LOOP:
C
      DO  4  J3=0, MAXQ-J1-J2
C
C INGRESSO NEL QUARTO PASSO DEL QUADRUPLO LOOP:
C
      DO  5  J4=0, MAXQ-J1-J2-J3
C
C NELL'ARRAY A4 E' CONTENUTO IL CONTRIBUTO DEL MONOMIO DI INDICE
C "J1,J2,J3,J4".
C ORA TUTTI GLI ELEMENTI DELL'ARRAY "A4" SONO MOLTIPLICATI PER IL COEFFICIENTE
C DI INDICI "J1,J2,J3,J4" DELL'ARRAY COMPONENDO "Q1,Q2,Q3,Q4"
C E ACCUMULATI NELL'ARRAY RISULTATO
C
      VAR=Q1(IPUNT(J1,J2,J3,J4))
      IF (DABS(VAR).GE.1.E-14) THEN
         DO   I1=1, NMAX
              C1(I1)=C1(I1)+
     .               A4(I1)*VAR
         ENDDO
      ENDIF
      VAR=Q2(IPUNT(J1,J2,J3,J4))
      IF (DABS(VAR).GE.1.E-14) THEN
         DO   I1=1, NMAX
              C2(I1)=C2(I1)+
     .               A4(I1)*VAR
         ENDDO
      ENDIF
      VAR=Q3(IPUNT(J1,J2,J3,J4))
      IF (DABS(VAR).GE.1.E-14) THEN
         DO   I1=1, NMAX
              C3(I1)=C3(I1)+
     .               A4(I1)*VAR
         ENDDO
      ENDIF
      VAR=Q4(IPUNT(J1,J2,J3,J4))
      IF (DABS(VAR).GE.1.E-14) THEN
         DO  I1=1, NMAX
              C4(I1)=C4(I1)+
     .               A4(I1)*VAR
         ENDDO
      ENDIF
C
C L'ARRAY SCRATCH "A4",PROVENIENTE DALL'ESTERNO DEL LOOP,E` MOLTIPLICATO
C PER L'ARRAY "B4", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
C
      IF (J4.LT.MAXQ-J3-J2-J1) THEN
C           CALL PRODO(A4, B4, D1, MAXD)
C
C DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "D1", E`
C RIMESSO IN "A4" CHE RISULTA COSI` AGGIORNATO AL VALORE
C B1**J1*B2**J2*B3**J3*B4**(J4+1)
C
           DO  59  I1=1, NMAX
 59        A4(I1) = B4(I1)
      ENDIF
 5    CONTINUE
C
C L'ARRAY SCRATCH "A3", PROVENIENTE DALL'ESTERNO DEL LOOP, E` MOLTIPLICATO
C PER L'ARRAY "B3", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
C
      IF (J3.LT.MAXQ-J1-J2) THEN
C           CALL PRODO(A3, B3, A4, MAXD)
C
C DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "A4", E`
C RIMESSO IN "A3" CHE RISULTA COSI` AGGIORNATO AL VALORE
C B1**J1*B2**J2*B3**(J3+1)
C
          DO  29  I1=1, NMAX
 29       A4(I1) = B3(I1)
      ENDIF
C
 4    CONTINUE
C
C L'ARRAY SCRATCH "A2" ,PROVENIENTE DALL'ESTERNO DEL LOOP,E` MOLTIPLICATO
C PER L'ARRAY "B2", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
C
      IF (J2.LT.MAXQ-J1) THEN
C           CALL PRODO(A2, B2, A4, MAXD)
C
C DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "D1", E`
C RIMESSO IN "A2" CHE RISULTA COSI` AGGIORNATO AL VALORE B1**J1*B2**(J2+1)
C
C L'ARRAY "A3" SCRATCH E` INIZIALIZZATO AL VALORE CORRENTE DELL'ARRAY "A2":
C
           DO  26  I1=1, NMAX
 26        A4(I1) = B2(I1)
      ENDIF
C
 3    CONTINUE
C
C L'ARRAY SCRATCH "A1" ,PROVENIENTE DALL'ESTERNO DEL LOOP,E` MOLTIPLICATO
C PER L'ARRAY "B1", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
C
      IF (J1.LT.MAXQ) THEN
C           CALL PRODO(A1, B1, A4, MAXD)
C
C DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "D1", E`
C RIMESSO IN "A1" CHE RISULTA COSI` AGGIORNATO AL VALORE B1**J1
C
C L'ARRAY "A2" SCRATCH E` INIZIALIZZATO AL VALORE CORRENTE DELL'ARRAY "A1":
C
           DO  23  I1=1, NMAX
 23        A4(I1) = B1(I1)
      ENDIF
C
  2   CONTINUE
C...............................................................
      RETURN
C...............................................................
      END
CDECK  ID>, SMAPDYN2.
C=================================================================
C COMPUTES STABLE POINTS INSIDE THE DYNAMIC APERTURE AND WRITES
C THEM ON UNIT 1.
C IMAP IS A PARAMETER SPECIFING THE MAP TYPE:
C IMAP = 1 ===> HENON MAP
C IMAP = 2 ===> ACCELERATOR MODEL
C NTURN IS THE TURN NUMBER USED TO COMPUTE THE STABILITY.
C IDIME IS THE PHASE SPACE DIMENSION.
C (XMIN,XMAX) REPRESENTS THE INTERVAL FOR X VALUES. NX INITIAL
C CONDITIONS ARE CONSIDERED.
C (YMIN,YMAX) REPRESENTS THE INTERVAL FOR Y VALUES. NY INITIAL
C CONDITIONS ARE CONSIDERED.
C N.B: IF IDIME = 2 THEN THE INITAL CONDITIONS ARE (X0,Y0), WHILE
C          IF IDIME = 4 THEN THE INITAL CONDITIONS ARE (X0,0,Y0,0).
C ITUNE IS A FLAG FOR TUNE CALCULATIONS:
C ITUNE > 0 THE TUNE IS COMPUTED. IT IS POSSIBLE TO CHOOSE THE
C TECHNIQUE USED FOR THE COMPUTATION OF THE TUNE:
C ITUNE=1 ===> TUNEAPA  IS USED
C ITUNE=2 ===> TUNEFFT  IS USED
C ITUNE=3 ===> TUNEFFTI IS USED
C ITUNE=4 ===> TUNELASK IS USED
C ITUNE=5 ===> TUNENEWT IS USED
C ITUNE=6 ===> TUNEFIT  IS USED
C ITUNE=7 ===> TUNEABT2 IS USED
C ITUNE = 0 ONLY THE STABILITY OF A POINT IS COMPUTED.
C XS,YS ARE ARRAYS CONTAINING THE COORDINATES OF THE STABLE POINTS.
C QX,QY ARE ARRAYS CONTAINING THE TUNES OF THE STABLE POINTS (IF ITUNE
C = 1 IS SPECIFIED).
C IPOINT REPRESENTS THE NUMBER OF STABLE INITIAL CONDITIONS.
C AMPL2 REPRESENTS THE AVERAGE RADIUS CORRESPONDING TO THE
C DISTRIBUTION OF STABLE POINTS.
C
C AUTHOR M. GIOVANNOZZI - CERN & BOLOGNA UNIVERSITY
C
 
      SUBROUTINE MAPDYN2(IMAP,NTURN,IDIME,XMIN,XMAX,NX,YMIN,YMAX,NY,
     .                   ITUNE,XS,YS,QX,QY,IPOINT,AMPL2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      PARAMETER (NMAX=100000,NPOINT=50000)
      DIMENSION X(4),X1(NMAX),XP1(NMAX),Y1(NMAX),YP1(NMAX)
      DIMENSION XS(NPOINT),YS(NPOINT),QX(NPOINT),QY(NPOINT)
C.............................................................
      EXTERNAL NLFUN
C.............................................................
      IF (IMAP.NE.1.AND.IMAP.NE.2) THEN
        WRITE(6,*) '***ERROR(MAPDYN2): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (NX*NY.GT.NPOINT) THEN
        WRITE(6,*) '***ERROR(MAPDYN2): TOO MANY INITIAL CONDITIONS'
        WRITE(6,*) '***MESSAGE(MAPDYN2): THE PRODUCT OF NX * NY',
     .             ' CONDITIONS MUST BE < ',NPOINT
        STOP
      ENDIF
C.............................................................
      IF (NTURN.GT.NMAX.OR.(NTURN+1.GT.NMAX.AND.ITUNE.GE.1)) THEN
        WRITE(6,*) '***ERROR(MAPDYN2): TOO MANY ITERATIONS'
        WRITE(6,*) '***MESSAGE(MAPDYN2): THE NUMBER OF ITERATIONS',
     .             ' MUST BE < ',NMAX
        STOP
      ENDIF
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MAPDYN2): FOURTH PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
      IF (ITUNE.LT.0.OR.ITUNE.GT.7) THEN
        WRITE(6,*) '***ERROR(MAPDYN2): TENTH PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
C.............................................................
      STEPX=(XMAX-XMIN)/DFLOAT(NX)     !..COMPUTES THE STEP IN THE GRID
      IF (STEPX.LT.0D0) THEN
        STEPX=-STEPX
        WRITE(6,*) '***MESSAGE(MAPDYN2): WRONG ORDER IN X INTERVAL'
        WRITE(6,*) '***MESSAGE(MAPDYN2): ORDER INVERTED'
      ENDIF
C.............................................................
      STEPY=(YMAX-YMIN)/DFLOAT(NY)
      IF (STEPY.LT.0D0) THEN
        STEPY=-STEPY
        WRITE(6,*) '***MESSAGE(MAPDYN2): WRONG ORDER IN Y INTERVAL'
        WRITE(6,*) '***MESSAGE(MAPDYN2): ORDER INVERTED'
      ENDIF
C.............................................................
      PI=4D0*DATAN(1D0)
C.............................................................
      IPOINT=0          !..INITIALIZES POINT COUNTER
      ISTABLE=0         !..INITIALIZES STABLE PARTICLE COUNTER
C.............................................................
      IF (ITUNE.EQ.0) THEN
C.............................................................
        DO ICOL=1,NX                   !..LOOP ON X
          XNEW=XMIN+STEPX*DFLOAT(ICOL)
C.............................................................
          DO IROW=1,NY                 !..LOOP ON Y
            YNEW=YMIN+STEPY*DFLOAT(IROW)
            IF (IDIME.EQ.2) THEN
C.............................................................
              X(1)=XNEW
              X(2)=YNEW
              X(3)=0D0
              X(4)=0D0
C.............................................................
            ELSE
C.............................................................
              X(1)=XNEW
              X(2)=0D0
              X(3)=YNEW
              X(4)=0D0
C.............................................................
            ENDIF
C.............................................................
            IF (IMAP.EQ.1) THEN
              ITURN=MANY_HENON(X,NTURN)
            ELSE
              ITURN=MANYTURN(X,NTURN,IDIME,NLFUN)  !..TRACKING
            ENDIF
C.............................................................
            IF (ITURN.EQ.NTURN) THEN
              IPOINT=IPOINT+1
              XS(IPOINT)=XNEW
              YS(IPOINT)=YNEW
              ISTABLE=ISTABLE+1
            ENDIF
C.............................................................
          ENDDO
C.............................................................
        ENDDO
C.............................................................
      ELSE
C.............................................................
        DO ICOL=1,NX                   !..LOOP ON X
          XNEW=XMIN+STEPX*DFLOAT(ICOL)
C.............................................................
          DO IROW=1,NY                 !..LOOP ON Y
            YNEW=YMIN+STEPY*DFLOAT(IROW)
C..SKIPS SOME POINTS TO AVOID PROBLEMS IN THE TUNE COMPUTATION
C...................................ORIGIN FOR EVERY STRUCTURE
            IF (DABS(XNEW)+DABS(YNEW).EQ.0D0) GOTO 10
C................................X AXIS FOR THE HENON MAP ONLY
            IF (IMAP.EQ.1.AND.YNEW.EQ.0D0) GOTO 10
C.............................................................
            IF (IDIME.EQ.2) THEN
C.............................................................
              X(1)=XNEW
              X(2)=YNEW
              X(3)=0D0
              X(4)=0D0
C.............................................................
              DO JITER=1,NTURN
C.............................................................
                IF (IMAP.EQ.1) THEN
                  ITURN=MANY_HENON(X,1)
                ELSE
                  ITURN=MANYTURN(X,1,IDIME,NLFUN)
                ENDIF
C.............................................................
                IF (ITURN.LT.0) GOTO 10
                X1(JITER)=X(1)
                XP1(JITER)=X(2)
              ENDDO
              IPOINT=IPOINT+1            !..NEW STABLE PARTICLE
              XS(IPOINT)=XNEW
              YS(IPOINT)=YNEW
C.............................................................
              IF (ITUNE.EQ.1) THEN
                QX(IPOINT)=TUNEAPA(X1,XP1,NTURN)  !..NEW TUNE
              ELSEIF (ITUNE.EQ.2) THEN
                QX(IPOINT)=TUNEFFT(X1,XP1,NTURN)  !..NEW TUNE
              ELSEIF (ITUNE.EQ.3) THEN
                QX(IPOINT)=TUNEFFTI(X1,XP1,NTURN) !..NEW TUNE
              ELSEIF (ITUNE.EQ.4) THEN
                QX(IPOINT)=TUNELASK(X1,XP1,NTURN) !..NEW TUNE
              ELSEIF (ITUNE.EQ.5) THEN
                QX(IPOINT)=TUNENEWT(X1,XP1,NTURN) !..NEW TUNE
              ELSEIF (ITUNE.EQ.6) THEN
                QX(IPOINT)=TUNEFIT(X1,XP1,NTURN) !..NEW TUNE
              ELSEIF (ITUNE.EQ.7) THEN
                QX(IPOINT)=TUNEABT2(X1,XP1,NTURN) !..NEW TUNE
              ENDIF
C.............................................................
              ISTABLE=ISTABLE+1
C.............................................................
            ELSE
C.............................................................
              X(1)=XNEW
              X(2)=0D0
              X(3)=YNEW
              X(4)=0D0
C.............................................................
              DO JITER=1,NTURN
C.............................................................
                IF (IMAP.EQ.1) THEN
                  ITURN=MANY_HENON(X,1)
                ELSE
                  ITURN=MANYTURN(X,1,IDIME,NLFUN)
                ENDIF
C.............................................................
                IF (ITURN.LT.0) GOTO 10
                X1(JITER)=X(1)
                XP1(JITER)=X(2)
                Y1(JITER)=X(3)
                YP1(JITER)=X(4)
              ENDDO
              IPOINT=IPOINT+1            !..NEW STABLE PARTICLE
              XS(IPOINT)=XNEW
              YS(IPOINT)=YNEW
C.............................................................
              IF (ITUNE.EQ.1) THEN
                QX(IPOINT)=TUNEAPA(X1,XP1,NTURN)  !..NEW TUNE
                QY(IPOINT)=TUNEAPA(Y1,YP1,NTURN)  !..NEW TUNE
              ELSEIF (ITUNE.EQ.2) THEN
                QX(IPOINT)=TUNEFFT(X1,XP1,NTURN)  !..NEW TUNE
                QY(IPOINT)=TUNEFFT(Y1,YP1,NTURN)  !..NEW TUNE
              ELSEIF (ITUNE.EQ.3) THEN
                QX(IPOINT)=TUNEFFTI(X1,XP1,NTURN) !..NEW TUNE
                QY(IPOINT)=TUNEFFTI(Y1,YP1,NTURN) !..NEW TUNE
              ELSEIF (ITUNE.EQ.4) THEN
                QX(IPOINT)=TUNELASK(X1,XP1,NTURN) !..NEW TUNE
                QY(IPOINT)=TUNELASK(Y1,YP1,NTURN) !..NEW TUNE
              ELSEIF (ITUNE.EQ.5) THEN
                QX(IPOINT)=TUNENEWT(X1,XP1,NTURN) !..NEW TUNE
                QY(IPOINT)=TUNENEWT(Y1,YP1,NTURN) !..NEW TUNE
              ELSEIF (ITUNE.EQ.6) THEN
                QX(IPOINT)=TUNEFIT(X1,XP1,NTURN) !..NEW TUNE
                QY(IPOINT)=TUNEFIT(Y1,YP1,NTURN) !..NEW TUNE
              ELSEIF (ITUNE.EQ.7) THEN
                QX(IPOINT)=TUNEABT2(X1,XP1,NTURN) !..NEW TUNE
                QY(IPOINT)=TUNEABT2(Y1,YP1,NTURN) !..NEW TUNE
              ENDIF
C.............................................................
              ISTABLE=ISTABLE+1
C.............................................................
            ENDIF
C.............................................................
 10         CONTINUE
          ENDDO
C.............................................................
        ENDDO
C.............................................................
      ENDIF
C.............................................................
      FRACTION=DFLOAT(ISTABLE)/DFLOAT(NX*NY) !..COMPUTES FRACTION
      AMPL2=SQRT(XMAX*YMAX*FRACTION*4D0/PI)
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, SMAPDYN3.
C============================================================
C COMPUTES STABLE POINTS INSIDE THE DYNAMIC APERTURE
C IMAP IS A PARAMETER SPECIFING THE MAP TYPE:
C IMAP = 1 ===> HENON MAP
C IMAP = 2 ===> ACCELERATOR MODEL
C NTURN IS THE TURN NUMBER USED TO COMPUTE THE STABILITY.
C (RMIN,RMAX) REPRESENTS THE INTERVAL FOR RADIUS VALUES.
C NR INITIAL CONDITIONS ARE CONSIDERED.
C (ANG1,ANG2) REPRESENTS THE INTERVAL FOR THE ANGLES IN UNIT
C OF TWOPI. NANG INITIAL DIRECTIONS ARE CONSIDERED.
C ITUNE IS A FLAG FOR TUNE CALCULATIONS:
C ITUNE > 0 THE TUNE IS COMPUTED. IT IS POSSIBLE TO CHOOSE THE
C TECHNIQUE USED FOR THE COMPUTATION OF THE TUNE:
C ITUNE=1 ===> TUNEAPA  IS USED
C ITUNE=2 ===> TUNEFFT  IS USED
C ITUNE=3 ===> TUNEFFTI IS USED
C ITUNE=4 ===> TUNELASK IS USED
C ITUNE=5 ===> TUNELASK IS USED
C ITUNE=6 ===> TUNEFIT  IS USED
C ITUNE=7 ===> TUNEABT2 IS USED
C ITUNE = 0 ONLY THE STABILITY OF A POINT IS COMPUTED.
C XS,YS ARE ARRAYS CONTAINING THE COORDINATES OF THE STABLE
C PARTICLES AS A FUNCTION OF THE ANGLE.
C XR,YR ARE THE ARRAYS CONTAINING THE COORDINATES OF THE BORDER
C OF THE STABILITY DOMAIN.
C IPOINT REPRESENTS THE NUMBER OF STABLE INITIAL CONDITIONS.
C AMPL IS  THE AVERAGE DISTANCE WHERE PARTICLES ARE LOST .
C
C AUTHOR : R. GRASSI (M. GIOVANNOZZI - CERN - HAS INTRODUCED SOME
C CHANGES).
C
      SUBROUTINE MAPDYN3(IMAP,NTURN,RMIN,RMAX,NR,ANG1,ANG2,NANG,
     .                   ITUNE,XS,YS,XR,YR,QX,QY,IPOINT,AMPL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      PARAMETER (NMAX=100000,NPOINT=50000)
      DIMENSION X(4),X1(NMAX),XP1(NMAX),Y1(NMAX),YP1(NMAX)
      DIMENSION XR(NPOINT),YR(NPOINT)
      DIMENSION XS(NPOINT),YS(NPOINT),QX(NPOINT),QY(NPOINT)
C.............................................................
      EXTERNAL NLFUN
C.............................................................
      IF (IMAP.NE.1.AND.IMAP.NE.2) THEN
        WRITE(6,*) '***ERROR(MAPDYN3): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (NR*NANG.GT.NPOINT) THEN
        WRITE(6,*) '***ERROR(MAPDYN3): TOO MANY INITIAL CONDITIONS'
        STOP
      ENDIF
C.............................................................
      IF (NTURN.GT.NMAX.OR.(NTURN+1.GT.NMAX.AND.ITUNE.GE.1)) THEN
        WRITE(6,*) '***ERROR(MAPDYN3): TOO MANY ITERATIONS'
        WRITE(6,*) '***MESSAGE(MAPDYN3): THE NUMBER OF ITERATIONS',
     .             ' MUST BE < ',NMAX
        STOP
      ENDIF
      IF (ITUNE.LT.0.OR.ITUNE.GT.7) THEN
        WRITE(6,*) '***ERROR(MAPDYN3): ITUNE PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
C.............................................................
      STEPR=(RMAX-RMIN)/DFLOAT(NR)   !..COMPUTES THE STEP IN THE RADIUS S
      IF (STEPR.LT.0D0) THEN
        STEPR=-STEPR
        WRITE(6,*) '***MESSAGE(MAPDYN3): WRONG ORDER IN R INTERVAL'
        WRITE(6,*) '***MESSAGE(MAPDYN3): ORDER INVERTED'
      ENDIF
C.............................................................
      STEPANG=(ANG2-ANG1)/DFLOAT(NANG+1)
      IF (STEPANG.LT.0D0) THEN
        STEPANG=-STEPANG
        WRITE(6,*) '***MESSAGE(MAPDYN3): WRONG ORDER IN ANGLE INTERVAL'
        WRITE(6,*) '***MESSAGE(MAPDYN3): ORDER INVERTED'
      ENDIF
C.............................................................
      TWOPI=8D0*DATAN(1D0)
      IDIME=4
      IPOINT=0          !..INITIALIZES POINT COUNTER
C...............................................................
      IF (ITUNE.EQ.0) THEN
C...............................................................
        DO IROW=1,NANG            !..LOOP ON ANGLE
          ANGNEW=ANG1+STEPANG*DFLOAT(IROW)
C...............................................................
          DO ICOL=1,NR           !..LOOP ON R
            RNEW=RMIN+STEPR*DFLOAT(ICOL)
C.............................................................
            X(1)=RNEW*COS(ANGNEW*TWOPI)
            X(2)=0D0
            X(3)=RNEW*SIN(ANGNEW*TWOPI)
            X(4)=0D0
C.............................................................
            IF (IMAP.EQ.1) THEN
              ITURN=MANY_HENON(X,NTURN)
            ELSE
              ITURN=MANYTURN(X,NTURN,IDIME,NLFUN)
            ENDIF
C.............................................................
            IF (ITURN.NE.NTURN) THEN
              IF (ICOL.EQ.1) THEN
                WRITE(6,*) '***WARNING(MAPDYN3): RMIN UNSTABLE'
              ENDIF
              XR(IROW)=(RNEW-STEPR)*COS(ANGNEW*TWOPI)
              YR(IROW)=(RNEW-STEPR)*SIN(ANGNEW*TWOPI)
              GOTO 10
            ELSE
              IPOINT=IPOINT+1
              XS(IPOINT)=RNEW*COS(ANGNEW*TWOPI)
              YS(IPOINT)=RNEW*SIN(ANGNEW*TWOPI)
            ENDIF
C.............................................................
          ENDDO
          XR(IROW)=RNEW*COS(ANGNEW*TWOPI)
          YR(IROW)=RNEW*SIN(ANGNEW*TWOPI)
          WRITE(6,*) '***WARNING(MAPDYN3): RMAX STABLE'
          WRITE(6,*) '***LAST PARTICLE SAVED!'
 10       CONTINUE
C.............................................................
        ENDDO
C.............................................................
      ELSE
C.............................................................
        DO IROW=1,NANG            !..LOOP ON ANGLE
          ANGNEW=ANG1+STEPANG*DFLOAT(IROW)
C...............................................................
          DO ICOL=1,NR           !..LOOP ON R
            RNEW=RMIN+STEPR*DFLOAT(ICOL)
C.............................................................
            X(1)=RNEW*COS(ANGNEW*TWOPI)
            X(2)=0D0
            X(3)=RNEW*SIN(ANGNEW*TWOPI)
            X(4)=0D0
C.............................................................
            DO JITER=1,NTURN
C.............................................................
              IF (IMAP.EQ.1) THEN
                ITURN=MANY_HENON(X,1)
              ELSE
                ITURN=MANYTURN(X,1,IDIME,NLFUN)
              ENDIF
C.............................................................
              IF (ITURN.LT.0) THEN
                IF (ICOL.EQ.1) THEN
                  WRITE(6,*) '***WARNING(MAPDYN3): RMIN UNSTABLE'
                ENDIF
                XR(IROW)=(RNEW-STEPR)*COS(ANGNEW*TWOPI)
                YR(IROW)=(RNEW-STEPR)*SIN(ANGNEW*TWOPI)
                GOTO 20
              ENDIF
              X1(JITER)=X(1)
              XP1(JITER)=X(2)
              Y1(JITER)=X(3)
              YP1(JITER)=X(4)
            ENDDO
            IPOINT=IPOINT+1            !..NEW STABLE PARTICLE
            XS(IPOINT)=RNEW*COS(ANGNEW*TWOPI)
            YS(IPOINT)=RNEW*SIN(ANGNEW*TWOPI)
C.............................................................
            IF (ITUNE.EQ.1) THEN
              QX(IPOINT)=TUNEAPA(X1,XP1,NTURN)  !..NEW TUNE
              QY(IPOINT)=TUNEAPA(Y1,YP1,NTURN)  !..NEW TUNE
            ELSEIF (ITUNE.EQ.2) THEN
              QX(IPOINT)=TUNEFFT(X1,XP1,NTURN)  !..NEW TUNE
              QY(IPOINT)=TUNEFFT(Y1,YP1,NTURN)  !..NEW TUNE
            ELSEIF (ITUNE.EQ.3) THEN
              QX(IPOINT)=TUNEFFTI(X1,XP1,NTURN) !..NEW TUNE
              QY(IPOINT)=TUNEFFTI(Y1,YP1,NTURN) !..NEW TUNE
            ELSEIF (ITUNE.EQ.4) THEN
              QX(IPOINT)=TUNELASK(X1,XP1,NTURN) !..NEW TUNE
              QY(IPOINT)=TUNELASK(Y1,YP1,NTURN) !..NEW TUNE
            ELSEIF (ITUNE.EQ.5) THEN
              QX(IPOINT)=TUNENEWT(X1,XP1,NTURN) !..NEW TUNE
              QY(IPOINT)=TUNENEWT(Y1,YP1,NTURN) !..NEW TUNE
            ELSEIF (ITUNE.EQ.6) THEN
              QX(IPOINT)=TUNEFIT(X1,XP1,NTURN) !..NEW TUNE
              QY(IPOINT)=TUNEFIT(Y1,YP1,NTURN) !..NEW TUNE
            ELSEIF (ITUNE.EQ.7) THEN
              QX(IPOINT)=TUNEABT2(X1,XP1,NTURN) !..NEW TUNE
              QY(IPOINT)=TUNEABT2(Y1,YP1,NTURN) !..NEW TUNE
            ENDIF
C.............................................................
          ENDDO
          XR(IROW)=RNEW*COS(ANGNEW*TWOPI)
          YR(IROW)=RNEW*SIN(ANGNEW*TWOPI)
          WRITE(6,*) '***WARNING(MAPDYN3): RMAX STABLE'
          WRITE(6,*) '***LAST PARTICLE SAVED!'
 20       CONTINUE
C.............................................................
        ENDDO
C.............................................................
      ENDIF
C.............................................................
      AMPL=0
      DO IROW=1,NANG
        AMPL=AMPL+SQRT(XR(IROW)**2+YR(IROW)**2)
      ENDDO
      AMPL=AMPL/NANG
C.............................................................
      RETURN
      END
CDECK  ID>, SMAPDYN4.
C============================================================
C COMPUTES STABLE POINTS INSIDE THE DYNAMIC APERTURE
C IMAP IS A PARAMETER SPECIFING THE MAP TYPE:
C IMAP = 1 ===> HENON MAP
C IMAP = 2 ===> ACCELERATOR MODEL
C NTURN IS THE TURN NUMBER USED TO COMPUTE THE STABILITY.
C (RMIN,RMAX) REPRESENTS THE INTERVAL FOR RADIUS VALUES.
C NR INITIAL CONDITIONS ARE CONSIDERED.
C (ANG1,ANG2) REPRESENTS THE INTERVAL FOR THE ANGLES IN UNIT
C OF TWOPI. NANG INITIAL DIRECTIONS ARE CONSIDERED.
C ITUNE IS A FLAG FOR TUNE CALCULATIONS:
C ITUNE > 0 THE TUNE IS COMPUTED. IT IS POSSIBLE TO CHOOSE THE
C TECHNIQUE USED FOR THE COMPUTATION OF THE TUNE:
C ITUNE=1 ===> TUNEAPA  IS USED
C ITUNE=2 ===> TUNEAPA + NORMAL FORMS ARE USED
C ITUNE=3 ===> TUNELASK IS USED
C ITUNE=4 ===> TUNENEWT IS USED
C ITUNE=5 ===> TUNEABT2 IS USED
C ITUNE= 0 ONLY THE STABILITY OF A POINT IS COMPUTED.
C IT PRODUCES SOME OUTPUT FILES:
C UNIT 25 ===> INITIAL COORDINATES AND STABILITY
C UNIT 26 ===> INITIAL COORDINATES, STABILITY AND TUNE DIFFERENCES
C              (SEE LONG WRITE UP)
C UNIT 27 ===> INITIAL COORDINATES, STABILITY AND TUNE DISTANCES (FIRST
C              TYPE) (SEE LONG WRITE UP)
C UNIT 28 ===> INITIAL COORDINATES, STABILITY AND TUNE DISTANCES (FIRST
C              TYPE) (SEE LONG WRITE UP)
C UNIT 29 ===> INITIAL COORDINATES, STABILITY AND TUNE DISTANCES (SECOND
C              TYPE) (SEE LONG WRITE UP)
C UNIT 30 ===> INITIAL COORDINATES, STABILITY AND TUNE DISTANCES (SECOND
C              TYPE) (SEE LONG WRITE UP)
C
C AUTHOR : M. GIOVANNOZZI - CERN
C
      SUBROUTINE MAPDYN4(IMAP,NTURN,RMIN,RMAX,NR,ANG1,ANG2,NANG,
     .                   ITUNE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      PARAMETER (NMAX=10000,NPOINT=50000)
      DIMENSION X(4),X1(NMAX),XP1(NMAX),Y1(NMAX),YP1(NMAX)
      DIMENSION X2(NMAX),XP2(NMAX),Y2(NMAX),YP2(NMAX)
      DIMENSION X3(NMAX),XP3(NMAX),Y3(NMAX),YP3(NMAX)
      DIMENSION DQX1(4,NMAX),DQY1(4,NMAX),DQX2(4,NMAX),DQY2(4,NMAX)
      DIMENSION DQ(4,NMAX),RHO(NMAX),ALPHA(NMAX)
      DIMENSION LOSS(NMAX),NT(4)
      COMPLEX*16 Z1(NMAX),Z2(NMAX),ZETA1(NMAX),ZETA2(NMAX)
      COMPLEX*16 AUTORES,ER1,ER2
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
      COMMON/RISONANZA/IRES1(2),IRES2(2),AUTORES(4),IRESON,ICASE
      COMMON/TUNEPAR/ADVSIG,ADVMIN,ADVMAX
C.............................................................
      EXTERNAL NLFUN
C.............................................................
      IF (IMAP.NE.1.AND.IMAP.NE.2) THEN
        WRITE(6,*) '***ERROR(MAPDYN4): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (NTURN.LE.0) THEN
        WRITE(6,*) '***ERROR(MAPDYN4): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (RMIN.LT.0D0.OR.RMAX.LT.0D0) THEN
        WRITE(6,*) '***ERROR(MAPDYN4): THIRD OR FOURTH PARAMETER OUT ',
     .             'OF BOUNDS'
        STOP
      ENDIF
      IF (NR.LT.0) THEN
        WRITE(6,*) '***ERROR(MAPDYN4): FIFTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (ANG1.LT.0D0.OR.ANG2.LT.0D0) THEN
        WRITE(6,*) '***ERROR(MAPDYN4): SIXTH OR SEVENTH PARAMETER OUT ',
     .             'OF BOUNDS'
        STOP
      ENDIF
      IF (NR*NANG.GT.NPOINT) THEN
        WRITE(6,*) '***ERROR(MAPDYN4): TOO MANY INITIAL CONDITIONS'
        STOP
      ENDIF
      IF (ITUNE.LT.0.OR.ITUNE.GT.5) THEN
        WRITE(6,*) '***ERROR(MAPDYN4): ITUNE PARAMETER OUT OF BOUNDS.'
        STOP
      ENDIF
C.............................................................
      STEPR=(RMAX-RMIN)/DFLOAT(NR)   !..COMPUTES THE STEP IN THE RADIUS S
      IF (STEPR.LT.0D0) THEN
        STEPR=-STEPR
        WRITE(6,*) '***MESSAGE(MAPDYN4): WRONG ORDER IN R INTERVAL'
        WRITE(6,*) '***MESSAGE(MAPDYN4): ORDER INVERTED'
      ENDIF
C.............................................................
      STEPANG=(ANG2-ANG1)/DFLOAT(NANG+1)
      IF (STEPANG.LT.0D0) THEN
        STEPANG=-STEPANG
        WRITE(6,*) '***MESSAGE(MAPDYN4): WRONG ORDER IN ANGLE INTERVAL'
        WRITE(6,*) '***MESSAGE(MAPDYN4): ORDER INVERTED'
      ENDIF
C.............................................................
      TWOPI=8D0*DATAN(1D0)
      IDIME=4
      IPOINT=0
C...............................................................
      NT(1)=64
      NT(2)=256
      NT(3)=1024
      NT(4)=4096
C...............................................................
      IF (ITUNE.EQ.2) THEN !..PREPARES EVERYTHING FOR NORMAL FORMS
C...............................................................
        IF (IMAP.EQ.1) THEN
C...............................................................
          CALL INIZ_PUNT(NORMAX)
          MORD=2
          IF (OCTUP.NE.0D0) MORD=3
C...............................................................
        ELSE
C...............................................................
          MORD=5
          CALL MTRANSF(MORD,IDIME)
C...............................................................
        ENDIF
C...............................................................
        IRESON=0
        ICASE=0
        OR1=0.
        OR2=0.
        J11=0
        J12=0
        J21=0
        J22=0
C...............................................................
        OMER1=OR1*TWOPI
        OMER2=OR2*TWOPI
        ER1=EXP(DCMPLX(0D0,OMER1))
        ER2=EXP(DCMPLX(0D0,OMER2))
C...............................................................
        AUTORES(1)=ER1
        AUTORES(2)=DCONJG(ER1)
        AUTORES(3)=ER2
        AUTORES(4)=DCONJG(ER2)
C...............................................................
        IRES1(1)=J11
        IRES2(1)=J12
        IRES1(2)=J21
        IRES2(2)=J22
C...............................................................
        NFORD=5
        INTERACT=0
        IINV=1
        IOUT=0
        IAZAN=0
        CALL GENERALE_4D(RMAP,MORD,NFORD,INTERACT,IINV,IOUT,IAZAN)
      ENDIF
C...............................................................
      IF (ITUNE.EQ.0) THEN
C...............................................................
        DO IROW=1,NANG            !..LOOP ON ANGLE
          ANGNEW=ANG1+STEPANG*DFLOAT(IROW)
C...............................................................
          DO ICOL=1,NR           !..LOOP ON R
            RNEW=RMIN+STEPR*DFLOAT(ICOL)
C.............................................................
            X(1)=RNEW*COS(ANGNEW*TWOPI)
            X(2)=0D0
            X(3)=RNEW*SIN(ANGNEW*TWOPI)
            X(4)=0D0
C.............................................................
            IPOINT=IPOINT+1
            RHO(IPOINT)=RNEW      !..DEFINES INITIAL CONDITIONS
            ALPHA(IPOINT)=ANGNEW
C.............................................................
            IF (IMAP.EQ.1) THEN
              ITURN=MANY_HENON(X,NTURN)
            ELSE
              ITURN=MANYTURN(X,NTURN,IDIME,NLFUN)
            ENDIF
C.............................................................
            IF (ITURN.NE.NTURN) THEN
              LOSS(IPOINT)=-ITURN
            ELSE
              LOSS(IPOINT)=ITURN+1
            ENDIF
C.............................................................
          ENDDO
C.............................................................
        ENDDO
C.............................................................
      ELSE
C.............................................................
        DO IROW=1,NANG            !..LOOP ON ANGLE
          ANGNEW=ANG1+STEPANG*DFLOAT(IROW)
C...............................................................
          DO ICOL=1,NR           !..LOOP ON R
            RNEW=RMIN+STEPR*DFLOAT(ICOL)
C.............................................................
            X(1)=RNEW*COS(ANGNEW*TWOPI)
            X(2)=0D0
            X(3)=RNEW*SIN(ANGNEW*TWOPI)
            X(4)=0D0
C.............................................................
            IPOINT=IPOINT+1
            RHO(IPOINT)=RNEW
            ALPHA(IPOINT)=ANGNEW
C.............................................................
            DO JITER=1,2*NT(4)
C.............................................................
              IF (IMAP.EQ.1) THEN
                ITURN=MANY_HENON(X,1)
              ELSE
                ITURN=MANYTURN(X,1,IDIME,NLFUN)
              ENDIF
C.............................................................
              IF (ITURN.LT.0) THEN
                LOSS(IPOINT)=JITER
                GOTO 20
              ENDIF
C.............................................................
              X1(JITER)=X(1)
              XP1(JITER)=X(2)
              Y1(JITER)=X(3)
              YP1(JITER)=X(4)
            ENDDO
C.............................................................
            LOSS(IPOINT)=2*NT(4)+1
C.............................................................
 20         DO IWIN=1,4
              IF (LOSS(IPOINT).GE.2*NT(IWIN)) THEN      
                IF (ITUNE.EQ.1) THEN      !..TUNE COMPUTATIONS
                  QX1=TUNEAPA(X1,XP1,NT(IWIN))  !..WINDOW1
C.............................................................
                  XMIN1=ADVMIN
                  XMAX1=ADVMAX
                  XSIG1=ADVSIG/DSQRT(DFLOAT(NT(IWIN)))
C.............................................................
                  QY1=TUNEAPA(Y1,YP1,NT(IWIN))
C.............................................................
                  YMIN1=ADVMIN
                  YMAX1=ADVMAX
                  YSIG1=ADVSIG/DSQRT(DFLOAT(NT(IWIN)))
C.............................................................
                ELSEIF (ITUNE.EQ.2) THEN
C.............................................................
                  DO IP=1,NT(IWIN)
                    Z1(IP)=DCMPLX(X1(IP),-XP1(IP))  !..COMPLEX COOR
                    Z2(IP)=DCMPLX(Y1(IP),-YP1(IP))
                    CALL EVAL_PSI(NFORD,Z1(IP),Z2(IP),
     .                            ZETA1(IP),ZETA2(IP))
                    X3(IP)=DREAL(ZETA1(IP))         !..BACK TO REAL
                    XP3(IP)=-DIMAG(ZETA1(IP))
                    Y3(IP)=DREAL(ZETA2(IP))
                    YP3(IP)=-DIMAG(ZETA2(IP))
                  ENDDO
C.............................................................
                  QX1=TUNEAPA(X3,XP3,NT(IWIN))  !..WINDOW1
C.............................................................
                  XMIN1=ADVMIN
                  XMAX1=ADVMAX
                  XSIG1=ADVSIG/DSQRT(DFLOAT(NT(IWIN)))
C.............................................................
                  QY1=TUNEAPA(Y3,YP3,NT(IWIN))
C.............................................................
                  YMIN1=ADVMIN
                  YMAX1=ADVMAX
                  YSIG1=ADVSIG/DSQRT(DFLOAT(NT(IWIN)))
C.............................................................
                ELSEIF (ITUNE.EQ.3) THEN
C.............................................................
                  QX1=TUNELASK(X1,XP1,NT(IWIN)) !..WINDOW1
C.............................................................
                  XMIN1=0D0
                  XMAX1=0D0
                  XSIG1=0D0
C.............................................................
                  QY1=TUNELASK(Y1,YP1,NT(IWIN))
C.............................................................
                  YMIN1=0D0
                  YMAX1=0D0
                  YSIG1=0D0
C.............................................................
                ELSEIF (ITUNE.EQ.4) THEN
C.............................................................
                  QX1=TUNENEWT(X1,XP1,NT(IWIN)) !..WINDOW1
C.............................................................
                  XMIN1=0D0
                  XMAX1=0D0
                  XSIG1=0D0
C.............................................................
                  QY1=TUNENEWT(Y1,YP1,NT(IWIN))
C.............................................................
                  YMIN1=0D0
                  YMAX1=0D0
                  YSIG1=0D0
C.............................................................
                ELSEIF (ITUNE.EQ.5) THEN
C.............................................................
                  QX1=TUNEABT2(X1,XP1,NT(IWIN)) !..WINDOW1
C.............................................................
                  XMIN1=0D0
                  XMAX1=0D0
                  XSIG1=0D0
C.............................................................
                  QY1=TUNEABT2(Y1,YP1,NT(IWIN))
C.............................................................
                  YMIN1=0D0
                  YMAX1=0D0
                  YSIG1=0D0
C.............................................................
                ENDIF
C.............................................................
                CALL MCOPYV(X1,NT(IWIN)+1,2*NT(IWIN),X2)     !..COPIES
                CALL MCOPYV(XP1,NT(IWIN)+1,2*NT(IWIN),XP2)
                CALL MCOPYV(Y1,NT(IWIN)+1,2*NT(IWIN),Y2)
                CALL MCOPYV(YP1,NT(IWIN)+1,2*NT(IWIN),YP2)
C.............................................................
                IF (ITUNE.EQ.1) THEN      !..TUNE COMPUTATIONS
                  QX2=TUNEAPA(X2,XP2,NT(IWIN))  !..WINDOW2
C.............................................................
                  XMIN2=ADVMIN
                  XMAX2=ADVMAX
                  XSIG2=ADVSIG/DSQRT(DFLOAT(NT(IWIN)))
C.............................................................
                  QY2=TUNEAPA(Y2,YP2,NT(IWIN))
C.............................................................
                  YMIN2=ADVMIN
                  YMAX2=ADVMAX
                  YSIG2=ADVSIG/DSQRT(DFLOAT(NT(IWIN)))
C.............................................................
                ELSEIF (ITUNE.EQ.2) THEN
C.............................................................
                  DO IP=1,NT(IWIN)
                    Z1(IP)=DCMPLX(X2(IP),-XP2(IP))  !..COMPLEX COOR
                    Z2(IP)=DCMPLX(Y2(IP),-YP2(IP))
                    CALL EVAL_PSI(NFORD,Z1(IP),Z2(IP),
     .                            ZETA1(IP),ZETA2(IP))
                    X3(IP)=DREAL(ZETA1(IP))         !..BACK TO REAL
                    XP3(IP)=-DIMAG(ZETA1(IP))
                    Y3(IP)=DREAL(ZETA2(IP))
                    YP3(IP)=-DIMAG(ZETA2(IP))
                  ENDDO
C.............................................................
                  QX2=TUNEAPA(X3,XP3,NT(IWIN))  !..WINDOW2
C.............................................................
                  XMIN2=ADVMIN
                  XMAX2=ADVMAX
                  XSIG2=ADVSIG/DSQRT(DFLOAT(NT(IWIN)))
C.............................................................
                  QY2=TUNEAPA(Y3,YP3,NT(IWIN))
C.............................................................
                  YMIN2=ADVMIN
                  YMAX2=ADVMAX
                  YSIG2=ADVSIG/DSQRT(DFLOAT(NT(IWIN)))
C.............................................................
                ELSEIF (ITUNE.EQ.3) THEN
C.............................................................
                  QX2=TUNELASK(X2,XP2,NT(IWIN)) !..WINDOW2
C.............................................................
                  XMIN2=0D0
                  XMAX2=0D0
                  XSIG2=0D0
C.............................................................
                  QY2=TUNELASK(Y2,YP2,NT(IWIN))
C.............................................................
                  YMIN2=0D0
                  YMAX2=0D0
                  YSIG2=0D0
C.............................................................
                ELSEIF (ITUNE.EQ.4) THEN
C.............................................................
                  QX2=TUNENEWT(X2,XP2,NT(IWIN)) !..WINDOW2
C.............................................................
                  XMIN2=0D0
                  XMAX2=0D0
                  XSIG2=0D0
C.............................................................
                  QY2=TUNENEWT(Y2,YP2,NT(IWIN))
C.............................................................
                  YMIN2=0D0
                  YMAX2=0D0
                  YSIG2=0D0
C.............................................................
                ELSEIF (ITUNE.EQ.5) THEN
C.............................................................
                  QX2=TUNEABT2(X2,XP2,NT(IWIN)) !..WINDOW2
C.............................................................
                  XMIN2=0D0
                  XMAX2=0D0
                  XSIG2=0D0
C.............................................................
                  QY2=TUNEABT2(Y2,YP2,NT(IWIN))
C.............................................................
                  YMIN2=0D0
                  YMAX2=0D0
                  YSIG2=0D0
C.............................................................
                ENDIF
C..........................................COMPUTES THE ERRORS
                DQX=QX1-QX2
                DQY=QY1-QY2
C.............................................................
                SUBSX1=3D0*(XSIG1+XSIG2)
                SUBSY1=3D0*(YSIG1+YSIG2)
C.............................................................
                IF (QX1.GE.QX2) THEN
                  SUBSX2=XMIN1+XMAX2
                ELSE
                  SUBSX2=XMAX1+XMIN2
                ENDIF
                IF (QY1.GE.QY2) THEN
                  SUBSY2=YMIN1+YMAX2
                ELSE
                  SUBSY2=YMAX1+YMIN2
                ENDIF
C.............................................................
                DTUNEX=DABS(DQX)
                DTUNEY=DABS(DQY)
C.............................................................
                DQX1(IWIN,IPOINT)=DTUNEX-SUBSX1
                DQY1(IWIN,IPOINT)=DTUNEY-SUBSY1
                DQX2(IWIN,IPOINT)=DTUNEX-SUBSX2
                DQY2(IWIN,IPOINT)=DTUNEY-SUBSY2
C....................................COMPUTES DISTANCE IN TUNE
                DQ(IWIN,IPOINT)=DSQRT(DQX*DQX+DQY*DQY)
C.............................................................
              ELSE
C.............................................................
                DQ(IWIN,IPOINT)=1D0
C.............................................................
                DQX1(IWIN,IPOINT)=0D0
                DQY1(IWIN,IPOINT)=0D0
                DQX2(IWIN,IPOINT)=0D0
                DQY2(IWIN,IPOINT)=0D0
C.............................................................
              ENDIF
C.............................................................
            ENDDO
C.............................................................
          ENDDO
C.............................................................
        ENDDO
C.............................................................
      ENDIF
C.............................................................
      DO IP=1,IPOINT  !..WRITES OUT THE I.C. AND THE STABILITY
        WRITE(25,40) ALPHA(IP),RHO(IP),LOSS(IP)
      ENDDO
C.............................................................
      IF (ITUNE.NE.0) THEN
        DO IP=1,IPOINT
          WRITE(26,50) ALPHA(IP),RHO(IP),(DQ(IWIN,IP),IWIN=1,4)
          WRITE(27,50) ALPHA(IP),RHO(IP),(DQX1(IWIN,IP),IWIN=1,4)
          WRITE(28,50) ALPHA(IP),RHO(IP),(DQY1(IWIN,IP),IWIN=1,4)
C.............................................................
          WRITE(29,50) ALPHA(IP),RHO(IP),(DQX2(IWIN,IP),IWIN=1,4)
          WRITE(30,50) ALPHA(IP),RHO(IP),(DQY2(IWIN,IP),IWIN=1,4)
        ENDDO
      ENDIF
C.............................................................
 40   FORMAT(1X,2F8.4,I12)
 50   FORMAT(1X,2F8.4,4F13.9)
C.............................................................
      RETURN
      END
CDECK  ID>, SMAPDYN5.
C============================================================
C COMPUTES THE DYNAMIC APERTURE AS THE CONNECTED VOLUME IN
C PHASE SPACE OF INITIAL CONDITIONS WHICH ARE STABLE UNDER NTURN
C ITERATIONS. THE INTEGRATION IS CARRIED OUT OVER THE FOUR VARIABLES
C OF PHASE SPACE. POLAR COORDINATES ARE USED:
C      X=R COS ALPHA COS THETA1
C      PX=R COS ALPHA SIN THETA1
C      Y=R SIN ALPHA COS THETA2
C      PY=R SIN ALPHA SIN THETA2
C
C IMAP IS A PARAMETER SPECIFING THE MAP TYPE:
C IMAP = 1 ===> HENON MAP
C IMAP = 2 ===> ACCELERATOR MODEL
C NANG1 IS THE NUMBER OF ANGLES CONSIDERED FOR THETA1 AND THETA2
C NANG2 IS THE NUMBER OF ANGLES CONSIDERED FOR ALPHA
C RMAX IS THE MAXIMUM DISTANCE TO THE ORIGIN CONSIDERED
C RSTEP IS THE STEP USED IN RHO
C IN OUTPUT, ONE HAS
C    AMPL: THE RADIUS OF THE HYPERSPHERE WHICH HAS THE SAME VOLUME
C          OF THE STABLE VOLUME IN PHASE SPACE
C  AMPLMIN: MINUMUM DISTANCE OF THE STABILITY BORDER TO THE ORIGIN
C
C     AUTHOR : E. TODESCO - BOLOGNA UNIVERSITY
C
      SUBROUTINE MAPDYN5(IMAP,NTURN,RMAX,RSTEP,NANG1,NANG2,AMPL,
     .                   AMPLMIN,AMPLMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      PARAMETER (NMAX=50000,NPOINT=50000)
      DIMENSION X(4)
C.............................................................
      EXTERNAL NLFUN
C.............................................................
      IF (IMAP.NE.1.AND.IMAP.NE.2) THEN
        WRITE(6,*) '***ERROR(MAPDYN5): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (NANG1.LT.1.OR.NANG2.LT.1) THEN
        WRITE(6,*) '***ERROR(MAPDYN5): THE NUMBER OF ANGLES IS <1 '
        STOP
      ENDIF
C.............................................................
      IF (RMAX.LT.0) THEN
        WRITE(6,*) '***ERROR(MAPDYN5): THE MAXIMUM RHO IS NEGATIVE'
        STOP
      ENDIF
C.............................................................
      IF (RSTEP.LT.0) THEN
        WRITE(6,*) '***ERROR(MAPDYN5): THE STEP IN RHO IS NEGATIVE'
        STOP
      ENDIF
C.............................................................
      IF (NTURN.GT.NMAX) THEN
        WRITE(6,*) '***ERROR(MAPDYN5): TOO MANY ITERATIONS'
        WRITE(6,*) '***MESSAGE(MAPDYN5): THE NUMBER OF ITERATIONS',
     .             ' MUST BE < ',NMAX
        STOP
      ENDIF
C.............................................................
      PI=4D0*DATAN(1D0)
      TWOPI=8D0*DATAN(1D0)
      IDIME=4
C.............................................................
      AMPLMIN=1D15
      AMPLMAX=0
      DYNAP=0
C.............................................................
      DO I1=1,NANG1                    ! LOOP IN THETA1
        A1=FLOAT(I1)/NANG1*2*PI
        CA1=COS(A1)
        SA1=SIN(A1)
C.............................................................
        DO I2=1,NANG1                  ! LOOP IN THETA2
          A2=FLOAT(I2)/NANG1*2*PI
          CA2=COS(A2)
          SA2=SIN(A2)
C.............................................................
          DO I3=1,NANG2                ! LOOP IN ALPHA
            A3=FLOAT(I3)/NANG2*PI/2
            CA3=COS(A3)
            SA3=SIN(A3)
C.............................................................
            DO RHO=RSTEP,RMAX,RSTEP    ! LOOP IN RHO
              X(1)=RHO*CA3*CA1
              X(2)=RHO*CA3*SA1
              X(3)=RHO*SA3*CA2
              X(4)=RHO*SA3*SA2
C.............................................................
              IF (IMAP.EQ.1) THEN
                ITURN=MANY_HENON(X,NTURN)
              ELSE
                ITURN=MANYTURN(X,NTURN,IDIME,NLFUN)
              ENDIF
C.............................................................
              IF (ITURN.NE.NTURN) THEN
                RHO0=RHO-RSTEP
                IF(RHO0.LT.AMPLMIN) AMPLMIN=RHO0
                IF(RHO0.GT.AMPLMAX) AMPLMAX=RHO0
                DYNAP=DYNAP+(RHO0**4)*SIN(2*A3)
                GOTO 11
              ENDIF
C.............................................................
            ENDDO
            WRITE(6,*) '***ERROR(MAPDYN5): RMAX STABLE'
            AMPLMIN=-1D0
            AMPLMAX=-1D0
            AMPL=-1D0
            RETURN
 11         CONTINUE
C.............................................................
          ENDDO
C.............................................................
        ENDDO
C.............................................................
      ENDDO
C.............................................................
      DYNAP=DYNAP*PI*PI*PI/4/NANG1/NANG1/NANG2
      AMPL=SQRT(SQRT(2*DYNAP/PI/PI))
C.............................................................
      RETURN
      END
CDECK  ID>, SMAPDYN6.
C============================================================
C COMPUTES THE DYNAMIC APERTURE AS THE CONNECTED VOLUME IN
C PHASE SPACE OF INITIAL CONDITIONS WHICH ARE STABLE UNDER NTURN
C ITERATIONS. THE INTEGRATION IS CARRIED OUT OVER TWO VARIABLES
C OF PHASE SPACE (R,ALPHA). PHASE SPACE COORDINATES:
C      X=R COS ALPHA COS THETA1
C      PX=R COS ALPHA SIN THETA1
C      Y=R SIN ALPHA COS THETA2
C      PY=R SIN ALPHA SIN THETA2
C THE INTEGRATION ON THETA1 AND THETA2 IS SUBSTITUTED
C BY THE INTEGRATION OVER THE ORBIT. 
C IMAP IS A PARAMETER SPECIFING THE MAP TYPE:
C IMAP = 1 ===> HENON MAP
C IMAP = 2 ===> ACCELERATOR MODEL
C NANG IS THE NUMBER OF ANGLES CONSIDERED FOR ALPHA
C RMAX IS THE MAXIMUM DISTANCE TO THE ORIGIN CONSIDERED
C RSTEP IS THE STEP USED IN RHO
C IN OUTPUT, ONE HAS
C    AMPL: THE RADIUS OF THE HYPERSPHERE WHICH HAS THE SAME VOLUME
C          OF THE STABLE VOLUME IN PHASE SPACE
C  AMPLMIN: MINUMUM DISTANCE OF THE STABILITY BORDER TO THE ORIGIN
C  AMPLMAX: MAXIMUM DISTANCE OF THE STABILITY BORDER TO THE ORIGIN
C
C     AUTHOR : E. TODESCO - BOLOGNA UNIVERSITY
C

      SUBROUTINE MAPDYN6(IMAP,NTURN,RMAX,RSTEP,NANG,AMPL,AMPLMIN,
     .                   AMPLMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=8000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      PARAMETER (NMAX=50000,NPOINT=50000)
      DIMENSION X(4),XX(NMAX),PXX(NMAX),YY(NMAX),PYY(NMAX)
C.............................................................
      EXTERNAL NLFUN
C.............................................................
      IF (IMAP.NE.1.AND.IMAP.NE.2) THEN
        WRITE(6,*) '***ERROR(MAPDYN6): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (NANG.LT.1) THEN
        WRITE(6,*) '***ERROR(MAPDYN6): THE NUMBER OF ANGLES IS <1 '
        STOP
      ENDIF
C.............................................................
      IF (RMAX.LT.0) THEN
        WRITE(6,*) '***ERROR(MAPDYN6): THE MAXIMUM RHO IS NEGATIVE'
        STOP
      ENDIF
C.............................................................
      IF (RSTEP.LT.0) THEN
        WRITE(6,*) '***ERROR(MAPDYN6): THE STEP IN RHO IS NEGATIVE'
        STOP
      ENDIF
C.............................................................
      IF (NTURN.GT.NMAX) THEN
        WRITE(6,*) '***ERROR(MAPDYN6): TOO MANY ITERATIONS'
        WRITE(6,*) '***MESSAGE(MAPDYN6): THE NUMBER OF ITERATIONS',
     .             ' MUST BE < ',NMAX
        STOP
      ENDIF
C.............................................................
      PI=4D0*DATAN(1D0)
      TWOPI=2D0*PI
      PIH=.5D0*PI
      IDIME=4
C.............................................................
      ASTEP1=PIH/DFLOAT(NANG+1)
C.............................................................
      AMPLMIN=1D15
      AMPLMAX=0D0
      AMPL=0D0
C.............................................................
      DO I1=1,NANG                 ! LOOP IN ALPHA
        A1=DFLOAT(I1)*ASTEP1
        CA1=DCOS(A1)
        SA1=DSIN(A1)
C.............................................................
        DO RHO=RSTEP,RMAX,RSTEP    ! LOOP TO EVALUATE LAST STABLE ORBIT
          X(1)=RHO*CA1
          X(2)=0D0
          X(3)=RHO*SA1
          X(4)=0D0
C.............................................................
          IF (IMAP.EQ.1) THEN
            ITURN=MANY_HENON(X,NTURN)
          ELSE
            ITURN=MANYTURN(X,NTURN,IDIME,NLFUN)
          ENDIF
C.............................................................
          IF (ITURN.NE.NTURN) THEN
            RHO0=RHO-RSTEP
            GOTO 11
          ENDIF
C.............................................................
        ENDDO
        WRITE(6,*) '***ERROR(MAPDYN6): RMAX STABLE'
        AMPLMIN=-1D0
        AMPLMAX=-1D0
        AMPL=-1D0
        RETURN
 11     CONTINUE
C.............................................................
        X(1)=RHO0*CA1          
        X(2)=0D0
        X(3)=RHO0*SA1
        X(4)=0D0
        XX(1)=X(1)
        PXX(1)=X(2)
        YY(1)=X(3)
        PYY(1)=X(4)
C.............................................................
        DO JJ=2,NTURN         ! STORES THE ITERATES OF LAST STABLE ORBIT
          IF (IMAP.EQ.1) THEN
            ITURN1=MANY_HENON(X,1)
          ELSE
            ITURN1=MANYTURN(X,1,IDIME,NLFUN)
          ENDIF
          XX(JJ)=X(1)
          PXX(JJ)=X(2)
          YY(JJ)=X(3)
          PYY(JJ)=X(4)
          R=DSQRT(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)+X(4)*X(4))
          IF(R.LT.AMPLMIN) AMPLMIN=R
          IF(R.GT.AMPLMAX) AMPLMAX=R
        END DO        
C.............................................................
        R=DINTEG(XX,PXX,YY,PYY,NTURN)  ! COMPUTES AVERAGE RADIUS OF ORBIT
        AMPL=AMPL+R*R*R*R*DSIN(2D0*A1)
      END DO
C.............................................................
      AMPL=AMPL*PI*PI*PI/4/NANG
      AMPL=DSQRT(DSQRT(2*AMPL/PI/PI))
C.............................................................
      RETURN
      END
CDECK  ID>, SMAPDYN7.
C============================================================
C COMPUTES THE LYAPUNOV EXPONENTS FOR A SET OF INITIAL 
C CONDITIONS. 
C IMAP IS A PARAMETER SPECIFING THE MAP TYPE:
C IMAP = 1 ===> HENON MAP
C IMAP = 2 ===> ACCELERATOR MODEL
C NTURN IS THE TURN NUMBER USED TO COMPUTE THE STABILITY.
C IF NTURN IS 0 THEN THE LYAPUNOV ARE COMPUTED ON THE FOUR 
C TIME WINDOWS DEFINED IN THE ARRAY NT(4).
C IDIME IS THE DIMENSION OF THE PHASE SPACE.
C (RMIN,RMAX) REPRESENTS THE INTERVAL FOR RADIUS VALUES.
C NR INITIAL CONDITIONS ARE CONSIDERED.
C (ANG1,ANG2) REPRESENTS THE INTERVAL FOR THE ANGLES IN UNIT
C OF TWOPI. NANG INITIAL DIRECTIONS ARE CONSIDERED.
C DELTA IS THE INITIAL DISTANCE BETWEEN NEARBY PARTICLES.
C IT PRODUCES SOME OUTPUT FILES:
C UNIT 25 ===> INITIAL COORDINATES AND STABILITY
C UNIT 26 ===> INITIAL COORDINATES AND LYAPUNOV
C              (SEE LONG WRITE UP)
C
C AUTHOR : M. GIOVANNOZZI - CERN 
C

      SUBROUTINE MAPDYN7(IMAP,NTURN,IDIME,RMIN,RMAX,NR,ANG1,ANG2,NANG,
     .                   DELTA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      PARAMETER (NMAX=100000,NPOINT=50000)
      DIMENSION X(4),ALYAP(4,NPOINT),BLYAP(NMAX)
      DIMENSION RHO(NPOINT),ALPHA(NPOINT)
      DIMENSION LOSS(NPOINT),NT(4)
C.............................................................
      IF (IMAP.NE.1.AND.IMAP.NE.2) THEN
        WRITE(6,*) '***ERROR(MAPDYN7): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (NR*NANG.GT.NPOINT) THEN
        WRITE(6,*) '***ERROR(MAPDYN7): TOO MANY INITIAL CONDITIONS'
        STOP
      ENDIF
C.............................................................
      IF (NTURN.LT.0) THEN
        WRITE(6,*) '***ERROR(MAPDYN7): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      IF (IDIME.NE.2.AND.IDIME.NE.4) THEN
        WRITE(6,*) '***ERROR(MAPDYN7): THIRD PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
C.............................................................
      STEPR=(RMAX-RMIN)/DFLOAT(NR)   !..COMPUTES THE STEP IN THE RADIUS S
      IF (STEPR.LT.0D0) THEN
        STEPR=-STEPR
        WRITE(6,*) '***MESSAGE(MAPDYN7): WRONG ORDER IN R INTERVAL'
        WRITE(6,*) '***MESSAGE(MAPDYN7): ORDER INVERTED'
      ENDIF
C.............................................................
      STEPANG=(ANG2-ANG1)/DFLOAT(NANG+1)
      IF (STEPANG.LT.0D0) THEN
        STEPANG=-STEPANG
        WRITE(6,*) '***MESSAGE(MAPDYN7): WRONG ORDER IN ANGLE INTERVAL'
        WRITE(6,*) '***MESSAGE(MAPDYN7): ORDER INVERTED'
      ENDIF
C.............................................................
      NT(1)=64     !..DEFINES TIME WINDOWS
      NT(2)=256
      NT(3)=1024
      NT(4)=4096
C.............................................................
      TWOPI=8D0*DATAN(1D0)
      IPOINT=0
C.............................................................
      DO IROW=1,NANG            !..LOOP ON ANGLE
        ANGNEW=ANG1+STEPANG*DFLOAT(IROW)
C.............................................................
        DO ICOL=1,NR           !..LOOP ON R
          RNEW=RMIN+STEPR*DFLOAT(ICOL)
C.............................................................
          IF(IDIME.EQ.2) THEN  !..PHASE SPACE DIMENSION
C.............................................................
            X(1)=RNEW*COS(ANGNEW*TWOPI)
            X(2)=RNEW*SIN(ANGNEW*TWOPI)
            X(3)=0D0
            X(4)=0D0
C.............................................................
          ELSE
C.............................................................
            X(1)=RNEW*COS(ANGNEW*TWOPI)
            X(2)=0D0
            X(3)=RNEW*SIN(ANGNEW*TWOPI)
            X(4)=0D0
C.............................................................
          ENDIF
C.............................................................
          IOUT=0                     !..COMPUTATION OF LYAPUNOV
          IF (NTURN.EQ.0) THEN
            ITURN=NT(4)
            IWMAX=4
          ELSE
            ITURN=NTURN
            NT(1)=NTURN
            IWMAX=1
          ENDIF
C.............................................................
          ISTEP=1
          CALL MLYAP(IMAP,ITURN,IDIME,X,DELTA,BLYAP,IOUT,ISTEP)
C.............................................................
          IPOINT=IPOINT+1            !..DEFINTION OF OUTPUT
          LOSS(IPOINT)=ITURN
          RHO(IPOINT)=RNEW
          ALPHA(IPOINT)=ANGNEW
C.............................................................
          DO IWIN=1,IWMAX
C.............................................................
            IF (LOSS(IPOINT).GE.NT(IWIN)) THEN
              ALYAP(IWIN,IPOINT)=BLYAP(NT(IWIN))
            ELSE
              ALYAP(IWIN,IPOINT)=1D0
            ENDIF
C.............................................................
          ENDDO
C.............................................................
        ENDDO
C.............................................................
      ENDDO
C.............................................................
      DO IP=1,IPOINT  !..WRITES OUT THE I.C. AND THE STABILITY
        WRITE(25,40) ALPHA(IP),RHO(IP),LOSS(IP)
      ENDDO
C.............................................................
      DO IP=1,IPOINT
        WRITE(26,50) ALPHA(IP),RHO(IP),(ALYAP(IWIN,IP),IWIN=1,IWMAX)
      ENDDO
C.............................................................
 40   FORMAT(1X,2F8.4,I12)
 50   FORMAT(1X,2F8.4,4F13.9)
C.............................................................
      RETURN
      END
CDECK  ID>, SMAPDYN8.
C=============================================================
C COMPUTES POINTS ON RESONANCE INSIDE THE DYNAMIC APERTURE
C (FOLLOWING E. TODESCO).
C IMAP IS A PARAMETER SPECIFING THE MAP TYPE:
C IMAP = 1 ===> HENON MAP
C IMAP = 2 ===> ACCELERATOR MODEL
C NTURN IS THE TURN NUMBER USED TO COMPUTE THE TUNE.
C (XMIN,XMAX) REPRESENTS THE INTERVAL FOR X VALUES.
C NXS INITIAL CONDITIONS ARE CONSIDERED.
C (YMIN,YMAX) REPRESENTS THE INTERVAL FOR Y VALUES.
C NYS INITIAL DIRECTIONS ARE CONSIDERED.
C NRMIN,NRMAX REPRESENT THE MINIMUM AND MAXIMUM ORDER OF THE
C RESONANCES CHECKED.
C EPS REPRESENTS THE DISTANCE FROM THE RESONANCE.
C IT PRODUCES SOME OUTPUT FILES:
C UNIT 25 ===> INITIAL COORDINATES ON RESONANCE AND RESONANCE
C              PARAMETERS
C
C AUTHOR : M. GIOVANNOZZI - CERN
C

      SUBROUTINE MAPDYN8(IMAP,NTURN,XMIN,XMAX,NXS,YMIN,YMAX,NYS,
     .                   NRMIN,NRMAX,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NELE=80000,NLUMP=NELE/10)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      PARAMETER (NMAX=10000)
      DIMENSION X(4),X1(NMAX),XP1(NMAX),Y1(NMAX),YP1(NMAX)
      REAL TMP
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
C.............................................................
      EXTERNAL NLFUN
C.............................................................
      IF (IMAP.NE.1.AND.IMAP.NE.2) THEN
        WRITE(6,*) '***ERROR(MAPDYN8): FIRST PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (NTURN.LE.0) THEN
        WRITE(6,*) '***ERROR(MAPDYN8): SECOND PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (XMIN.LT.0D0.OR.XMAX.LT.0D0) THEN
        WRITE(6,*) '***ERROR(MAPDYN8): THIRD OR FOURTH PARAMETER OUT ',
     .             'OF BOUNDS'
        STOP
      ENDIF
      IF (NXS.LT.0) THEN
        WRITE(6,*) '***ERROR(MAPDYN8): FIFTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (YMIN.LT.0D0.OR.YMAX.LT.0D0) THEN
        WRITE(6,*) '***ERROR(MAPDYN8): SIXTH OR SEVENTH PARAMETER OUT ',
     .             'OF BOUNDS'
        STOP
      ENDIF
      IF (NYS.LT.0) THEN
        WRITE(6,*) '***ERROR(MAPDYN8): EIGHTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (NRMIN.LT.0) THEN
        WRITE(6,*) '***ERROR(MAPDYN8): NINETH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (NRMAX.LT.0) THEN
        WRITE(6,*) '***ERROR(MAPDYN8): TENTH PARAMETER OUT OF BOUNDS'
        STOP
      ENDIF
      IF (NRMIN.GT.NRMAX) THEN
        WRITE(6,*) '***ERROR(MAPDN8): NINETH OR TENTH PARAMETER OUT ',
     .             'OF BOUNDS'
        STOP
      ENDIF
      IF (EPS.LT.0D0) THEN
        WRITE(6,*) '***ERROR(MAPDYN8): ELEVENTH PARAMETER OUT OF ',
     .             'BOUNDS'
        STOP
      ENDIF
C.............................................................
      STEPX=(XMAX-XMIN)/DFLOAT(NXS)  !..COMPUTES THE STEP IN X
      IF (STEPX.LT.0D0) THEN
        STEPX=-STEPX
        WRITE(6,*) '***MESSAGE(MAPDYN8): WRONG ORDER IN X INTERVAL'
        WRITE(6,*) '***MESSAGE(MAPDYN8): ORDER INVERTED'
      ENDIF
C.............................................................
      STEPY=(YMAX-YMIN)/DFLOAT(NYS)
      IF (STEPY.LT.0D0) THEN
        STEPY=-STEPY
        WRITE(6,*) '***MESSAGE(MAPDYN8): WRONG ORDER IN Y INTERVAL'
        WRITE(6,*) '***MESSAGE(MAPDYN8): ORDER INVERTED'
      ENDIF
C.............................................................
      IDIME=4
C.............................................................
      DO IROW=1,NYS           !..LOOP ON Y
        YNEW=YMIN+STEPY*DFLOAT(IROW)
C.............................................................
        DO ICOL=1,NXS         !..LOOP ON X
          XNEW=XMIN+STEPX*DFLOAT(ICOL)
C.............................................................
          X(1)=XNEW
          X(2)=0D0
          X(3)=YNEW
          X(4)=0D0
C.............................................................
          X0=X(1)
          Y0=X(3)
C.............................................................
          DO JITER=1,NTURN
C.............................................................
            IF (IMAP.EQ.1) THEN
              ITURN=MANY_HENON(X,1)
            ELSE
              ITURN=MANYTURN(X,1,IDIME,NLFUN)
            ENDIF
C.............................................................
            IF (ITURN.LT.0) THEN
              GOTO 10
            ENDIF
C.............................................................
            X1(JITER)=X(1)
            XP1(JITER)=X(2)
            Y1(JITER)=X(3)
            YP1(JITER)=X(4)
          ENDDO
C.............................................................
          QX=TUNEABT2(X1,XP1,NTURN) !..TUNE COMPUTATION
          QY=TUNEABT2(Y1,YP1,NTURN)
C.............................................................
          DO IRES=NRMIN,NRMAX
            DO NX=-IRES,IRES-1
              NY=IRES-ABS(NX)
              TMP=NX*QX+NY*QY
              IF (TMP.GT.0D0) ITUNE=INT(TMP+0.5)
              IF (TMP.LT.0D0) ITUNE=INT(TMP-0.5)
              FRAC=ABS(NX*QX+NY*QY-ITUNE)
              IF (FRAC.LT.EPS) THEN
                WRITE(25,40) X0,Y0,NX,NY
                GOTO 10
              ENDIF
            ENDDO
          ENDDO
 10       CONTINUE
        ENDDO
      ENDDO
C.............................................................
 40   FORMAT(1X,2F8.4,2I12)
C.............................................................
      RETURN
      END
CDECK  ID>, IPRIM.
C===============================================================================
C  RETURNS 1 IF THREE INTEGER NUMBERS HAVE A COMMON DIVISOR. IT IS USED TO 
C  AVOID TO COMPUTE UNUSEFUL RESONANCES 
C
C  AUTHOR: E. TODESCO - INFN
C
C
      INTEGER FUNCTION IPRIM(I1,I2,I3)

      IPRIM=0
      DO I=2,ABS(I1)+ABS(I2)
        IS1=I1/I
        IS2=I2/I
        IS3=I3/I
        IF(IS1*I.EQ.I1.AND.IS2*I.EQ.I2.AND.IS3*I.EQ.I3) IPRIM=1
      END DO

      END
CDECK  ID>, COMP_REALRES.
C===============================================================================
C  COMPUTES WHICH RESONANCES FALL INSIDE THE FOOTPRINT EVALUATED AT AMPLITUDE
C  AMPLMAX AND STORES THEM IN IRERES
C
C  AUTHOR: E. TODESCO - INFN
C
C
      SUBROUTINE COMP_REALRES(A,IRERES,MAXORD)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
C...............................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IR(NORMAX*NORMAX1,5),IRERES(NORMAX*NORMAX1,5)
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
C...............................................................................
      IORD=NMAXH
      NRES=0
C.............................................................
      PI2=8D0*DATAN(1D0)
      OM1=DACOS(RMAP(5,1))/PI2
      OM2=DACOS(RMAP(3,3))/PI2
      DO NN=5,MAXORD+1
        DO N1=NN,0,-1
          NRES=NRES+1
          N2=NN-N1
          SCRA=N1*OM1+N2*OM2          
          IF(SCRA.GT.0) L=INT(SCRA+0.5)
          IF(SCRA.EQ.0) L=0
          IF(SCRA.LT.0) L=INT(SCRA-0.5)
          IR(NRES,1)=N1
          IR(NRES,2)=N2
          IR(NRES,3)=L
          IF(N1*OM1+N2*OM2-L.GT.0) IR(NRES,4)=1
          IF(N1*OM1+N2*OM2-L.LT.0) IR(NRES,4)=-1
          IR(NRES,5)=0
        END DO
        DO N1=NN-1,1,-1
          NRES=NRES+1
          N2=N1-NN
          SCRA=N1*OM1+N2*OM2          
          IF(SCRA.GT.0) L=INT(SCRA+0.5)
          IF(SCRA.EQ.0) L=0
          IF(SCRA.LT.0) L=INT(SCRA-0.5)
          IR(NRES,1)=N1
          IR(NRES,2)=N2
          IR(NRES,3)=L
          IF(N1*OM1+N2*OM2-L.GT.0) IR(NRES,4)=1
          IF(N1*OM1+N2*OM2-L.LT.0) IR(NRES,4)=-1
          IR(NRES,5)=0
        END DO
      END DO
C...............................................................................
      NPUNT=10
      DO J1=0,NPUNT
        RH1=A*J1/DFLOAT(NPUNT)
        IF(RH1.EQ.0) RH1=1.D-10
        DO J2=0,NPUNT-J1
          RH2=A*J2/DFLOAT(NPUNT)
          IF(RH2.EQ.0) RH2=1.D-10
          T1=TUNE1(RH1,RH2,IORD)
          T2=TUNE2(RH1,RH2,IORD)
          DO JJ=1,NRES
            IF(IR(JJ,5).EQ.0.AND.
     .      IPRIM(IR(JJ,1),IR(JJ,2),IR(JJ,3)).NE.1.AND.
     .      (T1*IR(JJ,1)+T2*IR(JJ,2)-IR(JJ,3))*IR(JJ,4).LT.0) IR(JJ,5)=1
          END DO
        END DO
      END DO
C...............................................................................
      KRES=0
      DO JJ=1,NRES
        IF(IR(JJ,5).EQ.1) THEN
          KRES=KRES+1 
          IRERES(KRES,1)=IR(JJ,1)
          IRERES(KRES,2)=IR(JJ,2)
          IRERES(KRES,3)=IR(JJ,3)
        END IF
      END DO
      IRERES(1,5)=KRES
C...............................................................................
      END
CDECK  ID>, W_FLINES.
C===============================================================================
C  WRITES THE POSITION OF THE FIXED LINES 
C  IN INPUT Z1 AND Z2 IS THE POSITION OF THE FIXED LINES IN COMPLEX C.S. 
C  COORDINATES. I IS -1 IN THE HYPERBOLIC CASE AND +1 IN THE ELLIPTIC CASE.
C  STAB IS THE EIGENVALUE
C  IN UNIT 32 ONE HAS THE OUTPUT IN REAL C.S. COORDINATES 
C  CHANGE THE COMMENTS IF YOU WANT THE OUTPUT IN PHYSICAL COORDINATES 
C
C  AUTHOR: E. TODESCO - INFN
C
C
      SUBROUTINE W_FLINES(I,Z1,Z2,STAB)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 TYP
      COMPLEX*16 Z1,Z2
      COMMON/INVARIANT/R2
      COMMON/CSTRANS/CSM(4,4),CSMIN(4,4)
      IF(I.EQ.1) TYP='E'
      IF(I.EQ.-1) TYP='H'
C------------------output in c-s coordinates------------------------------------
      WRITE(32,990) TYP,R2,DREAL(Z1),-DIMAG(Z1),DREAL(Z2),
     .                          -DIMAG(Z2),SQRT(ABS(STAB))
C------------------output in physical coordinates-------------------------------
C      WRITE(42,990) TYP,R2,DREAL(Z1)*CSMIN(1,1),
C     .              DREAL(Z1)*CSMIN(2,1)-DIMAG(Z1)*CSMIN(2,2),
C     .              DREAL(Z2)*CSMIN(3,3),
C     .              DREAL(Z2)*CSMIN(4,3)-DIMAG(Z2)*CSMIN(4,4),
C     .              SQRT(ABS(STAB))

 990  FORMAT(1X,A1,E13.5,1X,4E13.5,1X,E13.5)
      END
CDECK  ID>, ANAL_FOOT.
C===============================================================================
C  COMPUTES THE ANALYTICAL FOOTPRINT USING NONRESONANT NORMAL FORMS
C  A IS THE MAXIMUM SUM OF THE NONLINEAR EMITTANCES
C  NPUNT INITIAL CONDITIONS ARE TAKES IN RHO1 AND IN RHO2, WITH RHO1+RHO2<A
C  THE CORRESPONDING NONLINEAR TUNES ARE WRITTEN ON FILE IUNIT
C
C  AUTHOR: E. TODESCO - INFN
C
C
      SUBROUTINE ANAL_FOOT(A,NPUNT,IUNIT) 

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      T1MAX=-1
      T1MIN=10
      T2MAX=-1
      T2MIN=10
      IORD=NMAXH
C...............................................................................
      DO J1=0,NPUNT
        RH1=A*J1/DFLOAT(NPUNT)
        IF(RH1.EQ.0) RH1=1.D-10
        DO J2=0,NPUNT-J1
          RH2=A*J2/DFLOAT(NPUNT)
          IF(RH2.EQ.0) RH2=1.D-10
          T1=TUNE1(RH1,RH2,IORD)
          T2=TUNE2(RH1,RH2,IORD)
          WRITE(IUNIT,'(1X,2F11.5)') T1,T2
          IF(T1.GT.T1MAX) T1MAX=T1
          IF(T1.LT.T1MIN) T1MIN=T1
          IF(T2.GT.T2MAX) T2MAX=T2
          IF(T2.LT.T2MIN) T2MIN=T2
        END DO
      END DO
      WRITE(IUNIT,*) T1MIN,T1MAX
      WRITE(IUNIT,*) T2MIN,T2MAX

      END
CDECK  ID>, HC0.
C===============================================================================
C  COMPUTES THE VALUE OF THE SINGLE RESONANT INTERPOLATING HAMILTONIAN 
C  AT THE NONLINEAR INVARIANTS RH10,RH20 NEGLECTING THE RESONANT PART
C  THIS IS NEVER USED BUT INSERTED FOR COMPLETENESS
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION HC0(RH10,RH20)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 RH10,RH20
C...............................................................................
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      HC0=0.D0
      DO N=1,NMAXH/2
        DO K=0,N
          HC0=HC0+A1(K,N-K,0)*DCOS(ALPHA1(K,N-K,0))*
     .           RH10**K*RH20**(N-K)
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, DHC0D10.
C===============================================================================
C  COMPUTES THE DERIVATIVE OF THE SINGLE RESONANT INTERPOLATING HAMILTONIAN 
C  WITH RESPECT TO THE FIRST INVARIANT AT THE NONLINEAR INVARIANTS RH10,RH20 
C  NEGLECTING THE RESONANT PART.
C  THIS IS USED FOR EVALUATING THE POSITION OF THE FIXED LINE.
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION DHC0D10(RH10,RH20)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 RH10,RH20
C...............................................................................
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      DHC0D10=0.D0
      DO N=1,NMAXH/2
        DO K=0,N
          DHC0D10=DHC0D10+A1(K,N-K,0)*DCOS(ALPHA1(K,N-K,0))*
     .           K*RH10**(K-1)*RH20**(N-K)
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, DHC0D01.
C===============================================================================
C  COMPUTES THE DERIVATIVE OF THE SINGLE RESONANT INTERPOLATING HAMILTONIAN 
C  WITH RESPECT TO THE SECOND INVARIANT AT THE NONLINEAR INVARIANTS RH10,RH20 
C  NEGLECTING THE RESONANT PART.
C  THIS IS USED FOR EVALUATING THE POSITION OF THE FIXED LINE.
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION DHC0D01(RH10,RH20)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 RH10,RH20
C...............................................................................
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      DHC0D01=0.D0
      DO N=1,NMAXH/2
        DO K=0,N
          DHC0D01=DHC0D01+A1(K,N-K,0)*DCOS(ALPHA1(K,N-K,0))*
     .           RH10**K*(N-K)*RH20**(N-K-1)
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, DHC0D20.
C===============================================================================
C  COMPUTES THE SECOND DERIVATIVE OF THE SINGLE RESONANT INTERPOLATING 
C  HAMILTONIAN WITH RESPECT TO THE FIRST INVARIANT AT THE NONLINEAR INVARIANTS 
C  RH10,RH20 NEGLECTING THE RESONANT PART.
C  THIS IS USED FOR EVALUATING THE POSITION OF THE FIXED LINE.
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION DHC0D20(RH10,RH20)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 RH10,RH20
C...............................................................................
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      DHC0D20=0.D0
      DO N=1,NMAXH/2
        DO K=0,N
          DHC0D20=DHC0D20+A1(K,N-K,0)*DCOS(ALPHA1(K,N-K,0))*
     .           K*(K-1)*RH10**(K-2)*RH20**(N-K)
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, DHC0D11.
C===============================================================================
C  COMPUTES THE SECOND DERIVATIVE OF THE SINGLE RESONANT INTERPOLATING 
C  HAMILTONIAN WITH RESPECT TO THE FIRST AND TO THE SECOND INVARIANT AT THE 
C  NONLINEAR INVARIANTS RH10,RH20 NEGLECTING THE RESONANT PART.
C  THIS IS USED FOR EVALUATING THE POSITION OF THE FIXED LINE.
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION DHC0D11(RH10,RH20)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 RH10,RH20
C...............................................................................
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      DHC0D11=0.D0
      DO N=1,NMAXH/2
        DO K=0,N
          DHC0D11=DHC0D11+A1(K,N-K,0)*DCOS(ALPHA1(K,N-K,0))*
     .           K*RH10**(K-1)*(N-K)*RH20**(N-K-1)
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, DHC0D02.
C===============================================================================
C  COMPUTES THE SECOND DERIVATIVE OF THE SINGLE RESONANT INTERPOLATING 
C  HAMILTONIAN WITH RESPECT TO THE SECOND INVARIANT AT THE NONLINEAR INVARIANTS 
C  RH10,RH20 NEGLECTING THE RESONANT PART.
C  THIS IS USED FOR EVALUATING THE POSITION OF THE FIXED LINE.
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION DHC0D02(RH10,RH20)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 RH10,RH20
C...............................................................................
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      DHC0D02=0.D0
      DO N=1,NMAXH/2
        DO K=0,N
          DHC0D02=DHC0D02+A1(K,N-K,0)*DCOS(ALPHA1(K,N-K,0))*
     .           RH10**K*(N-K)*(N-K-1)*RH20**(N-K-2)
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, HC1.
C===============================================================================
C  COMPUTES THE SINGLE RESONANT INTERPOLATING HAMILTONIAN 
C  AT PS,RH WITH THE SECOND INVARIANT R2 FIXED,
C  TRUNCATING AT THE FIRST SIGNIFICANT RESONANT ORDER.
C  THIS IS NEVER USED BUT INSERTED FOR COMPLETENESS.
C
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION HC1(PS,RH)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 PS,RH,R2
C...............................................................................
      COMMON/INVARIANT/R2
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      HC1=0
      N=KRES
      DO L=0,N/KSUM
        DO K2=0,(N-L*KSUM)/2
          JJ=(N-L*KSUM-2*K2)/2
          IF(2*JJ.EQ.N-L*KSUM-2*K2) THEN
            K1 = (N-L*KSUM-2*K2)/2
            HC1=HC1+A1(K1,K2,L)*RH**(K1+L*K11/2.D0)*
     .                      (PQ*RH+R2)**(K2+L*IAK12/2.D0)*
     .           DCOS(L*K11*PS+ALPHA1(K1,K2,L))
          END IF
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, DHC1D10.
C===============================================================================
C  COMPUTES THE FIRST DERIVATIVE OF THE SINGLE RESONANT INTERPOLATING 
C  HAMILTONIAN WITH RESPECT TO THE ANGLE PS
C  AT PS,RH WITH THE SECOND INVARIANT R2 FIXED,
C  TRUNCATING AT THE FIRST SIGNIFICANT RESONANT ORDER. 
C  THIS IS USED FOR THE EVALUATION OF THE FIRST GUESS OF THE ANGULAR SOLUTION 
C  OF THE FIXED LINES.
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION DHC1D10(PS,RH)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
C...............................................................................
      COMMON/INVARIANT/R2
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      DHC1D10=0
      N=KRES
      DO L=1,N/KSUM
        DO K2=0,(N-L*KSUM)/2
          JJ=(N-L*KSUM-2*K2)/2
          IF(2*JJ.EQ.N-L*KSUM-2*K2) THEN
            K1 = (N-L*KSUM-2*K2)/2
            DHC1D10=DHC1D10-A1(K1,K2,L)*RH**(K1+L*K11/2.D0)*
     .                      (PQ*RH+R2)**(K2+L*IAK12/2.D0)*
     .                      L*K11*DSIN(L*K11*PS+ALPHA1(K1,K2,L))
          END IF
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, DHC1D20.
C===============================================================================
C  COMPUTES THE SECOND DERIVATIVE OF THE SINGLE RESONANT INTERPOLATING 
C  HAMILTONIAN WITH RESPECT TO THE ANGLE PS
C  AT PS,RH WITH THE SECOND INVARIANT R2 FIXED,
C  TRUNCATING AT THE FIRST SIGNIFICANT RESONANT ORDER.
C  THIS IS USED FOR THE EVALUATION OF THE FIRST GUESS OF THE ANGULAR SOLUTION 
C  OF THE FIXED LINES.
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION DHC1D20(PS,RH)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 PS,RH,R2
C...............................................................................
      COMMON/INVARIANT/R2
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      DHC1D20=0
      N=KRES
      DO L=1,N/KSUM
        DO K2=0,(N-L*KSUM)/2
          JJ=(N-L*KSUM-2*K2)/2
          IF(2*JJ.EQ.N-L*KSUM-2*K2) THEN
            K1 = (N-L*KSUM-2*K2)/2
            DHC1D20=DHC1D20-A1(K1,K2,L)*RH**(K1+L*K11/2.D0)*
     .                   (PQ*RH+R2)**(K2+L*IAK12/2.D0)*
     .                   L*L*K11*K11*DCOS(L*K11*PS+ALPHA1(K1,K2,L))
          END IF
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, HC2.
C===============================================================================
C  COMPUTES THE SINGLE RESONANT INTERPOLATING 
C  HAMILTONIAN AT PS,RH WITH THE SECOND INVARIANT R2 FIXED,
C  TRUNCATING AT THE MAXIMUM ORDER NMAXH.
C  THIS IS USED FOR THE EVALUATION OF THE ISLAND WIDTH.
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION HC2(PS,RH)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 PS,RH,R2
C...............................................................................
      COMMON/INVARIANT/R2
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      HC2=0
      DO N=1,NMAXH
        DO L=0,N/KSUM
          DO K2=0,(N-L*KSUM)/2
            JJ=(N-L*KSUM-2*K2)/2
            IF(2*JJ.EQ.N-L*KSUM-2*K2) THEN
              K1 = (N-L*KSUM-2*K2)/2
              HC2=HC2+A1(K1,K2,L)*RH**(K1+L*K11/2.D0)*
     .                      (PQ*RH+R2)**(K2+L*IAK12/2.D0)*
     .           DCOS(L*K11*PS+ALPHA1(K1,K2,L))
            END IF
          END DO
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, DHC2D10.
C===============================================================================
C  COMPUTES THE DERIVATIVE OF THE SINGLE RESONANT INTERPOLATING 
C  HAMILTONIAN AT PS,RH WITH RESPECT TO PS 
C  WHERE THE SECOND INVARIANT R2 IS FIXED,
C  TRUNCATING AT THE MAXIMUM ORDER NMAXH.
C  THIS IS USED FOR THE EVALUATION OF THE FIXED LINE POSITION.
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION DHC2D10(PS,RH)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 PS,RH,R2
C...............................................................................
      COMMON/INVARIANT/R2
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      DHC2D10=0
      DO N=1,NMAXH
        DO L=0,N/KSUM
          DO K2=0,(N-L*KSUM)/2
            JJ=(N-L*KSUM-2*K2)/2
            IF(2*JJ.EQ.N-L*KSUM-2*K2) THEN
              K1 = (N-L*KSUM-2*K2)/2
              DHC2D10=DHC2D10-A1(K1,K2,L)*RH**(K1+L*K11/2.D0)*
     .                       (PQ*RH+R2)**(K2+L*IAK12/2.D0)*
     .                 L*K11*DSIN(L*K11*PS+ALPHA1(K1,K2,L))
            END IF
          END DO
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, DHC2D01.
C===============================================================================
C  COMPUTES THE DERIVATIVE OF THE SINGLE RESONANT INTERPOLATING 
C  HAMILTONIAN AT PS,RH WITH RESPECT TO RH
C  WHERE THE SECOND INVARIANT R2 IS FIXED,
C  TRUNCATING AT THE MAXIMUM ORDER NMAXH.
C  THIS IS USED FOR THE EVALUATION OF THE FIXED LINE POSITION.
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION DHC2D01(PS,RH)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 PS,RH,R2
C...............................................................................
      COMMON/INVARIANT/R2
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      DHC2D01=0
      DO N=1,NMAXH
        DO L=0,N/KSUM
          DO K2=0,(N-L*KSUM)/2
            JJ=(N-L*KSUM-2*K2)/2
            IF(2*JJ.EQ.N-L*KSUM-2*K2) THEN
              K1 = (N-L*KSUM-2*K2)/2
              DHC2D01=DHC2D01+A1(K1,K2,L)*RH**(K1+L*K11/2.D0-1)*
     .                      (PQ*RH+R2)**(K2+L*IAK12/2.D0-1)*
     .         ((K1+L*K11/2.D0)*(PQ*RH+R2)+PQ*(K2+L*IAK12/2.D0)*RH)*
     .           DCOS(L*K11*PS+ALPHA1(K1,K2,L))
            END IF
          END DO
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, DHC2D11.
C===============================================================================
C  COMPUTES THE SECOND DERIVATIVE OF THE SINGLE RESONANT INTERPOLATING 
C  HAMILTONIAN AT PS,RH WITH RESPECT TO PS AND RH
C  WHERE THE SECOND INVARIANT R2 IS FIXED,
C  TRUNCATING AT THE MAXIMUM ORDER NMAXH
C  THIS IS USED FOR THE EVALUATION OF THE FIXED LINE POSITION.
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION DHC2D11(PS,RH)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 PS,RH,R2
C...............................................................................
      COMMON/INVARIANT/R2
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      DHC2D11=0
      DO N=1,NMAXH
        DO L=0,N/KSUM
          DO K2=0,(N-L*KSUM)/2
            JJ=(N-L*KSUM-2*K2)/2
            IF(2*JJ.EQ.N-L*KSUM-2*K2) THEN
              K1 = (N-L*KSUM-2*K2)/2
              DHC2D11=DHC2D11-A1(K1,K2,L)*RH**(K1+L*K11/2.D0-1)*
     .                      (PQ*RH+R2)**(K2+L*IAK12/2.D0-1)*
     .         ((K1+L*K11/2.D0)*(PQ*RH+R2)+PQ*(K2+L*IAK12/2.D0)*RH)*
     .           L*K11*DSIN(L*K11*PS+ALPHA1(K1,K2,L))
            END IF
          END DO
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, DHC2D20.
C===============================================================================
C  COMPUTES THE SECOND DERIVATIVE OF THE SINGLE RESONANT INTERPOLATING 
C  HAMILTONIAN AT PS,RH WITH RESPECT TO PS 
C  WHERE THE SECOND INVARIANT R2 IS FIXED,
C  TRUNCATING AT THE MAXIMUM ORDER NMAXH
C  THIS IS USED FOR THE EVALUATION OF THE FIXED LINE POSITION.
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION DHC2D20(PS,RH)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 PS,RH,R2
C...............................................................................
      COMMON/INVARIANT/R2
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      DHC2D20=0
      DO N=1,NMAXH
        DO L=0,N/KSUM
          DO K2=0,(N-L*KSUM)/2
            JJ=(N-L*KSUM-2*K2)/2
            IF(2*JJ.EQ.N-L*KSUM-2*K2) THEN
              K1 = (N-L*KSUM-2*K2)/2
              DHC2D20=DHC2D20-A1(K1,K2,L)*RH**(K1+L*K11/2.D0)*
     .                       (PQ*RH+R2)**(K2+L*IAK12/2.D0)*
     .             L*K11*L*K11*DCOS(L*K11*PS+ALPHA1(K1,K2,L))
            END IF
          END DO
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, DHC2D02.
C===============================================================================
C  COMPUTES THE SECOND DERIVATIVE OF THE SINGLE RESONANT INTERPOLATING 
C  HAMILTONIAN AT PS,RH WITH RESPECT TO RH
C  WHERE THE SECOND INVARIANT R2 IS FIXED,
C  TRUNCATING AT THE MAXIMUM ORDER NMAXH
C  THIS IS USED FOR THE EVALUATION OF THE FIXED LINE POSITION.
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION DHC2D02(PS,RH)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 PS,RH,R2
C...............................................................................
      COMMON/INVARIANT/R2
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      DHC2D02=0
      DO N=1,NMAXH
        DO L=0,N/KSUM
          DO K2=0,(N-L*KSUM)/2
            JJ=(N-L*KSUM-2*K2)/2
            IF(2*JJ.EQ.N-L*KSUM-2*K2) THEN
              K1 = (N-L*KSUM-2*K2)/2
              DHC2D02=DHC2D02+A1(K1,K2,L)*RH**(K1+L*K11/2.D0-2)*
     .                      (PQ*RH+R2)**(K2+L*IAK12/2.D0-2)*
     .         ((K1+L*K11/2.D0)*(K1+L*K11/2.D0-1)*(PQ*RH+R2)**2+
     .           2*PQ*(K1+L*K11/2.D0)*(K2+L*IAK12/2.D0)*RH*(PQ*RH+R2)+
     .           PQ*PQ*(K2+L*IAK12/2.D0)*(K2+L*IAK12/2.D0-1)*RH**2)*
     .           DCOS(L*K11*PS+ALPHA1(K1,K2,L))
            END IF
          END DO
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, DHC3D2.
C===============================================================================
C  COMPUTES THE DERIVATIVE OF THE SINGLE RESONANT INTERPOLATING 
C  HAMILTONIAN AT RH WITH RESPECT TO RH
C  WHERE THE SECOND INVARIANT R2 IS FIXED, THE RESONANT PART IS NEGLECTED, 
C  AND TRUNCATING AT THE MAXIMUM ORDER NMAXH.
C  THIS IS USED FOR EVALUATING THE FIRST GUESS OF THE ISLAND WIDTH
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION DHC3D2(RH)
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 RH,R2
C...............................................................................
      COMMON/INVARIANT/R2
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C...............................................................................
      DHC3D2=0.D0
      DO N=2,NMAXH/2
        DO K=0,N
          DHC3D2=DHC3D2+A1(K,N-K,0)*DCOS(ALPHA1(K,N-K,0))*
     .           (K*(K-1)*RH**(K-2)*(PQ*RH+R2)**(N-K)+
     .           2*PQ*(N-K)*K*RH**(K-1)*(PQ*RH+R2)**(N-K-1)+
     .           PQ*PQ*(N-K)*(N-K-1)*RH**K*(PQ*RH+R2)**(N-K-2))
        END DO
      END DO
C...............................................................................
      END
CDECK  ID>, SWITCH.
C===============================================================================
C  SWITCHES THE X,PX WITH THE Y,PY PLANE AND CONSEQUENTLY THE RESONANT BASIS
C  K11 K12 FOR A SINGLE RESONANCE.
C  THIS IS USED FOR ANALYSING A SINGLE UNCOUPLED RESONANCE IN THE Y,PY PLANE
C  USING THE SAME CODE BUiLT FOR ANALYSING SINGLE UNCOUPLED RESONANCES IN THE 
C  X,PX PLANE.
C  MOREOVER, IT IS USED TO DEAL WITH RESONANCE [1,2] USING THE CODE OF RESONANCE
C  [2,1] 
C
C  AUTHOR: E. TODESCO - INFN
C

      SUBROUTINE SWITCH
C...............................................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...............................................................................
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
C...............................................................................
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C......................SWITCHES THE HAMILTONIAN COEFFICIENTS....................
      DO N=2,NMAXH
        DO L=0,N/(KSUM)
          DO K1=0,(N-L*KSUM)/2/2
            K2=(N-L*KSUM-2*K1)/2
            IF(K2*2.EQ.N-L*KSUM-2*K1) THEN
              SCR=A1(K1,K2,L)
              A1(K1,K2,L)=A1(K2,K1,L)
              A1(K2,K1,L)=SCR
              SCR=ALPHA1(K1,K2,L)
              ALPHA1(K1,K2,L)=ALPHA1(K2,K1,L)
              ALPHA1(K2,K1,L)=SCR
            END IF
          END DO
        END DO
      END DO
C......................SWITCHES THE RESONANT BASIS..............................
      ISCR=K11
      K11=K12
      K12=ISCR
      IF(K11.LT.0) THEN
        K11=ABS(K11)
        K12=-ABS(K12)
      END IF
      KSUM=K11+ABS(K12)
      IAK12=ABS(K12)
C...............................................................................
      END 
CDECK  ID>, FIXLIN.
C===============================================================================
C  MAKES A COMPLETE ANALYSIS OF A SINGLE RESONANCE, STARTING FROM THE   
C  COEFFICIENTS OF THE INTERPOLATING HAMILTONIAN IN ACTION ANGLE VARIABLE STORED
C  IN A1, ALPHA1 (COMMON BLOCK /HAMILT/).
C
C  IFLAG=1 MAKES THE PHASE SPACE ANALYSIS, I.E. THE EVALUATION OF THE POSITION
C          OF THE RESONANCES AND OF THEIR WIDTHS
C  IFLAG=0 ONLY EVALUATES THE NORMAL FORM COEFFICIENTS AND THE RELATIVE QUALITY 
C          FACTORS
C          
C  HYPVOL: HYPERVOLUME OF THE RESONANCE - IT GIVES 0.999 IF THE RESONANCE IS
C          EXCITED BUT THE COMPUTATION IS NOT POSSIBLE
C
C  OUTPUT FILES:
C  UNIT 31: HAMILTONIAN COEFFICIENTS, RESONANCE EXCITATION AND THE QUALITY 
C           FACTOR Q1 (MAP NORM) Q2 (TUNESHIFT NORM) AND Q3 (NORM OF THE 
C           RESONANT PART OF THE HAMILTONIAN)
C  UNIT 32: AMPLITUDES AND THE EIGENVALUES OF THE FIXED LINES (IFLAG=1)
C  UNIT 33: POSITION AND WIDTHS OF THE ISLANDS IN THE AMPLITUDE PLANE (IFLAG=1)
C
C  UNIT 38: THE DATA RELATIVE TO THE FIXED LINES THAT ARE NOT FOUND (I.E.,
C           VIRTUAL OR PATHOLOGICAL) ARE PRINTED OUT (IFLAG=1)
C  UNIT 39: ALL THE ERROR MESSAGES (I.E. NONCONVERGENT NEWTON, ETC) ARE 
C           PRINTED OUT (IFLAG=1)
C
C  AUTHOR: E. TODESCO - INFN
C

      SUBROUTINE FIXLIN(AMPLMAX,IFLAG,HYPVOL)
 
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMPLEX*16 Z1,Z2,ZETA1,ZETA2,AUTORES
      REAL*8 a0(0:NORMAX1,0:NORMAX1),a1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 a2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)

      REAL*8 RH10,DRH10,RH1P,RH1M,DRH1P,DRH1M
      REAL*8 PS1P,PS1M,DPS1P,DPS1M,PSEL,PSHY,PS,DPS,Q2i,Q3
      REAL*8 STABP,STABM,SEPEN,DETP,DETM,KAPPA,c
      REAL*8 RHIN,RHOUT,DRHIN,DRHOUT,WIDTH,DEL,RHINST,RHOUTST
      REAL*8 AMPLMAX,EPSNEWT,HYPVOL,STEPT,OVRFLW
      REAL*8 R2,R2OLD
      REAL*8 HC2,DHC2D01,DHC2D10,DHC2D02,DHC2D11,DHC2D20,DHC3D2
      REAL*8 HC1,DHC1D10,DHC1D20,SIGNOLD,SIGN
      REAL*8 PI

      REAL*8 RH20,DRH20,F1,F2
      REAL*8 P2Q2,DET,H20,H11,H02,H10,H01,T
      REAL*8 DHC0D01,DHC0D10,DHC0D20,DHC0D11,DHC0D02
      INTEGER FLAGSWITCH
      CHARACTER*1 RESTYP

      COMMON/INVARIANT/R2
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
      COMMON/RISONANZA/IRES1(2),IRES2(2),AUTORES(4),IRESON,ICASE
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
      COMMON/CSTRANS/CSM(4,4),CSMIN(4,4)
C-------------------------------------------------------------------------------
C                Resonance check
C-------------------------------------------------------------------------------
      IF(ICASE.NE.1) THEN
        WRITE(*,*) 
     .  '***ERROR(FIXED_LINES_4D): YOU DO NOT HAVE A SINGLE RESONANCE'
        STOP
      END IF
C-------------------------------------------------------------------------------
C                Defult settings
C-------------------------------------------------------------------------------
      MAXNEWT=30             !.......MAXIMUM ITERATIONS OF NEWTON
      EPSNEWT=1.D-10         !.......PRECISION REQUIRED BY NEWTON
      OVRFLW=1.D10           
      PI=4.D0*DATAN(1.D0)
      NTH=100                !.......STEPS IN THE SCANNING OF THE ANGLES
      NSTEPT=500             !.......STEPS IN THE SCANNING OF THE RESONANT LINE
      RESTYP='V'
C-------------------------------------------------------------------------------
C                Selection of switch
C-------------------------------------------------------------------------------
      IF(K11.NE.0) FLAGSWITCH=0
      IF(K11.EQ.0) FLAGSWITCH=1
C-------------------------------------------------------------------------------
c                analysis of resonance excitation:
C                KRES is the lowest order of a resonant cotribution
C                MULRES is the multiplicity of the resonance
C-------------------------------------------------------------------------------
      DO N=2,NMAXH
        DO L=0,N/(KSUM)
          DO K1=0,(N-L*KSUM)/2
            ISCR=(N-L*KSUM-2*K1)/2
            IF(ISCR*2.EQ.N-L*KSUM-2*K1) THEN
              K2=ISCR
              IF(L.GE.1.AND.A1(K1,K2,L).NE.0) THEN
                KRES=N
                GOTO 11
              END IF
            END IF
          END DO
        END DO
      END DO
      KRES=NMAXH+1
 11   CONTINUE

      N=KRES
      DO L=1,N/(KSUM)
        SUM=0
        DO K1=0,(N-L*KSUM)/2
          K2=(N-L*KSUM-2*K1)/2
          IF(K2*2.EQ.N-L*KSUM-2*K1) SUM=SUM+A1(K1,K2,L)
        END DO
        IF(SUM.NE.0) THEN
          MULRES=L
          GOTO 12
        END IF
      END DO
 12   CONTINUE
      IF(1.EQ.1) THEN
C-------------------------------------------------------------------------------
C                Defult settings
C-------------------------------------------------------------------------------
        IF(FLAGSWITCH.EQ.1) CALL SWITCH
        PQ=DFLOAT(K12)/DFLOAT(K11)
        STEPT=.4/NSTEPT
        IF(K12.EQ.0) A1(0,1,0)=0          !...SINGLE UNC. RESON. IMPLIES EPS2=0
        H10=DCOS(ALPHA1(1,0,0))*A1(1,0,0)
        H01=DCOS(ALPHA1(0,1,0))*A1(0,1,0)
        H20=DCOS(ALPHA1(2,0,0))*A1(2,0,0)
        H11=DCOS(ALPHA1(1,1,0))*A1(1,1,0)
        H02=DCOS(ALPHA1(0,2,0))*A1(0,2,0)
        P2Q2=SQRT(1.D0/(K11**2+K12**2))
        R2=0
        HYPVOL=0
        FLAGINIZ=0
        IF(KRES.EQ.3) THEN
C-------------------------------------------------------------------------------
C                Resonances excited at order three
C-------------------------------------------------------------------------------
          IF(K11.EQ.3.AND.K12.EQ.0) THEN
C-------------------------------------------------------------------------------
C                Resonance [3,0] or [0,3]
C-------------------------------------------------------------------------------
            DO R2=EPSNEWT,AMPLMAX+2*EPSNEWT,AMPLMAX/NSTEPT
C-------------------------------------------------------------------------------
C                First guess for the hyperbolic fixed point 
C-------------------------------------------------------------------------------
              IF(FLAGINIZ.EQ.0) THEN 
                RH10=4*(H10+R2*H11)**2/(9*A1(0,0,1)**2)
                IF(H10+R2*H11.GT.0) THEN
                  PS=(PI-ALPHA1(0,0,1))/3              
                ELSE IF(H10+R2*H11.LT.0) THEN
                  PS=-ALPHA1(0,0,1)/3              
                END IF
              END IF
C-------------------------------------------------------------------------------
C                Computation of the hyperbolic fixed point (arbitrary order) 
C-------------------------------------------------------------------------------
              DO JNEWT=1,MAXNEWT
                IF(RH10.GT.OVRFLW.OR.RH10.LT.0) THEN
                  IF(IFLAG.EQ.1) WRITE(34,992) 
     .                   R2,'Rh10 out of bounds     ',RH10,R2
                  FLAGINIZ=0
                  GOTO 199
                END IF
                DET=DHC2D20(PS,RH10)*DHC2D02(PS,RH10)-
     .              DHC2D11(PS,RH10)**2
                DRH10=(-DHC2D20(PS,RH10)*DHC2D01(PS,RH10)+
     .              DHC2D11(PS,RH10)*DHC2D10(PS,RH10))/DET
                DPS=(DHC2D11(PS,RH10)*DHC2D01(PS,RH10)-
     .              DHC2D02(PS,RH10)*DHC2D10(PS,RH10))/DET
                IF(ABS(DRH10).LT.EPSNEWT.AND.ABS(DPS).LT.EPSNEWT
     .              ) GOTO 200
                RH10=RH10+DRH10
                PS=PS+DPS
              END DO
              IF(IFLAG.EQ.1) WRITE(34,997) R2, 
     .            'Divergent Newton for hyperbolic fixed lines       '
              FLAGINIZ=0
 200          CONTINUE
C-------------------------------------------------------------------------------
C                Eigenvalues (arbitrary order)
C-------------------------------------------------------------------------------
              STAB=DHC2D02(PS,RH10)*DHC2D20(PS,RH10)-
     .            DHC2D11(PS,RH10)**2
              IF(STAB.GT.0) THEN
                IF(IFLAG.EQ.1) WRITE(34,992) 
     .                   R2,'Wrong stability        ',STAB
                FLAGINIZ=0
                GOTO 199
              END IF
              FLAGINIZ=1
C-------------------------------------------------------------------------------
C                Evaluation of the whole separatrix
C-------------------------------------------------------------------------------
              RHMIN=RH10
              ESEP=HC2(PS,RH10)
              DO PSI=PS+PI/100,PS+2*PI/3.d0,PI/50
                RHO=RH10*(DCOS(PI/3-ABS(PSI-PS-PI/3)))**2
                DO JNEWT=1,MAXNEWT
                  IF(RHO.GT.OVRFLW.OR.RHO.LE.0) THEN
                    IF(IFLAG.EQ.1) WRITE(34,992) 
     .                R2,'Separatrix out of bound',RHO,R2
                    GOTO 199
                  END IF
                  DRHO=-(HC2(PSI,RHO)-ESEP)/DHC2D01(PSI,RHO)
                  IF(ABS(DRHO).LT.EPSNEWT) GOTO 122
                  RHO=RHO+DRHO
                END DO
                IF(IFLAG.EQ.1) WRITE(34,997) R2, 
     .            'Divergent Newton for separatrix                   '
                GOTO 199
 122            CONTINUE
                IF(RHO.LT.RHMIN) RHMIN=RHO
 123            CONTINUE
              END DO
              IF(IFLAG.EQ.1) THEN
C-------------------------------------------------------------------------------
C                Output on for033.dat of resonant patterns
C-------------------------------------------------------------------------------
                IF(FLAGSWITCH.EQ.0) THEN 
                  WRITE(33,991) 0.,0.,RHMIN,R2,RH10,R2,K11,K12
                ELSE IF(FLAGSWITCH.EQ.1) THEN 
                  WRITE(33,991) 0.,0.,R2,RHMIN,R2,RH10,K12,K11
                END IF
C-------------------------------------------------------------------------------
C                Fixed points data output
C-------------------------------------------------------------------------------
                IF(FLAGSWITCH.EQ.0) THEN 
                  ZETA1=DSQRT(RH10)*CDEXP((0.D0,1.D0)*PS)
                  ZETA2=DSQRT(R2)
                ELSE IF(FLAGSWITCH.EQ.1) THEN 
                  ZETA2=DSQRT(RH10)*CDEXP((0.D0,1.D0)*PS)
                  ZETA1=DSQRT(R2)
                END IF
                CALL EVAL_PHI(NMAX,ZETA1,ZETA2,Z1,Z2)
                CALL W_FLINES(-1,Z1,Z2,STAB)
C-------------------------------------------------------------------------------
C                Evaluation resonant hypervolume
C-------------------------------------------------------------------------------
              END IF
              IF(RH10.LE.AMPLMAX) THEN 
                HYPVOL=HYPVOL+PI*AMPLMAX/NSTEPT*
     .                         (PI*(AMPLMAX-R2)-RH10*3*SQRT(3.)/4.)
                IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Separatrix succesfully computed                   '
              ELSE
                IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Separatrix outside amplitude bounds               '
               END IF
 199          CONTINUE
            END DO
C-------------------------------------------------------------------------------
C                Normalization resonant hypervolume
C-------------------------------------------------------------------------------
            HYPVOL=HYPVOL/PI/PI/AMPLMAX/AMPLMAX*2
          ELSE IF((K11.EQ.2.AND.K12.EQ.1).OR.
     .               (K11.EQ.1.AND.K12.EQ.2)) THEN
C-------------------------------------------------------------------------------
C                Resonance [1,2] or [2,1]
C-------------------------------------------------------------------------------
            IF(K11.EQ.2.AND.K12.EQ.1) THEN
              FLAGSWITCH=1
              CALL SWITCH 
            END IF
            EPS=H10+2*H01
            RR=EPS*EPS/A1(0,0,1)/A1(0,0,1)/6
            DO R2=AMPLMAX,-2*AMPLMAX,-3*AMPLMAX/NSTEPT
C-------------------------------------------------------------------------------
C                First guess for the hyperbolic fixed point 
C-------------------------------------------------------------------------------
              IF(R2.GT.-3*RR.AND.R2.LT.RR) THEN
                IF(FLAGINIZ.EQ.0) THEN 
                  RH10=-R2/6+RR/3*(1+SQRT(1-R2/RR))
                  IF(A1(0,0,1)/EPS*(RH10+R2/6).GT.0) THEN
                    PS=PI-ALPHA1(0,0,1)
                  ELSE IF(A1(0,0,1)/EPS*(RH10+R2/6).LT.0) THEN
                    PS=-ALPHA1(0,0,1)
                  END IF
                END IF
C-------------------------------------------------------------------------------
C                Computation of the hyperbolic fixed points (arbitrary order) 
C-------------------------------------------------------------------------------
                DO JNEWT=1,MAXNEWT
                  IF(RH10.GT.OVRFLW.OR.RH10.LE.0.OR.
     .               PQ*RH10+R2.LE.0) THEN
                     IF(IFLAG.EQ.1) WRITE(34,992)
     .                   R2,'Rh10 out of bounds     ',RH10,R2
                    FLAGINIZ=0
                    GOTO 299
                  END IF
                  DET=DHC2D20(PS,RH10)*DHC2D02(PS,RH10)-
     .              DHC2D11(PS,RH10)**2
                  DRH10=(-DHC2D20(PS,RH10)*DHC2D01(PS,RH10)+
     .              DHC2D11(PS,RH10)*DHC2D10(PS,RH10))/DET
                  DPS=(DHC2D11(PS,RH10)*DHC2D01(PS,RH10)-
     .              DHC2D02(PS,RH10)*DHC2D10(PS,RH10))/DET
                  IF(ABS(DRH10).LT.EPSNEWT.AND.ABS(DPS).LT.EPSNEWT
     .              ) GOTO 300
                  RH10=RH10+DRH10
                  PS=PS+DPS
                END DO
                IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Divergent Newton for hyperbolic fixed lines       '
                FLAGINIZ=0
 300            CONTINUE
C-------------------------------------------------------------------------------
C                Eigenvalues (arbitrary order)
C-------------------------------------------------------------------------------
                STAB=DHC2D02(PS,RH10)*DHC2D20(PS,RH10)-
     .            DHC2D11(PS,RH10)**2
                IF(STAB.GT.0) THEN
                  IF(IFLAG.EQ.1) WRITE(34,992)
     .                   R2,'Wrong stability        ',STAB
                  FLAGINIZ=0
                  GOTO 299
                END IF
                FLAGINIZ=1
C-------------------------------------------------------------------------------
C                Evaluation of the whole separatrix
C-------------------------------------------------------------------------------
                RHMIN=RH10
                ESEP=HC2(PS,RH10)
                RHOAVE=0
                DO PSI=0,2*PI-EPSNEWT,PI/20
                  DO RHO=EPSNEWT,RH10*1.5,RH10/20
                    IF(PQ*RHO+R2.GT.0) THEN
                      ENER=HC2(PS+PI,RHO)
                      IF(ESEP.LT.0.AND.ENER.LT.ESEP) GOTO 221
                      IF(ESEP.GT.0.AND.ENER.GT.ESEP) GOTO 221
                    END IF
                  END DO
                  GOTO 299
 221              CONTINUE
C-------------------------------------------------------------------------------
C                Newton for the whole separatrix
C-------------------------------------------------------------------------------
                  DO JNEWT=1,MAXNEWT
                    IF(RHO.GT.OVRFLW.OR.RHO.LE.0) THEN
                      IF(IFLAG.EQ.1) WRITE(34,992)
     .                  R2,'Separatrix out of bound',RHO,R2
                      GOTO 299
                    END IF
                    DRHO=-(HC2(PSI,RHO)-ESEP)/DHC2D01(PSI,RHO)
                    IF(ABS(DRHO).LT.EPSNEWT) GOTO 222
                    RHO=RHO+DRHO
                  END DO
                  IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .              'Divergent Newton for separatrix                   '
                  GOTO 299
 222              CONTINUE
                  RHOAVE=RHOAVE+RHO/20
                  IF(RHO.LT.RHMIN) RHMIN=RHO
 223              CONTINUE
                END DO
                IF(IFLAG.EQ.1) THEN
C-------------------------------------------------------------------------------
C                Output on for033.dat of resonant patterns
C-------------------------------------------------------------------------------
                  IF(FLAGSWITCH.EQ.0) WRITE(33,991) 
     .               0.,0.,RHMIN,2*RHMIN+R2,RH10,2*RH10+R2,K11,K12 
                  IF(FLAGSWITCH.EQ.1) WRITE(33,991) 
     .               0.,0.,2*RHMIN+R2,RHMIN,0.,2*RH10+R2,RH10,K12,K11
C-------------------------------------------------------------------------------
C                Fixed points data output
C-------------------------------------------------------------------------------
                  IF(FLAGSWITCH.EQ.0) THEN 
                    ZETA1=DSQRT(RH10)*CDEXP((0.D0,1.D0)*PS)
                    ZETA2=DSQRT(2*RH10+R2)
                  ELSE IF(FLAGSWITCH.EQ.1) THEN 
                    ZETA2=DSQRT(RH10)*CDEXP((0.D0,1.D0)*PS)
                    ZETA1=DSQRT(2*RH10+R2)
                  END IF
                  CALL EVAL_PHI(NMAX,ZETA1,ZETA2,Z1,Z2)
                  CALL W_FLINES(-1,Z1,Z2,STAB)
C-------------------------------------------------------------------------------
C                Hypervolume evaluation
C-------------------------------------------------------------------------------
                END IF
                IF(RHOAVE.LT.(AMPLMAX-R2)/3) THEN
                  IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Separatrix succesfully computed                   '
                  IF(R2.GE.0) THEN
                    HYPVOL=HYPVOL+PI*(RHOAVE)
                  ELSE IF(R2.LT.0) THEN
                    HYPVOL=HYPVOL+PI*(RHOAVE+R2/2)
                  END IF
                ELSE
                  IF(R2.GE.0) THEN
                    HYPVOL=HYPVOL+PI*((AMPLMAX-R2)/3)
                  ELSE IF(R2.LT.0) THEN
                    HYPVOL=HYPVOL+PI*((AMPLMAX-R2)/3+R2/2)
                  END IF
                  IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Separatrix outside amplitude bounds               '
                END IF
C-------------------------------------------------------------------------------
 299            CONTINUE
              END IF
            END DO
            HYPVOL=1-HYPVOL*3/NSTEPT/PI/AMPLMAX*2
          ELSE IF((K11.EQ.2.AND.K12.EQ.-1).OR.
     .               (K11.EQ.1.AND.K12.EQ.-2)) THEN
C-------------------------------------------------------------------------------
C                Resonance [2,-1] or [1,-2]
C-------------------------------------------------------------------------------
          END IF
 777      CONTINUE
C-------------------------------------------------------------------------------
C                End resonance three
C-------------------------------------------------------------------------------
        ELSE IF(KRES.EQ.4) THEN
C-------------------------------------------------------------------------------
C                Resonances excited at order four
C-------------------------------------------------------------------------------
          IF(K11.EQ.4.AND.K12.EQ.0) THEN
C-------------------------------------------------------------------------------
C                Resonance [4,0]
C-------------------------------------------------------------------------------
            DO R2=EPSNEWT,AMPLMAX+2*EPSNEWT,AMPLMAX/NSTEPT
C-------------------------------------------------------------------------------
C                Selection stable/unstable resonance
C-------------------------------------------------------------------------------
              IF(ABS(A1(0,0,1)).GT.ABS(H20)) THEN
C-------------------------------------------------------------------------------
C                First guess for the hyperbolic fixed point 
C-------------------------------------------------------------------------------
                IF(FLAGINIZ.EQ.0) THEN 
                  IF((H10+R2*H11)*(H20+A1(0,0,1)).LT.0) THEN
                    RH10=-(H10+R2*H11)/(H20+A1(0,0,1))/2
                    PS=-ALPHA1(0,0,1)/4              
                  ELSE 
                    RH10=-(H10+R2*H11)/(H20-A1(0,0,1))/2
                    PS=(PI-ALPHA1(0,0,1))/4              
                  END IF
                END IF
C-------------------------------------------------------------------------------
C                Computation of the hyperbolic fixed point (arbitrary order) 
C-------------------------------------------------------------------------------
                DO JNEWT=1,MAXNEWT
                  IF(RH10.GT.OVRFLW.OR.RH10.LT.0) THEN
                    IF(IFLAG.EQ.1) WRITE(38,992) 
     .                             R2,'RH10=',RH10,R2,JNEWT
                    FLAGINIZ=0
                    GOTO 499
                  END IF
                  DET=DHC2D20(PS,RH10)*DHC2D02(PS,RH10)-
     .              DHC2D11(PS,RH10)**2
                  DRH10=(-DHC2D20(PS,RH10)*DHC2D01(PS,RH10)+
     .              DHC2D11(PS,RH10)*DHC2D10(PS,RH10))/DET
                  DPS=(DHC2D11(PS,RH10)*DHC2D01(PS,RH10)-
     .              DHC2D02(PS,RH10)*DHC2D10(PS,RH10))/DET
                  IF(ABS(DRH10).LT.EPSNEWT.AND.ABS(DPS).LT.EPSNEWT
     .              ) GOTO 400
                  RH10=RH10+DRH10
                  PS=PS+DPS
                END DO
                IF(IFLAG.EQ.1) WRITE(39,997) 
     .                     R2,'Divergent Newton for hyp. fixed lines'
                FLAGINIZ=0
 400            CONTINUE
C-------------------------------------------------------------------------------
C                Eigenvalues (arbitrary order)
C-------------------------------------------------------------------------------
                STAB=DHC2D02(PS,RH10)*DHC2D20(PS,RH10)-
     .            DHC2D11(PS,RH10)**2
                IF(STAB.GT.0) THEN
                  IF(IFLAG.EQ.1) WRITE(*,*) 
     .                 '***ERROR: WRONG EIGENVALUES IN [4,0]'
                  FLAGINIZ=0
                  GOTO 499
                END IF
                FLAGINIZ=1
C-------------------------------------------------------------------------------
C                Evaluation of the whole separatrix
C-------------------------------------------------------------------------------
                RHMIN=RH10
                ESEP=HC2(PS,RH10)
                DO PSI=PS+PI/100,PS+2*PI/4.d0,PI/50
                  RHO=RH10*(DCOS(PI/4-ABS(PSI-PS-PI/4)))**2
                  DO JNEWT=1,MAXNEWT
                    IF(RHO.GT.OVRFLW.OR.RHO.LE.0) THEN
                      IF(IFLAG.EQ.1) WRITE(38,992) 
     .                               R2,'RHO=',RHO,R2,JNEWT
                      GOTO 423
                    END IF
                    DRHO=-(HC2(PSI,RHO)-ESEP)/DHC2D01(PSI,RHO)
                    IF(ABS(DRHO).LT.EPSNEWT) GOTO 422
                    RHO=RHO+DRHO
                  END DO
                  IF(IFLAG.EQ.1) WRITE(39,997) 
     .                     R2,'Divergent Newton for separatrix'
                  GOTO 499
 422              CONTINUE
                  IF(RHO.LT.RHMIN) RHMIN=RHO
 423              CONTINUE
                END DO
                IF(IFLAG.EQ.1) THEN
C-------------------------------------------------------------------------------
C                Output on for033.dat of resonant patterns
C-------------------------------------------------------------------------------
                  IF(FLAGSWITCH.EQ.0) THEN 
                    WRITE(33,991) 0.,0.,RHMIN,R2,RH10,R2,K11,K12
                  ELSE IF(FLAGSWITCH.EQ.1) THEN 
                    WRITE(33,991) 0.,0.,R2,RHMIN,R2,RH10,K12,K11
                  END IF
C-------------------------------------------------------------------------------
C                Fixed points data output
C-------------------------------------------------------------------------------
                  IF(FLAGSWITCH.EQ.0) THEN 
                    ZETA1=DSQRT(RH10)*CDEXP((0.D0,1.D0)*PS)
                    ZETA2=DSQRT(R2)
                  ELSE IF(FLAGSWITCH.EQ.1) THEN 
                    ZETA2=DSQRT(RH10)*CDEXP((0.D0,1.D0)*PS)
                    ZETA1=DSQRT(R2)
                  END IF
                  CALL EVAL_PHI(NMAX,ZETA1,ZETA2,Z1,Z2)
                  CALL W_FLINES(-1,Z1,Z2,STAB)
C-------------------------------------------------------------------------------
                END IF
 499            CONTINUE
              END IF
            END DO
            HYPVOL=1
          END IF
        ELSE
C-------------------------------------------------------------------------------
C                Resonances excited at order higher than four
C-------------------------------------------------------------------------------
          DET=4*H20*H02-H11*H11
          IF(ABS(DET).LT.EPSNEWT) THEN
            WRITE(*,*) 
     .        '***ERROR(FIXLIN): DEGENERATE HAMILTONIAN'
            RETURN
          END IF
C-------------------------------------------------------------------------------
C                Loop on the resonant tune: tune is scanned \pm 0.2 
C                                           from the resonant tune 
C-------------------------------------------------------------------------------
          DO T=-0.2,0.2,STEPT      
C-------------------------------------------------------------------------------
C                Computes rh10,rh20 (first order)
C-------------------------------------------------------------------------------
            RH10=(2*H02*(-K12*P2Q2*T-H10)-H11*(K11*P2Q2*T-H01))/DET
            RH20=(-H11*(-K12*P2Q2*T-H10)+2*H20*(K11*P2Q2*T-H01))/DET
C-------------------------------------------------------------------------------
C                Computes rh10,rh20 (arbitrary order)
C-------------------------------------------------------------------------------
            DO JNEWT=1,MAXNEWT
              DET=DHC0D20(RH10,RH20)*DHC0D02(RH10,RH20)-
     .             DHC0D11(RH10,RH20)**2
              F1=DHC0D10(RH10,RH20)+K12*P2Q2*T
              F2=DHC0D01(RH10,RH20)-K11*P2Q2*T
              DRH10=-(DHC0D02(RH10,RH20)*F1-
     .               DHC0D11(RH10,RH20)*F2)/DET
              DRH20=-(-DHC0D11(RH10,RH20)*F1+
     .               DHC0D20(RH10,RH20)*F2)/DET
              IF(ABS(DRH10).LT.EPSNEWT.AND.ABS(DRH20).LT.EPSNEWT
     .          ) GOTO 101
                RH10=RH10+DRH10
                RH20=RH20+DRH20
            END DO
            IF(IFLAG.EQ.1) WRITE(34,997) R2, 
     .            'Divergent Newton for rh10,rh20                    '
            GOTO 88
 101        CONTINUE
            R2OLD=R2
            R2=-PQ*RH10+RH20
 997        FORMAT(1H ,'R2= ',D12.4,2X,A50)
C-------------------------------------------------------------------------------
C                Selection of virtual or far resonances
C-------------------------------------------------------------------------------
            IF(RH10.LT.0.OR.RH20.LT.0) THEN
            ELSE IF(RH10+RH20.GT.AMPLMAX) THEN
 992        FORMAT(1H ,'R2= ',D12.4,2X,A23,E13.5,2X,E13.5,I3)
C-------------------------------------------------------------------------------
C                Selection of real resonances
C-------------------------------------------------------------------------------
            ELSE 
              RESTYP='R'
C-------------------------------------------------------------------------------
C                Resonances not excited 
C-------------------------------------------------------------------------------
              IF(KRES.GT.NMAXH) THEN
                IF(FLAGSWITCH.EQ.0) THEN 
                  IF(IFLAG.EQ.1) WRITE(33,991) RH10,PQ*RH10+R2,
     .                  RH10,PQ*RH10+R2,RH10,PQ*RH10+R2,K11,K12
                ELSE IF(FLAGSWITCH.EQ.1) THEN 
                  IF(IFLAG.EQ.1) WRITE(33,991) PQ*RH10+R2,RH10,
     .                  PQ*RH10+R2,RH10,PQ*RH10+R2,RH10,K12,K11
                END IF
                HYPVOL=0            
 991            FORMAT(1X,6E13.5,2I5,3X,2I3)
C-------------------------------------------------------------------------------
C                Fixed point guess for resonances excited at first order
C-------------------------------------------------------------------------------
              ELSE
                IF(KRES.EQ.KSUM) THEN
                  KAPPA=-ALPHA1(0,0,1)
                  PS1P=KAPPA/K11
                  PS1M=(PI+KAPPA)/K11
C-------------------------------------------------------------------------------
C                Fixed point guess for resonances excited at higher order
C-------------------------------------------------------------------------------
                ELSE IF(KRES.LE.NMAXH.AND.KRES.GT.KSUM) THEN
                  SIGNOLD=DHC1D10(0D0,RH10)
                  DO PS=0,2.D0*PI/K11,2.D0*PI/K11/NTH
                    SIGN=DHC1D10(PS,RH10)
                    IF(SIGN*SIGNOLD.LT.0) GOTO 120
                    SIGNOLD=SIGN
                  END DO
                  IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Wrong first order guess for fixed line angles     '
                  GOTO 88
 120              PS1P=PS
                  DO JNEWT=1,MAXNEWT
                    DPS1P=-DHC1D10(PS1P,RH10)/DHC1D20(PS1P,RH10)
                    IF(ABS(DPS1P).LT.EPSNEWT) GOTO 130
                    PS1P=PS1P+DPS1P
                  END DO
                  IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Divergent Newton for first order fixed line angles'
                  GOTO 88
 130              CONTINUE
                  SIGNOLD=DHC1D10(PS1P+EPSNEWT,RH10)
                  DO PS=PS1P+EPSNEWT,2.D0*PI/K11,2.D0*PI/K11/NTH
                    SIGN=DHC1D10(PS,RH10)
                    IF(SIGN*SIGNOLD.LT.0) GOTO 140
                    SIGNOLD=SIGN
                  END DO
                  IF(IFLAG.EQ.1) WRITE(34,997) 
     .            'Wrong first order guess for fixed line angles     '
                  GOTO 88
 140              PS1M=PS
                  DO JNEWT=1,MAXNEWT
                    DPS1M=-DHC1D10(PS1M,RH10)/DHC1D20(PS1M,RH10)
                    IF(ABS(DPS1M).LT.EPSNEWT) GOTO 150
                    PS1M=PS1M+DPS1M
                  END DO
                  IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Divergent Newton for first order fixed line angles'
                  GOTO 88
 150              CONTINUE
                  IF(DHC1D20(PS1P,RH10)*DHC1D20(PS1M,RH10).GE.0) THEN
                    IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Wrong first order guess for fixed line angles     '
                    GOTO 88
                  END IF                   
                END IF
C-------------------------------------------------------------------------------
C                Fixed points guess of the amplitude
C-------------------------------------------------------------------------------
                RH1P=RH10
                RH1M=RH10
C-------------------------------------------------------------------------------
C                Computation of the fixed points (arbitrary order) p solution
C-------------------------------------------------------------------------------
                DO JNEWT=1,MAXNEWT
                  IF(RH1P.LE.0) THEN
                    IF(IFLAG.EQ.1) WRITE(34,992) 
     .                   R2,'Rh1p negative          ',RH1P,PQ*RH1P+R2
                    GOTO 88              
                  ELSE IF(RH1P.GT.OVRFLW.OR.PQ*RH1P+R2.LE.0) THEN
                    IF(IFLAG.EQ.1) WRITE(34,992) 
     .                   R2,'Rh1p out of bounds     ',RH1P,PQ*RH1P+R2
                    GOTO 88
                  END IF
                  DETP=DHC2D20(PS1P,RH1P)*DHC2D02(PS1P,RH1P)-
     .              DHC2D11(PS1P,RH1P)**2
                  DRH1P=(-DHC2D20(PS1P,RH1P)*DHC2D01(PS1P,RH1P)+
     .              DHC2D11(PS1P,RH1P)*DHC2D10(PS1P,RH1P))/DETP
                  DPS1P=(DHC2D11(PS1P,RH1P)*DHC2D01(PS1P,RH1P)-
     .              DHC2D02(PS1P,RH1P)*DHC2D10(PS1P,RH1P))/DETP
                  IF(ABS(DRH1P).LT.EPSNEWT.AND.ABS(DPS1P).LT.EPSNEWT
     .              ) GOTO 100
                  RH1P=RH1P+DRH1P
                  PS1P=PS1P+DPS1P
                END DO
                IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Divergent Newton for + fixed lines                '
                GOTO 88
 100            CONTINUE
C-------------------------------------------------------------------------------
C                Computation of the fixed points (arbitrary order) m solution
C-------------------------------------------------------------------------------
                DO JNEWT=1,MAXNEWT
                  IF(RH1M.LE.0) THEN
                    IF(IFLAG.EQ.1) WRITE(34,992) 
     .                   R2,'Rh1m negative          ',RH1M,PQ*RH1M+R2
                    GOTO 88             
                  ELSE IF(RH1M.GT.OVRFLW.OR.PQ*RH1M+R2.LE.0) THEN
                    IF(IFLAG.EQ.1) WRITE(34,992) 
     .                   R2,'Rh1m out of bounds     ',RH1M,PQ*RH1M+R2
                    GOTO 88
                  END IF
                  DETM=DHC2D20(PS1M,RH1M)*DHC2D02(PS1M,RH1M)-
     .              DHC2D11(PS1M,RH1M)**2
                  DRH1M=(-DHC2D20(PS1M,RH1M)*DHC2D01(PS1M,RH1M)+
     .              DHC2D11(PS1M,RH1M)*DHC2D10(PS1M,RH1M))/DETM
                  DPS1M=(DHC2D11(PS1M,RH1M)*DHC2D01(PS1M,RH1M)-
     .              DHC2D02(PS1M,RH1M)*DHC2D10(PS1M,RH1M))/DETM
                  IF(ABS(DRH1M).LT.EPSNEWT.AND.ABS(DPS1M).LT.EPSNEWT
     .              ) GOTO 102
                  RH1M=RH1M+DRH1M
                  PS1M=PS1M+DPS1M
                END DO
                IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Divergent Newton for - fixed lines                '
                GOTO 88
 102            CONTINUE
C-------------------------------------------------------------------------------
C                Eigenvalues (arbitrary order)
C-------------------------------------------------------------------------------
                STABP=DHC2D02(PS1P,RH1P)*DHC2D20(PS1P,RH1P)-
     .            DHC2D11(PS1P,RH1P)**2
                IF(STABP.GT.0) ISTABP=1
                IF(STABP.LT.0) ISTABP=-1
                STABM=DHC2D02(PS1M,RH1M)*DHC2D20(PS1M,RH1M)-
     .            DHC2D11(PS1M,RH1M)**2
                IF(STABM.GT.0) ISTABM=1
                IF(STABM.LT.0) ISTABM=-1
                IF(IFLAG.EQ.1) THEN
C-------------------------------------------------------------------------------
C                Fixed points data output
C-------------------------------------------------------------------------------
                  IF(FLAGSWITCH.EQ.0) THEN 
                    ZETA1=DSQRT(RH1P)*CDEXP((0.D0,1.D0)*PS1P)
                    ZETA2=DSQRT(PQ*RH1P+R2)
                  ELSE IF(FLAGSWITCH.EQ.1) THEN 
                    ZETA2=DSQRT(RH1P)*CDEXP((0.D0,1.D0)*PS1P)
                    ZETA1=DSQRT(PQ*RH1P+R2)
                  END IF
                  CALL EVAL_PHI(NMAX,ZETA1,ZETA2,Z1,Z2)
                  CALL W_FLINES(ISTABP,Z1,Z2,STABP)
                  IF(FLAGSWITCH.EQ.0) THEN 
                    ZETA1=DSQRT(RH1M)*CDEXP((0.D0,1.D0)*PS1M)
                    ZETA2=DSQRT(PQ*RH1M+R2)
                  ELSE IF(FLAGSWITCH.EQ.1) THEN 
                    ZETA2=DSQRT(RH1M)*CDEXP((0.D0,1.D0)*PS1M)
                    ZETA1=DSQRT(PQ*RH1M+R2)
                  END IF
                  CALL EVAL_PHI(NMAX,ZETA1,ZETA2,Z1,Z2)
                  CALL W_FLINES(ISTABM,Z1,Z2,STABM)
C-------------------------------------------------------------------------------
C                Stability analysis
C-------------------------------------------------------------------------------
                END IF
                IF(STABP.GT.0.AND.STABM.LT.0) THEN
                  PSEL=PS1P
                  PSHY=PS1M
                  SEPEN=HC2(PS1M,RH1M)
                ELSE IF(STABP.LT.0.AND.STABM.GT.0) THEN
                  PSEL=PS1M
                  PSHY=PS1P
                  SEPEN=HC2(PS1P,RH1P)
                ELSE
                  IF(IFLAG.EQ.1) WRITE(34,992) 
     .                   r2,'Wrong stability        ',STABP,STABM
                  GOTO 88
                END IF
C-------------------------------------------------------------------------------
C                Computation first guess for resonances excited at first order
C-------------------------------------------------------------------------------
                IF(KRES.EQ.KSUM) THEN 
                  C=2*A1(0,0,1)*RH10**(K11/2.D0)*
     .                    (PQ*RH10+R2)**(ABS(K12)/2.D0)
C-------------------------------------------------------------------------------
C                Computation first guess for resonances excited at higher orders
C-------------------------------------------------------------------------------
                ELSE IF(KRES.GE.KSUM.AND.KRES.LE.NMAXH) THEN
                  C=HC1(PSEL,RH10)-HC1(PSHY,RH10)
                END IF
                DEL=SQRT(ABS(2*C/DHC3D2(RH10)))
                RHIN=RH10-DEL
                RHOUT=RH10+DEL
                IF(RHIN.LT.0) THEN 
                  RHIN=EPSNEWT
                END IF
C-------------------------------------------------------------------------------
C                Newton for the island width: inner part
C-------------------------------------------------------------------------------
                DO JNEWT=1,MAXNEWT
                  IF(RHIN.LT.0) THEN
                    IF(IFLAG.EQ.1) WRITE(34,992) 
     .                   R2,'Inner separatrix < 0   ',RHIN,PQ*RHIN+R2
                    GOTO 88             
                  ELSE IF(RHIN.GT.OVRFLW.OR.PQ*RHIN+R2.LE.0) THEN
                    IF(IFLAG.EQ.1) WRITE(34,992) 
     .                   R2,'Inner separatrix o.o.b.',RHIN,PQ*RHIN+R2
                    GOTO 88
                  END IF
                  DRHIN=-(HC2(PSEL,RHIN)-SEPEN)/DHC2D01(PSEL,RHIN)
                  IF(ABS(DRHIN).LT.EPSNEWT) GOTO 90
                  RHIN=RHIN+DRHIN
                END DO
                IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Divergent Newton for the innner width             '
                GOTO 88
 90             CONTINUE
                RHINST=RHIN
C-------------------------------------------------------------------------------
C                Newton for the island width: outer part
C-------------------------------------------------------------------------------
                DO JNEWT=1,MAXNEWT
                  IF(RHOUT.LT.0) THEN
                    IF(IFLAG.EQ.1) WRITE(34,992) 
     .                   R2,'Outer separatrix < 0   ',RHOUT,PQ*RHOUT+R2
                    GOTO 88              
                  ELSE IF(RHOUT.GT.OVRFLW.OR.PQ*RHOUT+R2.LE.0)THEN
                    IF(IFLAG.EQ.1) WRITE(34,992) 
     .                   R2,'Outer separatrix o.o.b.',RHOUT,PQ*RHOUT+R2
                    GOTO 88
                  END IF
                  DRHOUT=-(HC2(PSEL,RHOUT)-SEPEN)/DHC2D01(PSEL,RHOUT)
                  IF(ABS(DRHOUT).LT.EPSNEWT) GOTO 91
                  RHOUT=RHOUT+DRHOUT
                END DO
                IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Divergent Newton for the outer width              '
                GOTO 88
 91             CONTINUE
                RHOUTST=RHOUT
C-------------------------------------------------------------------------------
C                Output on for033.dat of resonant patterns
C-------------------------------------------------------------------------------
                WIDTH=RHOUTST-RHINST
                IF(FLAGSWITCH.EQ.0) THEN 
                  IF(IFLAG.EQ.1) WRITE(33,991) RH10,PQ*RH10+R2,
     .                RHINST,PQ*RHINST+R2,RHOUTST,PQ*RHOUTST+R2,K11,K12
                ELSE IF(FLAGSWITCH.EQ.1) THEN 
                  IF(IFLAG.EQ.1) WRITE(33,991) PQ*RH10+R2,RH10,
     .                PQ*RHINST+R2,RHINST,PQ*RHOUTST+R2,RHOUTST,K12,K11
                END IF
C-------------------------------------------------------------------------------
C                Computation of the resonant hypervolume
C-------------------------------------------------------------------------------
                IF(RH10*(1+PQ)+R2.LT.AMPLMAX) THEN
                  HYPVOL=HYPVOL+WIDTH*ABS(R2-R2OLD)
                  IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Island width succesfully computed                 '
                ELSE
                  IF(IFLAG.EQ.1) WRITE(34,997) R2,
     .            'Island width outside amplitude bounds             '
                END IF
              END IF
            END IF
 88         CONTINUE
          END DO
          HYPVOL=8*HYPVOL*PI
C-------------------------------------------------------------------------------
C                normalization of the resonant hypervolume
C-------------------------------------------------------------------------------
          HYPVOL=HYPVOL/(AMPLMAX*AMPLMAX*2*PI*PI)
C-------------------------------------------------------------------------------
C          
C-------------------------------------------------------------------------------
        END IF
        IF(RESTYP.EQ.'R'.AND.HYPVOL.EQ.0.AND.KRES.LE.NMAXH) 
     .    HYPVOL=0.999
C-------------------------------------------------------------------------------
C                        Output of island hypervolume
C-------------------------------------------------------------------------------
        IF(IFLAG.EQ.1) 
     .    WRITE(31,'(5x,''Island Hypervol.  '',''q4='',D12.5)') HYPVOL
      END IF

      END 
CDECK  ID>, PRINTUNE.
C===============================================================================
C  PRINTS OUT ON UNIT IUNIT THE TUNESHIFT COEFFICIENTS, I.E. THE DERIVATIVE
c  OF THE INTERPOLATING HAMILTONIAN OF THE NONRESONANT CASE WITH RESPECT TO 
C  THE AMPLITUDES, DIVIDED BY TWO PI. 
C  THE HAMILTONIAN MUST HAVE BEEN PREVIOUSLY COMPUTED (SUBROUTINE HAZANFILL)
C  AND STORED IN THE ARRAY A0, ALPHA0 IN THE COMMON BLOCK /HAMILT/
C
C  AUTHOR: E. TODESCO - INFN
C
C

      SUBROUTINE PRINTUNE(IUNIT)

      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 PI
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES

      PI=4*DATAN(1.D0)
 663  FORMAT(1H ,D24.16,4X,'(',2I2,')')
      WRITE(IUNIT,*) NMAXH
C-------------------------------------------------------------------------------
C                                  TUNE X
C-------------------------------------------------------------------------------
      WRITE(IUNIT,*) 
      WRITE(IUNIT,*) 'Tune x'
      DO N=2,NMAXH,2
        DO K1=1,N/2
          K2=N/2-K1
          WRITE(IUNIT,663) K1*A0(K1,K2)*DCOS(ALPHA0(K1,K2))/2/PI,K1-1,K2
        END DO
      END DO
C-------------------------------------------------------------------------------
C                                  TUNE Y
C-------------------------------------------------------------------------------
      WRITE(IUNIT,*) 
      WRITE(IUNIT,*) 'Tune y'
      DO N=2,NMAXH,2
        DO K1=0,N/2-1
          K2=N/2-K1
          WRITE(IUNIT,663) K2*A0(K1,K2)*DCOS(ALPHA0(K1,K2))/2/PI,K1,K2-1
        END DO
      END DO
      END
CDECK  ID>, Q1.
C=============================================================
C  COMPUTES THE NORM OF THE MAP (QUALITY FACTOR Q1) AT THE 
c  AMPLITUDE A USING THE MAP IN REAL COURANT-SNYDER 
c  COORDINATES STORED IN RMAP THE SUM OF THE ABSOLUTE VALUE OF 
c  THE COEFFICIENTS OF ORDER JORD IS STORED IN THE VECTOR 
C  QMAP(JORD), COMMON BLOCK /QFACT/
C
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION Q1(A,IORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C..............................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
C..............................................................
      COMMON/QFACT/QMAP(NORMAX1),QTUNE(NORMAX1),
     .             QRES1(0:NORMAX1,-NORMAX1:NORMAX1,NORMAX1),
     .             IFLRES(0:NORMAX1,-NORMAX1:NORMAX1),
     .             QRES2(0:NORMAX1,-NORMAX1:NORMAX1)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
C..............................................................
      Q1=0D0
      POT=A
      DO JORD=2,IORD
C..............................................................
        QMAP(JORD)=0D0
        DO I1=0,JORD
          DO I2=0,JORD-I1
            DO I3=0,JORD-I1-I2
              I4=JORD-I1-I2-I3
C..............................................................
              QMAP(JORD)=QMAP(JORD)+DABS(RMAP(IPUNT(I4,I3,I2,I1),1))
     .                             +DABS(RMAP(IPUNT(I4,I3,I2,I1),2))
     .                             +DABS(RMAP(IPUNT(I4,I3,I2,I1),3))
     .                             +DABS(RMAP(IPUNT(I4,I3,I2,I1),4))
C..............................................................
            END DO
          END DO
        END DO
C..............................................................
        POT=POT*A
        QMAP(JORD)=QMAP(JORD)*POT
        Q1=Q1+QMAP(JORD)
      END DO
      QMAP(1)=Q1
C..............................................................
      END
CDECK  ID>, Q2I.
C==============================================================
C  COMPUTES THE NORM OF THE TUNESHIFT (QUALITY FACTOR Q2) AT 
C  THE AMPLITUDE A USING THE INTERPOLATING HAMILTONIAN A0 
C  STORED IN THE COMMON BLOCK HAMILT.
C  IN QTUNE THE NORM OF THE TUNESHIFT COMPUTED AT LOWER ORDERS 
C  IS STORED:
C  QTUNE(4)   FIRST ORDER TUNESHIFT
C  QTUNE(6)   FIRST + SECOND ORDER TUNESHIFT
C  QTUNE(8)   FIRST + SECOND + THIRD ORDER TUNESHIFT
C  ....
C
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION Q2I(A)
C.............................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.............................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
C.............................................................
      COMMON/QFACT/QMAP(NORMAX1),QTUNE(NORMAX1),
     .             QRES1(0:NORMAX1,-NORMAX1:NORMAX1,NORMAX1),
     .             IFLRES(0:NORMAX1,-NORMAX1:NORMAX1),
     .             QRES2(0:NORMAX1,-NORMAX1:NORMAX1)
      COMMON/NLMAP/RMAT(NDIM,4),RKICK(NDIM,4),RMAP(NDIM,4),
     .             TEMP(NDIM,4),RDER(NDIM,16)
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C.............................................................
      PI2=8D0*DATAN(1D0)
      O1=DACOS(RMAP(5,1))/PI2
      O2=DACOS(RMAP(3,3))/PI2
C.............................................................
      DO JORD=4,NMAXH,2
        ACC=0D0
C.............................................................
        DO J=1,49
          RH1=A*J/50.D0
          RH2=A-RH1
          ACC=ACC+((TUNE1(RH1,RH2,JORD)-O1)*
     .             (TUNE1(RH1,RH2,JORD)-O1)
     .           + (TUNE2(RH1,RH2,JORD)-O2)*
     .             (TUNE2(RH1,RH2,JORD)-O2))/2
        END DO
C.............................................................
        QTUNE(JORD-1)=SQRT(ACC/49)
C.............................................................
      ENDDO
      QTUNE(1)=SQRT(ACC/49)
      Q2I=QTUNE(1)
C.............................................................
      RETURN
C.............................................................
      END
CDECK  ID>, Q3.
C=============================================================
C  COMPUTES THE NORM OF THE RESONANT PART OF THE HAMILTONIAN 
C  (QUALITY FACTOR Q3) AT THE AMPLITUDE A, ORDER IORD, USING 
C  THE INTERPOLATING HAMILTONIAN A1 STORED IN THE COMMON BLOCK 
C  HAMILT.
C  IN QRES1(K1,K2,N) ARE STORED THE COEFFICIENTS OF RESONANCE
C  K1,K2 AT ORDER N, COMMON BLOCK /QFACT/. 
C  IN QRES1(K1,K2,1) IS STORED THE SUM OF THE COEFFICIENTS OF 
C  RESONANCE K1,K2 ON DIFFERENT ORDERS.
C  IN QRES1(0,0,N) IS STORED THE SUM OF ALL THE RESONANT 
C  COEFFICIENTS OVER THE PERTURBATIVE ORDER.
C  IN QRES1(0,0,1) IS STORED THE GLOBAL SUM OVER RESONANCE AND
C  ORDER.
C
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION Q3(A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
C.............................................................
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
      COMMON/QFACT/QMAP(NORMAX1),QTUNE(NORMAX1),
     .             QRES1(0:NORMAX1,-NORMAX1:NORMAX1,NORMAX1),
     .             IFLRES(0:NORMAX1,-NORMAX1:NORMAX1),
     .             QRES2(0:NORMAX1,-NORMAX1:NORMAX1)
C.............................................................
      DO N=1,NORMAX1
        QRES1(K11,K12,N)=0D0
      END DO
      DO N=2,NMAXH
        SCRA=0
        DO L=0,N/(KSUM)
          DO K1=0,(N-L*KSUM)/2
            K2=(N-L*KSUM-2*K1)/2
            IF(K2*2.EQ.N-L*KSUM-2*K1) THEN
              IF(L.NE.0) QRES1(K11,K12,N)=QRES1(K11,K12,N)+A1(K1,K2,L)
            END IF
          END DO
        END DO
        QRES1(K11,K12,N)=QRES1(K11,K12,N)*A**(N/2.D0)
        QRES1(K11,K12,1)=QRES1(K11,K12,1)+QRES1(K11,K12,N)
      END DO
      Q3=QRES1(K11,K12,1)
      RETURN
      END 
CDECK  ID>, TUNE1.
C=============================================================
C  COMPUTES THE TUNE ON THE X,PX PLANE THROUGH NONRESONANT 
C  NORMAL FORMS AT THE NONLINEAR INVARIANTS RH1,RH2. THE 
C  INTERPOLATING HAMILTONIAN A0 IS USED.
C  THIS IS USED FOR EVALUATING THE ANALYTICAL TUNE FOOTPRINT.
C
C
C  AUTHOR: E. TODESCO - INFN
C

      REAL*8 FUNCTION TUNE1(RH1,RH2,JORD)
C..............................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 PI
C..............................................................
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C..............................................................
      PI=4*DATAN(1.D0)
      TUNE1=0 
      DO N=1,JORD/2
        DO K=0,N
          TUNE1=TUNE1+A0(K,N-K)*DCOS(ALPHA0(K,N-K))*
     .           K*RH1**(K-1)*RH2**(N-K)
        END DO
      END DO
      TUNE1=TUNE1/2/PI
      IF(TUNE1.LT.0) TUNE1=1+TUNE1
      IF(TUNE1.GT.1) TUNE1=TUNE1-1
C..............................................................
      END
CDECK  ID>, TUNE2.
C==============================================================
C  COMPUTES THE TUNE ON THE Y,PY PLANE THROUGH NONRESONANT 
C  NORMAL FORMS AT THE NONLINEAR INVARIANTS RH1,RH2. THE 
C  INTERPOLATING HAMILTONIAN A0 IS USED.
C  THIS IS USED FOR EVALUATING THE ANALYTICAL TUNE FOOTPRINT.
C
C
C  AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION TUNE2(RH1,RH2,JORD)
C..............................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 A0(0:NORMAX1,0:NORMAX1),A1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 A2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 ALPHA0(0:NORMAX1,0:NORMAX1),PQ,OR1,OR2
      REAL*8 ALPHA1(0:NORMAX1,0:NORMAX1,0:NORMAX1)
      REAL*8 ALPHA2(0:NORMAX1,0:NORMAX1,0:NORMAX1,-NORMAX1:NORMAX1)
      REAL*8 PI
C..............................................................
      COMMON/HAMILT/A0,A1,A2,ALPHA0,ALPHA1,ALPHA2,PQ,OR1,OR2
      COMMON/ORDER/NMAX,NMAXH,K11,K12,K21,K22,KSUM,IAK12,
     .             IQ1,IP1,IQ2,IP2,IAP1,IAP2,KRES
C..............................................................
      PI=4*DATAN(1.D0)
      TUNE2=0
      DO N=1,JORD/2
        DO K=0,N
          TUNE2=TUNE2+A0(K,N-K)*DCOS(ALPHA0(K,N-K))*
     .           RH1**K*(N-K)*RH2**(N-K-1)
        END DO
      END DO
      TUNE2=TUNE2/2/PI
      IF(TUNE2.LT.0) TUNE2=1+TUNE2
      IF(TUNE2.GT.1) TUNE2=TUNE2-1
C..............................................................
      END
CDECK  ID>, DIRES.
C============================================================
C COMPUTES THE DISTANCE OF A RESONANT LINE OF EQUATION 
C K1*X+K2*Y=INTEGER TO THE WORKING POINT (O1,O2)
C
C
C      AUTHOR: E. TODESCO - INFN
C

      DOUBLE PRECISION FUNCTION DIRES(O1,O2,K1,K2)
      DOUBLE PRECISION O1,O2,UGH,OR1,OR2
C.....................................SINGLE RESONANCE ON X,PX
      IF (K2.EQ.0) THEN
        OR1=DFLOAT(INT(O1*K1+0.5))/DFLOAT(K1)
        DIRES=ABS(OR1-O1)
C.....................................SINGLE RESONANCE ON Y,PY
      ELSE IF (K1.EQ.0) THEN
        OR2=DFLOAT(INT(O2*K2+0.5))/DFLOAT(K2)
        DIRES=ABS(OR2-O2)
C................................................SUM RESONANCE
      ELSE IF(K1*K2.GT.0) THEN
        DIRES=1000
        DO L=0,K1+K2
          UGH=ABS(K1*O1+K2*O2-L)/SQRT(FLOAT(K1*K1+K2*K2))
          IF(UGH.LT.DIRES) DIRES=UGH
        END DO
C.........................................DIFFERENCE RESONANCE
      ELSE IF(K1*K2.LT.0) THEN
        DIRES=1000
        DO L=K2,K1
          UGH=ABS(K1*O1+K2*O2-L)/SQRT(FLOAT(K1*K1+K2*K2))
          IF(UGH.LT.DIRES) DIRES=UGH
        END DO
      END IF
C.............................................................
      RETURN
C.............................................................
      END 
CDECK  ID>, PREP_NF_NONRES.
C=================================================================
C FILLS THE COMMON /RISONANZA/ WITH THE PARAMETERS NECESSARY TO 
C COMPUTE THE NONRESONANT NORMAL FORM, 
C
C AUTHOR: E. TODESCO - INFN - BOLOGNA
C
C
      SUBROUTINE PREP_NF_NONRES 
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      COMMON/RISONANZA/IRES1(2),IRES2(2),AUTORES(4),IRESON,ICASE

      ICASE=0
      IRESON=2
C......................RESONANT EIGENVALUES FOR THE LINEAR MAP
      AUTORES(1)=1
      AUTORES(2)=1
      AUTORES(3)=1
      AUTORES(4)=1
C..................BASIS FOR THE SPACE OF RESONANT FREQUENCIES
      IRES1(1)=0
      IRES2(1)=0
      IRES1(2)=0
      IRES2(2)=0
C.............................................................

      END 
CDECK  ID>, PREP_NF_SINGRES.
C=================================================================
C FILLS THE COMMON /RISONANZA/ WITH THE PARAMETERS NECESSARY TO 
C COMPUTE THE RESONANT NORMAL FORM, SINGLE RESONANCE (K1,K2)
C
C AUTHOR: E. TODESCO - INFN - BOLOGNA
C
C
      SUBROUTINE PREP_NF_SINGRES(O1,O2,K1,K2)
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      PARAMETER(MINEPS=1D-5)
      DOUBLE PRECISION PI,OR1,OR2,O1,O2,EPS,UGH
      COMMON/RISONANZA/IRES1(2),IRES2(2),AUTORES(4),IRESON,ICASE
C............................................................
      PI=4.D0*DATAN(1.D0)
      J11=0
      J12=0
      J21=0
      J22=0
      ICASE=1
C....................................SINGLE RESONANCE ON X,PX
      IF (K2.EQ.0) THEN
        OR1=DFLOAT(INT(O1*K1+0.5))/DFLOAT(K1)
        OR2=0.
        EPS=DABS(OR1-O1)
        J11=K1
C....................................SINGLE RESONANCE ON Y,PY
      ELSEIF (K1.EQ.0) THEN
        OR1=0.
        OR2=DFLOAT(INT(O2*K2+0.5))/DFLOAT(K2)
        EPS=DABS(OR2-O2)
        J12=K2
C...........................................COUPLED RESONANCE
      ELSE 
        IF(K1*K2.GT.0) LMIN=0
        IF(K1*K2.GT.0) LMAX=K1+K2
        IF(K1*K2.LT.0) LMIN=K2
        IF(K1*K2.LT.0) LMAX=K1
        EPS=1000
        DO L=LMIN,LMAX
          UGH=ABS(K1*O1+K2*O2-L)/SQRT(FLOAT(K1*K1+K2*K2))
          IF(UGH.LT.EPS) THEN
            EPS=UGH
            LL=L
          END IF
        END DO
        OR1=(LL*K1+K2*K2*O1-K1*K2*O2)/(K1*K1+K2*K2)
        OR2=(LL*K2+K1*K1*O2-K1*K2*O1)/(K1*K1+K2*K2)
        J11=K1
        J12=K2
      END IF
C.............................................................
      IF (DABS(EPS).LT.MINEPS) THEN
        IRESON=1
      ELSE
        IRESON=2
      ENDIF
C.............................................................
      ER1 = EXP( DCMPLX(0.D0,2*OR1*PI) )
      ER2 = EXP( DCMPLX(0.D0,2*OR2*PI) )
C......................RESONANT EIGENVALUES FOR THE LINEAR MAP
      AUTORES(1)=ER1
      AUTORES(2)=DCONJG(ER1)
      AUTORES(3)=ER2
      AUTORES(4)=DCONJG(ER2)
C..................BASIS FOR THE SPACE OF RESONANT FREQUENCIES
      IRES1(1)=J11
      IRES2(1)=J12
      IRES1(2)=J21
      IRES2(2)=J22
C.............................................................

      END 
CDECK  ID>, GENERALE_4D.
C=================================================================
C MAIN SUBROUTINE TO COMPUTE NORMAL FORMS.
C REALMAP IS A TRANSFER MAP STORED IN THE ARRAY REALMAP(NDIM,4).
C IT REPRESENTS A MAP IN REAL COURANT-SNYDER COORDINATES.
C MAXMAP IS THE ORDER OF THE TRANSFER MAP.
C MAXORD IS THE ORDER OF THE NORMAL FORM.
C INTERACT IS A FLAG TO SUPPRESS INTERACTIVE INPUT:
C INTERACT = 0 INTERACTIVE INPUT SUPPRESSED.
C INTERACT = 1 INTERACTIVE INPUT ENABLED.
C IINV IS A FLAG TO SUPPRESS THE COMPUTATION OF PSI IN ORDER TO SAVE
C TIME:
C IINV = 0 PSI COMPUTATION SUPPRESSED.
C IINV = 1 PSI COMPUTATION ENABLED.
C IOUT IS A FLAG FOR THE OUTPUT:
C IOUT = 0 NO OUTPUT.
C IOUT > 1 OUTPUT ENABLED.
C IN THE LAST CASE THE OUTPUT CONSISTS OF THREE OR FOUR FILES
C (ACCORDING TO THE VALUES OF IOUT & IINV):
C
C PHI.DAT    ===> CONJUGATION FUNCTION.
C                 IF IOUT = 1 IT IS EXPRESSED AS A REAL
C                 POLYNOMIAL IN REAL C-S COORDINATES
C                 IF IOUT = 2 IT IS EXPRESSED AS A REAL
C                 POLYNOMIAL IN REAL PHYSICAL COORDINATES
C                 IF IOUT = 3 IT IS EXPRESSED AS A COMPLEX
C                 POLYNOMIAL IN COMPLEX C-S COORDINATES
C PSI.DAT    ===> INVERSE OF THE CONJUGATION FUNCTION IN REAL
C                 COORDINATES (SEE PHI.DAT FOR MORE DETAILS ON
C                 THE CONVENTION)(IN THE CASE IINV = 1)
C HAMIL.DAT  ===> REAL HAMILTONIAN IN THE INVARIANTS.
C HAZANG.DAT ===> REAL HAMILTONIAN IN ACTION ANGLE VARIABLES.
C
C IAZAN IS A FLAG TO FILL THE COMMONS HAZAN* WITH THE HAMILTONIAN IN
C ACTION ANGLE VARIABLES (MORE DETAILS IN SUBROUTINE MHAZANFILL):
C IAZAN = O THE COMMON BLOCK IS LEFT EMPTY.
C IAZAN = 1 THE COMMON BLOCK IS FILLED.
C
C AUTHOR: A. BAZZANI (M. GIOVANNOZZI HAS CREATED A SUBROUTINE OUT OF
C IT).
C
 
      SUBROUTINE GENERALE_4D(REALMAP,MAXMAP,MAXORD,INTERACT,IINV,IOUT,
     .                       IAZAN)
      IMPLICIT COMPLEX*16 (A-H,O-Z)
C.............................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      PARAMETER(MINEPS=1D-5)
      DOUBLE PRECISION REALMAP(NDIM,4),O1,O2,OR1,OR2,OMER1,OMER2,EPS
      DOUBLE PRECISION DIST,DISTMIN,PI2
      DIMENSION H(NDIMPLUS)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/AUTOVALORI/AUTO(4)
      COMMON/RESTI1/RESTO1(NDIM,4)
      COMMON/RESTI2/RESTO2(NDIM,4)
      COMMON/ IDENTITY/ K(4,2)
      COMMON/CAMPO/ A(NDIM,4)
      COMMON/RISONANZA/IRES1(2),IRES2(2),AUTORES(4),IRESON,ICASE
      COMMON/DIRETTA/PHI1(NDIM,4)
      COMMON/INVERSA/PSI(NDIM,4)
      COMMON/MAPPA/F(NDIM,4)
      COMMON/SCRATCH1/Q(NDIM,4)
      COMMON/SCRATCH2/D1(NDIM)
      COMMON/SCRATCH3/R(NDIM,4)
      COMMON/SCRATCH4/D(NDIM,4)
C.............................................................
      IF (MAXMAP.LT.0.OR.MAXMAP.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(GENERALE_4D): SECOND PARAMETER OUT OF ',
     .             'BOUNDS'
        STOP
      ENDIF
      IF (MAXORD.LT.0.OR.MAXORD.GT.NORMAX) THEN
        WRITE(6,*) '***ERROR(GENERALE_4D): THIRD PARAMETER OUT OF ',
     .             'BOUNDS'
        STOP
      ENDIF
      IF (INTERACT.NE.0.AND.INTERACT.NE.1) THEN
        WRITE(6,*) '***ERROR(GENERALE_4D): FOURTH PARAMETER OUT OF ',
     .             'BOUNDS'
        STOP
      ENDIF
      IF (IINV.NE.0.AND.IINV.NE.1) THEN
        WRITE(6,*) '***ERROR(GENERALE_4D): FIVETH PARAMETER OUT OF ',
     .             'BOUNDS'
        STOP
      ENDIF
      IF (IOUT.LT.0.OR.IOUT.GT.3) THEN
        WRITE(6,*) '***ERROR(GENERALE_4D): SIXTH PARAMETER OUT OF ',
     .             'BOUNDS'
        STOP
      ENDIF
      IF (IAZAN.NE.0.AND.IAZAN.NE.1) THEN
        WRITE(6,*) '***ERROR(GENERALE_4D): SEVENTH PARAMETER OUT OF ',
     .             'BOUNDS'
        STOP
      ENDIF
C.............................................................
      NFOUR=4
      IDIME=NFOUR
      CALL MRTOC(REALMAP,F,MAXMAP,NFOUR) !..CONVERSION FROM R TO C
C................................COMPUTES THE LINEAR FREQUENCY
      PI2=8D0*DATAN(1D0)
      O1=DACOS(REALMAP(5,1))/PI2
      O2=DACOS(REALMAP(3,3))/PI2
C.............................................READS.PARAMETERS
      IF (INTERACT.EQ.1) THEN
C.............................................................
        IF (O1.LT.0) O1=1+O1
        IF (O2.LT.0) O2=1+O2
        WRITE(*,*) '0: NONRESONANT CASE'
        WRITE(*,*) '1: SINGLE OR COUPLED RESONANCE CASE'
        WRITE(*,*) '2: DOUBLE RESONANCE CASE'
        READ(*,*) ICASE
C.............................................................
        IF (ICASE.EQ.0) THEN
          CALL PREP_NF_NONRES
          WRITE(6,*) 'ICASE = ',ICASE,'IRESON = ',IRESON
C.............................................................
        ELSEIF (ICASE.EQ.1) THEN
C.............................................................
          WRITE(*,*) 'PLEASE INPUT RESONANCE K1,K2   (K1*NU1+K2*NU2=L)'
          READ(*,*) K1,K2
          CALL PREP_NF_SINGRES(O1,O2,K1,K2)
C.............................................................
          IF (IRESON.EQ.1) THEN
            WRITE(6,*) 'RESONANT CASE'
          ELSE
            WRITE(6,*) 'QUASI RESONANT CASE'
          ENDIF
          WRITE(6,*) 'ICASE = ',ICASE,'IRESON = ',IRESON
          WRITE(*,*) 'NON RESONANT FREQUENCIES',O1,O2
          WRITE(*,*) 'RESONANT FREQUENCIES',
     .                CDLOG(AUTORES(1))/PI2,CDLOG(AUTORES(3))/PI2
C.............................................................
        ELSEIF (ICASE.EQ.2) THEN
          WRITE(*,*) 'PLEASE INPUT P1,Q1,P2,Q2'
          READ(*,*) IP1,IQ1,IP2,IQ2
          OR1=DFLOAT(IP1)/DFLOAT(IQ1)
          OR2=DFLOAT(IP2)/DFLOAT(IQ2)
          IF ((DABS(OR1-O1).LT.MINEPS).AND.
     .       (DABS(OR2-O2).LT.MINEPS)) THEN
            IRESON=1
          ELSE
            IRESON=2
          ENDIF
C..................BASIS FOR THE SPACE OF RESONANT FREQUENCIES
          IRES1(1)=IQ1
          IRES2(1)=0
          IRES1(2)=0
          IRES2(2)=IQ2
          WRITE(6,*) 'ICASE = ',ICASE,'IRESON = ',IRESON
          WRITE(6,*) 'IP1 = ',IP1,'IQ1 = ',IQ1,
     .               'IP2 = ',IP2,'IQ2 = ',IQ2
C......................RESONANT EIGENVALUES FOR THE LINEAR MAP
          AUTORES(1)=EXP(DCMPLX(0.D0,OR1*PI2))
          AUTORES(2)=DCONJG(AUTORES(1))
          AUTORES(3)=EXP(DCMPLX(0.D0,OR2*PI2))
          AUTORES(4)=DCONJG(AUTORES(3))
C.............................................................
        ENDIF
C.............................................................
      ENDIF
C....................................................END INPUT
      CALL NORMAL_FORM(H,MAXMAP,MAXORD,IINV,IOUT)
C...............................FILLS THE HAZAN* COMMON BLOCKS
      IF (IAZAN.EQ.1) CALL MHAZANFILL(H,MAXORD+1,IDIME)
C..........................OUTPUT HAMILTONIAN (ZETA VARIABLES)
      IF (IOUT.GT.0) THEN
        IUNIT=12
        CALL MPRINTHC(H,MAXORD+1,IDIME,IUNIT)
C..................OUTPUT HAMILTONIAN (ACTION-ANGLE VARIABLES)
        IUNIT=13
        CALL MPRINTAZAN(IUNIT)
      ENDIF
C.............................................................
      END
CDECK  ID>, NORMAL_FORM.
C =================================================================
C SUBROUTINE TO COMPUTE NORMAL FORMS.
C H IS AN ARRAY CONTAINING THE HAMILTONIAN.
C MAXMAP REPRESENTS THE ORDER OF THE TRANSFER MAP.
C MAXORD IS THE ORDER OF THE NORMAL FORM.
C IINV IS A FLAG FOR THE COMPUTATION OF THE INVERSE OF PHI:
C IINV = 0 PSI COMPUTATION SUPPRESSED.
C IINV = 1 PSI COMPUTATION ENABLED.
C IOUT IS A FLAG TO SELECT THE OUTPUT (SEE COMMENTS FOR ROUTINE
C GENERALE_4D).
C
C AUTHOR: A. BAZZANI (M. GIOVANNOZZI - CERN - HAS INTRODUCED SOME
C CHANGES).
C
 
      SUBROUTINE NORMAL_FORM(H,MAXMAP,MAXORD,IINV,IOUT)
C.............................................................
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION  PHI(NDIM,4), U(NDIM,4), FOLD(NDIM,4)
      DIMENSION H(NDIMPLUS)
      DOUBLE PRECISION RPHI(NDIM,4),RPSI(NDIM,4)
      INTEGER NPHI,NU
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/AUTOVALORI/AUTO(4)
      COMMON/RESTI1/RESTO1(NDIM,4)
      COMMON/RESTI2/RESTO2(NDIM,4)
      COMMON/ IDENTITY/ K(4,2)
      COMMON/CAMPO/ A(NDIM,4)
      COMMON/RISONANZA/IRES1(2),IRES2(2),AUTORES(4),IRESON,ICASE
      COMMON/DIRETTA/PHI1(NDIM,4)
      COMMON/INVERSA/PSI(NDIM,4)
      COMMON/MAPPA/F(NDIM,4)
      COMMON/SCRATCH1/Q(NDIM,4)
      COMMON/SCRATCH2/D1(NDIM)
      COMMON/SCRATCH3/R(NDIM,4)
      COMMON/SCRATCH4/D(NDIM,4)
C.............................................................
      NFOUR=4
C.............................................................
      NCOEF=IPUNT(MAXORD,0,0,0)
C........................SAVES THE RESONANT BASIS AND F.......
      IR11=IRES1(1)
      IR12=IRES1(2)
      IR21=IRES2(1)
      IR22=IRES2(2)
      DO J=1,4
        DO I=1,NDIM
          FOLD(I,J)=F(I,J)
        END DO
      END DO
C..............................THE ARRAY K(*,*) IS INITIALIZED
      K(1,1)=1
      K(1,2)=0
      K(2,1)=-1
      K(2,2)=0
      K(3,1)=0
      K(3,2)=1
      K(4,1)=0
      K(4,2)=-1
C........................................PHI AND U SET TO ZERO
      CALL MZEROC(A,NFOUR)
      CALL MZEROC(U,NFOUR)
C.............................................................
      DO NCONT=2, 5
          AUTO(NCONT-1) = F(7-NCONT,NCONT-1)
       ENDDO
C..........................................PHI SET TO IDENTITY
       CALL MONEC(PHI,NFOUR)
C....................................................MAIN LOOP
       DO 100 NORD=2, MAXORD
C........................................COMPUTE THE REMAINDER
          NF=MIN(NORD,MAXMAP)
          NPHI=NORD-1
          NU=NORD-1
C..............................COMPUTES PHI AS A LIE TRANSFORM
          NTERM=IPUNT(NORD,0,0,0)
          CALL MONEC(PSI,NFOUR)
C..............................COMPUTES THE LIE TRANSFORM OF A
          CALL LIE_TRASFC(A,PSI,PHI,NPHI,NPHI,NORD)
C.............................................................
          CALL LIE_TRASFC(A, F, RESTO1, NPHI, NF, NORD)
C.............................................................
          CALL COMPLIN(PHI,AUTO,PSI,NORD)
          CALL LIE_TRASFC(U, PSI, RESTO2, NU, NORD, NORD)
          NTERM1= NORD*(NORD+1)*(NORD+2)*(NORD+3)/24
          NTERM2= (NORD+1)*(NORD+2)*(NORD+3)*(NORD+4)/24
          DO ICONT=NTERM1+1, NTERM2
            DO J=1, 4
               RESTO1(ICONT,J)= RESTO1(ICONT,J)-RESTO2(ICONT,J)
            ENDDO
          ENDDO
C.................................SOLVES THE OMOLOGIC EQUATION
          CALL SOL_EQ_OM(PHI,U,RESTO1,NORD)
100    CONTINUE
C............................COMPUTES THE INVERSE OF PHI (PSI)
       DO J=1,4
         DO ICONT=1, NCOEF
           A(ICONT,J) = -A(ICONT,J)
         ENDDO
       ENDDO
C...............STORES PHI IN PHI1 FOR SUCCESSIVE COMPUTATIONS
       DO ICOMP=1,4
         DO IEL=1,NDIM
           PHI1(IEL,ICOMP)=PHI(IEL,ICOMP)
         ENDDO
       ENDDO
C.............................................................
       IF (IINV.EQ.1) CALL INVERS(PHI,PSI,MAXORD)
C...................................................OUTPUT PHI
       IF (IOUT.GT.0) THEN
         IUNIT=14
         IDIME=NFOUR
C.............................................................
         IF (IOUT.EQ.3) THEN
           CALL MPRINTC(PHI,MAXORD,IDIME,IUNIT) !..PRINTS COMPLEX PHI
         ELSE
           CALL MCTOR(PHI,RPHI,MAXORD,IDIME)
           IF (IOUT.EQ.2) CALL MNORTOPHY(RPHI,MAXORD,IDIME)
           CALL MPRINT(RPHI,MAXORD,IDIME,IUNIT)
         ENDIF
C...................................................OUTPUT PSI
         IF (IINV.EQ.1) THEN
           IUNIT=15
           IDIME=NFOUR
C.............................................................
           IF (IOUT.EQ.3) THEN
             CALL MPRINTC(PSI,MAXORD,IDIME,IUNIT) !..PRINTS COMPLEX PSI
           ELSE
             CALL MCTOR(PSI,RPSI,MAXORD,IDIME)
             IF (IOUT.EQ.2) CALL MNORTOPHY(RPSI,MAXORD,IDIME)
             CALL MPRINT(RPSI,MAXORD,IDIME,IUNIT)
           ENDIF
C.............................................................
         ENDIF
       ENDIF
C..........................................QUASI-RESONANT PART
       IF (IRESON.EQ.2) THEN
C..............COMPUTES THE NORMAL FORM U AND STORES IT IN PHI
          CALL MONEC(PSI,NFOUR)
C..............................COMPUTES THE LIE TRANSFORM OF U
         CALL LIE_TRASFC(U,PSI,PHI,MAXORD,MAXORD,MAXORD)
C..................................................REDEFINES F
          DO J=1,4
             RATIO=AUTO(J)/AUTORES(J)
             DO ICONT=1, NCOEF
                F(ICONT,J)=PHI(ICONT,J)*RATIO
             ENDDO
           ENDDO
           IRES1(1)=0
           IRES2(1)=0
           IRES1(2)=0
           IRES2(2)=0
CC.......................................INITIALIZES PHI AND U
          CALL MONEC(PHI,NFOUR)
          CALL MZEROC(A,NFOUR)
          CALL MZEROC(U,NFOUR)
          DO NCONT=2, 5
             AUTO(NCONT-1) = F(7-NCONT,NCONT-1)
          ENDDO
C.............................................................
          DO 200 NORD=2, MAXORD
            NF=NORD
            NPHI=NORD-1
            NU=NORD-1
C..............................COMPUTES PHI AS A LIE TRANSFORM
            CALL MONEC(PSI,NFOUR)
C..............................COMPUTES THE LIE TRANSFORM OF A
           CALL LIE_TRASFC(A,PSI,PHI,NPHI,NPHI,NORD)
C.............................................................
            CALL LIE_TRASFC(A, F, RESTO1, NPHI, NF, NORD)
            CALL COMPLIN(PHI,AUTO,PSI,NORD)
            CALL LIE_TRASFC(U, PSI, RESTO2, NU, NORD, NORD)
            NTERM1= NORD*(NORD+1)*(NORD+2)*(NORD+3)/24
            NTERM2= (NORD+1)*(NORD+2)*(NORD+3)*(NORD+4)/24
            DO ICONT=NTERM1+1, NTERM2
              DO J=1, 4
                 RESTO1(ICONT,J)= RESTO1(ICONT,J)-RESTO2(ICONT,J)
              ENDDO
            ENDDO
C................................SOLVES THE OMOLOGIC EQUATION
            CALL SOL_EQ_OM(PHI,U,RESTO1,NORD)
200       CONTINUE
C..................................INVERTS PHI (COMPUTES PSI)
          DO J=1,4
             DO ICONT=1, NCOEF
                A(ICONT,J)=-A(ICONT,J)
             ENDDO
          ENDDO
C...................................ADDS THE LINEAR PART TO U
          U(2,4)=LOG(AUTO(4))
          U(3,3)=LOG(AUTO(3))
          U(4,2)=LOG(AUTO(2))
          U(5,1)=LOG(AUTO(1))
C...................................BACK TO RESONANT VARIABLES
          CALL LIE_DERIVAC(U, PHI, PSI, MAXORD, MAXORD, MAXORD)
          CALL LIE_TRASFC(A, PSI, U, MAXORD, MAXORD, MAXORD)
       ENDIF
C.......................COMPUTES THE INTERPOLATING HAMILTONIAN
      CALL HAMIL(U,H,MAXORD)
C.......................RESTORES THE RESONANT BASIS AND F.....
      IRES1(1)=IR11
      IRES1(2)=IR12
      IRES2(1)=IR21
      IRES2(2)=IR22
      DO J=1,4
        DO I=1,NDIM
          F(I,J)=FOLD(I,J)
        END DO
      END DO
C.............................................................
      RETURN
      END
CDECK  ID>, COMPLIN.
C=================================================================
C SUBROUTINE TO COMPOSE THE LINEAR PARTS OF THE MAP.
C P IS AN ARRAY CONTAINING THE INPUT.
C VEC CONTAINS THE LINEAR PARTS.
C Q IS AN ARRAY CONTAINING THE OUTPUT (P*VEC).
C MAXG IS THE ORDER OF THE INPUT ARRAY.
C
C AUTHOR: A. BAZZANI (M. GIOVANNOZZI - CERN - HAS INTRODUCED SOME
C CHANGES).
C
 
      SUBROUTINE COMPLIN(P, VEC,Q, MAXG)
C--------------------------------------------------------
C---- PRECISIONE IMPLICITA
C--------------------------------------------------------
      IMPLICIT COMPLEX*16 (A-H,O-Z)
C--------------------------------------------------------
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
C------------------------------------------------
      DIMENSION P(NDIM,4),Q(NDIM,4),VEC(4)
C------------------------------------------------
      RAP14=VEC(1)/VEC(4)
      RAP24=VEC(2)/VEC(4)
      RAP34=VEC(3)/VEC(4)
      ICONT=1
      EVAL=1.D0
      DO NORD=1, MAXG
         EVAL=VEC(4)*EVAL
         EVAL1=EVAL
         DO N1=0, NORD
            EVAL2=EVAL1
            DO N2=0, NORD-N1
               EVAL3=EVAL2
               DO N3=0, NORD-N1-N2
                  N4=NORD-N1-N2-N3
                  ICONT=ICONT+1
                  DO J=1, 4
                     Q(ICONT,J)=P(ICONT,J)*EVAL3
                  ENDDO
                  EVAL3=EVAL3*RAP34
               ENDDO
               EVAL2=EVAL2*RAP24
            ENDDO
            EVAL1=EVAL1*RAP14
         ENDDO
      ENDDO
      RETURN
      END
CDECK  ID>, COMPOSEC.
C=================================================================
C SUBROUTINE TO COMPOSE COMPLEX POLYNOMIALS.
C THE MEANING OF THE PARAMETERS IS THE SAME AS IN THE ROUTINE
C COMPOSE (SEE ABOVE FOR MORE DETAILS).
C
C AUTHOR: A. BAZZANI (M. GIOVANNOZZI - CERN - HAS INTRODUCED SOME
C CHANGES).
C
 
      SUBROUTINE COMPOSEC(Q, B, C, NQ, NB, MAXG)
C--------------------------------------------------------
C---- PRECISIONE IMPLICITA
C--------------------------------------------------------
      IMPLICIT COMPLEX*16 (A-H,O-Z)
C--------------------------------------------------------
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION Q(NDIM,4)
      DIMENSION B(NDIM,4)
      DIMENSION C(NDIM,4)
*
*     NQ ED NB CONTENGONO RISPETTIVAMENTE I GRADI DEI POLINOMI Q E B
*
      INTEGER NQ,NB
C---------------------------------------------------------
      CALL COMPONIC(Q(1,1),Q(1,2),Q(1,3),Q(1,4),B(1,1),B(1,2),
     .             B(1,3),B(1,4),C(1,1),C(1,2),C(1,3),C(1,4),
     .             NQ,NB,MAXG)
      END
CDECK  ID>, COMPONIC.
C=================================================================
C SUBROUTINE TO PERFORM THE COMPOSITION OF TWO COMPLEX POLYNOMIALS
C STORED IN THE ARRAYS Q*, B*, C*.
C
C AUTHOR: A. BAZZANI (M. GIOVANNOZZI - CERN - HAS INTRODUCED SOME
C CHANGES).
C
 
      SUBROUTINE COMPONIC(Q1, Q2, Q3, Q4, B1, B2, B3, B4,
     .    C1, C2, C3, C4, NQ, NB, MAXG)
C---------------------------------------------------------------
C----- DA TOGLIERE SE SI USA IL CRAY
C----- DA CAMBIARE PER L'USO DEI COMPLESSI
C-----------------------------------------------------------------
      IMPLICIT COMPLEX*16 (A-H,O-Z)
C----------------------------------------------------------------
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION Q1(NDIM),Q2(NDIM),Q3(NDIM),Q4(NDIM)
      DIMENSION B1(NDIM),B2(NDIM),B3(NDIM),B4(NDIM)
      DIMENSION C1(NDIM),C2(NDIM),C3(NDIM),C4(NDIM)
*
*     NQ ED NB CONTENGONO RISPETTIVAMENTE I GRADI DEI POLINOMI Q E B
*     SI CONSIDERANO I B TUTTI DELLO STESSO GRADO PER NON INTERFERIRE
*     CON LA VETTORIZZAZIONE
*
      INTEGER NQ,NB
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/SCRATCH1/A1(NDIM),A2(NDIM),A3(NDIM),A4(NDIM)
      COMMON/SCRATCH2/D1(NDIM)
*
      NCOPPIE(0)=0
      MAXQ=NQ
      MAXP=NB
*
*     CALCOLA GLI INDIRIZZI PER IL PRODOTTO P*Q DI GRADO MAXP E MAXQ
*     TENENDO I TERMINI FINO AL GRADO MAX
*
      CALL INIZPRODO4(MAXG,MAXP,MAXG)
*
*      M = MAXG
*
*     DA' IL NUMERO TOTALE DI TERMINI IN UN POLINOMIO DI GRADO MAXG
*
      NMAX=(MAXG+4)*(MAXG+3)*(MAXG+2)*(MAXG+1)/24
*
*
* AZZERAMENTO INIZIALE COMPLETO DELL'ARRAY RISULTATO:
*
      DO  1  I1=1, NMAX
      A1 (I1)=0.D0
      A2 (I1)=0.D0
      A3 (I1)=0.D0
      A4 (I1)=0.D0
      C1 (I1)=0.D0
      C2 (I1)=0.D0
      C3 (I1)=0.D0
      C4 (I1)=0.D0
  1   D1 (I1)=0.D0
*
* L'ARRAY "A1" SCRATCH E` INIZIALIZZATO ALLA COSTANTE 1.
      D1(1) = 1.D0
      A1(1) = 1.D0
      A2(1) = 1.D0
      A3(1) = 1.D0
      A4(1) = 1.D0
*
* INGRESSO NEL QUADRUPLO LOOP CHE NUMERA I CONTRIBUTI ALLA COMPOSIZIONE:
*
      DO  2  J1=0, MAXQ
*
* INGRESSO NEL SECONDO PASSO DEL QUADRUPLO LOOP:
*
      DO  3  J2=0, MAXQ-J1
*
* INGRESSO NEL TERZO PASSO DEL QUADRUPLO LOOP:
      DO  4  J3=0, MAXQ-J1-J2
*
* INGRESSO NEL QUARTO PASSO DEL QUADRUPLO LOOP:
      DO  5  J4=0, MAXQ-J1-J2-J3
*
* NELL'ARRAY A4 E' CONTENUTO IL CONTRIBUTO DEL MONOMIO DI INDICE
* "J1,J2,J3,J4".
* ORA TUTTI GLI ELEMENTI DELL'ARRAY "A4" SONO MOLTIPLICATI PER IL COEFFICIENTE
* DI INDICI "J1,J2,J3,J4" DELL'ARRAY COMPONENDO "Q1,Q2,Q3,Q4"
* E ACCUMULATI NELL'ARRAY RISULTATO
*
      VAR=Q1(IPUNT(J1,J2,J3,J4))
      IF (CDABS(VAR).GE.1.E-14) THEN
         DO   I1=1, NMAX
              C1(I1)=C1(I1)+
     &               A4(I1)*VAR
         ENDDO
      ENDIF
      VAR=Q2(IPUNT(J1,J2,J3,J4))
      IF (CDABS(VAR).GE.1.E-14) THEN
         DO   I1=1, NMAX
              C2(I1)=C2(I1)+
     &               A4(I1)*VAR
         ENDDO
      ENDIF
      VAR=Q3(IPUNT(J1,J2,J3,J4))
      IF (CDABS(VAR).GE.1.E-14) THEN
         DO   I1=1, NMAX
              C3(I1)=C3(I1)+
     &               A4(I1)*VAR
         ENDDO
      ENDIF
      VAR=Q4(IPUNT(J1,J2,J3,J4))
      IF (CDABS(VAR).GE.1.E-14) THEN
         DO  I1=1, NMAX
              C4(I1)=C4(I1)+
     &               A4(I1)*VAR
         ENDDO
      ENDIF
*
* L'ARRAY SCRATCH "A4",PROVENIENTE DALL'ESTERNO DEL LOOP,E` MOLTIPLICATO
* PER L'ARRAY "B4", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
*
      IF (J4.LT.MAXQ-J3-J2-J1) THEN
           CALL PRODOC(A4, B4, D1, MAXG)
*
* DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "D1", E`
* RIMESSO IN "A4" CHE RISULTA COSI` AGGIORNATO AL VALORE
* B1**J1*B2**J2*B3**J3*B4**(J4+1)
*
           DO  59  I1=1, NMAX
 59        A4(I1) = D1(I1)
      ENDIF
 5    CONTINUE
*
* L'ARRAY SCRATCH "A3", PROVENIENTE DALL'ESTERNO DEL LOOP, E` MOLTIPLICATO
* PER L'ARRAY "B3", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
*
      IF (J3.LT.MAXQ-J1-J2) THEN
           CALL PRODOC(A3, B3, A4, MAXG)
*
* DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "A4", E`
* RIMESSO IN "A3" CHE RISULTA COSI` AGGIORNATO AL VALORE
* B1**J1*B2**J2*B3**(J3+1)
*
          DO  29  I1=1, NMAX
 29       A3(I1) = A4(I1)
      ENDIF
*
 4    CONTINUE
*
* L'ARRAY SCRATCH "A2" ,PROVENIENTE DALL'ESTERNO DEL LOOP,E` MOLTIPLICATO
* PER L'ARRAY "B2", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
      IF (J2.LT.MAXQ-J1) THEN
           CALL PRODOC(A2, B2, A4, MAXG)
*
* DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "D1", E`
* RIMESSO IN "A2" CHE RISULTA COSI` AGGIORNATO AL VALORE B1**J1*B2**(J2+1)
*
* L'ARRAY "A3" SCRATCH E` INIZIALIZZATO AL VALORE CORRENTE DELL'ARRAY "A2":
*
           DO  26  I1=1, NMAX
           A3(I1) = A4(I1)
 26        A2(I1) = A3(I1)
      ENDIF
*
 3    CONTINUE
*
* L'ARRAY SCRATCH "A1" ,PROVENIENTE DALL'ESTERNO DEL LOOP,E` MOLTIPLICATO
* PER L'ARRAY "B1", IN MODO DA ESTRARRE TUTTI I CONTRIBUTI DI GRADO "M":
      IF (J1.LT.MAXQ) THEN
           CALL PRODOC(A1, B1, A4, MAXG)
*
* DOPO DI CHE IL RISULTATO DEI PRODOTTI, POSTO NELL'ARRAY SCRATCH "D1", E`
* RIMESSO IN "A1" CHE RISULTA COSI` AGGIORNATO AL VALORE B1**J1
*
* L'ARRAY "A2" SCRATCH E` INIZIALIZZATO AL VALORE CORRENTE DELL'ARRAY "A1":
*
           DO  23  I1=1, NMAX
           A3(I1) = A4(I1)
           A2(I1) = A3(I1)
 23        A1(I1) = A2(I1)
      ENDIF
*
  2   CONTINUE
*
*     SCRITTURE
*
*      ICONT=0
*      DO 13 NORD=0, MAXG
*       NORD=MAXG
*       ICONT=(MAXG+1)*MAXG/2
*        DO 13 N1=0, NORD
*        DO 13 N2=0,NORD-N1
*        DO 13 N3=0,NORD-N1-N2
*        N4=NORD-N1-N2-N3
*          ICONT=ICONT+1
      END
CDECK  ID>, PRODOC.
C=================================================================
C SUBROUTINE TO COMPUTE THE COMPOSITION OF THE POLYNOMIALS.
C
C AUTHOR: A. BAZZANI (M. GIOVANNOZZI - CERN - HAS INTRODUCED SOME
C CHANGES).
C
 
      SUBROUTINE PRODOC(Q, P, R, MAXG)
      IMPLICIT COMPLEX*16 (A-H,O-Z)
C------------------------------------------------------------
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION P(NDIM)
      DIMENSION Q(NDIM)
      DIMENSION R(NDIM)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      MAXNT=IPUNT(MAXG,0,0,0)
      DO 1 J=1, MAXNT
1      R(J)=0.D0
      DO 2 J=1, MAXNT
      DO 2 NJ=NCOPPIE(J-1)+1, NCOPPIE(J)
2           R(J) = R(J) + P(INDP(NJ))*Q(INDQ(NJ))
      RETURN
      END
CDECK  ID>, HAMIL.
C=================================================================
C SUBROUTINE TO COMPUTE A HAMILTONIAN STARTING FROM THE VECTOR
C FIELD.
C H IS AN ARRAY CONTAINING THE HAMILTONIAN.
C A IS THE FIELD.
C MAXG IS THE ORDER FOR THE COMPUTATIONS.
C
C AUTHOR: A. BAZZANI (M. GIOVANNOZZI - CERN - HAS INTRODUCED SOME
C CHANGES).
C
 
      SUBROUTINE HAMIL(A,H,MAXG)
      IMPLICIT COMPLEX*16 (A-H,O-Z)
C------------------------------------------------------------
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION A(NDIM,4),H(NDIMPLUS)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/RESTI1/C(NDIM,4)
      COMMON/RESTI2/D(NDIM,4)
      COMMON/AUTOVALORI/AUTO(4)
C----------------------------------------------------------
      MAXNT=IPUNT(MAXG,0,0,0)
C-----------------------------------------------------------
C---- INIZIALIZZO LA PARTE LINEARE
C----------------------------------------------------------
      A(2,4)=LOG(AUTO(4))
      A(3,3)=LOG(AUTO(3))
      A(4,2)=LOG(AUTO(2))
      A(5,1)=LOG(AUTO(1))
      ICONT=0
      DO 100 N=0, MAXG
         DO N1=0, N
            DO N2=0, N-N1
               DO N3=0, N-N1-N2
                  N4=N-N1-N2-N3
                  ICONT=ICONT+1
                  H(IPUNT(N1,N2,N3,N4+1))=A(ICONT,3)/(N4+1)
               ENDDO
            ENDDO
         ENDDO
         DO N1=0, N
            DO N2=0, N-N1
                 N3=N-N1-N2
                 H(IPUNT(N1,N2,N3+1,0))=-A(IPUNT(N1,N2,N3,0),4)/(N3+1)
            ENDDO
         ENDDO
         DO N1=0, N
            N2=N-N1
            H(IPUNT(N1,N2+1,0,0))=A(IPUNT(N1,N2,0,0),1)/(N2+1)
         ENDDO
           DO N1=0, N
              H(IPUNT(N1+1,0,0,0))=-A(IPUNT(N1,0,0,0),2)/(N1+1)
           ENDDO
100        CONTINUE
           RETURN
           END
CDECK  ID>, LIE_TRASFC.
C=================================================================
C SUBROUTINE TO COMPUTE LIE TRANSFORMS OF COMPLEX POLYNOMIALS.
C A IS AN ARRAY CONTAINING THE OPERATOR
C S IS AN ARRAY CONTAINING THE POLYNOMIAL TO BE TRANSFORMED.
C C IS AN ARRAY CONTAINING THE LIE TRANSFORM.
C MAXA REPRESENTS THE ORDER OF A.
C MAXQ REPRESENTS THE ORDER OF S.
C MAXG REPRESENTS THE ORDER OF C
C
C AUTHOR: A. BAZZANI (M. GIOVANNOZZI - CERN - HAS INTRODUCED SOME
C CHANGES).
C
 
      SUBROUTINE LIE_TRASFC(A,S,C,MAXA,MAXQ,MAXG)
      IMPLICIT COMPLEX*16 (A-H,O-Z)
C------------------------------------------------------------
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION A(NDIM,4),C(NDIM,4),S(NDIM,4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/SCRATCH1/Q(NDIM,4)
      COMMON/SCRATCH3/R(NDIM,4)
C------------------------------------------------
      MAXNT=IPUNT(MAXG,0,0,0)
C     LE SEGUENTE DEFINIZIONE SERVE AD EVITARE WARNING IN 
C     COMPILAZIONE.
      PIPPO=MAXQ
C---------------------------------------------------------
C----- INIZIALIZZO C
C---------------------------------------------------------
      DO ITERM=1, MAXNT
        Q(ITERM,1)=S(ITERM,1)
        Q(ITERM,2)=S(ITERM,2)
        Q(ITERM,3)=S(ITERM,3)
        Q(ITERM,4)=S(ITERM,4)
        C(ITERM,1)=Q(ITERM,1)
        C(ITERM,2)=Q(ITERM,2)
        C(ITERM,3)=Q(ITERM,3)
        C(ITERM,4)=Q(ITERM,4)
      ENDDO
C-----------------------------------------------------------
C----- MAIN LOOP
C---------------------------------------------------------
      DO 1 NORD=1, MAXG-1
         CALL LIE_DERIVAC(A,Q,R,MAXA,MAXG,MAXG)
         DO ITERM=1, MAXNT
            Q(ITERM,1)=R(ITERM,1)/NORD
            Q(ITERM,2)=R(ITERM,2)/NORD
            Q(ITERM,3)=R(ITERM,3)/NORD
            Q(ITERM,4)=R(ITERM,4)/NORD
         ENDDO
         DO ITERM=1, MAXNT
            C(ITERM,1)=C(ITERM,1)+Q(ITERM,1)
            C(ITERM,2)=C(ITERM,2)+Q(ITERM,2)
            C(ITERM,3)=C(ITERM,3)+Q(ITERM,3)
            C(ITERM,4)=C(ITERM,4)+Q(ITERM,4)
          ENDDO
1     CONTINUE
      RETURN
      END
CDECK  ID>, LIE_DERIVAC.
C=================================================================
C SUBROUTINE TO COMPUTE LIE DERIVATIVES OF COMPLEX POLYNOMIALS.
C A IS AN ARRAY CONTAINING THE OPERATOR
C C IS AN ARRAY CONTAINING THE POLYNOMIAL TO BE TRANSFORMED.
C R IS AN ARRAY CONTAINING THE LIE DERIVATIVE.
C NA REPRESENTS THE ORDER OF A.
C NC REPRESENTS THE ORDER OF C.
C MAXG REPRESENTS THE ORDER OF R
C
C AUTHOR: A. BAZZANI (M. GIOVANNOZZI - CERN - HAS INTRODUCED SOME
C CHANGES).
C
 
      SUBROUTINE LIE_DERIVAC(A,C,R,NA,NC,MAXG)
      IMPLICIT COMPLEX*16 (A-H,O-Z)
C------------------------------------------------------------
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION A(NDIM,4),C(NDIM,4),R(NDIM,4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/SCRATCH4/D(NDIM,4)
      COMMON/SCRATCH2/B(NDIM)
C------------------------------------------------
      MAXNT=IPUNT(MAXG,0,0,0)
      CALL INIZPRODO4(NA,NC-1,MAXG)
      DO J=1,4
        ICONT=0
        DO 1 NORD=0, MAXG
           DO 1 N1=0, NORD
              DO 1 N2=0, NORD-N1
                 DO 1 N3=0, NORD-N1-N2
                    N4=NORD-N1-N2-N3
                    ICONT=ICONT+1
                    IF (N1.GE.1) D(IDERIV(ICONT,1),1)=N1*C(ICONT,J)
                    IF (N2.GE.1) D(IDERIV(ICONT,2),2)=N2*C(ICONT,J)
                    IF (N3.GE.1) D(IDERIV(ICONT,3),3)=N3*C(ICONT,J)
1                   IF (N4.GE.1) D(IDERIV(ICONT,4),4)=N4*C(ICONT,J)
         CALL PRODOC(A(1,1), D(1,1) , B, MAXG)
         DO ICONT=1, MAXNT
            R(ICONT,J)=B(ICONT)
         ENDDO
         CALL PRODOC(A(1,2), D(1,2) , B, MAXG)
         DO ICONT=1, MAXNT
            R(ICONT,J)=R(ICONT,J)+B(ICONT)
         ENDDO
         CALL PRODOC(A(1,3), D(1,3) , B, MAXG)
         DO ICONT=1, MAXNT
            R(ICONT,J)=R(ICONT,J)+B(ICONT)
         ENDDO
         CALL PRODOC(A(1,4), D(1,4) , B, MAXG)
         DO ICONT=1, MAXNT
            R(ICONT,J)=R(ICONT,J)+B(ICONT)
         ENDDO
       ENDDO
       RETURN
       END
CDECK  ID>, SOL_EQ_OM.
C=================================================================
C SUBROUTINE TO SOLVE THE OMOLOGIC EQUATION AND TO PERFORM THE
C SIMPLECTIFICATION OF THE SOLUTION.
C PHI IS A COMPLEX ARRAY CONTANING THE NORMAL FORM TRANSFORMATION.
C U IS A COMPLEX ARRAY CONTANING THE NORMAL FORM.
C RESTO IS AN ARRAY CONTAING THE REMAINDER OF THE OMOLOGIC EQUATION.
C MAXG IS THE ORDER OF THE COMPOSITION.
C THE REST IS TOP SECRET.
C
C AUTHOR: A. BAZZANI
C
 
      SUBROUTINE SOL_EQ_OM(PHI, U, RESTO, MAXG)
C-----------------------------------------------------------------------------
C---- PRECISIONE IMPLICITA
C-----------------------------------------------------------------------------
      IMPLICIT COMPLEX*16 (A-H,O-Z)
C-----------------------------------------------------------------------------
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION PHI(NDIM,4),U(NDIM,4),RESTO(NDIM,4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/AUTOVALORI/AUTO(4)
      COMMON/CAMPO/ A(NDIM,4)
      COMMON/ IDENTITY/ K(4,2)
      COMMON/RISONANZA/IRES1(2),IRES2(2),AUTORES(4),IRESON,ICASE
      COMMON/INVERSA/C(NDIM,4)
      COMMON/RESTI2/D(NDIM,4)
      COMMON/SCRATCH1/Q(NDIM,4)
C-----------------------------------------------------------------------------
C---- INIZIALIZZO
C-----------------------------------------------------------------------------
      NTERM1=IPUNT(MAXG-1,0,0,0)
      NTERM2=IPUNT(MAXG,0,0,0)
C-----------------------------------------------------------------------------
C---- FACCIO UN LOOP SU TUTTI I TEMINI DI ORDINE N DISCRIMINANDO I TERMINI
C---- CHE SONO COMBINAZIONI INTERE DEI VETTORI RISONANTI
C-----------------------------------------------------------------------------
      INORM1=IRES1(1)*IRES1(1)+IRES2(1)*IRES2(1)
      INORM2=IRES1(2)*IRES1(2)+IRES2(2)*IRES2(2)
      IPRODV=IRES1(1)*IRES2(2)-IRES2(1)*IRES1(2)
      ICONT=NTERM1
      PD1=1.D0
      PD4=AUTO(4)**MAXG
      RAP14=AUTO(1)/AUTO(4)
      RAP24=AUTO(2)/AUTO(4)
      RAP34=AUTO(3)/AUTO(4)
      DO N1=0, MAXG
         PD2=PD1
         DO N2=0, MAXG-N1
            PD3=PD2
            DO N3=0, MAXG-N1-N2
               N4=MAXG-N1-N2-N3
               ICONT=ICONT+1
               PD=PD3*PD4
               DO J=1, 4
                  PDEN=PD-AUTO(J)
                  K1=N1-N2-K(J,1)
                  K2=N3-N4-K(J,2)
                  IMOD1=1
                  IMOD2=1
                  IF (INORM2.NE.0) THEN
                      IPRODV1=K1*IRES2(1)-K2*IRES1(1)
                      IPRODV2=K1*IRES2(2)-K2*IRES1(2)
                      IMOD2=MOD(-IPRODV1,IPRODV)
                      IMOD1=MOD(IPRODV2,IPRODV)
                  ELSE
                      IF (INORM1.NE.0) THEN
                         ISCAL1=K1*IRES1(1)+K2*IRES2(1)
                         ISCAL2=K1*IRES2(1)-K2*IRES1(1)
                         IMOD1=MOD(ISCAL1,INORM1)
                         IMOD2=ISCAL2
                      ENDIF
                  ENDIF
                  IF ((K1.EQ.0.AND.K2.EQ.0).OR.
     &                (IMOD1.EQ.0.AND.IMOD2.EQ.0)) THEN
                      A(ICONT,J)=0.D0
                      U(ICONT,J)=RESTO(ICONT,J)/AUTO(J)
C-----------------------------------------------------------------------------
C-------- LA SCELTA NON SIMPLETTICA E'  PHI(ICONT,J)=0.D0
C-----------------------------------------------------------------------------
                  ELSE
                    IF (CDABS(RESTO(ICONT,J)).GT.1D-18) THEN
                      A(ICONT,J)=RESTO(ICONT,J)/PDEN
                      PHI(ICONT,J)=PHI(ICONT,J)+A(ICONT,J)
                    ELSE
                      A(ICONT,J)=0D0
                    ENDIF
                    U(ICONT,J)=0.D0
                  ENDIF
                ENDDO
               PD3=PD3*RAP34
            ENDDO
            PD2=PD2*RAP24
         ENDDO
         PD1=PD1*RAP14
      ENDDO
      RETURN
      END
CDECK  ID>, INVERS.
C=================================================================
C SUBROUTINE TO INVERT THE TRANSFORMATION PHI.
C PHI IS A COMPLEX ARRAY CONTAINING THE TRANSFORMATION.
C PSI IS A COMPLEX ARRAY CONTAINING THE INVERSE.
C MAXG REPRESENTS THE ORDER OF THE POLYNOMIALS.
C
C AUTHOR: A. BAZZANI (M. GIOVANNOZZI - CERN - HAS INTRODUCED SOME
C CHANGES).
C
 
      SUBROUTINE INVERS(PHI,PSI,MAXG)
C-------------------------------------------------------
C---- PRECISIONE IMPLICITA
C-------------------------------------------------------
      IMPLICIT COMPLEX*16 (A-H,O-Z)
C-------------------------------------------------------
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION PHI(NDIM,4),PSI(NDIM,4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/CAMPO/ A(NDIM,4)
      COMMON/ IDENTITY/ K(4,2)
      COMMON/RESTI1/C(NDIM,4)
      COMMON/SCRATCH3/R(NDIM,4)
C---------------------------------------------------
C---- INIZIALIZZAZIONE
C--------------------------------------------------
      NFOUR=4
C     LA SEGUENTE DEFINIZIONE SERVE PER EVITARE WARNING IN 
C     COMPILAZIONE.
      PIPPO=PHI(1,1)
C--------------------------------------------------
      CALL MONEC(C,NFOUR)
      CALL MONEC(PSI,NFOUR)
C--------------------------------------------------
      NTERM=IPUNT(MAXG,0,0,0)
C-----------------------------------------------------------
C----- MAIN LOOP
C---------------------------------------------------------
      DO 1 NORD=1, MAXG
         CALL LIE_DERIVAC(A,C,R,MAXG,MAXG,MAXG)
         DO ITERM=1, NTERM
            C(ITERM,1)=R(ITERM,1)/NORD
            C(ITERM,2)=R(ITERM,2)/NORD
            C(ITERM,3)=R(ITERM,3)/NORD
            C(ITERM,4)=R(ITERM,4)/NORD
         ENDDO
         DO ITERM=6, NTERM
            PSI(ITERM,1)=PSI(ITERM,1)+C(ITERM,1)
            PSI(ITERM,2)=PSI(ITERM,2)+C(ITERM,2)
            PSI(ITERM,3)=PSI(ITERM,3)+C(ITERM,3)
            PSI(ITERM,4)=PSI(ITERM,4)+C(ITERM,4)
          ENDDO
1     CONTINUE
      RETURN
      END
CDECK  ID>, LIE_TRASF.
C=================================================================
C SUBROUTINE TO COMPUTE LIE TRANSFORMS OF COMPLEX POLYNOMIALS.
C A IS AN ARRAY CONTAINING THE OPERATOR
C S IS AN ARRAY CONTAINING THE POLYNOMIAL TO BE TRANSFORMED.
C C IS AN ARRAY CONTAINING THE LIE TRANSFORM.
C MAXA REPRESENTS THE ORDER OF A.
C MAXQ REPRESENTS THE ORDER OF S.
C MAXD REPRESENTS THE ORDER OF C
C
C AUTHOR: A. BAZZANI (M. GIOVANNOZZI - CERN - HAS INTRODUCED SOME
C CHANGES).
C
 
      SUBROUTINE LIE_TRASF(A,S,C,MAXA,MAXQ,MAXD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.............................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION A(NDIM,4),C(NDIM,4),S(NDIM,4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/SCRATCH11/Q(NDIM,4)
      COMMON/SCRATCH33/R(NDIM,4)
C.............................................................
      MAXNT=IPUNT(MAXD,0,0,0)
C     LA SEGUENTE DEFINIZIONE SERVE PER EVITARE WARNING IN 
C     COMPILAZIONE
      PIPPO=MAXQ
C.............................................................
C INIZIALIZZO C
C.............................................................
      DO ITERM=1, MAXNT
        Q(ITERM,1)=S(ITERM,1)
        Q(ITERM,2)=S(ITERM,2)
        Q(ITERM,3)=S(ITERM,3)
        Q(ITERM,4)=S(ITERM,4)
        C(ITERM,1)=Q(ITERM,1)
        C(ITERM,2)=Q(ITERM,2)
        C(ITERM,3)=Q(ITERM,3)
        C(ITERM,4)=Q(ITERM,4)
      ENDDO
C.............................................................
C MAIN LOOP
C.............................................................
      CALL INIZPRODO4(MAXA,MAXD-1,MAXD)
      DO 1 NORD=1, MAXD-1
         CALL LIE_DERIVA(A,Q,R,MAXA,MAXD,MAXD)
         DO ITERM=1, MAXNT
            Q(ITERM,1)=R(ITERM,1)/NORD
            Q(ITERM,2)=R(ITERM,2)/NORD
            Q(ITERM,3)=R(ITERM,3)/NORD
            Q(ITERM,4)=R(ITERM,4)/NORD
         ENDDO
         DO ITERM=1, MAXNT
            C(ITERM,1)=C(ITERM,1)+Q(ITERM,1)
            C(ITERM,2)=C(ITERM,2)+Q(ITERM,2)
            C(ITERM,3)=C(ITERM,3)+Q(ITERM,3)
            C(ITERM,4)=C(ITERM,4)+Q(ITERM,4)
          ENDDO
1     CONTINUE
      RETURN
      END
CDECK  ID>, LIE_DERIVA.
C=================================================================
C SUBROUTINE TO COMPUTE LIE DERIVATIVES OF REAL POLYNOMIALS.
C A IS AN ARRAY CONTAINING THE OPERATOR
C C IS AN ARRAY CONTAINING THE POLYNOMIAL TO BE TRANSFORMED.
C R IS AN ARRAY CONTAINING THE LIE DERIVATIVE.
C NA REPRESENTS THE ORDER OF A.
C NC REPRESENTS THE ORDER OF C.
C MAXD REPRESENTS THE ORDER OF R
C IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!
C THE FIRST AND THE THIRD COMPONENT ARE SKIPPED IN ORDER TO
C IMPROVE THE PERFORMANCE.
C
C AUTHOR: A. BAZZANI (M. GIOVANNOZZI - CERN - HAS INTRODUCED SOME
C CHANGES).
C
 
      SUBROUTINE LIE_DERIVA(A,C,R,NA,NC,MAXD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.............................................................
      PARAMETER (NORMAX=10,NORMAX1=NORMAX+1,
     . NDIM=(NORMAX+4)*(NORMAX+3)*(NORMAX+2)*(NORMAX+1)/24,
     . NDIMPLUS=(NORMAX+5)*(NORMAX+4)*(NORMAX+3)*(NORMAX+2)/24,
     . NCOP=(NORMAX+8)*(NORMAX+7)*(NORMAX+6)*(NORMAX+5)*NDIM/1680)
      DIMENSION A(NDIM,4),C(NDIM,4),R(NDIM,4)
      COMMON/INDIRIZZI/IPUNT(0:NORMAX1,0:NORMAX1,0:NORMAX1,0:NORMAX1),
     .   INDP(NCOP),INDQ(NCOP),NCOPPIE(0:NDIM),IDERIV(NDIMPLUS,4)
      COMMON/SCRATCH44/D(NDIM,4)
      COMMON/SCRATCH22/B(NDIM)
C.............................................................
      MAXNT=IPUNT(MAXD,0,0,0)
C     LA SEGUENTE DEFINIZIONE SERVE PER EVITARE WARNING IN 
C     COMPILAZIONE
      PIPPO=NA
      PLUTO=NC
C      CALL INIZPRODO4(NA,NC-1,MAXD)
      DO J=1,4
        ICONT=0
        DO 1 NORD=0, MAXD
           DO 1 N1=0, NORD
              DO 1 N2=0, NORD-N1
                 DO 1 N3=0, NORD-N1-N2
                    N4=NORD-N1-N2-N3
                    ICONT=ICONT+1
                    IF (N2.GE.1) D(IDERIV(ICONT,2),2)=N2*C(ICONT,J)
1                   IF (N4.GE.1) D(IDERIV(ICONT,4),4)=N4*C(ICONT,J)
         CALL PRODO(A(1,2), D(1,2) , B, MAXD)
         DO ICONT=1, MAXNT
            R(ICONT,J)=B(ICONT)
         ENDDO
         CALL PRODO(A(1,4), D(1,4) , B, MAXD)
         DO ICONT=1, MAXNT
            R(ICONT,J)=R(ICONT,J)+B(ICONT)
         ENDDO
       ENDDO
       RETURN
       END
