
  !---- Zero the arrays

normal = zero
skew = zero
pnl = zero
psl = zero

 

    irrtr = irm_rf(ix)
  !  NORMALDA = nor_rf_amp(irrtr,:)
  !  SKEWDA = skew_rf_amp(irrtr,:)
  !  PNLDA = nor_rf_ph(irrtr,:)
  !  PSLDA = skew_rf_ph(irrtr,:)
    nordm=nmu_rf(irrtr)
    crabfreq = freq_rfm(irrtr)
  

!FOX  KCRABDA=(SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI) ;
    
! FOX  DXI=XL*C1M3;
! FOX  DYI=ZL*C1M3;

 do iord = 0, nordm
! FOX  FIELDCDA(1, IORD) = (NORMAL(IORD) * COS(PNL(IORD)*2D0*PI  - KCRABDA))
! FOX  FIELDSDA(1, IORD) = (NORMAL(IORD) * SIN(PNL(IORD)*2D0*PI  - KCRABDA))
! FOX  FIELDCDA(2, IORD) = (SKEW(IORD)   * COS(PSL(IORD)*2D0*PI  - KCRABDA))
! FOX  FIELDSDA(2, IORD) = (SKEW(IORD)   * SIN(PSL(IORD)*2D0*PI  - KCRABDA))
enddo


!    do iord = nordm, 1, -1

! FOX CP_RE =CP_RE*DXI - CP_IM*DYI  + FIELDC(1, IORD);
! FOX CP_IM =CP_RE*DYI + CP_IM*DXI  + FIELDC(2, IORD);

! FOX SP_RE =SP_RE*DXI - SP_IM*DYI  + FIELDS(1, IORD);
! FOX SP_IM =SP_RE*DYI + SP_IM*DXI  + FIELDS(2, IORD);

!    enddo
!  FOX SP_RE = SP_RE*DXI;
!  FOX SP_IM = SP_IM*DYI;


!  FOX Y(1)=Y(1) - CP_RE*C1E3*MTCDA/(ONE+DPDA);
!  FOX Y(2)=Y(2) + CP_IM*C1E3*MTCDA/(ONE+DPDA);
!  FOX EJ1 = EJ1 - C1E3*E0F*SP_RE*KCRABDA/SIGMDA;

   



