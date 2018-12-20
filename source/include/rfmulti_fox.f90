
  !---- Zero the arrays



 

irrtr = irm_rf(ix)
nordm=nmu_rf(irrtr)
crabfreq = freq_rfm(irrtr)
  

!FOX  KCRABDA=(SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI) ;
    
!FOX  DXRF=XL*C1M3;
!FOX  DYRF=ZL*C1M3;
print *, "crabdaaaaa", dare(kcrabda)
print *, "DXRFaaaa", dare(dxrf)
print *, "DYRFaaaa", dare(dyrf)

 do iord = 1, nordm
!FOX  FSIDA(1, IORD)=(NORRFAMP(IORD, IRRTR)* 
!FOX  COS(NORRFPH(IORD, IRRTR)*2D0*PI  - KCRABDA)) ;

!FOX  FSIDA(1, IORD)=(NORRFAMP(IORD, IRRTR)*
!FOX  SIN(NORRFPH(IORD, IRRTR)*2D0*PI  - KCRABDA)) ;

!FOX  FCODA(2, IORD)=(SKRFAMP(IORD, IRRTR)*
!FOX  COS(SKRFPH(IORD, IRRTR)*2D0*PI  - KCRABDA)) ;

!FOX  FSIDA(2, IORD)=(SKRFAMP(IORD, IRRTR)*
!FOX  SIN(SKRFPH(IORD, IRRTR)*2D0*PI  - KCRABDA)) ;
print *, "sida1daaaaa", dare(fsida(1,iord))
print *, "sida1daaaaa", dare(fsida(2,iord))
print *, "sida1daaaaa", dare(fcoda(1,iord))
print *, "sida1daaaaa", dare(fcoda(2,iord))
enddo


print *, "sid2aaaaa", dare(fsida(1,2))
print *, "cod1aaaaa", dare(fcoda(1,1))

do iord = nordm, 1, -1

!FOX  CP_RE=(CP_RE*DXRF - CP_IM*DYRF)/(IORD)  + FCODA(1, IORD) ;
!FOX  CP_IM=(CP_RE*DYRF + CP_IM*DXRF)/(IORD)  + FCODA(2, IORD) ;

!FOX  SP_RE=(SP_RE*DXRF - SP_IM*DYRF)/(IORD+1)  + FSIDA(1, IORD) ;
!FOX  SP_IM=(SP_RE*DYRF + SP_IM*DXRF)/(IORD+1)  + FSIDA(2, IORD) ;
print *, "sida1daaaaapppppsp_re", dare(sp_re)
print *, "sida1daaaaapppppsp_re", dare(sp_im)

print *, "sida1daaaaapppppcp_re", dare(cp_re)
print *, "sida1daaaaapppppcp_re", dare(cp_im)
enddo

!FOX  SP_RE = SP_RE*DXRF;
!FOX  SP_IM = SP_IM*DYRF;

print *, "sida1daaaaapppppsp_re_after", dare(sp_re)
print *, "sida1daaaaapppppsp_re_after", dare(sp_im)

print *, "sida1daaaaapppppsp_re_aftercp", dare(cp_re)
print *, "sida1daaaaapppppsp_re_aftercp", dare(cp_im)
!FOX  Y(1)=Y(1) - CP_RE*C1E3*MTCDA/(ONE+DPDA);
!FOX  Y(2)=Y(2) + CP_IM*C1E3*MTCDA/(ONE+DPDA);
print *, "sida1daaaaapppppsp_y1", dare(y(1))
print *, "sida1daaaaapppppsp_y2", dare(y(2))
!FOX  EJ1 = EJ1 - C1E3*E0F*SP_RE/((CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI);
print *, "sida1daaaaapppppsp_y1", dare(y(1))
print *, "sida1daaaaapppppsp_y2", dare(y(2))



   



