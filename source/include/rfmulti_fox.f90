
  !---- Zero the arrays



 

irrtr = irm_rf(ix)
nordm=nmu_rf(irrtr)
crabfreq = freq_rfm(irrtr)
  

!FOX  KCRABDA=(SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI) ;
    
!FOX  DXRF=ONE*C1M3 ;
!FOX  DYRF=3D0*C1M3 ;
print *, "crabdaaaaa", dare(kcrabda)
print *, "DXRFaaaa", dare(dxrf)
print *, "DYRFaaaa", dare(dyrf)

 do iord = 1, nordm
!FOX  FCODA(1, IORD)=(NORRFAMP(IRRTR, IORD)* 
!FOX  COS(NORRFPH(IRRTR, IORD)*2D0*PI  - KCRABDA)) ;

!FOX  FSIDA(1, IORD)=(NORRFAMP(IRRTR, IORD)*
!FOX  SIN(NORRFPH(IRRTR, IORD)*2D0*PI  - KCRABDA)) ;

!FOX  FCODA(2, IORD)=(SKRFAMP(IRRTR, IORD)*
!FOX  COS(SKRFPH(IRRTR, IORD)*2D0*PI  - KCRABDA)) ;

!FOX  FSIDA(2, IORD)=(SKRFAMP(IRRTR, IORD)*
!FOX  SIN(SKRFPH(IRRTR, IORD)*2D0*PI  - KCRABDA)) ;
print * , "irrrtr", irrtr
print *, "sida1daaaaa", iord, dare(fsida(1,iord))
print *, "sida1daaaaa",iord, dare(fsida(2,iord))
print *, "sida1daaaaa",iord, dare(fcoda(1,iord))
print *, "sida1daaaaa",iord, dare(fcoda(2,iord))

enddo

!FOX  CP_RE=ZERO ;
!FOX  CP_IM=ZERO ;
!FOX  SP_RE=ZERO ;
!FOX  SP_IM=ZERO ;
do iord = nordm, 1, -1
! Måste fixaaaa temp variablens!!!
!FOX  CP_RETP = CP_RE ;
!FOX  CP_RE=(CP_RETP*DXRF - CP_IM*DYRF)/(IORD)  + FCODA(1, IORD) ;
!FOX  CP_IM=(CP_RETP*DYRF + CP_IM*DXRF)/(IORD)  + FCODA(2, IORD) ;

!FOX  SP_RETP = SP_RE ;
!FOX  SP_RE=(SP_RETP*DXRF - SP_IM*DYRF)/(IORD+1)  + FSIDA(1, IORD) ;
!FOX  SP_IM=(SP_RETP*DYRF + SP_IM*DXRF)/(IORD+1)  + FSIDA(2, IORD) ;
print *, "sida1daaaaapppppsp_resp", dare(sp_re)
print *, "sida1daaaaapppppsp_resp", dare(sp_im)

print *, "sida1daaaaapppppcp_recp", dare(cp_re)
print *, "sida1daaaaapppppcp_recp", dare(cp_im)
enddo

!FOX  SP_RETP = SP_RE ;
!FOX  SP_RE=(SP_RETP*DXRF - SP_IM*DYRF) ;
!FOX  SP_IM=(SP_RETP*DYRF + SP_IM*DXRF) ;

print *, "sida1daaaaapppppsp_re_after", dare(sp_re)
print *, "sida1daaaaapppppsp_re_after", dare(sp_im)

print *, "sida1daaaaapppppsp_re_aftercp", dare(cp_re)
print *, "sida1daaaaapppppsp_re_aftercp", dare(cp_im)
!FOX  Y(1)=Y(1) - CP_RE*C1E3*MTCDA/(ONE+DPDA);
!FOX  Y(2)=Y(2) + CP_IM*C1E3*MTCDA/(ONE+DPDA);
print *, "sida1daaaaapppppsp_y1", dare(y(1))
print *, "sida1daaaaapppppsp_y2", dare(y(2))

print *, "ejjjjj111before", dare(ej1)
!FOX  EJ1 = EJ1 - (C1E3*E0F*SP_RE)/((CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI);



print *, "ejjjjj111", dare(ej1)

print *, "sida1daaaaapppppsp_y1", dare(y(1))
print *, "sida1daaaaapppppsp_y2", dare(y(2))



   



