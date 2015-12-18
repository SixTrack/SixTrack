      subroutine envada
!-----------------------------------------------------------------------
!  CALCULATION OF : MOMENTUM-DEPENDING ELEMENT-MATRICES AND
!                   CHANGE OF PATH LENGTHS FOR EACH PARTICLE.
!      SPECIALLY PREPARED FOR NEW D.A. (SIX-DIMENSIONAL VERSION)
!-----------------------------------------------------------------------
      implicit none
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
      integer i,ien,ih,ip,kz1,l,idaa
      double precision dare,result
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,    &
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,    &
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele=1200,nblo=600,nper=16,nelb=140,nblz=20000,         &
     &nzfz = 300000,mmul = 20)
      parameter(nran = 2000000,ncom = 100,mran = 500,mpa = 6,nrco = 5,  &
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 99)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 350)
      double precision c180e0,c1e1,c1e12,c1e13,c1e15,c1e16,c1e2,c1e3,   &
     &c1e4,c1e6,c1m1,c1m7,c1m10,c1m12,c1m13,c1m15,c1m18,c1m2,c1m21,     &
     &c1m24,c1m3,c1m36,c1m38,c1m6,c1m9,c2e3,c4e3,crade,clight,four,half,&
     &one,pieni,pmae,pmap,three,two,zero
      parameter(pieni = 1d-38)
      parameter(zero = 0.0d0,half = 0.5d0,one = 1.0d0)
      parameter(two = 2.0d0,three = 3.0d0,four = 4.0d0)
      parameter(c1e1 = 1.0d1,c1e2 = 1.0d2,c1m2 = 1.0d-2)
      parameter(c1e3 = 1.0d3,c2e3 = 2.0d3,c4e3 = 4.0d3,c1e4 = 1.0d4)
      parameter(c1e12 = 1.0d12,c1e13 = 1.0d13,c1e15 = 1.0d15,c1e16 =    &
     &1.0d16)
      parameter(c180e0 = 180.0d0,c1e6 = 1.0d6)
      parameter(c1m1 = 1.0d-1,c1m3 = 1.0d-3,c1m6 = 1.0d-6,c1m7 = 1.0d-7)
      parameter(c1m9 = 1.0d-9,c1m10 = 1.0d-10,c1m12 = 1.0d-12)
      parameter(c1m13 = 1.0d-13,c1m15 = 1.0d-15)
      parameter(c1m18 = 1.0d-18,c1m21 = 1.0d-21,c1m24 = 1.0d-24)
      parameter(c1m36 = 1.0d-36,c1m38 = 1.0d-38)
!     electron mass from PDG, 2002
      parameter(pmap = 938.271998d0,pmae = .510998902d0)
      parameter(crade = 2.817940285d-15, clight = 2.99792458d8)
      integer iav,ibb6d,ibbc,ibeco,ibidu,ibtyp,ic,icext,icextal,iclo,   &
     &iclo6,iclo6r,icode,icoe,icomb,icomb0,iconv,icow,icr,idam,idfor,   &
     &idis,idp,ierro,iffw,ifh,iicav,il,ilin,imad,imbb,                  &
     &imc,imtr,iorg,iout,                                               &
     &ipos,ipr,iprint,ipt,iq,iqmod,iqmod6,iratioe,ird,ire,ires,         &
     &irew,irm,irmod2,ise,ise1,ise2,ise3,isea,iskew,iskip,istw,         &
     &isub,itco,itcro,itf,ithick,ition,itionc,itqv,its6d,iu,iver,ivox,  &
     &ivoz,iwg,ixcav,izu0,kanf,kp,kpa,kwtype,kz,lhc,m21,m22,m23,mblo,   &
     &mbloz,mcut,mel,mesa,mmac,mout2,mp,mper,mstr,msym,mtyp,mzu,napx,   &
     &napxo,nbeam,nch,ncororb,ncorrep,ncorru,ncy,ndafi,nde,nhcorr,      &
     &nhmoni,niu,nlin,nmu,npp,nprint,nqc,nre,nrr,nskew,                 &
     &nstart,nstop,nt,nta,ntco,nte,ntwin,nu,numl,numlr,nur,nvcorr,      &
     &nvmoni,nwr, nturn1, nturn2, nturn3, nturn4,numlcp,numlmax,nnuml
      double precision a,ak0,aka,alfx,alfz,amp0,aper,apx,apz,ape,bbcu,  &
     &bclorb,beamoff,benkc,benki,betac,betam,betx,betz,bk0,bka,bl1,bl2, &
     &clo6,clobeam,clop6,cma1,cma2,cotr,crad,de0,dech,ded,dfft,         &
     &di0,dip0,dki,dkq,dma,dmap,dphix,dphiz,dppoff,dpscor,dqq,dres,dsi, &
     &dsm0,dtr,e0,ed,ej,ejf,ek,el,elbe,emitx,emity,emitz,extalign,      &
     &exterr,eui,euii,gammar,hsy,hsyc,pac,pam,parbe,parbe14,partnum,    &
     &phas,phas0,phasc,pi,pi2,pisqrt,pma,ptnfac,qs,qw0,qwsk,qx0,qxt,qz0,&
     &qzt,r00,rad,rat,ratio,ratioe,rrtr,rtc,rts,rvf,                    &
     &sigcor,sige,sigma0,sigman,sigman2,sigmanq,sigmoff,sigz,sm,ta,tam1,&
     &tam2,tiltc,tilts,tlen,totl,track6d,xpl,xrms,zfz,zpl,zrms,wirel,   &
     &acdipph, crabph, bbbx, bbby, bbbs,                                &
     &crabph2, crabph3, crabph4
      real hmal
      character*16 bez,bezb,bezr,erbez,bezl
      character*80 toptit,sixtit,commen
      common/erro/ierro,erbez
      common/kons/pi,pi2,pisqrt,rad
      common/str /il,mper,mblo,mbloz,msym(nper),kanf,iu,ic(nblz)
      common/ell /ed(nele),el(nele),ek(nele),sm(nele),kz(nele),kp(nele)
      common/bbb /bbbx(nele), bbby(nele), bbbs(nele)
      common/pla /xpl(nele),xrms(nele),zpl(nele),zrms(nele)
      common/str2 /mel(nblo),mtyp(nblo,nelb),mstr(nblo)
      common/mat/a(nele,2,6),bl1(nblo,2,6),bl2(nblo,2,6)
      common/syos2/rvf(mpa)
      common/tra1/rat,idfor,napx,napxo,numl,niu(2),numlr,nde(2),nwr(4), &
     &ird,imc,irew,ntwin,iclo6,iclo6r,iver,ibidu,numlcp,numlmax,nnuml
      common/syn/qs,e0,pma,ej(mpa),ejf(mpa),phas0,phas,hsy(3),          &
     &crad,hsyc(nele),phasc(nele),dppoff,sigmoff(nblz),tlen,            &
     &iicav,itionc(nele),ition,idp,ncy,ixcav
      common/corcom/dpscor,sigcor,icode,idam,its6d
      common/multi/bk0(nele,mmul),ak0(nele,mmul),                       &
     &bka(nele,mmul),aka(nele,mmul)
      common/mult1/benki,benkc(nele),r00(nele),irm(nele),nmu(nele)
      common/rand0/zfz(nzfz),iorg,mzu(nblz),bezr(3,nele),izu0,mmac,mcut
      common/rand1/exterr(nblz,40),extalign(nblz,3),tiltc(nblz),        &
     &tilts(nblz),mout2,icext(nblz),icextal(nblz)
      common/beo /aper(2),di0(2),dip0(2),ta(6,6)
      common/clo/dma,dmap,dkq,dqq,de0,ded,dsi,dech,dsm0,itco,itcro,itqv,&
     &iout
      common/qmodi/qw0(3),amp0,iq(3),iqmod,kpa(nele),iqmod6
      common/linop/bez(nele),elbe(nblo),bezb(nblo),ilin,nt,iprint,      &
     &ntco,eui,euii,nlin,bezl(nele)
      common/cororb/betam(nmon1,2),pam(nmon1,2),betac(ncor1,2),         &
     &pac(ncor1,2),bclorb(nmon1,2),nhmoni,nhcorr,nvmoni,nvcorr,         &
     &ncororb(nele)
      common/apert/apx(nele),apz(nele),ape(3,nele)
      common/clos/sigma0(2),iclo,ncorru,ncorrep
      common/combin/icomb0(20),icomb(ncom,20),ratio(ncom,20),           &
     &ratioe(nele),iratioe(nele),icoe
      common/seacom/ise,mesa,mp,m21,m22,m23,ise1,ise2,ise3,isea(nele)
      common/subres/qxt,qzt,tam1,tam2,isub,nta,nte,ipt,totl
      common/secom/rtc(9,18,10,5),rts(9,18,10,5),ire(12),ipr(5),irmod2
      common/secom1/dtr(10),nre,nur,nch,nqc,npp,nrr(5),nu(5)
      common/postr/dphix,dphiz,qx0,qz0,dres,dfft,cma1,cma2,             &
     &nstart,nstop,iskip,iconv,imad
      common/posti1/ipos,iav,iwg,ivox,ivoz,ires,ifh,toptit(5)
      common/posti2/kwtype,itf,icr,idis,icow,istw,iffw,nprint,ndafi
      common/skew/qwsk(2),betx(2),betz(2),alfx(2),alfz(2),iskew,nskew(6)
      common/pawc/hmal(nplo)
      common/tit/sixtit,commen,ithick
      common/co6d/clo6(3),clop6(3)
      common/dkic/dki(nele,3)
      common/beam/sigman(2,nbb),sigman2(2,nbb),sigmanq(2,nbb),          &
     &clobeam(6,nbb),beamoff(6,nbb),parbe(nele,5),track6d(6,npart),     &
     &ptnfac(nele),sigz,sige,partnum,parbe14,emitx,emity,emitz,gammar,  &
     &nbeam,ibbc,ibeco,ibtyp,lhc
      common/trom/ cotr(ntr,6),rrtr(ntr,6,6),imtr(nele)
      common/bb6d/ bbcu(nbb,12),ibb6d,imbb(nblz)
      common/wireco/ wirel(nele)
      common/acdipco/ acdipph(nele), nturn1(nele), nturn2(nele),        &
     &nturn3(nele), nturn4(nele)
      common/crabco/ crabph(nele),crabph2(nele),                        &
     &crabph3(nele),crabph4(nele)
      integer idz,itra
      double precision a2,al,as,at,chi0,chid,dp1,dps,exz,sigm
      common/syos/as(6,2,2,nele),at(6,2,2,nele),a2(6,2,2,nele),         &
     &al(6,2,2,nele),sigm(mpa),dps(mpa),idz(2)
      common/anf/chi0,chid,exz(2,6),dp1,itra
      integer ichrom,is
      double precision alf0,amp,bet0,clo,clop,cro,x,y
      common/tra/x(mpa,2),y(mpa,2),amp(2),bet0(2),alf0(2),clo(2),clop(2)
      common/chrom/cro(2),is(2),ichrom
      integer icorr,idial,idptr,imod1,imod2,inorm,ipar,namp,ncor,nctype,&
     &ndimf,nmom,nmom1,nmom2,nord,nord1,nordf,nsix,nvar,nvar2,nvarf
      double precision dpmax,preda,weig1,weig2
      character*16 coel
      common/dial/preda,idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
      common/norf/nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
      common/tcorr/icorr,nctype,namp,nmom,nmom1,nmom2,weig1,weig2,dpmax,&
     &coel(10)
      common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda, &
     &ej1,ejf1,rv
      double precision ald6,asd6
      common/dael6/ald6(nele,2,6,nema),asd6(nele,2,6,nema)
      save
!-----------------------------------------------------------------------
*FOX  B D ;
*FOX  D V DA EXT SIGMDA NORD NVAR ; D V DA EXT DPDA NORD NVAR ;
*FOX  D V DA EXT DPDA1 NORD NVAR ; D V DA EXT RV NORD NVAR ;
*FOX  D V DA EXT XX NORD NVAR 2 ; D V DA EXT YY NORD NVAR 2 ;
*FOX  D V DA EXT EJ1 NORD NVAR ; D V DA EXT EJF1 NORD NVAR ;
*FOX  D V DA EXT ALDA NORD NVAR 2 6 ; D V DA EXT ASDA NORD NVAR 2 6 ;
*FOX  D V DA EXT ALDAQ NORD NVAR 2 6 ; D V DA EXT ASDAQ NORD NVAR 2 6 ;
*FOX  D V DA EXT SMIDA NORD NVAR MCOR ;
*FOX  D V DA INT FOKQ NORD NVAR ; D V DA INT WFHI NORD NVAR ;
*FOX  D V DA INT DPD NORD NVAR ; D V DA INT DPSQ NORD NVAR ;
*FOX  D V DA INT FOK NORD NVAR ; D V DA INT RHO NORD NVAR ;
*FOX  D V DA INT FOK1 NORD NVAR ; D V DA INT SM1 NORD NVAR ;
*FOX  D V DA INT SM2 NORD NVAR ; D V DA INT SM3 NORD NVAR ;
*FOX  D V DA INT SM4 NORD NVAR ; D V DA INT SM5 NORD NVAR ;
*FOX  D V DA INT SM6 NORD NVAR ; D V DA INT SM12 NORD NVAR ;
*FOX  D V DA INT SM23 NORD NVAR ; D V DA INT AS3 NORD NVAR ;
*FOX  D V DA INT AS4 NORD NVAR ; D V DA INT AS6 NORD NVAR ;
*FOX  D V DA INT SI NORD NVAR ; D V DA INT CO NORD NVAR ;
*FOX  D V DA INT G NORD NVAR ; D V DA INT GL NORD NVAR ;
*FOX  D V DA INT SIQ NORD NVAR ; D V DA INT RHOC NORD NVAR ;
*FOX  D V DA INT HI NORD NVAR ; D V DA INT FI NORD NVAR ;
*FOX  D V DA INT AEK NORD NVAR ; D V DA INT HI1 NORD NVAR ;
*FOX  D V DA INT HP NORD NVAR ; D V DA INT HM NORD NVAR ;
*FOX  D V DA INT HC NORD NVAR ; D V DA INT HS NORD NVAR ;
*FOX  D V DA INT FOKC NORD NVAR ; D V DA INT WF NORD NVAR ;
*FOX  D V DA INT AFOK NORD NVAR ; D V DA INT RHOI NORD NVAR ;
*FOX  D V DA INT WFA NORD NVAR ; D V RE INT RATIOE NELE ;
*FOX  D V RE INT EL NELE ; D V RE INT EK NELE ; D V RE INT ED NELE ;
*FOX  D V RE INT ONE ; D V RE INT ZERO ; D V RE INT TWO ;
*FOX  D V RE INT HALF ; D V RE INT FOUR ; D V RE INT C1E3 ;
*FOX  D V RE INT C2E3 ; D V RE INT C4E3 ;
*FOX  D V IN INT I ; D V IN INT L ; D V IN INT IH ; D V IN INT NE ;
*FOX  D V IN INT NA ; D V IN INT IP ; D V IN INT IPCH ;
*FOX  D F RE DARE 1 ;
*FOX  E D ;
*FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
*FOX  DPD=ONE+DPDA ;
*FOX  DPSQ=SQRT(DPD) ;
      do 220 i=1,il
        do 10 ih=1,2
          do 10 ip=1,6
*FOX  ALDA(IH,IP)=ZERO ;
*FOX  ASDA(IH,IP)=ZERO ;
   10   continue
        if(abs(el(i)).le.pieni) goto 190
        kz1=kz(i)+1
!       goto(20,40,100,60,80,90,130,170,180),kz1
        if (kz1.eq.1) goto 20
        if (kz1.eq.2) goto 40
        if (kz1.eq.3) goto 100
        if (kz1.eq.4) goto 60
        if (kz1.eq.5) goto 80
        if (kz1.eq.6) goto 90
        if (kz1.eq.7) goto 130
        if (kz1.eq.8) goto 170
        if (kz1.eq.9) goto 180
        goto 220
!-----------------------------------------------------------------------
!  DRIFTLENGTH
!-----------------------------------------------------------------------
   20   do 30 l=1,2
*FOX  ALDA(L,1)=ONE  ;
*FOX  ALDA(L,2)=EL(I) ;
*FOX  ALDA(L,3)=ZERO ;
*FOX  ALDA(L,4)=ONE ;
*FOX  ASDA(L,6)=-RV*ALDA(L,2)/C2E3 ;
   30   continue
*FOX  ASDA(1,1)=EL(I)*(ONE-RV)*C1E3 ;
        goto 190
!-----------------------------------------------------------------------
!  RECTANGULAR MAGNET
!  HORIZONTAL
!-----------------------------------------------------------------------
   40   ih=1
   50   continue
        if(abs(ed(i)).le.pieni) goto 20
*FOX  FOK=EL(I)*ED(I)/DPSQ ;
*FOX  RHO=ONE/ED(I)*DPSQ ;
*FOX  FOK1=SIN(FOK*HALF)/COS(FOK*HALF)/RHO ;
*FOX  SI=SIN(FOK) ;
*FOX  CO=COS(FOK) ;
*FOX  ALDA(IH,1)=ONE ;
*FOX  ALDA(IH,2)=RHO*SI ;
*FOX  ALDA(IH,3)=ZERO ;
*FOX  ALDA(IH,4)=ONE ;
*FOX  ALDA(IH,5)=-DPDA*RHO*(ONE-CO)/DPSQ*C1E3 ;
*FOX  ALDA(IH,6)=-DPDA*TWO*SIN(FOK*HALF)/COS(FOK*HALF)/DPSQ*C1E3 ;
*FOX  SM1=COS(FOK) ;
*FOX  SM2=SIN(FOK)*RHO ;
*FOX  SM3=-SIN(FOK)/RHO ;
*FOX  SM5=-RHO*DPSQ*(ONE-SM1) ;
*FOX  SM6=-SM2*DPSQ/RHO ;
*FOX  SM12=EL(I)-SM1*SM2 ;
*FOX  SM23=SM2*SM3 ;
*FOX  AS3=-RV*(DPDA*RHO/(TWO*DPSQ)*SM23+SM5) ;
*FOX  AS4=-RV*SM23/C2E3 ;
*FOX  AS6=-RV*(EL(I)+SM1*SM2)/C4E3 ;
*FOX  ASDA(IH,1)=(-RV*(DPDA*DPDA/(FOUR*DPD)*SM12+DPDA*(EL(I)-SM2))
*FOX  +EL(I)*(ONE-RV))*C1E3 ;
*FOX  ASDA(IH,2)=-RV*(DPDA/(TWO*RHO*DPSQ)*SM12+SM6)+FOK1*AS3 ;
*FOX  ASDA(IH,3)=AS3 ;
*FOX  ASDA(IH,4)=AS4+TWO*AS6*FOK1 ;
*FOX  ASDA(IH,5)=-RV*SM12/(C4E3*RHO*RHO)+AS6*FOK1*FOK1+FOK1*AS4  ;
*FOX  ASDA(IH,6)=AS6 ;
!--VERTIKAL
        ih=ih+1
        if(ih.gt.2) ih=1
*FOX  G=SIN(FOK*HALF)/COS(FOK*HALF)/RHO ;
*FOX  GL=EL(I)*G ;
*FOX  ALDA(IH,1)=ONE-GL ;
*FOX  ALDA(IH,2)=EL(I) ;
*FOX  ALDA(IH,3)=-G*(TWO-GL) ;
*FOX  ALDA(IH,4)=ALDA(IH,1) ;
*FOX  AS6=-RV*ALDA(IH,2)/C2E3 ;
*FOX  ASDA(IH,4)=-TWO*AS6*FOK1 ;
*FOX  ASDA(IH,5)=AS6*FOK1*FOK1 ;
*FOX  ASDA(IH,6)=AS6 ;
        goto 190
!-----------------------------------------------------------------------
!  SEKTORMAGNET
!  HORIZONTAL
!-----------------------------------------------------------------------
   60   ih=1
   70   continue
        if(abs(ed(i)).le.pieni) goto 20
*FOX  FOK=EL(I)*ED(I)/DPSQ ;
*FOX  RHO=(ONE/ED(I))*DPSQ ;
*FOX  SI=SIN(FOK) ;
*FOX  CO=COS(FOK) ;
*FOX  RHOC=RHO*(ONE-CO)/DPSQ ;
*FOX  SIQ=SI/DPSQ ;
*FOX  ALDA(IH,1)=CO ;
*FOX  ALDA(IH,2)=RHO*SI ;
*FOX  ALDA(IH,3)=-SI/RHO ;
*FOX  ALDA(IH,4)=CO ;
*FOX  ALDA(IH,5)=-DPDA*RHOC*C1E3 ;
*FOX  ALDA(IH,6)=-DPDA*SIQ*C1E3 ;
*FOX  SM12=EL(I)-ALDA(IH,1)*ALDA(IH,2) ;
*FOX  SM23=ALDA(IH,2)*ALDA(IH,3) ;
*FOX  ASDA(IH,1)=(-RV*(DPDA*DPDA/(FOUR*DPD)*SM12
*FOX  +DPDA*(EL(I)-ALDA(IH,2)))+EL(I)*(ONE-RV))*C1E3 ;
*FOX  ASDA(IH,2)=-RV*(DPDA/(TWO*RHO*DPSQ)*SM12-DPD*SIQ) ;
*FOX  ASDA(IH,3)=-RV*(DPDA*RHO/(TWO*DPSQ)*SM23-DPD*RHOC) ;
*FOX  ASDA(IH,4)=-RV*SM23/C2E3 ;
*FOX  ASDA(IH,5)=-RV*SM12/(C4E3*RHO*RHO) ;
*FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;
!--VERTIKAL
        ih=ih+1
        if(ih.gt.2) ih=1
*FOX  ALDA(IH,1)=ONE ;
*FOX  ALDA(IH,2)=EL(I) ;
*FOX  ALDA(IH,3)=ZERO ;
*FOX  ALDA(IH,4)=ONE ;
*FOX  ASDA(IH,6)=-RV*ALDA(IH,2)/C2E3 ;
        goto 190
!-----------------------------------------------------------------------
!  RECTANGULAR MAGNET VERTIKAL
!-----------------------------------------------------------------------
   80   ih=2
        goto 50
!-----------------------------------------------------------------------
!  SEKTORMAGNET VERTIKAL
!-----------------------------------------------------------------------
   90   ih=2
        goto 70
!-----------------------------------------------------------------------
!  QUADRUPOLE
!  FOCUSSING
!-----------------------------------------------------------------------
  100   continue
        if(abs(ek(i)).le.pieni) goto 20
*FOX  FOK=EK(I)/(ONE+DPDA) ;
*FOX  AEK=FOK ;
        if(dare(aek).lt.zero) then
*FOX  AEK=-AEK ;
        endif
        ih=0
*FOX  HI=SQRT(AEK) ;
*FOX  FI=EL(I)*HI ;
        if(ek(i).gt.zero) goto 120
  110   ih=ih+1
*FOX  ALDA(IH,1)=COS(FI) ;
*FOX  HI1=SIN(FI) ;
*FOX  ALDA(IH,2)=HI1/HI ;
*FOX  ALDA(IH,3)=-HI1*HI ;
*FOX  ALDA(IH,4)=ALDA(IH,1) ;
*FOX  ASDA(IH,1)=EL(I)*(ONE-RV)*C1E3 ;
*FOX  ASDA(IH,4)=-RV*ALDA(IH,2)*ALDA(IH,3)/C2E3 ;
*FOX  ASDA(IH,5)=-RV*(EL(I)-ALDA(IH,1)*ALDA(IH,2))*AEK/C4E3 ;
*FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;
        if(ih.eq.2) goto 190
!--DEFOCUSSING
  120   ih=ih+1
*FOX  HP=EXP(FI) ;
*FOX  HM=ONE/HP ;
*FOX  HC=(HP+HM)*HALF ;
*FOX  HS=(HP-HM)*HALF ;
*FOX  ALDA(IH,1)=HC ;
*FOX  ALDA(IH,2)=HS/HI ;
*FOX  ALDA(IH,3)=HS*HI ;
*FOX  ALDA(IH,4)=HC ;
*FOX  ASDA(IH,4)=-RV*ALDA(IH,2)*ALDA(IH,3)/C2E3 ;
*FOX  ASDA(IH,5)=+RV*(EL(I)-ALDA(IH,1)*ALDA(IH,2))*AEK/C4E3 ;
*FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;
        if(ih.eq.1) goto 110
        goto 190
!-----------------------------------------------------------------------
!  COMBINED FUNCTION MAGNET HORIZONTAL
!  FOCUSSING
!-----------------------------------------------------------------------
  130   ih=0
*FOX  FOKQ=EK(I) ;
  140   continue
        if(abs(ek(i)).le.pieni) goto 60
        if(abs(ed(i)).le.pieni) goto 100
!hr08   if(abs(ek(i)-ed(i)*ed(i)).le.pieni) goto 20
        if(abs(ek(i)-ed(i)**2).le.pieni) goto 20                         !hr08
*FOX  WF=ED(I)/DPSQ ;
*FOX  FOK=FOKQ/DPD-WF*WF ;
*FOX  AFOK=FOK ;
      if(dare(afok).lt.zero) then
*FOX  AFOK=-AFOK ;
      endif
*FOX  HI=SQRT(AFOK) ;
*FOX  FI=HI*EL(I) ;
        if(dare(fok).gt.zero) goto 160
  150   ih=ih+1
*FOX  SI=SIN(FI) ;
*FOX  CO=COS(FI) ;
*FOX  WFA=WF/AFOK*(ONE-CO)/DPSQ ;
*FOX  WFHI=WF/HI*SI/DPSQ ;
*FOX  ALDA(IH,1)=CO ;
*FOX  ALDA(IH,2)=SI/HI ;
*FOX  ALDA(IH,3)=-SI*HI ;
*FOX  ALDA(IH,4)=CO ;
*FOX  ALDA(IH,5)=-WFA*DPDA*C1E3 ;
*FOX  ALDA(IH,6)=-WFHI*DPDA*C1E3 ;
*FOX  SM12=EL(I)-ALDA(IH,1)*ALDA(IH,2) ;
*FOX  SM23=ALDA(IH,2)*ALDA(IH,3) ;
*FOX  ASDA(IH,1)=(-RV*(DPDA*DPDA/(FOUR*DPD)*SM12+DPDA
*FOX  *(EL(I)-ALDA(IH,2)))/AFOK*WF*WF+EL(I)*(ONE-RV))*C1E3 ;
*FOX  ASDA(IH,2)=-RV*(DPDA*WF/(TWO*DPSQ)*SM12-DPD*WFHI) ;
*FOX  ASDA(IH,3)=-RV*(DPDA*HALF/AFOK/DPD*ED(I)*SM23-DPD*WFA) ;
*FOX  ASDA(IH,4)=-RV*SM23/C2E3 ;
*FOX  ASDA(IH,5)=-RV*SM12*AFOK/C4E3 ;
*FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;
        ih=ih+1
        if(ih.gt.2) ih=1
*FOX  AEK=EK(I)/DPD ;
      if(dare(aek).lt.zero) then
*FOX  AEK=-AEK ;
      endif
*FOX  HI=SQRT(AEK) ;
*FOX  FI=HI*EL(I) ;
*FOX  HP=EXP(FI) ;
*FOX  HM=ONE/HP ;
*FOX  HC=(HP+HM)*HALF ;
*FOX  HS=(HP-HM)*HALF ;
*FOX  ALDA(IH,1)=HC ;
*FOX  ALDA(IH,2)=EL(I) ;
*FOX  ALDA(IH,2)=HS/HI ;
*FOX  ALDA(IH,3)=HS*HI ;
*FOX  ALDA(IH,4)=HC ;
*FOX  ASDA(IH,4)=-RV*ALDA(IH,2)*ALDA(IH,3)/C2E3 ;
*FOX  ASDA(IH,5)=+RV*(EL(I)-ALDA(IH,1)*ALDA(IH,2))*AEK/C4E3 ;
*FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;
        goto 190
!--DEFOCUSSING
  160   ih=ih+1
*FOX  HP=EXP(FI) ;
*FOX  HM=ONE/HP ;
*FOX  HC=(HP+HM)*HALF ;
*FOX  HS=(HP-HM)*HALF ;
*FOX  ALDA(IH,1)=HC ;
*FOX  ALDA(IH,2)=HS/HI ;
*FOX  ALDA(IH,3)=HS*HI ;
*FOX  ALDA(IH,4)=HC ;
*FOX  WFA=WF/AFOK*(ONE-HC)/DPSQ ;
*FOX  WFHI=WF/HI*HS/DPSQ ;
*FOX  ALDA(IH,5)= WFA*DPDA*C1E3 ;
*FOX  ALDA(IH,6)=-WFHI*DPDA*C1E3 ;
*FOX  SM12=EL(I)-ALDA(IH,1)*ALDA(IH,2) ;
*FOX  SM23=ALDA(IH,2)*ALDA(IH,3) ;
*FOX  ASDA(IH,1)=(RV*(DPDA*DPDA/(FOUR*DPD)*SM12
*FOX  +DPDA*(EL(I)-ALDA(IH,2)))/AFOK*WF*WF+EL(I)*(ONE-RV))*C1E3 ;
*FOX  ASDA(IH,2)=-RV*(DPDA*WF/(TWO*DPSQ)*SM12-DPD*WFHI) ;
*FOX  ASDA(IH,3)=RV*(DPDA*HALF/AFOK/DPD*ED(I)*SM23-DPD*WFA) ;
*FOX  ASDA(IH,4)=-RV*SM23/C2E3 ;
*FOX  ASDA(IH,5)=+RV*SM12*AFOK/C4E3 ;
*FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;
        ih=ih+1
        if(ih.gt.2) ih=1
*FOX  AEK=EK(I)/DPD ;
      if(dare(aek).lt.zero) then
*FOX  AEK=-AEK ;
      endif
*FOX  HI=SQRT(AEK) ;
*FOX  FI=HI*EL(I) ;
*FOX  SI=SIN(FI) ;
*FOX  CO=COS(FI) ;
*FOX  ALDA(IH,1)=CO ;
*FOX  ALDA(IH,2)=SI/HI ;
*FOX  ALDA(IH,3)=-SI*HI ;
*FOX  ALDA(IH,4)=CO ;
*FOX  ASDA(IH,4)=-RV*ALDA(IH,2)*ALDA(IH,3)/C2E3 ;
*FOX  ASDA(IH,5)=-RV*(EL(I)-ALDA(IH,1)*ALDA(IH,2))*AEK/C4E3 ;
*FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;
        goto 190
!-----------------------------------------------------------------------
!  COMBINED FUNCTION MAGNET VERTICAL
!-----------------------------------------------------------------------
  170   ih=1
*FOX  FOKQ=-EK(I) ;
        goto 140
!-----------------------------------------------------------------------
!  EDGE FOCUSSING
!-----------------------------------------------------------------------
  180   continue
*FOX  RHOI=ED(I)/DPSQ ;
*FOX  FOK=RHOI*SIN(EL(I)*RHOI*HALF)/COS(EL(I)*RHOI*HALF) ;
*FOX  ALDA(1,1)=ONE ;
*FOX  ALDA(1,2)=ZERO ;
*FOX  ALDA(1,3)=FOK ;
*FOX  ALDA(1,4)=ONE ;
*FOX  ALDA(2,1)=ONE ;
*FOX  ALDA(2,2)=ZERO ;
*FOX  ALDA(2,3)=-FOK ;
*FOX  ALDA(2,4)=ONE ;
        goto 190
!-----------------------------------------------------------------------
!   NONLINEAR INSERTION
!-----------------------------------------------------------------------
  190   continue
        do 210 ih=1,2
          do 210 ip=1,6
            do 200 ien=1,nord+1
              if (nvar2.eq.5) then
                call dapri6(alda(ih,ip),result,ien,5)
                ald6(i,ih,ip,ien) = result
                call dapri6(asda(ih,ip),result,ien,5)
                asd6(i,ih,ip,ien) = result
              else if (nvar2.eq.6) then
                call dapri6(alda(ih,ip),result,ien,6)
                ald6(i,ih,ip,ien) = result
                call dapri6(asda(ih,ip),result,ien,6)
                asd6(i,ih,ip,ien) = result
              else if (nvar2.eq.4) then
                call dapri6(alda(ih,ip),result,ien,4)
                ald6(i,ih,ip,ien) = result
                call dapri6(asda(ih,ip),result,ien,4)
                asd6(i,ih,ip,ien) = result
              endif
  200       continue
  210   continue
  220 continue
!     DADAL AUTOMATIC INCLUSION
      return
      end
      subroutine envquad(i,ipch)
!-----------------------------------------------------------------------
!  CALCULATION OF : MOMENTUM-DEPENDING ELEMENT-MATRICES AND
!                   CHANGE OF PATH LENGTHS FOR EACH PARTICLE.
!      SPECIALLY PREPARED FOR NEW D.A. (SIX-DIMENSIONAL VERSION)
!-----------------------------------------------------------------------
      implicit none
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
      integer i,ih,ipch,idaa
      double precision dare
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,    &
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,    &
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele=1200,nblo=600,nper=16,nelb=140,nblz=20000,         &
     &nzfz = 300000,mmul = 20)
      parameter(nran = 2000000,ncom = 100,mran = 500,mpa = 6,nrco = 5,  &
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 99)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 350)
      double precision c180e0,c1e1,c1e12,c1e13,c1e15,c1e16,c1e2,c1e3,   &
     &c1e4,c1e6,c1m1,c1m7,c1m10,c1m12,c1m13,c1m15,c1m18,c1m2,c1m21,     &
     &c1m24,c1m3,c1m36,c1m38,c1m6,c1m9,c2e3,c4e3,crade,clight,four,half,&
     &one,pieni,pmae,pmap,three,two,zero
      parameter(pieni = 1d-38)
      parameter(zero = 0.0d0,half = 0.5d0,one = 1.0d0)
      parameter(two = 2.0d0,three = 3.0d0,four = 4.0d0)
      parameter(c1e1 = 1.0d1,c1e2 = 1.0d2,c1m2 = 1.0d-2)
      parameter(c1e3 = 1.0d3,c2e3 = 2.0d3,c4e3 = 4.0d3,c1e4 = 1.0d4)
      parameter(c1e12 = 1.0d12,c1e13 = 1.0d13,c1e15 = 1.0d15,c1e16 =    &
     &1.0d16)
      parameter(c180e0 = 180.0d0,c1e6 = 1.0d6)
      parameter(c1m1 = 1.0d-1,c1m3 = 1.0d-3,c1m6 = 1.0d-6,c1m7 = 1.0d-7)
      parameter(c1m9 = 1.0d-9,c1m10 = 1.0d-10,c1m12 = 1.0d-12)
      parameter(c1m13 = 1.0d-13,c1m15 = 1.0d-15)
      parameter(c1m18 = 1.0d-18,c1m21 = 1.0d-21,c1m24 = 1.0d-24)
      parameter(c1m36 = 1.0d-36,c1m38 = 1.0d-38)
!     electron mass from PDG, 2002
      parameter(pmap = 938.271998d0,pmae = .510998902d0)
      parameter(crade = 2.817940285d-15, clight = 2.99792458d8)
      integer iav,ibb6d,ibbc,ibeco,ibidu,ibtyp,ic,icext,icextal,iclo,   &
     &iclo6,iclo6r,icode,icoe,icomb,icomb0,iconv,icow,icr,idam,idfor,   &
     &idis,idp,ierro,iffw,ifh,iicav,il,ilin,imad,imbb,                  &
     &imc,imtr,iorg,iout,                                               &
     &ipos,ipr,iprint,ipt,iq,iqmod,iqmod6,iratioe,ird,ire,ires,         &
     &irew,irm,irmod2,ise,ise1,ise2,ise3,isea,iskew,iskip,istw,         &
     &isub,itco,itcro,itf,ithick,ition,itionc,itqv,its6d,iu,iver,ivox,  &
     &ivoz,iwg,ixcav,izu0,kanf,kp,kpa,kwtype,kz,lhc,m21,m22,m23,mblo,   &
     &mbloz,mcut,mel,mesa,mmac,mout2,mp,mper,mstr,msym,mtyp,mzu,napx,   &
     &napxo,nbeam,nch,ncororb,ncorrep,ncorru,ncy,ndafi,nde,nhcorr,      &
     &nhmoni,niu,nlin,nmu,npp,nprint,nqc,nre,nrr,nskew,                 &
     &nstart,nstop,nt,nta,ntco,nte,ntwin,nu,numl,numlr,nur,nvcorr,      &
     &nvmoni,nwr, nturn1, nturn2, nturn3, nturn4,numlcp,numlmax,nnuml
      double precision a,ak0,aka,alfx,alfz,amp0,aper,apx,apz,ape,bbcu,  &
     &bclorb,beamoff,benkc,benki,betac,betam,betx,betz,bk0,bka,bl1,bl2, &
     &clo6,clobeam,clop6,cma1,cma2,cotr,crad,de0,dech,ded,dfft,         &
     &di0,dip0,dki,dkq,dma,dmap,dphix,dphiz,dppoff,dpscor,dqq,dres,dsi, &
     &dsm0,dtr,e0,ed,ej,ejf,ek,el,elbe,emitx,emity,emitz,extalign,      &
     &exterr,eui,euii,gammar,hsy,hsyc,pac,pam,parbe,parbe14,partnum,    &
     &phas,phas0,phasc,pi,pi2,pisqrt,pma,ptnfac,qs,qw0,qwsk,qx0,qxt,qz0,&
     &qzt,r00,rad,rat,ratio,ratioe,rrtr,rtc,rts,rvf,                    &
     &sigcor,sige,sigma0,sigman,sigman2,sigmanq,sigmoff,sigz,sm,ta,tam1,&
     &tam2,tiltc,tilts,tlen,totl,track6d,xpl,xrms,zfz,zpl,zrms,wirel,   &
     &acdipph, crabph, bbbx, bbby, bbbs,                                &
     &crabph2, crabph3, crabph4
      real hmal
      character*16 bez,bezb,bezr,erbez,bezl
      character*80 toptit,sixtit,commen
      common/erro/ierro,erbez
      common/kons/pi,pi2,pisqrt,rad
      common/str /il,mper,mblo,mbloz,msym(nper),kanf,iu,ic(nblz)
      common/ell /ed(nele),el(nele),ek(nele),sm(nele),kz(nele),kp(nele)
      common/bbb /bbbx(nele), bbby(nele), bbbs(nele)
      common/pla /xpl(nele),xrms(nele),zpl(nele),zrms(nele)
      common/str2 /mel(nblo),mtyp(nblo,nelb),mstr(nblo)
      common/mat/a(nele,2,6),bl1(nblo,2,6),bl2(nblo,2,6)
      common/syos2/rvf(mpa)
      common/tra1/rat,idfor,napx,napxo,numl,niu(2),numlr,nde(2),nwr(4), &
     &ird,imc,irew,ntwin,iclo6,iclo6r,iver,ibidu,numlcp,numlmax,nnuml
      common/syn/qs,e0,pma,ej(mpa),ejf(mpa),phas0,phas,hsy(3),          &
     &crad,hsyc(nele),phasc(nele),dppoff,sigmoff(nblz),tlen,            &
     &iicav,itionc(nele),ition,idp,ncy,ixcav
      common/corcom/dpscor,sigcor,icode,idam,its6d
      common/multi/bk0(nele,mmul),ak0(nele,mmul),                       &
     &bka(nele,mmul),aka(nele,mmul)
      common/mult1/benki,benkc(nele),r00(nele),irm(nele),nmu(nele)
      common/rand0/zfz(nzfz),iorg,mzu(nblz),bezr(3,nele),izu0,mmac,mcut
      common/rand1/exterr(nblz,40),extalign(nblz,3),tiltc(nblz),        &
     &tilts(nblz),mout2,icext(nblz),icextal(nblz)
      common/beo /aper(2),di0(2),dip0(2),ta(6,6)
      common/clo/dma,dmap,dkq,dqq,de0,ded,dsi,dech,dsm0,itco,itcro,itqv,&
     &iout
      common/qmodi/qw0(3),amp0,iq(3),iqmod,kpa(nele),iqmod6
      common/linop/bez(nele),elbe(nblo),bezb(nblo),ilin,nt,iprint,      &
     &ntco,eui,euii,nlin,bezl(nele)
      common/cororb/betam(nmon1,2),pam(nmon1,2),betac(ncor1,2),         &
     &pac(ncor1,2),bclorb(nmon1,2),nhmoni,nhcorr,nvmoni,nvcorr,         &
     &ncororb(nele)
      common/apert/apx(nele),apz(nele),ape(3,nele)
      common/clos/sigma0(2),iclo,ncorru,ncorrep
      common/combin/icomb0(20),icomb(ncom,20),ratio(ncom,20),           &
     &ratioe(nele),iratioe(nele),icoe
      common/seacom/ise,mesa,mp,m21,m22,m23,ise1,ise2,ise3,isea(nele)
      common/subres/qxt,qzt,tam1,tam2,isub,nta,nte,ipt,totl
      common/secom/rtc(9,18,10,5),rts(9,18,10,5),ire(12),ipr(5),irmod2
      common/secom1/dtr(10),nre,nur,nch,nqc,npp,nrr(5),nu(5)
      common/postr/dphix,dphiz,qx0,qz0,dres,dfft,cma1,cma2,             &
     &nstart,nstop,iskip,iconv,imad
      common/posti1/ipos,iav,iwg,ivox,ivoz,ires,ifh,toptit(5)
      common/posti2/kwtype,itf,icr,idis,icow,istw,iffw,nprint,ndafi
      common/skew/qwsk(2),betx(2),betz(2),alfx(2),alfz(2),iskew,nskew(6)
      common/pawc/hmal(nplo)
      common/tit/sixtit,commen,ithick
      common/co6d/clo6(3),clop6(3)
      common/dkic/dki(nele,3)
      common/beam/sigman(2,nbb),sigman2(2,nbb),sigmanq(2,nbb),          &
     &clobeam(6,nbb),beamoff(6,nbb),parbe(nele,5),track6d(6,npart),     &
     &ptnfac(nele),sigz,sige,partnum,parbe14,emitx,emity,emitz,gammar,  &
     &nbeam,ibbc,ibeco,ibtyp,lhc
      common/trom/ cotr(ntr,6),rrtr(ntr,6,6),imtr(nele)
      common/bb6d/ bbcu(nbb,12),ibb6d,imbb(nblz)
      common/wireco/ wirel(nele)
      common/acdipco/ acdipph(nele), nturn1(nele), nturn2(nele),        &
     &nturn3(nele), nturn4(nele)
      common/crabco/ crabph(nele),crabph2(nele),                        &
     &crabph3(nele),crabph4(nele)
      integer icorr,idial,idptr,imod1,imod2,inorm,ipar,namp,ncor,nctype,&
     &ndimf,nmom,nmom1,nmom2,nord,nord1,nordf,nsix,nvar,nvar2,nvarf
      double precision dpmax,preda,weig1,weig2
      character*16 coel
      common/dial/preda,idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
      common/norf/nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
      common/tcorr/icorr,nctype,namp,nmom,nmom1,nmom2,weig1,weig2,dpmax,&
     &coel(10)
      common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda, &
     &ej1,ejf1,rv
      save
!-----------------------------------------------------------------------
*FOX  B D ;
*FOX  D V DA EXT SIGMDA NORD NVAR ; D V DA EXT DPDA NORD NVAR ;
*FOX  D V DA EXT DPDA1 NORD NVAR ; D V DA EXT RV NORD NVAR ;
*FOX  D V DA EXT XX NORD NVAR 2 ; D V DA EXT YY NORD NVAR 2 ;
*FOX  D V DA EXT EJ1 NORD NVAR ; D V DA EXT EJF1 NORD NVAR ;
*FOX  D V DA EXT ALDA NORD NVAR 2 6 ; D V DA EXT ASDA NORD NVAR 2 6 ;
*FOX  D V DA EXT ALDAQ NORD NVAR 2 6 ; D V DA EXT ASDAQ NORD NVAR 2 6 ;
*FOX  D V DA EXT SMIDA NORD NVAR MCOR ;
*FOX  D V DA INT FOKQ NORD NVAR ; D V DA INT WFHI NORD NVAR ;
*FOX  D V DA INT DPD NORD NVAR ; D V DA INT DPSQ NORD NVAR ;
*FOX  D V DA INT FOK NORD NVAR ; D V DA INT RHO NORD NVAR ;
*FOX  D V DA INT FOK1 NORD NVAR ; D V DA INT SM1 NORD NVAR ;
*FOX  D V DA INT SM2 NORD NVAR ; D V DA INT SM3 NORD NVAR ;
*FOX  D V DA INT SM4 NORD NVAR ; D V DA INT SM5 NORD NVAR ;
*FOX  D V DA INT SM6 NORD NVAR ; D V DA INT SM12 NORD NVAR ;
*FOX  D V DA INT SM23 NORD NVAR ; D V DA INT AS3 NORD NVAR ;
*FOX  D V DA INT AS4 NORD NVAR ; D V DA INT AS6 NORD NVAR ;
*FOX  D V DA INT SI NORD NVAR ; D V DA INT CO NORD NVAR ;
*FOX  D V DA INT G NORD NVAR ; D V DA INT GL NORD NVAR ;
*FOX  D V DA INT SIQ NORD NVAR ; D V DA INT RHOC NORD NVAR ;
*FOX  D V DA INT HI NORD NVAR ; D V DA INT FI NORD NVAR ;
*FOX  D V DA INT AEK NORD NVAR ; D V DA INT HI1 NORD NVAR ;
*FOX  D V DA INT HP NORD NVAR ; D V DA INT HM NORD NVAR ;
*FOX  D V DA INT HC NORD NVAR ; D V DA INT HS NORD NVAR ;
*FOX  D V DA INT FOKC NORD NVAR ; D V DA INT WF NORD NVAR ;
*FOX  D V DA INT AFOK NORD NVAR ; D V DA INT RHOI NORD NVAR ;
*FOX  D V DA INT WFA NORD NVAR ; D V RE INT RATIOE NELE ;
*FOX  D V RE INT EL NELE ; D V RE INT EK NELE ; D V RE INT ED NELE ;
*FOX  D V RE INT ONE ; D V RE INT ZERO ; D V RE INT TWO ;
*FOX  D V RE INT HALF ; D V RE INT FOUR ; D V RE INT C1E3 ;
*FOX  D V RE INT C2E3 ; D V RE INT C4E3 ;
*FOX  D V IN INT I ; D V IN INT L ; D V IN INT IH ; D V IN INT NE ;
*FOX  D V IN INT NA ; D V IN INT IP ; D V IN INT IPCH ;
*FOX  D F RE DARE 1 ;
*FOX  E D ;
*FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
*FOX  DPD=ONE+DPDA ;
*FOX  DPSQ=SQRT(DPD) ;
!-----------------------------------------------------------------------
!  QUADRUPOLE
!  FOCUSSING
!-----------------------------------------------------------------------
      if(abs(ek(i)).le.pieni) goto 100
*FOX  FOK=(SMIDA(IPCH)*RATIOE(I))/(ONE+DPDA) ;
*FOX  AEK=FOK ;
      if(dare(aek).lt.zero) then
*FOX  AEK=-AEK ;
      endif
      ih=0
*FOX  HI=SQRT(AEK) ;
*FOX  FI=EL(I)*HI ;
      if(ek(i).gt.zero) goto 30
   20 ih=ih+1
*FOX  ALDAQ(IH,1)=COS(FI) ;
*FOX  HI1=SIN(FI) ;
*FOX  ALDAQ(IH,2)=HI1/HI ;
*FOX  ALDAQ(IH,3)=-HI1*HI ;
*FOX  ALDAQ(IH,4)=ALDAQ(IH,1) ;
*FOX  ASDAQ(IH,1)=EL(I)*(ONE-RV)*C1E3 ;
*FOX  ASDAQ(IH,4)=-RV*ALDAQ(IH,2)*ALDAQ(IH,3)/C2E3 ;
*FOX  ASDAQ(IH,5)=-RV*(EL(I)-ALDAQ(IH,1)*ALDAQ(IH,2))*AEK/C4E3 ;
*FOX  ASDAQ(IH,6)=-RV*(EL(I)+ALDAQ(IH,1)*ALDAQ(IH,2))/C4E3 ;
      if(ih.eq.2) goto 100
!--DEFOCUSSING
   30 ih=ih+1
*FOX  HP=EXP(FI) ;
*FOX  HM=ONE/HP ;
*FOX  HC=(HP+HM)*HALF ;
*FOX  HS=(HP-HM)*HALF ;
*FOX  ALDAQ(IH,1)=HC ;
*FOX  ALDAQ(IH,2)=HS/HI ;
*FOX  ALDAQ(IH,3)=HS*HI ;
*FOX  ALDAQ(IH,4)=HC ;
*FOX  ASDAQ(IH,4)=-RV*ALDAQ(IH,2)*ALDAQ(IH,3)/C2E3 ;
*FOX  ASDAQ(IH,5)=+RV*(EL(I)-ALDAQ(IH,1)*ALDAQ(IH,2))*AEK/C4E3 ;
*FOX  ASDAQ(IH,6)=-RV*(EL(I)+ALDAQ(IH,1)*ALDAQ(IH,2))/C4E3 ;
      if(ih.eq.1) goto 20
  100 continue
!     DADAL AUTOMATIC INCLUSION
      return
      end
      subroutine umlauda
!-----------------------------------------------------------------------
!  CENTRAL LOOP FOR 6-DIMENSIONAL CLOSED ORBIT
!-----------------------------------------------------------------------
      implicit none
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
      integer i,ibb,iii,i2,i3,i4,icav,icoonly,ien,iflag,iflag1,iflag2,  &
     &ii,ii2,ip,ipch,irrtr,ivar,ivar1,iwrite,ix,j,j1,jb,jj,jmel,jx,k,   &
     &kkk,kpz,kzz,mfile,nd2,nmz,idaa,angno,damap,damapi,damap1,f,aa2,   &
     &aa2r,a1,a1r,xy,h,df
      double precision al1,al2,al3,angp,angnoe,au,aui,b1,b2,b3,beamoff1,&
     &beamoff2,beamoff4,beamoff5,beamoff6,betr0,c,c5m4,cbxb,cbzb,coefh1,&
     &cik,coefh2,coefv1,coefv2,cp,crk,crxb,crzb,cx,d,dicu,dare,det1,dp, &
     &dpdav,dpdav2,dphi,dps1,dps11,dummy,ed1,ed2,g1,g2,g3,ox,oxp,oxp1,  &
     &oz,ozp,ozp1,phi,r0,r2b,r2bf,rb,rbf,rdd,rho2b,rkb,rkbf,rrad,       &
     &scikveb,scrkveb,sfac1,sfac2,sfac2s,sfac3,sfac4,sfac5,sigm1,       &
     &sigmdac,startco,sx,tas,tkb,tl,x2pi,xbb,xrb,xs,zbb,zfeld1,zfeld2,  &
     &zrb,zs,  crabfreq, crabpht, crabpht2, crabpht3, crabpht4
      character*16 typ
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,    &
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,    &
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele=1200,nblo=600,nper=16,nelb=140,nblz=20000,         &
     &nzfz = 300000,mmul = 20)
      parameter(nran = 2000000,ncom = 100,mran = 500,mpa = 6,nrco = 5,  &
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 99)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 350)
      double precision c180e0,c1e1,c1e12,c1e13,c1e15,c1e16,c1e2,c1e3,   &
     &c1e4,c1e6,c1m1,c1m7,c1m10,c1m12,c1m13,c1m15,c1m18,c1m2,c1m21,     &
     &c1m24,c1m3,c1m36,c1m38,c1m6,c1m9,c2e3,c4e3,crade,clight,four,half,&
     &one,pieni,pmae,pmap,three,two,zero
      parameter(pieni = 1d-38)
      parameter(zero = 0.0d0,half = 0.5d0,one = 1.0d0)
      parameter(two = 2.0d0,three = 3.0d0,four = 4.0d0)
      parameter(c1e1 = 1.0d1,c1e2 = 1.0d2,c1m2 = 1.0d-2)
      parameter(c1e3 = 1.0d3,c2e3 = 2.0d3,c4e3 = 4.0d3,c1e4 = 1.0d4)
      parameter(c1e12 = 1.0d12,c1e13 = 1.0d13,c1e15 = 1.0d15,c1e16 =    &
     &1.0d16)
      parameter(c180e0 = 180.0d0,c1e6 = 1.0d6)
      parameter(c1m1 = 1.0d-1,c1m3 = 1.0d-3,c1m6 = 1.0d-6,c1m7 = 1.0d-7)
      parameter(c1m9 = 1.0d-9,c1m10 = 1.0d-10,c1m12 = 1.0d-12)
      parameter(c1m13 = 1.0d-13,c1m15 = 1.0d-15)
      parameter(c1m18 = 1.0d-18,c1m21 = 1.0d-21,c1m24 = 1.0d-24)
      parameter(c1m36 = 1.0d-36,c1m38 = 1.0d-38)
!     electron mass from PDG, 2002
      parameter(pmap = 938.271998d0,pmae = .510998902d0)
      parameter(crade = 2.817940285d-15, clight = 2.99792458d8)
!-----------------------------------------------------------------------
!  COMMON FOR EXACT VERSION
!-----------------------------------------------------------------------
      integer iexact
      common/exact/iexact
      integer iav,ibb6d,ibbc,ibeco,ibidu,ibtyp,ic,icext,icextal,iclo,   &
     &iclo6,iclo6r,icode,icoe,icomb,icomb0,iconv,icow,icr,idam,idfor,   &
     &idis,idp,ierro,iffw,ifh,iicav,il,ilin,imad,imbb,                  &
     &imc,imtr,iorg,iout,                                               &
     &ipos,ipr,iprint,ipt,iq,iqmod,iqmod6,iratioe,ird,ire,ires,         &
     &irew,irm,irmod2,ise,ise1,ise2,ise3,isea,iskew,iskip,istw,         &
     &isub,itco,itcro,itf,ithick,ition,itionc,itqv,its6d,iu,iver,ivox,  &
     &ivoz,iwg,ixcav,izu0,kanf,kp,kpa,kwtype,kz,lhc,m21,m22,m23,mblo,   &
     &mbloz,mcut,mel,mesa,mmac,mout2,mp,mper,mstr,msym,mtyp,mzu,napx,   &
     &napxo,nbeam,nch,ncororb,ncorrep,ncorru,ncy,ndafi,nde,nhcorr,      &
     &nhmoni,niu,nlin,nmu,npp,nprint,nqc,nre,nrr,nskew,                 &
     &nstart,nstop,nt,nta,ntco,nte,ntwin,nu,numl,numlr,nur,nvcorr,      &
     &nvmoni,nwr, nturn1, nturn2, nturn3, nturn4,numlcp,numlmax,nnuml
      double precision a,ak0,aka,alfx,alfz,amp0,aper,apx,apz,ape,bbcu,  &
     &bclorb,beamoff,benkc,benki,betac,betam,betx,betz,bk0,bka,bl1,bl2, &
     &clo6,clobeam,clop6,cma1,cma2,cotr,crad,de0,dech,ded,dfft,         &
     &di0,dip0,dki,dkq,dma,dmap,dphix,dphiz,dppoff,dpscor,dqq,dres,dsi, &
     &dsm0,dtr,e0,ed,ej,ejf,ek,el,elbe,emitx,emity,emitz,extalign,      &
     &exterr,eui,euii,gammar,hsy,hsyc,pac,pam,parbe,parbe14,partnum,    &
     &phas,phas0,phasc,pi,pi2,pisqrt,pma,ptnfac,qs,qw0,qwsk,qx0,qxt,qz0,&
     &qzt,r00,rad,rat,ratio,ratioe,rrtr,rtc,rts,rvf,                    &
     &sigcor,sige,sigma0,sigman,sigman2,sigmanq,sigmoff,sigz,sm,ta,tam1,&
     &tam2,tiltc,tilts,tlen,totl,track6d,xpl,xrms,zfz,zpl,zrms,wirel,   &
     &acdipph, crabph, bbbx, bbby, bbbs,                                &
     &crabph2, crabph3, crabph4
      real hmal
      character*16 bez,bezb,bezr,erbez,bezl
      character*80 toptit,sixtit,commen
      common/erro/ierro,erbez
      common/kons/pi,pi2,pisqrt,rad
      common/str /il,mper,mblo,mbloz,msym(nper),kanf,iu,ic(nblz)
      common/ell /ed(nele),el(nele),ek(nele),sm(nele),kz(nele),kp(nele)
      common/bbb /bbbx(nele), bbby(nele), bbbs(nele)
      common/pla /xpl(nele),xrms(nele),zpl(nele),zrms(nele)
      common/str2 /mel(nblo),mtyp(nblo,nelb),mstr(nblo)
      common/mat/a(nele,2,6),bl1(nblo,2,6),bl2(nblo,2,6)
      common/syos2/rvf(mpa)
      common/tra1/rat,idfor,napx,napxo,numl,niu(2),numlr,nde(2),nwr(4), &
     &ird,imc,irew,ntwin,iclo6,iclo6r,iver,ibidu,numlcp,numlmax,nnuml
      common/syn/qs,e0,pma,ej(mpa),ejf(mpa),phas0,phas,hsy(3),          &
     &crad,hsyc(nele),phasc(nele),dppoff,sigmoff(nblz),tlen,            &
     &iicav,itionc(nele),ition,idp,ncy,ixcav
      common/corcom/dpscor,sigcor,icode,idam,its6d
      common/multi/bk0(nele,mmul),ak0(nele,mmul),                       &
     &bka(nele,mmul),aka(nele,mmul)
      common/mult1/benki,benkc(nele),r00(nele),irm(nele),nmu(nele)
      common/rand0/zfz(nzfz),iorg,mzu(nblz),bezr(3,nele),izu0,mmac,mcut
      common/rand1/exterr(nblz,40),extalign(nblz,3),tiltc(nblz),        &
     &tilts(nblz),mout2,icext(nblz),icextal(nblz)
      common/beo /aper(2),di0(2),dip0(2),ta(6,6)
      common/clo/dma,dmap,dkq,dqq,de0,ded,dsi,dech,dsm0,itco,itcro,itqv,&
     &iout
      common/qmodi/qw0(3),amp0,iq(3),iqmod,kpa(nele),iqmod6
      common/linop/bez(nele),elbe(nblo),bezb(nblo),ilin,nt,iprint,      &
     &ntco,eui,euii,nlin,bezl(nele)
      common/cororb/betam(nmon1,2),pam(nmon1,2),betac(ncor1,2),         &
     &pac(ncor1,2),bclorb(nmon1,2),nhmoni,nhcorr,nvmoni,nvcorr,         &
     &ncororb(nele)
      common/apert/apx(nele),apz(nele),ape(3,nele)
      common/clos/sigma0(2),iclo,ncorru,ncorrep
      common/combin/icomb0(20),icomb(ncom,20),ratio(ncom,20),           &
     &ratioe(nele),iratioe(nele),icoe
      common/seacom/ise,mesa,mp,m21,m22,m23,ise1,ise2,ise3,isea(nele)
      common/subres/qxt,qzt,tam1,tam2,isub,nta,nte,ipt,totl
      common/secom/rtc(9,18,10,5),rts(9,18,10,5),ire(12),ipr(5),irmod2
      common/secom1/dtr(10),nre,nur,nch,nqc,npp,nrr(5),nu(5)
      common/postr/dphix,dphiz,qx0,qz0,dres,dfft,cma1,cma2,             &
     &nstart,nstop,iskip,iconv,imad
      common/posti1/ipos,iav,iwg,ivox,ivoz,ires,ifh,toptit(5)
      common/posti2/kwtype,itf,icr,idis,icow,istw,iffw,nprint,ndafi
      common/skew/qwsk(2),betx(2),betz(2),alfx(2),alfz(2),iskew,nskew(6)
      common/pawc/hmal(nplo)
      common/tit/sixtit,commen,ithick
      common/co6d/clo6(3),clop6(3)
      common/dkic/dki(nele,3)
      common/beam/sigman(2,nbb),sigman2(2,nbb),sigmanq(2,nbb),          &
     &clobeam(6,nbb),beamoff(6,nbb),parbe(nele,5),track6d(6,npart),     &
     &ptnfac(nele),sigz,sige,partnum,parbe14,emitx,emity,emitz,gammar,  &
     &nbeam,ibbc,ibeco,ibtyp,lhc
      common/trom/ cotr(ntr,6),rrtr(ntr,6,6),imtr(nele)
      common/bb6d/ bbcu(nbb,12),ibb6d,imbb(nblz)
      common/wireco/ wirel(nele)
      common/acdipco/ acdipph(nele), nturn1(nele), nturn2(nele),        &
     &nturn3(nele), nturn4(nele)
      common/crabco/ crabph(nele),crabph2(nele),                        &
     &crabph3(nele),crabph4(nele)
      integer idz,itra
      double precision a2,al,as,at,chi0,chid,dp1,dps,exz,sigm
      common/syos/as(6,2,2,nele),at(6,2,2,nele),a2(6,2,2,nele),         &
     &al(6,2,2,nele),sigm(mpa),dps(mpa),idz(2)
      common/anf/chi0,chid,exz(2,6),dp1,itra
      integer ichrom,issss
      double precision alf0,amp,bet0,clo,clop,cro,xxtr,yytr
      common/tra/xxtr(mpa,2),yytr(mpa,2),amp(2),                        &
     &bet0(2),alf0(2),clo(2),clop(2)
      common/chrom/cro(2),issss(2),ichrom
      integer icorr,idial,idptr,imod1,imod2,inorm,ipar,namp,ncor,nctype,&
     &ndimf,nmom,nmom1,nmom2,nord,nord1,nordf,nsix,nvar,nvar2,nvarf
      double precision dpmax,preda,weig1,weig2
      character*16 coel
      common/dial/preda,idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
      common/norf/nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
      common/tcorr/icorr,nctype,namp,nmom,nmom1,nmom2,weig1,weig2,dpmax,&
     &coel(10)
      double precision aml6,edcor
      common/sixdim/aml6(6,6),edcor(2)
      double precision aai,ampt,bbi,damp,smi,smizf,xsi,                 &
     &zsi
      integer napxto
      real tlim,time0,time1,time2,time3,trtime
! fixes for CPU time (for all versions, not just crlibm).
      real pretime,posttime,tottime
      common/xz/xsi(nblz),zsi(nblz),smi(nblz),smizf(nblz),              &
     &aai(nblz,mmul),bbi(nblz,mmul)
      common/damp/damp,ampt
      common/ttime/tlim,time0,time1,time2,time3,trtime,napxto,          &
     &pretime,posttime,tottime
      integer numx
      double precision e0f
      common/main4/ e0f,numx
      common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda, &
     &ej1,ejf1,rv
      double precision ald6,asd6
      common/dael6/ald6(nele,2,6,nema),asd6(nele,2,6,nema)
      integer ichromc,ilinc,iqmodc
      double precision clon,chromc,corr,wxys
      common/correct/ corr(3,3),chromc(2),wxys(3),clon(6),iqmodc,       &
     &ichromc,ilinc
      double precision tasm
      common/tasm/tasm(6,6)
      dimension damap(6),damapi(6),damap1(6)
      dimension aa2(6),aa2r(6),a1(6),a1r(6),xy(6),df(6)
      dimension zfeld1(100),zfeld2(100)
      dimension jj(100),dpdav2(6),rrad(3),rdd(6,6),dicu(20)
      dimension angnoe(3),angp(2,6),phi(3),dphi(3)
      dimension b1(3),b2(3),b3(3),al1(3),al2(3),al3(3),g1(3),g2(3),g3(3)
      dimension d(3),dp(3),c(3),cp(3),au(6,6),aui(2)
      dimension i4(10,2)
      save
!-----------------------------------------------------------------------
*FOX  B D ;
*FOX  D V DA COM SIGMDA NORD NVAR ; D V DA COM DPDA NORD NVAR ;
*FOX  D V DA COM DPDA1 NORD NVAR ; D V DA COM RV NORD NVAR ;
*FOX  D V DA COM XX NORD NVAR 2 ; D V DA COM YY NORD NVAR 2 ;
*FOX  D V DA COM EJ1 NORD NVAR ; D V DA COM EJF1 NORD NVAR ;
*FOX  D V DA COM ALDA NORD NVAR 2 6 ; D V DA COM ASDA NORD NVAR 2 6 ;
*FOX  D V DA COM ALDAQ NORD NVAR 2 6 ; D V DA COM ASDAQ NORD NVAR 2 6 ;
*FOX  D V DA COM SMIDA NORD NVAR MCOR ;
*FOX  D V DA INT X NORD NVAR 2 ; D V DA INT Y NORD NVAR 2 ;
*FOX  D V DA INT YP NORD NVAR 2 ; D V DA INT DKIP NORD NVAR ;
*FOX  D V DA INT CORROLD NORD NVAR MCOP ;
*FOX  D V DA INT CORRNEW NORD NVAR MCOP ;
*FOX  D V DA INT CORRAU1 NORD NVAR MCOP ;
*FOX  D V DA INT CORRAU2 NORD NVAR MCOP ;
*FOX  D V DA INT AA NORD NVAR 11 ;  D V DA INT BB NORD NVAR 11 ;
*FOX  D V DA INT TRACKI NORD NVAR 6 ;
*FOX  D V DA INT PUX NORD NVAR ; D V DA INT PUZ NORD NVAR ;
*FOX  D V DA INT EJF0 NORD NVAR ; D V DA INT EKK NORD NVAR ;
*FOX  D V DA INT XL NORD NVAR ; D V DA INT ZL NORD NVAR ;
*FOX  D V DA INT CRKVE NORD NVAR ; D V DA INT CIKVE NORD NVAR ;
*FOX  D V DA INT CRKVEUK NORD NVAR ; D V DA INT CBZBF NORD NVAR ;
*FOX  D V DA INT YV1J NORD NVAR ; D V DA INT YV2J NORD NVAR ;
*FOX  D V DA INT CRKVEBF NORD NVAR ; D V DA INT CIKVEBF NORD NVAR ;
*FOX  D V DA INT RHO2BF NORD NVAR ; D V DA INT TKBF NORD NVAR ;
*FOX  D V DA INT XRBF NORD NVAR ; D V DA INT CCCC NORD NVAR ;
*FOX  D V DA INT ZRBF NORD NVAR ; D V DA INT XBBF NORD NVAR ;
*FOX  D V DA INT ZBBF NORD NVAR ; D V DA INT CRXBF NORD NVAR ;
*FOX  D V DA INT CBXBF NORD NVAR ; D V DA INT CRZBF NORD NVAR ;
*FOX  D V DA INT WX NORD NVAR ; D V DA INT WY NORD NVAR ;
*FOX  D V DA INT CRABAMP NORD NVAR ;
*FOX  D V DA INT CRABAMP2 NORD NVAR ;
*FOX  D V DA INT CRABAMP3 NORD NVAR ;
*FOX  D V DA INT CRABAMP4 NORD NVAR ;
*FOX  D V RE INT AAI NBLZ MMUL ; D V RE INT BBI NBLZ MMUL ;
*FOX  D V RE INT TILTC NBLZ ; D V RE INT TILTS NBLZ ;
*FOX  D V RE INT DPS MPA ; D V RE INT SIGM MPA ;
*FOX  D V RE INT DKI NELE 3 ; D V RE INT BL1 NBLO 2 6 ;
*FOX  D V RE INT EL NELE ; D V RE INT EJ MPA ; D V RE INT EJF MPA ;
*FOX  D V RE INT SMI NBLZ ; D V RE INT SMIZF NBLZ ;
*FOX  D V RE INT ED1 ; D V RE INT ED2 ;
*FOX  D V RE INT DPDAV2 6 ; D V RE INT RRTR NTR 6 6 ;
*FOX  D V RE INT COTR NTR 6 ; D V RE INT DPDAV ; D V RE INT BETR0 ;
*FOX  D V RE INT E0 ; D V RE INT E0F ; D V RE INT PMA ;
*FOX  D V RE INT XS ; D V RE INT ZS ; D V RE INT OX ; D V RE INT OXP ;
*FOX  D V RE INT OZ ; D V RE INT OZP ; D V RE INT SIGM1 ;
*FOX  D V RE INT BEAMOFF1 ; D V RE INT BEAMOFF2 ;
*FOX  D V RE INT BEAMOFF4 ; D V RE INT BEAMOFF5 ; D V RE INT BEAMOFF6 ;
*FOX  D V RE INT DPS1 ; D V RE INT RKBF ; D V RE INT RBF ;
*FOX  D V RE INT R2BF ; D V RE INT BBCU NBB 12 ;
*FOX  D V RE INT SIGMAN 2 NBB ; D V RE INT PTNFAC NELE ;
*FOX  D V RE INT CRAD ; D V RE INT GAMMAR ;
*FOX  D V RE INT PARTNUM ; D V RE INT PISQRT ; D V RE INT SCRKVEB ;
*FOX  D V RE INT SCIKVEB ; D V RE INT STARTCO ; D V RE INT RATIOE NELE ;
*FOX  D V RE INT PARBE14 ; D V RE INT PI ;
*FOX  D V RE INT SIGMDAC ; D V RE INT DUMMY ;
*FOX  D V RE INT ED NELE ; D V RE INT EK NELE ;
*FOX  D V RE INT C5M4 ; D V RE INT C2E3 ; D V RE INT C1E6 ;
*FOX  D V RE INT C1E3 ; D V RE INT C1M3 ; D V RE INT C1M6 ;
*FOX  D V RE INT C1M9 ; D V RE INT C1M12 ; D V RE INT C1M15 ;
*FOX  D V RE INT C1M18 ; D V RE INT C1M21 ; D V RE INT C1M24 ;
*FOX  D V RE INT ONE ; D V RE INT TWO ; D V RE INT THREE ;
*FOX  D V RE INT FOUR ; D V RE INT ZERO ; D V RE INT HALF ;
*FOX  D V RE INT CRABFREQ ; D V RE INT CRABPHT ;
*FOX  D V RE INT CRABPHT2 ; D V RE INT CRABPHT3 ;
*FOX  D V RE INT CRABPHT4 ;
*FOX  D V RE INT CLIGHT ;
*FOX  D V IN INT IDZ 2 ; D V IN INT KX ; D V IN INT IX ; D V IN INT JX ;
*FOX  D V IN INT I ; D V IN INT IPCH ; D V IN INT K ; D V IN INT KKK ;
*FOX  D V IN INT IVAR ; D V IN INT IRRTR ; D V IN INT KK ;
*FOX  D V IN INT IMBB NBLZ ;
*FOX  D F RE DARE 1 ;
*FOX  D V DA INT PZ NORD NVAR ;
*FOX  E D ;
*FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
      nd2=ndimf*2
      call etall(damap,6)
      call etall(damapi,6)
      call etall(damap1,6)
      call etall(angno,1)
      call etall(f,1)
      call etall(aa2,6)
      call etall(aa2r,6)
      call etall(a1,6)
      call etall(a1r,6)
      call etall(xy,6)
      call etall(h,1)
      call etall(df,6)
      if(iqmodc.eq.1) call danot(2)
      if(iqmodc.eq.3) call danot(1)
      if(ichromc.eq.1) call danot(3)
      icoonly=0
      if(iqmodc.eq.2.or.iqmodc.eq.4.or.ichromc.eq.2) icoonly=1
      c5m4=5.0d-4
      do j=1,2
        angnoe(j)=zero
        do i=1,6
          angp(j,i)=zero
        enddo
      enddo
      do i=1,100
        jj(i)=0
      enddo
      x2pi=atan_rn(one)*8d0
      i4(1,1)=1
      i4(1,2)=1
      i4(2,1)=3
      i4(2,2)=3
      i4(3,1)=1
      i4(3,2)=3
      i4(4,1)=1
      i4(4,2)=2
      i4(5,1)=1
      i4(5,2)=4
      i4(6,1)=2
      i4(6,2)=2
      i4(7,1)=2
      i4(7,2)=3
      i4(8,1)=2
      i4(8,2)=4
      i4(9,1)=3
      i4(9,2)=4
      i4(10,1)=4
      i4(10,2)=4
!hr05 e0f=sqrt(e0*e0-pma*pma)
      e0f=sqrt(e0**2-pma**2)                                             !hr05
      betr0=sqrt(one-(pma/e0)**2)
      ox=xxtr(1,1)
      oxp=yytr(1,1)
      oz=xxtr(1,2)
      ozp=yytr(1,2)
      sigm1=sigm(1)
      dps1=dps(1)
      if(iqmodc.eq.1) then
        if(el(iq(1)).le.pieni) then
          ed1=ed(iq(1))
        else
          ed1=ek(iq(1))
        endif
        if(el(iq(2)).le.pieni) then
          ed2=ed(iq(2))
        else
          ed2=ek(iq(2))
        endif
      endif
      if(ichromc.eq.1) then
        ed1=ed(issss(1))
        ed2=ed(issss(2))
      endif
      call davar(x(1),ox,1)
      oxp1=oxp*(one+dps1)
      call davar(yp(1),oxp1,2)
      ivar=2
      if(nvar2.ge.4) then
        call davar(x(2),oz,3)
        ozp1=ozp*(one+dps1)
        call davar(yp(2),ozp1,4)
        ivar=4
      else
*FOX  X(2)=OZ ;
*FOX  YP(2)=OZP*(ONE+DPS1) ;
      endif
      dps11=dps1*c1e3
      if(nvar2.eq.3) then
        call davar(dpda1,dps11,3)
        ivar=ivar+1
      elseif(nvar2.eq.5) then
        call davar(dpda1,dps11,5)
        ivar=ivar+1
      elseif(nvar2.eq.6) then
        call davar(sigmda,sigm1,5)
        call davar(dpda1,dps11,6)
        ivar=ivar+2
      else
*FOX  SIGMDA=SIGM1 ;
*FOX  DPDA1=DPS1*C1E3 ;
      endif
      ivar1=ivar
      if(iqmodc.eq.1.or.ichromc.eq.1) then
        call davar(smida(1),ed1,ivar+1)
        call davar(smida(2),ed2,ivar+2)
        ivar=ivar+2
      endif
!--Normal Form Analysis for calculation of linear lattice functions
      if(ilinc.eq.1.or.ilinc.eq.2) then
        mfile=18
!Eric
        rewind mfile
        rewind 111
!ERIC HERE
        call daread(damap,nvar,mfile,1.d0)
        call mapnorm(damap,f,aa2,a1,xy,h,nord1)
        do j=1,nvar
          call dacop(damap(j),damap1(j))
          dummy=dare(damap1(j))
          call dacsu(damap1(j),dummy,damap1(j))
        enddo
        if(ndimf.eq.3) then
          call damul(damap1(5),damap1(5),angno)
          call averaged(angno,damap1,.true.,angno)
          jj(5)=1
          jj(6)=1
          call dapek(angno,jj,emitz)
          jj(5)=0
          jj(6)=0
          if(abs(emitz).le.pieni) then
            emitz=zero
          else
!hr05       emitz=sigz*sigz/emitz*half*c1e6
            emitz=((sigz**2/emitz)*half)*c1e6                            !hr05
          endif
        endif
        jj(5)=1
        do j=1,nd2
          call dapek(a1(j),jj,dicu(j))
        enddo
        jj(5)=0
      endif
*FOX  DPDA=DPDA1*C1M3 ;
*FOX  CORROLD(1)=X(1) ;
*FOX  CORROLD(2)=YP(1) ;
*FOX  CORROLD(3)=X(2) ;
*FOX  CORROLD(4)=YP(2) ;
*FOX  CORROLD(5)=SIGMDA ;
*FOX  CORROLD(6)=DPDA1 ;
            do 5 kkk=1,6
              dpdav=dare(corrold(kkk))
*FOX  CORROLD(KKK)=CORROLD(KKK)-DPDAV ;
    5       continue
*FOX  Y(1)=YP(1)/(ONE+DPDA) ;
*FOX  Y(2)=YP(2)/(ONE+DPDA) ;
      iflag=0
      iflag1=0
      iflag2=0
*FOX  EJF1=E0F*(ONE+DPDA) ;
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;
*FOX  RV=EJ1/E0*E0F/EJF1 ;
      if(ithick.eq.1) call envada
      icav=0
      typ='START'
      phi(1)=zero
      phi(2)=zero
      phi(3)=zero
      ibb=0
      do 430 i=1,iu
        if(iqmodc.eq.2.or.iqmodc.eq.4) then
          if(i.eq.niu(1)) then
            do ii=1,2
              ii2=2*ii
              clon(ii2-1)=dare(x(ii))
              clon(ii2)=dare(y(ii))
            enddo
            clon(5)=dare(sigmda)
            clon(6)=dare(dpda)
          endif
        endif
        if(ilinc.eq.1.and.i.eq.1) then
          write(*,10000) nd2
          if(iprint.eq.1) write(*,10130)
          write(*,10010)
          write(*,10020)
          write(*,10010)
          tl=zero
          iwrite=0
          if(nlin.eq.0) then
            iwrite=1
          else
            do ii=1,nlin
              if(typ.eq.bezl(ii)) iwrite=1
            enddo
          endif
*FOX  YP(1)=Y(1)*(ONE+DPDA) ;
*FOX  YP(2)=Y(2)*(ONE+DPDA) ;
*FOX  DPDA1=DPDA*C1E3 ;
          call dacop(x(1),damap(1))
          call dacop(yp(1),damap(2))
          call dacop(x(2),damap(3))
          call dacop(yp(2),damap(4))
          do j=1,2
            ii=2*j
            call dapek(damap(ii-1),jj,c(j))
            call dapek(damap(ii),jj,cp(j))
          enddo
          call dacsu(damap(1),c(1),damap(1))
          call dacsu(damap(2),cp(1),damap(2))
          call dacsu(damap(3),c(2),damap(3))
          call dacsu(damap(4),cp(2),damap(4))
          if(ndimf.eq.3) then
            call dacop(sigmda,damap(5))
            call dacop(dpda1,damap(6))
            call dapek(damap(5),jj,c(3))
            call dapek(damap(6),jj,cp(3))
            call dacsu(damap(5),c(3),damap(5))
            call dacsu(damap(6),cp(3),damap(6))
            if(iflag1.eq.1.and.ithick.eq.1) then
              call dacct(damap,nvar,corrnew,nvar,damap,nvar)
            endif
          else
            call dacop(dpda1,damap(5))
            do j1=1,4
              do ii=1,4
                jj(ii)=1
                call dapek(damap(j1),jj,rdd(j1,ii))
                jj(ii)=0
              enddo
            enddo
            jj(5)=1
            do j1=1,4
              call dapek(damap(j1),jj,rdd(j1,5))
            enddo
            jj(5)=0
            do j1=1,2
              ii=2*j1
!hr03         d(j1)=rdd(ii-1,1)*dicu(1)+rdd(ii-1,2)*dicu(2)+            &
!hr03&rdd(ii-1,3)*dicu(3)+rdd(ii-1,4)*dicu(4)+rdd(ii-1,5)
              d(j1)=(((rdd(ii-1,1)*dicu(1)+rdd(ii-1,2)*dicu(2))+        &!hr03
     &rdd(ii-1,3)*dicu(3))+rdd(ii-1,4)*dicu(4))+rdd(ii-1,5)              !hr03
!hr03         dp(j1)=rdd(ii,1)*dicu(1)+rdd(ii,2)*dicu(2)+               &
!hr03&rdd(ii,3)*dicu(3)+rdd(ii,4)*dicu(4)+rdd(ii,5)
              dp(j1)=(((rdd(ii,1)*dicu(1)+rdd(ii,2)*dicu(2))+           &!hr03
     &rdd(ii,3)*dicu(3))+rdd(ii,4)*dicu(4))+rdd(ii,5)                    !hr03
            enddo
          endif
          call dacct(damap,nvar,aa2,nvar,damap,nvar)
          do j=1,ndimf
            ii=2*j
            if(j.eq.1) then
              i2=4
              i3=6
            elseif(j.eq.2) then
              i2=2
              i3=6
            elseif(j.eq.3) then
              i2=2
              i3=4
            endif
            jj(ii-1)=1
            call dapek(damap(ii-1),jj,angp(1,ii-1))
            call dapek(damap(ii),jj,au(ii,ii-1))
            jj(ii-1)=0
            jj(ii)=1
            call dapek(damap(ii-1),jj,angp(1,ii))
            call dapek(damap(ii),jj,au(ii,ii))
            jj(ii)=0
            jj(i2-1)=1
            call dapek(damap(ii-1),jj,au(i2-1,i2-1))
            call dapek(damap(ii),jj,au(i2,i2-1))
            jj(i2-1)=0
            jj(i2)=1
            call dapek(damap(ii-1),jj,au(i2-1,i2))
            call dapek(damap(ii),jj,au(i2,i2))
            jj(i2)=0
            jj(i3-1)=1
            call dapek(damap(ii-1),jj,au(i3-1,i3-1))
            call dapek(damap(ii),jj,au(i3,i3-1))
            jj(i3-1)=0
            jj(i3)=1
            call dapek(damap(ii-1),jj,au(i3-1,i3))
            call dapek(damap(ii),jj,au(i3,i3))
            jj(i3)=0
!hr08       b1(j)=angp(1,ii-1)*angp(1,ii-1)+angp(1,ii)*angp(1,ii)
            b1(j)=angp(1,ii-1)**2+angp(1,ii)**2                          !hr08
!hr08       b2(j)=au(i2-1,i2-1)*au(i2-1,i2-1)+au(i2-1,i2)*au(i2-1,i2)
            b2(j)=au(i2-1,i2-1)**2+au(i2-1,i2)**2                        !hr08
!hr08       b3(j)=au(i3-1,i3-1)*au(i3-1,i3-1)+au(i3-1,i3)*au(i3-1,i3)
            b3(j)=au(i3-1,i3-1)**2+au(i3-1,i3)**2                        !hr08
!hr03       al1(j)=-(angp(1,ii-1)*au(ii,ii-1)+angp(1,ii)*au(ii,ii))
            al1(j)=-1d0*(angp(1,ii-1)*au(ii,ii-1)+angp(1,ii)*au(ii,ii))  !hr03
!hr03       al2(j)=-(au(i2-1,i2-1)*au(i2,i2-1)+au(i2-1,i2)*au(i2,i2))
           al2(j)=-1d0*(au(i2-1,i2-1)*au(i2,i2-1)+au(i2-1,i2)*au(i2,i2)) !hr03
!hr03       al3(j)=-(au(i3-1,i3-1)*au(i3,i3-1)+au(i3-1,i3)*au(i3,i3))
           al3(j)=-1d0*(au(i3-1,i3-1)*au(i3,i3-1)+au(i3-1,i3)*au(i3,i3)) !hr03
!hr04       g1(j)=au(ii,ii-1)*au(ii,ii-1)+au(ii,ii)*au(ii,ii)
!hr04       g2(j)=au(i2,i2-1)*au(i2,i2-1)+au(i2,i2)*au(i2,i2)
!hr04       g3(j)=au(i3,i3-1)*au(i3,i3-1)+au(i3,i3)*au(i3,i3)
            g1(j)=au(ii,ii-1)**2+au(ii,ii)**2                            !hr04
            g2(j)=au(i2,i2-1)**2+au(i2,i2)**2                            !hr04
            g3(j)=au(i3,i3-1)**2+au(i3,i3)**2                            !hr04
            if(ndimf.eq.3) then
              call dainv(damap,nvar,damapi,nvar)
              jj(6)=1
              call dapek(damapi(5),jj,aui(1))
              call dapek(damapi(6),jj,aui(2))
              jj(6)=0
              if(j.lt.3) then
                d(j)=au(i3-1,i3-1)*aui(1)+au(i3-1,i3)*aui(2)
                dp(j)=au(i3,i3-1)*aui(1)+au(i3,i3)*aui(2)
              else
                d(j)=angp(1,ii-1)*aui(1)+angp(1,ii)*aui(2)
                dp(j)=au(ii,ii-1)*aui(1)+au(ii,ii)*aui(2)
              endif
            endif
            sx=angp(2,ii-1)*angp(1,ii)-angp(1,ii-1)*angp(2,ii)
            cx=angp(1,ii-1)*angp(2,ii-1)+angp(1,ii)*angp(2,ii)
            if(abs(sx).gt.c1m15.or.abs(cx).gt.c1m15) then
              dphi(j)=atan2_rn(sx,cx)/x2pi
            else
              dphi(j)=zero
            endif
            phi(j)=phi(j)+dphi(j)
          enddo
          do j=1,ndimf
            ii=2*j
            angp(2,ii-1)=angp(1,ii-1)
            angp(2,ii)=angp(1,ii)
          enddo
          if(iwrite.eq.1) then
            iii=i
            if(typ(:8).eq.'START   ') iii=0
            write(*,10030) iii,typ(:8),tl,phi(1),b1(1),al1(1),g1(1),    &
     &d(1),dp(1),c(1),cp(1)
            if(ndimf.eq.3) then
              write(*,10040) b2(1),al2(1),g2(1)
              write(*,10050) typ(9:16),b3(1),al3(1),g3(1)
            else
              write(*,10055) typ(9:16),b2(1),al2(1),g2(1)
            endif
            write(*,10060)
            write(*,10070) phi(2),b1(2),al1(2),g1(2),d(2),dp(2),c(2),   &
     &cp(2)
            write(*,10080) b2(2),al2(2),g2(2)
            if(ndimf.eq.3) then
              write(*,10090) b3(2),al3(2),g3(2)
              write(*,10060)
              write(*,10100) -phi(3),b1(3),al1(3),g1(3),d(3),dp(3),c(3),&
     &cp(3)*c1m3
              write(*,10080) b2(3),al2(3),g2(3)
              write(*,10040) b3(3),al3(3),g3(3)
            endif
            write(*,10010)
          endif
        endif
        if(iflag.eq.1) then
*FOX  EJF1=E0F*(ONE+DPDA) ;
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;
*FOX  RV=EJ1/E0*E0F/EJF1 ;
          if(ithick.eq.1) then
*FOX  YP(1)=Y(1)*(ONE+DPDA) ;
*FOX  YP(2)=Y(2)*(ONE+DPDA) ;
            if(icav.eq.0) then
*FOX  CORRNEW(1)=X(1) ;
*FOX  CORRNEW(2)=YP(1) ;
*FOX  CORRNEW(3)=X(2) ;
*FOX  CORRNEW(4)=YP(2) ;
*FOX  CORRNEW(5)=SIGMDA ;
*FOX  CORRNEW(6)=DPDA1 ;
              do 24 kkk=1,6
                dpdav=dare(corrnew(kkk))
*FOX  CORRNEW(KKK)=CORRNEW(KKK)-DPDAV ;
   24         continue
            else
*FOX  CORRAU2(1)=X(1) ;
*FOX  CORRAU2(2)=YP(1) ;
*FOX  CORRAU2(3)=X(2) ;
*FOX  CORRAU2(4)=YP(2) ;
*FOX  CORRAU2(5)=SIGMDA ;
*FOX  CORRAU2(6)=DPDA1 ;
              do 25 kkk=1,6
*FOX  CORRAU1(KKK)=CORRNEW(KKK) ;
                dpdav=dare(corrau2(kkk))
*FOX  CORRAU2(KKK)=CORRAU2(KKK)-DPDAV ;
   25         continue
              if(ivar.gt.ivar1) then
*FOX  CORRAU2(7)=SMIDA(1) ;
*FOX  CORRAU2(8)=SMIDA(2) ;
                dpdav=dare(smida(1))
*FOX  CORRAU1(7)=SMIDA(1)-DPDAV ;
                dpdav=dare(smida(2))
*FOX  CORRAU1(8)=SMIDA(2)-DPDAV ;
              endif
              call dacct(corrau2,nvar,corrau1,nvar,corrnew,nvar)
            endif
            dpdav=dare(x(1))
*FOX  X(1)=CORROLD(1)+DPDAV ;
            dpdav=dare(yp(1))
*FOX  YP(1)=CORROLD(2)+DPDAV ;
            dpdav=dare(x(2))
*FOX  X(2)=CORROLD(3)+DPDAV ;
            dpdav=dare(yp(2))
*FOX  YP(2)=CORROLD(4)+DPDAV ;
            dpdav=dare(sigmda)
*FOX  SIGMDA=CORROLD(5)+DPDAV ;
            dpdav=dare(dpda1)
*FOX  DPDA1=CORROLD(6)+DPDAV ;
*FOX  DPDA=DPDA1*C1M3 ;
*FOX  Y(1)=YP(1)/(ONE+DPDA) ;
*FOX  Y(2)=YP(2)/(ONE+DPDA) ;
*FOX  EJF1=E0F*(ONE+DPDA) ;
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;
*FOX  RV=EJ1/E0*E0F/EJF1 ;
            icav=icav+1
            call envada
          endif
          iflag=0
        endif
        ix=ic(i)
        if(ix.gt.nblo) goto 50
        if(ix.le.0) then
          call prror(93)
        endif
        jmel=mel(ix)
        if(idp.eq.0.or.ition.eq.0) then
          if(ithick.eq.1) then
            do jb=1,jmel
              jx=mtyp(ix,jb)
              do ip=1,6
                do ien=1,nord+1
                  zfeld1(ien)=ald6(jx,1,ip,ien)
                enddo
                if(nvar2.eq.4) then
                  call darea6(alda(1,ip),zfeld1,4)
                else if(nvar2.eq.5) then
                  call darea6(alda(1,ip),zfeld1,5)
                endif
                do ien=1,nord+1
                  zfeld1(ien)=ald6(jx,2,ip,ien)
                enddo
                if(nvar2.eq.4) then
                  call darea6(alda(2,ip),zfeld1,4)
                else if(nvar2.eq.5) then
                  call darea6(alda(2,ip),zfeld1,5)
                endif
              enddo
              ipch=0
              if(iqmodc.eq.1.and.kz(jx).eq.2) then
                if(jx.eq.iq(1).or.iratioe(jx).eq.iq(1)) then
                  ipch=1
                else if(jx.eq.iq(2).or.iratioe(jx).eq.iq(2)) then
                  ipch=2
                endif
              endif
              if(ipch.ne.0) then
                call envquad(jx,ipch)
*FOX  PUX=X(1) ;
*FOX  PUZ=Y(1) ;
*FOX  X(1)=ALDAQ(1,1)*PUX+ALDAQ(1,2)*PUZ+ALDAQ(1,5)*IDZ(1) ;
*FOX  Y(1)=ALDAQ(1,3)*PUX+ALDAQ(1,4)*PUZ+ALDAQ(1,6)*IDZ(1) ;
*FOX  PUX=X(2) ;
*FOX  PUZ=Y(2) ;
*FOX  X(2)=ALDAQ(2,1)*PUX+ALDAQ(2,2)*PUZ+ALDAQ(2,5)*IDZ(2) ;
*FOX  Y(2)=ALDAQ(2,3)*PUX+ALDAQ(2,4)*PUZ+ALDAQ(2,6)*IDZ(2) ;
              else
*FOX  PUX=X(1) ;
*FOX  PUZ=Y(1) ;
*FOX  X(1)=ALDA(1,1)*PUX+ALDA(1,2)*PUZ+ALDA(1,5)*IDZ(1) ;
*FOX  Y(1)=ALDA(1,3)*PUX+ALDA(1,4)*PUZ+ALDA(1,6)*IDZ(1) ;
*FOX  PUX=X(2) ;
*FOX  PUZ=Y(2) ;
*FOX  X(2)=ALDA(2,1)*PUX+ALDA(2,2)*PUZ+ALDA(2,5)*IDZ(2) ;
*FOX  Y(2)=ALDA(2,3)*PUX+ALDA(2,4)*PUZ+ALDA(2,6)*IDZ(2) ;
              endif
            enddo
          else
*FOX  X(1)=X(1)+BL1(IX,1,2)*Y(1) ;
*FOX  X(2)=X(2)+BL1(IX,2,2)*Y(2) ;
          endif
          if(ilinc.eq.1) then
            do jb=1,jmel
              jx=mtyp(ix,jb)
              typ=bez(jx)
              tl=tl+el(jx)
          iwrite=0
          if(nlin.eq.0) then
            iwrite=1
          else
            do ii=1,nlin
              if(typ.eq.bezl(ii)) iwrite=1
            enddo
          endif
*FOX  YP(1)=Y(1)*(ONE+DPDA) ;
*FOX  YP(2)=Y(2)*(ONE+DPDA) ;
*FOX  DPDA1=DPDA*C1E3 ;
          call dacop(x(1),damap(1))
          call dacop(yp(1),damap(2))
          call dacop(x(2),damap(3))
          call dacop(yp(2),damap(4))
          do j=1,2
            ii=2*j
            call dapek(damap(ii-1),jj,c(j))
            call dapek(damap(ii),jj,cp(j))
          enddo
          call dacsu(damap(1),c(1),damap(1))
          call dacsu(damap(2),cp(1),damap(2))
          call dacsu(damap(3),c(2),damap(3))
          call dacsu(damap(4),cp(2),damap(4))
          if(ndimf.eq.3) then
            call dacop(sigmda,damap(5))
            call dacop(dpda1,damap(6))
            call dapek(damap(5),jj,c(3))
            call dapek(damap(6),jj,cp(3))
            call dacsu(damap(5),c(3),damap(5))
            call dacsu(damap(6),cp(3),damap(6))
            if(iflag1.eq.1.and.ithick.eq.1) then
              call dacct(damap,nvar,corrnew,nvar,damap,nvar)
            endif
          else
            call dacop(dpda1,damap(5))
            do j1=1,4
              do ii=1,4
                jj(ii)=1
                call dapek(damap(j1),jj,rdd(j1,ii))
                jj(ii)=0
              enddo
            enddo
            jj(5)=1
            do j1=1,4
              call dapek(damap(j1),jj,rdd(j1,5))
            enddo
            jj(5)=0
            do j1=1,2
              ii=2*j1
!hr03         d(j1)=rdd(ii-1,1)*dicu(1)+rdd(ii-1,2)*dicu(2)+            &
!hr03&rdd(ii-1,3)*dicu(3)+rdd(ii-1,4)*dicu(4)+rdd(ii-1,5)
              d(j1)=(((rdd(ii-1,1)*dicu(1)+rdd(ii-1,2)*dicu(2))+        &!hr03
     &rdd(ii-1,3)*dicu(3))+rdd(ii-1,4)*dicu(4))+rdd(ii-1,5)              !hr03
!hr03         dp(j1)=rdd(ii,1)*dicu(1)+rdd(ii,2)*dicu(2)+               &
!hr03&rdd(ii,3)*dicu(3)+rdd(ii,4)*dicu(4)+rdd(ii,5)
              dp(j1)=(((rdd(ii,1)*dicu(1)+rdd(ii,2)*dicu(2))+           &!hr03
     &rdd(ii,3)*dicu(3))+rdd(ii,4)*dicu(4))+rdd(ii,5)                    !hr03
            enddo
          endif
          call dacct(damap,nvar,aa2,nvar,damap,nvar)
          do j=1,ndimf
            ii=2*j
            if(j.eq.1) then
              i2=4
              i3=6
            elseif(j.eq.2) then
              i2=2
              i3=6
            elseif(j.eq.3) then
              i2=2
              i3=4
            endif
            jj(ii-1)=1
            call dapek(damap(ii-1),jj,angp(1,ii-1))
            call dapek(damap(ii),jj,au(ii,ii-1))
            jj(ii-1)=0
            jj(ii)=1
            call dapek(damap(ii-1),jj,angp(1,ii))
            call dapek(damap(ii),jj,au(ii,ii))
            jj(ii)=0
            jj(i2-1)=1
            call dapek(damap(ii-1),jj,au(i2-1,i2-1))
            call dapek(damap(ii),jj,au(i2,i2-1))
            jj(i2-1)=0
            jj(i2)=1
            call dapek(damap(ii-1),jj,au(i2-1,i2))
            call dapek(damap(ii),jj,au(i2,i2))
            jj(i2)=0
            jj(i3-1)=1
            call dapek(damap(ii-1),jj,au(i3-1,i3-1))
            call dapek(damap(ii),jj,au(i3,i3-1))
            jj(i3-1)=0
            jj(i3)=1
            call dapek(damap(ii-1),jj,au(i3-1,i3))
            call dapek(damap(ii),jj,au(i3,i3))
            jj(i3)=0
!hr08       b1(j)=angp(1,ii-1)*angp(1,ii-1)+angp(1,ii)*angp(1,ii)
            b1(j)=angp(1,ii-1)**2+angp(1,ii)**2                          !hr08
!hr08       b2(j)=au(i2-1,i2-1)*au(i2-1,i2-1)+au(i2-1,i2)*au(i2-1,i2)
            b2(j)=au(i2-1,i2-1)**2+au(i2-1,i2)**2                        !hr08
!hr08       b3(j)=au(i3-1,i3-1)*au(i3-1,i3-1)+au(i3-1,i3)*au(i3-1,i3)
            b3(j)=au(i3-1,i3-1)**2+au(i3-1,i3)**2                        !hr08
!hr03       al1(j)=-(angp(1,ii-1)*au(ii,ii-1)+angp(1,ii)*au(ii,ii))
            al1(j)=-1d0*(angp(1,ii-1)*au(ii,ii-1)+angp(1,ii)*au(ii,ii))  !hr03
!hr03       al2(j)=-(au(i2-1,i2-1)*au(i2,i2-1)+au(i2-1,i2)*au(i2,i2))
           al2(j)=-1d0*(au(i2-1,i2-1)*au(i2,i2-1)+au(i2-1,i2)*au(i2,i2)) !hr03
!hr03       al3(j)=-(au(i3-1,i3-1)*au(i3,i3-1)+au(i3-1,i3)*au(i3,i3))
           al3(j)=-1d0*(au(i3-1,i3-1)*au(i3,i3-1)+au(i3-1,i3)*au(i3,i3)) !hr03
!hr04       g1(j)=au(ii,ii-1)*au(ii,ii-1)+au(ii,ii)*au(ii,ii)
!hr04       g2(j)=au(i2,i2-1)*au(i2,i2-1)+au(i2,i2)*au(i2,i2)
!hr04       g3(j)=au(i3,i3-1)*au(i3,i3-1)+au(i3,i3)*au(i3,i3)
            g1(j)=au(ii,ii-1)**2+au(ii,ii)**2                            !hr04
            g2(j)=au(i2,i2-1)**2+au(i2,i2)**2                            !hr04
            g3(j)=au(i3,i3-1)**2+au(i3,i3)**2                            !hr04
            if(ndimf.eq.3) then
              call dainv(damap,nvar,damapi,nvar)
              jj(6)=1
              call dapek(damapi(5),jj,aui(1))
              call dapek(damapi(6),jj,aui(2))
              jj(6)=0
              if(j.lt.3) then
                d(j)=au(i3-1,i3-1)*aui(1)+au(i3-1,i3)*aui(2)
                dp(j)=au(i3,i3-1)*aui(1)+au(i3,i3)*aui(2)
              else
                d(j)=angp(1,ii-1)*aui(1)+angp(1,ii)*aui(2)
                dp(j)=au(ii,ii-1)*aui(1)+au(ii,ii)*aui(2)
              endif
            endif
            sx=angp(2,ii-1)*angp(1,ii)-angp(1,ii-1)*angp(2,ii)
            cx=angp(1,ii-1)*angp(2,ii-1)+angp(1,ii)*angp(2,ii)
            if(abs(sx).gt.c1m15.or.abs(cx).gt.c1m15) then
              dphi(j)=atan2_rn(sx,cx)/x2pi
            else
              dphi(j)=zero
            endif
            phi(j)=phi(j)+dphi(j)
          enddo
          do j=1,ndimf
            ii=2*j
            angp(2,ii-1)=angp(1,ii-1)
            angp(2,ii)=angp(1,ii)
          enddo
          if(iwrite.eq.1) then
            iii=i
            if(typ(:8).eq.'START   ') iii=0
            write(*,10030) iii,typ(:8),tl,phi(1),b1(1),al1(1),g1(1),    &
     &d(1),dp(1),c(1),cp(1)
            if(ndimf.eq.3) then
              write(*,10040) b2(1),al2(1),g2(1)
              write(*,10050) typ(9:16),b3(1),al3(1),g3(1)
            else
              write(*,10055) typ(9:16),b2(1),al2(1),g2(1)
            endif
            write(*,10060)
            write(*,10070) phi(2),b1(2),al1(2),g1(2),d(2),dp(2),c(2),   &
     &cp(2)
            write(*,10080) b2(2),al2(2),g2(2)
            if(ndimf.eq.3) then
              write(*,10090) b3(2),al3(2),g3(2)
              write(*,10060)
              write(*,10100) -phi(3),b1(3),al1(3),g1(3),d(3),dp(3),c(3),&
     &cp(3)*c1m3
              write(*,10080) b2(3),al2(3),g2(3)
              write(*,10040) b3(3),al3(3),g3(3)
            endif
            write(*,10010)
          endif
              if(i.eq.nt) goto 470
            enddo
          endif
        else
          do jb=1,jmel
            jx=mtyp(ix,jb)
            if(ithick.eq.1) then
              do ip=1,6
                do ien=1,nord+1
                  zfeld1(ien)=ald6(jx,1,ip,ien)
                  zfeld2(ien)=asd6(jx,1,ip,ien)
                enddo
                if(nvar2.eq.5) then
                  call darea6(alda(1,ip),zfeld1,5)
                  call darea6(asda(1,ip),zfeld2,5)
                else if(nvar2.eq.6) then
                  call darea6(alda(1,ip),zfeld1,6)
                  call darea6(asda(1,ip),zfeld2,6)
                endif
                do ien=1,nord+1
                  zfeld1(ien)=ald6(jx,2,ip,ien)
                  zfeld2(ien)=asd6(jx,2,ip,ien)
                enddo
                if(nvar2.eq.5) then
                  call darea6(alda(2,ip),zfeld1,5)
                  call darea6(asda(2,ip),zfeld2,5)
                else if(nvar2.eq.6) then
                  call darea6(alda(2,ip),zfeld1,6)
                  call darea6(asda(2,ip),zfeld2,6)
                endif
              enddo
              ipch=0
              if(iqmodc.eq.1.and.kz(jx).eq.2) then
                if(jx.eq.iq(1).or.iratioe(jx).eq.iq(1)) then
                  ipch=1
                else if(jx.eq.iq(2).or.iratioe(jx).eq.iq(2)) then
                  ipch=2
                endif
              endif
              if(ipch.ne.0) then
                call envquad(jx,ipch)
*FOX  PUX=X(1) ;
*FOX  PUZ=Y(1) ;
*FOX  SIGMDA=SIGMDA+ASDAQ(1,1)+ASDAQ(1,2)*PUX+
*FOX  ASDAQ(1,3)*PUZ+ASDAQ(1,4)*PUX*PUZ+ASDAQ(1,5)*PUX*PUX+
*FOX  ASDAQ(1,6)*PUZ*PUZ ;
*FOX  X(1)=ALDAQ(1,1)*PUX+ALDAQ(1,2)*PUZ+ALDAQ(1,5)*IDZ(1) ;
*FOX  Y(1)=ALDAQ(1,3)*PUX+ALDAQ(1,4)*PUZ+ALDAQ(1,6)*IDZ(1) ;
*FOX  PUX=X(2) ;
*FOX  PUZ=Y(2) ;
*FOX  SIGMDA=SIGMDA+ASDAQ(2,1)+ASDAQ(2,2)*PUX+
*FOX  ASDAQ(2,3)*PUZ+ASDAQ(2,4)*PUX*PUZ+ASDAQ(2,5)*PUX*PUX+
*FOX  ASDAQ(2,6)*PUZ*PUZ ;
*FOX  X(2)=ALDAQ(2,1)*PUX+ALDAQ(2,2)*PUZ+ALDAQ(2,5)*IDZ(2) ;
*FOX  Y(2)=ALDAQ(2,3)*PUX+ALDAQ(2,4)*PUZ+ALDAQ(2,6)*IDZ(2) ;
              else
*FOX  PUX=X(1) ;
*FOX  PUZ=Y(1) ;
*FOX  SIGMDA=SIGMDA+ASDA(1,1)+ASDA(1,2)*PUX+
*FOX  ASDA(1,3)*PUZ+ASDA(1,4)*PUX*PUZ+ASDA(1,5)*PUX*PUX+
*FOX  ASDA(1,6)*PUZ*PUZ ;
*FOX  X(1)=ALDA(1,1)*PUX+ALDA(1,2)*PUZ+ALDA(1,5)*IDZ(1) ;
*FOX  Y(1)=ALDA(1,3)*PUX+ALDA(1,4)*PUZ+ALDA(1,6)*IDZ(1) ;
*FOX  PUX=X(2) ;
*FOX  PUZ=Y(2) ;
*FOX  SIGMDA=SIGMDA+ASDA(2,1)+ASDA(2,2)*PUX+
*FOX  ASDA(2,3)*PUZ+ASDA(2,4)*PUX*PUZ+ASDA(2,5)*PUX*PUX+
*FOX  ASDA(2,6)*PUZ*PUZ ;
*FOX  X(2)=ALDA(2,1)*PUX+ALDA(2,2)*PUZ+ALDA(2,5)*IDZ(2) ;
*FOX  Y(2)=ALDA(2,3)*PUX+ALDA(2,4)*PUZ+ALDA(2,6)*IDZ(2) ;
              endif
            else
              if(iexact.eq.1) then
!-----------------------------------------------------------------------
!  EXACT DRIFT
!-----------------------------------------------------------------------
*FOX  X(1)=X(1)*C1M3 ;
*FOX  X(2)=X(2)*C1M3 ;
*FOX  Y(1)=Y(1)*C1M3 ;
*FOX  Y(2)=Y(2)*C1M3 ;
*FOX  SIGMDA=SIGMDA*C1M3 ;
*FOX  PZ=SQRT(ONE-Y(1)*Y(1)-Y(2)*Y(2)) ;
*FOX  X(1)=X(1)+EL(JX)*(Y(1)/PZ) ;
*FOX  X(2)=X(2)+EL(JX)*(Y(2)/PZ) ;
*FOX  SIGMDA=SIGMDA+(ONE-(RV/PZ))*EL(JX) ;
*FOX  X(1)=X(1)*C1E3 ;
*FOX  X(2)=X(2)*C1E3 ;
*FOX  Y(1)=Y(1)*C1E3 ;
*FOX  Y(2)=Y(2)*C1E3 ;
*FOX  SIGMDA=SIGMDA*C1E3 ;
!-----------------------------------------------------------------------
              else
! Regular drift
!            else !moved outside of dalin6 /Mattias
*FOX  X(1)=X(1)+EL(JX)*Y(1) ;
*FOX  X(2)=X(2)+EL(JX)*Y(2) ;
!
*FOX  SIGMDA=SIGMDA+
*FOX  EL(JX)*(C1E3-RV*(C1E3+(Y(1)*Y(1)+Y(2)*Y(2))*C5M4)) ;
              endif
            endif
            if(ilinc.eq.1) then
              typ=bez(jx)
              tl=tl+el(jx)
          iwrite=0
          if(nlin.eq.0) then
            iwrite=1
          else
            do ii=1,nlin
              if(typ.eq.bezl(ii)) iwrite=1
            enddo
          endif
*FOX  YP(1)=Y(1)*(ONE+DPDA) ;
*FOX  YP(2)=Y(2)*(ONE+DPDA) ;
*FOX  DPDA1=DPDA*C1E3 ;
          call dacop(x(1),damap(1))
          call dacop(yp(1),damap(2))
          call dacop(x(2),damap(3))
          call dacop(yp(2),damap(4))
          do j=1,2
            ii=2*j
            call dapek(damap(ii-1),jj,c(j))
            call dapek(damap(ii),jj,cp(j))
          enddo
          call dacsu(damap(1),c(1),damap(1))
          call dacsu(damap(2),cp(1),damap(2))
          call dacsu(damap(3),c(2),damap(3))
          call dacsu(damap(4),cp(2),damap(4))
          if(ndimf.eq.3) then
            call dacop(sigmda,damap(5))
            call dacop(dpda1,damap(6))
            call dapek(damap(5),jj,c(3))
            call dapek(damap(6),jj,cp(3))
            call dacsu(damap(5),c(3),damap(5))
            call dacsu(damap(6),cp(3),damap(6))
            if(iflag1.eq.1.and.ithick.eq.1) then
              call dacct(damap,nvar,corrnew,nvar,damap,nvar)
            endif
          else
            call dacop(dpda1,damap(5))
            do j1=1,4
              do ii=1,4
                jj(ii)=1
                call dapek(damap(j1),jj,rdd(j1,ii))
                jj(ii)=0
              enddo
            enddo
            jj(5)=1
            do j1=1,4
              call dapek(damap(j1),jj,rdd(j1,5))
            enddo
            jj(5)=0
            do j1=1,2
              ii=2*j1
!hr03         d(j1)=rdd(ii-1,1)*dicu(1)+rdd(ii-1,2)*dicu(2)+            &
!hr03&rdd(ii-1,3)*dicu(3)+rdd(ii-1,4)*dicu(4)+rdd(ii-1,5)
              d(j1)=(((rdd(ii-1,1)*dicu(1)+rdd(ii-1,2)*dicu(2))+        &!hr03
     &rdd(ii-1,3)*dicu(3))+rdd(ii-1,4)*dicu(4))+rdd(ii-1,5)              !hr03
!hr03         dp(j1)=rdd(ii,1)*dicu(1)+rdd(ii,2)*dicu(2)+               &
!hr03&rdd(ii,3)*dicu(3)+rdd(ii,4)*dicu(4)+rdd(ii,5)
              dp(j1)=(((rdd(ii,1)*dicu(1)+rdd(ii,2)*dicu(2))+           &!hr03
     &rdd(ii,3)*dicu(3))+rdd(ii,4)*dicu(4))+rdd(ii,5)                    !hr03
            enddo
          endif
          call dacct(damap,nvar,aa2,nvar,damap,nvar)
          do j=1,ndimf
            ii=2*j
            if(j.eq.1) then
              i2=4
              i3=6
            elseif(j.eq.2) then
              i2=2
              i3=6
            elseif(j.eq.3) then
              i2=2
              i3=4
            endif
            jj(ii-1)=1
            call dapek(damap(ii-1),jj,angp(1,ii-1))
            call dapek(damap(ii),jj,au(ii,ii-1))
            jj(ii-1)=0
            jj(ii)=1
            call dapek(damap(ii-1),jj,angp(1,ii))
            call dapek(damap(ii),jj,au(ii,ii))
            jj(ii)=0
            jj(i2-1)=1
            call dapek(damap(ii-1),jj,au(i2-1,i2-1))
            call dapek(damap(ii),jj,au(i2,i2-1))
            jj(i2-1)=0
            jj(i2)=1
            call dapek(damap(ii-1),jj,au(i2-1,i2))
            call dapek(damap(ii),jj,au(i2,i2))
            jj(i2)=0
            jj(i3-1)=1
            call dapek(damap(ii-1),jj,au(i3-1,i3-1))
            call dapek(damap(ii),jj,au(i3,i3-1))
            jj(i3-1)=0
            jj(i3)=1
            call dapek(damap(ii-1),jj,au(i3-1,i3))
            call dapek(damap(ii),jj,au(i3,i3))
            jj(i3)=0
!hr08       b1(j)=angp(1,ii-1)*angp(1,ii-1)+angp(1,ii)*angp(1,ii)
            b1(j)=angp(1,ii-1)**2+angp(1,ii)**2                          !hr08
!hr08       b2(j)=au(i2-1,i2-1)*au(i2-1,i2-1)+au(i2-1,i2)*au(i2-1,i2)
            b2(j)=au(i2-1,i2-1)**2+au(i2-1,i2)**2                        !hr08
!hr08       b3(j)=au(i3-1,i3-1)*au(i3-1,i3-1)+au(i3-1,i3)*au(i3-1,i3)
            b3(j)=au(i3-1,i3-1)**2+au(i3-1,i3)**2                        !hr08
!hr03       al1(j)=-(angp(1,ii-1)*au(ii,ii-1)+angp(1,ii)*au(ii,ii))
            al1(j)=-1d0*(angp(1,ii-1)*au(ii,ii-1)+angp(1,ii)*au(ii,ii))  !hr03
!hr03       al2(j)=-(au(i2-1,i2-1)*au(i2,i2-1)+au(i2-1,i2)*au(i2,i2))
           al2(j)=-1d0*(au(i2-1,i2-1)*au(i2,i2-1)+au(i2-1,i2)*au(i2,i2)) !hr03
!hr03       al3(j)=-(au(i3-1,i3-1)*au(i3,i3-1)+au(i3-1,i3)*au(i3,i3))
           al3(j)=-1d0*(au(i3-1,i3-1)*au(i3,i3-1)+au(i3-1,i3)*au(i3,i3)) !hr03
!hr04       g1(j)=au(ii,ii-1)*au(ii,ii-1)+au(ii,ii)*au(ii,ii)
!hr04       g2(j)=au(i2,i2-1)*au(i2,i2-1)+au(i2,i2)*au(i2,i2)
!hr04       g3(j)=au(i3,i3-1)*au(i3,i3-1)+au(i3,i3)*au(i3,i3)
            g1(j)=au(ii,ii-1)**2+au(ii,ii)**2                            !hr04
            g2(j)=au(i2,i2-1)**2+au(i2,i2)**2                            !hr04
            g3(j)=au(i3,i3-1)**2+au(i3,i3)**2                            !hr04
            if(ndimf.eq.3) then
              call dainv(damap,nvar,damapi,nvar)
              jj(6)=1
              call dapek(damapi(5),jj,aui(1))
              call dapek(damapi(6),jj,aui(2))
              jj(6)=0
              if(j.lt.3) then
                d(j)=au(i3-1,i3-1)*aui(1)+au(i3-1,i3)*aui(2)
                dp(j)=au(i3,i3-1)*aui(1)+au(i3,i3)*aui(2)
              else
                d(j)=angp(1,ii-1)*aui(1)+angp(1,ii)*aui(2)
                dp(j)=au(ii,ii-1)*aui(1)+au(ii,ii)*aui(2)
              endif
            endif
            sx=angp(2,ii-1)*angp(1,ii)-angp(1,ii-1)*angp(2,ii)
            cx=angp(1,ii-1)*angp(2,ii-1)+angp(1,ii)*angp(2,ii)
            if(abs(sx).gt.c1m15.or.abs(cx).gt.c1m15) then
              dphi(j)=atan2_rn(sx,cx)/x2pi
            else
              dphi(j)=zero
            endif
            phi(j)=phi(j)+dphi(j)
          enddo
          do j=1,ndimf
            ii=2*j
            angp(2,ii-1)=angp(1,ii-1)
            angp(2,ii)=angp(1,ii)
          enddo
          if(iwrite.eq.1) then
            iii=i
            if(typ(:8).eq.'START   ') iii=0
            write(*,10030) iii,typ(:8),tl,phi(1),b1(1),al1(1),g1(1),    &
     &d(1),dp(1),c(1),cp(1)
            if(ndimf.eq.3) then
              write(*,10040) b2(1),al2(1),g2(1)
              write(*,10050) typ(9:16),b3(1),al3(1),g3(1)
            else
              write(*,10055) typ(9:16),b2(1),al2(1),g2(1)
            endif
            write(*,10060)
            write(*,10070) phi(2),b1(2),al1(2),g1(2),d(2),dp(2),c(2),   &
     &cp(2)
            write(*,10080) b2(2),al2(2),g2(2)
            if(ndimf.eq.3) then
              write(*,10090) b3(2),al3(2),g3(2)
              write(*,10060)
              write(*,10100) -phi(3),b1(3),al1(3),g1(3),d(3),dp(3),c(3),&
     &cp(3)*c1m3
              write(*,10080) b2(3),al2(3),g2(3)
              write(*,10040) b3(3),al3(3),g3(3)
            endif
            write(*,10010)
          endif
              if(i.eq.nt) goto 470
            endif
          enddo
        endif
        goto 430
   50   ix=ix-nblo
        if(abs(dare(x(1))).gt.aint(aper(1)).or.                         &
     &abs(dare(x(2))).gt.aint(aper(2))) then
          write(*,10120) j,i,dare(x(1)),aper(1),dare(x(2)),aper(2),ix,  &
     &kz(ix),bez(ix)
          call prror(97)
        endif
        kpz=abs(kp(ix))
        if(kpz.eq.0) goto 80
        goto(80,80,80,80,80,70),kpz
        goto 430
   70   continue
        if(ition.ne.0) then
*FOX  EJF0=EJF1 ;
          ixcav=ix
          if(abs(dppoff).gt.pieni) then
            sigmdac=dare(sigmda)
            sigmoff(i)=sigmdac
*FOX  SIGMDA=SIGMDA-SIGMDAC ;
          endif
          call synoda
*FOX  DPDA=DPDA1*C1M3 ;
*FOX  Y(1)=EJF0/EJF1*Y(1) ;
*FOX  Y(2)=EJF0/EJF1*Y(2) ;
          if(nvar2.eq.6.and.nsix.ne.2) then
            iflag=1
            iflag1=1
            iflag2=1
          endif
        endif
        goto 440
   80   kzz=kz(ix)
        if(kzz.eq.15) then
          ixcav=ix
*FOX  XX(1)=X(1) ;
*FOX  XX(2)=X(2) ;
*FOX  YY(1)=Y(1) ;
*FOX  YY(2)=Y(2) ;
          call wireda
*FOX  X(1)=XX(1) ;
*FOX  X(2)=XX(2) ;
*FOX  Y(1)=YY(1) ;
*FOX  Y(2)=YY(2) ;
          goto 440
        endif
        if(ilinc.eq.2.and.kzz.eq.20) then
          if(nbeam.ge.1) then
          ibb=ibb+1
          if(ibb.gt.nbb) call prror(102)
          imbb(i)=ibb
*FOX  YP(1)=Y(1)*(ONE+DPDA) ;
*FOX  YP(2)=Y(2)*(ONE+DPDA) ;
*FOX  DPDA1=DPDA*C1E3 ;
          call dacop(x(1),damap(1))
          call dacop(yp(1),damap(2))
          call dacop(x(2),damap(3))
          call dacop(yp(2),damap(4))
          do j=1,2
            ii=2*j
            call dapek(damap(ii-1),jj,c(j))
            call dapek(damap(ii),jj,cp(j))
          enddo
          call dacsu(damap(1),c(1),damap(1))
          call dacsu(damap(2),cp(1),damap(2))
          call dacsu(damap(3),c(2),damap(3))
          call dacsu(damap(4),cp(2),damap(4))
          if(ndimf.eq.3) then
            call dacop(sigmda,damap(5))
            call dacop(dpda1,damap(6))
            call dapek(damap(5),jj,c(3))
            call dapek(damap(6),jj,cp(3))
            call dacsu(damap(5),c(3),damap(5))
            call dacsu(damap(6),cp(3),damap(6))
            if(iflag2.eq.1.and.ithick.eq.1) then
              call dacct(damap,nvar,corrnew,nvar,damap,nvar)
            endif
          endif
          call dainv(damap,nvar,damapi,nvar)
          call dacct(damap,nvar,aa2,nvar,aa2r,nvar)
          call dacct(damap,nvar,damap1,nvar,damap,nvar)
          call dacct(damap,nvar,damapi,nvar,damap,nvar)
            do ii=1,3
              call damul(damap(i4(ii,1)),damap(i4(ii,2)),angno)
              call averaged(angno,aa2r,.false.,angno)
              do j=1,ndimf
                j1=2*j
                jj(j1-1)=1
                jj(j1)=1
                call dapek(angno,jj,angnoe(j))
                jj(j1-1)=0
                jj(j1)=0
              enddo
              if(ndimf.eq.3) then
!hr03           bbcu(ibb,ii)=two*(emitx*angnoe(1)+emity*angnoe(2)+      &
                bbcu(ibb,ii)=two*((emitx*angnoe(1)+emity*angnoe(2))+    &!hr03
     &emitz*angnoe(3))
              else
                bbcu(ibb,ii)=two*(emitx*angnoe(1)+emity*angnoe(2))
              endif
            enddo
!hr12       if(parbe(ix,2).gt.0) then
            if(parbe(ix,2).gt.0d0) then                                  !hr12
              do ii=4,10
                call damul(damap(i4(ii,1)),damap(i4(ii,2)),angno)
                call averaged(angno,aa2r,.false.,angno)
                do j=1,ndimf
                  j1=2*j
                  jj(j1-1)=1
                  jj(j1)=1
                  call dapek(angno,jj,angnoe(j))
                  jj(j1-1)=0
                  jj(j1)=0
                enddo
                if(ndimf.eq.3) then
!hr03           bbcu(ibb,ii)=two*(emitx*angnoe(1)+emity*angnoe(2)+      &
                bbcu(ibb,ii)=two*((emitx*angnoe(1)+emity*angnoe(2))+    &!hr03
     &emitz*angnoe(3))
                else
                  bbcu(ibb,ii)=two*(emitx*angnoe(1)+emity*angnoe(2))
                endif
              enddo
            endif
          if(lhc.eq.1) then
            dummy=bbcu(ibb,1)
            bbcu(ibb,1)=bbcu(ibb,2)
            bbcu(ibb,2)=dummy
            dummy=bbcu(ibb,4)
            bbcu(ibb,4)=bbcu(ibb,9)
            bbcu(ibb,9)=dummy
            dummy=bbcu(ibb,5)
            bbcu(ibb,5)=bbcu(ibb,7)
            bbcu(ibb,7)=dummy
            dummy=bbcu(ibb,6)
            bbcu(ibb,6)=bbcu(ibb,10)
            bbcu(ibb,10)=dummy
          endif
          if(lhc.eq.2) then
            bbcu(ibb,1)=bbbx(ix)
            bbcu(ibb,2)=bbby(ix)
            bbcu(ibb,3)=bbbs(ix)
          endif
          if((bbcu(ibb,1).le.pieni).or.(bbcu(ibb,2).le.pieni)) then
            call prror(88)
          endif
          if(ibbc.eq.1) then
            sfac1=bbcu(ibb,1)+bbcu(ibb,2)
            sfac2=bbcu(ibb,1)-bbcu(ibb,2)
            sfac2s=one
!hr08       if(sfac2.lt.zero) sfac2s=-one
            if(sfac2.lt.zero) sfac2s=-1d0*one                            !hr08
!hr03       sfac3=sqrt(sfac2*sfac2+four*bbcu(ibb,3)*bbcu(ibb,3))
            sfac3=sqrt(sfac2**2+(four*bbcu(ibb,3))*bbcu(ibb,3))          !hr03
            if(sfac3.gt.sfac1) call prror(103)
!hr03       sfac4=sfac2s*sfac2/sfac3
            sfac4=(sfac2s*sfac2)/sfac3                                   !hr03
!hr03       sfac5=-sfac2s*two*bbcu(ibb,3)/sfac3
            sfac5=(((-1d0*sfac2s)*two)*bbcu(ibb,3))/sfac3                !hr03
!hr03       sigman(1,ibb)=sqrt((sfac1+sfac2*sfac4+                      &
!hr03&two*bbcu(ibb,3)*sfac5)*half)
            sigman(1,ibb)=sqrt(((sfac1+sfac2*sfac4)+                    &!hr03
     &(two*bbcu(ibb,3))*sfac5)*half)                                     !hr03
!hr03       sigman(2,ibb)=sqrt((sfac1-sfac2*sfac4-                      &
!hr03&two*bbcu(ibb,3)*sfac5)*half)
            sigman(2,ibb)=sqrt(((sfac1-sfac2*sfac4)-                    &!hr03
     &(two*bbcu(ibb,3))*sfac5)*half)                                     !hr03
            bbcu(ibb,11)=sqrt(half*(one+sfac4))
!hr03       bbcu(ibb,12)=-sfac2s*sqrt(half*(one-sfac4))
            bbcu(ibb,12)=(-1d0*sfac2s)*sqrt(half*(one-sfac4))            !hr03
!hr03       if(bbcu(ibb,3).lt.zero) bbcu(ibb,12)=-bbcu(ibb,12)
            if(bbcu(ibb,3).lt.zero) bbcu(ibb,12)=-1d0*bbcu(ibb,12)       !hr03
          else
            bbcu(ibb,11)=one
            sigman(1,ibb)=sqrt(bbcu(ibb,1))
            sigman(2,ibb)=sqrt(bbcu(ibb,2))
          endif
!hr08     if(parbe(ix,2).gt.0) then
          if(parbe(ix,2).gt.0d0) then                                    !hr08
            do ii=1,10
              bbcu(ibb,ii)=bbcu(ibb,ii)*c1m6
            enddo
          endif
          endif
          goto 440
        endif
        if(kzz.eq.20.and.iqmodc.eq.4) goto 440
!hr12   if(kzz.eq.20.and.parbe(ix,2).eq.0) then
        if(kzz.eq.20.and.parbe(ix,2).eq.0d0) then                        !hr12
          if(nbeam.ge.1) then
            if(ilinc.eq.0) then
              clobeam(1,imbb(i))=dare(x(1))
              clobeam(2,imbb(i))=dare(x(2))
              clobeam(4,imbb(i))=dare(y(1))*(one+dare(dpda))
              clobeam(5,imbb(i))=dare(y(2))*(one+dare(dpda))
              if(ndimf.eq.3) then
                clobeam(3,imbb(i))=dare(sigmda)
                clobeam(6,imbb(i))=dare(dpda)
              endif
            endif
            if(sigman(1,imbb(i)).eq.sigman(2,imbb(i))) then
              if(ibeco.eq.1) then                                       &
            if(ibbc.eq.0) then
              crk=ed(ix)
              cik=ek(ix)
            else
              crk=ed(ix)*bbcu(imbb(i),11)+ek(ix)*bbcu(imbb(i),12)
!hr03         cik=-ed(ix)*bbcu(imbb(i),12)+ek(ix)*bbcu(imbb(i),11)
              cik=ek(ix)*bbcu(imbb(i),11)-ed(ix)*bbcu(imbb(i),12)        !hr03
            endif
!hr03       rho2b=crk*crk+cik*cik
            rho2b=crk**2+cik**2                                          !hr03
            if(rho2b.gt.pieni)
     &then
            if(abs(sigman(1,imbb(i))).lt.pieni) call prror(88)
!hr03       tkb=rho2b/(two*sigman(1,imbb(i))*sigman(1,imbb(i)))
            tkb=rho2b/((two*sigman(1,imbb(i)))*sigman(1,imbb(i)))        !hr03
!hr03       beamoff4=crad*ptnfac(ix)*crk/                               &
!hr03       beamoff4=crad*ptnfac(ix)*crk/                               &
!hr03&rho2b*(one-exp_rn(-tkb))
            beamoff4=(((crad*ptnfac(ix))*crk)/                          &!hr03
     &rho2b)*(one-exp_rn(-1d0*tkb))                                      !hr03
!hr03       beamoff5=crad*ptnfac(ix)*cik/                               &
!hr03       beamoff5=crad*ptnfac(ix)*cik/                               &
!hr03&rho2b*(one-exp_rn(-tkb))
            beamoff5=(((crad*ptnfac(ix))*cik)/                          &!hr03
     &rho2b)*(one-exp_rn(-1d0*tkb))                                      !hr03
                endif
              endif
*FOX  CRKVEBF=X(1) ;
*FOX  CIKVEBF=X(2) ;
!hr03       startco=dare(x(1))-clobeam(1,imbb(i))+ed(ix)
            startco=(dare(x(1))-clobeam(1,imbb(i)))+ed(ix)               !hr03
            call dapok(crkvebf,jj,startco)
!hr03       startco=dare(x(2))-clobeam(2,imbb(i))+ek(ix)
            startco=(dare(x(2))-clobeam(2,imbb(i)))+ek(ix)
            call dapok(cikvebf,jj,startco)
            if(ibbc.eq.1) then
*FOX  CCCC=CRKVEBF ;
*FOX  CRKVEBF=CCCC*BBCU(IMBB(I),11)+CIKVEBF*BBCU(IMBB(I),12) ;
*FOX  CIKVEBF=-CCCC*BBCU(IMBB(I),12)+CIKVEBF*BBCU(IMBB(I),11) ;
            endif
*FOX  RHO2BF=CRKVEBF*CRKVEBF+CIKVEBF*CIKVEBF ;
              if(abs(dare(rho2bf)).gt.pieni) then
      if(abs(sigman(1,imbb(i))).lt.pieni) call prror(88)
*FOX  TKBF=RHO2BF/(TWO*SIGMAN(1,IMBB(I))*SIGMAN(1,IMBB(I))) ;
      if(ibbc.eq.0) then
*FOX   Y(1)=Y(1)+(CRAD*CRKVEBF/RHO2BF*
*FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF4)/(ONE+DPDA) ;
*FOX   Y(2)=Y(2)+(CRAD*CIKVEBF/RHO2BF*
*FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF5)/(ONE+DPDA) ;
      else
*FOX   CCCC=(CRAD*CRKVEBF/RHO2BF*
*FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF4)*BBCU(IMBB(I),11)-
*FOX   (CRAD*CIKVEBF/RHO2BF*
*FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF5)*BBCU(IMBB(I),12) ;
*FOX   Y(1)=Y(1)+CCCC/(ONE+DPDA) ;
*FOX   CCCC=(CRAD*CRKVEBF/RHO2BF*
*FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF4)*BBCU(IMBB(I),12)+
*FOX   (CRAD*CIKVEBF/RHO2BF*
*FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF5)*BBCU(IMBB(I),11) ;
*FOX   Y(2)=Y(2)+CCCC/(ONE+DPDA) ;
      endif
              endif
            else if(sigman(1,imbb(i)).gt.sigman(2,imbb(i))) then
              if(ibeco.eq.1) then
            if(abs(sigman(1,imbb(i))).lt.pieni.or.                      &
     &abs(sigman(2,imbb(i))).lt.pieni) call prror(88)
!hr08       r2b=two*(sigman(1,imbb(i))*sigman(1,imbb(i))-               &
!hr08&sigman(2,imbb(i))*sigman(2,imbb(i)))
            r2b=two*(sigman(1,imbb(i))**2-                              &!hr08
     &sigman(2,imbb(i))**2)                                              !hr08
            rb=sqrt(r2b)
!hr03       rkb=crad*ptnfac(ix)*pisqrt/rb
            rkb=((crad*ptnfac(ix))*pisqrt)/rb                            !hr03
            if(ibbc.eq.0) then
              crk=ed(ix)
              cik=ek(ix)
            else
              crk=ed(ix)*bbcu(imbb(i),11)+ek(ix)*bbcu(imbb(i),12)
!hr03         cik=-ed(ix)*bbcu(imbb(i),12)+ek(ix)*bbcu(imbb(i),11)
              cik=ek(ix)*bbcu(imbb(i),11)-ed(ix)*bbcu(imbb(i),12)        !hr03
            endif
            xrb=abs(crk)/rb
            zrb=abs(cik)/rb
            call errf(xrb,zrb,crxb,crzb)
            if(abs(sigman(1,imbb(i))).lt.pieni.or.                      &
     &abs(sigman(2,imbb(i))).lt.pieni) call prror(88)
!hr03       tkb=(crk*crk/(sigman(1,imbb(i))*sigman(1,imbb(i)))+         &
!hr03&cik*cik/(sigman(2,imbb(i))*sigman(2,imbb(i))))*half
            tkb=(crk**2/sigman(1,imbb(i))**2+                           &!hr03
     &cik**2/sigman(2,imbb(i))**2)*half                                  !hr03
!hr03       xbb=sigman(2,imbb(i))/sigman(1,imbb(i))*xrb
            xbb=(sigman(2,imbb(i))/sigman(1,imbb(i)))*xrb                !hr03
!hr03       zbb=sigman(1,imbb(i))/sigman(2,imbb(i))*zrb
            zbb=(sigman(1,imbb(i))/sigman(2,imbb(i)))*zrb                !hr03
            call errf(xbb,zbb,cbxb,cbzb)
!hr03         beamoff4=rkb*(crzb-exp_rn(-tkb)*cbzb)*
!hr03&sign(one,crk)
              beamoff4=(rkb*(crzb-exp_rn(-1d0*tkb)*cbzb))*              &!hr03
     &sign(one,crk)                                                      !hr03
!hr03&sign(one,crk)
!hr03         beamoff5=rkb*(crxb-exp_rn(-tkb)*cbxb)*
!hr03&sign(one,cik)
              beamoff5=(rkb*(crxb-exp_rn(-1d0*tkb)*cbxb))*              &!hr03
     &sign(one,cik)                                                      !hr03
!hr03&sign(one,cik)
              endif
            if(abs(sigman(1,imbb(i))).lt.pieni.or.                      &
     &abs(sigman(2,imbb(i))).lt.pieni) call prror(88)
!hr08       r2bf=two*(sigman(1,imbb(i))*sigman(1,imbb(i))-              &
!hr08&sigman(2,imbb(i))*sigman(2,imbb(i)))
            r2bf=two*(sigman(1,imbb(i))**2-                             &!hr08
     &sigman(2,imbb(i))**2)                                              !hr08
            rbf=sqrt(r2bf)
!hr03       rkbf=crad*ptnfac(ix)*pisqrt/rbf
            rkbf=((crad*ptnfac(ix))*pisqrt)/rbf                          !hr03
*FOX  CRKVEBF=X(1) ;
*FOX  CIKVEBF=X(2) ;
!hr03       startco=dare(x(1))-clobeam(1,imbb(i))+ed(ix)
            startco=(dare(x(1))-clobeam(1,imbb(i)))+ed(ix)               !hr03
            call dapok(crkvebf,jj,startco)
!hr03       startco=dare(x(2))-clobeam(2,imbb(i))+ek(ix)
            startco=(dare(x(2))-clobeam(2,imbb(i)))+ek(ix)
            call dapok(cikvebf,jj,startco)
            if(ibbc.eq.1) then
*FOX  CCCC=CRKVEBF ;
*FOX  CRKVEBF=CCCC*BBCU(IMBB(I),11)+CIKVEBF*BBCU(IMBB(I),12) ;
*FOX  CIKVEBF=-CCCC*BBCU(IMBB(I),12)+CIKVEBF*BBCU(IMBB(I),11) ;
            endif
*FOX  XRBF=CRKVEBF/RBF ;
      if(dare(xrbf).lt.zero) then
*FOX  XRBF=-XRBF ;
      endif
*FOX  ZRBF=CIKVEBF/RBF ;
      if(dare(zrbf).lt.zero) then
*FOX  ZRBF=-ZRBF ;
      endif
            call errff(xrbf,zrbf,crxbf,crzbf)
      if(abs(sigman(1,imbb(i))).lt.pieni.or.                            &
     &abs(sigman(2,imbb(i))).lt.pieni) call prror(88)
*FOX  TKBF=(CRKVEBF*CRKVEBF/(SIGMAN(1,IMBB(I))*SIGMAN(1,IMBB(I)))+
*FOX  CIKVEBF*CIKVEBF/(SIGMAN(2,IMBB(I))*SIGMAN(2,IMBB(I))))*HALF ;
*FOX  XBBF=SIGMAN(2,IMBB(I))/SIGMAN(1,IMBB(I))*XRBF ;
*FOX  ZBBF=SIGMAN(1,IMBB(I))/SIGMAN(2,IMBB(I))*ZRBF ;
            call errff(xbbf,zbbf,cbxbf,cbzbf)
      scrkveb=sign(one,dare(crkvebf))
      scikveb=sign(one,dare(cikvebf))
      if(ibbc.eq.0) then
*FOX  Y(1)=Y(1)+(RKBF*(CRZBF-EXP(-TKBF)*
*FOX  CBZBF)*SCRKVEB-BEAMOFF4)/(ONE+DPDA) ;
*FOX  Y(2)=Y(2)+(RKBF*(CRXBF-EXP(-TKBF)*
*FOX  CBXBF)*SCIKVEB-BEAMOFF5)/(ONE+DPDA) ;
      else
*FOX  CCCC=(RKBF*(CRZBF-EXP(-TKBF)*
*FOX  CBZBF)*SCRKVEB-BEAMOFF4)*BBCU(IMBB(I),11)-
*FOX  (RKBF*(CRXBF-EXP(-TKBF)*
*FOX  CBXBF)*SCIKVEB-BEAMOFF5)*BBCU(IMBB(I),12) ;
*FOX   Y(1)=Y(1)+CCCC/(ONE+DPDA) ;
*FOX  CCCC=(RKBF*(CRZBF-EXP(-TKBF)*
*FOX  CBZBF)*SCRKVEB-BEAMOFF4)*BBCU(IMBB(I),12)+
*FOX  (RKBF*(CRXBF-EXP(-TKBF)*
*FOX  CBXBF)*SCIKVEB-BEAMOFF5)*BBCU(IMBB(I),11) ;
*FOX   Y(2)=Y(2)+CCCC/(ONE+DPDA) ;
      endif
            else if(sigman(1,imbb(i)).lt.sigman(2,imbb(i))) then
              if(ibeco.eq.1) then
            if(abs(sigman(1,imbb(i))).lt.pieni.or.                      &
     &abs(sigman(2,imbb(i))).lt.pieni) call prror(88)
!hr08       r2b=two*(sigman(2,imbb(i))*sigman(2,imbb(i))-               &
!hr08&sigman(1,imbb(i))*sigman(1,imbb(i)))
            r2b=two*(sigman(2,imbb(i))**2-                              &!hr08
     &sigman(1,imbb(i))**2)                                              !hr08
            rb=sqrt(r2b)
!hr03       rkb=crad*ptnfac(ix)*pisqrt/rb
            rkb=((crad*ptnfac(ix))*pisqrt)/rb                            !hr03
            if(ibbc.eq.0) then
              crk=ed(ix)
              cik=ek(ix)
            else
              crk=ed(ix)*bbcu(imbb(i),11)+ek(ix)*bbcu(imbb(i),12)
!hr03         cik=-ed(ix)*bbcu(imbb(i),12)+ek(ix)*bbcu(imbb(i),11)
              cik=ek(ix)*bbcu(imbb(i),11)-ed(ix)*bbcu(imbb(i),12)        !hr03
            endif
            xrb=abs(crk)/rb
            zrb=abs(cik)/rb
            call errf(zrb,xrb,crzb,crxb)
            if(abs(sigman(1,imbb(i))).lt.pieni.or.                      &
     &abs(sigman(2,imbb(i))).lt.pieni) call prror(88)
!hr03       tkb=(crk*crk/(sigman(1,imbb(i))*sigman(1,imbb(i)))+         &
!hr03&cik*cik/(sigman(2,imbb(i))*sigman(2,imbb(i))))*half
            tkb=(crk**2/sigman(1,imbb(i))**2+                           &!hr03
     &cik**2/sigman(2,imbb(i))**2)*half                                  !hr03
!hr03       xbb=sigman(2,imbb(i))/sigman(1,imbb(i))*xrb
            xbb=(sigman(2,imbb(i))/sigman(1,imbb(i)))*xrb                !hr03
!hr03       zbb=sigman(1,imbb(i))/sigman(2,imbb(i))*zrb
            zbb=(sigman(1,imbb(i))/sigman(2,imbb(i)))*zrb                !hr03
            call errf(zbb,xbb,cbzb,cbxb)
!hr03         beamoff4=rkb*(crzb-exp_rn(-tkb)*cbzb)*
!hr03&sign(one,crk)
              beamoff4=(rkb*(crzb-exp_rn(-1d0*tkb)*cbzb))*              &!hr03
     &sign(one,crk)                                                      !hr03
!hr03&sign(one,crk)
!hr03         beamoff5=rkb*(crxb-exp_rn(-tkb)*cbxb)*
!hr03&sign(one,cik)
              beamoff5=(rkb*(crxb-exp_rn(-1d0*tkb)*cbxb))*              &!hr03
     &sign(one,cik)                                                      !hr03
!hr03&sign(one,cik)
              endif
            if(abs(sigman(1,imbb(i))).lt.pieni.or.                      &
     &abs(sigman(2,imbb(i))).lt.pieni) call prror(88)
!hr08       r2bf=two*(sigman(2,imbb(i))*sigman(2,imbb(i))-              &
!hr08&sigman(1,imbb(i))*sigman(1,imbb(i)))
            r2bf=two*(sigman(2,imbb(i))**2-                             &!hr08
     &sigman(1,imbb(i))**2)                                              !hr08
            rbf=sqrt(r2bf)
!hr03       rkbf=crad*ptnfac(ix)*pisqrt/rbf
            rkbf=((crad*ptnfac(ix))*pisqrt)/rbf                          !hr03
*FOX  CRKVEBF=X(1) ;
*FOX  CIKVEBF=X(2) ;
!hr03       startco=dare(x(1))-clobeam(1,imbb(i))+ed(ix)
            startco=(dare(x(1))-clobeam(1,imbb(i)))+ed(ix)               !hr03
            call dapok(crkvebf,jj,startco)
!hr03       startco=dare(x(2))-clobeam(2,imbb(i))+ek(ix)
            startco=(dare(x(2))-clobeam(2,imbb(i)))+ek(ix)
            call dapok(cikvebf,jj,startco)
            if(ibbc.eq.1) then
*FOX  CCCC=CRKVEBF ;
*FOX  CRKVEBF=CCCC*BBCU(IMBB(I),11)+CIKVEBF*BBCU(IMBB(I),12) ;
*FOX  CIKVEBF=-CCCC*BBCU(IMBB(I),12)+CIKVEBF*BBCU(IMBB(I),11) ;
            endif
*FOX  XRBF=CRKVEBF/RBF ;
      if(dare(xrbf).lt.zero) then
*FOX  XRBF=-XRBF ;
      endif
*FOX  ZRBF=CIKVEBF/RBF ;
      if(dare(zrbf).lt.zero) then
*FOX  ZRBF=-ZRBF ;
      endif
            call errff(zrbf,xrbf,crzbf,crxbf)
      if(abs(sigman(1,imbb(i))).lt.pieni.or.                            &
     &abs(sigman(2,imbb(i))).lt.pieni) call prror(88)
*FOX  TKBF=(CRKVEBF*CRKVEBF/(SIGMAN(1,IMBB(I))*SIGMAN(1,IMBB(I)))+
*FOX  CIKVEBF*CIKVEBF/(SIGMAN(2,IMBB(I))*SIGMAN(2,IMBB(I))))*HALF ;
*FOX  XBBF=SIGMAN(2,IMBB(I))/SIGMAN(1,IMBB(I))*XRBF ;
*FOX  ZBBF=SIGMAN(1,IMBB(I))/SIGMAN(2,IMBB(I))*ZRBF ;
            call errff(zbbf,xbbf,cbzbf,cbxbf)
      scrkveb=sign(one,dare(crkvebf))
      scikveb=sign(one,dare(cikvebf))
      if(ibbc.eq.0) then
*FOX  Y(1)=Y(1)+(RKBF*(CRZBF-EXP(-TKBF)*
*FOX  CBZBF)*SCRKVEB-BEAMOFF4)/(ONE+DPDA) ;
*FOX  Y(2)=Y(2)+(RKBF*(CRXBF-EXP(-TKBF)*
*FOX  CBXBF)*SCIKVEB-BEAMOFF5)/(ONE+DPDA) ;
      else
*FOX  CCCC=(RKBF*(CRZBF-EXP(-TKBF)*
*FOX  CBZBF)*SCRKVEB-BEAMOFF4)*BBCU(IMBB(I),11)-
*FOX  (RKBF*(CRXBF-EXP(-TKBF)*
*FOX  CBXBF)*SCIKVEB-BEAMOFF5)*BBCU(IMBB(I),12) ;
*FOX   Y(1)=Y(1)+CCCC/(ONE+DPDA) ;
*FOX  CCCC=(RKBF*(CRZBF-EXP(-TKBF)*
*FOX  CBZBF)*SCRKVEB-BEAMOFF4)*BBCU(IMBB(I),12)+
*FOX  (RKBF*(CRXBF-EXP(-TKBF)*
*FOX  CBXBF)*SCIKVEB-BEAMOFF5)*BBCU(IMBB(I),11) ;
*FOX   Y(2)=Y(2)+CCCC/(ONE+DPDA) ;
      endif
            endif
            goto 440
          endif
          goto 440
        endif
!hr12   if(kzz.eq.20.and.parbe(ix,2).gt.0) then
        if(kzz.eq.20.and.parbe(ix,2).gt.0d0) then                        !hr12
          if(ilinc.eq.0)then
            clobeam(1,imbb(i))=dare(x(1))
            clobeam(2,imbb(i))=dare(x(2))
            clobeam(4,imbb(i))=dare(y(1))*(one+dare(dpda))
            clobeam(5,imbb(i))=dare(y(2))*(one+dare(dpda))
            if(ndimf.eq.3)then
              clobeam(3,imbb(i))=dare(sigmda)
              clobeam(6,imbb(i))=dare(dpda)
            endif
          endif
!hr08     parbe(ix,4)=-crad*ptnfac(ix)*half*c1m6
          parbe(ix,4)=(((-1d0*crad)*ptnfac(ix))*half)*c1m6               !hr08
!--Hirata's 6D beam-beam kick
          dummy=dare(x(1))
*FOX      TRACKI(1)=(X(1)+ED(IX)-DUMMY)*C1M3 ;
*FOX      YP(1)=Y(1)*(ONE+DPDA) ;
          dummy=dare(yp(1))
*FOX      TRACKI(2)=(YP(1)-DUMMY)*C1M3 ;
          dummy=dare(x(2))
*FOX      TRACKI(3)=(X(2)+EK(IX)-DUMMY)*C1M3 ;
*FOX      YP(2)=Y(2)*(ONE+DPDA) ;
          dummy=dare(yp(2))
*FOX      TRACKI(4)=(YP(2)-DUMMY)*C1M3 ;
          dummy=dare(sigmda)
*FOX      TRACKI(5)=(SIGMDA-DUMMY)*C1M3 ;
          dummy=dare(dpda)
*FOX      TRACKI(6)=DPDA-DUMMY ;
          call beaminf(tracki,parbe,sigz,bbcu,imbb(i),ix,ibbc)
          if(ibeco.eq.1) then
            beamoff1=dare(tracki(1))*c1e3
            beamoff2=dare(tracki(3))*c1e3
            beamoff4=dare(tracki(2))*c1e3
            beamoff5=dare(tracki(4))*c1e3
            beamoff6=dare(tracki(6))
          else
            beamoff1=zero
            beamoff2=zero
            beamoff4=zero
            beamoff5=zero
            beamoff6=zero
          endif
          dummy=dare(x(1))
*FOX      X(1)=TRACKI(1)*C1E3+DUMMY-BEAMOFF1 ;
          dummy=dare(x(2))
*FOX      X(2)=TRACKI(3)*C1E3+DUMMY-BEAMOFF2 ;
          dummy=dare(yp(1))
*FOX      YP(1)=TRACKI(2)*C1E3+DUMMY-BEAMOFF4 ;
          dummy=dare(yp(2))
*FOX      YP(2)=TRACKI(4)*C1E3+DUMMY-BEAMOFF5 ;
          dummy=dare(dpda)
*FOX      DPDA=TRACKI(6)+DUMMY-BEAMOFF6 ;
*FOX      DPDA1=DPDA*C1E3 ;
*FOX      Y(1)=YP(1)/(ONE+DPDA) ;
*FOX      Y(2)=YP(2)/(ONE+DPDA) ;
*FOX      EJF1=E0F*(ONE+DPDA) ;
*FOX      EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;
*FOX      RV=EJ1/E0*E0F/EJF1 ;
          if(ithick.eq.1) call envada
          goto 440
        endif
          pi=4d0*atan_rn(1d0)
        if(kzz.eq.23) then
*FOX  CRABAMP=ED(IX)/(EJF1) ;
!       call dapri(EJF1,234)
!       write(*,*) crabamp, EJF1, EJF0,clight, "HELLO"
        crabfreq=ek(ix)*c1e3
        crabpht=crabph(ix)
*FOX  Y(1)=Y(1) - CRABAMP*C1E3*
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT) ;
*FOX  DPDA1=DPDA1 - CRABAMP*CRABFREQ*2D0*PI/CLIGHT*X(1)*
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT) ;
*FOX  EJF0=EJF1 ;
*FOX  DPDA=DPDA1*C1M3 ;
*FOX  EJF1=E0F*(ONE+DPDA) ;
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;
*FOX  Y(1)=EJF0/EJF1*Y(1) ;
*FOX  Y(2)=EJF0/EJF1*Y(2) ;
          goto 440
      endif
        if(kzz.eq.-23) then
*FOX  CRABAMP=ED(IX)/(EJF1) ;
           crabfreq=ek(ix)*c1e3
           crabpht=crabph(ix)
*FOX  Y(2)=Y(2) - CRABAMP*C1E3*
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT) ;
*FOX  DPDA1=DPDA1 - CRABAMP*CRABFREQ*2D0*PI/CLIGHT*X(2)*
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT) ;
*FOX  EJF0=EJF1 ;
*FOX  DPDA=DPDA1*C1M3 ;
*FOX  EJF1=E0F*(ONE+DPDA) ;
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;
*FOX  Y(1)=EJF0/EJF1*Y(1) ;
*FOX  Y(2)=EJF0/EJF1*Y(2) ;
          goto 440
      endif
 
! JBG RF CC Multipoles
        if(kzz.eq.26) then
            ! JBG bypass this element if 4D/5D case
            if(iclo6.eq.0) then
!                write(*,*)'Bypassing RF mult 4D or 5D case'
                goto 440
            endif
          xs=xsi(i) ! JBG change of variables for misal calculations
          zs=zsi(i)
*FOX  XL=(X(1)-XS)*TILTC(I)+(X(2)-ZS)*TILTS(I) ;
*FOX  ZL=-(X(1)-XS)*TILTS(I)+(X(2)-ZS)*TILTC(I) ;
*FOX  CRKVE=XL ;
*FOX  CIKVE=ZL ;
*FOX  CRABAMP2=ED(IX)/(ONE+DPDA) ;
!       call dapri(EJF1,234)
!       write(*,*) crabamp, EJF1, EJF0,clight, "HELLO"
        crabfreq=ek(ix)*c1e3 !JBG Input in MHz changed to kHz
        crabpht2=crabph2(ix)
*FOX  Y(1)=Y(1) + (CRABAMP2*CRKVE)*
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT2);
*FOX  Y(2)=Y(2) - (CRABAMP2*CIKVE)*
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT2);
*FOX  DPDA1=DPDA1 - (1/2.)*(CRABAMP2)*(CRKVE*CRKVE-
*FOX  CIKVE*CIKVE)*(((CRABFREQ*2D0)*PI)/CLIGHT)*C1M3*
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI+CRABPHT2) ;
*FOX  EJF0=EJF1 ;
*FOX  DPDA=DPDA1*C1M3 ;
*FOX  EJF1=E0F*(ONE+DPDA) ;
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;
*FOX  Y(1)=EJF0/EJF1*Y(1) ;
*FOX  Y(2)=EJF0/EJF1*Y(2) ;
          goto 440
      endif
          if(kzz.eq.-26) then
            ! JBG bypass this element if 4D/5D case
            if(iclo6.eq.0) then
!                write(*,*)'Bypassing RF mult 4D or 5D case'
                goto 440
            endif
          xs=xsi(i) ! JBG change of variables for misal calculations
          zs=zsi(i)
*FOX  CRABAMP2=ED(IX)/(ONE+DPDA) ;
             crabfreq=ek(ix)*c1e3
             crabpht2=crabph2(ix)
*FOX  Y(2)=Y(2) + (CRABAMP2*CRKVE)*
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT2);
*FOX  Y(1)=Y(1) + (CRABAMP2*CIKVE)*
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT2);
*FOX  DPDA1=DPDA1 - (CRABAMP2)*(CIKVE*CRKVE)
*FOX  *(((CRABFREQ*2D0)*PI)/CLIGHT)*C1M3*
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI+CRABPHT2) ;
*FOX  EJF0=EJF1 ;
*FOX  DPDA=DPDA1*C1M3 ;
*FOX  EJF1=E0F*(ONE+DPDA) ;
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;
*FOX  Y(1)=EJF0/EJF1*Y(1) ;
*FOX  Y(2)=EJF0/EJF1*Y(2) ;
          endif
          if(kzz.eq.27) then
            ! JBG bypass this element if 4D/5D case
            if(iclo6.eq.0) then
!                write(*,*)'Bypassing RF mult 4D or 5D case'
                goto 440
            endif
          xs=xsi(i)
          zs=zsi(i)
*FOX  XL=(X(1)-XS)*TILTC(I)+(X(2)-ZS)*TILTS(I) ;
*FOX  ZL=-(X(1)-XS)*TILTS(I)+(X(2)-ZS)*TILTC(I) ;
*FOX  CRKVE=XL ;
*FOX  CIKVE=ZL ;
*FOX  CRABAMP3=ED(IX)/(ONE+DPDA) ;
             crabfreq=ek(ix)*c1e3
             crabpht3=crabph3(ix)
*FOX  Y(1)=Y(1) + 2*(1/2.)*CRABAMP3*((CRKVE*CRKVE)-
*FOX  (CIKVE*CIKVE))*C1M3*
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT3);
*FOX  Y(2)=Y(2) - 2*CRABAMP3*(CRKVE*CIKVE)*C1M3*
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT3);
*FOX  DPDA1=DPDA1 - 2*(1/6.)*(CRABAMP3)*(CRKVE*CRKVE*CRKVE-
*FOX  3*CIKVE*CIKVE*CRKVE)*(((CRABFREQ*2D0)*PI)/CLIGHT)*
*FOX  C1M6*
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI+CRABPHT3) ;
*FOX  EJF0=EJF1 ;
*FOX  DPDA=DPDA1*C1M3 ;
*FOX  EJF1=E0F*(ONE+DPDA) ;
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;
*FOX  Y(1)=EJF0/EJF1*Y(1) ;
*FOX  Y(2)=EJF0/EJF1*Y(2) ;
          goto 440
          endif
          if(kzz.eq.-27) then
            ! JBG bypass this element if 4D/5D case
            if(iclo6.eq.0) then
!                write(*,*)'Bypassing RF mult 4D or 5D case'
                goto 440
            endif
          xs=xsi(i)
          zs=zsi(i)
*FOX  XL=(X(1)-XS)*TILTC(I)+(X(2)-ZS)*TILTS(I) ;
*FOX  ZL=-(X(1)-XS)*TILTS(I)+(X(2)-ZS)*TILTC(I) ;
*FOX  CRKVE=XL ;
*FOX  CIKVE=ZL ;
*FOX  CRABAMP3=ED(IX)/(ONE+DPDA) ;
             crabfreq=ek(ix)*c1e3
             crabpht3=crabph3(ix)
*FOX  Y(2)=Y(2) - 2*(1/2.)*CRABAMP3*((CIKVE*CIKVE)-
*FOX  (CRKVE*CRKVE))*C1M3*
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT3);
*FOX  Y(1)=Y(1) + 2*CRABAMP3*(CRKVE*CIKVE)*C1M3*
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT3);
*FOX  DPDA1=DPDA1 + 2*(1/6.)*(CRABAMP3)*(CIKVE*CIKVE*CIKVE-
*FOX  3*CIKVE*CRKVE*CRKVE)*(((CRABFREQ*2D0)*PI)/CLIGHT)*
*FOX  C1M6*
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI+CRABPHT3) ;
*FOX  EJF0=EJF1 ;
*FOX  DPDA=DPDA1*C1M3 ;
*FOX  EJF1=E0F*(ONE+DPDA) ;
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;
*FOX  Y(1)=EJF0/EJF1*Y(1) ;
*FOX  Y(2)=EJF0/EJF1*Y(2) ;
          endif
          if(kzz.eq.28) then
            ! JBG bypass this element if 4D/5D case
            if(iclo6.eq.0) then
!                write(*,*)'Bypassing RF mult 4D or 5D case'
                goto 440
            endif
          xs=xsi(i)
          zs=zsi(i)
*FOX  XL=(X(1)-XS)*TILTC(I)+(X(2)-ZS)*TILTS(I) ;
*FOX  ZL=-(X(1)-XS)*TILTS(I)+(X(2)-ZS)*TILTC(I) ;
*FOX  CRKVE=XL ;
*FOX  CIKVE=ZL ;
*FOX  CRABAMP4=ED(IX)/(ONE+DPDA) ;
             crabfreq=ek(ix)*c1e3
             crabpht4=crabph4(ix)
*FOX  Y(1)=Y(1) + 6*(1/6.)*(CRABAMP4)*
*FOX  (CRKVE*CRKVE*CRKVE-(3*CRKVE*CIKVE*CIKVE))*C1M3*C1M3*
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT4);
*FOX  Y(2)=Y(2) - 6*(1/6.)*(CRABAMP4)*
*FOX  (3*CIKVE*CRKVE*CRKVE-CIKVE*CIKVE*CIKVE)*C1M3*C1M3*
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT4);
*FOX  DPDA1=DPDA1-6*(1/24.)*(CRABAMP4)*(CRKVE*CRKVE*CRKVE*CRKVE-
*FOX  6*CRKVE*CRKVE*CIKVE*CIKVE+CIKVE*CIKVE*CIKVE*CIKVE)*
*FOX  C1M9*(((CRABFREQ*2D0)*PI)/CLIGHT)*
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI+CRABPHT4) ;
*FOX  EJF0=EJF1 ;
*FOX  DPDA=DPDA1*C1M3 ;
*FOX  EJF1=E0F*(ONE+DPDA) ;
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;
*FOX  Y(1)=EJF0/EJF1*Y(1) ;
*FOX  Y(2)=EJF0/EJF1*Y(2) ;
          goto 440
          endif
          if(kzz.eq.-28) then
            ! JBG bypass this element if 4D/5D case
            if(iclo6.eq.0) then
!                write(*,*)'Bypassing RF mult 4D or 5D case'
                goto 440
            endif
          xs=xsi(i)
          zs=zsi(i)
*FOX  XL=(X(1)-XS)*TILTC(I)+(X(2)-ZS)*TILTS(I) ;
*FOX  ZL=-(X(1)-XS)*TILTS(I)+(X(2)-ZS)*TILTC(I) ;
*FOX  CRKVE=XL ;
*FOX  CIKVE=ZL ;
*FOX  CRABAMP4=ED(IX)/(ONE+DPDA) ;
             crabfreq=ek(ix)*c1e3
             crabpht4=crabph4(ix)
*FOX  Y(1)=Y(1) + 6*(1/6.)*(CRABAMP4)*
*FOX  (CIKVE*CIKVE*CIKVE-(3*CIKVE*CRKVE*CRKVE))*C1M3*C1M3*
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT4);
*FOX  Y(2)=Y(2) + 6*(1/6.)*(CRABAMP4)*
*FOX  (3*CRKVE*CIKVE*CIKVE-CRKVE*CRKVE*CRKVE)*C1M3*C1M3*
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT4);
*FOX  DPDA1=DPDA1+6*(1/6.)*(CRABAMP4)*(CRKVE*CRKVE*CRKVE*CIKVE-
*FOX  CIKVE*CIKVE*CIKVE*CRKVE)*
*FOX  C1M9*(((CRABFREQ*2D0)*PI)/CLIGHT)*
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI+CRABPHT4) ;
*FOX  EJF0=EJF1 ;
*FOX  DPDA=DPDA1*C1M3 ;
*FOX  EJF1=E0F*(ONE+DPDA) ;
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;
*FOX  Y(1)=EJF0/EJF1*Y(1) ;
*FOX  Y(2)=EJF0/EJF1*Y(2) ;
          endif
        if(kzz.eq.22) then
          irrtr=imtr(ix)
*FOX  SIGMDA=SIGMDA+COTR(IRRTR,5)+RRTR(IRRTR,5,1)*X(1)+
*FOX  RRTR(IRRTR,5,2)*Y(1)+RRTR(IRRTR,5,3)*X(2)+RRTR(IRRTR,5,4)*Y(2) ;
*FOX  PUX=X(1) ;
*FOX  PUZ=Y(1) ;
*FOX  X(1)=COTR(IRRTR,1)+RRTR(IRRTR,1,1)*PUX+RRTR(IRRTR,1,2)*PUZ+
*FOX  RRTR(IRRTR,1,6)*DPDA1 ;
*FOX  Y(1)=COTR(IRRTR,2)+RRTR(IRRTR,2,1)*PUX+RRTR(IRRTR,2,2)*PUZ+
*FOX  RRTR(IRRTR,2,6)*DPDA1 ;
*FOX  PUX=X(2) ;
*FOX  PUZ=Y(2) ;
*FOX  X(2)=COTR(IRRTR,3)+RRTR(IRRTR,3,3)*PUX+RRTR(IRRTR,3,4)*PUZ+
*FOX  RRTR(IRRTR,3,6)*DPDA1 ;
*FOX  Y(2)=COTR(IRRTR,4)+RRTR(IRRTR,4,3)*PUX+RRTR(IRRTR,4,4)*PUZ+
*FOX  RRTR(IRRTR,4,6)*DPDA1 ;
      endif
        if(kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22) goto 440
        ipch=0
        if(iqmodc.eq.1) then
          if(ix.eq.iq(1).or.iratioe(ix).eq.iq(1)) then
            ipch=1
          else if(ix.eq.iq(2).or.iratioe(ix).eq.iq(2)) then
            ipch=2
          endif
        endif
        if(ichromc.eq.1) then
          if(ix.eq.issss(1).or.iratioe(ix).eq.issss(1)) then
            ipch=1
          else if(ix.eq.issss(2).or.iratioe(ix).eq.issss(2)) then
            ipch=2
          endif
        endif
        if(ipch.ne.0) then
*FOX  EKK=(SMIDA(IPCH)*RATIOE(IX)+SMIZF(I))/(ONE+DPDA) ;
        else
*FOX  EKK=SMI(I)/(ONE+DPDA) ;
        endif
        xs=xsi(i)
        zs=zsi(i)
*FOX  XL=(X(1)-XS)*TILTC(I)+(X(2)-ZS)*TILTS(I) ;
*FOX  ZL=-(X(1)-XS)*TILTS(I)+(X(2)-ZS)*TILTC(I) ;
*FOX  CRKVE=XL ;
*FOX  CIKVE=ZL ;
        if(kzz.lt.0) goto 320
        goto(90,100,110,120,130,140,150,160,170,180,190,440,440,440,    &
     &       440,440,440,440,440,440,440,440,440,185,186),kzz
        goto 440
!--HORIZONTAL DIPOLE
   90   continue
*FOX  EKK=EKK*C1E3 ;
*FOX  Y(1)=Y(1)+EKK*TILTC(I) ;
*FOX  Y(2)=Y(2)+EKK*TILTS(I) ;
        goto 440
!--NORMAL QUADRUPOLE
  100   continue
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;
        goto 440
!---NORMAL SEXTUPOLE
  110   continue
*FOX  EKK=EKK*C1M3 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;
        goto 440
!--NORMAL OCTUPOLE
  120   continue
*FOX  EKK=EKK*C1M6 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;
        goto 440
!--NORMAL DECAPOLE
  130   continue
*FOX  EKK=EKK*C1M9 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;
        goto 440
!---NORMAL DODECAPOL
  140   continue
*FOX  EKK=EKK*C1M12 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;
        goto 440
!---NORMAL 14-POL
  150   continue
*FOX  EKK=EKK*C1M15 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;
        goto 440
!---NORMAL 16-POL
  160   continue
*FOX  EKK=EKK*C1M18 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;
        goto 440
!---NORMAL 18-POL
  170   continue
*FOX  EKK=EKK*C1M21 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;
        goto 440
!---NORMAL 20-POL
  180   continue
*FOX  EKK=EKK*C1M24 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;
        goto 440
!--DIPEDGE ELEMENT
  185   continue
*FOX  Y(1)=Y(1)+(ED(IX)*TILTC(I)*CRKVE-EK(IX)*TILTS(I)*CIKVE)/
*FOX  (ONE+DPDA) ;
*FOX  Y(2)=Y(2)+(EK(IX)*TILTC(I)*CIKVE+ED(IX)*TILTS(I)*CRKVE)/
*FOX  (ONE+DPDA) ;
        goto 440
!--solenoid
  186   continue
*FOX  Y(1)=Y(1)-X(2)*ED(IX) ;
*FOX  Y(2)=Y(2)+X(1)*ED(IX) ;
*FOX  CRKVE=Y(1)-X(1)*ED(IX)*EK(IX) ;
*FOX  CIKVE=Y(2)-X(2)*ED(IX)*EK(IX) ;
*FOX  Y(1)=CRKVE*COS(EK(IX))+CIKVE*SIN(EK(IX)) ;
*FOX  Y(2)=-CRKVE*SIN(EK(IX))+CIKVE*COS(EK(IX)) ;
*FOX  CRKVE=X(1)*COS(EK(IX))+X(2)*SIN(EK(IX)) ;
*FOX  CIKVE=-X(1)*SIN(EK(IX))+X(2)*COS(EK(IX)) ;
*FOX  X(1)=CRKVE ;
*FOX  X(2)=CIKVE ;
*FOX  Y(1)=Y(1)+X(2)*ED(IX) ;
*FOX  Y(2)=Y(2)-X(1)*ED(IX) ;
        goto 440
  190   r0=ek(ix)
        nmz=nmu(ix)
          if(abs(dki(ix,1)).gt.pieni) then
            if(abs(dki(ix,3)).gt.pieni) then
*FOX  DKIP=DKI(IX,1)/(ONE+DPDA) ;
*FOX  Y(1)=Y(1)-(DKI(IX,1)/DKI(IX,3)*XL+DPDA*C1E3)*
*FOX  TILTC(I)*DKIP
*FOX  +C1E3*DKI(IX,1)/(ONE+DPDA)*(ONE-TILTC(I)) ;
*FOX  Y(2)=Y(2)-(DKI(IX,1)/DKI(IX,3)*XL+DPDA*C1E3)*
*FOX  TILTS(I)*DKIP
*FOX  +C1E3*DKI(IX,1)/(ONE+DPDA)*TILTS(I) ;
            else
*FOX  Y(1)=Y(1)-DKI(IX,1)*DPDA*C1E3/(ONE+DPDA)*TILTC(I)
*FOX  +C1E3*DKI(IX,1)/(ONE+DPDA)*(ONE-TILTC(I)) ;
*FOX  Y(2)=Y(2)-DKI(IX,1)*DPDA*C1E3/(ONE+DPDA)*TILTS(I)
*FOX  +C1E3*DKI(IX,1)/(ONE+DPDA)*TILTS(I) ;
            endif
            if(idp.eq.1.and.iabs(ition).eq.1) then
*FOX  SIGMDA=SIGMDA+RV*DKI(IX,1)*XL ;
            endif
          endif
          if(abs(dki(ix,2)).gt.pieni) then
            if(abs(dki(ix,3)).gt.pieni) then
*FOX  DKIP=DKI(IX,2)/(ONE+DPDA) ;
*FOX  Y(1)=Y(1)+(DKI(IX,2)/DKI(IX,3)*ZL-DPDA*C1E3)*
*FOX  TILTS(I)*DKIP
*FOX  +C1E3*DKI(IX,2)/(ONE+DPDA)*TILTS(I) ;
*FOX  Y(2)=Y(2)-(DKI(IX,2)/DKI(IX,3)*ZL-DPDA*C1E3)*
*FOX  TILTC(I)*DKIP
*FOX  -C1E3*DKI(IX,2)/(ONE+DPDA)*(ONE-TILTC(I)) ;
            else
*FOX  Y(1)=Y(1)-DKI(IX,2)*DPDA*C1E3/(ONE+DPDA)*TILTS(I)
*FOX  +C1E3*DKI(IX,2)/(ONE+DPDA)*TILTS(I) ;
*FOX  Y(2)=Y(2)+DKI(IX,2)*DPDA*C1E3/(ONE+DPDA)*TILTC(I)
*FOX  -C1E3*DKI(IX,2)/(ONE+DPDA)*(ONE-TILTC(I)) ;
            endif
            if(idp.eq.1.and.iabs(ition).eq.1) then
*FOX  SIGMDA=SIGMDA-RV*DKI(IX,2)*ZL ;
            endif
          endif
        if(abs(r0).le.pieni.or.nmz.eq.0) goto 440
        if(nmz.ge.2) then
*FOX  YV1J=BBI(I,1)+BBI(I,2)*XL+AAI(I,2)*ZL ;
*FOX  YV2J=AAI(I,1)-BBI(I,2)*ZL+AAI(I,2)*XL ;
*FOX  CRKVE=XL ;
*FOX  CIKVE=ZL ;
          do 200 k=3,nmz
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  YV1J=YV1J+BBI(I,K)*CRKVE+AAI(I,K)*CIKVE ;
*FOX  YV2J=YV2J-BBI(I,K)*CIKVE+AAI(I,K)*CRKVE ;
  200     continue
*FOX  Y(1)=Y(1)+(TILTC(I)*YV1J-TILTS(I)*YV2J)/(ONE+DPDA) ;
*FOX  Y(2)=Y(2)+(TILTC(I)*YV2J+TILTS(I)*YV1J)/(ONE+DPDA) ;
        else
*FOX  Y(1)=Y(1)+(TILTC(I)*BBI(I,1)-TILTS(I)*AAI(I,1))/(ONE+DPDA) ;
*FOX  Y(2)=Y(2)+(TILTC(I)*AAI(I,1)+TILTS(I)*BBI(I,1))/(ONE+DPDA) ;
        endif
        goto 440
!--SKEW ELEMENTS
  320   kzz=-kzz
        goto(330,340,350,360,370,380,390,400,410,420),kzz
        goto 440
!---VERTICAL DIPOLE
  330   continue
*FOX  EKK=EKK*C1E3 ;
*FOX  Y(1)=Y(1)-EKK*TILTS(I) ;
*FOX  Y(2)=Y(2)+EKK*TILTC(I) ;
        goto 440
!---SKEW QUADRUPOLE
  340   continue
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
        goto 440
!---SKEW SEXTUPOLE
  350   continue
*FOX  EKK=EKK*C1M3 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
        goto 440
!---SKEW OCTUPOLE
  360   continue
*FOX  EKK=EKK*C1M6 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
        goto 440
!---SKEW DECAPOLE
  370   continue
*FOX  EKK=EKK*C1M9 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
        goto 440
!---SKEW DODECAPOL
  380   continue
*FOX  EKK=EKK*C1M12 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
        goto 440
!---SKEW 14-POL
  390   continue
*FOX  EKK=EKK*C1M15 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
        goto 440
!---SKEW 16-POL
  400   continue
*FOX  EKK=EKK*C1M18 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
        goto 440
!---SKEW 18-POL
  410   continue
*FOX  EKK=EKK*C1M21 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
        goto 440
!---SKEW 20-POL
  420   continue
*FOX  EKK=EKK*C1M24 ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;
*FOX  CRKVE=CRKVEUK ;
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;
 440  continue
      if(ilinc.eq.1) then
        typ=bez(ix)
          iwrite=0
          if(nlin.eq.0) then
            iwrite=1
          else
            do ii=1,nlin
              if(typ.eq.bezl(ii)) iwrite=1
            enddo
          endif
*FOX  YP(1)=Y(1)*(ONE+DPDA) ;
*FOX  YP(2)=Y(2)*(ONE+DPDA) ;
*FOX  DPDA1=DPDA*C1E3 ;
          call dacop(x(1),damap(1))
          call dacop(yp(1),damap(2))
          call dacop(x(2),damap(3))
          call dacop(yp(2),damap(4))
          do j=1,2
            ii=2*j
            call dapek(damap(ii-1),jj,c(j))
            call dapek(damap(ii),jj,cp(j))
          enddo
          call dacsu(damap(1),c(1),damap(1))
          call dacsu(damap(2),cp(1),damap(2))
          call dacsu(damap(3),c(2),damap(3))
          call dacsu(damap(4),cp(2),damap(4))
          if(ndimf.eq.3) then
            call dacop(sigmda,damap(5))
            call dacop(dpda1,damap(6))
            call dapek(damap(5),jj,c(3))
            call dapek(damap(6),jj,cp(3))
            call dacsu(damap(5),c(3),damap(5))
            call dacsu(damap(6),cp(3),damap(6))
            if(iflag1.eq.1.and.ithick.eq.1) then
              call dacct(damap,nvar,corrnew,nvar,damap,nvar)
            endif
          else
            call dacop(dpda1,damap(5))
            do j1=1,4
              do ii=1,4
                jj(ii)=1
                call dapek(damap(j1),jj,rdd(j1,ii))
                jj(ii)=0
              enddo
            enddo
            jj(5)=1
            do j1=1,4
              call dapek(damap(j1),jj,rdd(j1,5))
            enddo
            jj(5)=0
            do j1=1,2
              ii=2*j1
!hr03         d(j1)=rdd(ii-1,1)*dicu(1)+rdd(ii-1,2)*dicu(2)+            &
!hr03&rdd(ii-1,3)*dicu(3)+rdd(ii-1,4)*dicu(4)+rdd(ii-1,5)
              d(j1)=(((rdd(ii-1,1)*dicu(1)+rdd(ii-1,2)*dicu(2))+        &!hr03
     &rdd(ii-1,3)*dicu(3))+rdd(ii-1,4)*dicu(4))+rdd(ii-1,5)              !hr03
!hr03         dp(j1)=rdd(ii,1)*dicu(1)+rdd(ii,2)*dicu(2)+               &
!hr03&rdd(ii,3)*dicu(3)+rdd(ii,4)*dicu(4)+rdd(ii,5)
              dp(j1)=(((rdd(ii,1)*dicu(1)+rdd(ii,2)*dicu(2))+           &!hr03
     &rdd(ii,3)*dicu(3))+rdd(ii,4)*dicu(4))+rdd(ii,5)                    !hr03
            enddo
          endif
          call dacct(damap,nvar,aa2,nvar,damap,nvar)
          do j=1,ndimf
            ii=2*j
            if(j.eq.1) then
              i2=4
              i3=6
            elseif(j.eq.2) then
              i2=2
              i3=6
            elseif(j.eq.3) then
              i2=2
              i3=4
            endif
            jj(ii-1)=1
            call dapek(damap(ii-1),jj,angp(1,ii-1))
            call dapek(damap(ii),jj,au(ii,ii-1))
            jj(ii-1)=0
            jj(ii)=1
            call dapek(damap(ii-1),jj,angp(1,ii))
            call dapek(damap(ii),jj,au(ii,ii))
            jj(ii)=0
            jj(i2-1)=1
            call dapek(damap(ii-1),jj,au(i2-1,i2-1))
            call dapek(damap(ii),jj,au(i2,i2-1))
            jj(i2-1)=0
            jj(i2)=1
            call dapek(damap(ii-1),jj,au(i2-1,i2))
            call dapek(damap(ii),jj,au(i2,i2))
            jj(i2)=0
            jj(i3-1)=1
            call dapek(damap(ii-1),jj,au(i3-1,i3-1))
            call dapek(damap(ii),jj,au(i3,i3-1))
            jj(i3-1)=0
            jj(i3)=1
            call dapek(damap(ii-1),jj,au(i3-1,i3))
            call dapek(damap(ii),jj,au(i3,i3))
            jj(i3)=0
!hr08       b1(j)=angp(1,ii-1)*angp(1,ii-1)+angp(1,ii)*angp(1,ii)
            b1(j)=angp(1,ii-1)**2+angp(1,ii)**2                          !hr08
!hr08       b2(j)=au(i2-1,i2-1)*au(i2-1,i2-1)+au(i2-1,i2)*au(i2-1,i2)
            b2(j)=au(i2-1,i2-1)**2+au(i2-1,i2)**2                        !hr08
!hr08       b3(j)=au(i3-1,i3-1)*au(i3-1,i3-1)+au(i3-1,i3)*au(i3-1,i3)
            b3(j)=au(i3-1,i3-1)**2+au(i3-1,i3)**2                        !hr08
!hr03       al1(j)=-(angp(1,ii-1)*au(ii,ii-1)+angp(1,ii)*au(ii,ii))
            al1(j)=-1d0*(angp(1,ii-1)*au(ii,ii-1)+angp(1,ii)*au(ii,ii))  !hr03
!hr03       al2(j)=-(au(i2-1,i2-1)*au(i2,i2-1)+au(i2-1,i2)*au(i2,i2))
           al2(j)=-1d0*(au(i2-1,i2-1)*au(i2,i2-1)+au(i2-1,i2)*au(i2,i2)) !hr03
!hr03       al3(j)=-(au(i3-1,i3-1)*au(i3,i3-1)+au(i3-1,i3)*au(i3,i3))
           al3(j)=-1d0*(au(i3-1,i3-1)*au(i3,i3-1)+au(i3-1,i3)*au(i3,i3)) !hr03
!hr04       g1(j)=au(ii,ii-1)*au(ii,ii-1)+au(ii,ii)*au(ii,ii)
!hr04       g2(j)=au(i2,i2-1)*au(i2,i2-1)+au(i2,i2)*au(i2,i2)
!hr04       g3(j)=au(i3,i3-1)*au(i3,i3-1)+au(i3,i3)*au(i3,i3)
            g1(j)=au(ii,ii-1)**2+au(ii,ii)**2                            !hr04
            g2(j)=au(i2,i2-1)**2+au(i2,i2)**2                            !hr04
            g3(j)=au(i3,i3-1)**2+au(i3,i3)**2                            !hr04
            if(ndimf.eq.3) then
              call dainv(damap,nvar,damapi,nvar)
              jj(6)=1
              call dapek(damapi(5),jj,aui(1))
              call dapek(damapi(6),jj,aui(2))
              jj(6)=0
              if(j.lt.3) then
                d(j)=au(i3-1,i3-1)*aui(1)+au(i3-1,i3)*aui(2)
                dp(j)=au(i3,i3-1)*aui(1)+au(i3,i3)*aui(2)
              else
                d(j)=angp(1,ii-1)*aui(1)+angp(1,ii)*aui(2)
                dp(j)=au(ii,ii-1)*aui(1)+au(ii,ii)*aui(2)
              endif
            endif
            sx=angp(2,ii-1)*angp(1,ii)-angp(1,ii-1)*angp(2,ii)
            cx=angp(1,ii-1)*angp(2,ii-1)+angp(1,ii)*angp(2,ii)
            if(abs(sx).gt.c1m15.or.abs(cx).gt.c1m15) then
              dphi(j)=atan2_rn(sx,cx)/x2pi
            else
              dphi(j)=zero
            endif
            phi(j)=phi(j)+dphi(j)
          enddo
          do j=1,ndimf
            ii=2*j
            angp(2,ii-1)=angp(1,ii-1)
            angp(2,ii)=angp(1,ii)
          enddo
          if(iwrite.eq.1) then
            iii=i
            if(typ(:8).eq.'START   ') iii=0
            write(*,10030) iii,typ(:8),tl,phi(1),b1(1),al1(1),g1(1),    &
     &d(1),dp(1),c(1),cp(1)
            if(ndimf.eq.3) then
              write(*,10040) b2(1),al2(1),g2(1)
              write(*,10050) typ(9:16),b3(1),al3(1),g3(1)
            else
              write(*,10055) typ(9:16),b2(1),al2(1),g2(1)
            endif
            write(*,10060)
            write(*,10070) phi(2),b1(2),al1(2),g1(2),d(2),dp(2),c(2),   &
     &cp(2)
            write(*,10080) b2(2),al2(2),g2(2)
            if(ndimf.eq.3) then
              write(*,10090) b3(2),al3(2),g3(2)
              write(*,10060)
              write(*,10100) -phi(3),b1(3),al1(3),g1(3),d(3),dp(3),c(3),&
     &cp(3)*c1m3
              write(*,10080) b2(3),al2(3),g2(3)
              write(*,10040) b3(3),al3(3),g3(3)
            endif
            write(*,10010)
          endif
        if(i.eq.nt) goto 470
      endif
 430  continue
*FOX  YP(1)=Y(1)*(ONE+DPDA) ;
*FOX  YP(2)=Y(2)*(ONE+DPDA) ;
      if(icav.eq.0.or.ithick.ne.1) then
        if(icoonly.eq.1) then
          xxtr(1,1) = dare(x(1))
          yytr(1,1) = dare(y(1))
          xxtr(1,2) = dare(x(2))
          yytr(1,2) = dare(y(2))
          sigm(1) = dare(sigmda)
          dps(1) = dare(dpda)
        endif
      else
*FOX  CORRAU1(1)=X(1) ;
*FOX  CORRAU1(2)=YP(1) ;
*FOX  CORRAU1(3)=X(2) ;
*FOX  CORRAU1(4)=YP(2) ;
*FOX  CORRAU1(5)=SIGMDA ;
*FOX  CORRAU1(6)=DPDA1 ;
        do 435 kkk=1,6
          dpdav2(kkk)=dare(corrau1(kkk))
*FOX  CORRAU1(KKK)=CORRAU1(KKK)-DPDAV2(KKK) ;
  435   continue
        if(ivar.gt.ivar1) then
*FOX  CORRAU1(7)=SMIDA(1) ;
*FOX  CORRAU1(8)=SMIDA(2) ;
          dpdav=dare(smida(1))
*FOX  CORRNEW(7)=SMIDA(1)-DPDAV ;
          dpdav=dare(smida(2))
*FOX  CORRNEW(8)=SMIDA(2)-DPDAV ;
        endif
        call dacct(corrau1,nvar,corrnew,nvar,corrau2,nvar)
        do 436 kkk=1,6
*FOX  CORRAU2(KKK)=CORRAU2(KKK)+DPDAV2(KKK) ;
  436   continue
*FOX  CORRAU1(2)=CORRAU2(2)/(ONE+CORRAU2(6)) ;
*FOX  CORRAU1(4)=CORRAU2(4)/(ONE+CORRAU2(6)) ;
*FOX  X(1)=CORRAU2(1) ;
*FOX  YP(1)=CORRAU2(2) ;
*FOX  X(2)=CORRAU2(3) ;
*FOX  YP(2)=CORRAU2(4) ;
*FOX  SIGMDA=CORRAU2(5) ;
*FOX  DPDA1=CORRAU2(6) ;
*FOX  DPDA=DPDA1*C1M3 ;
*FOX  Y(1)=YP(1)/(ONE+DPDA) ;
*FOX  Y(2)=YP(2)/(ONE+DPDA) ;
        if(icoonly.eq.1) then
          xxtr(1,1) = dare(x(1))
          yytr(1,1) = dare(y(1))
          xxtr(1,2) = dare(x(2))
          yytr(1,2) = dare(y(2))
          sigm(1) = dare(sigmda)
          dps(1) = dare(dpda)
        endif
      endif
      call dacop(x(1),damap(1))
      call dacop(x(2),damap(3))
      if(ndimf.eq.3) call dacop(sigmda,damap(5))
      if(icoonly.eq.1.or.iqmodc.eq.3) then
        call dacop(y(1),damap(2))
        call dacop(y(2),damap(4))
        if(ndimf.eq.3) call dacop(dpda,damap(6))
        do i=1,nd2
          jj(i)=1
          do ii=1,nd2
            call dapek(damap(ii),jj,aml6(ii,i))
            if(i.eq.6) aml6(ii,i)=aml6(ii,i)*c1e3
          enddo
          jj(i)=0
        enddo
        do i=1,nd2
          aml6(i,i)=aml6(i,i)-one
        enddo
      endif
      call dacop(yp(1),damap(2))
      call dacop(yp(2),damap(4))
      if(ndimf.eq.3) then
        call dacop(sigmda,damap(5))
        call dacop(dpda1,damap(6))
      else
        call dacop(dpda1,damap(5))
      endif
      if(iqmodc.eq.2.or.iqmodc.eq.4.or.ilin.ge.2) then
        rewind 18
!Eric
        rewind 111
        call daprid(damap,1,nvar,18)
      endif
!--now do the output
      if(iqmodc.eq.1) call danot(3)
      if(iqmodc.eq.3) call danot(2)
      if(ichromc.eq.1) call danot(4)
      if(ilinc.eq.1.or.ilinc.eq.2.or.iqmodc.eq.1.or.iqmodc.eq.3.or.     &
     &ichromc.eq.1) then
        call mapnorm(damap,f,aa2,a1,xy,h,nord1)
      endif
      if(iqmodc.eq.1.or.iqmodc.eq.3) then
        call gettura(wxys,rrad)
        wxys(3)=abs(wxys(3))
        write(*,*) (wxys(i),i=1,ndimf)
        do i=1,nd2
          jj(i)=1
          do ii=1,nd2
            call dapek(aa2(ii),jj,tas)
            if(i.eq.6.and.ii.ne.6) tas=tas*c1e3
            if(ii.eq.6.and.i.ne.6) tas=tas*c1m3
            tasm(ii,i)=tas
          enddo
          jj(i)=0
        enddo
      endif
      if(iqmodc.eq.1) then
        call dhdj(h,df)
        do i=1,ndimf
          call dapek(df(ndimf+i),jj,corr(1,i))
        enddo
        corr(1,3)=abs(corr(1,3))
        jj(nd2+1)=1
        call dapek(df(ndimf+1),jj,coefh1)
        call dapek(df(ndimf+2),jj,coefv1)
        jj(nd2+1)=0
        jj(nd2+2)=1
        call dapek(df(ndimf+1),jj,coefh2)
        call dapek(df(ndimf+2),jj,coefv2)
        jj(nd2+2)=0
        det1=coefh1*coefv2-coefv1*coefh2
        if(abs(det1).le.pieni) call prror(90)
        corr(2,1)=coefv2/det1
!hr05   corr(2,2)=-coefh2/det1
        corr(2,2)=(-1d0*coefh2)/det1                                     !hr05
!hr05   corr(3,1)=-coefv1/det1
        corr(3,1)=(-1d0*coefv1)/det1                                     !hr05
        corr(3,2)=coefh1/det1
      endif
      if(ichromc.eq.1) then
        call dhdj(h,df)
        jj(nd2+1)=1
        call dapek(df(ndimf+1),jj,corr(1,1))
        call dapek(df(ndimf+2),jj,corr(1,2))
        jj(nd2+2)=1
        call dapek(df(ndimf+1),jj,coefh1)
        call dapek(df(ndimf+2),jj,coefv1)
        jj(nd2+2)=0
        jj(nd2+3)=1
        call dapek(df(ndimf+1),jj,coefh2)
        call dapek(df(ndimf+2),jj,coefv2)
        jj(nd2+3)=0
        jj(nd2+1)=0
        det1=coefh1*coefv2-coefv1*coefh2
        if(abs(det1).le.pieni) call prror(96)
        corr(2,1)=coefv2/det1
!hr05   corr(2,2)=-coefh2/det1
        corr(2,2)=(-1d0*coefh2)/det1                                     !hr05
!hr05   corr(3,1)=-coefv1/det1
        corr(3,1)=(-1d0*coefv1)/det1                                     !hr05
        corr(3,2)=coefh1/det1
      endif
 470  continue
      call dadal(damap,6)
      call dadal(damapi,6)
      call dadal(damap1,6)
      call dadal(angno,1)
      call dadal(f,1)
      call dadal(aa2,6)
      call dadal(aa2r,6)
      call dadal(a1,6)
      call dadal(a1r,6)
      call dadal(xy,6)
      call dadal(h,1)
      call dadal(df,6)
!     DADAL AUTOMATIC INCLUSION
      return
!-----------------------------------------------------------------------
10000 format(/t5 ,'---- ENTRY ',i1,'D LINOPT ----')
10010 format(132('-'))
10020 format('  NR     TYP      L-TOTAL    P     PHI          ',        &
     &'BETA         ALFA         GAMMA        DIS        DISP         ',&
     &'CLO        CLOP'/ 1x,                                            &
     &'                    (M)           (2*PI)        ',               &
     &'(M)          (RAD)         (M)         (M)        (RAD)        ',&
     &'(MM)       (MRAD)')
10030 format('|',i6,'|',a8,'|',f12.5,'|','X','|',f12.7,'|',f12.6,'|',   &
     &f13.7,'|',f11.6,'|',f11.7,'|',f11.7,'|',f11.7,'|',f11.7,'|')
10040 format('|',6x,'|',8x,'|',12x,'|','Y','|',12x,'|',f12.6,'|', f13.7,&
     &'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10050 format('|',6x,'|',a8,'|',12x,'|','S','|',12x,'|',f12.6,'|', f13.  &
     &7,'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10055 format('|',6x,'|',a8,'|',12x,'|','Y','|',12x,'|',f12.6,'|', f13.  &
     &7,'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10060 format('|',6x,'|',8x,'|',12x,'|',102('-'))
10070 format('|',6x,'|',8x,'|',12x,'|','Y','|',f12.7,'|',f12.6,'|', f13.&
     &7,'|',f11.6,'|',f11.7,'|',f11.7,'|',f11.7,'|',f11.7,'|')
10080 format('|',6x,'|',8x,'|',12x,'|','X','|',12x,'|',f12.6,'|', f13.7,&
     &'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10090 format('|',6x,'|',8x,'|',12x,'|','S','|',12x,'|',f12.6,'|', f13.7,&
     &'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10100 format('|',6x,'|',8x,'|',12x,'|','S','|',f12.7,'|',f12.6,'|', f13.&
     &7,'|',f11.6,'|',f11.7,'|',f11.7,'|',f11.7,'|',f11.7,'|')
10110 format(/t10,'CO-TRACKING ENDED ABNORMALLY'/t10, 'PARTICLE NO. '   &
     &,i3,' AT ELEMENT ',i4/ t10,'HORIZ:  AMPLITUDE = ',f15.3,          &
     &'   APERTURE = ',f15.3/ t10,'VERT:   AMPLITUDE = ',f15.3,         &
     &'   APERTURE = ',f15.3/ t10,'ELEMENT - LIST NUMBER ',i4,          &
     &' TYP NUMBER ',i4,' NAME ',a16/)
10120 format(/t10,'CO-TRACKING ENDED ABNORMALLY'/t10, 'PARTICLE NO. '   &
     &,i3,' AT ELEMENT ',i4/ t10,'HORIZ:  AMPLITUDE = ',f15.3,          &
     &'   APERTURE = ',f15.3/ t10,'VERT:   AMPLITUDE = ',f15.3,         &
     &'   APERTURE = ',f15.3/ t10,'ELEMENT - LIST NUMBER ',i4,          &
     &' TYP NUMBER ',i4,' NAME ',a16/)
10130 format('  LINEAR OPTICS CALCULATION WITH PRINTOUT ',              &
     &'AFTER EACH BLOCK'/                                               &
     &'   A T T E N T I O N : BETATRON PHASE CALCULATION MIGHT BE WRONG'&
     &,' BY A MULTIPLE OF 0.5 FOR EACH LARGE BLOCK'/)
      end
      subroutine wireda
!     This program sends a particle with coordinates (x,a,y,b,d)
!       through
!     the map of a straight current wire.
!     a=p_x/p_0, b=p_y/p_0; p_0 is the reference particle's momentum;
!     d - relative momentum deviation.
!     The wire may have arbitrary orientation (tilt), and the map is
!       computed
!     in the thin lens approximation (<=>second order symplectic
!       integration).
!     The current wire is specified by the following parameters:
!       l -length
!       cur - current (positive if it flows in the +z direction)
!       (dx,dy) - horizontal and vertical distance between wire midpoint
!       and closed orbit
!       (tx,ty) - tilt of the wire w.r.t the closed orbit in the
!       horizontal and vertical planes (in degrees)
!       embl - length of the embedding drift. Without loss of generality
!       it is assumed that the wire midpoint and the embedding
!       drift midpoint coincide.
!     Current is in Amperes and all distances are in meters.
      implicit none
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
      integer ix,idaa
      double precision l,cur,dx,dy,tx,ty,embl,leff,rx,ry,lin,chi
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,    &
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,    &
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele=1200,nblo=600,nper=16,nelb=140,nblz=20000,         &
     &nzfz = 300000,mmul = 20)
      parameter(nran = 2000000,ncom = 100,mran = 500,mpa = 6,nrco = 5,  &
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 99)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 350)
      double precision c180e0,c1e1,c1e12,c1e13,c1e15,c1e16,c1e2,c1e3,   &
     &c1e4,c1e6,c1m1,c1m7,c1m10,c1m12,c1m13,c1m15,c1m18,c1m2,c1m21,     &
     &c1m24,c1m3,c1m36,c1m38,c1m6,c1m9,c2e3,c4e3,crade,clight,four,half,&
     &one,pieni,pmae,pmap,three,two,zero
      parameter(pieni = 1d-38)
      parameter(zero = 0.0d0,half = 0.5d0,one = 1.0d0)
      parameter(two = 2.0d0,three = 3.0d0,four = 4.0d0)
      parameter(c1e1 = 1.0d1,c1e2 = 1.0d2,c1m2 = 1.0d-2)
      parameter(c1e3 = 1.0d3,c2e3 = 2.0d3,c4e3 = 4.0d3,c1e4 = 1.0d4)
      parameter(c1e12 = 1.0d12,c1e13 = 1.0d13,c1e15 = 1.0d15,c1e16 =    &
     &1.0d16)
      parameter(c180e0 = 180.0d0,c1e6 = 1.0d6)
      parameter(c1m1 = 1.0d-1,c1m3 = 1.0d-3,c1m6 = 1.0d-6,c1m7 = 1.0d-7)
      parameter(c1m9 = 1.0d-9,c1m10 = 1.0d-10,c1m12 = 1.0d-12)
      parameter(c1m13 = 1.0d-13,c1m15 = 1.0d-15)
      parameter(c1m18 = 1.0d-18,c1m21 = 1.0d-21,c1m24 = 1.0d-24)
      parameter(c1m36 = 1.0d-36,c1m38 = 1.0d-38)
!     electron mass from PDG, 2002
      parameter(pmap = 938.271998d0,pmae = .510998902d0)
      parameter(crade = 2.817940285d-15, clight = 2.99792458d8)
      integer iav,ibb6d,ibbc,ibeco,ibidu,ibtyp,ic,icext,icextal,iclo,   &
     &iclo6,iclo6r,icode,icoe,icomb,icomb0,iconv,icow,icr,idam,idfor,   &
     &idis,idp,ierro,iffw,ifh,iicav,il,ilin,imad,imbb,                  &
     &imc,imtr,iorg,iout,                                               &
     &ipos,ipr,iprint,ipt,iq,iqmod,iqmod6,iratioe,ird,ire,ires,         &
     &irew,irm,irmod2,ise,ise1,ise2,ise3,isea,iskew,iskip,istw,         &
     &isub,itco,itcro,itf,ithick,ition,itionc,itqv,its6d,iu,iver,ivox,  &
     &ivoz,iwg,ixcav,izu0,kanf,kp,kpa,kwtype,kz,lhc,m21,m22,m23,mblo,   &
     &mbloz,mcut,mel,mesa,mmac,mout2,mp,mper,mstr,msym,mtyp,mzu,napx,   &
     &napxo,nbeam,nch,ncororb,ncorrep,ncorru,ncy,ndafi,nde,nhcorr,      &
     &nhmoni,niu,nlin,nmu,npp,nprint,nqc,nre,nrr,nskew,                 &
     &nstart,nstop,nt,nta,ntco,nte,ntwin,nu,numl,numlr,nur,nvcorr,      &
     &nvmoni,nwr, nturn1, nturn2, nturn3, nturn4,numlcp,numlmax,nnuml
      double precision a,ak0,aka,alfx,alfz,amp0,aper,apx,apz,ape,bbcu,  &
     &bclorb,beamoff,benkc,benki,betac,betam,betx,betz,bk0,bka,bl1,bl2, &
     &clo6,clobeam,clop6,cma1,cma2,cotr,crad,de0,dech,ded,dfft,         &
     &di0,dip0,dki,dkq,dma,dmap,dphix,dphiz,dppoff,dpscor,dqq,dres,dsi, &
     &dsm0,dtr,e0,ed,ej,ejf,ek,el,elbe,emitx,emity,emitz,extalign,      &
     &exterr,eui,euii,gammar,hsy,hsyc,pac,pam,parbe,parbe14,partnum,    &
     &phas,phas0,phasc,pi,pi2,pisqrt,pma,ptnfac,qs,qw0,qwsk,qx0,qxt,qz0,&
     &qzt,r00,rad,rat,ratio,ratioe,rrtr,rtc,rts,rvf,                    &
     &sigcor,sige,sigma0,sigman,sigman2,sigmanq,sigmoff,sigz,sm,ta,tam1,&
     &tam2,tiltc,tilts,tlen,totl,track6d,xpl,xrms,zfz,zpl,zrms,wirel,   &
     &acdipph, crabph, bbbx, bbby, bbbs,                                &
     &crabph2, crabph3, crabph4
      real hmal
      character*16 bez,bezb,bezr,erbez,bezl
      character*80 toptit,sixtit,commen
      common/erro/ierro,erbez
      common/kons/pi,pi2,pisqrt,rad
      common/str /il,mper,mblo,mbloz,msym(nper),kanf,iu,ic(nblz)
      common/ell /ed(nele),el(nele),ek(nele),sm(nele),kz(nele),kp(nele)
      common/bbb /bbbx(nele), bbby(nele), bbbs(nele)
      common/pla /xpl(nele),xrms(nele),zpl(nele),zrms(nele)
      common/str2 /mel(nblo),mtyp(nblo,nelb),mstr(nblo)
      common/mat/a(nele,2,6),bl1(nblo,2,6),bl2(nblo,2,6)
      common/syos2/rvf(mpa)
      common/tra1/rat,idfor,napx,napxo,numl,niu(2),numlr,nde(2),nwr(4), &
     &ird,imc,irew,ntwin,iclo6,iclo6r,iver,ibidu,numlcp,numlmax,nnuml
      common/syn/qs,e0,pma,ej(mpa),ejf(mpa),phas0,phas,hsy(3),          &
     &crad,hsyc(nele),phasc(nele),dppoff,sigmoff(nblz),tlen,            &
     &iicav,itionc(nele),ition,idp,ncy,ixcav
      common/corcom/dpscor,sigcor,icode,idam,its6d
      common/multi/bk0(nele,mmul),ak0(nele,mmul),                       &
     &bka(nele,mmul),aka(nele,mmul)
      common/mult1/benki,benkc(nele),r00(nele),irm(nele),nmu(nele)
      common/rand0/zfz(nzfz),iorg,mzu(nblz),bezr(3,nele),izu0,mmac,mcut
      common/rand1/exterr(nblz,40),extalign(nblz,3),tiltc(nblz),        &
     &tilts(nblz),mout2,icext(nblz),icextal(nblz)
      common/beo /aper(2),di0(2),dip0(2),ta(6,6)
      common/clo/dma,dmap,dkq,dqq,de0,ded,dsi,dech,dsm0,itco,itcro,itqv,&
     &iout
      common/qmodi/qw0(3),amp0,iq(3),iqmod,kpa(nele),iqmod6
      common/linop/bez(nele),elbe(nblo),bezb(nblo),ilin,nt,iprint,      &
     &ntco,eui,euii,nlin,bezl(nele)
      common/cororb/betam(nmon1,2),pam(nmon1,2),betac(ncor1,2),         &
     &pac(ncor1,2),bclorb(nmon1,2),nhmoni,nhcorr,nvmoni,nvcorr,         &
     &ncororb(nele)
      common/apert/apx(nele),apz(nele),ape(3,nele)
      common/clos/sigma0(2),iclo,ncorru,ncorrep
      common/combin/icomb0(20),icomb(ncom,20),ratio(ncom,20),           &
     &ratioe(nele),iratioe(nele),icoe
      common/seacom/ise,mesa,mp,m21,m22,m23,ise1,ise2,ise3,isea(nele)
      common/subres/qxt,qzt,tam1,tam2,isub,nta,nte,ipt,totl
      common/secom/rtc(9,18,10,5),rts(9,18,10,5),ire(12),ipr(5),irmod2
      common/secom1/dtr(10),nre,nur,nch,nqc,npp,nrr(5),nu(5)
      common/postr/dphix,dphiz,qx0,qz0,dres,dfft,cma1,cma2,             &
     &nstart,nstop,iskip,iconv,imad
      common/posti1/ipos,iav,iwg,ivox,ivoz,ires,ifh,toptit(5)
      common/posti2/kwtype,itf,icr,idis,icow,istw,iffw,nprint,ndafi
      common/skew/qwsk(2),betx(2),betz(2),alfx(2),alfz(2),iskew,nskew(6)
      common/pawc/hmal(nplo)
      common/tit/sixtit,commen,ithick
      common/co6d/clo6(3),clop6(3)
      common/dkic/dki(nele,3)
      common/beam/sigman(2,nbb),sigman2(2,nbb),sigmanq(2,nbb),          &
     &clobeam(6,nbb),beamoff(6,nbb),parbe(nele,5),track6d(6,npart),     &
     &ptnfac(nele),sigz,sige,partnum,parbe14,emitx,emity,emitz,gammar,  &
     &nbeam,ibbc,ibeco,ibtyp,lhc
      common/trom/ cotr(ntr,6),rrtr(ntr,6,6),imtr(nele)
      common/bb6d/ bbcu(nbb,12),ibb6d,imbb(nblz)
      common/wireco/ wirel(nele)
      common/acdipco/ acdipph(nele), nturn1(nele), nturn2(nele),        &
     &nturn3(nele), nturn4(nele)
      common/crabco/ crabph(nele),crabph2(nele),                        &
     &crabph3(nele),crabph4(nele)
      integer idz,itra
      double precision a2,al,as,at,chi0,chid,dp1,dps,exz,sigm
      common/syos/as(6,2,2,nele),at(6,2,2,nele),a2(6,2,2,nele),         &
     &al(6,2,2,nele),sigm(mpa),dps(mpa),idz(2)
      common/anf/chi0,chid,exz(2,6),dp1,itra
      integer ichrom,issss
      double precision alf0,amp,bet0,clo,clop,cro,xxtr,yytr
      common/tra/xxtr(mpa,2),yytr(mpa,2),amp(2),                        &
     &bet0(2),alf0(2),clo(2),clop(2)
      common/chrom/cro(2),issss(2),ichrom
      common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda, &
     &ej1,ejf1,rv
      integer icorr,idial,idptr,imod1,imod2,inorm,ipar,namp,ncor,nctype,&
     &ndimf,nmom,nmom1,nmom2,nord,nord1,nordf,nsix,nvar,nvar2,nvarf
      double precision dpmax,preda,weig1,weig2
      character*16 coel
      common/dial/preda,idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
      common/norf/nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
      common/tcorr/icorr,nctype,namp,nmom,nmom1,nmom2,weig1,weig2,dpmax,&
     &coel(10)
      double precision aml6,edcor
      common/sixdim/aml6(6,6),edcor(2)
      double precision aai,ampt,bbi,damp,smi,smizf,xsi,                 &
     &zsi
      integer napxto
      real tlim,time0,time1,time2,time3,trtime
! fixes for CPU time (for all versions, not just crlibm).
      real pretime,posttime,tottime
      common/xz/xsi(nblz),zsi(nblz),smi(nblz),smizf(nblz),              &
     &aai(nblz,mmul),bbi(nblz,mmul)
      common/damp/damp,ampt
      common/ttime/tlim,time0,time1,time2,time3,trtime,napxto,          &
     &pretime,posttime,tottime
      save
!-----------------------------------------------------------------------
*FOX  B D ;
*FOX  D V DA EXT SIGMDA NORD NVAR ; D V DA EXT DPDA NORD NVAR ;
*FOX  D V DA EXT DPDA1 NORD NVAR ; D V DA EXT RV NORD NVAR ;
*FOX  D V DA EXT XX NORD NVAR 2 ; D V DA EXT YY NORD NVAR 2 ;
*FOX  D V DA EXT EJ1 NORD NVAR ; D V DA EXT EJF1 NORD NVAR ;
*FOX  D V DA EXT ALDA NORD NVAR 2 6 ; D V DA EXT ASDA NORD NVAR 2 6 ;
*FOX  D V DA EXT ALDAQ NORD NVAR 2 6 ; D V DA EXT ASDAQ NORD NVAR 2 6 ;
*FOX  D V DA EXT SMIDA NORD NVAR MCOR ;
*FOX  D V DA INT XI NORD NVAR ; D V DA INT YI NORD NVAR ;
*FOX  D V RE INT EMBL ; D V RE INT TX ; D V RE INT TY ;
*FOX  D V RE INT DX ; D V RE INT DY ; D V RE INT CHI ;
*FOX  D V RE INT LEFF ; D V RE INT LIN ; D V RE INT RX ;
*FOX  D V RE INT RY ; D V RE INT CUR ;
*FOX  D V RE INT L ; D V RE INT ONE ; D V RE INT TWO ;
*FOX  D V RE INT C1M7 ;
*FOX  D V RE INT C1E3 ; D V RE INT C1M3 ;
*FOX  E D ;
*FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
 
!     magnetic rigidity
!hr05 chi = sqrt(e0*e0-pmap*pmap)*c1e6/clight
      chi = (sqrt(e0**2-pmap**2)*c1e6)/clight                            !hr05
 
!     The wire map consists of the following sequence of elements:
!     drift, tilt, drift, kick, drift, tilt, shift, and drift.
!     This way it is a zero length insertion. For details see upcoming paper.
 
      ix = ixcav
      tx = xrms(ix)
      ty = zrms(ix)
      dx = xpl(ix)
      dy = zpl(ix)
      embl = ek(ix)
      l = wirel(ix)
      cur = ed(ix)
 
!hr05 leff = embl/cos_rn(tx)/cos_rn(ty)
      leff = (embl/cos_rn(tx))/cos_rn(ty)                                !hr05
!hr05 rx = dx *cos_rn(tx)-embl*sin_rn(tx)/two
      rx = dx *cos_rn(tx)-(embl*sin_rn(tx))/two                          !hr05
!hr05 lin= dx *sin_rn(tx)+embl*cos_rn(tx)/two
      lin= dx *sin_rn(tx)+(embl*cos_rn(tx))/two                          !hr05
      ry = dy *cos_rn(ty)-lin *sin_rn(ty)
      lin= lin*cos_rn(ty)+dy  *sin_rn(ty)
 
 
*FOX  XX(1)=XX(1)*C1M3;
*FOX  XX(2)=XX(2)*C1M3;
*FOX  YY(1)=YY(1)*C1M3;
*FOX  YY(2)=YY(2)*C1M3;
 
!      CALL DRIFT(-EMBL/2)
*FOX  XX(1)=XX(1)-EMBL/TWO*YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-
*FOX  YY(1)*YY(1)-YY(2)*YY(2)) ;
*FOX  XX(2)=XX(2)-EMBL/TWO*YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-
*FOX  YY(1)*YY(1)-YY(2)*YY(2)) ;
!      CALL TILT(TX,TY)
*FOX  XX(2)=XX(2)-XX(1)*SIN(TX)*YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-
*FOX  YY(2)*YY(2))/COS(ATAN(YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-
*FOX  YY(1)*YY(1)-YY(2)*YY(2)))-TX) ;
*FOX  XX(1)=XX(1)*(COS(TX)-SIN(TX)*TAN(ATAN(YY(1)/SQRT((ONE+DPDA)*
*FOX  (ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TX)) ;
*FOX  YY(1)=SQRT((ONE+DPDA)*(ONE+DPDA)-YY(2)*YY(2))*SIN(ATAN(YY(1)/
*FOX  SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TX) ;
*FOX  XX(1)=XX(1)-XX(2)*SIN(TY)*YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-
*FOX  YY(1)*YY(1))/COS(ATAN(YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-
*FOX  YY(1)*YY(1)-YY(2)*YY(2)))-TY) ;
*FOX  XX(2)=XX(2)*(COS(TY)-SIN(TY)*TAN(ATAN(YY(2)/SQRT((ONE+DPDA)*
*FOX  (ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TY)) ;
*FOX  YY(2)=SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1))*SIN(ATAN(YY(2)/
*FOX  SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TY) ;
!      CALL DRIFT(LIN)
*FOX  XX(1)=XX(1)+LIN*YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-
*FOX  YY(2)*YY(2)) ;
*FOX  XX(2)=XX(2)+LIN*YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-
*FOX  YY(2)*YY(2)) ;
!      CALL KICK(L,CUR,LIN,RX,RY)
*FOX  XI=XX(1)-RX ;
*FOX  YI=XX(2)-RY ;
*FOX  YY(1)=YY(1)-C1M7*CUR/CHI*XI/(XI*XI+YI*YI)*(SQRT((LIN+L)*(LIN+L)+
*FOX  XI*XI+YI*YI)-SQRT((LIN-L)*(LIN-L)+XI*XI+YI*YI)) ;
*FOX  YY(2)=YY(2)-C1M7*CUR/CHI*YI/(XI*XI+YI*YI)*(SQRT((LIN+L)*(LIN+L)+
*FOX  XI*XI+YI*YI)-SQRT((LIN-L)*(LIN-L)+XI*XI+YI*YI)) ;
!     A = A-1E-7*2*L*CUR/CHI*XI/(XI*XI+YI*YI)
!     B = B-1E-7*2*L*CUR/CHI*YI/(XI*XI+YI*YI)
!
!     THIS HAPPENS WHEN THE EMBEDDING DRIFT LENGTH IS PRACTICALLY INFINITE.
!     CALL DRIFT(LEFF-LIN)
*FOX  XX(1)=XX(1)+(LEFF-LIN)*YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-
*FOX  YY(1)*YY(1)-YY(2)*YY(2)) ;
*FOX  XX(2)=XX(2)+(LEFF-LIN)*YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-
*FOX  YY(1)*YY(1)-YY(2)*YY(2)) ;
!      CALL INVTILT(TX,TY)
*FOX  XX(1)=XX(1)+XX(2)*SIN(TY)*YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-
*FOX  YY(1)*YY(1))/COS(ATAN(YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-
*FOX  YY(1)*YY(1)-YY(2)*YY(2)))+TY) ;
*FOX  XX(2)=XX(2)*(COS(TY)+SIN(TY)*TAN(ATAN(YY(2)/SQRT((ONE+DPDA)*
*FOX  (ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))+TY)) ;
*FOX  YY(2)=SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1))*SIN(ATAN(YY(2)/
*FOX  SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))+TY) ;
*FOX  XX(2)=XX(2)+XX(1)*SIN(TX)*YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-
*FOX  YY(2)*YY(2))/COS(ATAN(YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-
*FOX  YY(1)*YY(1)-YY(2)*YY(2)))+TX) ;
*FOX  XX(1)=XX(1)*(COS(TX)+SIN(TX)*TAN(ATAN(YY(1)/SQRT((ONE+DPDA)*
*FOX  (ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))+TX)) ;
*FOX  YY(1)=SQRT((ONE+DPDA)*(ONE+DPDA)-YY(2)*YY(2))*SIN(ATAN(YY(1)/
*FOX  SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))+TX) ;
!     CALL SHIFT(-EMBL*TAN(TX),-EMBL*TAN(TY)/COS(TX))
*FOX  XX(1)=XX(1)+EMBL*TAN(TX) ;
*FOX  XX(2)=XX(2)+EMBL*TAN(TY)/COS(TX) ;
!     CALL DRIFT(-EMBL/2)
*FOX  XX(1)=XX(1)-EMBL/TWO*YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-
*FOX  YY(1)*YY(1)-YY(2)*YY(2)) ;
*FOX  XX(2)=XX(2)-EMBL/TWO*YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-
*FOX  YY(1)*YY(1)-YY(2)*YY(2)) ;
 
*FOX  XX(1)=XX(1)*C1E3;
*FOX  XX(2)=XX(2)*C1E3;
*FOX  YY(1)=YY(1)*C1E3;
*FOX  YY(2)=YY(2)*C1E3;
 
!     DADAL AUTOMATIC INCLUSION
      end
      subroutine synoda
!-----------------------------------------------------------------------
!  SYNCHROTRON OSCILLATIONS
!        SPECIALLY PREPARED FOR NEW D.A.
!-----------------------------------------------------------------------
      implicit none
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
      integer ix,idaa
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,    &
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,    &
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele=1200,nblo=600,nper=16,nelb=140,nblz=20000,         &
     &nzfz = 300000,mmul = 20)
      parameter(nran = 2000000,ncom = 100,mran = 500,mpa = 6,nrco = 5,  &
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 99)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 350)
      double precision c180e0,c1e1,c1e12,c1e13,c1e15,c1e16,c1e2,c1e3,   &
     &c1e4,c1e6,c1m1,c1m7,c1m10,c1m12,c1m13,c1m15,c1m18,c1m2,c1m21,     &
     &c1m24,c1m3,c1m36,c1m38,c1m6,c1m9,c2e3,c4e3,crade,clight,four,half,&
     &one,pieni,pmae,pmap,three,two,zero
      parameter(pieni = 1d-38)
      parameter(zero = 0.0d0,half = 0.5d0,one = 1.0d0)
      parameter(two = 2.0d0,three = 3.0d0,four = 4.0d0)
      parameter(c1e1 = 1.0d1,c1e2 = 1.0d2,c1m2 = 1.0d-2)
      parameter(c1e3 = 1.0d3,c2e3 = 2.0d3,c4e3 = 4.0d3,c1e4 = 1.0d4)
      parameter(c1e12 = 1.0d12,c1e13 = 1.0d13,c1e15 = 1.0d15,c1e16 =    &
     &1.0d16)
      parameter(c180e0 = 180.0d0,c1e6 = 1.0d6)
      parameter(c1m1 = 1.0d-1,c1m3 = 1.0d-3,c1m6 = 1.0d-6,c1m7 = 1.0d-7)
      parameter(c1m9 = 1.0d-9,c1m10 = 1.0d-10,c1m12 = 1.0d-12)
      parameter(c1m13 = 1.0d-13,c1m15 = 1.0d-15)
      parameter(c1m18 = 1.0d-18,c1m21 = 1.0d-21,c1m24 = 1.0d-24)
      parameter(c1m36 = 1.0d-36,c1m38 = 1.0d-38)
!     electron mass from PDG, 2002
      parameter(pmap = 938.271998d0,pmae = .510998902d0)
      parameter(crade = 2.817940285d-15, clight = 2.99792458d8)
      integer iav,ibb6d,ibbc,ibeco,ibidu,ibtyp,ic,icext,icextal,iclo,   &
     &iclo6,iclo6r,icode,icoe,icomb,icomb0,iconv,icow,icr,idam,idfor,   &
     &idis,idp,ierro,iffw,ifh,iicav,il,ilin,imad,imbb,                  &
     &imc,imtr,iorg,iout,                                               &
     &ipos,ipr,iprint,ipt,iq,iqmod,iqmod6,iratioe,ird,ire,ires,         &
     &irew,irm,irmod2,ise,ise1,ise2,ise3,isea,iskew,iskip,istw,         &
     &isub,itco,itcro,itf,ithick,ition,itionc,itqv,its6d,iu,iver,ivox,  &
     &ivoz,iwg,ixcav,izu0,kanf,kp,kpa,kwtype,kz,lhc,m21,m22,m23,mblo,   &
     &mbloz,mcut,mel,mesa,mmac,mout2,mp,mper,mstr,msym,mtyp,mzu,napx,   &
     &napxo,nbeam,nch,ncororb,ncorrep,ncorru,ncy,ndafi,nde,nhcorr,      &
     &nhmoni,niu,nlin,nmu,npp,nprint,nqc,nre,nrr,nskew,                 &
     &nstart,nstop,nt,nta,ntco,nte,ntwin,nu,numl,numlr,nur,nvcorr,      &
     &nvmoni,nwr, nturn1, nturn2, nturn3, nturn4,numlcp,numlmax,nnuml
      double precision a,ak0,aka,alfx,alfz,amp0,aper,apx,apz,ape,bbcu,  &
     &bclorb,beamoff,benkc,benki,betac,betam,betx,betz,bk0,bka,bl1,bl2, &
     &clo6,clobeam,clop6,cma1,cma2,cotr,crad,de0,dech,ded,dfft,         &
     &di0,dip0,dki,dkq,dma,dmap,dphix,dphiz,dppoff,dpscor,dqq,dres,dsi, &
     &dsm0,dtr,e0,ed,ej,ejf,ek,el,elbe,emitx,emity,emitz,extalign,      &
     &exterr,eui,euii,gammar,hsy,hsyc,pac,pam,parbe,parbe14,partnum,    &
     &phas,phas0,phasc,pi,pi2,pisqrt,pma,ptnfac,qs,qw0,qwsk,qx0,qxt,qz0,&
     &qzt,r00,rad,rat,ratio,ratioe,rrtr,rtc,rts,rvf,                    &
     &sigcor,sige,sigma0,sigman,sigman2,sigmanq,sigmoff,sigz,sm,ta,tam1,&
     &tam2,tiltc,tilts,tlen,totl,track6d,xpl,xrms,zfz,zpl,zrms,wirel,   &
     &acdipph, crabph, bbbx, bbby, bbbs,                                &
     &crabph2, crabph3, crabph4
      real hmal
      character*16 bez,bezb,bezr,erbez,bezl
      character*80 toptit,sixtit,commen
      common/erro/ierro,erbez
      common/kons/pi,pi2,pisqrt,rad
      common/str /il,mper,mblo,mbloz,msym(nper),kanf,iu,ic(nblz)
      common/ell /ed(nele),el(nele),ek(nele),sm(nele),kz(nele),kp(nele)
      common/bbb /bbbx(nele), bbby(nele), bbbs(nele)
      common/pla /xpl(nele),xrms(nele),zpl(nele),zrms(nele)
      common/str2 /mel(nblo),mtyp(nblo,nelb),mstr(nblo)
      common/mat/a(nele,2,6),bl1(nblo,2,6),bl2(nblo,2,6)
      common/syos2/rvf(mpa)
      common/tra1/rat,idfor,napx,napxo,numl,niu(2),numlr,nde(2),nwr(4), &
     &ird,imc,irew,ntwin,iclo6,iclo6r,iver,ibidu,numlcp,numlmax,nnuml
      common/syn/qs,e0,pma,ej(mpa),ejf(mpa),phas0,phas,hsy(3),          &
     &crad,hsyc(nele),phasc(nele),dppoff,sigmoff(nblz),tlen,            &
     &iicav,itionc(nele),ition,idp,ncy,ixcav
      common/corcom/dpscor,sigcor,icode,idam,its6d
      common/multi/bk0(nele,mmul),ak0(nele,mmul),                       &
     &bka(nele,mmul),aka(nele,mmul)
      common/mult1/benki,benkc(nele),r00(nele),irm(nele),nmu(nele)
      common/rand0/zfz(nzfz),iorg,mzu(nblz),bezr(3,nele),izu0,mmac,mcut
      common/rand1/exterr(nblz,40),extalign(nblz,3),tiltc(nblz),        &
     &tilts(nblz),mout2,icext(nblz),icextal(nblz)
      common/beo /aper(2),di0(2),dip0(2),ta(6,6)
      common/clo/dma,dmap,dkq,dqq,de0,ded,dsi,dech,dsm0,itco,itcro,itqv,&
     &iout
      common/qmodi/qw0(3),amp0,iq(3),iqmod,kpa(nele),iqmod6
      common/linop/bez(nele),elbe(nblo),bezb(nblo),ilin,nt,iprint,      &
     &ntco,eui,euii,nlin,bezl(nele)
      common/cororb/betam(nmon1,2),pam(nmon1,2),betac(ncor1,2),         &
     &pac(ncor1,2),bclorb(nmon1,2),nhmoni,nhcorr,nvmoni,nvcorr,         &
     &ncororb(nele)
      common/apert/apx(nele),apz(nele),ape(3,nele)
      common/clos/sigma0(2),iclo,ncorru,ncorrep
      common/combin/icomb0(20),icomb(ncom,20),ratio(ncom,20),           &
     &ratioe(nele),iratioe(nele),icoe
      common/seacom/ise,mesa,mp,m21,m22,m23,ise1,ise2,ise3,isea(nele)
      common/subres/qxt,qzt,tam1,tam2,isub,nta,nte,ipt,totl
      common/secom/rtc(9,18,10,5),rts(9,18,10,5),ire(12),ipr(5),irmod2
      common/secom1/dtr(10),nre,nur,nch,nqc,npp,nrr(5),nu(5)
      common/postr/dphix,dphiz,qx0,qz0,dres,dfft,cma1,cma2,             &
     &nstart,nstop,iskip,iconv,imad
      common/posti1/ipos,iav,iwg,ivox,ivoz,ires,ifh,toptit(5)
      common/posti2/kwtype,itf,icr,idis,icow,istw,iffw,nprint,ndafi
      common/skew/qwsk(2),betx(2),betz(2),alfx(2),alfz(2),iskew,nskew(6)
      common/pawc/hmal(nplo)
      common/tit/sixtit,commen,ithick
      common/co6d/clo6(3),clop6(3)
      common/dkic/dki(nele,3)
      common/beam/sigman(2,nbb),sigman2(2,nbb),sigmanq(2,nbb),          &
     &clobeam(6,nbb),beamoff(6,nbb),parbe(nele,5),track6d(6,npart),     &
     &ptnfac(nele),sigz,sige,partnum,parbe14,emitx,emity,emitz,gammar,  &
     &nbeam,ibbc,ibeco,ibtyp,lhc
      common/trom/ cotr(ntr,6),rrtr(ntr,6,6),imtr(nele)
      common/bb6d/ bbcu(nbb,12),ibb6d,imbb(nblz)
      common/wireco/ wirel(nele)
      common/acdipco/ acdipph(nele), nturn1(nele), nturn2(nele),        &
     &nturn3(nele), nturn4(nele)
      common/crabco/ crabph(nele),crabph2(nele),                        &
     &crabph3(nele),crabph4(nele)
      integer idz,itra
      double precision a2,al,as,at,chi0,chid,dp1,dps,exz,sigm
      common/syos/as(6,2,2,nele),at(6,2,2,nele),a2(6,2,2,nele),         &
     &al(6,2,2,nele),sigm(mpa),dps(mpa),idz(2)
      common/anf/chi0,chid,exz(2,6),dp1,itra
      integer ichrom,is
      double precision alf0,amp,bet0,clo,clop,cro,x,y
      common/tra/x(mpa,2),y(mpa,2),amp(2),bet0(2),alf0(2),clo(2),clop(2)
      common/chrom/cro(2),is(2),ichrom
      integer icorr,idial,idptr,imod1,imod2,inorm,ipar,namp,ncor,nctype,&
     &ndimf,nmom,nmom1,nmom2,nord,nord1,nordf,nsix,nvar,nvar2,nvarf
      double precision dpmax,preda,weig1,weig2
      character*16 coel
      common/dial/preda,idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
      common/norf/nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
      common/tcorr/icorr,nctype,namp,nmom,nmom1,nmom2,weig1,weig2,dpmax,&
     &coel(10)
      common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda, &
     &ej1,ejf1,rv
      integer numx
      double precision e0f
      common/main4/ e0f,numx
      save
!-----------------------------------------------------------------------
*FOX  B D ;
*FOX  D V DA EXT SIGMDA NORD NVAR ; D V DA EXT DPDA NORD NVAR ;
*FOX  D V DA EXT DPDA1 NORD NVAR ; D V DA EXT RV NORD NVAR ;
*FOX  D V DA EXT XX NORD NVAR 2 ; D V DA EXT YY NORD NVAR 2 ;
*FOX  D V DA EXT EJ1 NORD NVAR ; D V DA EXT EJF1 NORD NVAR ;
*FOX  D V DA EXT ALDA NORD NVAR 2 6 ; D V DA EXT ASDA NORD NVAR 2 6 ;
*FOX  D V DA EXT ALDAQ NORD NVAR 2 6 ; D V DA EXT ASDAQ NORD NVAR 2 6 ;
*FOX  D V DA EXT SMIDA NORD NVAR MCOR ;
*FOX  D V RE INT E0 ; D V RE INT PMA ; D V RE EXT E0F ;
*FOX  D V RE INT HSY 3 ; D V RE INT PHAS ;
*FOX  D V RE EXT ED NELE ; D V RE EXT HSYC NELE ;
*FOX  D V RE EXT PHASC NELE ;
*FOX  D V RE INT C1E3 ; D V RE INT ONE ;
*FOX  D V IN EXT ITIONC NELE ; D V IN INT ITION ; D V IN INT IX ;
*FOX  E D ;
*FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
      ix=ixcav
      if(kz(ix).eq.12) then
*FOX  EJ1=EJ1+ED(IX)*SIN(HSYC(IX)*SIGMDA/C1E3*
*FOX  ITIONC(IX)+PHASC(IX)) ;
      else
*FOX  EJ1=EJ1+HSY(1)*SIN(HSY(3)*SIGMDA/C1E3*ITION+PHAS) ;
      endif
*FOX  EJF1=SQRT(EJ1*EJ1-PMA*PMA) ;
*FOX  DPDA1=(EJF1-E0F)/E0F*C1E3 ;
      return
      end
      subroutine errff(xx,yy,wx,wy)
!----------------------------------------------------------------------*
! PURPOSE:                                                             *
!   MODIFICATION OF WWERF, DOUBLE PRECISION COMPLEX ERROR FUNCTION,    *
!   WRITTEN AT CERN BY K. KOELBIG.                                     *
!   TAKEN FROM MAD8                                                    *
!   VERSION FOR MAP PRODUCTION USING BERZ'S DA PACKAGE                 *
! INPUT:                                                               *
!   XX, YY    (REAL)    ARGUMENT TO CERF.                              *
! OUTPUT:                                                              *
!   WX, WY    (REAL)    FUNCTION RESULT.                               *
!----------------------------------------------------------------------*
      implicit none
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
      integer n,n1,nc,nuu,nuu1,idaa
      double precision cc,dare,dum,xlim,ylim
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,    &
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,    &
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele=1200,nblo=600,nper=16,nelb=140,nblz=20000,         &
     &nzfz = 300000,mmul = 20)
      parameter(nran = 2000000,ncom = 100,mran = 500,mpa = 6,nrco = 5,  &
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 99)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 350)
      double precision c180e0,c1e1,c1e12,c1e13,c1e15,c1e16,c1e2,c1e3,   &
     &c1e4,c1e6,c1m1,c1m7,c1m10,c1m12,c1m13,c1m15,c1m18,c1m2,c1m21,     &
     &c1m24,c1m3,c1m36,c1m38,c1m6,c1m9,c2e3,c4e3,crade,clight,four,half,&
     &one,pieni,pmae,pmap,three,two,zero
      parameter(pieni = 1d-38)
      parameter(zero = 0.0d0,half = 0.5d0,one = 1.0d0)
      parameter(two = 2.0d0,three = 3.0d0,four = 4.0d0)
      parameter(c1e1 = 1.0d1,c1e2 = 1.0d2,c1m2 = 1.0d-2)
      parameter(c1e3 = 1.0d3,c2e3 = 2.0d3,c4e3 = 4.0d3,c1e4 = 1.0d4)
      parameter(c1e12 = 1.0d12,c1e13 = 1.0d13,c1e15 = 1.0d15,c1e16 =    &
     &1.0d16)
      parameter(c180e0 = 180.0d0,c1e6 = 1.0d6)
      parameter(c1m1 = 1.0d-1,c1m3 = 1.0d-3,c1m6 = 1.0d-6,c1m7 = 1.0d-7)
      parameter(c1m9 = 1.0d-9,c1m10 = 1.0d-10,c1m12 = 1.0d-12)
      parameter(c1m13 = 1.0d-13,c1m15 = 1.0d-15)
      parameter(c1m18 = 1.0d-18,c1m21 = 1.0d-21,c1m24 = 1.0d-24)
      parameter(c1m36 = 1.0d-36,c1m38 = 1.0d-38)
!     electron mass from PDG, 2002
      parameter(pmap = 938.271998d0,pmae = .510998902d0)
      parameter(crade = 2.817940285d-15, clight = 2.99792458d8)
      parameter(cc = 1.12837916709551d0)
      parameter(xlim = 5.33d0)
      parameter(ylim = 4.29d0)
      integer iav,ibb6d,ibbc,ibeco,ibidu,ibtyp,ic,icext,icextal,iclo,   &
     &iclo6,iclo6r,icode,icoe,icomb,icomb0,iconv,icow,icr,idam,idfor,   &
     &idis,idp,ierro,iffw,ifh,iicav,il,ilin,imad,imbb,                  &
     &imc,imtr,iorg,iout,                                               &
     &ipos,ipr,iprint,ipt,iq,iqmod,iqmod6,iratioe,ird,ire,ires,         &
     &irew,irm,irmod2,ise,ise1,ise2,ise3,isea,iskew,iskip,istw,         &
     &isub,itco,itcro,itf,ithick,ition,itionc,itqv,its6d,iu,iver,ivox,  &
     &ivoz,iwg,ixcav,izu0,kanf,kp,kpa,kwtype,kz,lhc,m21,m22,m23,mblo,   &
     &mbloz,mcut,mel,mesa,mmac,mout2,mp,mper,mstr,msym,mtyp,mzu,napx,   &
     &napxo,nbeam,nch,ncororb,ncorrep,ncorru,ncy,ndafi,nde,nhcorr,      &
     &nhmoni,niu,nlin,nmu,npp,nprint,nqc,nre,nrr,nskew,                 &
     &nstart,nstop,nt,nta,ntco,nte,ntwin,nu,numl,numlr,nur,nvcorr,      &
     &nvmoni,nwr, nturn1, nturn2, nturn3, nturn4,numlcp,numlmax,nnuml
      double precision a,ak0,aka,alfx,alfz,amp0,aper,apx,apz,ape,bbcu,  &
     &bclorb,beamoff,benkc,benki,betac,betam,betx,betz,bk0,bka,bl1,bl2, &
     &clo6,clobeam,clop6,cma1,cma2,cotr,crad,de0,dech,ded,dfft,         &
     &di0,dip0,dki,dkq,dma,dmap,dphix,dphiz,dppoff,dpscor,dqq,dres,dsi, &
     &dsm0,dtr,e0,ed,ej,ejf,ek,el,elbe,emitx,emity,emitz,extalign,      &
     &exterr,eui,euii,gammar,hsy,hsyc,pac,pam,parbe,parbe14,partnum,    &
     &phas,phas0,phasc,pi,pi2,pisqrt,pma,ptnfac,qs,qw0,qwsk,qx0,qxt,qz0,&
     &qzt,r00,rad,rat,ratio,ratioe,rrtr,rtc,rts,rvf,                    &
     &sigcor,sige,sigma0,sigman,sigman2,sigmanq,sigmoff,sigz,sm,ta,tam1,&
     &tam2,tiltc,tilts,tlen,totl,track6d,xpl,xrms,zfz,zpl,zrms,wirel,   &
     &acdipph, crabph, bbbx, bbby, bbbs,                                &
     &crabph2, crabph3, crabph4
      real hmal
      character*16 bez,bezb,bezr,erbez,bezl
      character*80 toptit,sixtit,commen
      common/erro/ierro,erbez
      common/kons/pi,pi2,pisqrt,rad
      common/str /il,mper,mblo,mbloz,msym(nper),kanf,iu,ic(nblz)
      common/ell /ed(nele),el(nele),ek(nele),sm(nele),kz(nele),kp(nele)
      common/bbb /bbbx(nele), bbby(nele), bbbs(nele)
      common/pla /xpl(nele),xrms(nele),zpl(nele),zrms(nele)
      common/str2 /mel(nblo),mtyp(nblo,nelb),mstr(nblo)
      common/mat/a(nele,2,6),bl1(nblo,2,6),bl2(nblo,2,6)
      common/syos2/rvf(mpa)
      common/tra1/rat,idfor,napx,napxo,numl,niu(2),numlr,nde(2),nwr(4), &
     &ird,imc,irew,ntwin,iclo6,iclo6r,iver,ibidu,numlcp,numlmax,nnuml
      common/syn/qs,e0,pma,ej(mpa),ejf(mpa),phas0,phas,hsy(3),          &
     &crad,hsyc(nele),phasc(nele),dppoff,sigmoff(nblz),tlen,            &
     &iicav,itionc(nele),ition,idp,ncy,ixcav
      common/corcom/dpscor,sigcor,icode,idam,its6d
      common/multi/bk0(nele,mmul),ak0(nele,mmul),                       &
     &bka(nele,mmul),aka(nele,mmul)
      common/mult1/benki,benkc(nele),r00(nele),irm(nele),nmu(nele)
      common/rand0/zfz(nzfz),iorg,mzu(nblz),bezr(3,nele),izu0,mmac,mcut
      common/rand1/exterr(nblz,40),extalign(nblz,3),tiltc(nblz),        &
     &tilts(nblz),mout2,icext(nblz),icextal(nblz)
      common/beo /aper(2),di0(2),dip0(2),ta(6,6)
      common/clo/dma,dmap,dkq,dqq,de0,ded,dsi,dech,dsm0,itco,itcro,itqv,&
     &iout
      common/qmodi/qw0(3),amp0,iq(3),iqmod,kpa(nele),iqmod6
      common/linop/bez(nele),elbe(nblo),bezb(nblo),ilin,nt,iprint,      &
     &ntco,eui,euii,nlin,bezl(nele)
      common/cororb/betam(nmon1,2),pam(nmon1,2),betac(ncor1,2),         &
     &pac(ncor1,2),bclorb(nmon1,2),nhmoni,nhcorr,nvmoni,nvcorr,         &
     &ncororb(nele)
      common/apert/apx(nele),apz(nele),ape(3,nele)
      common/clos/sigma0(2),iclo,ncorru,ncorrep
      common/combin/icomb0(20),icomb(ncom,20),ratio(ncom,20),           &
     &ratioe(nele),iratioe(nele),icoe
      common/seacom/ise,mesa,mp,m21,m22,m23,ise1,ise2,ise3,isea(nele)
      common/subres/qxt,qzt,tam1,tam2,isub,nta,nte,ipt,totl
      common/secom/rtc(9,18,10,5),rts(9,18,10,5),ire(12),ipr(5),irmod2
      common/secom1/dtr(10),nre,nur,nch,nqc,npp,nrr(5),nu(5)
      common/postr/dphix,dphiz,qx0,qz0,dres,dfft,cma1,cma2,             &
     &nstart,nstop,iskip,iconv,imad
      common/posti1/ipos,iav,iwg,ivox,ivoz,ires,ifh,toptit(5)
      common/posti2/kwtype,itf,icr,idis,icow,istw,iffw,nprint,ndafi
      common/skew/qwsk(2),betx(2),betz(2),alfx(2),alfz(2),iskew,nskew(6)
      common/pawc/hmal(nplo)
      common/tit/sixtit,commen,ithick
      common/co6d/clo6(3),clop6(3)
      common/dkic/dki(nele,3)
      common/beam/sigman(2,nbb),sigman2(2,nbb),sigmanq(2,nbb),          &
     &clobeam(6,nbb),beamoff(6,nbb),parbe(nele,5),track6d(6,npart),     &
     &ptnfac(nele),sigz,sige,partnum,parbe14,emitx,emity,emitz,gammar,  &
     &nbeam,ibbc,ibeco,ibtyp,lhc
      common/trom/ cotr(ntr,6),rrtr(ntr,6,6),imtr(nele)
      common/bb6d/ bbcu(nbb,12),ibb6d,imbb(nblz)
      common/wireco/ wirel(nele)
      common/acdipco/ acdipph(nele), nturn1(nele), nturn2(nele),        &
     &nturn3(nele), nturn4(nele)
      common/crabco/ crabph(nele),crabph2(nele),                        &
     &crabph3(nele),crabph4(nele)
      integer idz,itra
      double precision a2,al,as,at,chi0,chid,dp1,dps,exz,sigm
      common/syos/as(6,2,2,nele),at(6,2,2,nele),a2(6,2,2,nele),         &
     &al(6,2,2,nele),sigm(mpa),dps(mpa),idz(2)
      common/anf/chi0,chid,exz(2,6),dp1,itra
      integer ichrom,issss
      double precision alf0,amp,bet0,clo,clop,cro,xxtr,yytr
      common/tra/xxtr(mpa,2),yytr(mpa,2),amp(2),                        &
     &bet0(2),alf0(2),clo(2),clop(2)
      common/chrom/cro(2),issss(2),ichrom
      integer icorr,idial,idptr,imod1,imod2,inorm,ipar,namp,ncor,nctype,&
     &ndimf,nmom,nmom1,nmom2,nord,nord1,nordf,nsix,nvar,nvar2,nvarf
      double precision dpmax,preda,weig1,weig2
      character*16 coel
      common/dial/preda,idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
      common/norf/nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
      common/tcorr/icorr,nctype,namp,nmom,nmom1,nmom2,weig1,weig2,dpmax,&
     &coel(10)
      double precision aml6,edcor
      common/sixdim/aml6(6,6),edcor(2)
      double precision aai,ampt,bbi,damp,smi,smizf,xsi,                 &
     &zsi
      integer napxto
      real tlim,time0,time1,time2,time3,trtime
! fixes for CPU time (for all versions, not just crlibm).
      real pretime,posttime,tottime
      common/xz/xsi(nblz),zsi(nblz),smi(nblz),smizf(nblz),              &
     &aai(nblz,mmul),bbi(nblz,mmul)
      common/damp/damp,ampt
      common/ttime/tlim,time0,time1,time2,time3,trtime,napxto,          &
     &pretime,posttime,tottime
      save
!-----------------------------------------------------------------------
*FOX  B D ;
*FOX  D V DA EXT XX NORD NVAR ; D V DA EXT YY NORD NVAR ;
*FOX  D V DA EXT WX NORD NVAR ; D V DA EXT WY NORD NVAR ;
*FOX  D V DA INT X NORD NVAR ; D V DA INT Y NORD NVAR ;
*FOX  D V DA INT Q NORD NVAR ; D V DA INT H NORD NVAR ;
*FOX  D V DA INT XH NORD NVAR ; D V DA INT YH NORD NVAR ;
*FOX  D V DA INT RX NORD NVAR 33 ; D V DA INT RY NORD NVAR 33 ;
*FOX  D V DA INT TX NORD NVAR ; D V DA INT TN NORD NVAR ;
*FOX  D V DA INT TY NORD NVAR ;D V DA INT SAUX NORD NVAR ;
*FOX  D V DA INT SX NORD NVAR ; D V DA INT SY NORD NVAR ;
*FOX  D V DA INT XL NORD NVAR ;
*FOX  D V RE INT XLIM ; D V RE INT YLIM ; D V RE INT TWO ;
*FOX  D V RE INT ONE ; D V RE INT ZERO ; D V RE INT HALF ;
*FOX  D V RE INT CC ; D V RE INT DUM ;
*FOX  D V IN INT NC ; D V IN INT N ; D V IN INT N1 ; D V IN INT NUU ;
*FOX  D V IN INT NUU1 ; D V IN INT NCC ;
*FOX  D F RE DARE 1 ;
*FOX  E D ;
*FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
*FOX  X=XX ;
*FOX  Y=YY ;
      if(dare(x).lt.zero) then
        write(*,*) ' Problem in DA complex error function: dare(x) < 0'
*FOX    X=-X ;
      endif
      if(dare(y).lt.zero) then
        write(*,*) ' Problem in DA complex error function: dare(y) < 0'
*FOX    Y=-Y ;
      endif
      if(dare(y).lt.ylim.and.dare(x).lt.xlim) then
*FOX    Q=(ONE-Y/YLIM)*SQRT(ONE-X*X/XLIM/XLIM) ;
*FOX    DUM=3.2D0 ;
*FOX    H=ONE/(DUM*Q) ;
        nc=7+int(23.0d0*dare(q))
*FOX    XL=EXP((1-NC)*LOG(H)) ;
*FOX    XH=Y+HALF/H ;
*FOX    YH=X ;
        nuu=10+int(21.0d0*dare(q))
        nuu1=nuu+1
*FOX    RX(NUU1)=ZERO ;
*FOX    RY(NUU1)=ZERO ;
        do 10 n=nuu,1,-1
          n1=n+1
*FOX      TX=XH+N*RX(N1) ;
*FOX      TY=YH-N*RY(N1) ;
*FOX      TN=TX*TX+TY*TY ;
*FOX      RX(N)=HALF*TX/TN ;
*FOX      RY(N)=HALF*TY/TN ;
   10   continue
*FOX    SX=ZERO ;
*FOX    SY=ZERO ;
        do 20 n=nc,1,-1
*FOX      SAUX=SX+XL ;
*FOX      SX=RX(N)*SAUX-RY(N)*SY ;
*FOX      SY=RX(N)*SY+RY(N)*SAUX ;
*FOX      XL=H*XL ;
   20   continue
*FOX    WX=CC*SX ;
*FOX    WY=CC*SY ;
      else
*FOX    XH=Y ;
*FOX    YH=X ;
*FOX    RX(1)=ZERO ;
*FOX    RY(1)=ZERO ;
        do 30 n=9,1,-1
*FOX      TX=XH+N*RX(1) ;
*FOX      TY=YH-N*RY(1) ;
*FOX      TN=TX*TX+TY*TY ;
*FOX      RX(1)=HALF*TX/TN ;
*FOX      RY(1)=HALF*TY/TN ;
   30   continue
*FOX    WX=CC*RX(1) ;
*FOX    WY=CC*RY(1) ;
      endif
!      if(dare(y).eq.0.) then
!*FOX    WX=EXP(-X*X) ;
!      endif
!hr05 if(dare(yy).lt.0.) then
      if(dare(yy).lt.0.d0) then                                          !hr05
*FOX    WX=TWO*EXP(Y*Y-X*X)*COS(TWO*X*Y)-WX ;
*FOX    WY=-TWO*EXP(Y*Y-X*X)*SIN(TWO*X*Y)-WY ;
!hr05   if(dare(xx).gt.0.) then
        if(dare(xx).gt.0.d0) then                                        !hr05
*FOX      WY=-WY ;
        endif
      else
!hr05   if(dare(xx).lt.0.) then
        if(dare(xx).lt.0.d0) then                                        !hr05
*FOX      WY=-WY ;
        endif
      endif
!     DADAL AUTOMATIC INCLUSION
      return
      end
      subroutine beaminf(track,param,sigzs,bcu,ibb,ne,ibbc)
!-----------------------------------------------------------------------
!
!   Hirata's 6d beam-beam from BBC
!   DA VERSION for SIXTRACK courtesy Peter Leunissen
!   January 1999
!
!-----------------------------------------------------------------------
      implicit none
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
      integer ibb,ibbc,ne,nsli,idaa
      double precision alpha,bcu,calpha,cphi,f,param,phi,salpha,sigzs,  &
     &sphi,star,tphi,phi2,cphi2,sphi2,tphi2
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,    &
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,    &
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele=1200,nblo=600,nper=16,nelb=140,nblz=20000,         &
     &nzfz = 300000,mmul = 20)
      parameter(nran = 2000000,ncom = 100,mran = 500,mpa = 6,nrco = 5,  &
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 99)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 350)
      double precision c180e0,c1e1,c1e12,c1e13,c1e15,c1e16,c1e2,c1e3,   &
     &c1e4,c1e6,c1m1,c1m7,c1m10,c1m12,c1m13,c1m15,c1m18,c1m2,c1m21,     &
     &c1m24,c1m3,c1m36,c1m38,c1m6,c1m9,c2e3,c4e3,crade,clight,four,half,&
     &one,pieni,pmae,pmap,three,two,zero
      parameter(pieni = 1d-38)
      parameter(zero = 0.0d0,half = 0.5d0,one = 1.0d0)
      parameter(two = 2.0d0,three = 3.0d0,four = 4.0d0)
      parameter(c1e1 = 1.0d1,c1e2 = 1.0d2,c1m2 = 1.0d-2)
      parameter(c1e3 = 1.0d3,c2e3 = 2.0d3,c4e3 = 4.0d3,c1e4 = 1.0d4)
      parameter(c1e12 = 1.0d12,c1e13 = 1.0d13,c1e15 = 1.0d15,c1e16 =    &
     &1.0d16)
      parameter(c180e0 = 180.0d0,c1e6 = 1.0d6)
      parameter(c1m1 = 1.0d-1,c1m3 = 1.0d-3,c1m6 = 1.0d-6,c1m7 = 1.0d-7)
      parameter(c1m9 = 1.0d-9,c1m10 = 1.0d-10,c1m12 = 1.0d-12)
      parameter(c1m13 = 1.0d-13,c1m15 = 1.0d-15)
      parameter(c1m18 = 1.0d-18,c1m21 = 1.0d-21,c1m24 = 1.0d-24)
      parameter(c1m36 = 1.0d-36,c1m38 = 1.0d-38)
!     electron mass from PDG, 2002
      parameter(pmap = 938.271998d0,pmae = .510998902d0)
      parameter(crade = 2.817940285d-15, clight = 2.99792458d8)
      integer icorr,idial,idptr,imod1,imod2,inorm,ipar,namp,ncor,nctype,&
     &ndimf,nmom,nmom1,nmom2,nord,nord1,nordf,nsix,nvar,nvar2,nvarf
      double precision dpmax,preda,weig1,weig2
      character*16 coel
      common/dial/preda,idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
      common/norf/nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
      common/tcorr/icorr,nctype,namp,nmom,nmom1,nmom2,weig1,weig2,dpmax,&
     &coel(10)
      ! JBG Increaseing param to dimension 5 for xstr
      dimension param(nele,5),bcu(nbb,12),star(3,mbea)
      save
!-----------------------------------------------------------------------
*FOX  B D ;
*FOX  D V DA EXT TRACK NORD NVAR 6 ;
*FOX  E D ;
*FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
      phi=param(ne,1)
      nsli=param(ne,2)
      alpha=param(ne,3)
      phi2=param(ne,5)
!hr05 f=param(ne,4)/nsli
      f=param(ne,4)/dble(nsli)                                           !hr05
      sphi=sin_rn(phi)
      sphi2=sin_rn(phi2)
      cphi=cos_rn(phi)
      cphi2=cos_rn(phi2)
      tphi=tan_rn(phi)
      tphi2=tan_rn(phi2)
      salpha=sin_rn(alpha)
      calpha=cos_rn(alpha)
!     define slices
      call stsld(star,cphi2,sphi2,sigzs,nsli,calpha,salpha)
      call boostf(sphi,cphi,tphi,salpha,calpha,track)
      call sbcf(star,cphi,cphi2,nsli,f,ibb,bcu,track,ibbc)
      call boostif(sphi,cphi,tphi,salpha,calpha,track)
!     DADAL AUTOMATIC INCLUSION
      return
      end
      subroutine boostf(sphi,cphi,tphi,salpha,calpha,track)
!-----------------------------------------------------------------------
!
!   Hirata's 6d beam-beam from BBC
!   DA VERSION for SIXTRACK courtesy Peter Leunissen
!   January 1999
!
! BOOSTF Boost Operation *******************************************
!    P,Q,E are all normalized by P0
!-----------------------------------------------------------------------
      implicit none
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
      integer idaa
      double precision calpha,cphi,salpha,sphi,tphi,cphi2,sphi2,tphi2    &
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,    &
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,    &
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele=1200,nblo=600,nper=16,nelb=140,nblz=20000,         &
     &nzfz = 300000,mmul = 20)
      parameter(nran = 2000000,ncom = 100,mran = 500,mpa = 6,nrco = 5,  &
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 99)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 350)
      double precision c180e0,c1e1,c1e12,c1e13,c1e15,c1e16,c1e2,c1e3,   &
     &c1e4,c1e6,c1m1,c1m7,c1m10,c1m12,c1m13,c1m15,c1m18,c1m2,c1m21,     &
     &c1m24,c1m3,c1m36,c1m38,c1m6,c1m9,c2e3,c4e3,crade,clight,four,half,&
     &one,pieni,pmae,pmap,three,two,zero
      parameter(pieni = 1d-38)
      parameter(zero = 0.0d0,half = 0.5d0,one = 1.0d0)
      parameter(two = 2.0d0,three = 3.0d0,four = 4.0d0)
      parameter(c1e1 = 1.0d1,c1e2 = 1.0d2,c1m2 = 1.0d-2)
      parameter(c1e3 = 1.0d3,c2e3 = 2.0d3,c4e3 = 4.0d3,c1e4 = 1.0d4)
      parameter(c1e12 = 1.0d12,c1e13 = 1.0d13,c1e15 = 1.0d15,c1e16 =    &
     &1.0d16)
      parameter(c180e0 = 180.0d0,c1e6 = 1.0d6)
      parameter(c1m1 = 1.0d-1,c1m3 = 1.0d-3,c1m6 = 1.0d-6,c1m7 = 1.0d-7)
      parameter(c1m9 = 1.0d-9,c1m10 = 1.0d-10,c1m12 = 1.0d-12)
      parameter(c1m13 = 1.0d-13,c1m15 = 1.0d-15)
      parameter(c1m18 = 1.0d-18,c1m21 = 1.0d-21,c1m24 = 1.0d-24)
      parameter(c1m36 = 1.0d-36,c1m38 = 1.0d-38)
!     electron mass from PDG, 2002
      parameter(pmap = 938.271998d0,pmae = .510998902d0)
      parameter(crade = 2.817940285d-15, clight = 2.99792458d8)
      integer icorr,idial,idptr,imod1,imod2,inorm,ipar,namp,ncor,nctype,&
     &ndimf,nmom,nmom1,nmom2,nord,nord1,nordf,nsix,nvar,nvar2,nvarf
      double precision dpmax,preda,weig1,weig2
      character*16 coel
      common/dial/preda,idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
      common/norf/nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
      common/tcorr/icorr,nctype,namp,nmom,nmom1,nmom2,weig1,weig2,dpmax,&
     &coel(10)
      save
!-----------------------------------------------------------------------
*FOX  B D ;
*FOX  D V DA EXT TRACK NORD NVAR 6 ; D V DA INT A NORD NVAR ;
*FOX  D V DA INT H NORD NVAR ; D V DA INT SQR1A NORD NVAR ;
*FOX  D V DA INT A1 NORD NVAR ; D V DA INT HD1 NORD NVAR ;
*FOX  D V DA INT H1X NORD NVAR ; D V DA INT H1Y NORD NVAR ;
*FOX  D V DA INT H1Z NORD NVAR ; D V DA INT X1 NORD NVAR ;
*FOX  D V DA INT Y1 NORD NVAR ;
*FOX  D V RE EXT SPHI ; D V RE EXT CPHI ; D V RE EXT TPHI ;
*FOX  D V RE EXT SPHI2 ; D V RE EXT CPHI2 ; D V RE EXT TPHI2 ;
*FOX  D V RE EXT SALPHA ; D V RE EXT CALPHA ;
*FOX  D V RE INT ONE ; D V RE INT C1E3 ;
*FOX  D V DA INT DET NORD NVAR ; D V DA INT H1 NORD NVAR ;
*FOX  E D ;
*FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
*FOX    H=TRACK(6)+ONE-SQRT((ONE+TRACK(6))*(ONE+TRACK(6))-
*FOX    TRACK(2)*TRACK(2)-TRACK(4)*TRACK(4)) ;
*FOX    TRACK(6)=TRACK(6)-CALPHA*TPHI*TRACK(2)
*FOX              -TRACK(4)*SALPHA*TPHI+H*TPHI*TPHI ;
*FOX    TRACK(2)=(TRACK(2)-TPHI*H*CALPHA)/CPHI ;
*FOX    TRACK(4)=(TRACK(4)-TPHI*H*SALPHA)/CPHI ;
*FOX    HD1=SQRT((ONE+TRACK(6))*(ONE+TRACK(6))-TRACK(2)*TRACK(2)-
*FOX    TRACK(4)*TRACK(4)) ;
*FOX    H1X=TRACK(2)/HD1 ;
*FOX    H1Y=TRACK(4)/HD1 ;
*FOX    H1Z=ONE-(ONE+TRACK(6))/HD1 ;
*FOX    X1=CALPHA*TPHI*TRACK(5)+(ONE+CALPHA*SPHI*H1X)*TRACK(1)
*FOX       +TRACK(3)*SALPHA*SPHI*H1X ;
*FOX    Y1=SALPHA*TPHI*TRACK(5)+(ONE+SALPHA*SPHI*H1Y)*TRACK(3)
*FOX       +TRACK(1)*CALPHA*SPHI*H1Y ;
*FOX    TRACK(5)=TRACK(5)/CPHI+H1Z*(SPHI*CALPHA*TRACK(1)
*FOX       +SPHI*SALPHA*TRACK(3)) ;
*FOX    TRACK(1)=X1 ;
*FOX    TRACK(3)=Y1 ;
!     DADAL AUTOMATIC INCLUSION
      return
      end
      subroutine sbcf(star,cphi,cphi2,nsli,f,ibb,bcu,track,ibbc)
!-----------------------------------------------------------------------
!
!   Hirata's 6d beam-beam from BBC
!   DA VERSION for SIXTRACK courtesy Peter Leunissen
!   January 1999
!
!-----------------------------------------------------------------------
!**SBCF ***Synchro-Beam for headon collision*********************
!****************************************************************
      implicit none
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
      integer ibb,ibbc,ibbc1,jsli,nsli,idaa
      double precision bcu,cphi,cphi2,dare,f,sfac,star
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,    &
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,    &
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele=1200,nblo=600,nper=16,nelb=140,nblz=20000,         &
     &nzfz = 300000,mmul = 20)
      parameter(nran = 2000000,ncom = 100,mran = 500,mpa = 6,nrco = 5,  &
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 99)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 350)
      double precision c180e0,c1e1,c1e12,c1e13,c1e15,c1e16,c1e2,c1e3,   &
     &c1e4,c1e6,c1m1,c1m7,c1m10,c1m12,c1m13,c1m15,c1m18,c1m2,c1m21,     &
     &c1m24,c1m3,c1m36,c1m38,c1m6,c1m9,c2e3,c4e3,crade,clight,four,half,&
     &one,pieni,pmae,pmap,three,two,zero
      parameter(pieni = 1d-38)
      parameter(zero = 0.0d0,half = 0.5d0,one = 1.0d0)
      parameter(two = 2.0d0,three = 3.0d0,four = 4.0d0)
      parameter(c1e1 = 1.0d1,c1e2 = 1.0d2,c1m2 = 1.0d-2)
      parameter(c1e3 = 1.0d3,c2e3 = 2.0d3,c4e3 = 4.0d3,c1e4 = 1.0d4)
      parameter(c1e12 = 1.0d12,c1e13 = 1.0d13,c1e15 = 1.0d15,c1e16 =    &
     &1.0d16)
      parameter(c180e0 = 180.0d0,c1e6 = 1.0d6)
      parameter(c1m1 = 1.0d-1,c1m3 = 1.0d-3,c1m6 = 1.0d-6,c1m7 = 1.0d-7)
      parameter(c1m9 = 1.0d-9,c1m10 = 1.0d-10,c1m12 = 1.0d-12)
      parameter(c1m13 = 1.0d-13,c1m15 = 1.0d-15)
      parameter(c1m18 = 1.0d-18,c1m21 = 1.0d-21,c1m24 = 1.0d-24)
      parameter(c1m36 = 1.0d-36,c1m38 = 1.0d-38)
!     electron mass from PDG, 2002
      parameter(pmap = 938.271998d0,pmae = .510998902d0)
      parameter(crade = 2.817940285d-15, clight = 2.99792458d8)
      integer icorr,idial,idptr,imod1,imod2,inorm,ipar,namp,ncor,nctype,&
     &ndimf,nmom,nmom1,nmom2,nord,nord1,nordf,nsix,nvar,nvar2,nvarf
      double precision dpmax,preda,weig1,weig2
      character*16 coel
      common/dial/preda,idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
      common/norf/nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
      common/tcorr/icorr,nctype,namp,nmom,nmom1,nmom2,weig1,weig2,dpmax,&
     &coel(10)
      dimension star(3,mbea),bcu(nbb,12)
      save
!-----------------------------------------------------------------------
*FOX  B D ;
*FOX  D V DA EXT TRACK NORD NVAR 6 ;
*FOX  D V DA INT S NORD NVAR ; D V DA INT SP NORD NVAR ;
*FOX  D V DA INT DD NORD NVAR ;
*FOX  D V DA INT DUM NORD NVAR 13 ;
*FOX  D V DA INT SX NORD NVAR ; D V DA INT SY NORD NVAR ;
*FOX  D V DA INT SEPX NORD NVAR ; D V DA INT SEPY NORD NVAR ;
*FOX  D V DA INT SEPX0 NORD NVAR ; D V DA INT SEPY0 NORD NVAR ;
*FOX  D V DA INT COSTH NORD NVAR ; D V DA INT SINTH NORD NVAR ;
*FOX  D V DA INT COSTHP NORD NVAR ; D V DA INT SINTHP NORD NVAR ;
*FOX  D V DA INT BBFX NORD NVAR ; D V DA INT BBFY NORD NVAR ;
*FOX  D V DA INT BBF0 NORD NVAR ; D V DA INT BBGX NORD NVAR ;
*FOX  D V DA INT BBGY NORD NVAR ;
*FOX  D V RE EXT STAR 3 MBEA ; D V RE EXT F ; D V RE EXT BCU NBB 12 ;
*FOX  D V IN EXT IBB ;
*FOX  D V RE INT HALF ; D V RE INT TWO ; D V RE INT FOUR ;
*FOX  D V RE INT ZERO ; D V RE INT ONE ; D V RE INT C1E3 ;
*FOX  D V RE INT SFAC ; D V RE INT CPHI ; D V RE INT CPHI2 ;
*FOX  D V IN INT JSLI ;
*FOX  D F RE DARE 1 ;
*FOX  E D ;
*FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
      do 2000 jsli=1,nsli
*FOX    S=(TRACK(5)-STAR(3,JSLI))*HALF ;
*FOX    SP=S/CPHI2 ;
*FOX    DUM(1)=BCU(IBB,1)+TWO*BCU(IBB,4)*SP+BCU(IBB,6)*SP*SP ;
*FOX    DUM(2)=BCU(IBB,2)+TWO*BCU(IBB,9)*SP+BCU(IBB,10)*SP*SP ;
*FOX    DUM(3)=BCU(IBB,3)+(BCU(IBB,5)+BCU(IBB,7))*SP+
*FOX    BCU(IBB,8)*SP*SP ;
*FOX    DUM(4)=DUM(1)-DUM(2) ;
*FOX    DUM(5)=DUM(4)*DUM(4)+FOUR*DUM(3)*DUM(3) ;
        if(ibbc.eq.1.and.(abs(dare(dum(4))).gt.pieni.and.               &
     &abs(dare(dum(5))).gt.pieni)) then
          ibbc1=1
*FOX    DUM(5)=SQRT(DUM(5)) ;
        else
          ibbc1=0
        endif
*FOX    SEPX0=TRACK(1)+TRACK(2)*S-STAR(1,JSLI) ;
*FOX    SEPY0=TRACK(3)+TRACK(4)*S-STAR(2,JSLI) ;
        if(ibbc1.eq.1) then
          sfac=one
!hr05     if(dare(dum(4)).lt.zero) sfac=-one
          if(dare(dum(4)).lt.zero) sfac=(-1d0*one)                       !hr05
*FOX    DUM(6)=SFAC*DUM(4)/DUM(5) ;
*FOX    DUM(7)=DUM(1)+DUM(2) ;
*FOX    COSTH=HALF*(ONE+DUM(6)) ;
          if(abs(dare(costh)).gt.pieni) then
*FOX    COSTH=SQRT(COSTH) ;
          else
*FOX    COSTH=ZERO ;
          endif
*FOX    SINTH=HALF*(ONE-DUM(6)) ;
          if(abs(dare(sinth)).gt.pieni) then
*FOX    SINTH=-SFAC*SQRT(SINTH) ;
          else
*FOX    SINTH=ZERO ;
          endif
          if(dare(dum(3)).lt.zero) then
*FOX    SINTH=-SINTH ;
          endif
*FOX    SY=SFAC*DUM(5) ;
*FOX    SX=(DUM(7)+SY)*HALF ;
*FOX    SY=(DUM(7)-SY)*HALF ;
*FOX    SEPX=SEPX0*COSTH+SEPY0*SINTH ;
*FOX    SEPY=-SEPX0*SINTH+SEPY0*COSTH ;
        else
*FOX    SX=DUM(1) ;
*FOX    SY=DUM(2) ;
*FOX    SEPX=SEPX0 ;
*FOX    SEPY=SEPY0 ;
        endif
        if(dare(sx).gt.dare(sy)) then
          call bbff(sepx,sepy,sx,sy,bbfx,bbfy,bbgx,bbgy)
        else
          call bbff(sepy,sepx,sy,sx,bbfy,bbfx,bbgy,bbgx)
        endif
*FOX    BBFX=F*BBFX ;
*FOX    BBFY=F*BBFY ;
*FOX    BBGX=F*BBGX ;
*FOX    BBGY=F*BBGY ;
        if(ibbc1.eq.1) then
*FOX    DUM(8)=TWO*(BCU(IBB,4)-BCU(IBB,9)+(BCU(IBB,6)-BCU(IBB,10))*SP) ;
*FOX    DUM(9)=BCU(IBB,5)+BCU(IBB,7)+TWO*BCU(IBB,8)*SP ;
*FOX    DUM(10)=(DUM(4)*DUM(8)+FOUR*DUM(3)*DUM(9))/
*FOX    DUM(5)/DUM(5)/DUM(5) ;
*FOX    DUM(11)=SFAC*(DUM(8)/DUM(5)-DUM(4)*DUM(10)) ;
*FOX    DUM(12)=BCU(IBB,4)+BCU(IBB,9)+(BCU(IBB,6)+BCU(IBB,10))*SP ;
*FOX    DUM(13)=SFAC*(DUM(4)*DUM(8)*HALF+TWO*DUM(3)*DUM(9))/DUM(5) ;
          if(abs(dare(costh)).gt.pieni) then
*FOX    COSTHP=DUM(11)/FOUR/COSTH ;
          else
*FOX    COSTHP=ZERO ;
          endif
          if(abs(dare(sinth)).gt.pieni) then
*FOX    SINTHP=-DUM(11)/FOUR/SINTH ;
          else
*FOX    SINTHP=ZERO ;
          endif
*FOX    TRACK(6)=TRACK(6)-
*FOX    (BBFX*(COSTHP*SEPX0+SINTHP*SEPY0)+
*FOX    BBFY*(-SINTHP*SEPX0+COSTHP*SEPY0)+
*FOX    BBGX*(DUM(12)+DUM(13))+BBGY*(DUM(12)-DUM(13)))/
*FOX    CPHI*HALF ;
*FOX    BBF0=BBFX ;
*FOX    BBFX=BBF0*COSTH-BBFY*SINTH ;
*FOX    BBFY=BBF0*SINTH+BBFY*COSTH ;
        else
*FOX    TRACK(6)=TRACK(6)-
*FOX    (BBGX*(BCU(IBB,4)+BCU(IBB,6)*SP)+
*FOX    BBGY*(BCU(IBB,9)+BCU(IBB,10)*SP))/CPHI ;
        endif
*FOX    TRACK(6)=TRACK(6)-(BBFX*(TRACK(2)-BBFX*HALF)+
*FOX    BBFY*(TRACK(4)-BBFY*HALF))*HALF ;
*FOX    TRACK(1)=TRACK(1)+S*BBFX ;
*FOX    TRACK(2)=TRACK(2)-BBFX ;
*FOX    TRACK(3)=TRACK(3)+S*BBFY ;
*FOX    TRACK(4)=TRACK(4)-BBFY ;
 2000 continue
!     DADAL AUTOMATIC INCLUSION
      return
      end
      subroutine boostif(sphi,cphi,tphi,salpha,calpha,track)
!-----------------------------------------------------------------------
!
!   Hirata's 6d beam-beam from BBC
!   DA VERSION for SIXTRACK courtesy Peter Leunissen
!   January 1999
!
!-----------------------------------------------------------------------
! BOOSTIF **************inverse boost ****************
      implicit none
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
      integer idaa
      double precision calpha,cphi,salpha,sphi,tphi
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,    &
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,    &
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele=1200,nblo=600,nper=16,nelb=140,nblz=20000,         &
     &nzfz = 300000,mmul = 20)
      parameter(nran = 2000000,ncom = 100,mran = 500,mpa = 6,nrco = 5,  &
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 99)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 350)
      double precision c180e0,c1e1,c1e12,c1e13,c1e15,c1e16,c1e2,c1e3,   &
     &c1e4,c1e6,c1m1,c1m7,c1m10,c1m12,c1m13,c1m15,c1m18,c1m2,c1m21,     &
     &c1m24,c1m3,c1m36,c1m38,c1m6,c1m9,c2e3,c4e3,crade,clight,four,half,&
     &one,pieni,pmae,pmap,three,two,zero
      parameter(pieni = 1d-38)
      parameter(zero = 0.0d0,half = 0.5d0,one = 1.0d0)
      parameter(two = 2.0d0,three = 3.0d0,four = 4.0d0)
      parameter(c1e1 = 1.0d1,c1e2 = 1.0d2,c1m2 = 1.0d-2)
      parameter(c1e3 = 1.0d3,c2e3 = 2.0d3,c4e3 = 4.0d3,c1e4 = 1.0d4)
      parameter(c1e12 = 1.0d12,c1e13 = 1.0d13,c1e15 = 1.0d15,c1e16 =    &
     &1.0d16)
      parameter(c180e0 = 180.0d0,c1e6 = 1.0d6)
      parameter(c1m1 = 1.0d-1,c1m3 = 1.0d-3,c1m6 = 1.0d-6,c1m7 = 1.0d-7)
      parameter(c1m9 = 1.0d-9,c1m10 = 1.0d-10,c1m12 = 1.0d-12)
      parameter(c1m13 = 1.0d-13,c1m15 = 1.0d-15)
      parameter(c1m18 = 1.0d-18,c1m21 = 1.0d-21,c1m24 = 1.0d-24)
      parameter(c1m36 = 1.0d-36,c1m38 = 1.0d-38)
!     electron mass from PDG, 2002
      parameter(pmap = 938.271998d0,pmae = .510998902d0)
      parameter(crade = 2.817940285d-15, clight = 2.99792458d8)
      integer icorr,idial,idptr,imod1,imod2,inorm,ipar,namp,ncor,nctype,&
     &ndimf,nmom,nmom1,nmom2,nord,nord1,nordf,nsix,nvar,nvar2,nvarf
      double precision dpmax,preda,weig1,weig2
      character*16 coel
      common/dial/preda,idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
      common/norf/nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
      common/tcorr/icorr,nctype,namp,nmom,nmom1,nmom2,weig1,weig2,dpmax,&
     &coel(10)
      save
!-----------------------------------------------------------------------
*FOX  B D ;
*FOX  D V DA EXT TRACK NORD NVAR 6 ; D V DA INT A1 NORD NVAR ;
*FOX  D V DA INT H1 NORD NVAR ; D V DA INT SQR1A NORD NVAR ;
*FOX  D V DA INT H1D NORD NVAR ; D V DA INT H1X NORD NVAR ;
*FOX  D V DA INT H1Y NORD NVAR ; D V DA INT H1Z NORD NVAR ;
*FOX  D V DA INT DET NORD NVAR ; D V DA INT X1 NORD NVAR ;
*FOX  D V DA INT Y1 NORD NVAR ; D V DA INT Z1 NORD NVAR ;
*FOX  D V RE INT ONE ; D V RE INT TWO ;
*FOX  D V RE EXT SPHI ; D V RE EXT CPHI ; D V RE EXT TPHI ;
*FOX  D V RE EXT SPHI2 ; D V RE EXT CPHI2 ; D V RE EXT TPHI2 ;
*FOX  D V RE EXT SALPHA ; D V RE EXT CALPHA ;
*FOX  E D ;
*FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
*FOX    H1D=SQRT((ONE+TRACK(6))*(ONE+TRACK(6))-
*FOX    TRACK(2)*TRACK(2)-TRACK(4)*TRACK(4)) ;
*FOX    H1X=TRACK(2)/H1D ;
*FOX    H1Y=TRACK(4)/H1D ;
*FOX    H1Z=ONE-(ONE+TRACK(6))/H1D ;
*FOX    H1=(TRACK(6)+ONE-SQRT((ONE+TRACK(6))*(ONE+TRACK(6))-
*FOX    TRACK(2)*TRACK(2)-TRACK(4)*TRACK(4)))*CPHI*CPHI ;
*FOX    DET=ONE/CPHI+TPHI*(H1X*CALPHA+H1Y*SALPHA-H1Z*SPHI) ;
*FOX    X1= TRACK(1)*(ONE/CPHI+SALPHA*(H1Y-H1Z*SALPHA*SPHI)*TPHI)
*FOX       +TRACK(3)*SALPHA*TPHI*(-H1X+H1Z*CALPHA*SPHI)
*FOX       -TRACK(5)*(CALPHA+H1Y*CALPHA*SALPHA*SPHI
*FOX       -H1X*SALPHA*SALPHA*SPHI)*TPHI ;
*FOX    Y1= TRACK(1)*CALPHA*TPHI*(-H1Y+H1Z*SALPHA*SPHI)
*FOX       +TRACK(3)*(ONE/CPHI+CALPHA*(H1X-H1Z*CALPHA*SPHI)*TPHI)
*FOX       -TRACK(5)*(SALPHA-H1Y*CALPHA*CALPHA*SPHI
*FOX       +H1X*CALPHA*SALPHA*SPHI)*TPHI ;
*FOX    Z1=-TRACK(1)*H1Z*CALPHA*SPHI
*FOX       -TRACK(3)*H1Z*SALPHA*SPHI
*FOX       +TRACK(5)*(ONE+H1X*CALPHA*SPHI+H1Y*SALPHA*SPHI) ;
*FOX    TRACK(1)=X1/DET ;
*FOX    TRACK(3)=Y1/DET ;
*FOX    TRACK(5)=Z1/DET ;
*FOX    TRACK(6)=TRACK(6)+CALPHA*SPHI*TRACK(2)
*FOX            +SALPHA*SPHI*TRACK(4) ;
*FOX    TRACK(2)=(TRACK(2)+CALPHA*SPHI*H1)*CPHI ;
*FOX    TRACK(4)=(TRACK(4)+SALPHA*SPHI*H1)*CPHI ;
!     DADAL AUTOMATIC INCLUSION
      return
      end
      subroutine bbff(sepx,sepy,sigxx,sigyy,bbfx,bbfy,bbgx,bbgy)
!-----------------------------------------------------------------------
!
!   Hirata's 6d beam-beam from BBC
!   DA VERSION for SIXTRACK courtesy Peter Leunissen
!   January 1999
!
!-----------------------------------------------------------------------
!**********************************************************************
! BBFF gives transverse (f_x and f_y) and longitudinal(g_x and g_y)
! beam-beam kicks except for the kinematical term (nr_e/\gamma)
! SIGXX is \Sigma
!**********************************************************************
      implicit none
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
      integer idaa
      double precision dare,hundred,sqrpi2
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,    &
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,    &
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele=1200,nblo=600,nper=16,nelb=140,nblz=20000,         &
     &nzfz = 300000,mmul = 20)
      parameter(nran = 2000000,ncom = 100,mran = 500,mpa = 6,nrco = 5,  &
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 99)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 350)
      double precision c180e0,c1e1,c1e12,c1e13,c1e15,c1e16,c1e2,c1e3,   &
     &c1e4,c1e6,c1m1,c1m7,c1m10,c1m12,c1m13,c1m15,c1m18,c1m2,c1m21,     &
     &c1m24,c1m3,c1m36,c1m38,c1m6,c1m9,c2e3,c4e3,crade,clight,four,half,&
     &one,pieni,pmae,pmap,three,two,zero
      parameter(pieni = 1d-38)
      parameter(zero = 0.0d0,half = 0.5d0,one = 1.0d0)
      parameter(two = 2.0d0,three = 3.0d0,four = 4.0d0)
      parameter(c1e1 = 1.0d1,c1e2 = 1.0d2,c1m2 = 1.0d-2)
      parameter(c1e3 = 1.0d3,c2e3 = 2.0d3,c4e3 = 4.0d3,c1e4 = 1.0d4)
      parameter(c1e12 = 1.0d12,c1e13 = 1.0d13,c1e15 = 1.0d15,c1e16 =    &
     &1.0d16)
      parameter(c180e0 = 180.0d0,c1e6 = 1.0d6)
      parameter(c1m1 = 1.0d-1,c1m3 = 1.0d-3,c1m6 = 1.0d-6,c1m7 = 1.0d-7)
      parameter(c1m9 = 1.0d-9,c1m10 = 1.0d-10,c1m12 = 1.0d-12)
      parameter(c1m13 = 1.0d-13,c1m15 = 1.0d-15)
      parameter(c1m18 = 1.0d-18,c1m21 = 1.0d-21,c1m24 = 1.0d-24)
      parameter(c1m36 = 1.0d-36,c1m38 = 1.0d-38)
!     electron mass from PDG, 2002
      parameter(pmap = 938.271998d0,pmae = .510998902d0)
      parameter(crade = 2.817940285d-15, clight = 2.99792458d8)
      parameter(sqrpi2 = 3.544907701811032d0,hundred = 100d0)
      integer icorr,idial,idptr,imod1,imod2,inorm,ipar,namp,ncor,nctype,&
     &ndimf,nmom,nmom1,nmom2,nord,nord1,nordf,nsix,nvar,nvar2,nvarf
      double precision dpmax,preda,weig1,weig2
      character*16 coel
      common/dial/preda,idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
      common/norf/nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
      common/tcorr/icorr,nctype,namp,nmom,nmom1,nmom2,weig1,weig2,dpmax,&
     &coel(10)
      save
!-----------------------------------------------------------------------
*FOX  B D ;
*FOX  D V DA EXT SEPX NORD NVAR ; D V DA EXT SEPY NORD NVAR ;
*FOX  D V DA EXT SIGXX NORD NVAR ; D V DA EXT SIGYY NORD NVAR ;
*FOX  D V DA EXT BBFX NORD NVAR ; D V DA EXT BBFY NORD NVAR ;
*FOX  D V DA EXT BBGX NORD NVAR ; D V DA EXT BBGY NORD NVAR ;
*FOX  D V DA INT SIGXY NORD NVAR ; D V DA INT EXPFAC NORD NVAR ;
*FOX  D V DA INT X NORD NVAR ; D V DA INT FAC NORD NVAR ;
*FOX  D V DA INT FAC2 NORD NVAR ; D V DA INT CONST NORD NVAR ;
*FOX  D V DA INT ARG1X NORD NVAR ; D V DA INT ARG1Y NORD NVAR ;
*FOX  D V DA INT ARG2X NORD NVAR ; D V DA INT ARG2Y NORD NVAR ;
*FOX  D V DA INT WX1 NORD NVAR ; D V DA INT WY1 NORD NVAR ;
*FOX  D V DA INT WX2 NORD NVAR ; D V DA INT WY2 NORD NVAR ;
*FOX  D V DA INT COMFAC NORD NVAR ; D V DA INT COMFAC2 NORD NVAR ;
*FOX  D V RE INT ZERO ; D V RE INT HALF ; D V RE INT ONE ;
*FOX  D V RE INT TWO ; D V RE INT FOUR ; D V RE INT HUNDRED ;
*FOX  D V RE INT SQRPI2 ;
*FOX  D F RE DARE 1 ;
*FOX  E D ;
*FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
      if(dare(sigxx).eq.dare(sigyy)) then
*FOX    X=SEPX*SEPX+SEPY*SEPY ;
        if(abs(dare(sigxx)+dare(sigyy)).gt.pieni) then
*FOX      CONST=X/(SIGXX+SIGYY) ;
        else
*FOX      CONST=ZERO ;
        endif
*FOX    EXPFAC=EXP(-CONST) ;
        if(abs(dare(x)).gt.pieni) then
*FOX      BBFX=TWO*SEPX*(ONE-EXPFAC)/X ;
*FOX      BBFY=TWO*SEPY*(ONE-EXPFAC)/X ;
*FOX      COMFAC=-SEPX*BBFX+SEPY*BBFY ;
          if(dare(sigxx).lt.zero) then
*FOX        SIGXX=-SIGXX ;
          endif
          if(dare(sigyy).lt.zero) then
*FOX        SIGYY=-SIGYY ;
          endif
*FOX      COMFAC2=(SIGXX+SIGYY)*(SIGXX+SIGYY) ;
*FOX      BBGX=(COMFAC+FOUR*SEPX*SEPX*CONST/X*EXPFAC)/(TWO*X) ;
*FOX      BBGY=(-COMFAC+FOUR*SEPY*SEPY*CONST/X*EXPFAC)/(TWO*X) ;
        else
*FOX      BBFX=ZERO ;
*FOX      BBFY=ZERO ;
*FOX      BBGX=ZERO ;
*FOX      BBGY=ZERO ;
        endif
      else
*FOX    X=SEPX*SEPX/SIGXX+SEPY*SEPY/SIGYY ;
*FOX    FAC2=TWO*(SIGXX-SIGYY) ;
        if(dare(sigxx).lt.dare(sigyy)) then
*FOX      FAC2=TWO*(SIGYY-SIGXX) ;
        endif
*FOX    FAC=SQRT(FAC2) ;
*FOX    CONST=SQRPI2/FAC ;
*FOX    SIGXY=SQRT(SIGXX/SIGYY) ;
*FOX    ARG1X=(SEPX/FAC) ;
        if(dare(sepx).lt.zero) then
*FOX      ARG1X=-(SEPX/FAC) ;
        endif
*FOX    ARG1Y=(SEPY/FAC) ;
        if(dare(sepy).lt.zero) then
*FOX      ARG1Y=-(SEPY/FAC) ;
        endif
        call errff(arg1x,arg1y,wy1,wx1)
        if(dare(x).lt.hundred) then
*FOX      EXPFAC=EXP(-X*HALF) ;
*FOX      ARG2X=ARG1X/SIGXY ;
*FOX      ARG2Y=ARG1Y*SIGXY ;
          call errff(arg2x,arg2y,wy2,wx2)
*FOX      BBFX=CONST*(WX1-EXPFAC*WX2) ;
*FOX      BBFY=CONST*(WY1-EXPFAC*WY2) ;
          if(dare(sepx).lt.zero) then
*FOX        BBFX=-BBFX ;
          endif
          if(dare(sepy).lt.zero) then
*FOX        BBFY=-BBFY ;
          endif
*FOX      COMFAC=SEPX*BBFX+SEPY*BBFY ;
*FOX      BBGX=-(COMFAC+TWO*(EXPFAC/SIGXY-ONE))/FAC2 ;
*FOX      BBGY= (COMFAC+TWO*(EXPFAC*SIGXY-ONE))/FAC2 ;
        else
*FOX      BBFX=CONST*WX1 ;
*FOX      BBFY=CONST*WY1 ;
          if(dare(sepx).lt.zero) then
*FOX        BBFX=-BBFX ;
          endif
          if(dare(sepy).lt.zero) then
*FOX        BBFY=-BBFY ;
          endif
*FOX      COMFAC=SEPX*BBFX+SEPY*BBFY ;
*FOX      BBGX=-(COMFAC-TWO)/FAC2 ;
*FOX      BBGY= -BBGX ;
        endif
      endif
!     DADAL AUTOMATIC INCLUSION
      return
      end
      subroutine clorda(nn,idummy,am)
!-----------------------------------------------------------------------
!  CALCULATION OF THE SIX-DIMENSIONAL CLOSED ORBIT
!-----------------------------------------------------------------------
      implicit none
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
      integer i,i4,icheck,ii,j,j4,k,l,ll,nd2,nn
      double precision am,cloc,cor,coro,dc,dd,dlo,xx
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,    &
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,    &
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele=1200,nblo=600,nper=16,nelb=140,nblz=20000,         &
     &nzfz = 300000,mmul = 20)
      parameter(nran = 2000000,ncom = 100,mran = 500,mpa = 6,nrco = 5,  &
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 99)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 350)
      double precision c180e0,c1e1,c1e12,c1e13,c1e15,c1e16,c1e2,c1e3,   &
     &c1e4,c1e6,c1m1,c1m7,c1m10,c1m12,c1m13,c1m15,c1m18,c1m2,c1m21,     &
     &c1m24,c1m3,c1m36,c1m38,c1m6,c1m9,c2e3,c4e3,crade,clight,four,half,&
     &one,pieni,pmae,pmap,three,two,zero
      parameter(pieni = 1d-38)
      parameter(zero = 0.0d0,half = 0.5d0,one = 1.0d0)
      parameter(two = 2.0d0,three = 3.0d0,four = 4.0d0)
      parameter(c1e1 = 1.0d1,c1e2 = 1.0d2,c1m2 = 1.0d-2)
      parameter(c1e3 = 1.0d3,c2e3 = 2.0d3,c4e3 = 4.0d3,c1e4 = 1.0d4)
      parameter(c1e12 = 1.0d12,c1e13 = 1.0d13,c1e15 = 1.0d15,c1e16 =    &
     &1.0d16)
      parameter(c180e0 = 180.0d0,c1e6 = 1.0d6)
      parameter(c1m1 = 1.0d-1,c1m3 = 1.0d-3,c1m6 = 1.0d-6,c1m7 = 1.0d-7)
      parameter(c1m9 = 1.0d-9,c1m10 = 1.0d-10,c1m12 = 1.0d-12)
      parameter(c1m13 = 1.0d-13,c1m15 = 1.0d-15)
      parameter(c1m18 = 1.0d-18,c1m21 = 1.0d-21,c1m24 = 1.0d-24)
      parameter(c1m36 = 1.0d-36,c1m38 = 1.0d-38)
!     electron mass from PDG, 2002
      parameter(pmap = 938.271998d0,pmae = .510998902d0)
      parameter(crade = 2.817940285d-15, clight = 2.99792458d8)
      integer iav,ibb6d,ibbc,ibeco,ibidu,ibtyp,ic,icext,icextal,iclo,   &
     &iclo6,iclo6r,icode,icoe,icomb,icomb0,iconv,icow,icr,idam,idfor,   &
     &idis,idp,ierro,iffw,ifh,iicav,il,ilin,imad,imbb,                  &
     &imc,imtr,iorg,iout,                                               &
     &ipos,ipr,iprint,ipt,iq,iqmod,iqmod6,iratioe,ird,ire,ires,         &
     &irew,irm,irmod2,ise,ise1,ise2,ise3,isea,iskew,iskip,istw,         &
     &isub,itco,itcro,itf,ithick,ition,itionc,itqv,its6d,iu,iver,ivox,  &
     &ivoz,iwg,ixcav,izu0,kanf,kp,kpa,kwtype,kz,lhc,m21,m22,m23,mblo,   &
     &mbloz,mcut,mel,mesa,mmac,mout2,mp,mper,mstr,msym,mtyp,mzu,napx,   &
     &napxo,nbeam,nch,ncororb,ncorrep,ncorru,ncy,ndafi,nde,nhcorr,      &
     &nhmoni,niu,nlin,nmu,npp,nprint,nqc,nre,nrr,nskew,                 &
     &nstart,nstop,nt,nta,ntco,nte,ntwin,nu,numl,numlr,nur,nvcorr,      &
     &nvmoni,nwr, nturn1, nturn2, nturn3, nturn4,numlcp,numlmax,nnuml
      double precision a,ak0,aka,alfx,alfz,amp0,aper,apx,apz,ape,bbcu,  &
     &bclorb,beamoff,benkc,benki,betac,betam,betx,betz,bk0,bka,bl1,bl2, &
     &clo6,clobeam,clop6,cma1,cma2,cotr,crad,de0,dech,ded,dfft,         &
     &di0,dip0,dki,dkq,dma,dmap,dphix,dphiz,dppoff,dpscor,dqq,dres,dsi, &
     &dsm0,dtr,e0,ed,ej,ejf,ek,el,elbe,emitx,emity,emitz,extalign,      &
     &exterr,eui,euii,gammar,hsy,hsyc,pac,pam,parbe,parbe14,partnum,    &
     &phas,phas0,phasc,pi,pi2,pisqrt,pma,ptnfac,qs,qw0,qwsk,qx0,qxt,qz0,&
     &qzt,r00,rad,rat,ratio,ratioe,rrtr,rtc,rts,rvf,                    &
     &sigcor,sige,sigma0,sigman,sigman2,sigmanq,sigmoff,sigz,sm,ta,tam1,&
     &tam2,tiltc,tilts,tlen,totl,track6d,xpl,xrms,zfz,zpl,zrms,wirel,   &
     &acdipph, crabph, bbbx, bbby, bbbs,                                &
     &crabph2, crabph3, crabph4
      real hmal
      character*16 bez,bezb,bezr,erbez,bezl
      character*80 toptit,sixtit,commen
      common/erro/ierro,erbez
      common/kons/pi,pi2,pisqrt,rad
      common/str /il,mper,mblo,mbloz,msym(nper),kanf,iu,ic(nblz)
      common/ell /ed(nele),el(nele),ek(nele),sm(nele),kz(nele),kp(nele)
      common/bbb /bbbx(nele), bbby(nele), bbbs(nele)
      common/pla /xpl(nele),xrms(nele),zpl(nele),zrms(nele)
      common/str2 /mel(nblo),mtyp(nblo,nelb),mstr(nblo)
      common/mat/a(nele,2,6),bl1(nblo,2,6),bl2(nblo,2,6)
      common/syos2/rvf(mpa)
      common/tra1/rat,idfor,napx,napxo,numl,niu(2),numlr,nde(2),nwr(4), &
     &ird,imc,irew,ntwin,iclo6,iclo6r,iver,ibidu,numlcp,numlmax,nnuml
      common/syn/qs,e0,pma,ej(mpa),ejf(mpa),phas0,phas,hsy(3),          &
     &crad,hsyc(nele),phasc(nele),dppoff,sigmoff(nblz),tlen,            &
     &iicav,itionc(nele),ition,idp,ncy,ixcav
      common/corcom/dpscor,sigcor,icode,idam,its6d
      common/multi/bk0(nele,mmul),ak0(nele,mmul),                       &
     &bka(nele,mmul),aka(nele,mmul)
      common/mult1/benki,benkc(nele),r00(nele),irm(nele),nmu(nele)
      common/rand0/zfz(nzfz),iorg,mzu(nblz),bezr(3,nele),izu0,mmac,mcut
      common/rand1/exterr(nblz,40),extalign(nblz,3),tiltc(nblz),        &
     &tilts(nblz),mout2,icext(nblz),icextal(nblz)
      common/beo /aper(2),di0(2),dip0(2),ta(6,6)
      common/clo/dma,dmap,dkq,dqq,de0,ded,dsi,dech,dsm0,itco,itcro,itqv,&
     &iout
      common/qmodi/qw0(3),amp0,iq(3),iqmod,kpa(nele),iqmod6
      common/linop/bez(nele),elbe(nblo),bezb(nblo),ilin,nt,iprint,      &
     &ntco,eui,euii,nlin,bezl(nele)
      common/cororb/betam(nmon1,2),pam(nmon1,2),betac(ncor1,2),         &
     &pac(ncor1,2),bclorb(nmon1,2),nhmoni,nhcorr,nvmoni,nvcorr,         &
     &ncororb(nele)
      common/apert/apx(nele),apz(nele),ape(3,nele)
      common/clos/sigma0(2),iclo,ncorru,ncorrep
      common/combin/icomb0(20),icomb(ncom,20),ratio(ncom,20),           &
     &ratioe(nele),iratioe(nele),icoe
      common/seacom/ise,mesa,mp,m21,m22,m23,ise1,ise2,ise3,isea(nele)
      common/subres/qxt,qzt,tam1,tam2,isub,nta,nte,ipt,totl
      common/secom/rtc(9,18,10,5),rts(9,18,10,5),ire(12),ipr(5),irmod2
      common/secom1/dtr(10),nre,nur,nch,nqc,npp,nrr(5),nu(5)
      common/postr/dphix,dphiz,qx0,qz0,dres,dfft,cma1,cma2,             &
     &nstart,nstop,iskip,iconv,imad
      common/posti1/ipos,iav,iwg,ivox,ivoz,ires,ifh,toptit(5)
      common/posti2/kwtype,itf,icr,idis,icow,istw,iffw,nprint,ndafi
      common/skew/qwsk(2),betx(2),betz(2),alfx(2),alfz(2),iskew,nskew(6)
      common/pawc/hmal(nplo)
      common/tit/sixtit,commen,ithick
      common/co6d/clo6(3),clop6(3)
      common/dkic/dki(nele,3)
      common/beam/sigman(2,nbb),sigman2(2,nbb),sigmanq(2,nbb),          &
     &clobeam(6,nbb),beamoff(6,nbb),parbe(nele,5),track6d(6,npart),     &
     &ptnfac(nele),sigz,sige,partnum,parbe14,emitx,emity,emitz,gammar,  &
     &nbeam,ibbc,ibeco,ibtyp,lhc
      common/trom/ cotr(ntr,6),rrtr(ntr,6,6),imtr(nele)
      common/bb6d/ bbcu(nbb,12),ibb6d,imbb(nblz)
      common/wireco/ wirel(nele)
      common/acdipco/ acdipph(nele), nturn1(nele), nturn2(nele),        &
     &nturn3(nele), nturn4(nele)
      common/crabco/ crabph(nele),crabph2(nele),                        &
     &crabph3(nele),crabph4(nele)
      integer idz,itra
      double precision a2,al,as,at,chi0,chid,dp1,dps,exz,sigm
      common/syos/as(6,2,2,nele),at(6,2,2,nele),a2(6,2,2,nele),         &
     &al(6,2,2,nele),sigm(mpa),dps(mpa),idz(2)
      common/anf/chi0,chid,exz(2,6),dp1,itra
      integer ichrom,is
      double precision alf0,amp,bet0,clo,clop,cro,x,y
      common/tra/x(mpa,2),y(mpa,2),amp(2),bet0(2),alf0(2),clo(2),clop(2)
      common/chrom/cro(2),is(2),ichrom
      double precision aml6,edcor
      common/sixdim/aml6(6,6),edcor(2)
      integer icorr,idial,idptr,imod1,imod2,inorm,ipar,namp,ncor,nctype,&
     &ndimf,nmom,nmom1,nmom2,nord,nord1,nordf,nsix,nvar,nvar2,nvarf
      double precision dpmax,preda,weig1,weig2
      character*16 coel
      common/dial/preda,idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
      common/norf/nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
      common/tcorr/icorr,nctype,namp,nmom,nmom1,nmom2,weig1,weig2,dpmax,&
     &coel(10)
      integer idummy(nn)
      character*6 chp(3),chd(3)
      dimension xx(6),dlo(6),cloc(6),dd(6),dc(6),am(nn,nn)
      integer nerror
      save
!-----------------------------------------------------------------------
      nd2=2*ndimf
      write(*,10010) nd2
      do l=1,nd2
        xx(l)=zero
        cloc(l)=zero
        dd(l)=zero
        dc(l)=zero
        dlo(l)=zero
        do i=1,nd2
          am(l,i)=zero
        enddo
      enddo
      chp(1)=' CLOX '
      chp(2)=' CLOY '
      chp(3)=' CLOS '
      chd(1)='  D-X '
      chd(2)='  D-Y '
      chd(3)='  D-S '
      cor=zero
      coro=1d38
      ii=0
      if(ndimf.eq.3) then
        do l=1,2
          ll=2*l
          cloc(ll-1)=clo6(l)
          cloc(ll)=clop6(l)
        enddo
        cloc(5)=clo6(3)
        cloc(6)=clop6(3)
        if(abs(dppoff).gt.pieni) cloc(6)=dppoff
      else
        do l=1,ndimf
          ll=2*l
          cloc(ll-1)=clo(l)
          cloc(ll)=clop(l)
        enddo
        do l=ndimf+1,3
          ll=2*l
          cloc(ll-1)=zero
          cloc(ll)=zero
        enddo
        cloc(6)=dps(1)
      endif
      do 80 ii=1,itco
        do l=1,2
          ll=2*l
          x(1,l)=cloc(ll-1)
          y(1,l)=cloc(ll)
        enddo
        sigm(1)=cloc(5)
        dps(1)=cloc(6)
        call umlauda
        do i4=1,nd2
          do j4=1,nd2
            am(i4,j4)=aml6(i4,j4)
          enddo
        enddo
        call dinv(nd2,am,nd2,idummy,nerror)
        if(nerror.ne.0) write(*,*) ' ATTENTION, MATRIX SINGULAR '
        if(ndimf.eq.3) then
          do l=1,2
            ll=2*l
            xx(ll-1)=x(1,l)
            xx(ll)=y(1,l)
          enddo
          xx(5)=sigm(1)
          xx(6)=dps(1)
        else
          do l=1,ndimf
            ll=2*l
            xx(ll-1)=x(1,l)
            xx(ll)=y(1,l)
          enddo
          do l=ndimf+1,3
            ll=2*l
            xx(ll-1)=zero
            xx(ll)=zero
          enddo
        endif
        do l=1,nd2
          dd(l)=cloc(l)-xx(l)
          dc(l)=abs(dd(l))
          if(l.eq.5) dc(5)=dc(5)*c1m2
        enddo
        icheck=0
        do l=1,ndimf
          ll=2*l
          if(dc(ll-1).gt.dma) icheck=1
          if(dc(ll).gt.dmap) icheck=1
        enddo
        if(icheck.eq.0) goto 90
        do k=1,nd2
          dlo(k)=zero
          do j=1,nd2
            dlo(k)=am(k,j)*dd(j)+dlo(k)
          enddo
          if(abs(dppoff).gt.pieni) dlo(6)=zero
        enddo
        write(*,10020)
        cor=zero
        do l=1,ndimf
          ll=2*l
          write(*,10060) chp(l),cloc(ll-1),cloc(ll)
!hr06     cor=cor+dc(ll-1)*dc(ll-1)
          cor=cor+dc(ll-1)**2                                            !hr06
        enddo
        cor=sqrt(cor)
        if(ii.eq.1.or.cor.lt.coro) then
          coro=cor
          do l=1,nd2
            cloc(l)=cloc(l)+dlo(l)
          enddo
          if(ii.ne.itco) then
            write(*,10030)
            do l=1,ndimf
              ll=2*l
              write(*,10060) chp(l),cloc(ll-1),cloc(ll)
            enddo
            write(*,10080) ii,cor
          endif
        else
          write(*,10040) nd2,ii
          goto 91
        endif
 80   continue
      write(*,10000) itco
      ii=itco
 90   continue
      if(ii.ne.itco) then
        do k=1,nd2
          dlo(k)=zero
          do j=1,nd2
            dlo(k)=am(k,j)*dd(j)+dlo(k)
          enddo
          if(abs(dppoff).gt.pieni) dlo(6)=zero
        enddo
        write(*,10020)
        cor=zero
        do l=1,ndimf
          ll=2*l
          write(*,10060) chp(l),cloc(ll-1),cloc(ll)
!hr06     cor=cor+dc(ll-1)*dc(ll-1)
          cor=cor+dc(ll-1)**2                                            !hr06
        enddo
        cor=sqrt(cor)
        if(cor.lt.coro) then
          coro=cor
          do l=1,nd2
            cloc(l)=cloc(l)+dlo(l)
          enddo
          write(*,10030)
          do l=1,ndimf
            ll=2*l
            write(*,10060) chp(l),cloc(ll-1),cloc(ll)
          enddo
          write(*,10080) ii,cor
        else
          write(*,10040) nd2,ii
          goto 91
        endif
        do l=1,2
          ll=2*l
          x(1,l)=cloc(ll-1)
          y(1,l)=cloc(ll)
        enddo
        sigm(1)=cloc(5)
        dps(1)=cloc(6)
        call umlauda
        do i4=1,nd2
          do j4=1,nd2
            am(i4,j4)=aml6(i4,j4)
          enddo
        enddo
        call dinv(nd2,am,nd2,idummy,nerror)
        if(nerror.ne.0) write(*,*) ' ATTENTION, MATRIX SINGULAR '
        if(ndimf.eq.3) then
          do l=1,2
            ll=2*l
            xx(ll-1)=x(1,l)
            xx(ll)=y(1,l)
          enddo
          xx(5)=sigm(1)
          xx(6)=dps(1)
        else
          do l=1,ndimf
            ll=2*l
            xx(ll-1)=x(1,l)
            xx(ll)=y(1,l)
          enddo
          do l=ndimf+1,3
            ll=2*l
            xx(ll-1)=zero
            xx(ll)=zero
          enddo
        endif
        do l=1,nd2
          dc(l)=abs(cloc(l)-xx(l))
          if(l.eq.5) dc(5)=dc(5)*c1m2
        enddo
      endif
      write(*,10050) nd2,ii
      cor=zero
      do l=1,ndimf
        ll=2*l
        write(*,10070) chp(l),cloc(ll-1),cloc(ll),                      &
     &chd(l),dc(ll-1),dc(ll)
!hr06   cor=cor+dc(ll-1)*dc(ll-1)
        cor=cor+dc(ll-1)**2                                              !hr06
      enddo
      cor=sqrt(cor)
      write(*,10080) ii,cor
 91   continue
      if(ndimf.eq.3) then
        do l=1,2
          ll=2*l
          clo6(l)=cloc(ll-1)
          clop6(l)=cloc(ll)
        enddo
        clo6(3)=cloc(5)
        clop6(3)=cloc(6)
      else
        do l=1,ndimf
          ll=2*l
          clo(l)=cloc(ll-1)
          clop(l)=cloc(ll)
        enddo
      endif
!-----------------------------------------------------------------------
      return
10000 format(t10,'DA CLOSED ORBIT CALCULATION'/ t10,                    &
     &'MAXIMUM NUMBER OF ITERATIONS ACHIEVED--->',2x,i4/ t10,           &
     &'PROCEDURE MAY NOT HAVE CONVERGED')
10010 format(/131('-')/t10,'ENTERING ',i1,                              &
     &'-D DA CLOSED ORBIT CALCULATION'/)
10020 format(5x,'---- closed orbit before correction----')
10030 format(5x,'---- after DA correction----')
10040 format(/5x,'NO IMPROVEMENT OF ',i1,'-D DA CLOSED ORBIT ',         &
     &'CALCULATION IN ITERATION: ',i4/)
10050 format(t5,'SUCCESSFULL END OF ',i1,'-D DA CLOSED ORBIT ',         &
     &'CALCULATION IN ITERATION: ',i4/)
10060 format(5x,a6,1p,2(1x,g16.9))
10070 format(5x,a6,1p,2(1x,g16.9)/5x,a6,1p,2(1x,g16.9))
10080 format(5x,' ITERAT.=',i3,' ACCURACY=',d13.6/)
      end
      subroutine mydaini(ncase,nnord,nnvar,nndim,nnvar2,nnord1)
!-----------------------------------------------------------------------
!  CALCULATION OF THE 4-DIMENSIONAL CLOSED ORBIT INCLUDING DELTA
!-----------------------------------------------------------------------
      implicit none
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
      integer idummy,ncase,ndimfo,ndpt,nis,nndim,                       &
     &nnord,nnord1,nnvar,nnvar2,nord1o,nordo,nvar2o,nvaro
      double precision am
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,    &
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,    &
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele=1200,nblo=600,nper=16,nelb=140,nblz=20000,         &
     &nzfz = 300000,mmul = 20)
      parameter(nran = 2000000,ncom = 100,mran = 500,mpa = 6,nrco = 5,  &
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 99)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 350)
      integer icorr,idial,idptr,imod1,imod2,inorm,ipar,namp,ncor,nctype,&
     &ndimf,nmom,nmom1,nmom2,nord,nord1,nordf,nsix,nvar,nvar2,nvarf
      double precision dpmax,preda,weig1,weig2
      character*16 coel
      common/dial/preda,idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
      common/norf/nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
      common/tcorr/icorr,nctype,namp,nmom,nmom1,nmom2,weig1,weig2,dpmax,&
     &coel(10)
      integer idao,iscrri
      integer          iscrda
      double precision rscrri
      common/dascr/iscrda(100),rscrri(100),iscrri(100),idao
      integer ichromc,ilinc,iqmodc
      double precision clon,chromc,corr,wxys
      common/correct/ corr(3,3),chromc(2),wxys(3),clon(6),iqmodc,       &
     &ichromc,ilinc
      dimension am(6,6),idummy(6)
      save
!-----------------------------------------------------------------------
      if(nndim.lt.2.or.nndim.gt.3) call prror(95)
!--------------------
      nordo=nord
      nvaro=nvar
      ndimfo=ndimf
      nvar2o=nvar2
      nord1o=nord1
!--------------------
      nord=nnord
      nvar=nnvar
      ndimf=nndim
      nvar2=nnvar2
      nord1=nnord1
!--------------------
      ndpt=0
      nis=0
!--------------------
      call daeps(preda)
      call idprset(-102)
      call lieinit(nord,nvar,ndimf,ndpt,0,nis)
      write(*,10000) nord,nvar,ndimf
      call daall(iscrda,100,'$$IS      ',nord,nvar)
!--closed orbit
      if(ncase.eq.1) call clorda(2*ndimf,idummy,am)
!--tune variation
      if(ncase.eq.2) call umlauda
      iqmodc=0
      ichromc=0
      ilinc=0
      call dadal(iscrda,100)
!--------------------
      nord=nordo
      nvar=nvaro
      nvar2=nvar2o
      ndimf=ndimfo
      nord1=nord1o
!-----------------------------------------------------------------------
10000 format(/131('-')/10x,'DA INITIALIZATION: ORDER = ',i2,            &
     &', # of VARIABLES = ',i2,', DIMENSION = ',i2/)
      return
      end
