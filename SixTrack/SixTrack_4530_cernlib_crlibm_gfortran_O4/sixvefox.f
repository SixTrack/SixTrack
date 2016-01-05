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
      double precision al,as,chi0,chid,dp1,dps,exz,sigm
      common/syos/as(6,2,npart,nele),al(6,2,npart,nele),sigm(mpa),      &
     &dps(mpa),idz(2)
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
*FOX{
      INTEGER SIGMDA  
      INTEGER DPDA    
      INTEGER DPDA1   
      INTEGER RV      
      INTEGER XX      (2)
      INTEGER YY      (2)
      INTEGER EJ1     
      INTEGER EJF1    
      INTEGER ALDA    (2,6)
      INTEGER ASDA    (2,6)
      INTEGER ALDAQ   (2,6)
      INTEGER ASDAQ   (2,6)
      INTEGER SMIDA   (MCOR)
      INTEGER FOKQ    
      INTEGER WFHI    
      INTEGER DPD     
      INTEGER DPSQ    
      INTEGER FOK     
      INTEGER RHO     
      INTEGER FOK1    
      INTEGER SM1     
      INTEGER SM2     
      INTEGER SM3     
      INTEGER SM4     
      INTEGER SM5     
      INTEGER SM6     
      INTEGER SM12    
      INTEGER SM23    
      INTEGER AS3     
      INTEGER AS4     
      INTEGER AS6     
      INTEGER SI      
      INTEGER CO      
      INTEGER G       
      INTEGER GL      
      INTEGER SIQ     
      INTEGER RHOC    
      INTEGER HI      
      INTEGER FI      
      INTEGER AEK     
      INTEGER HI1     
      INTEGER HP      
      INTEGER HM      
      INTEGER HC      
      INTEGER HS      
      INTEGER FOKC    
      INTEGER WF      
      INTEGER AFOK    
      INTEGER RHOI    
      INTEGER WFA     
      INTEGER ISCRDA, ISCRRI,IDAO
      DOUBLE PRECISION RSCRRI
      COMMON/DASCR/ISCRDA(100),RSCRRI(100),ISCRRI(100),IDAO
      if(1.eq.1) then                                                 
         CALL DAKEY('FOX V2.1')
         CALL DAALL(FOKQ    ,1,'FOKQ      ',NORD,NVAR)
         CALL DAALL(WFHI    ,1,'WFHI      ',NORD,NVAR)
         CALL DAALL(DPD     ,1,'DPD       ',NORD,NVAR)
         CALL DAALL(DPSQ    ,1,'DPSQ      ',NORD,NVAR)
         CALL DAALL(FOK     ,1,'FOK       ',NORD,NVAR)
         CALL DAALL(RHO     ,1,'RHO       ',NORD,NVAR)
         CALL DAALL(FOK1    ,1,'FOK1      ',NORD,NVAR)
         CALL DAALL(SM1     ,1,'SM1       ',NORD,NVAR)
         CALL DAALL(SM2     ,1,'SM2       ',NORD,NVAR)
         CALL DAALL(SM3     ,1,'SM3       ',NORD,NVAR)
         CALL DAALL(SM4     ,1,'SM4       ',NORD,NVAR)
         CALL DAALL(SM5     ,1,'SM5       ',NORD,NVAR)
         CALL DAALL(SM6     ,1,'SM6       ',NORD,NVAR)
         CALL DAALL(SM12    ,1,'SM12      ',NORD,NVAR)
         CALL DAALL(SM23    ,1,'SM23      ',NORD,NVAR)
         CALL DAALL(AS3     ,1,'AS3       ',NORD,NVAR)
         CALL DAALL(AS4     ,1,'AS4       ',NORD,NVAR)
         CALL DAALL(AS6     ,1,'AS6       ',NORD,NVAR)
         CALL DAALL(SI      ,1,'SI        ',NORD,NVAR)
         CALL DAALL(CO      ,1,'CO        ',NORD,NVAR)
         CALL DAALL(G       ,1,'G         ',NORD,NVAR)
         CALL DAALL(GL      ,1,'GL        ',NORD,NVAR)
         CALL DAALL(SIQ     ,1,'SIQ       ',NORD,NVAR)
         CALL DAALL(RHOC    ,1,'RHOC      ',NORD,NVAR)
         CALL DAALL(HI      ,1,'HI        ',NORD,NVAR)
         CALL DAALL(FI      ,1,'FI        ',NORD,NVAR)
         CALL DAALL(AEK     ,1,'AEK       ',NORD,NVAR)
         CALL DAALL(HI1     ,1,'HI1       ',NORD,NVAR)
         CALL DAALL(HP      ,1,'HP        ',NORD,NVAR)
         CALL DAALL(HM      ,1,'HM        ',NORD,NVAR)
         CALL DAALL(HC      ,1,'HC        ',NORD,NVAR)
         CALL DAALL(HS      ,1,'HS        ',NORD,NVAR)
         CALL DAALL(FOKC    ,1,'FOKC      ',NORD,NVAR)
         CALL DAALL(WF      ,1,'WF        ',NORD,NVAR)
         CALL DAALL(AFOK    ,1,'AFOK      ',NORD,NVAR)
         CALL DAALL(RHOI    ,1,'RHOI      ',NORD,NVAR)
         CALL DAALL(WFA     ,1,'WFA       ',NORD,NVAR)
      ENDIF
      IDAA = IDAO
*FOX}
!-----------------------------------------------------------------------
*FOX  DPD=ONE+DPDA ;                                                    *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPD        )                          
*FOX  DPSQ=SQRT(DPD) ;                                                  *FOX
      CALL DAFUN('SQRT',DPD        ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),DPSQ       )                          
      do 220 i=1,il
        do 10 ih=1,2
          do 10 ip=1,6
*FOX  ALDA(IH,IP)=ZERO ;                                                *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(ALDA       (IH         ,IP         ),RSCRRI(100))      
*FOX  ASDA(IH,IP)=ZERO ;                                                *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(ASDA       (IH         ,IP         ),RSCRRI(100))      
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
*FOX  ALDA(L,1)=ONE  ;                                                  *FOX
      RSCRRI(100) = ONE                                                 
      CALL DACON(ALDA       (L          ,(1)),RSCRRI(100))              
*FOX  ALDA(L,2)=EL(I) ;                                                 *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      RSCRRI(100) = RSCRRI(  1+IDAA)                                    
      CALL DACON(ALDA       (L          ,(2)),RSCRRI(100))              
*FOX  ALDA(L,3)=ZERO ;                                                  *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(ALDA       (L          ,(3)),RSCRRI(100))              
*FOX  ALDA(L,4)=ONE ;                                                   *FOX
      RSCRRI(100) = ONE                                                 
      CALL DACON(ALDA       (L          ,(4)),RSCRRI(100))              
*FOX  ASDA(L,6)=-RV*ALDA(L,2)/C2E3 ;                                    *FOX
      CALL DACOP(ALDA       (L          ,(2)),ISCRDA(  1+IDAA))         
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  2+IDAA))        
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACDI(ISCRDA(  3+IDAA),ONE*C2E3       ,ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),ASDA       (L          ,(6)))         
   30   continue
*FOX  ASDA(1,1)=EL(I)*(ONE-RV)*C1E3 ;                                   *FOX
      CALL DASUC(RV         ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      RSCRRI(  2+IDAA) = EL         (I          )                       
      CALL DACMU(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  3+IDAA),ONE*C1E3       ,ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),ASDA       ((1),(1)))                 
        goto 190
!-----------------------------------------------------------------------
!  RECTANGULAR MAGNET
!  HORIZONTAL
!-----------------------------------------------------------------------
   40   ih=1
   50   continue
        if(abs(ed(i)).le.pieni) goto 20
*FOX  FOK=EL(I)*ED(I)/DPSQ ;                                            *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      RSCRRI(  2+IDAA) = ED         (I          )                       
      RSCRRI(  3+IDAA) = RSCRRI(  1+IDAA) * RSCRRI(  2+IDAA)            
      CALL DADIC(DPSQ       ,ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),FOK        )                          
*FOX  RHO=ONE/ED(I)*DPSQ ;                                              *FOX
      RSCRRI(  1+IDAA) = ED         (I          )                       
      RSCRRI(  2+IDAA) = ONE         / RSCRRI(  1+IDAA)                 
      CALL DACMU(DPSQ       ,ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),RHO        )                          
*FOX  FOK1=SIN(FOK*HALF)/COS(FOK*HALF)/RHO ;                            *FOX
      CALL DACMU(FOK        ,ONE*HALF       ,ISCRDA(  1+IDAA))          
      CALL DACMU(FOK        ,ONE*HALF       ,ISCRDA(  2+IDAA))          
      CALL DAFUN('SIN ',ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))              
      CALL DAFUN('COS ',ISCRDA(  2+IDAA),ISCRDA(  4+IDAA))              
      CALL DADIV(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DADIV(ISCRDA(  5+IDAA),RHO        ,ISCRDA(  6+IDAA))         
      CALL DACOP(ISCRDA(  6+IDAA),FOK1       )                          
*FOX  SI=SIN(FOK) ;                                                     *FOX
      CALL DAFUN('SIN ',FOK        ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),SI         )                          
*FOX  CO=COS(FOK) ;                                                     *FOX
      CALL DAFUN('COS ',FOK        ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),CO         )                          
*FOX  ALDA(IH,1)=ONE ;                                                  *FOX
      RSCRRI(100) = ONE                                                 
      CALL DACON(ALDA       (IH         ,(1)),RSCRRI(100))              
*FOX  ALDA(IH,2)=RHO*SI ;                                               *FOX
      CALL DAMUL(RHO        ,SI         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(2)))         
*FOX  ALDA(IH,3)=ZERO ;                                                 *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(ALDA       (IH         ,(3)),RSCRRI(100))              
*FOX  ALDA(IH,4)=ONE ;                                                  *FOX
      RSCRRI(100) = ONE                                                 
      CALL DACON(ALDA       (IH         ,(4)),RSCRRI(100))              
*FOX  ALDA(IH,5)=-DPDA*RHO*(ONE-CO)/DPSQ*C1E3 ;                         *FOX
      CALL DASUC(CO         ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(DPDA       ,ONE*(-ONE       ),ISCRDA(  2+IDAA))        
      CALL DAMUL(ISCRDA(  2+IDAA),RHO        ,ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  1+IDAA),ISCRDA(  4+IDAA))    
      CALL DADIV(ISCRDA(  4+IDAA),DPSQ       ,ISCRDA(  5+IDAA))         
      CALL DACMU(ISCRDA(  5+IDAA),ONE*C1E3       ,ISCRDA(  6+IDAA))     
      CALL DACOP(ISCRDA(  6+IDAA),ALDA       (IH         ,(5)))         
*FOX  ALDA(IH,6)=-DPDA*TWO*SIN(FOK*HALF)/COS(FOK*HALF)/DPSQ*C1E3 ;      *FOX
      CALL DACMU(FOK        ,ONE*HALF       ,ISCRDA(  1+IDAA))          
      CALL DACMU(FOK        ,ONE*HALF       ,ISCRDA(  2+IDAA))          
      CALL DAFUN('SIN ',ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))              
      CALL DAFUN('COS ',ISCRDA(  2+IDAA),ISCRDA(  4+IDAA))              
      CALL DACMU(DPDA       ,ONE*(-ONE       ),ISCRDA(  5+IDAA))        
      CALL DACMU(ISCRDA(  5+IDAA),ONE*TWO        ,ISCRDA(  6+IDAA))     
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  3+IDAA),ISCRDA(  7+IDAA))    
      CALL DADIV(ISCRDA(  7+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DADIV(ISCRDA(  8+IDAA),DPSQ       ,ISCRDA(  9+IDAA))         
      CALL DACMU(ISCRDA(  9+IDAA),ONE*C1E3       ,ISCRDA( 10+IDAA))     
      CALL DACOP(ISCRDA( 10+IDAA),ALDA       (IH         ,(6)))         
*FOX  SM1=COS(FOK) ;                                                    *FOX
      CALL DAFUN('COS ',FOK        ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),SM1        )                          
*FOX  SM2=SIN(FOK)*RHO ;                                                *FOX
      CALL DAFUN('SIN ',FOK        ,ISCRDA(  1+IDAA))                   
      CALL DAMUL(ISCRDA(  1+IDAA),RHO        ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),SM2        )                          
*FOX  SM3=-SIN(FOK)/RHO ;                                               *FOX
      CALL DAFUN('SIN ',FOK        ,ISCRDA(  1+IDAA))                   
      CALL DACMU(ISCRDA(  1+IDAA),ONE*(-ONE       ),ISCRDA(  2+IDAA))   
      CALL DADIV(ISCRDA(  2+IDAA),RHO        ,ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),SM3        )                          
*FOX  SM5=-RHO*DPSQ*(ONE-SM1) ;                                         *FOX
      CALL DASUC(SM1        ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(RHO        ,ONE*(-ONE       ),ISCRDA(  2+IDAA))        
      CALL DAMUL(ISCRDA(  2+IDAA),DPSQ       ,ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  1+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),SM5        )                          
*FOX  SM6=-SM2*DPSQ/RHO ;                                               *FOX
      CALL DACMU(SM2        ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),DPSQ       ,ISCRDA(  2+IDAA))         
      CALL DADIV(ISCRDA(  2+IDAA),RHO        ,ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),SM6        )                          
*FOX  SM12=EL(I)-SM1*SM2 ;                                              *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DAMUL(SM1        ,SM2        ,ISCRDA(  2+IDAA))              
      CALL DASUC(ISCRDA(  2+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  3+IDAA),SM12       )                          
*FOX  SM23=SM2*SM3 ;                                                    *FOX
      CALL DAMUL(SM2        ,SM3        ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),SM23       )                          
*FOX  AS3=-RV*(DPDA*RHO/(TWO*DPSQ)*SM23+SM5) ;                          *FOX
      CALL DACMU(DPSQ       ,ONE*TWO        ,ISCRDA(  1+IDAA))          
      CALL DAMUL(DPDA       ,RHO        ,ISCRDA(  2+IDAA))              
      CALL DADIV(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),SM23       ,ISCRDA(  4+IDAA))         
      CALL DAADD(ISCRDA(  4+IDAA),SM5        ,ISCRDA(  5+IDAA))         
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DACOP(ISCRDA(  7+IDAA),AS3        )                          
*FOX  AS4=-RV*SM23/C2E3 ;                                               *FOX
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),SM23       ,ISCRDA(  2+IDAA))         
      CALL DACDI(ISCRDA(  2+IDAA),ONE*C2E3       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),AS4        )                          
*FOX  AS6=-RV*(EL(I)+SM1*SM2)/C4E3 ;                                    *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DAMUL(SM1        ,SM2        ,ISCRDA(  2+IDAA))              
      CALL DACAD(ISCRDA(  2+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  4+IDAA))        
      CALL DAMUL(ISCRDA(  4+IDAA),ISCRDA(  3+IDAA),ISCRDA(  5+IDAA))    
      CALL DACDI(ISCRDA(  5+IDAA),ONE*C4E3       ,ISCRDA(  6+IDAA))     
      CALL DACOP(ISCRDA(  6+IDAA),AS6        )                          
*FOX  ASDA(IH,1)=(-RV*(DPDA*DPDA/(FOUR*DPD)*SM12+DPDA*(EL(I)-SM2))      *FOX
*FOX  +EL(I)*(ONE-RV))*C1E3 ;                                           *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACMU(DPD        ,ONE*FOUR       ,ISCRDA(  2+IDAA))          
      CALL DASUC(SM2        ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DAMUL(DPDA       ,DPDA       ,ISCRDA(  4+IDAA))              
      CALL DADIV(ISCRDA(  4+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),SM12       ,ISCRDA(  6+IDAA))         
      CALL DAMUL(DPDA       ,ISCRDA(  3+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DASUC(RV         ,ONE*ONE        ,ISCRDA(  9+IDAA))          
      RSCRRI( 10+IDAA) = EL         (I          )                       
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA( 11+IDAA))        
      CALL DAMUL(ISCRDA( 11+IDAA),ISCRDA(  8+IDAA),ISCRDA( 12+IDAA))    
      CALL DACMU(ISCRDA(  9+IDAA),ONE*RSCRRI( 10+IDAA),ISCRDA( 13+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA( 12+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DACMU(ISCRDA( 14+IDAA),ONE*C1E3       ,ISCRDA( 15+IDAA))     
      CALL DACOP(ISCRDA( 15+IDAA),ASDA       (IH         ,(1)))         
*FOX  ASDA(IH,2)=-RV*(DPDA/(TWO*RHO*DPSQ)*SM12+SM6)+FOK1*AS3 ;          *FOX
      CALL DACMU(RHO        ,ONE*TWO        ,ISCRDA(  1+IDAA))          
      CALL DAMUL(ISCRDA(  1+IDAA),DPSQ       ,ISCRDA(  2+IDAA))         
      CALL DADIV(DPDA       ,ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  3+IDAA),SM12       ,ISCRDA(  4+IDAA))         
      CALL DAADD(ISCRDA(  4+IDAA),SM6        ,ISCRDA(  5+IDAA))         
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(FOK1       ,AS3        ,ISCRDA(  8+IDAA))              
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),ASDA       (IH         ,(2)))         
*FOX  ASDA(IH,3)=AS3 ;                                                  *FOX
      CALL DACOP(AS3        ,ASDA       (IH         ,(3)))              
*FOX  ASDA(IH,4)=AS4+TWO*AS6*FOK1 ;                                     *FOX
      CALL DACMU(AS6        ,ONE*TWO        ,ISCRDA(  1+IDAA))          
      CALL DAMUL(ISCRDA(  1+IDAA),FOK1       ,ISCRDA(  2+IDAA))         
      CALL DAADD(AS4        ,ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),ASDA       (IH         ,(4)))         
*FOX  ASDA(IH,5)=-RV*SM12/(C4E3*RHO*RHO)+AS6*FOK1*FOK1+FOK1*AS4  ;      *FOX
      CALL DACMU(RHO        ,ONE*C4E3       ,ISCRDA(  1+IDAA))          
      CALL DAMUL(ISCRDA(  1+IDAA),RHO        ,ISCRDA(  2+IDAA))         
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  3+IDAA))        
      CALL DAMUL(ISCRDA(  3+IDAA),SM12       ,ISCRDA(  4+IDAA))         
      CALL DADIV(ISCRDA(  4+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DAMUL(AS6        ,FOK1       ,ISCRDA(  6+IDAA))              
      CALL DAMUL(ISCRDA(  6+IDAA),FOK1       ,ISCRDA(  7+IDAA))         
      CALL DAMUL(FOK1       ,AS4        ,ISCRDA(  8+IDAA))              
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DAADD(ISCRDA(  9+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),ASDA       (IH         ,(5)))         
*FOX  ASDA(IH,6)=AS6 ;                                                  *FOX
      CALL DACOP(AS6        ,ASDA       (IH         ,(6)))              
!--VERTIKAL
        ih=ih+1
        if(ih.gt.2) ih=1
*FOX  G=SIN(FOK*HALF)/COS(FOK*HALF)/RHO ;                               *FOX
      CALL DACMU(FOK        ,ONE*HALF       ,ISCRDA(  1+IDAA))          
      CALL DACMU(FOK        ,ONE*HALF       ,ISCRDA(  2+IDAA))          
      CALL DAFUN('SIN ',ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))              
      CALL DAFUN('COS ',ISCRDA(  2+IDAA),ISCRDA(  4+IDAA))              
      CALL DADIV(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DADIV(ISCRDA(  5+IDAA),RHO        ,ISCRDA(  6+IDAA))         
      CALL DACOP(ISCRDA(  6+IDAA),G          )                          
*FOX  GL=EL(I)*G ;                                                      *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACMU(G          ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),GL         )                          
*FOX  ALDA(IH,1)=ONE-GL ;                                               *FOX
      CALL DASUC(GL         ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(1)))         
*FOX  ALDA(IH,2)=EL(I) ;                                                *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      RSCRRI(100) = RSCRRI(  1+IDAA)                                    
      CALL DACON(ALDA       (IH         ,(2)),RSCRRI(100))              
*FOX  ALDA(IH,3)=-G*(TWO-GL) ;                                          *FOX
      CALL DASUC(GL         ,ONE*TWO        ,ISCRDA(  1+IDAA))          
      CALL DACMU(G          ,ONE*(-ONE       ),ISCRDA(  2+IDAA))        
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),ALDA       (IH         ,(3)))         
*FOX  ALDA(IH,4)=ALDA(IH,1) ;                                           *FOX
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  1+IDAA))         
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(4)))         
*FOX  AS6=-RV*ALDA(IH,2)/C2E3 ;                                         *FOX
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  1+IDAA))         
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  2+IDAA))        
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACDI(ISCRDA(  3+IDAA),ONE*C2E3       ,ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),AS6        )                          
*FOX  ASDA(IH,4)=-TWO*AS6*FOK1 ;                                        *FOX
      RSCRRI(  1+IDAA) = (-ONE       ) * TWO                            
      CALL DACMU(AS6        ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DAMUL(ISCRDA(  2+IDAA),FOK1       ,ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),ASDA       (IH         ,(4)))         
*FOX  ASDA(IH,5)=AS6*FOK1*FOK1 ;                                        *FOX
      CALL DAMUL(AS6        ,FOK1       ,ISCRDA(  1+IDAA))              
      CALL DAMUL(ISCRDA(  1+IDAA),FOK1       ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),ASDA       (IH         ,(5)))         
*FOX  ASDA(IH,6)=AS6 ;                                                  *FOX
      CALL DACOP(AS6        ,ASDA       (IH         ,(6)))              
        goto 190
!-----------------------------------------------------------------------
!  SEKTORMAGNET
!  HORIZONTAL
!-----------------------------------------------------------------------
   60   ih=1
   70   continue
        if(abs(ed(i)).le.pieni) goto 20
*FOX  FOK=EL(I)*ED(I)/DPSQ ;                                            *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      RSCRRI(  2+IDAA) = ED         (I          )                       
      RSCRRI(  3+IDAA) = RSCRRI(  1+IDAA) * RSCRRI(  2+IDAA)            
      CALL DADIC(DPSQ       ,ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),FOK        )                          
*FOX  RHO=(ONE/ED(I))*DPSQ ;                                            *FOX
      RSCRRI(  1+IDAA) = ED         (I          )                       
      RSCRRI(  2+IDAA) = ONE         / RSCRRI(  1+IDAA)                 
      CALL DACMU(DPSQ       ,ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),RHO        )                          
*FOX  SI=SIN(FOK) ;                                                     *FOX
      CALL DAFUN('SIN ',FOK        ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),SI         )                          
*FOX  CO=COS(FOK) ;                                                     *FOX
      CALL DAFUN('COS ',FOK        ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),CO         )                          
*FOX  RHOC=RHO*(ONE-CO)/DPSQ ;                                          *FOX
      CALL DASUC(CO         ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DAMUL(RHO        ,ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))         
      CALL DADIV(ISCRDA(  2+IDAA),DPSQ       ,ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),RHOC       )                          
*FOX  SIQ=SI/DPSQ ;                                                     *FOX
      CALL DADIV(SI         ,DPSQ       ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),SIQ        )                          
*FOX  ALDA(IH,1)=CO ;                                                   *FOX
      CALL DACOP(CO         ,ALDA       (IH         ,(1)))              
*FOX  ALDA(IH,2)=RHO*SI ;                                               *FOX
      CALL DAMUL(RHO        ,SI         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(2)))         
*FOX  ALDA(IH,3)=-SI/RHO ;                                              *FOX
      CALL DACMU(SI         ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DADIV(ISCRDA(  1+IDAA),RHO        ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),ALDA       (IH         ,(3)))         
*FOX  ALDA(IH,4)=CO ;                                                   *FOX
      CALL DACOP(CO         ,ALDA       (IH         ,(4)))              
*FOX  ALDA(IH,5)=-DPDA*RHOC*C1E3 ;                                      *FOX
      CALL DACMU(DPDA       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),RHOC       ,ISCRDA(  2+IDAA))         
      CALL DACMU(ISCRDA(  2+IDAA),ONE*C1E3       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),ALDA       (IH         ,(5)))         
*FOX  ALDA(IH,6)=-DPDA*SIQ*C1E3 ;                                       *FOX
      CALL DACMU(DPDA       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),SIQ        ,ISCRDA(  2+IDAA))         
      CALL DACMU(ISCRDA(  2+IDAA),ONE*C1E3       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),ALDA       (IH         ,(6)))         
*FOX  SM12=EL(I)-ALDA(IH,1)*ALDA(IH,2) ;                                *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DASUC(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  5+IDAA),SM12       )                          
*FOX  SM23=ALDA(IH,2)*ALDA(IH,3) ;                                      *FOX
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  1+IDAA))         
      CALL DACOP(ALDA       (IH         ,(3)),ISCRDA(  2+IDAA))         
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),SM23       )                          
*FOX  ASDA(IH,1)=(-RV*(DPDA*DPDA/(FOUR*DPD)*SM12                        *FOX
*FOX  +DPDA*(EL(I)-ALDA(IH,2)))+EL(I)*(ONE-RV))*C1E3 ;                  *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  2+IDAA))         
      CALL DACMU(DPD        ,ONE*FOUR       ,ISCRDA(  3+IDAA))          
      CALL DASUC(ISCRDA(  2+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      CALL DAMUL(DPDA       ,DPDA       ,ISCRDA(  5+IDAA))              
      CALL DADIV(ISCRDA(  5+IDAA),ISCRDA(  3+IDAA),ISCRDA(  6+IDAA))    
      CALL DAMUL(ISCRDA(  6+IDAA),SM12       ,ISCRDA(  7+IDAA))         
      CALL DAMUL(DPDA       ,ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))         
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUC(RV         ,ONE*ONE        ,ISCRDA( 10+IDAA))          
      RSCRRI( 11+IDAA) = EL         (I          )                       
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA( 12+IDAA))        
      CALL DAMUL(ISCRDA( 12+IDAA),ISCRDA(  9+IDAA),ISCRDA( 13+IDAA))    
      CALL DACMU(ISCRDA( 10+IDAA),ONE*RSCRRI( 11+IDAA),ISCRDA( 14+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA( 13+IDAA),ISCRDA( 14+IDAA),ISCRDA( 15+IDAA))    
      CALL DACMU(ISCRDA( 15+IDAA),ONE*C1E3       ,ISCRDA( 16+IDAA))     
      CALL DACOP(ISCRDA( 16+IDAA),ASDA       (IH         ,(1)))         
*FOX  ASDA(IH,2)=-RV*(DPDA/(TWO*RHO*DPSQ)*SM12-DPD*SIQ) ;               *FOX
      CALL DACMU(RHO        ,ONE*TWO        ,ISCRDA(  1+IDAA))          
      CALL DAMUL(ISCRDA(  1+IDAA),DPSQ       ,ISCRDA(  2+IDAA))         
      CALL DADIV(DPDA       ,ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  3+IDAA),SM12       ,ISCRDA(  4+IDAA))         
      CALL DAMUL(DPD        ,SIQ        ,ISCRDA(  5+IDAA))              
      CALL DASUB(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  7+IDAA))        
      CALL DAMUL(ISCRDA(  7+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),ASDA       (IH         ,(2)))         
*FOX  ASDA(IH,3)=-RV*(DPDA*RHO/(TWO*DPSQ)*SM23-DPD*RHOC) ;              *FOX
      CALL DACMU(DPSQ       ,ONE*TWO        ,ISCRDA(  1+IDAA))          
      CALL DAMUL(DPDA       ,RHO        ,ISCRDA(  2+IDAA))              
      CALL DADIV(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),SM23       ,ISCRDA(  4+IDAA))         
      CALL DAMUL(DPD        ,RHOC       ,ISCRDA(  5+IDAA))              
      CALL DASUB(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  7+IDAA))        
      CALL DAMUL(ISCRDA(  7+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),ASDA       (IH         ,(3)))         
*FOX  ASDA(IH,4)=-RV*SM23/C2E3 ;                                        *FOX
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),SM23       ,ISCRDA(  2+IDAA))         
      CALL DACDI(ISCRDA(  2+IDAA),ONE*C2E3       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),ASDA       (IH         ,(4)))         
*FOX  ASDA(IH,5)=-RV*SM12/(C4E3*RHO*RHO) ;                              *FOX
      CALL DACMU(RHO        ,ONE*C4E3       ,ISCRDA(  1+IDAA))          
      CALL DAMUL(ISCRDA(  1+IDAA),RHO        ,ISCRDA(  2+IDAA))         
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  3+IDAA))        
      CALL DAMUL(ISCRDA(  3+IDAA),SM12       ,ISCRDA(  4+IDAA))         
      CALL DADIV(ISCRDA(  4+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(ISCRDA(  5+IDAA),ASDA       (IH         ,(5)))         
*FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;               *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACAD(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DACDI(ISCRDA(  7+IDAA),ONE*C4E3       ,ISCRDA(  8+IDAA))     
      CALL DACOP(ISCRDA(  8+IDAA),ASDA       (IH         ,(6)))         
!--VERTIKAL
        ih=ih+1
        if(ih.gt.2) ih=1
*FOX  ALDA(IH,1)=ONE ;                                                  *FOX
      RSCRRI(100) = ONE                                                 
      CALL DACON(ALDA       (IH         ,(1)),RSCRRI(100))              
*FOX  ALDA(IH,2)=EL(I) ;                                                *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      RSCRRI(100) = RSCRRI(  1+IDAA)                                    
      CALL DACON(ALDA       (IH         ,(2)),RSCRRI(100))              
*FOX  ALDA(IH,3)=ZERO ;                                                 *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(ALDA       (IH         ,(3)),RSCRRI(100))              
*FOX  ALDA(IH,4)=ONE ;                                                  *FOX
      RSCRRI(100) = ONE                                                 
      CALL DACON(ALDA       (IH         ,(4)),RSCRRI(100))              
*FOX  ASDA(IH,6)=-RV*ALDA(IH,2)/C2E3 ;                                  *FOX
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  1+IDAA))         
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  2+IDAA))        
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACDI(ISCRDA(  3+IDAA),ONE*C2E3       ,ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),ASDA       (IH         ,(6)))         
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
*FOX  FOK=EK(I)/(ONE+DPDA) ;                                            *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      RSCRRI(  2+IDAA) = EK         (I          )                       
      CALL DADIC(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  3+IDAA),FOK        )                          
*FOX  AEK=FOK ;                                                         *FOX
      CALL DACOP(FOK        ,AEK        )                               
        if(dare(aek).lt.zero) then
*FOX  AEK=-AEK ;                                                        *FOX
      CALL DACMU(AEK        ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),AEK        )                          
        endif
        ih=0
*FOX  HI=SQRT(AEK) ;                                                    *FOX
      CALL DAFUN('SQRT',AEK        ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),HI         )                          
*FOX  FI=EL(I)*HI ;                                                     *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACMU(HI         ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),FI         )                          
        if(ek(i).gt.zero) goto 120
  110   ih=ih+1
*FOX  ALDA(IH,1)=COS(FI) ;                                              *FOX
      CALL DAFUN('COS ',FI         ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(1)))         
*FOX  HI1=SIN(FI) ;                                                     *FOX
      CALL DAFUN('SIN ',FI         ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),HI1        )                          
*FOX  ALDA(IH,2)=HI1/HI ;                                               *FOX
      CALL DADIV(HI1        ,HI         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(2)))         
*FOX  ALDA(IH,3)=-HI1*HI ;                                              *FOX
      CALL DACMU(HI1        ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),HI         ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),ALDA       (IH         ,(3)))         
*FOX  ALDA(IH,4)=ALDA(IH,1) ;                                           *FOX
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  1+IDAA))         
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(4)))         
*FOX  ASDA(IH,1)=EL(I)*(ONE-RV)*C1E3 ;                                  *FOX
      CALL DASUC(RV         ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      RSCRRI(  2+IDAA) = EL         (I          )                       
      CALL DACMU(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  3+IDAA),ONE*C1E3       ,ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),ASDA       (IH         ,(1)))         
*FOX  ASDA(IH,4)=-RV*ALDA(IH,2)*ALDA(IH,3)/C2E3 ;                       *FOX
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  1+IDAA))         
      CALL DACOP(ALDA       (IH         ,(3)),ISCRDA(  2+IDAA))         
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  3+IDAA))        
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  1+IDAA),ISCRDA(  4+IDAA))    
      CALL DAMUL(ISCRDA(  4+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DACDI(ISCRDA(  5+IDAA),ONE*C2E3       ,ISCRDA(  6+IDAA))     
      CALL DACOP(ISCRDA(  6+IDAA),ASDA       (IH         ,(4)))         
*FOX  ASDA(IH,5)=-RV*(EL(I)-ALDA(IH,1)*ALDA(IH,2))*AEK/C4E3 ;           *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DASUC(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  7+IDAA),AEK        ,ISCRDA(  8+IDAA))         
      CALL DACDI(ISCRDA(  8+IDAA),ONE*C4E3       ,ISCRDA(  9+IDAA))     
      CALL DACOP(ISCRDA(  9+IDAA),ASDA       (IH         ,(5)))         
*FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;               *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACAD(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DACDI(ISCRDA(  7+IDAA),ONE*C4E3       ,ISCRDA(  8+IDAA))     
      CALL DACOP(ISCRDA(  8+IDAA),ASDA       (IH         ,(6)))         
        if(ih.eq.2) goto 190
!--DEFOCUSSING
  120   ih=ih+1
*FOX  HP=EXP(FI) ;                                                      *FOX
      CALL DAFUN('EXP ',FI         ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),HP         )                          
*FOX  HM=ONE/HP ;                                                       *FOX
      CALL DADIC(HP         ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),HM         )                          
*FOX  HC=(HP+HM)*HALF ;                                                 *FOX
      CALL DAADD(HP         ,HM         ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*HALF       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),HC         )                          
*FOX  HS=(HP-HM)*HALF ;                                                 *FOX
      CALL DASUB(HP         ,HM         ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*HALF       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),HS         )                          
*FOX  ALDA(IH,1)=HC ;                                                   *FOX
      CALL DACOP(HC         ,ALDA       (IH         ,(1)))              
*FOX  ALDA(IH,2)=HS/HI ;                                                *FOX
      CALL DADIV(HS         ,HI         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(2)))         
*FOX  ALDA(IH,3)=HS*HI ;                                                *FOX
      CALL DAMUL(HS         ,HI         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(3)))         
*FOX  ALDA(IH,4)=HC ;                                                   *FOX
      CALL DACOP(HC         ,ALDA       (IH         ,(4)))              
*FOX  ASDA(IH,4)=-RV*ALDA(IH,2)*ALDA(IH,3)/C2E3 ;                       *FOX
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  1+IDAA))         
      CALL DACOP(ALDA       (IH         ,(3)),ISCRDA(  2+IDAA))         
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  3+IDAA))        
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  1+IDAA),ISCRDA(  4+IDAA))    
      CALL DAMUL(ISCRDA(  4+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DACDI(ISCRDA(  5+IDAA),ONE*C2E3       ,ISCRDA(  6+IDAA))     
      CALL DACOP(ISCRDA(  6+IDAA),ASDA       (IH         ,(4)))         
*FOX  ASDA(IH,5)=+RV*(EL(I)-ALDA(IH,1)*ALDA(IH,2))*AEK/C4E3 ;           *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DASUC(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DAMUL(RV         ,ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))         
      CALL DAMUL(ISCRDA(  6+IDAA),AEK        ,ISCRDA(  7+IDAA))         
      CALL DACDI(ISCRDA(  7+IDAA),ONE*C4E3       ,ISCRDA(  8+IDAA))     
      CALL DACOP(ISCRDA(  8+IDAA),ASDA       (IH         ,(5)))         
*FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;               *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACAD(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DACDI(ISCRDA(  7+IDAA),ONE*C4E3       ,ISCRDA(  8+IDAA))     
      CALL DACOP(ISCRDA(  8+IDAA),ASDA       (IH         ,(6)))         
        if(ih.eq.1) goto 110
        goto 190
!-----------------------------------------------------------------------
!  COMBINED FUNCTION MAGNET HORIZONTAL
!  FOCUSSING
!-----------------------------------------------------------------------
  130   ih=0
*FOX  FOKQ=EK(I) ;                                                      *FOX
      RSCRRI(  1+IDAA) = EK         (I          )                       
      RSCRRI(100) = RSCRRI(  1+IDAA)                                    
      CALL DACON(FOKQ       ,RSCRRI(100))                               
  140   continue
        if(abs(ek(i)).le.pieni) goto 60
        if(abs(ed(i)).le.pieni) goto 100
!hr08   if(abs(ek(i)-ed(i)*ed(i)).le.pieni) goto 20
        if(abs(ek(i)-ed(i)**2).le.pieni) goto 20                         !hr08
*FOX  WF=ED(I)/DPSQ ;                                                   *FOX
      RSCRRI(  1+IDAA) = ED         (I          )                       
      CALL DADIC(DPSQ       ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),WF         )                          
*FOX  FOK=FOKQ/DPD-WF*WF ;                                              *FOX
      CALL DADIV(FOKQ       ,DPD        ,ISCRDA(  1+IDAA))              
      CALL DAMUL(WF         ,WF         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),FOK        )                          
*FOX  AFOK=FOK ;                                                        *FOX
      CALL DACOP(FOK        ,AFOK       )                               
      if(dare(afok).lt.zero) then
*FOX  AFOK=-AFOK ;                                                      *FOX
      CALL DACMU(AFOK       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),AFOK       )                          
      endif
*FOX  HI=SQRT(AFOK) ;                                                   *FOX
      CALL DAFUN('SQRT',AFOK       ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),HI         )                          
*FOX  FI=HI*EL(I) ;                                                     *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACMU(HI         ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),FI         )                          
        if(dare(fok).gt.zero) goto 160
  150   ih=ih+1
*FOX  SI=SIN(FI) ;                                                      *FOX
      CALL DAFUN('SIN ',FI         ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),SI         )                          
*FOX  CO=COS(FI) ;                                                      *FOX
      CALL DAFUN('COS ',FI         ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),CO         )                          
*FOX  WFA=WF/AFOK*(ONE-CO)/DPSQ ;                                       *FOX
      CALL DASUC(CO         ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DADIV(WF         ,AFOK       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DADIV(ISCRDA(  3+IDAA),DPSQ       ,ISCRDA(  4+IDAA))         
      CALL DACOP(ISCRDA(  4+IDAA),WFA        )                          
*FOX  WFHI=WF/HI*SI/DPSQ ;                                              *FOX
      CALL DADIV(WF         ,HI         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(ISCRDA(  1+IDAA),SI         ,ISCRDA(  2+IDAA))         
      CALL DADIV(ISCRDA(  2+IDAA),DPSQ       ,ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),WFHI       )                          
*FOX  ALDA(IH,1)=CO ;                                                   *FOX
      CALL DACOP(CO         ,ALDA       (IH         ,(1)))              
*FOX  ALDA(IH,2)=SI/HI ;                                                *FOX
      CALL DADIV(SI         ,HI         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(2)))         
*FOX  ALDA(IH,3)=-SI*HI ;                                               *FOX
      CALL DACMU(SI         ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),HI         ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),ALDA       (IH         ,(3)))         
*FOX  ALDA(IH,4)=CO ;                                                   *FOX
      CALL DACOP(CO         ,ALDA       (IH         ,(4)))              
*FOX  ALDA(IH,5)=-WFA*DPDA*C1E3 ;                                       *FOX
      CALL DACMU(WFA        ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),DPDA       ,ISCRDA(  2+IDAA))         
      CALL DACMU(ISCRDA(  2+IDAA),ONE*C1E3       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),ALDA       (IH         ,(5)))         
*FOX  ALDA(IH,6)=-WFHI*DPDA*C1E3 ;                                      *FOX
      CALL DACMU(WFHI       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),DPDA       ,ISCRDA(  2+IDAA))         
      CALL DACMU(ISCRDA(  2+IDAA),ONE*C1E3       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),ALDA       (IH         ,(6)))         
*FOX  SM12=EL(I)-ALDA(IH,1)*ALDA(IH,2) ;                                *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DASUC(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  5+IDAA),SM12       )                          
*FOX  SM23=ALDA(IH,2)*ALDA(IH,3) ;                                      *FOX
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  1+IDAA))         
      CALL DACOP(ALDA       (IH         ,(3)),ISCRDA(  2+IDAA))         
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),SM23       )                          
*FOX  ASDA(IH,1)=(-RV*(DPDA*DPDA/(FOUR*DPD)*SM12+DPDA                   *FOX
*FOX  *(EL(I)-ALDA(IH,2)))/AFOK*WF*WF+EL(I)*(ONE-RV))*C1E3 ;            *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  2+IDAA))         
      CALL DACMU(DPD        ,ONE*FOUR       ,ISCRDA(  3+IDAA))          
      CALL DASUC(ISCRDA(  2+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      CALL DAMUL(DPDA       ,DPDA       ,ISCRDA(  5+IDAA))              
      CALL DADIV(ISCRDA(  5+IDAA),ISCRDA(  3+IDAA),ISCRDA(  6+IDAA))    
      CALL DAMUL(ISCRDA(  6+IDAA),SM12       ,ISCRDA(  7+IDAA))         
      CALL DAMUL(DPDA       ,ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))         
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUC(RV         ,ONE*ONE        ,ISCRDA( 10+IDAA))          
      RSCRRI( 11+IDAA) = EL         (I          )                       
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA( 12+IDAA))        
      CALL DAMUL(ISCRDA( 12+IDAA),ISCRDA(  9+IDAA),ISCRDA( 13+IDAA))    
      CALL DADIV(ISCRDA( 13+IDAA),AFOK       ,ISCRDA( 14+IDAA))         
      CALL DAMUL(ISCRDA( 14+IDAA),WF         ,ISCRDA( 15+IDAA))         
      CALL DAMUL(ISCRDA( 15+IDAA),WF         ,ISCRDA( 16+IDAA))         
      CALL DACMU(ISCRDA( 10+IDAA),ONE*RSCRRI( 11+IDAA),ISCRDA( 17+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA( 16+IDAA),ISCRDA( 17+IDAA),ISCRDA( 18+IDAA))    
      CALL DACMU(ISCRDA( 18+IDAA),ONE*C1E3       ,ISCRDA( 19+IDAA))     
      CALL DACOP(ISCRDA( 19+IDAA),ASDA       (IH         ,(1)))         
*FOX  ASDA(IH,2)=-RV*(DPDA*WF/(TWO*DPSQ)*SM12-DPD*WFHI) ;               *FOX
      CALL DACMU(DPSQ       ,ONE*TWO        ,ISCRDA(  1+IDAA))          
      CALL DAMUL(DPDA       ,WF         ,ISCRDA(  2+IDAA))              
      CALL DADIV(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),SM12       ,ISCRDA(  4+IDAA))         
      CALL DAMUL(DPD        ,WFHI       ,ISCRDA(  5+IDAA))              
      CALL DASUB(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  7+IDAA))        
      CALL DAMUL(ISCRDA(  7+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),ASDA       (IH         ,(2)))         
*FOX  ASDA(IH,3)=-RV*(DPDA*HALF/AFOK/DPD*ED(I)*SM23-DPD*WFA) ;          *FOX
      RSCRRI(  1+IDAA) = ED         (I          )                       
      CALL DACMU(DPDA       ,ONE*HALF       ,ISCRDA(  2+IDAA))          
      CALL DADIV(ISCRDA(  2+IDAA),AFOK       ,ISCRDA(  3+IDAA))         
      CALL DADIV(ISCRDA(  3+IDAA),DPD        ,ISCRDA(  4+IDAA))         
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA(  5+IDAA),SM23       ,ISCRDA(  6+IDAA))         
      CALL DAMUL(DPD        ,WFA        ,ISCRDA(  7+IDAA))              
      CALL DASUB(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  9+IDAA))        
      CALL DAMUL(ISCRDA(  9+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),ASDA       (IH         ,(3)))         
*FOX  ASDA(IH,4)=-RV*SM23/C2E3 ;                                        *FOX
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),SM23       ,ISCRDA(  2+IDAA))         
      CALL DACDI(ISCRDA(  2+IDAA),ONE*C2E3       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),ASDA       (IH         ,(4)))         
*FOX  ASDA(IH,5)=-RV*SM12*AFOK/C4E3 ;                                   *FOX
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),SM12       ,ISCRDA(  2+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),AFOK       ,ISCRDA(  3+IDAA))         
      CALL DACDI(ISCRDA(  3+IDAA),ONE*C4E3       ,ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),ASDA       (IH         ,(5)))         
*FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;               *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACAD(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DACDI(ISCRDA(  7+IDAA),ONE*C4E3       ,ISCRDA(  8+IDAA))     
      CALL DACOP(ISCRDA(  8+IDAA),ASDA       (IH         ,(6)))         
        ih=ih+1
        if(ih.gt.2) ih=1
*FOX  AEK=EK(I)/DPD ;                                                   *FOX
      RSCRRI(  1+IDAA) = EK         (I          )                       
      CALL DADIC(DPD        ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),AEK        )                          
      if(dare(aek).lt.zero) then
*FOX  AEK=-AEK ;                                                        *FOX
      CALL DACMU(AEK        ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),AEK        )                          
      endif
*FOX  HI=SQRT(AEK) ;                                                    *FOX
      CALL DAFUN('SQRT',AEK        ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),HI         )                          
*FOX  FI=HI*EL(I) ;                                                     *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACMU(HI         ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),FI         )                          
*FOX  HP=EXP(FI) ;                                                      *FOX
      CALL DAFUN('EXP ',FI         ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),HP         )                          
*FOX  HM=ONE/HP ;                                                       *FOX
      CALL DADIC(HP         ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),HM         )                          
*FOX  HC=(HP+HM)*HALF ;                                                 *FOX
      CALL DAADD(HP         ,HM         ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*HALF       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),HC         )                          
*FOX  HS=(HP-HM)*HALF ;                                                 *FOX
      CALL DASUB(HP         ,HM         ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*HALF       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),HS         )                          
*FOX  ALDA(IH,1)=HC ;                                                   *FOX
      CALL DACOP(HC         ,ALDA       (IH         ,(1)))              
*FOX  ALDA(IH,2)=EL(I) ;                                                *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      RSCRRI(100) = RSCRRI(  1+IDAA)                                    
      CALL DACON(ALDA       (IH         ,(2)),RSCRRI(100))              
*FOX  ALDA(IH,2)=HS/HI ;                                                *FOX
      CALL DADIV(HS         ,HI         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(2)))         
*FOX  ALDA(IH,3)=HS*HI ;                                                *FOX
      CALL DAMUL(HS         ,HI         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(3)))         
*FOX  ALDA(IH,4)=HC ;                                                   *FOX
      CALL DACOP(HC         ,ALDA       (IH         ,(4)))              
*FOX  ASDA(IH,4)=-RV*ALDA(IH,2)*ALDA(IH,3)/C2E3 ;                       *FOX
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  1+IDAA))         
      CALL DACOP(ALDA       (IH         ,(3)),ISCRDA(  2+IDAA))         
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  3+IDAA))        
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  1+IDAA),ISCRDA(  4+IDAA))    
      CALL DAMUL(ISCRDA(  4+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DACDI(ISCRDA(  5+IDAA),ONE*C2E3       ,ISCRDA(  6+IDAA))     
      CALL DACOP(ISCRDA(  6+IDAA),ASDA       (IH         ,(4)))         
*FOX  ASDA(IH,5)=+RV*(EL(I)-ALDA(IH,1)*ALDA(IH,2))*AEK/C4E3 ;           *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DASUC(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DAMUL(RV         ,ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))         
      CALL DAMUL(ISCRDA(  6+IDAA),AEK        ,ISCRDA(  7+IDAA))         
      CALL DACDI(ISCRDA(  7+IDAA),ONE*C4E3       ,ISCRDA(  8+IDAA))     
      CALL DACOP(ISCRDA(  8+IDAA),ASDA       (IH         ,(5)))         
*FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;               *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACAD(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DACDI(ISCRDA(  7+IDAA),ONE*C4E3       ,ISCRDA(  8+IDAA))     
      CALL DACOP(ISCRDA(  8+IDAA),ASDA       (IH         ,(6)))         
        goto 190
!--DEFOCUSSING
  160   ih=ih+1
*FOX  HP=EXP(FI) ;                                                      *FOX
      CALL DAFUN('EXP ',FI         ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),HP         )                          
*FOX  HM=ONE/HP ;                                                       *FOX
      CALL DADIC(HP         ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),HM         )                          
*FOX  HC=(HP+HM)*HALF ;                                                 *FOX
      CALL DAADD(HP         ,HM         ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*HALF       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),HC         )                          
*FOX  HS=(HP-HM)*HALF ;                                                 *FOX
      CALL DASUB(HP         ,HM         ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*HALF       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),HS         )                          
*FOX  ALDA(IH,1)=HC ;                                                   *FOX
      CALL DACOP(HC         ,ALDA       (IH         ,(1)))              
*FOX  ALDA(IH,2)=HS/HI ;                                                *FOX
      CALL DADIV(HS         ,HI         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(2)))         
*FOX  ALDA(IH,3)=HS*HI ;                                                *FOX
      CALL DAMUL(HS         ,HI         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(3)))         
*FOX  ALDA(IH,4)=HC ;                                                   *FOX
      CALL DACOP(HC         ,ALDA       (IH         ,(4)))              
*FOX  WFA=WF/AFOK*(ONE-HC)/DPSQ ;                                       *FOX
      CALL DASUC(HC         ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DADIV(WF         ,AFOK       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DADIV(ISCRDA(  3+IDAA),DPSQ       ,ISCRDA(  4+IDAA))         
      CALL DACOP(ISCRDA(  4+IDAA),WFA        )                          
*FOX  WFHI=WF/HI*HS/DPSQ ;                                              *FOX
      CALL DADIV(WF         ,HI         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(ISCRDA(  1+IDAA),HS         ,ISCRDA(  2+IDAA))         
      CALL DADIV(ISCRDA(  2+IDAA),DPSQ       ,ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),WFHI       )                          
*FOX  ALDA(IH,5)= WFA*DPDA*C1E3 ;                                       *FOX
      CALL DAMUL(WFA        ,DPDA       ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1E3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),ALDA       (IH         ,(5)))         
*FOX  ALDA(IH,6)=-WFHI*DPDA*C1E3 ;                                      *FOX
      CALL DACMU(WFHI       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),DPDA       ,ISCRDA(  2+IDAA))         
      CALL DACMU(ISCRDA(  2+IDAA),ONE*C1E3       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),ALDA       (IH         ,(6)))         
*FOX  SM12=EL(I)-ALDA(IH,1)*ALDA(IH,2) ;                                *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DASUC(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  5+IDAA),SM12       )                          
*FOX  SM23=ALDA(IH,2)*ALDA(IH,3) ;                                      *FOX
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  1+IDAA))         
      CALL DACOP(ALDA       (IH         ,(3)),ISCRDA(  2+IDAA))         
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),SM23       )                          
*FOX  ASDA(IH,1)=(RV*(DPDA*DPDA/(FOUR*DPD)*SM12                         *FOX
*FOX  +DPDA*(EL(I)-ALDA(IH,2)))/AFOK*WF*WF+EL(I)*(ONE-RV))*C1E3 ;       *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  2+IDAA))         
      CALL DACMU(DPD        ,ONE*FOUR       ,ISCRDA(  3+IDAA))          
      CALL DASUC(ISCRDA(  2+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      CALL DAMUL(DPDA       ,DPDA       ,ISCRDA(  5+IDAA))              
      CALL DADIV(ISCRDA(  5+IDAA),ISCRDA(  3+IDAA),ISCRDA(  6+IDAA))    
      CALL DAMUL(ISCRDA(  6+IDAA),SM12       ,ISCRDA(  7+IDAA))         
      CALL DAMUL(DPDA       ,ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))         
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUC(RV         ,ONE*ONE        ,ISCRDA( 10+IDAA))          
      RSCRRI( 11+IDAA) = EL         (I          )                       
      CALL DAMUL(RV         ,ISCRDA(  9+IDAA),ISCRDA( 12+IDAA))         
      CALL DADIV(ISCRDA( 12+IDAA),AFOK       ,ISCRDA( 13+IDAA))         
      CALL DAMUL(ISCRDA( 13+IDAA),WF         ,ISCRDA( 14+IDAA))         
      CALL DAMUL(ISCRDA( 14+IDAA),WF         ,ISCRDA( 15+IDAA))         
      CALL DACMU(ISCRDA( 10+IDAA),ONE*RSCRRI( 11+IDAA),ISCRDA( 16+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA( 15+IDAA),ISCRDA( 16+IDAA),ISCRDA( 17+IDAA))    
      CALL DACMU(ISCRDA( 17+IDAA),ONE*C1E3       ,ISCRDA( 18+IDAA))     
      CALL DACOP(ISCRDA( 18+IDAA),ASDA       (IH         ,(1)))         
*FOX  ASDA(IH,2)=-RV*(DPDA*WF/(TWO*DPSQ)*SM12-DPD*WFHI) ;               *FOX
      CALL DACMU(DPSQ       ,ONE*TWO        ,ISCRDA(  1+IDAA))          
      CALL DAMUL(DPDA       ,WF         ,ISCRDA(  2+IDAA))              
      CALL DADIV(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),SM12       ,ISCRDA(  4+IDAA))         
      CALL DAMUL(DPD        ,WFHI       ,ISCRDA(  5+IDAA))              
      CALL DASUB(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  7+IDAA))        
      CALL DAMUL(ISCRDA(  7+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),ASDA       (IH         ,(2)))         
*FOX  ASDA(IH,3)=RV*(DPDA*HALF/AFOK/DPD*ED(I)*SM23-DPD*WFA) ;           *FOX
      RSCRRI(  1+IDAA) = ED         (I          )                       
      CALL DACMU(DPDA       ,ONE*HALF       ,ISCRDA(  2+IDAA))          
      CALL DADIV(ISCRDA(  2+IDAA),AFOK       ,ISCRDA(  3+IDAA))         
      CALL DADIV(ISCRDA(  3+IDAA),DPD        ,ISCRDA(  4+IDAA))         
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA(  5+IDAA),SM23       ,ISCRDA(  6+IDAA))         
      CALL DAMUL(DPD        ,WFA        ,ISCRDA(  7+IDAA))              
      CALL DASUB(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(RV         ,ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))         
      CALL DACOP(ISCRDA(  9+IDAA),ASDA       (IH         ,(3)))         
*FOX  ASDA(IH,4)=-RV*SM23/C2E3 ;                                        *FOX
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),SM23       ,ISCRDA(  2+IDAA))         
      CALL DACDI(ISCRDA(  2+IDAA),ONE*C2E3       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),ASDA       (IH         ,(4)))         
*FOX  ASDA(IH,5)=+RV*SM12*AFOK/C4E3 ;                                   *FOX
      CALL DAMUL(RV         ,SM12       ,ISCRDA(  1+IDAA))              
      CALL DAMUL(ISCRDA(  1+IDAA),AFOK       ,ISCRDA(  2+IDAA))         
      CALL DACDI(ISCRDA(  2+IDAA),ONE*C4E3       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),ASDA       (IH         ,(5)))         
*FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;               *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACAD(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DACDI(ISCRDA(  7+IDAA),ONE*C4E3       ,ISCRDA(  8+IDAA))     
      CALL DACOP(ISCRDA(  8+IDAA),ASDA       (IH         ,(6)))         
        ih=ih+1
        if(ih.gt.2) ih=1
*FOX  AEK=EK(I)/DPD ;                                                   *FOX
      RSCRRI(  1+IDAA) = EK         (I          )                       
      CALL DADIC(DPD        ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),AEK        )                          
      if(dare(aek).lt.zero) then
*FOX  AEK=-AEK ;                                                        *FOX
      CALL DACMU(AEK        ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),AEK        )                          
      endif
*FOX  HI=SQRT(AEK) ;                                                    *FOX
      CALL DAFUN('SQRT',AEK        ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),HI         )                          
*FOX  FI=HI*EL(I) ;                                                     *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACMU(HI         ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),FI         )                          
*FOX  SI=SIN(FI) ;                                                      *FOX
      CALL DAFUN('SIN ',FI         ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),SI         )                          
*FOX  CO=COS(FI) ;                                                      *FOX
      CALL DAFUN('COS ',FI         ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),CO         )                          
*FOX  ALDA(IH,1)=CO ;                                                   *FOX
      CALL DACOP(CO         ,ALDA       (IH         ,(1)))              
*FOX  ALDA(IH,2)=SI/HI ;                                                *FOX
      CALL DADIV(SI         ,HI         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       (IH         ,(2)))         
*FOX  ALDA(IH,3)=-SI*HI ;                                               *FOX
      CALL DACMU(SI         ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),HI         ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),ALDA       (IH         ,(3)))         
*FOX  ALDA(IH,4)=CO ;                                                   *FOX
      CALL DACOP(CO         ,ALDA       (IH         ,(4)))              
*FOX  ASDA(IH,4)=-RV*ALDA(IH,2)*ALDA(IH,3)/C2E3 ;                       *FOX
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  1+IDAA))         
      CALL DACOP(ALDA       (IH         ,(3)),ISCRDA(  2+IDAA))         
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  3+IDAA))        
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  1+IDAA),ISCRDA(  4+IDAA))    
      CALL DAMUL(ISCRDA(  4+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DACDI(ISCRDA(  5+IDAA),ONE*C2E3       ,ISCRDA(  6+IDAA))     
      CALL DACOP(ISCRDA(  6+IDAA),ASDA       (IH         ,(4)))         
*FOX  ASDA(IH,5)=-RV*(EL(I)-ALDA(IH,1)*ALDA(IH,2))*AEK/C4E3 ;           *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DASUC(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  7+IDAA),AEK        ,ISCRDA(  8+IDAA))         
      CALL DACDI(ISCRDA(  8+IDAA),ONE*C4E3       ,ISCRDA(  9+IDAA))     
      CALL DACOP(ISCRDA(  9+IDAA),ASDA       (IH         ,(5)))         
*FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;               *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDA       (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDA       (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACAD(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DACDI(ISCRDA(  7+IDAA),ONE*C4E3       ,ISCRDA(  8+IDAA))     
      CALL DACOP(ISCRDA(  8+IDAA),ASDA       (IH         ,(6)))         
        goto 190
!-----------------------------------------------------------------------
!  COMBINED FUNCTION MAGNET VERTICAL
!-----------------------------------------------------------------------
  170   ih=1
*FOX  FOKQ=-EK(I) ;                                                     *FOX
      RSCRRI(  1+IDAA) = EK         (I          )                       
      RSCRRI(  2+IDAA) = (-ONE       ) * RSCRRI(  1+IDAA)               
      RSCRRI(100) = RSCRRI(  2+IDAA)                                    
      CALL DACON(FOKQ       ,RSCRRI(100))                               
        goto 140
!-----------------------------------------------------------------------
!  EDGE FOCUSSING
!-----------------------------------------------------------------------
  180   continue
*FOX  RHOI=ED(I)/DPSQ ;                                                 *FOX
      RSCRRI(  1+IDAA) = ED         (I          )                       
      CALL DADIC(DPSQ       ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),RHOI       )                          
*FOX  FOK=RHOI*SIN(EL(I)*RHOI*HALF)/COS(EL(I)*RHOI*HALF) ;              *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      RSCRRI(  2+IDAA) = EL         (I          )                       
      CALL DACMU(RHOI       ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*HALF       ,ISCRDA(  4+IDAA))     
      CALL DACMU(RHOI       ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DACMU(ISCRDA(  5+IDAA),ONE*HALF       ,ISCRDA(  6+IDAA))     
      CALL DAFUN('SIN ',ISCRDA(  4+IDAA),ISCRDA(  7+IDAA))              
      CALL DAFUN('COS ',ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))              
      CALL DAMUL(RHOI       ,ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))         
      CALL DADIV(ISCRDA(  9+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),FOK        )                          
*FOX  ALDA(1,1)=ONE ;                                                   *FOX
      RSCRRI(100) = ONE                                                 
      CALL DACON(ALDA       ((1),(1)),RSCRRI(100))                      
*FOX  ALDA(1,2)=ZERO ;                                                  *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(ALDA       ((1),(2)),RSCRRI(100))                      
*FOX  ALDA(1,3)=FOK ;                                                   *FOX
      CALL DACOP(FOK        ,ALDA       ((1),(3)))                      
*FOX  ALDA(1,4)=ONE ;                                                   *FOX
      RSCRRI(100) = ONE                                                 
      CALL DACON(ALDA       ((1),(4)),RSCRRI(100))                      
*FOX  ALDA(2,1)=ONE ;                                                   *FOX
      RSCRRI(100) = ONE                                                 
      CALL DACON(ALDA       ((2),(1)),RSCRRI(100))                      
*FOX  ALDA(2,2)=ZERO ;                                                  *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(ALDA       ((2),(2)),RSCRRI(100))                      
*FOX  ALDA(2,3)=-FOK ;                                                  *FOX
      CALL DACMU(FOK        ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),ALDA       ((2),(3)))                 
*FOX  ALDA(2,4)=ONE ;                                                   *FOX
      RSCRRI(100) = ONE                                                 
      CALL DACON(ALDA       ((2),(4)),RSCRRI(100))                      
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
        CALL DADAL(WFA     ,1)                                                  
        CALL DADAL(RHOI    ,1)                                                  
        CALL DADAL(AFOK    ,1)                                                  
        CALL DADAL(WF      ,1)                                                  
        CALL DADAL(FOKC    ,1)                                                  
        CALL DADAL(HS      ,1)                                                  
        CALL DADAL(HC      ,1)                                                  
        CALL DADAL(HM      ,1)                                                  
        CALL DADAL(HP      ,1)                                                  
        CALL DADAL(HI1     ,1)                                                  
        CALL DADAL(AEK     ,1)                                                  
        CALL DADAL(FI      ,1)                                                  
        CALL DADAL(HI      ,1)                                                  
        CALL DADAL(RHOC    ,1)                                                  
        CALL DADAL(SIQ     ,1)                                                  
        CALL DADAL(GL      ,1)                                                  
        CALL DADAL(G       ,1)                                                  
        CALL DADAL(CO      ,1)                                                  
        CALL DADAL(SI      ,1)                                                  
        CALL DADAL(AS6     ,1)                                                  
        CALL DADAL(AS4     ,1)                                                  
        CALL DADAL(AS3     ,1)                                                  
        CALL DADAL(SM23    ,1)                                                  
        CALL DADAL(SM12    ,1)                                                  
        CALL DADAL(SM6     ,1)                                                  
        CALL DADAL(SM5     ,1)                                                  
        CALL DADAL(SM4     ,1)                                                  
        CALL DADAL(SM3     ,1)                                                  
        CALL DADAL(SM2     ,1)                                                  
        CALL DADAL(SM1     ,1)                                                  
        CALL DADAL(FOK1    ,1)                                                  
        CALL DADAL(RHO     ,1)                                                  
        CALL DADAL(FOK     ,1)                                                  
        CALL DADAL(DPSQ    ,1)                                                  
        CALL DADAL(DPD     ,1)                                                  
        CALL DADAL(WFHI    ,1)                                                  
        CALL DADAL(FOKQ    ,1)                                                  
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
*FOX{
      INTEGER SIGMDA  
      INTEGER DPDA    
      INTEGER DPDA1   
      INTEGER RV      
      INTEGER XX      (2)
      INTEGER YY      (2)
      INTEGER EJ1     
      INTEGER EJF1    
      INTEGER ALDA    (2,6)
      INTEGER ASDA    (2,6)
      INTEGER ALDAQ   (2,6)
      INTEGER ASDAQ   (2,6)
      INTEGER SMIDA   (MCOR)
      INTEGER FOKQ    
      INTEGER WFHI    
      INTEGER DPD     
      INTEGER DPSQ    
      INTEGER FOK     
      INTEGER RHO     
      INTEGER FOK1    
      INTEGER SM1     
      INTEGER SM2     
      INTEGER SM3     
      INTEGER SM4     
      INTEGER SM5     
      INTEGER SM6     
      INTEGER SM12    
      INTEGER SM23    
      INTEGER AS3     
      INTEGER AS4     
      INTEGER AS6     
      INTEGER SI      
      INTEGER CO      
      INTEGER G       
      INTEGER GL      
      INTEGER SIQ     
      INTEGER RHOC    
      INTEGER HI      
      INTEGER FI      
      INTEGER AEK     
      INTEGER HI1     
      INTEGER HP      
      INTEGER HM      
      INTEGER HC      
      INTEGER HS      
      INTEGER FOKC    
      INTEGER WF      
      INTEGER AFOK    
      INTEGER RHOI    
      INTEGER WFA     
      INTEGER ISCRDA, ISCRRI,IDAO
      DOUBLE PRECISION RSCRRI
      COMMON/DASCR/ISCRDA(100),RSCRRI(100),ISCRRI(100),IDAO
      if(1.eq.1) then                                                 
         CALL DAKEY('FOX V2.1')
         CALL DAALL(FOKQ    ,1,'FOKQ      ',NORD,NVAR)
         CALL DAALL(WFHI    ,1,'WFHI      ',NORD,NVAR)
         CALL DAALL(DPD     ,1,'DPD       ',NORD,NVAR)
         CALL DAALL(DPSQ    ,1,'DPSQ      ',NORD,NVAR)
         CALL DAALL(FOK     ,1,'FOK       ',NORD,NVAR)
         CALL DAALL(RHO     ,1,'RHO       ',NORD,NVAR)
         CALL DAALL(FOK1    ,1,'FOK1      ',NORD,NVAR)
         CALL DAALL(SM1     ,1,'SM1       ',NORD,NVAR)
         CALL DAALL(SM2     ,1,'SM2       ',NORD,NVAR)
         CALL DAALL(SM3     ,1,'SM3       ',NORD,NVAR)
         CALL DAALL(SM4     ,1,'SM4       ',NORD,NVAR)
         CALL DAALL(SM5     ,1,'SM5       ',NORD,NVAR)
         CALL DAALL(SM6     ,1,'SM6       ',NORD,NVAR)
         CALL DAALL(SM12    ,1,'SM12      ',NORD,NVAR)
         CALL DAALL(SM23    ,1,'SM23      ',NORD,NVAR)
         CALL DAALL(AS3     ,1,'AS3       ',NORD,NVAR)
         CALL DAALL(AS4     ,1,'AS4       ',NORD,NVAR)
         CALL DAALL(AS6     ,1,'AS6       ',NORD,NVAR)
         CALL DAALL(SI      ,1,'SI        ',NORD,NVAR)
         CALL DAALL(CO      ,1,'CO        ',NORD,NVAR)
         CALL DAALL(G       ,1,'G         ',NORD,NVAR)
         CALL DAALL(GL      ,1,'GL        ',NORD,NVAR)
         CALL DAALL(SIQ     ,1,'SIQ       ',NORD,NVAR)
         CALL DAALL(RHOC    ,1,'RHOC      ',NORD,NVAR)
         CALL DAALL(HI      ,1,'HI        ',NORD,NVAR)
         CALL DAALL(FI      ,1,'FI        ',NORD,NVAR)
         CALL DAALL(AEK     ,1,'AEK       ',NORD,NVAR)
         CALL DAALL(HI1     ,1,'HI1       ',NORD,NVAR)
         CALL DAALL(HP      ,1,'HP        ',NORD,NVAR)
         CALL DAALL(HM      ,1,'HM        ',NORD,NVAR)
         CALL DAALL(HC      ,1,'HC        ',NORD,NVAR)
         CALL DAALL(HS      ,1,'HS        ',NORD,NVAR)
         CALL DAALL(FOKC    ,1,'FOKC      ',NORD,NVAR)
         CALL DAALL(WF      ,1,'WF        ',NORD,NVAR)
         CALL DAALL(AFOK    ,1,'AFOK      ',NORD,NVAR)
         CALL DAALL(RHOI    ,1,'RHOI      ',NORD,NVAR)
         CALL DAALL(WFA     ,1,'WFA       ',NORD,NVAR)
      ENDIF
      IDAA = IDAO
*FOX}
!-----------------------------------------------------------------------
*FOX  DPD=ONE+DPDA ;                                                    *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPD        )                          
*FOX  DPSQ=SQRT(DPD) ;                                                  *FOX
      CALL DAFUN('SQRT',DPD        ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),DPSQ       )                          
!-----------------------------------------------------------------------
!  QUADRUPOLE
!  FOCUSSING
!-----------------------------------------------------------------------
      if(abs(ek(i)).le.pieni) goto 100
*FOX  FOK=(SMIDA(IPCH)*RATIOE(I))/(ONE+DPDA) ;                          *FOX
      CALL DACOP(SMIDA      (IPCH       ),ISCRDA(  1+IDAA))             
      RSCRRI(  2+IDAA) = RATIOE     (I          )                       
      CALL DACMU(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  4+IDAA))          
      CALL DADIV(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(ISCRDA(  5+IDAA),FOK        )                          
*FOX  AEK=FOK ;                                                         *FOX
      CALL DACOP(FOK        ,AEK        )                               
      if(dare(aek).lt.zero) then
*FOX  AEK=-AEK ;                                                        *FOX
      CALL DACMU(AEK        ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),AEK        )                          
      endif
      ih=0
*FOX  HI=SQRT(AEK) ;                                                    *FOX
      CALL DAFUN('SQRT',AEK        ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),HI         )                          
*FOX  FI=EL(I)*HI ;                                                     *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACMU(HI         ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),FI         )                          
      if(ek(i).gt.zero) goto 30
   20 ih=ih+1
*FOX  ALDAQ(IH,1)=COS(FI) ;                                             *FOX
      CALL DAFUN('COS ',FI         ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),ALDAQ      (IH         ,(1)))         
*FOX  HI1=SIN(FI) ;                                                     *FOX
      CALL DAFUN('SIN ',FI         ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),HI1        )                          
*FOX  ALDAQ(IH,2)=HI1/HI ;                                              *FOX
      CALL DADIV(HI1        ,HI         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ALDAQ      (IH         ,(2)))         
*FOX  ALDAQ(IH,3)=-HI1*HI ;                                             *FOX
      CALL DACMU(HI1        ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),HI         ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),ALDAQ      (IH         ,(3)))         
*FOX  ALDAQ(IH,4)=ALDAQ(IH,1) ;                                         *FOX
      CALL DACOP(ALDAQ      (IH         ,(1)),ISCRDA(  1+IDAA))         
      CALL DACOP(ISCRDA(  1+IDAA),ALDAQ      (IH         ,(4)))         
*FOX  ASDAQ(IH,1)=EL(I)*(ONE-RV)*C1E3 ;                                 *FOX
      CALL DASUC(RV         ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      RSCRRI(  2+IDAA) = EL         (I          )                       
      CALL DACMU(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  3+IDAA),ONE*C1E3       ,ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),ASDAQ      (IH         ,(1)))         
*FOX  ASDAQ(IH,4)=-RV*ALDAQ(IH,2)*ALDAQ(IH,3)/C2E3 ;                    *FOX
      CALL DACOP(ALDAQ      (IH         ,(2)),ISCRDA(  1+IDAA))         
      CALL DACOP(ALDAQ      (IH         ,(3)),ISCRDA(  2+IDAA))         
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  3+IDAA))        
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  1+IDAA),ISCRDA(  4+IDAA))    
      CALL DAMUL(ISCRDA(  4+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DACDI(ISCRDA(  5+IDAA),ONE*C2E3       ,ISCRDA(  6+IDAA))     
      CALL DACOP(ISCRDA(  6+IDAA),ASDAQ      (IH         ,(4)))         
*FOX  ASDAQ(IH,5)=-RV*(EL(I)-ALDAQ(IH,1)*ALDAQ(IH,2))*AEK/C4E3 ;        *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDAQ      (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDAQ      (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DASUC(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  7+IDAA),AEK        ,ISCRDA(  8+IDAA))         
      CALL DACDI(ISCRDA(  8+IDAA),ONE*C4E3       ,ISCRDA(  9+IDAA))     
      CALL DACOP(ISCRDA(  9+IDAA),ASDAQ      (IH         ,(5)))         
*FOX  ASDAQ(IH,6)=-RV*(EL(I)+ALDAQ(IH,1)*ALDAQ(IH,2))/C4E3 ;            *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDAQ      (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDAQ      (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACAD(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DACDI(ISCRDA(  7+IDAA),ONE*C4E3       ,ISCRDA(  8+IDAA))     
      CALL DACOP(ISCRDA(  8+IDAA),ASDAQ      (IH         ,(6)))         
      if(ih.eq.2) goto 100
!--DEFOCUSSING
   30 ih=ih+1
*FOX  HP=EXP(FI) ;                                                      *FOX
      CALL DAFUN('EXP ',FI         ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),HP         )                          
*FOX  HM=ONE/HP ;                                                       *FOX
      CALL DADIC(HP         ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),HM         )                          
*FOX  HC=(HP+HM)*HALF ;                                                 *FOX
      CALL DAADD(HP         ,HM         ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*HALF       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),HC         )                          
*FOX  HS=(HP-HM)*HALF ;                                                 *FOX
      CALL DASUB(HP         ,HM         ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*HALF       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),HS         )                          
*FOX  ALDAQ(IH,1)=HC ;                                                  *FOX
      CALL DACOP(HC         ,ALDAQ      (IH         ,(1)))              
*FOX  ALDAQ(IH,2)=HS/HI ;                                               *FOX
      CALL DADIV(HS         ,HI         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ALDAQ      (IH         ,(2)))         
*FOX  ALDAQ(IH,3)=HS*HI ;                                               *FOX
      CALL DAMUL(HS         ,HI         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ALDAQ      (IH         ,(3)))         
*FOX  ALDAQ(IH,4)=HC ;                                                  *FOX
      CALL DACOP(HC         ,ALDAQ      (IH         ,(4)))              
*FOX  ASDAQ(IH,4)=-RV*ALDAQ(IH,2)*ALDAQ(IH,3)/C2E3 ;                    *FOX
      CALL DACOP(ALDAQ      (IH         ,(2)),ISCRDA(  1+IDAA))         
      CALL DACOP(ALDAQ      (IH         ,(3)),ISCRDA(  2+IDAA))         
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  3+IDAA))        
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  1+IDAA),ISCRDA(  4+IDAA))    
      CALL DAMUL(ISCRDA(  4+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DACDI(ISCRDA(  5+IDAA),ONE*C2E3       ,ISCRDA(  6+IDAA))     
      CALL DACOP(ISCRDA(  6+IDAA),ASDAQ      (IH         ,(4)))         
*FOX  ASDAQ(IH,5)=+RV*(EL(I)-ALDAQ(IH,1)*ALDAQ(IH,2))*AEK/C4E3 ;        *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDAQ      (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDAQ      (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DASUC(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DAMUL(RV         ,ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))         
      CALL DAMUL(ISCRDA(  6+IDAA),AEK        ,ISCRDA(  7+IDAA))         
      CALL DACDI(ISCRDA(  7+IDAA),ONE*C4E3       ,ISCRDA(  8+IDAA))     
      CALL DACOP(ISCRDA(  8+IDAA),ASDAQ      (IH         ,(5)))         
*FOX  ASDAQ(IH,6)=-RV*(EL(I)+ALDAQ(IH,1)*ALDAQ(IH,2))/C4E3 ;            *FOX
      RSCRRI(  1+IDAA) = EL         (I          )                       
      CALL DACOP(ALDAQ      (IH         ,(1)),ISCRDA(  2+IDAA))         
      CALL DACOP(ALDAQ      (IH         ,(2)),ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACAD(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(RV         ,ONE*(-ONE       ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DACDI(ISCRDA(  7+IDAA),ONE*C4E3       ,ISCRDA(  8+IDAA))     
      CALL DACOP(ISCRDA(  8+IDAA),ASDAQ      (IH         ,(6)))         
      if(ih.eq.1) goto 20
  100 continue
        CALL DADAL(WFA     ,1)                                                  
        CALL DADAL(RHOI    ,1)                                                  
        CALL DADAL(AFOK    ,1)                                                  
        CALL DADAL(WF      ,1)                                                  
        CALL DADAL(FOKC    ,1)                                                  
        CALL DADAL(HS      ,1)                                                  
        CALL DADAL(HC      ,1)                                                  
        CALL DADAL(HM      ,1)                                                  
        CALL DADAL(HP      ,1)                                                  
        CALL DADAL(HI1     ,1)                                                  
        CALL DADAL(AEK     ,1)                                                  
        CALL DADAL(FI      ,1)                                                  
        CALL DADAL(HI      ,1)                                                  
        CALL DADAL(RHOC    ,1)                                                  
        CALL DADAL(SIQ     ,1)                                                  
        CALL DADAL(GL      ,1)                                                  
        CALL DADAL(G       ,1)                                                  
        CALL DADAL(CO      ,1)                                                  
        CALL DADAL(SI      ,1)                                                  
        CALL DADAL(AS6     ,1)                                                  
        CALL DADAL(AS4     ,1)                                                  
        CALL DADAL(AS3     ,1)                                                  
        CALL DADAL(SM23    ,1)                                                  
        CALL DADAL(SM12    ,1)                                                  
        CALL DADAL(SM6     ,1)                                                  
        CALL DADAL(SM5     ,1)                                                  
        CALL DADAL(SM4     ,1)                                                  
        CALL DADAL(SM3     ,1)                                                  
        CALL DADAL(SM2     ,1)                                                  
        CALL DADAL(SM1     ,1)                                                  
        CALL DADAL(FOK1    ,1)                                                  
        CALL DADAL(RHO     ,1)                                                  
        CALL DADAL(FOK     ,1)                                                  
        CALL DADAL(DPSQ    ,1)                                                  
        CALL DADAL(DPD     ,1)                                                  
        CALL DADAL(WFHI    ,1)                                                  
        CALL DADAL(FOKQ    ,1)                                                  
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
      double precision al,as,chi0,chid,dp1,dps,exz,sigm
      common/syos/as(6,2,npart,nele),al(6,2,npart,nele),sigm(mpa),      &
     &dps(mpa),idz(2)
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
*FOX{
      INTEGER SIGMDA  
      INTEGER DPDA    
      INTEGER DPDA1   
      INTEGER RV      
      INTEGER XX      (2)
      INTEGER YY      (2)
      INTEGER EJ1     
      INTEGER EJF1    
      INTEGER ALDA    (2,6)
      INTEGER ASDA    (2,6)
      INTEGER ALDAQ   (2,6)
      INTEGER ASDAQ   (2,6)
      INTEGER SMIDA   (MCOR)
      INTEGER X       (2)
      INTEGER Y       (2)
      INTEGER YP      (2)
      INTEGER DKIP    
      INTEGER CORROLD (MCOP)
      INTEGER CORRNEW (MCOP)
      INTEGER CORRAU1 (MCOP)
      INTEGER CORRAU2 (MCOP)
      INTEGER AA      (11)
      INTEGER BB      (11)
      INTEGER TRACKI  (6)
      INTEGER PUX     
      INTEGER PUZ     
      INTEGER EJF0    
      INTEGER EKK     
      INTEGER XL      
      INTEGER ZL      
      INTEGER CRKVE   
      INTEGER CIKVE   
      INTEGER CRKVEUK 
      INTEGER CBZBF   
      INTEGER YV1J    
      INTEGER YV2J    
      INTEGER CRKVEBF 
      INTEGER CIKVEBF 
      INTEGER RHO2BF  
      INTEGER TKBF    
      INTEGER XRBF    
      INTEGER CCCC    
      INTEGER ZRBF    
      INTEGER XBBF    
      INTEGER ZBBF    
      INTEGER CRXBF   
      INTEGER CBXBF   
      INTEGER CRZBF   
      INTEGER WX      
      INTEGER WY      
      INTEGER CRABAMP 
      INTEGER CRABAMP2
      INTEGER CRABAMP3
      INTEGER CRABAMP4
      INTEGER PZ      
      INTEGER ISCRDA, ISCRRI,IDAO
      DOUBLE PRECISION RSCRRI
      COMMON/DASCR/ISCRDA(100),RSCRRI(100),ISCRRI(100),IDAO
      if(1.eq.1) then                                                 
         CALL DAKEY('FOX V2.1')
         CALL DAALL(SIGMDA  ,1,'SIGMDA    ',NORD,NVAR)
         CALL DAALL(DPDA    ,1,'DPDA      ',NORD,NVAR)
         CALL DAALL(DPDA1   ,1,'DPDA1     ',NORD,NVAR)
         CALL DAALL(RV      ,1,'RV        ',NORD,NVAR)
         CALL DAALL(XX      ,1*(2),'XX        ',NORD,NVAR)
         CALL DAALL(YY      ,1*(2),'YY        ',NORD,NVAR)
         CALL DAALL(EJ1     ,1,'EJ1       ',NORD,NVAR)
         CALL DAALL(EJF1    ,1,'EJF1      ',NORD,NVAR)
         CALL DAALL(ALDA    ,1*(2)*(6),'ALDA      ',NORD,NVAR)
         CALL DAALL(ASDA    ,1*(2)*(6),'ASDA      ',NORD,NVAR)
         CALL DAALL(ALDAQ   ,1*(2)*(6),'ALDAQ     ',NORD,NVAR)
         CALL DAALL(ASDAQ   ,1*(2)*(6),'ASDAQ     ',NORD,NVAR)
         CALL DAALL(SMIDA   ,1*(MCOR),'SMIDA     ',NORD,NVAR)
         CALL DAALL(X       ,1*(2),'X         ',NORD,NVAR)
         CALL DAALL(Y       ,1*(2),'Y         ',NORD,NVAR)
         CALL DAALL(YP      ,1*(2),'YP        ',NORD,NVAR)
         CALL DAALL(DKIP    ,1,'DKIP      ',NORD,NVAR)
         CALL DAALL(CORROLD ,1*(MCOP),'CORROLD   ',NORD,NVAR)
         CALL DAALL(CORRNEW ,1*(MCOP),'CORRNEW   ',NORD,NVAR)
         CALL DAALL(CORRAU1 ,1*(MCOP),'CORRAU1   ',NORD,NVAR)
         CALL DAALL(CORRAU2 ,1*(MCOP),'CORRAU2   ',NORD,NVAR)
         CALL DAALL(AA      ,1*(11),'AA        ',NORD,NVAR)
         CALL DAALL(BB      ,1*(11),'BB        ',NORD,NVAR)
         CALL DAALL(TRACKI  ,1*(6),'TRACKI    ',NORD,NVAR)
         CALL DAALL(PUX     ,1,'PUX       ',NORD,NVAR)
         CALL DAALL(PUZ     ,1,'PUZ       ',NORD,NVAR)
         CALL DAALL(EJF0    ,1,'EJF0      ',NORD,NVAR)
         CALL DAALL(EKK     ,1,'EKK       ',NORD,NVAR)
         CALL DAALL(XL      ,1,'XL        ',NORD,NVAR)
         CALL DAALL(ZL      ,1,'ZL        ',NORD,NVAR)
         CALL DAALL(CRKVE   ,1,'CRKVE     ',NORD,NVAR)
         CALL DAALL(CIKVE   ,1,'CIKVE     ',NORD,NVAR)
         CALL DAALL(CRKVEUK ,1,'CRKVEUK   ',NORD,NVAR)
         CALL DAALL(CBZBF   ,1,'CBZBF     ',NORD,NVAR)
         CALL DAALL(YV1J    ,1,'YV1J      ',NORD,NVAR)
         CALL DAALL(YV2J    ,1,'YV2J      ',NORD,NVAR)
         CALL DAALL(CRKVEBF ,1,'CRKVEBF   ',NORD,NVAR)
         CALL DAALL(CIKVEBF ,1,'CIKVEBF   ',NORD,NVAR)
         CALL DAALL(RHO2BF  ,1,'RHO2BF    ',NORD,NVAR)
         CALL DAALL(TKBF    ,1,'TKBF      ',NORD,NVAR)
         CALL DAALL(XRBF    ,1,'XRBF      ',NORD,NVAR)
         CALL DAALL(CCCC    ,1,'CCCC      ',NORD,NVAR)
         CALL DAALL(ZRBF    ,1,'ZRBF      ',NORD,NVAR)
         CALL DAALL(XBBF    ,1,'XBBF      ',NORD,NVAR)
         CALL DAALL(ZBBF    ,1,'ZBBF      ',NORD,NVAR)
         CALL DAALL(CRXBF   ,1,'CRXBF     ',NORD,NVAR)
         CALL DAALL(CBXBF   ,1,'CBXBF     ',NORD,NVAR)
         CALL DAALL(CRZBF   ,1,'CRZBF     ',NORD,NVAR)
         CALL DAALL(WX      ,1,'WX        ',NORD,NVAR)
         CALL DAALL(WY      ,1,'WY        ',NORD,NVAR)
         CALL DAALL(CRABAMP ,1,'CRABAMP   ',NORD,NVAR)
         CALL DAALL(CRABAMP2,1,'CRABAMP2  ',NORD,NVAR)
         CALL DAALL(CRABAMP3,1,'CRABAMP3  ',NORD,NVAR)
         CALL DAALL(CRABAMP4,1,'CRABAMP4  ',NORD,NVAR)
         CALL DAALL(PZ      ,1,'PZ        ',NORD,NVAR)
      ENDIF
      IDAA = IDAO
*FOX}
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
*FOX  X(2)=OZ ;                                                         *FOX
      RSCRRI(100) = OZ                                                  
      CALL DACON(X          ((2)),RSCRRI(100))                          
*FOX  YP(2)=OZP*(ONE+DPS1) ;                                            *FOX
      RSCRRI(  1+IDAA) = ONE         + DPS1                             
      RSCRRI(  2+IDAA) = OZP         * RSCRRI(  1+IDAA)                 
      RSCRRI(100) = RSCRRI(  2+IDAA)                                    
      CALL DACON(YP         ((2)),RSCRRI(100))                          
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
*FOX  SIGMDA=SIGM1 ;                                                    *FOX
      RSCRRI(100) = SIGM1                                               
      CALL DACON(SIGMDA     ,RSCRRI(100))                               
*FOX  DPDA1=DPS1*C1E3 ;                                                 *FOX
      RSCRRI(  1+IDAA) = DPS1        * C1E3                             
      RSCRRI(100) = RSCRRI(  1+IDAA)                                    
      CALL DACON(DPDA1      ,RSCRRI(100))                               
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
*FOX  DPDA=DPDA1*C1M3 ;                                                 *FOX
      CALL DACMU(DPDA1      ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA       )                          
*FOX  CORROLD(1)=X(1) ;                                                 *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORROLD    ((1)))                     
*FOX  CORROLD(2)=YP(1) ;                                                *FOX
      CALL DACOP(YP         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORROLD    ((2)))                     
*FOX  CORROLD(3)=X(2) ;                                                 *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORROLD    ((3)))                     
*FOX  CORROLD(4)=YP(2) ;                                                *FOX
      CALL DACOP(YP         ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORROLD    ((4)))                     
*FOX  CORROLD(5)=SIGMDA ;                                               *FOX
      CALL DACOP(SIGMDA     ,CORROLD    ((5)))                          
*FOX  CORROLD(6)=DPDA1 ;                                                *FOX
      CALL DACOP(DPDA1      ,CORROLD    ((6)))                          
            do 5 kkk=1,6
              dpdav=dare(corrold(kkk))
*FOX  CORROLD(KKK)=CORROLD(KKK)-DPDAV ;                                 *FOX
      CALL DACOP(CORROLD    (KKK        ),ISCRDA(  1+IDAA))             
      CALL DACSU(ISCRDA(  1+IDAA),ONE*DPDAV      ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),CORROLD    (KKK        ))             
    5       continue
*FOX  Y(1)=YP(1)/(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(YP         ((1)),ISCRDA(  2+IDAA))                     
      CALL DADIV(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((1)))                     
*FOX  Y(2)=YP(2)/(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(YP         ((2)),ISCRDA(  2+IDAA))                     
      CALL DADIV(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((2)))                     
      iflag=0
      iflag1=0
      iflag2=0
*FOX  EJF1=E0F*(ONE+DPDA) ;                                             *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),EJF1       )                          
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;                                     *FOX
      CALL DAMUL(EJF1       ,EJF1       ,ISCRDA(  1+IDAA))              
      RSCRRI(  2+IDAA) = PMA         * PMA                              
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DAFUN('SQRT',ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))              
      CALL DACOP(ISCRDA(  4+IDAA),EJ1        )                          
*FOX  RV=EJ1/E0*E0F/EJF1 ;                                              *FOX
      CALL DACDI(EJ1        ,ONE*E0         ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DADIV(ISCRDA(  2+IDAA),EJF1       ,ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),RV         )                          
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
*FOX  YP(1)=Y(1)*(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((1)))                     
*FOX  YP(2)=Y(2)*(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((2)))                     
*FOX  DPDA1=DPDA*C1E3 ;                                                 *FOX
      CALL DACMU(DPDA       ,ONE*C1E3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA1      )                          
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
*FOX  EJF1=E0F*(ONE+DPDA) ;                                             *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),EJF1       )                          
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;                                     *FOX
      CALL DAMUL(EJF1       ,EJF1       ,ISCRDA(  1+IDAA))              
      RSCRRI(  2+IDAA) = PMA         * PMA                              
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DAFUN('SQRT',ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))              
      CALL DACOP(ISCRDA(  4+IDAA),EJ1        )                          
*FOX  RV=EJ1/E0*E0F/EJF1 ;                                              *FOX
      CALL DACDI(EJ1        ,ONE*E0         ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DADIV(ISCRDA(  2+IDAA),EJF1       ,ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),RV         )                          
          if(ithick.eq.1) then
*FOX  YP(1)=Y(1)*(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((1)))                     
*FOX  YP(2)=Y(2)*(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((2)))                     
            if(icav.eq.0) then
*FOX  CORRNEW(1)=X(1) ;                                                 *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRNEW    ((1)))                     
*FOX  CORRNEW(2)=YP(1) ;                                                *FOX
      CALL DACOP(YP         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRNEW    ((2)))                     
*FOX  CORRNEW(3)=X(2) ;                                                 *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRNEW    ((3)))                     
*FOX  CORRNEW(4)=YP(2) ;                                                *FOX
      CALL DACOP(YP         ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRNEW    ((4)))                     
*FOX  CORRNEW(5)=SIGMDA ;                                               *FOX
      CALL DACOP(SIGMDA     ,CORRNEW    ((5)))                          
*FOX  CORRNEW(6)=DPDA1 ;                                                *FOX
      CALL DACOP(DPDA1      ,CORRNEW    ((6)))                          
              do 24 kkk=1,6
                dpdav=dare(corrnew(kkk))
*FOX  CORRNEW(KKK)=CORRNEW(KKK)-DPDAV ;                                 *FOX
      CALL DACOP(CORRNEW    (KKK        ),ISCRDA(  1+IDAA))             
      CALL DACSU(ISCRDA(  1+IDAA),ONE*DPDAV      ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),CORRNEW    (KKK        ))             
   24         continue
            else
*FOX  CORRAU2(1)=X(1) ;                                                 *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRAU2    ((1)))                     
*FOX  CORRAU2(2)=YP(1) ;                                                *FOX
      CALL DACOP(YP         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRAU2    ((2)))                     
*FOX  CORRAU2(3)=X(2) ;                                                 *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRAU2    ((3)))                     
*FOX  CORRAU2(4)=YP(2) ;                                                *FOX
      CALL DACOP(YP         ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRAU2    ((4)))                     
*FOX  CORRAU2(5)=SIGMDA ;                                               *FOX
      CALL DACOP(SIGMDA     ,CORRAU2    ((5)))                          
*FOX  CORRAU2(6)=DPDA1 ;                                                *FOX
      CALL DACOP(DPDA1      ,CORRAU2    ((6)))                          
              do 25 kkk=1,6
*FOX  CORRAU1(KKK)=CORRNEW(KKK) ;                                       *FOX
      CALL DACOP(CORRNEW    (KKK        ),ISCRDA(  1+IDAA))             
      CALL DACOP(ISCRDA(  1+IDAA),CORRAU1    (KKK        ))             
                dpdav=dare(corrau2(kkk))
*FOX  CORRAU2(KKK)=CORRAU2(KKK)-DPDAV ;                                 *FOX
      CALL DACOP(CORRAU2    (KKK        ),ISCRDA(  1+IDAA))             
      CALL DACSU(ISCRDA(  1+IDAA),ONE*DPDAV      ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),CORRAU2    (KKK        ))             
   25         continue
              if(ivar.gt.ivar1) then
*FOX  CORRAU2(7)=SMIDA(1) ;                                             *FOX
      CALL DACOP(SMIDA      ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRAU2    ((7)))                     
*FOX  CORRAU2(8)=SMIDA(2) ;                                             *FOX
      CALL DACOP(SMIDA      ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRAU2    ((8)))                     
                dpdav=dare(smida(1))
*FOX  CORRAU1(7)=SMIDA(1)-DPDAV ;                                       *FOX
      CALL DACOP(SMIDA      ((1)),ISCRDA(  1+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*DPDAV      ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),CORRAU1    ((7)))                     
                dpdav=dare(smida(2))
*FOX  CORRAU1(8)=SMIDA(2)-DPDAV ;                                       *FOX
      CALL DACOP(SMIDA      ((2)),ISCRDA(  1+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*DPDAV      ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),CORRAU1    ((8)))                     
              endif
              call dacct(corrau2,nvar,corrau1,nvar,corrnew,nvar)
            endif
            dpdav=dare(x(1))
*FOX  X(1)=CORROLD(1)+DPDAV ;                                           *FOX
      CALL DACOP(CORROLD    ((1)),ISCRDA(  1+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*DPDAV      ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),X          ((1)))                     
            dpdav=dare(yp(1))
*FOX  YP(1)=CORROLD(2)+DPDAV ;                                          *FOX
      CALL DACOP(CORROLD    ((2)),ISCRDA(  1+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*DPDAV      ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),YP         ((1)))                     
            dpdav=dare(x(2))
*FOX  X(2)=CORROLD(3)+DPDAV ;                                           *FOX
      CALL DACOP(CORROLD    ((3)),ISCRDA(  1+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*DPDAV      ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),X          ((2)))                     
            dpdav=dare(yp(2))
*FOX  YP(2)=CORROLD(4)+DPDAV ;                                          *FOX
      CALL DACOP(CORROLD    ((4)),ISCRDA(  1+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*DPDAV      ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),YP         ((2)))                     
            dpdav=dare(sigmda)
*FOX  SIGMDA=CORROLD(5)+DPDAV ;                                         *FOX
      CALL DACOP(CORROLD    ((5)),ISCRDA(  1+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*DPDAV      ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),SIGMDA     )                          
            dpdav=dare(dpda1)
*FOX  DPDA1=CORROLD(6)+DPDAV ;                                          *FOX
      CALL DACOP(CORROLD    ((6)),ISCRDA(  1+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*DPDAV      ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),DPDA1      )                          
*FOX  DPDA=DPDA1*C1M3 ;                                                 *FOX
      CALL DACMU(DPDA1      ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA       )                          
*FOX  Y(1)=YP(1)/(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(YP         ((1)),ISCRDA(  2+IDAA))                     
      CALL DADIV(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((1)))                     
*FOX  Y(2)=YP(2)/(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(YP         ((2)),ISCRDA(  2+IDAA))                     
      CALL DADIV(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((2)))                     
*FOX  EJF1=E0F*(ONE+DPDA) ;                                             *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),EJF1       )                          
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;                                     *FOX
      CALL DAMUL(EJF1       ,EJF1       ,ISCRDA(  1+IDAA))              
      RSCRRI(  2+IDAA) = PMA         * PMA                              
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DAFUN('SQRT',ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))              
      CALL DACOP(ISCRDA(  4+IDAA),EJ1        )                          
*FOX  RV=EJ1/E0*E0F/EJF1 ;                                              *FOX
      CALL DACDI(EJ1        ,ONE*E0         ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DADIV(ISCRDA(  2+IDAA),EJF1       ,ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),RV         )                          
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
*FOX  PUX=X(1) ;                                                        *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUX        )                          
*FOX  PUZ=Y(1) ;                                                        *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUZ        )                          
*FOX  X(1)=ALDAQ(1,1)*PUX+ALDAQ(1,2)*PUZ+ALDAQ(1,5)*IDZ(1) ;            *FOX
      CALL DACOP(ALDAQ      ((1),(1)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDAQ      ((1),(2)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDAQ      ((1),(5)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((1))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),X          ((1)))                     
*FOX  Y(1)=ALDAQ(1,3)*PUX+ALDAQ(1,4)*PUZ+ALDAQ(1,6)*IDZ(1) ;            *FOX
      CALL DACOP(ALDAQ      ((1),(3)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDAQ      ((1),(4)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDAQ      ((1),(6)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((1))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((1)))                     
*FOX  PUX=X(2) ;                                                        *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUX        )                          
*FOX  PUZ=Y(2) ;                                                        *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUZ        )                          
*FOX  X(2)=ALDAQ(2,1)*PUX+ALDAQ(2,2)*PUZ+ALDAQ(2,5)*IDZ(2) ;            *FOX
      CALL DACOP(ALDAQ      ((2),(1)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDAQ      ((2),(2)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDAQ      ((2),(5)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((2))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),X          ((2)))                     
*FOX  Y(2)=ALDAQ(2,3)*PUX+ALDAQ(2,4)*PUZ+ALDAQ(2,6)*IDZ(2) ;            *FOX
      CALL DACOP(ALDAQ      ((2),(3)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDAQ      ((2),(4)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDAQ      ((2),(6)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((2))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((2)))                     
              else
*FOX  PUX=X(1) ;                                                        *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUX        )                          
*FOX  PUZ=Y(1) ;                                                        *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUZ        )                          
*FOX  X(1)=ALDA(1,1)*PUX+ALDA(1,2)*PUZ+ALDA(1,5)*IDZ(1) ;               *FOX
      CALL DACOP(ALDA       ((1),(1)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDA       ((1),(2)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDA       ((1),(5)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((1))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),X          ((1)))                     
*FOX  Y(1)=ALDA(1,3)*PUX+ALDA(1,4)*PUZ+ALDA(1,6)*IDZ(1) ;               *FOX
      CALL DACOP(ALDA       ((1),(3)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDA       ((1),(4)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDA       ((1),(6)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((1))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((1)))                     
*FOX  PUX=X(2) ;                                                        *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUX        )                          
*FOX  PUZ=Y(2) ;                                                        *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUZ        )                          
*FOX  X(2)=ALDA(2,1)*PUX+ALDA(2,2)*PUZ+ALDA(2,5)*IDZ(2) ;               *FOX
      CALL DACOP(ALDA       ((2),(1)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDA       ((2),(2)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDA       ((2),(5)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((2))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),X          ((2)))                     
*FOX  Y(2)=ALDA(2,3)*PUX+ALDA(2,4)*PUZ+ALDA(2,6)*IDZ(2) ;               *FOX
      CALL DACOP(ALDA       ((2),(3)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDA       ((2),(4)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDA       ((2),(6)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((2))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((2)))                     
              endif
            enddo
          else
*FOX  X(1)=X(1)+BL1(IX,1,2)*Y(1) ;                                      *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = BL1        (IX         ,(1),(2))               
      CALL DACOP(Y          ((1)),ISCRDA(  3+IDAA))                     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(ISCRDA(  5+IDAA),X          ((1)))                     
*FOX  X(2)=X(2)+BL1(IX,2,2)*Y(2) ;                                      *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = BL1        (IX         ,(2),(2))               
      CALL DACOP(Y          ((2)),ISCRDA(  3+IDAA))                     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(ISCRDA(  5+IDAA),X          ((2)))                     
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
*FOX  YP(1)=Y(1)*(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((1)))                     
*FOX  YP(2)=Y(2)*(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((2)))                     
*FOX  DPDA1=DPDA*C1E3 ;                                                 *FOX
      CALL DACMU(DPDA       ,ONE*C1E3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA1      )                          
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
*FOX  PUX=X(1) ;                                                        *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUX        )                          
*FOX  PUZ=Y(1) ;                                                        *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUZ        )                          
*FOX  SIGMDA=SIGMDA+ASDAQ(1,1)+ASDAQ(1,2)*PUX+                          *FOX
*FOX  ASDAQ(1,3)*PUZ+ASDAQ(1,4)*PUX*PUZ+ASDAQ(1,5)*PUX*PUX+             *FOX
*FOX  ASDAQ(1,6)*PUZ*PUZ ;                                              *FOX
      CALL DACOP(ASDAQ      ((1),(1)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ASDAQ      ((1),(2)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ASDAQ      ((1),(3)),ISCRDA(  3+IDAA))                 
      CALL DACOP(ASDAQ      ((1),(4)),ISCRDA(  4+IDAA))                 
      CALL DACOP(ASDAQ      ((1),(5)),ISCRDA(  5+IDAA))                 
      CALL DACOP(ASDAQ      ((1),(6)),ISCRDA(  6+IDAA))                 
      CALL DAMUL(ISCRDA(  2+IDAA),PUX        ,ISCRDA(  7+IDAA))         
      CALL DAMUL(ISCRDA(  3+IDAA),PUZ        ,ISCRDA(  8+IDAA))         
      CALL DAMUL(ISCRDA(  4+IDAA),PUX        ,ISCRDA(  9+IDAA))         
      CALL DAMUL(ISCRDA(  9+IDAA),PUZ        ,ISCRDA( 10+IDAA))         
      CALL DAMUL(ISCRDA(  5+IDAA),PUX        ,ISCRDA( 11+IDAA))         
      CALL DAMUL(ISCRDA( 11+IDAA),PUX        ,ISCRDA( 12+IDAA))         
      CALL DAMUL(ISCRDA(  6+IDAA),PUZ        ,ISCRDA( 13+IDAA))         
      CALL DAMUL(ISCRDA( 13+IDAA),PUZ        ,ISCRDA( 14+IDAA))         
      CALL DAADD(SIGMDA     ,ISCRDA(  1+IDAA),ISCRDA( 15+IDAA))         
      CALL DAADD(ISCRDA( 15+IDAA),ISCRDA(  7+IDAA),ISCRDA( 16+IDAA))    
      CALL DAADD(ISCRDA( 16+IDAA),ISCRDA(  8+IDAA),ISCRDA( 17+IDAA))    
      CALL DAADD(ISCRDA( 17+IDAA),ISCRDA( 10+IDAA),ISCRDA( 18+IDAA))    
      CALL DAADD(ISCRDA( 18+IDAA),ISCRDA( 12+IDAA),ISCRDA( 19+IDAA))    
      CALL DAADD(ISCRDA( 19+IDAA),ISCRDA( 14+IDAA),ISCRDA( 20+IDAA))    
      CALL DACOP(ISCRDA( 20+IDAA),SIGMDA     )                          
*FOX  X(1)=ALDAQ(1,1)*PUX+ALDAQ(1,2)*PUZ+ALDAQ(1,5)*IDZ(1) ;            *FOX
      CALL DACOP(ALDAQ      ((1),(1)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDAQ      ((1),(2)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDAQ      ((1),(5)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((1))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),X          ((1)))                     
*FOX  Y(1)=ALDAQ(1,3)*PUX+ALDAQ(1,4)*PUZ+ALDAQ(1,6)*IDZ(1) ;            *FOX
      CALL DACOP(ALDAQ      ((1),(3)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDAQ      ((1),(4)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDAQ      ((1),(6)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((1))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((1)))                     
*FOX  PUX=X(2) ;                                                        *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUX        )                          
*FOX  PUZ=Y(2) ;                                                        *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUZ        )                          
*FOX  SIGMDA=SIGMDA+ASDAQ(2,1)+ASDAQ(2,2)*PUX+                          *FOX
*FOX  ASDAQ(2,3)*PUZ+ASDAQ(2,4)*PUX*PUZ+ASDAQ(2,5)*PUX*PUX+             *FOX
*FOX  ASDAQ(2,6)*PUZ*PUZ ;                                              *FOX
      CALL DACOP(ASDAQ      ((2),(1)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ASDAQ      ((2),(2)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ASDAQ      ((2),(3)),ISCRDA(  3+IDAA))                 
      CALL DACOP(ASDAQ      ((2),(4)),ISCRDA(  4+IDAA))                 
      CALL DACOP(ASDAQ      ((2),(5)),ISCRDA(  5+IDAA))                 
      CALL DACOP(ASDAQ      ((2),(6)),ISCRDA(  6+IDAA))                 
      CALL DAMUL(ISCRDA(  2+IDAA),PUX        ,ISCRDA(  7+IDAA))         
      CALL DAMUL(ISCRDA(  3+IDAA),PUZ        ,ISCRDA(  8+IDAA))         
      CALL DAMUL(ISCRDA(  4+IDAA),PUX        ,ISCRDA(  9+IDAA))         
      CALL DAMUL(ISCRDA(  9+IDAA),PUZ        ,ISCRDA( 10+IDAA))         
      CALL DAMUL(ISCRDA(  5+IDAA),PUX        ,ISCRDA( 11+IDAA))         
      CALL DAMUL(ISCRDA( 11+IDAA),PUX        ,ISCRDA( 12+IDAA))         
      CALL DAMUL(ISCRDA(  6+IDAA),PUZ        ,ISCRDA( 13+IDAA))         
      CALL DAMUL(ISCRDA( 13+IDAA),PUZ        ,ISCRDA( 14+IDAA))         
      CALL DAADD(SIGMDA     ,ISCRDA(  1+IDAA),ISCRDA( 15+IDAA))         
      CALL DAADD(ISCRDA( 15+IDAA),ISCRDA(  7+IDAA),ISCRDA( 16+IDAA))    
      CALL DAADD(ISCRDA( 16+IDAA),ISCRDA(  8+IDAA),ISCRDA( 17+IDAA))    
      CALL DAADD(ISCRDA( 17+IDAA),ISCRDA( 10+IDAA),ISCRDA( 18+IDAA))    
      CALL DAADD(ISCRDA( 18+IDAA),ISCRDA( 12+IDAA),ISCRDA( 19+IDAA))    
      CALL DAADD(ISCRDA( 19+IDAA),ISCRDA( 14+IDAA),ISCRDA( 20+IDAA))    
      CALL DACOP(ISCRDA( 20+IDAA),SIGMDA     )                          
*FOX  X(2)=ALDAQ(2,1)*PUX+ALDAQ(2,2)*PUZ+ALDAQ(2,5)*IDZ(2) ;            *FOX
      CALL DACOP(ALDAQ      ((2),(1)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDAQ      ((2),(2)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDAQ      ((2),(5)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((2))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),X          ((2)))                     
*FOX  Y(2)=ALDAQ(2,3)*PUX+ALDAQ(2,4)*PUZ+ALDAQ(2,6)*IDZ(2) ;            *FOX
      CALL DACOP(ALDAQ      ((2),(3)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDAQ      ((2),(4)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDAQ      ((2),(6)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((2))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((2)))                     
              else
*FOX  PUX=X(1) ;                                                        *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUX        )                          
*FOX  PUZ=Y(1) ;                                                        *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUZ        )                          
*FOX  SIGMDA=SIGMDA+ASDA(1,1)+ASDA(1,2)*PUX+                            *FOX
*FOX  ASDA(1,3)*PUZ+ASDA(1,4)*PUX*PUZ+ASDA(1,5)*PUX*PUX+                *FOX
*FOX  ASDA(1,6)*PUZ*PUZ ;                                               *FOX
      CALL DACOP(ASDA       ((1),(1)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ASDA       ((1),(2)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ASDA       ((1),(3)),ISCRDA(  3+IDAA))                 
      CALL DACOP(ASDA       ((1),(4)),ISCRDA(  4+IDAA))                 
      CALL DACOP(ASDA       ((1),(5)),ISCRDA(  5+IDAA))                 
      CALL DACOP(ASDA       ((1),(6)),ISCRDA(  6+IDAA))                 
      CALL DAMUL(ISCRDA(  2+IDAA),PUX        ,ISCRDA(  7+IDAA))         
      CALL DAMUL(ISCRDA(  3+IDAA),PUZ        ,ISCRDA(  8+IDAA))         
      CALL DAMUL(ISCRDA(  4+IDAA),PUX        ,ISCRDA(  9+IDAA))         
      CALL DAMUL(ISCRDA(  9+IDAA),PUZ        ,ISCRDA( 10+IDAA))         
      CALL DAMUL(ISCRDA(  5+IDAA),PUX        ,ISCRDA( 11+IDAA))         
      CALL DAMUL(ISCRDA( 11+IDAA),PUX        ,ISCRDA( 12+IDAA))         
      CALL DAMUL(ISCRDA(  6+IDAA),PUZ        ,ISCRDA( 13+IDAA))         
      CALL DAMUL(ISCRDA( 13+IDAA),PUZ        ,ISCRDA( 14+IDAA))         
      CALL DAADD(SIGMDA     ,ISCRDA(  1+IDAA),ISCRDA( 15+IDAA))         
      CALL DAADD(ISCRDA( 15+IDAA),ISCRDA(  7+IDAA),ISCRDA( 16+IDAA))    
      CALL DAADD(ISCRDA( 16+IDAA),ISCRDA(  8+IDAA),ISCRDA( 17+IDAA))    
      CALL DAADD(ISCRDA( 17+IDAA),ISCRDA( 10+IDAA),ISCRDA( 18+IDAA))    
      CALL DAADD(ISCRDA( 18+IDAA),ISCRDA( 12+IDAA),ISCRDA( 19+IDAA))    
      CALL DAADD(ISCRDA( 19+IDAA),ISCRDA( 14+IDAA),ISCRDA( 20+IDAA))    
      CALL DACOP(ISCRDA( 20+IDAA),SIGMDA     )                          
*FOX  X(1)=ALDA(1,1)*PUX+ALDA(1,2)*PUZ+ALDA(1,5)*IDZ(1) ;               *FOX
      CALL DACOP(ALDA       ((1),(1)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDA       ((1),(2)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDA       ((1),(5)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((1))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),X          ((1)))                     
*FOX  Y(1)=ALDA(1,3)*PUX+ALDA(1,4)*PUZ+ALDA(1,6)*IDZ(1) ;               *FOX
      CALL DACOP(ALDA       ((1),(3)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDA       ((1),(4)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDA       ((1),(6)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((1))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((1)))                     
*FOX  PUX=X(2) ;                                                        *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUX        )                          
*FOX  PUZ=Y(2) ;                                                        *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUZ        )                          
*FOX  SIGMDA=SIGMDA+ASDA(2,1)+ASDA(2,2)*PUX+                            *FOX
*FOX  ASDA(2,3)*PUZ+ASDA(2,4)*PUX*PUZ+ASDA(2,5)*PUX*PUX+                *FOX
*FOX  ASDA(2,6)*PUZ*PUZ ;                                               *FOX
      CALL DACOP(ASDA       ((2),(1)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ASDA       ((2),(2)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ASDA       ((2),(3)),ISCRDA(  3+IDAA))                 
      CALL DACOP(ASDA       ((2),(4)),ISCRDA(  4+IDAA))                 
      CALL DACOP(ASDA       ((2),(5)),ISCRDA(  5+IDAA))                 
      CALL DACOP(ASDA       ((2),(6)),ISCRDA(  6+IDAA))                 
      CALL DAMUL(ISCRDA(  2+IDAA),PUX        ,ISCRDA(  7+IDAA))         
      CALL DAMUL(ISCRDA(  3+IDAA),PUZ        ,ISCRDA(  8+IDAA))         
      CALL DAMUL(ISCRDA(  4+IDAA),PUX        ,ISCRDA(  9+IDAA))         
      CALL DAMUL(ISCRDA(  9+IDAA),PUZ        ,ISCRDA( 10+IDAA))         
      CALL DAMUL(ISCRDA(  5+IDAA),PUX        ,ISCRDA( 11+IDAA))         
      CALL DAMUL(ISCRDA( 11+IDAA),PUX        ,ISCRDA( 12+IDAA))         
      CALL DAMUL(ISCRDA(  6+IDAA),PUZ        ,ISCRDA( 13+IDAA))         
      CALL DAMUL(ISCRDA( 13+IDAA),PUZ        ,ISCRDA( 14+IDAA))         
      CALL DAADD(SIGMDA     ,ISCRDA(  1+IDAA),ISCRDA( 15+IDAA))         
      CALL DAADD(ISCRDA( 15+IDAA),ISCRDA(  7+IDAA),ISCRDA( 16+IDAA))    
      CALL DAADD(ISCRDA( 16+IDAA),ISCRDA(  8+IDAA),ISCRDA( 17+IDAA))    
      CALL DAADD(ISCRDA( 17+IDAA),ISCRDA( 10+IDAA),ISCRDA( 18+IDAA))    
      CALL DAADD(ISCRDA( 18+IDAA),ISCRDA( 12+IDAA),ISCRDA( 19+IDAA))    
      CALL DAADD(ISCRDA( 19+IDAA),ISCRDA( 14+IDAA),ISCRDA( 20+IDAA))    
      CALL DACOP(ISCRDA( 20+IDAA),SIGMDA     )                          
*FOX  X(2)=ALDA(2,1)*PUX+ALDA(2,2)*PUZ+ALDA(2,5)*IDZ(2) ;               *FOX
      CALL DACOP(ALDA       ((2),(1)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDA       ((2),(2)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDA       ((2),(5)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((2))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),X          ((2)))                     
*FOX  Y(2)=ALDA(2,3)*PUX+ALDA(2,4)*PUZ+ALDA(2,6)*IDZ(2) ;               *FOX
      CALL DACOP(ALDA       ((2),(3)),ISCRDA(  1+IDAA))                 
      CALL DACOP(ALDA       ((2),(4)),ISCRDA(  2+IDAA))                 
      CALL DACOP(ALDA       ((2),(6)),ISCRDA(  3+IDAA))                 
      ISCRRI(  4+IDAA) = IDZ        ((2))                               
      CALL DAMUL(ISCRDA(  1+IDAA),PUX        ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),PUZ        ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ISCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((2)))                     
              endif
            else
              if(iexact.eq.1) then
!-----------------------------------------------------------------------
!  EXACT DRIFT
!-----------------------------------------------------------------------
*FOX  X(1)=X(1)*C1M3 ;                                                  *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1M3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),X          ((1)))                     
*FOX  X(2)=X(2)*C1M3 ;                                                  *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1M3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),X          ((2)))                     
*FOX  Y(1)=Y(1)*C1M3 ;                                                  *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1M3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)*C1M3 ;                                                  *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1M3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),Y          ((2)))                     
*FOX  SIGMDA=SIGMDA*C1M3 ;                                              *FOX
      CALL DACMU(SIGMDA     ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),SIGMDA     )                          
*FOX  PZ=SQRT(ONE-Y(1)*Y(1)-Y(2)*Y(2)) ;                                *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(Y          ((1)),ISCRDA(  2+IDAA))                     
      CALL DACOP(Y          ((2)),ISCRDA(  3+IDAA))                     
      CALL DACOP(Y          ((2)),ISCRDA(  4+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  6+IDAA))    
      CALL DASUC(ISCRDA(  5+IDAA),ONE*ONE        ,ISCRDA(  7+IDAA))     
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAFUN('SQRT',ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))              
      CALL DACOP(ISCRDA(  9+IDAA),PZ         )                          
*FOX  X(1)=X(1)+EL(JX)*(Y(1)/PZ) ;                                      *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DADIV(ISCRDA(  1+IDAA),PZ         ,ISCRDA(  2+IDAA))         
      CALL DACOP(X          ((1)),ISCRDA(  3+IDAA))                     
      RSCRRI(  4+IDAA) = EL         (JX         )                       
      CALL DACMU(ISCRDA(  2+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACOP(ISCRDA(  6+IDAA),X          ((1)))                     
*FOX  X(2)=X(2)+EL(JX)*(Y(2)/PZ) ;                                      *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DADIV(ISCRDA(  1+IDAA),PZ         ,ISCRDA(  2+IDAA))         
      CALL DACOP(X          ((2)),ISCRDA(  3+IDAA))                     
      RSCRRI(  4+IDAA) = EL         (JX         )                       
      CALL DACMU(ISCRDA(  2+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACOP(ISCRDA(  6+IDAA),X          ((2)))                     
*FOX  SIGMDA=SIGMDA+(ONE-(RV/PZ))*EL(JX) ;                              *FOX
      CALL DADIV(RV         ,PZ         ,ISCRDA(  1+IDAA))              
      CALL DASUC(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  2+IDAA))     
      RSCRRI(  3+IDAA) = EL         (JX         )                       
      CALL DACMU(ISCRDA(  2+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      CALL DAADD(SIGMDA     ,ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))         
      CALL DACOP(ISCRDA(  5+IDAA),SIGMDA     )                          
*FOX  X(1)=X(1)*C1E3 ;                                                  *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1E3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),X          ((1)))                     
*FOX  X(2)=X(2)*C1E3 ;                                                  *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1E3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),X          ((2)))                     
*FOX  Y(1)=Y(1)*C1E3 ;                                                  *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1E3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)*C1E3 ;                                                  *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1E3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),Y          ((2)))                     
*FOX  SIGMDA=SIGMDA*C1E3 ;                                              *FOX
      CALL DACMU(SIGMDA     ,ONE*C1E3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),SIGMDA     )                          
!-----------------------------------------------------------------------
              else
! Regular drift
!            else !moved outside of dalin6 /Mattias
*FOX  X(1)=X(1)+EL(JX)*Y(1) ;                                           *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = EL         (JX         )                       
      CALL DACOP(Y          ((1)),ISCRDA(  3+IDAA))                     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(ISCRDA(  5+IDAA),X          ((1)))                     
*FOX  X(2)=X(2)+EL(JX)*Y(2) ;                                           *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = EL         (JX         )                       
      CALL DACOP(Y          ((2)),ISCRDA(  3+IDAA))                     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(ISCRDA(  5+IDAA),X          ((2)))                     
!      
*FOX  SIGMDA=SIGMDA+                                                    *FOX
*FOX  EL(JX)*(C1E3-RV*(C1E3+(Y(1)*Y(1)+Y(2)*Y(2))*C5M4)) ;              *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(Y          ((1)),ISCRDA(  2+IDAA))                     
      CALL DACOP(Y          ((2)),ISCRDA(  3+IDAA))                     
      CALL DACOP(Y          ((2)),ISCRDA(  4+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  6+IDAA))    
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  7+IDAA))    
      CALL DACMU(ISCRDA(  7+IDAA),ONE*C5M4       ,ISCRDA(  8+IDAA))     
      CALL DACAD(ISCRDA(  8+IDAA),ONE*C1E3       ,ISCRDA(  9+IDAA))     
      CALL DAMUL(RV         ,ISCRDA(  9+IDAA),ISCRDA( 10+IDAA))         
      CALL DASUC(ISCRDA( 10+IDAA),ONE*C1E3       ,ISCRDA( 11+IDAA))     
      RSCRRI( 12+IDAA) = EL         (JX         )                       
      CALL DACMU(ISCRDA( 11+IDAA),ONE*RSCRRI( 12+IDAA),ISCRDA( 13+IDAA))
     *                                                                  
      CALL DAADD(SIGMDA     ,ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))         
      CALL DACOP(ISCRDA( 14+IDAA),SIGMDA     )                          
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
*FOX  YP(1)=Y(1)*(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((1)))                     
*FOX  YP(2)=Y(2)*(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((2)))                     
*FOX  DPDA1=DPDA*C1E3 ;                                                 *FOX
      CALL DACMU(DPDA       ,ONE*C1E3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA1      )                          
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
*FOX  EJF0=EJF1 ;                                                       *FOX
      CALL DACOP(EJF1       ,EJF0       )                               
          ixcav=ix
          if(abs(dppoff).gt.pieni) then
            sigmdac=dare(sigmda)
            sigmoff(i)=sigmdac
*FOX  SIGMDA=SIGMDA-SIGMDAC ;                                           *FOX
      CALL DACSU(SIGMDA     ,ONE*SIGMDAC    ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),SIGMDA     )                          
          endif
          call synoda
*FOX  DPDA=DPDA1*C1M3 ;                                                 *FOX
      CALL DACMU(DPDA1      ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA       )                          
*FOX  Y(1)=EJF0/EJF1*Y(1) ;                                             *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((1)))                     
*FOX  Y(2)=EJF0/EJF1*Y(2) ;                                             *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((2)))                     
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
*FOX  XX(1)=X(1) ;                                                      *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),XX         ((1)))                     
*FOX  XX(2)=X(2) ;                                                      *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),XX         ((2)))                     
*FOX  YY(1)=Y(1) ;                                                      *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),YY         ((1)))                     
*FOX  YY(2)=Y(2) ;                                                      *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),YY         ((2)))                     
          call wireda
*FOX  X(1)=XX(1) ;                                                      *FOX
      CALL DACOP(XX         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),X          ((1)))                     
*FOX  X(2)=XX(2) ;                                                      *FOX
      CALL DACOP(XX         ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),X          ((2)))                     
*FOX  Y(1)=YY(1) ;                                                      *FOX
      CALL DACOP(YY         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),Y          ((1)))                     
*FOX  Y(2)=YY(2) ;                                                      *FOX
      CALL DACOP(YY         ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),Y          ((2)))                     
          goto 440
        endif
        if(ilinc.eq.2.and.kzz.eq.20) then
          if(nbeam.ge.1) then
          ibb=ibb+1
          if(ibb.gt.nbb) call prror(102)
          imbb(i)=ibb
*FOX  YP(1)=Y(1)*(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((1)))                     
*FOX  YP(2)=Y(2)*(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((2)))                     
*FOX  DPDA1=DPDA*C1E3 ;                                                 *FOX
      CALL DACMU(DPDA       ,ONE*C1E3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA1      )                          
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
*FOX  CRKVEBF=X(1) ;                                                    *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CRKVEBF    )                          
*FOX  CIKVEBF=X(2) ;                                                    *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CIKVEBF    )                          
!hr03       startco=dare(x(1))-clobeam(1,imbb(i))+ed(ix)
            startco=(dare(x(1))-clobeam(1,imbb(i)))+ed(ix)               !hr03
            call dapok(crkvebf,jj,startco)
!hr03       startco=dare(x(2))-clobeam(2,imbb(i))+ek(ix)
            startco=(dare(x(2))-clobeam(2,imbb(i)))+ek(ix)
            call dapok(cikvebf,jj,startco)
            if(ibbc.eq.1) then
*FOX  CCCC=CRKVEBF ;                                                    *FOX
      CALL DACOP(CRKVEBF    ,CCCC       )                               
*FOX  CRKVEBF=CCCC*BBCU(IMBB(I),11)+CIKVEBF*BBCU(IMBB(I),12) ;          *FOX
      ISCRRI(  1+IDAA) = IMBB       (I          )                       
      ISCRRI(  2+IDAA) = IMBB       (I          )                       
      RSCRRI(  3+IDAA) = BBCU       (ISCRRI(  1+IDAA),(11))             
      RSCRRI(  4+IDAA) = BBCU       (ISCRRI(  2+IDAA),(12))             
      CALL DACMU(CCCC       ,ONE*RSCRRI(  3+IDAA),ISCRDA(  5+IDAA))     
      CALL DACMU(CIKVEBF    ,ONE*RSCRRI(  4+IDAA),ISCRDA(  6+IDAA))     
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  7+IDAA))    
      CALL DACOP(ISCRDA(  7+IDAA),CRKVEBF    )                          
*FOX  CIKVEBF=-CCCC*BBCU(IMBB(I),12)+CIKVEBF*BBCU(IMBB(I),11) ;         *FOX
      ISCRRI(  1+IDAA) = IMBB       (I          )                       
      ISCRRI(  2+IDAA) = IMBB       (I          )                       
      RSCRRI(  3+IDAA) = BBCU       (ISCRRI(  1+IDAA),(12))             
      RSCRRI(  4+IDAA) = BBCU       (ISCRRI(  2+IDAA),(11))             
      CALL DACMU(CCCC       ,ONE*(-ONE       ),ISCRDA(  5+IDAA))        
      CALL DACMU(ISCRDA(  5+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))
     *                                                                  
      CALL DACMU(CIKVEBF    ,ONE*RSCRRI(  4+IDAA),ISCRDA(  7+IDAA))     
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),CIKVEBF    )                          
            endif
*FOX  RHO2BF=CRKVEBF*CRKVEBF+CIKVEBF*CIKVEBF ;                          *FOX
      CALL DAMUL(CRKVEBF    ,CRKVEBF    ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVEBF    ,CIKVEBF    ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),RHO2BF     )                          
              if(abs(dare(rho2bf)).gt.pieni) then
      if(abs(sigman(1,imbb(i))).lt.pieni) call prror(88)
*FOX  TKBF=RHO2BF/(TWO*SIGMAN(1,IMBB(I))*SIGMAN(1,IMBB(I))) ;           *FOX
      ISCRRI(  1+IDAA) = IMBB       (I          )                       
      ISCRRI(  2+IDAA) = IMBB       (I          )                       
      RSCRRI(  3+IDAA) = SIGMAN     ((1),ISCRRI(  1+IDAA))              
      RSCRRI(  4+IDAA) = SIGMAN     ((1),ISCRRI(  2+IDAA))              
      RSCRRI(  5+IDAA) = TWO         * RSCRRI(  3+IDAA)                 
      RSCRRI(  6+IDAA) = RSCRRI(  5+IDAA) * RSCRRI(  4+IDAA)            
      CALL DACDI(RHO2BF     ,ONE*RSCRRI(  6+IDAA),ISCRDA(  7+IDAA))     
      CALL DACOP(ISCRDA(  7+IDAA),TKBF       )                          
      if(ibbc.eq.0) then
*FOX   Y(1)=Y(1)+(CRAD*CRKVEBF/RHO2BF*                                  *FOX
*FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF4)/(ONE+DPDA) ;               *FOX
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAFUN('EXP ',ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))              
      CALL DASUC(ISCRDA(  2+IDAA),ONE*ONE        ,ISCRDA(  3+IDAA))     
      RSCRRI(  4+IDAA) = PTNFAC     (IX         )                       
      CALL DACMU(CRKVEBF    ,ONE*CRAD       ,ISCRDA(  5+IDAA))          
      CALL DADIV(ISCRDA(  5+IDAA),RHO2BF     ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  6+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA(  7+IDAA),ISCRDA(  3+IDAA),ISCRDA(  8+IDAA))    
      CALL DACSU(ISCRDA(  8+IDAA),ONE*BEAMOFF4   ,ISCRDA(  9+IDAA))     
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 10+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA( 11+IDAA))                     
      CALL DADIV(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 12+IDAA))    
      CALL DAADD(ISCRDA( 11+IDAA),ISCRDA( 12+IDAA),ISCRDA( 13+IDAA))    
      CALL DACOP(ISCRDA( 13+IDAA),Y          ((1)))                     
*FOX   Y(2)=Y(2)+(CRAD*CIKVEBF/RHO2BF*                                  *FOX
*FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF5)/(ONE+DPDA) ;               *FOX
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAFUN('EXP ',ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))              
      CALL DASUC(ISCRDA(  2+IDAA),ONE*ONE        ,ISCRDA(  3+IDAA))     
      RSCRRI(  4+IDAA) = PTNFAC     (IX         )                       
      CALL DACMU(CIKVEBF    ,ONE*CRAD       ,ISCRDA(  5+IDAA))          
      CALL DADIV(ISCRDA(  5+IDAA),RHO2BF     ,ISCRDA(  6+IDAA))         
      CALL DACMU(ISCRDA(  6+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA(  7+IDAA),ISCRDA(  3+IDAA),ISCRDA(  8+IDAA))    
      CALL DACSU(ISCRDA(  8+IDAA),ONE*BEAMOFF5   ,ISCRDA(  9+IDAA))     
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 10+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA( 11+IDAA))                     
      CALL DADIV(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 12+IDAA))    
      CALL DAADD(ISCRDA( 11+IDAA),ISCRDA( 12+IDAA),ISCRDA( 13+IDAA))    
      CALL DACOP(ISCRDA( 13+IDAA),Y          ((2)))                     
      else
*FOX   CCCC=(CRAD*CRKVEBF/RHO2BF*                                       *FOX
*FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF4)*BBCU(IMBB(I),11)-          *FOX
*FOX   (CRAD*CIKVEBF/RHO2BF*                                            *FOX
*FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF5)*BBCU(IMBB(I),12) ;         *FOX
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  2+IDAA))        
      CALL DAFUN('EXP ',ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))              
      CALL DAFUN('EXP ',ISCRDA(  2+IDAA),ISCRDA(  4+IDAA))              
      CALL DASUC(ISCRDA(  3+IDAA),ONE*ONE        ,ISCRDA(  5+IDAA))     
      CALL DASUC(ISCRDA(  4+IDAA),ONE*ONE        ,ISCRDA(  6+IDAA))     
      RSCRRI(  7+IDAA) = PTNFAC     (IX         )                       
      ISCRRI(  8+IDAA) = IMBB       (I          )                       
      RSCRRI(  9+IDAA) = PTNFAC     (IX         )                       
      ISCRRI( 10+IDAA) = IMBB       (I          )                       
      CALL DACMU(CRKVEBF    ,ONE*CRAD       ,ISCRDA( 11+IDAA))          
      CALL DADIV(ISCRDA( 11+IDAA),RHO2BF     ,ISCRDA( 12+IDAA))         
      CALL DACMU(ISCRDA( 12+IDAA),ONE*RSCRRI(  7+IDAA),ISCRDA( 13+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA( 13+IDAA),ISCRDA(  5+IDAA),ISCRDA( 14+IDAA))    
      CALL DACMU(CIKVEBF    ,ONE*CRAD       ,ISCRDA( 15+IDAA))          
      CALL DADIV(ISCRDA( 15+IDAA),RHO2BF     ,ISCRDA( 16+IDAA))         
      CALL DACMU(ISCRDA( 16+IDAA),ONE*RSCRRI(  9+IDAA),ISCRDA( 17+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA( 17+IDAA),ISCRDA(  6+IDAA),ISCRDA( 18+IDAA))    
      CALL DACSU(ISCRDA( 14+IDAA),ONE*BEAMOFF4   ,ISCRDA( 19+IDAA))     
      CALL DACSU(ISCRDA( 18+IDAA),ONE*BEAMOFF5   ,ISCRDA( 20+IDAA))     
      RSCRRI( 21+IDAA) = BBCU       (ISCRRI(  8+IDAA),(11))             
      RSCRRI( 22+IDAA) = BBCU       (ISCRRI( 10+IDAA),(12))             
      CALL DACMU(ISCRDA( 19+IDAA),ONE*RSCRRI( 21+IDAA),ISCRDA( 23+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 20+IDAA),ONE*RSCRRI( 22+IDAA),ISCRDA( 24+IDAA))
     *                                                                  
      CALL DASUB(ISCRDA( 23+IDAA),ISCRDA( 24+IDAA),ISCRDA( 25+IDAA))    
      CALL DACOP(ISCRDA( 25+IDAA),CCCC       )                          
*FOX   Y(1)=Y(1)+CCCC/(ONE+DPDA) ;                                      *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  2+IDAA))                     
      CALL DADIV(CCCC       ,ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))         
      CALL DAADD(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),Y          ((1)))                     
*FOX   CCCC=(CRAD*CRKVEBF/RHO2BF*                                       *FOX
*FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF4)*BBCU(IMBB(I),12)+          *FOX
*FOX   (CRAD*CIKVEBF/RHO2BF*                                            *FOX
*FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF5)*BBCU(IMBB(I),11) ;         *FOX
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  2+IDAA))        
      CALL DAFUN('EXP ',ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))              
      CALL DAFUN('EXP ',ISCRDA(  2+IDAA),ISCRDA(  4+IDAA))              
      CALL DASUC(ISCRDA(  3+IDAA),ONE*ONE        ,ISCRDA(  5+IDAA))     
      CALL DASUC(ISCRDA(  4+IDAA),ONE*ONE        ,ISCRDA(  6+IDAA))     
      RSCRRI(  7+IDAA) = PTNFAC     (IX         )                       
      ISCRRI(  8+IDAA) = IMBB       (I          )                       
      RSCRRI(  9+IDAA) = PTNFAC     (IX         )                       
      ISCRRI( 10+IDAA) = IMBB       (I          )                       
      CALL DACMU(CRKVEBF    ,ONE*CRAD       ,ISCRDA( 11+IDAA))          
      CALL DADIV(ISCRDA( 11+IDAA),RHO2BF     ,ISCRDA( 12+IDAA))         
      CALL DACMU(ISCRDA( 12+IDAA),ONE*RSCRRI(  7+IDAA),ISCRDA( 13+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA( 13+IDAA),ISCRDA(  5+IDAA),ISCRDA( 14+IDAA))    
      CALL DACMU(CIKVEBF    ,ONE*CRAD       ,ISCRDA( 15+IDAA))          
      CALL DADIV(ISCRDA( 15+IDAA),RHO2BF     ,ISCRDA( 16+IDAA))         
      CALL DACMU(ISCRDA( 16+IDAA),ONE*RSCRRI(  9+IDAA),ISCRDA( 17+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA( 17+IDAA),ISCRDA(  6+IDAA),ISCRDA( 18+IDAA))    
      CALL DACSU(ISCRDA( 14+IDAA),ONE*BEAMOFF4   ,ISCRDA( 19+IDAA))     
      CALL DACSU(ISCRDA( 18+IDAA),ONE*BEAMOFF5   ,ISCRDA( 20+IDAA))     
      RSCRRI( 21+IDAA) = BBCU       (ISCRRI(  8+IDAA),(12))             
      RSCRRI( 22+IDAA) = BBCU       (ISCRRI( 10+IDAA),(11))             
      CALL DACMU(ISCRDA( 19+IDAA),ONE*RSCRRI( 21+IDAA),ISCRDA( 23+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 20+IDAA),ONE*RSCRRI( 22+IDAA),ISCRDA( 24+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA( 23+IDAA),ISCRDA( 24+IDAA),ISCRDA( 25+IDAA))    
      CALL DACOP(ISCRDA( 25+IDAA),CCCC       )                          
*FOX   Y(2)=Y(2)+CCCC/(ONE+DPDA) ;                                      *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  2+IDAA))                     
      CALL DADIV(CCCC       ,ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))         
      CALL DAADD(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),Y          ((2)))                     
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
*FOX  CRKVEBF=X(1) ;                                                    *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CRKVEBF    )                          
*FOX  CIKVEBF=X(2) ;                                                    *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CIKVEBF    )                          
!hr03       startco=dare(x(1))-clobeam(1,imbb(i))+ed(ix)
            startco=(dare(x(1))-clobeam(1,imbb(i)))+ed(ix)               !hr03
            call dapok(crkvebf,jj,startco)
!hr03       startco=dare(x(2))-clobeam(2,imbb(i))+ek(ix)
            startco=(dare(x(2))-clobeam(2,imbb(i)))+ek(ix)
            call dapok(cikvebf,jj,startco)
            if(ibbc.eq.1) then
*FOX  CCCC=CRKVEBF ;                                                    *FOX
      CALL DACOP(CRKVEBF    ,CCCC       )                               
*FOX  CRKVEBF=CCCC*BBCU(IMBB(I),11)+CIKVEBF*BBCU(IMBB(I),12) ;          *FOX
      ISCRRI(  1+IDAA) = IMBB       (I          )                       
      ISCRRI(  2+IDAA) = IMBB       (I          )                       
      RSCRRI(  3+IDAA) = BBCU       (ISCRRI(  1+IDAA),(11))             
      RSCRRI(  4+IDAA) = BBCU       (ISCRRI(  2+IDAA),(12))             
      CALL DACMU(CCCC       ,ONE*RSCRRI(  3+IDAA),ISCRDA(  5+IDAA))     
      CALL DACMU(CIKVEBF    ,ONE*RSCRRI(  4+IDAA),ISCRDA(  6+IDAA))     
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  7+IDAA))    
      CALL DACOP(ISCRDA(  7+IDAA),CRKVEBF    )                          
*FOX  CIKVEBF=-CCCC*BBCU(IMBB(I),12)+CIKVEBF*BBCU(IMBB(I),11) ;         *FOX
      ISCRRI(  1+IDAA) = IMBB       (I          )                       
      ISCRRI(  2+IDAA) = IMBB       (I          )                       
      RSCRRI(  3+IDAA) = BBCU       (ISCRRI(  1+IDAA),(12))             
      RSCRRI(  4+IDAA) = BBCU       (ISCRRI(  2+IDAA),(11))             
      CALL DACMU(CCCC       ,ONE*(-ONE       ),ISCRDA(  5+IDAA))        
      CALL DACMU(ISCRDA(  5+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))
     *                                                                  
      CALL DACMU(CIKVEBF    ,ONE*RSCRRI(  4+IDAA),ISCRDA(  7+IDAA))     
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),CIKVEBF    )                          
            endif
*FOX  XRBF=CRKVEBF/RBF ;                                                *FOX
      CALL DACDI(CRKVEBF    ,ONE*RBF        ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),XRBF       )                          
      if(dare(xrbf).lt.zero) then
*FOX  XRBF=-XRBF ;                                                      *FOX
      CALL DACMU(XRBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),XRBF       )                          
      endif
*FOX  ZRBF=CIKVEBF/RBF ;                                                *FOX
      CALL DACDI(CIKVEBF    ,ONE*RBF        ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),ZRBF       )                          
      if(dare(zrbf).lt.zero) then
*FOX  ZRBF=-ZRBF ;                                                      *FOX
      CALL DACMU(ZRBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),ZRBF       )                          
      endif
            call errff(xrbf,zrbf,crxbf,crzbf)
      if(abs(sigman(1,imbb(i))).lt.pieni.or.                            &
     &abs(sigman(2,imbb(i))).lt.pieni) call prror(88)
*FOX  TKBF=(CRKVEBF*CRKVEBF/(SIGMAN(1,IMBB(I))*SIGMAN(1,IMBB(I)))+      *FOX
*FOX  CIKVEBF*CIKVEBF/(SIGMAN(2,IMBB(I))*SIGMAN(2,IMBB(I))))*HALF ;     *FOX
      ISCRRI(  1+IDAA) = IMBB       (I          )                       
      ISCRRI(  2+IDAA) = IMBB       (I          )                       
      ISCRRI(  3+IDAA) = IMBB       (I          )                       
      ISCRRI(  4+IDAA) = IMBB       (I          )                       
      RSCRRI(  5+IDAA) = SIGMAN     ((1),ISCRRI(  1+IDAA))              
      RSCRRI(  6+IDAA) = SIGMAN     ((1),ISCRRI(  2+IDAA))              
      RSCRRI(  7+IDAA) = SIGMAN     ((2),ISCRRI(  3+IDAA))              
      RSCRRI(  8+IDAA) = SIGMAN     ((2),ISCRRI(  4+IDAA))              
      RSCRRI(  9+IDAA) = RSCRRI(  5+IDAA) * RSCRRI(  6+IDAA)            
      RSCRRI( 10+IDAA) = RSCRRI(  7+IDAA) * RSCRRI(  8+IDAA)            
      CALL DAMUL(CRKVEBF    ,CRKVEBF    ,ISCRDA( 11+IDAA))              
      CALL DACDI(ISCRDA( 11+IDAA),ONE*RSCRRI(  9+IDAA),ISCRDA( 12+IDAA))
     *                                                                  
      CALL DAMUL(CIKVEBF    ,CIKVEBF    ,ISCRDA( 13+IDAA))              
      CALL DACDI(ISCRDA( 13+IDAA),ONE*RSCRRI( 10+IDAA),ISCRDA( 14+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA( 12+IDAA),ISCRDA( 14+IDAA),ISCRDA( 15+IDAA))    
      CALL DACMU(ISCRDA( 15+IDAA),ONE*HALF       ,ISCRDA( 16+IDAA))     
      CALL DACOP(ISCRDA( 16+IDAA),TKBF       )                          
*FOX  XBBF=SIGMAN(2,IMBB(I))/SIGMAN(1,IMBB(I))*XRBF ;                   *FOX
      ISCRRI(  1+IDAA) = IMBB       (I          )                       
      ISCRRI(  2+IDAA) = IMBB       (I          )                       
      RSCRRI(  3+IDAA) = SIGMAN     ((2),ISCRRI(  1+IDAA))              
      RSCRRI(  4+IDAA) = SIGMAN     ((1),ISCRRI(  2+IDAA))              
      RSCRRI(  5+IDAA) = RSCRRI(  3+IDAA) / RSCRRI(  4+IDAA)            
      CALL DACMU(XRBF       ,ONE*RSCRRI(  5+IDAA),ISCRDA(  6+IDAA))     
      CALL DACOP(ISCRDA(  6+IDAA),XBBF       )                          
*FOX  ZBBF=SIGMAN(1,IMBB(I))/SIGMAN(2,IMBB(I))*ZRBF ;                   *FOX
      ISCRRI(  1+IDAA) = IMBB       (I          )                       
      ISCRRI(  2+IDAA) = IMBB       (I          )                       
      RSCRRI(  3+IDAA) = SIGMAN     ((1),ISCRRI(  1+IDAA))              
      RSCRRI(  4+IDAA) = SIGMAN     ((2),ISCRRI(  2+IDAA))              
      RSCRRI(  5+IDAA) = RSCRRI(  3+IDAA) / RSCRRI(  4+IDAA)            
      CALL DACMU(ZRBF       ,ONE*RSCRRI(  5+IDAA),ISCRDA(  6+IDAA))     
      CALL DACOP(ISCRDA(  6+IDAA),ZBBF       )                          
            call errff(xbbf,zbbf,cbxbf,cbzbf)
      scrkveb=sign(one,dare(crkvebf))
      scikveb=sign(one,dare(cikvebf))
      if(ibbc.eq.0) then
*FOX  Y(1)=Y(1)+(RKBF*(CRZBF-EXP(-TKBF)*                                *FOX
*FOX  CBZBF)*SCRKVEB-BEAMOFF4)/(ONE+DPDA) ;                             *FOX
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAFUN('EXP ',ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),CBZBF      ,ISCRDA(  3+IDAA))         
      CALL DASUB(CRZBF      ,ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))         
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RKBF       ,ISCRDA(  5+IDAA))     
      CALL DACMU(ISCRDA(  5+IDAA),ONE*SCRKVEB    ,ISCRDA(  6+IDAA))     
      CALL DACSU(ISCRDA(  6+IDAA),ONE*BEAMOFF4   ,ISCRDA(  7+IDAA))     
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  8+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  9+IDAA))                     
      CALL DADIV(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DAADD(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(ISCRDA( 11+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+(RKBF*(CRXBF-EXP(-TKBF)*                                *FOX
*FOX  CBXBF)*SCIKVEB-BEAMOFF5)/(ONE+DPDA) ;                             *FOX
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAFUN('EXP ',ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),CBXBF      ,ISCRDA(  3+IDAA))         
      CALL DASUB(CRXBF      ,ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))         
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RKBF       ,ISCRDA(  5+IDAA))     
      CALL DACMU(ISCRDA(  5+IDAA),ONE*SCIKVEB    ,ISCRDA(  6+IDAA))     
      CALL DACSU(ISCRDA(  6+IDAA),ONE*BEAMOFF5   ,ISCRDA(  7+IDAA))     
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  8+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  9+IDAA))                     
      CALL DADIV(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DAADD(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(ISCRDA( 11+IDAA),Y          ((2)))                     
      else
*FOX  CCCC=(RKBF*(CRZBF-EXP(-TKBF)*                                     *FOX
*FOX  CBZBF)*SCRKVEB-BEAMOFF4)*BBCU(IMBB(I),11)-                        *FOX
*FOX  (RKBF*(CRXBF-EXP(-TKBF)*                                          *FOX
*FOX  CBXBF)*SCIKVEB-BEAMOFF5)*BBCU(IMBB(I),12) ;                       *FOX
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  2+IDAA))        
      CALL DAFUN('EXP ',ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))              
      CALL DAFUN('EXP ',ISCRDA(  2+IDAA),ISCRDA(  4+IDAA))              
      CALL DAMUL(ISCRDA(  3+IDAA),CBZBF      ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  4+IDAA),CBXBF      ,ISCRDA(  6+IDAA))         
      CALL DASUB(CRZBF      ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DASUB(CRXBF      ,ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))         
      ISCRRI(  9+IDAA) = IMBB       (I          )                       
      ISCRRI( 10+IDAA) = IMBB       (I          )                       
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RKBF       ,ISCRDA( 11+IDAA))     
      CALL DACMU(ISCRDA( 11+IDAA),ONE*SCRKVEB    ,ISCRDA( 12+IDAA))     
      CALL DACMU(ISCRDA(  8+IDAA),ONE*RKBF       ,ISCRDA( 13+IDAA))     
      CALL DACMU(ISCRDA( 13+IDAA),ONE*SCIKVEB    ,ISCRDA( 14+IDAA))     
      CALL DACSU(ISCRDA( 12+IDAA),ONE*BEAMOFF4   ,ISCRDA( 15+IDAA))     
      CALL DACSU(ISCRDA( 14+IDAA),ONE*BEAMOFF5   ,ISCRDA( 16+IDAA))     
      RSCRRI( 17+IDAA) = BBCU       (ISCRRI(  9+IDAA),(11))             
      RSCRRI( 18+IDAA) = BBCU       (ISCRRI( 10+IDAA),(12))             
      CALL DACMU(ISCRDA( 15+IDAA),ONE*RSCRRI( 17+IDAA),ISCRDA( 19+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 16+IDAA),ONE*RSCRRI( 18+IDAA),ISCRDA( 20+IDAA))
     *                                                                  
      CALL DASUB(ISCRDA( 19+IDAA),ISCRDA( 20+IDAA),ISCRDA( 21+IDAA))    
      CALL DACOP(ISCRDA( 21+IDAA),CCCC       )                          
*FOX   Y(1)=Y(1)+CCCC/(ONE+DPDA) ;                                      *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  2+IDAA))                     
      CALL DADIV(CCCC       ,ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))         
      CALL DAADD(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),Y          ((1)))                     
*FOX  CCCC=(RKBF*(CRZBF-EXP(-TKBF)*                                     *FOX
*FOX  CBZBF)*SCRKVEB-BEAMOFF4)*BBCU(IMBB(I),12)+                        *FOX
*FOX  (RKBF*(CRXBF-EXP(-TKBF)*                                          *FOX
*FOX  CBXBF)*SCIKVEB-BEAMOFF5)*BBCU(IMBB(I),11) ;                       *FOX
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  2+IDAA))        
      CALL DAFUN('EXP ',ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))              
      CALL DAFUN('EXP ',ISCRDA(  2+IDAA),ISCRDA(  4+IDAA))              
      CALL DAMUL(ISCRDA(  3+IDAA),CBZBF      ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  4+IDAA),CBXBF      ,ISCRDA(  6+IDAA))         
      CALL DASUB(CRZBF      ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DASUB(CRXBF      ,ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))         
      ISCRRI(  9+IDAA) = IMBB       (I          )                       
      ISCRRI( 10+IDAA) = IMBB       (I          )                       
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RKBF       ,ISCRDA( 11+IDAA))     
      CALL DACMU(ISCRDA( 11+IDAA),ONE*SCRKVEB    ,ISCRDA( 12+IDAA))     
      CALL DACMU(ISCRDA(  8+IDAA),ONE*RKBF       ,ISCRDA( 13+IDAA))     
      CALL DACMU(ISCRDA( 13+IDAA),ONE*SCIKVEB    ,ISCRDA( 14+IDAA))     
      CALL DACSU(ISCRDA( 12+IDAA),ONE*BEAMOFF4   ,ISCRDA( 15+IDAA))     
      CALL DACSU(ISCRDA( 14+IDAA),ONE*BEAMOFF5   ,ISCRDA( 16+IDAA))     
      RSCRRI( 17+IDAA) = BBCU       (ISCRRI(  9+IDAA),(12))             
      RSCRRI( 18+IDAA) = BBCU       (ISCRRI( 10+IDAA),(11))             
      CALL DACMU(ISCRDA( 15+IDAA),ONE*RSCRRI( 17+IDAA),ISCRDA( 19+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 16+IDAA),ONE*RSCRRI( 18+IDAA),ISCRDA( 20+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA( 19+IDAA),ISCRDA( 20+IDAA),ISCRDA( 21+IDAA))    
      CALL DACOP(ISCRDA( 21+IDAA),CCCC       )                          
*FOX   Y(2)=Y(2)+CCCC/(ONE+DPDA) ;                                      *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  2+IDAA))                     
      CALL DADIV(CCCC       ,ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))         
      CALL DAADD(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),Y          ((2)))                     
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
*FOX  CRKVEBF=X(1) ;                                                    *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CRKVEBF    )                          
*FOX  CIKVEBF=X(2) ;                                                    *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CIKVEBF    )                          
!hr03       startco=dare(x(1))-clobeam(1,imbb(i))+ed(ix)
            startco=(dare(x(1))-clobeam(1,imbb(i)))+ed(ix)               !hr03
            call dapok(crkvebf,jj,startco)
!hr03       startco=dare(x(2))-clobeam(2,imbb(i))+ek(ix)
            startco=(dare(x(2))-clobeam(2,imbb(i)))+ek(ix)
            call dapok(cikvebf,jj,startco)
            if(ibbc.eq.1) then
*FOX  CCCC=CRKVEBF ;                                                    *FOX
      CALL DACOP(CRKVEBF    ,CCCC       )                               
*FOX  CRKVEBF=CCCC*BBCU(IMBB(I),11)+CIKVEBF*BBCU(IMBB(I),12) ;          *FOX
      ISCRRI(  1+IDAA) = IMBB       (I          )                       
      ISCRRI(  2+IDAA) = IMBB       (I          )                       
      RSCRRI(  3+IDAA) = BBCU       (ISCRRI(  1+IDAA),(11))             
      RSCRRI(  4+IDAA) = BBCU       (ISCRRI(  2+IDAA),(12))             
      CALL DACMU(CCCC       ,ONE*RSCRRI(  3+IDAA),ISCRDA(  5+IDAA))     
      CALL DACMU(CIKVEBF    ,ONE*RSCRRI(  4+IDAA),ISCRDA(  6+IDAA))     
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  7+IDAA))    
      CALL DACOP(ISCRDA(  7+IDAA),CRKVEBF    )                          
*FOX  CIKVEBF=-CCCC*BBCU(IMBB(I),12)+CIKVEBF*BBCU(IMBB(I),11) ;         *FOX
      ISCRRI(  1+IDAA) = IMBB       (I          )                       
      ISCRRI(  2+IDAA) = IMBB       (I          )                       
      RSCRRI(  3+IDAA) = BBCU       (ISCRRI(  1+IDAA),(12))             
      RSCRRI(  4+IDAA) = BBCU       (ISCRRI(  2+IDAA),(11))             
      CALL DACMU(CCCC       ,ONE*(-ONE       ),ISCRDA(  5+IDAA))        
      CALL DACMU(ISCRDA(  5+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))
     *                                                                  
      CALL DACMU(CIKVEBF    ,ONE*RSCRRI(  4+IDAA),ISCRDA(  7+IDAA))     
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),CIKVEBF    )                          
            endif
*FOX  XRBF=CRKVEBF/RBF ;                                                *FOX
      CALL DACDI(CRKVEBF    ,ONE*RBF        ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),XRBF       )                          
      if(dare(xrbf).lt.zero) then
*FOX  XRBF=-XRBF ;                                                      *FOX
      CALL DACMU(XRBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),XRBF       )                          
      endif
*FOX  ZRBF=CIKVEBF/RBF ;                                                *FOX
      CALL DACDI(CIKVEBF    ,ONE*RBF        ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),ZRBF       )                          
      if(dare(zrbf).lt.zero) then
*FOX  ZRBF=-ZRBF ;                                                      *FOX
      CALL DACMU(ZRBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),ZRBF       )                          
      endif
            call errff(zrbf,xrbf,crzbf,crxbf)
      if(abs(sigman(1,imbb(i))).lt.pieni.or.                            &
     &abs(sigman(2,imbb(i))).lt.pieni) call prror(88)
*FOX  TKBF=(CRKVEBF*CRKVEBF/(SIGMAN(1,IMBB(I))*SIGMAN(1,IMBB(I)))+      *FOX
*FOX  CIKVEBF*CIKVEBF/(SIGMAN(2,IMBB(I))*SIGMAN(2,IMBB(I))))*HALF ;     *FOX
      ISCRRI(  1+IDAA) = IMBB       (I          )                       
      ISCRRI(  2+IDAA) = IMBB       (I          )                       
      ISCRRI(  3+IDAA) = IMBB       (I          )                       
      ISCRRI(  4+IDAA) = IMBB       (I          )                       
      RSCRRI(  5+IDAA) = SIGMAN     ((1),ISCRRI(  1+IDAA))              
      RSCRRI(  6+IDAA) = SIGMAN     ((1),ISCRRI(  2+IDAA))              
      RSCRRI(  7+IDAA) = SIGMAN     ((2),ISCRRI(  3+IDAA))              
      RSCRRI(  8+IDAA) = SIGMAN     ((2),ISCRRI(  4+IDAA))              
      RSCRRI(  9+IDAA) = RSCRRI(  5+IDAA) * RSCRRI(  6+IDAA)            
      RSCRRI( 10+IDAA) = RSCRRI(  7+IDAA) * RSCRRI(  8+IDAA)            
      CALL DAMUL(CRKVEBF    ,CRKVEBF    ,ISCRDA( 11+IDAA))              
      CALL DACDI(ISCRDA( 11+IDAA),ONE*RSCRRI(  9+IDAA),ISCRDA( 12+IDAA))
     *                                                                  
      CALL DAMUL(CIKVEBF    ,CIKVEBF    ,ISCRDA( 13+IDAA))              
      CALL DACDI(ISCRDA( 13+IDAA),ONE*RSCRRI( 10+IDAA),ISCRDA( 14+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA( 12+IDAA),ISCRDA( 14+IDAA),ISCRDA( 15+IDAA))    
      CALL DACMU(ISCRDA( 15+IDAA),ONE*HALF       ,ISCRDA( 16+IDAA))     
      CALL DACOP(ISCRDA( 16+IDAA),TKBF       )                          
*FOX  XBBF=SIGMAN(2,IMBB(I))/SIGMAN(1,IMBB(I))*XRBF ;                   *FOX
      ISCRRI(  1+IDAA) = IMBB       (I          )                       
      ISCRRI(  2+IDAA) = IMBB       (I          )                       
      RSCRRI(  3+IDAA) = SIGMAN     ((2),ISCRRI(  1+IDAA))              
      RSCRRI(  4+IDAA) = SIGMAN     ((1),ISCRRI(  2+IDAA))              
      RSCRRI(  5+IDAA) = RSCRRI(  3+IDAA) / RSCRRI(  4+IDAA)            
      CALL DACMU(XRBF       ,ONE*RSCRRI(  5+IDAA),ISCRDA(  6+IDAA))     
      CALL DACOP(ISCRDA(  6+IDAA),XBBF       )                          
*FOX  ZBBF=SIGMAN(1,IMBB(I))/SIGMAN(2,IMBB(I))*ZRBF ;                   *FOX
      ISCRRI(  1+IDAA) = IMBB       (I          )                       
      ISCRRI(  2+IDAA) = IMBB       (I          )                       
      RSCRRI(  3+IDAA) = SIGMAN     ((1),ISCRRI(  1+IDAA))              
      RSCRRI(  4+IDAA) = SIGMAN     ((2),ISCRRI(  2+IDAA))              
      RSCRRI(  5+IDAA) = RSCRRI(  3+IDAA) / RSCRRI(  4+IDAA)            
      CALL DACMU(ZRBF       ,ONE*RSCRRI(  5+IDAA),ISCRDA(  6+IDAA))     
      CALL DACOP(ISCRDA(  6+IDAA),ZBBF       )                          
            call errff(zbbf,xbbf,cbzbf,cbxbf)
      scrkveb=sign(one,dare(crkvebf))
      scikveb=sign(one,dare(cikvebf))
      if(ibbc.eq.0) then
*FOX  Y(1)=Y(1)+(RKBF*(CRZBF-EXP(-TKBF)*                                *FOX
*FOX  CBZBF)*SCRKVEB-BEAMOFF4)/(ONE+DPDA) ;                             *FOX
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAFUN('EXP ',ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),CBZBF      ,ISCRDA(  3+IDAA))         
      CALL DASUB(CRZBF      ,ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))         
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RKBF       ,ISCRDA(  5+IDAA))     
      CALL DACMU(ISCRDA(  5+IDAA),ONE*SCRKVEB    ,ISCRDA(  6+IDAA))     
      CALL DACSU(ISCRDA(  6+IDAA),ONE*BEAMOFF4   ,ISCRDA(  7+IDAA))     
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  8+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  9+IDAA))                     
      CALL DADIV(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DAADD(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(ISCRDA( 11+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+(RKBF*(CRXBF-EXP(-TKBF)*                                *FOX
*FOX  CBXBF)*SCIKVEB-BEAMOFF5)/(ONE+DPDA) ;                             *FOX
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAFUN('EXP ',ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),CBXBF      ,ISCRDA(  3+IDAA))         
      CALL DASUB(CRXBF      ,ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))         
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RKBF       ,ISCRDA(  5+IDAA))     
      CALL DACMU(ISCRDA(  5+IDAA),ONE*SCIKVEB    ,ISCRDA(  6+IDAA))     
      CALL DACSU(ISCRDA(  6+IDAA),ONE*BEAMOFF5   ,ISCRDA(  7+IDAA))     
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  8+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  9+IDAA))                     
      CALL DADIV(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DAADD(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(ISCRDA( 11+IDAA),Y          ((2)))                     
      else
*FOX  CCCC=(RKBF*(CRZBF-EXP(-TKBF)*                                     *FOX
*FOX  CBZBF)*SCRKVEB-BEAMOFF4)*BBCU(IMBB(I),11)-                        *FOX
*FOX  (RKBF*(CRXBF-EXP(-TKBF)*                                          *FOX
*FOX  CBXBF)*SCIKVEB-BEAMOFF5)*BBCU(IMBB(I),12) ;                       *FOX
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  2+IDAA))        
      CALL DAFUN('EXP ',ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))              
      CALL DAFUN('EXP ',ISCRDA(  2+IDAA),ISCRDA(  4+IDAA))              
      CALL DAMUL(ISCRDA(  3+IDAA),CBZBF      ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  4+IDAA),CBXBF      ,ISCRDA(  6+IDAA))         
      CALL DASUB(CRZBF      ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DASUB(CRXBF      ,ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))         
      ISCRRI(  9+IDAA) = IMBB       (I          )                       
      ISCRRI( 10+IDAA) = IMBB       (I          )                       
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RKBF       ,ISCRDA( 11+IDAA))     
      CALL DACMU(ISCRDA( 11+IDAA),ONE*SCRKVEB    ,ISCRDA( 12+IDAA))     
      CALL DACMU(ISCRDA(  8+IDAA),ONE*RKBF       ,ISCRDA( 13+IDAA))     
      CALL DACMU(ISCRDA( 13+IDAA),ONE*SCIKVEB    ,ISCRDA( 14+IDAA))     
      CALL DACSU(ISCRDA( 12+IDAA),ONE*BEAMOFF4   ,ISCRDA( 15+IDAA))     
      CALL DACSU(ISCRDA( 14+IDAA),ONE*BEAMOFF5   ,ISCRDA( 16+IDAA))     
      RSCRRI( 17+IDAA) = BBCU       (ISCRRI(  9+IDAA),(11))             
      RSCRRI( 18+IDAA) = BBCU       (ISCRRI( 10+IDAA),(12))             
      CALL DACMU(ISCRDA( 15+IDAA),ONE*RSCRRI( 17+IDAA),ISCRDA( 19+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 16+IDAA),ONE*RSCRRI( 18+IDAA),ISCRDA( 20+IDAA))
     *                                                                  
      CALL DASUB(ISCRDA( 19+IDAA),ISCRDA( 20+IDAA),ISCRDA( 21+IDAA))    
      CALL DACOP(ISCRDA( 21+IDAA),CCCC       )                          
*FOX   Y(1)=Y(1)+CCCC/(ONE+DPDA) ;                                      *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  2+IDAA))                     
      CALL DADIV(CCCC       ,ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))         
      CALL DAADD(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),Y          ((1)))                     
*FOX  CCCC=(RKBF*(CRZBF-EXP(-TKBF)*                                     *FOX
*FOX  CBZBF)*SCRKVEB-BEAMOFF4)*BBCU(IMBB(I),12)+                        *FOX
*FOX  (RKBF*(CRXBF-EXP(-TKBF)*                                          *FOX
*FOX  CBXBF)*SCIKVEB-BEAMOFF5)*BBCU(IMBB(I),11) ;                       *FOX
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACMU(TKBF       ,ONE*(-ONE       ),ISCRDA(  2+IDAA))        
      CALL DAFUN('EXP ',ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))              
      CALL DAFUN('EXP ',ISCRDA(  2+IDAA),ISCRDA(  4+IDAA))              
      CALL DAMUL(ISCRDA(  3+IDAA),CBZBF      ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  4+IDAA),CBXBF      ,ISCRDA(  6+IDAA))         
      CALL DASUB(CRZBF      ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DASUB(CRXBF      ,ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))         
      ISCRRI(  9+IDAA) = IMBB       (I          )                       
      ISCRRI( 10+IDAA) = IMBB       (I          )                       
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RKBF       ,ISCRDA( 11+IDAA))     
      CALL DACMU(ISCRDA( 11+IDAA),ONE*SCRKVEB    ,ISCRDA( 12+IDAA))     
      CALL DACMU(ISCRDA(  8+IDAA),ONE*RKBF       ,ISCRDA( 13+IDAA))     
      CALL DACMU(ISCRDA( 13+IDAA),ONE*SCIKVEB    ,ISCRDA( 14+IDAA))     
      CALL DACSU(ISCRDA( 12+IDAA),ONE*BEAMOFF4   ,ISCRDA( 15+IDAA))     
      CALL DACSU(ISCRDA( 14+IDAA),ONE*BEAMOFF5   ,ISCRDA( 16+IDAA))     
      RSCRRI( 17+IDAA) = BBCU       (ISCRRI(  9+IDAA),(12))             
      RSCRRI( 18+IDAA) = BBCU       (ISCRRI( 10+IDAA),(11))             
      CALL DACMU(ISCRDA( 15+IDAA),ONE*RSCRRI( 17+IDAA),ISCRDA( 19+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 16+IDAA),ONE*RSCRRI( 18+IDAA),ISCRDA( 20+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA( 19+IDAA),ISCRDA( 20+IDAA),ISCRDA( 21+IDAA))    
      CALL DACOP(ISCRDA( 21+IDAA),CCCC       )                          
*FOX   Y(2)=Y(2)+CCCC/(ONE+DPDA) ;                                      *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  2+IDAA))                     
      CALL DADIV(CCCC       ,ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))         
      CALL DAADD(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),Y          ((2)))                     
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
*FOX      TRACKI(1)=(X(1)+ED(IX)-DUMMY)*C1M3 ;                          *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = ED         (IX         )                       
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACSU(ISCRDA(  3+IDAA),ONE*DUMMY      ,ISCRDA(  4+IDAA))     
      CALL DACMU(ISCRDA(  4+IDAA),ONE*C1M3       ,ISCRDA(  5+IDAA))     
      CALL DACOP(ISCRDA(  5+IDAA),TRACKI     ((1)))                     
*FOX      YP(1)=Y(1)*(ONE+DPDA) ;                                       *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((1)))                     
          dummy=dare(yp(1))
*FOX      TRACKI(2)=(YP(1)-DUMMY)*C1M3 ;                                *FOX
      CALL DACOP(YP         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*DUMMY      ,ISCRDA(  2+IDAA))     
      CALL DACMU(ISCRDA(  2+IDAA),ONE*C1M3       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),TRACKI     ((2)))                     
          dummy=dare(x(2))
*FOX      TRACKI(3)=(X(2)+EK(IX)-DUMMY)*C1M3 ;                          *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = EK         (IX         )                       
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACSU(ISCRDA(  3+IDAA),ONE*DUMMY      ,ISCRDA(  4+IDAA))     
      CALL DACMU(ISCRDA(  4+IDAA),ONE*C1M3       ,ISCRDA(  5+IDAA))     
      CALL DACOP(ISCRDA(  5+IDAA),TRACKI     ((3)))                     
*FOX      YP(2)=Y(2)*(ONE+DPDA) ;                                       *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((2)))                     
          dummy=dare(yp(2))
*FOX      TRACKI(4)=(YP(2)-DUMMY)*C1M3 ;                                *FOX
      CALL DACOP(YP         ((2)),ISCRDA(  1+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*DUMMY      ,ISCRDA(  2+IDAA))     
      CALL DACMU(ISCRDA(  2+IDAA),ONE*C1M3       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),TRACKI     ((4)))                     
          dummy=dare(sigmda)
*FOX      TRACKI(5)=(SIGMDA-DUMMY)*C1M3 ;                               *FOX
      CALL DACSU(SIGMDA     ,ONE*DUMMY      ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1M3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),TRACKI     ((5)))                     
          dummy=dare(dpda)
*FOX      TRACKI(6)=DPDA-DUMMY ;                                        *FOX
      CALL DACSU(DPDA       ,ONE*DUMMY      ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),TRACKI     ((6)))                     
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
*FOX      X(1)=TRACKI(1)*C1E3+DUMMY-BEAMOFF1 ;                          *FOX
      CALL DACOP(TRACKI     ((1)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1E3       ,ISCRDA(  2+IDAA))     
      CALL DACAD(ISCRDA(  2+IDAA),ONE*DUMMY      ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  3+IDAA),ONE*BEAMOFF1   ,ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),X          ((1)))                     
          dummy=dare(x(2))
*FOX      X(2)=TRACKI(3)*C1E3+DUMMY-BEAMOFF2 ;                          *FOX
      CALL DACOP(TRACKI     ((3)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1E3       ,ISCRDA(  2+IDAA))     
      CALL DACAD(ISCRDA(  2+IDAA),ONE*DUMMY      ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  3+IDAA),ONE*BEAMOFF2   ,ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),X          ((2)))                     
          dummy=dare(yp(1))
*FOX      YP(1)=TRACKI(2)*C1E3+DUMMY-BEAMOFF4 ;                         *FOX
      CALL DACOP(TRACKI     ((2)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1E3       ,ISCRDA(  2+IDAA))     
      CALL DACAD(ISCRDA(  2+IDAA),ONE*DUMMY      ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  3+IDAA),ONE*BEAMOFF4   ,ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),YP         ((1)))                     
          dummy=dare(yp(2))
*FOX      YP(2)=TRACKI(4)*C1E3+DUMMY-BEAMOFF5 ;                         *FOX
      CALL DACOP(TRACKI     ((4)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1E3       ,ISCRDA(  2+IDAA))     
      CALL DACAD(ISCRDA(  2+IDAA),ONE*DUMMY      ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  3+IDAA),ONE*BEAMOFF5   ,ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),YP         ((2)))                     
          dummy=dare(dpda)
*FOX      DPDA=TRACKI(6)+DUMMY-BEAMOFF6 ;                               *FOX
      CALL DACOP(TRACKI     ((6)),ISCRDA(  1+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*DUMMY      ,ISCRDA(  2+IDAA))     
      CALL DACSU(ISCRDA(  2+IDAA),ONE*BEAMOFF6   ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),DPDA       )                          
*FOX      DPDA1=DPDA*C1E3 ;                                             *FOX
      CALL DACMU(DPDA       ,ONE*C1E3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA1      )                          
*FOX      Y(1)=YP(1)/(ONE+DPDA) ;                                       *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(YP         ((1)),ISCRDA(  2+IDAA))                     
      CALL DADIV(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((1)))                     
*FOX      Y(2)=YP(2)/(ONE+DPDA) ;                                       *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(YP         ((2)),ISCRDA(  2+IDAA))                     
      CALL DADIV(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((2)))                     
*FOX      EJF1=E0F*(ONE+DPDA) ;                                         *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),EJF1       )                          
*FOX      EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;                                 *FOX
      CALL DAMUL(EJF1       ,EJF1       ,ISCRDA(  1+IDAA))              
      RSCRRI(  2+IDAA) = PMA         * PMA                              
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DAFUN('SQRT',ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))              
      CALL DACOP(ISCRDA(  4+IDAA),EJ1        )                          
*FOX      RV=EJ1/E0*E0F/EJF1 ;                                          *FOX
      CALL DACDI(EJ1        ,ONE*E0         ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DADIV(ISCRDA(  2+IDAA),EJF1       ,ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),RV         )                          
          if(ithick.eq.1) call envada
          goto 440
        endif
          pi=4d0*atan_rn(1d0)
        if(kzz.eq.23) then
*FOX  CRABAMP=ED(IX)/(EJF1) ;                                           *FOX
      RSCRRI(  1+IDAA) = ED         (IX         )                       
      CALL DADIC(EJF1       ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),CRABAMP    )                          
!       call dapri(EJF1,234)                                            *FOX
!       write(*,*) crabamp, EJF1, EJF0,clight, "HELLO"                  *FOX
        crabfreq=ek(ix)*c1e3
        crabpht=crabph(ix)
*FOX  Y(1)=Y(1) - CRABAMP*C1E3*                                         *FOX
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT) ;               *FOX
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  1+IDAA))          
      CALL DACDI(ISCRDA(  1+IDAA),ONE*CLIGHT     ,ISCRDA(  2+IDAA))     
      CALL DACMU(ISCRDA(  2+IDAA),ONE*CRABFREQ   ,ISCRDA(  3+IDAA))     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*(2D0        ),ISCRDA(  4+IDAA))   
      CALL DACMU(ISCRDA(  4+IDAA),ONE*PI         ,ISCRDA(  5+IDAA))     
      CALL DACAD(ISCRDA(  5+IDAA),ONE*CRABPHT    ,ISCRDA(  6+IDAA))     
      CALL DACOP(Y          ((1)),ISCRDA(  7+IDAA))                     
      CALL DAFUN('SIN ',ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))              
      CALL DACMU(CRABAMP    ,ONE*C1E3       ,ISCRDA(  9+IDAA))          
      CALL DAMUL(ISCRDA(  9+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(ISCRDA( 11+IDAA),Y          ((1)))                     
*FOX  DPDA1=DPDA1 - CRABAMP*CRABFREQ*2D0*PI/CLIGHT*X(1)*                *FOX
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT) ;               *FOX
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  1+IDAA))          
      CALL DACDI(ISCRDA(  1+IDAA),ONE*CLIGHT     ,ISCRDA(  2+IDAA))     
      CALL DACMU(ISCRDA(  2+IDAA),ONE*CRABFREQ   ,ISCRDA(  3+IDAA))     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*(2D0        ),ISCRDA(  4+IDAA))   
      CALL DACMU(ISCRDA(  4+IDAA),ONE*PI         ,ISCRDA(  5+IDAA))     
      CALL DACAD(ISCRDA(  5+IDAA),ONE*CRABPHT    ,ISCRDA(  6+IDAA))     
      CALL DACOP(X          ((1)),ISCRDA(  7+IDAA))                     
      CALL DAFUN('COS ',ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))              
      CALL DACMU(CRABAMP    ,ONE*CRABFREQ   ,ISCRDA(  9+IDAA))          
      CALL DACMU(ISCRDA(  9+IDAA),ONE*(2D0        ),ISCRDA( 10+IDAA))   
      CALL DACMU(ISCRDA( 10+IDAA),ONE*PI         ,ISCRDA( 11+IDAA))     
      CALL DACDI(ISCRDA( 11+IDAA),ONE*CLIGHT     ,ISCRDA( 12+IDAA))     
      CALL DAMUL(ISCRDA( 12+IDAA),ISCRDA(  7+IDAA),ISCRDA( 13+IDAA))    
      CALL DAMUL(ISCRDA( 13+IDAA),ISCRDA(  8+IDAA),ISCRDA( 14+IDAA))    
      CALL DASUB(DPDA1      ,ISCRDA( 14+IDAA),ISCRDA( 15+IDAA))         
      CALL DACOP(ISCRDA( 15+IDAA),DPDA1      )                          
*FOX  EJF0=EJF1 ;                                                       *FOX
      CALL DACOP(EJF1       ,EJF0       )                               
*FOX  DPDA=DPDA1*C1M3 ;                                                 *FOX
      CALL DACMU(DPDA1      ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA       )                          
*FOX  EJF1=E0F*(ONE+DPDA) ;                                             *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),EJF1       )                          
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;                                     *FOX
      CALL DAMUL(EJF1       ,EJF1       ,ISCRDA(  1+IDAA))              
      RSCRRI(  2+IDAA) = PMA         * PMA                              
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DAFUN('SQRT',ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))              
      CALL DACOP(ISCRDA(  4+IDAA),EJ1        )                          
*FOX  Y(1)=EJF0/EJF1*Y(1) ;                                             *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((1)))                     
*FOX  Y(2)=EJF0/EJF1*Y(2) ;                                             *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((2)))                     
          goto 440
      endif
        if(kzz.eq.-23) then
*FOX  CRABAMP=ED(IX)/(EJF1) ;                                           *FOX
      RSCRRI(  1+IDAA) = ED         (IX         )                       
      CALL DADIC(EJF1       ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),CRABAMP    )                          
           crabfreq=ek(ix)*c1e3
           crabpht=crabph(ix)
*FOX  Y(2)=Y(2) - CRABAMP*C1E3*                                         *FOX
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT) ;               *FOX
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  1+IDAA))          
      CALL DACDI(ISCRDA(  1+IDAA),ONE*CLIGHT     ,ISCRDA(  2+IDAA))     
      CALL DACMU(ISCRDA(  2+IDAA),ONE*CRABFREQ   ,ISCRDA(  3+IDAA))     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*(2D0        ),ISCRDA(  4+IDAA))   
      CALL DACMU(ISCRDA(  4+IDAA),ONE*PI         ,ISCRDA(  5+IDAA))     
      CALL DACAD(ISCRDA(  5+IDAA),ONE*CRABPHT    ,ISCRDA(  6+IDAA))     
      CALL DACOP(Y          ((2)),ISCRDA(  7+IDAA))                     
      CALL DAFUN('SIN ',ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))              
      CALL DACMU(CRABAMP    ,ONE*C1E3       ,ISCRDA(  9+IDAA))          
      CALL DAMUL(ISCRDA(  9+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(ISCRDA( 11+IDAA),Y          ((2)))                     
*FOX  DPDA1=DPDA1 - CRABAMP*CRABFREQ*2D0*PI/CLIGHT*X(2)*                *FOX
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT) ;               *FOX
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  1+IDAA))          
      CALL DACDI(ISCRDA(  1+IDAA),ONE*CLIGHT     ,ISCRDA(  2+IDAA))     
      CALL DACMU(ISCRDA(  2+IDAA),ONE*CRABFREQ   ,ISCRDA(  3+IDAA))     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*(2D0        ),ISCRDA(  4+IDAA))   
      CALL DACMU(ISCRDA(  4+IDAA),ONE*PI         ,ISCRDA(  5+IDAA))     
      CALL DACAD(ISCRDA(  5+IDAA),ONE*CRABPHT    ,ISCRDA(  6+IDAA))     
      CALL DACOP(X          ((2)),ISCRDA(  7+IDAA))                     
      CALL DAFUN('COS ',ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))              
      CALL DACMU(CRABAMP    ,ONE*CRABFREQ   ,ISCRDA(  9+IDAA))          
      CALL DACMU(ISCRDA(  9+IDAA),ONE*(2D0        ),ISCRDA( 10+IDAA))   
      CALL DACMU(ISCRDA( 10+IDAA),ONE*PI         ,ISCRDA( 11+IDAA))     
      CALL DACDI(ISCRDA( 11+IDAA),ONE*CLIGHT     ,ISCRDA( 12+IDAA))     
      CALL DAMUL(ISCRDA( 12+IDAA),ISCRDA(  7+IDAA),ISCRDA( 13+IDAA))    
      CALL DAMUL(ISCRDA( 13+IDAA),ISCRDA(  8+IDAA),ISCRDA( 14+IDAA))    
      CALL DASUB(DPDA1      ,ISCRDA( 14+IDAA),ISCRDA( 15+IDAA))         
      CALL DACOP(ISCRDA( 15+IDAA),DPDA1      )                          
*FOX  EJF0=EJF1 ;                                                       *FOX
      CALL DACOP(EJF1       ,EJF0       )                               
*FOX  DPDA=DPDA1*C1M3 ;                                                 *FOX
      CALL DACMU(DPDA1      ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA       )                          
*FOX  EJF1=E0F*(ONE+DPDA) ;                                             *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),EJF1       )                          
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;                                     *FOX
      CALL DAMUL(EJF1       ,EJF1       ,ISCRDA(  1+IDAA))              
      RSCRRI(  2+IDAA) = PMA         * PMA                              
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DAFUN('SQRT',ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))              
      CALL DACOP(ISCRDA(  4+IDAA),EJ1        )                          
*FOX  Y(1)=EJF0/EJF1*Y(1) ;                                             *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((1)))                     
*FOX  Y(2)=EJF0/EJF1*Y(2) ;                                             *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((2)))                     
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
*FOX  XL=(X(1)-XS)*TILTC(I)+(X(2)-ZS)*TILTS(I) ;                        *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*XS         ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  2+IDAA),ONE*ZS         ,ISCRDA(  4+IDAA))     
      RSCRRI(  5+IDAA) = TILTC      (I          )                       
      RSCRRI(  6+IDAA) = TILTS      (I          )                       
      CALL DACMU(ISCRDA(  3+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),XL         )                          
*FOX  ZL=-(X(1)-XS)*TILTS(I)+(X(2)-ZS)*TILTC(I) ;                       *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*XS         ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  2+IDAA),ONE*ZS         ,ISCRDA(  4+IDAA))     
      RSCRRI(  5+IDAA) = TILTS      (I          )                       
      RSCRRI(  6+IDAA) = TILTC      (I          )                       
      CALL DACMU(ISCRDA(  3+IDAA),ONE*(-ONE       ),ISCRDA(  7+IDAA))   
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  9+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  9+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),ZL         )                          
*FOX  CRKVE=XL ;                                                        *FOX
      CALL DACOP(XL         ,CRKVE      )                               
*FOX  CIKVE=ZL ;                                                        *FOX
      CALL DACOP(ZL         ,CIKVE      )                               
*FOX  CRABAMP2=ED(IX)/(ONE+DPDA) ;                                      *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      RSCRRI(  2+IDAA) = ED         (IX         )                       
      CALL DADIC(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  3+IDAA),CRABAMP2   )                          
!       call dapri(EJF1,234)                                            *FOX
!       write(*,*) crabamp, EJF1, EJF0,clight, "HELLO"                  *FOX
        crabfreq=ek(ix)*c1e3 !JBG Input in MHz changed to kHz
        crabpht2=crabph2(ix)
*FOX  Y(1)=Y(1) + (CRABAMP2*CRKVE)*                                     *FOX
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT2);               *FOX
      CALL DAMUL(CRABAMP2   ,CRKVE      ,ISCRDA(  1+IDAA))              
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  2+IDAA))          
      CALL DACDI(ISCRDA(  2+IDAA),ONE*CLIGHT     ,ISCRDA(  3+IDAA))     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*CRABFREQ   ,ISCRDA(  4+IDAA))     
      CALL DACMU(ISCRDA(  4+IDAA),ONE*(2D0        ),ISCRDA(  5+IDAA))   
      CALL DACMU(ISCRDA(  5+IDAA),ONE*PI         ,ISCRDA(  6+IDAA))     
      CALL DACAD(ISCRDA(  6+IDAA),ONE*CRABPHT2   ,ISCRDA(  7+IDAA))     
      CALL DACOP(Y          ((1)),ISCRDA(  8+IDAA))                     
      CALL DAFUN('COS ',ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))              
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  9+IDAA),ISCRDA( 10+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(ISCRDA( 11+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2) - (CRABAMP2*CIKVE)*                                     *FOX
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT2);               *FOX
      CALL DAMUL(CRABAMP2   ,CIKVE      ,ISCRDA(  1+IDAA))              
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  2+IDAA))          
      CALL DACDI(ISCRDA(  2+IDAA),ONE*CLIGHT     ,ISCRDA(  3+IDAA))     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*CRABFREQ   ,ISCRDA(  4+IDAA))     
      CALL DACMU(ISCRDA(  4+IDAA),ONE*(2D0        ),ISCRDA(  5+IDAA))   
      CALL DACMU(ISCRDA(  5+IDAA),ONE*PI         ,ISCRDA(  6+IDAA))     
      CALL DACAD(ISCRDA(  6+IDAA),ONE*CRABPHT2   ,ISCRDA(  7+IDAA))     
      CALL DACOP(Y          ((2)),ISCRDA(  8+IDAA))                     
      CALL DAFUN('COS ',ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))              
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  9+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA(  8+IDAA),ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(ISCRDA( 11+IDAA),Y          ((2)))                     
*FOX  DPDA1=DPDA1 - (1/2.)*(CRABAMP2)*(CRKVE*CRKVE-                     *FOX
*FOX  CIKVE*CIKVE)*(((CRABFREQ*2D0)*PI)/CLIGHT)*C1M3*                   *FOX
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI+CRABPHT2) ;                *FOX
      RSCRRI(  1+IDAA) = CRABFREQ    * (2D0        )                    
      RSCRRI(  2+IDAA) = RSCRRI(  1+IDAA) * PI                          
      RSCRRI(  3+IDAA) = (1          ) / (2.         )                  
      CALL DAMUL(CRKVE      ,CRKVE      ,ISCRDA(  4+IDAA))              
      CALL DAMUL(CIKVE      ,CIKVE      ,ISCRDA(  5+IDAA))              
      RSCRRI(  6+IDAA) = RSCRRI(  2+IDAA) / CLIGHT                      
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  7+IDAA))          
      CALL DACDI(ISCRDA(  7+IDAA),ONE*CLIGHT     ,ISCRDA(  8+IDAA))     
      CALL DACMU(ISCRDA(  8+IDAA),ONE*CRABFREQ   ,ISCRDA(  9+IDAA))     
      CALL DACMU(ISCRDA(  9+IDAA),ONE*(2D0        ),ISCRDA( 10+IDAA))   
      CALL DACMU(ISCRDA( 10+IDAA),ONE*PI         ,ISCRDA( 11+IDAA))     
      CALL DASUB(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA( 12+IDAA))    
      CALL DACAD(ISCRDA( 11+IDAA),ONE*CRABPHT2   ,ISCRDA( 13+IDAA))     
      CALL DAFUN('SIN ',ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))              
      CALL DACMU(CRABAMP2   ,ONE*RSCRRI(  3+IDAA),ISCRDA( 15+IDAA))     
      CALL DAMUL(ISCRDA( 15+IDAA),ISCRDA( 12+IDAA),ISCRDA( 16+IDAA))    
      CALL DACMU(ISCRDA( 16+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA( 17+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 17+IDAA),ONE*C1M3       ,ISCRDA( 18+IDAA))     
      CALL DAMUL(ISCRDA( 18+IDAA),ISCRDA( 14+IDAA),ISCRDA( 19+IDAA))    
      CALL DASUB(DPDA1      ,ISCRDA( 19+IDAA),ISCRDA( 20+IDAA))         
      CALL DACOP(ISCRDA( 20+IDAA),DPDA1      )                          
*FOX  EJF0=EJF1 ;                                                       *FOX
      CALL DACOP(EJF1       ,EJF0       )                               
*FOX  DPDA=DPDA1*C1M3 ;                                                 *FOX
      CALL DACMU(DPDA1      ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA       )                          
*FOX  EJF1=E0F*(ONE+DPDA) ;                                             *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),EJF1       )                          
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;                                     *FOX
      CALL DAMUL(EJF1       ,EJF1       ,ISCRDA(  1+IDAA))              
      RSCRRI(  2+IDAA) = PMA         * PMA                              
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DAFUN('SQRT',ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))              
      CALL DACOP(ISCRDA(  4+IDAA),EJ1        )                          
*FOX  Y(1)=EJF0/EJF1*Y(1) ;                                             *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((1)))                     
*FOX  Y(2)=EJF0/EJF1*Y(2) ;                                             *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((2)))                     
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
*FOX  CRABAMP2=ED(IX)/(ONE+DPDA) ;                                      *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      RSCRRI(  2+IDAA) = ED         (IX         )                       
      CALL DADIC(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  3+IDAA),CRABAMP2   )                          
             crabfreq=ek(ix)*c1e3
             crabpht2=crabph2(ix)
*FOX  Y(2)=Y(2) + (CRABAMP2*CRKVE)*                                     *FOX
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT2);               *FOX
      CALL DAMUL(CRABAMP2   ,CRKVE      ,ISCRDA(  1+IDAA))              
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  2+IDAA))          
      CALL DACDI(ISCRDA(  2+IDAA),ONE*CLIGHT     ,ISCRDA(  3+IDAA))     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*CRABFREQ   ,ISCRDA(  4+IDAA))     
      CALL DACMU(ISCRDA(  4+IDAA),ONE*(2D0        ),ISCRDA(  5+IDAA))   
      CALL DACMU(ISCRDA(  5+IDAA),ONE*PI         ,ISCRDA(  6+IDAA))     
      CALL DACAD(ISCRDA(  6+IDAA),ONE*CRABPHT2   ,ISCRDA(  7+IDAA))     
      CALL DACOP(Y          ((2)),ISCRDA(  8+IDAA))                     
      CALL DAFUN('COS ',ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))              
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  9+IDAA),ISCRDA( 10+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(ISCRDA( 11+IDAA),Y          ((2)))                     
*FOX  Y(1)=Y(1) + (CRABAMP2*CIKVE)*                                     *FOX
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT2);               *FOX
      CALL DAMUL(CRABAMP2   ,CIKVE      ,ISCRDA(  1+IDAA))              
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  2+IDAA))          
      CALL DACDI(ISCRDA(  2+IDAA),ONE*CLIGHT     ,ISCRDA(  3+IDAA))     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*CRABFREQ   ,ISCRDA(  4+IDAA))     
      CALL DACMU(ISCRDA(  4+IDAA),ONE*(2D0        ),ISCRDA(  5+IDAA))   
      CALL DACMU(ISCRDA(  5+IDAA),ONE*PI         ,ISCRDA(  6+IDAA))     
      CALL DACAD(ISCRDA(  6+IDAA),ONE*CRABPHT2   ,ISCRDA(  7+IDAA))     
      CALL DACOP(Y          ((1)),ISCRDA(  8+IDAA))                     
      CALL DAFUN('COS ',ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))              
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  9+IDAA),ISCRDA( 10+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(ISCRDA( 11+IDAA),Y          ((1)))                     
*FOX  DPDA1=DPDA1 - (CRABAMP2)*(CIKVE*CRKVE)                            *FOX
*FOX  *(((CRABFREQ*2D0)*PI)/CLIGHT)*C1M3*
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI+CRABPHT2) ;                *FOX
      RSCRRI(  1+IDAA) = CRABFREQ    * (2D0        )                    
      RSCRRI(  2+IDAA) = RSCRRI(  1+IDAA) * PI                          
      CALL DAMUL(CIKVE      ,CRKVE      ,ISCRDA(  3+IDAA))              
      RSCRRI(  4+IDAA) = RSCRRI(  2+IDAA) / CLIGHT                      
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  5+IDAA))          
      CALL DACDI(ISCRDA(  5+IDAA),ONE*CLIGHT     ,ISCRDA(  6+IDAA))     
      CALL DACMU(ISCRDA(  6+IDAA),ONE*CRABFREQ   ,ISCRDA(  7+IDAA))     
      CALL DACMU(ISCRDA(  7+IDAA),ONE*(2D0        ),ISCRDA(  8+IDAA))   
      CALL DACMU(ISCRDA(  8+IDAA),ONE*PI         ,ISCRDA(  9+IDAA))     
      CALL DACAD(ISCRDA(  9+IDAA),ONE*CRABPHT2   ,ISCRDA( 10+IDAA))     
      CALL DAFUN('SIN ',ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))              
      CALL DAMUL(CRABAMP2   ,ISCRDA(  3+IDAA),ISCRDA( 12+IDAA))         
      CALL DACMU(ISCRDA( 12+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA( 13+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 13+IDAA),ONE*C1M3       ,ISCRDA( 14+IDAA))     
      CALL DAMUL(ISCRDA( 14+IDAA),ISCRDA( 11+IDAA),ISCRDA( 15+IDAA))    
      CALL DASUB(DPDA1      ,ISCRDA( 15+IDAA),ISCRDA( 16+IDAA))         
      CALL DACOP(ISCRDA( 16+IDAA),DPDA1      )                          
*FOX  EJF0=EJF1 ;                                                       *FOX
      CALL DACOP(EJF1       ,EJF0       )                               
*FOX  DPDA=DPDA1*C1M3 ;                                                 *FOX
      CALL DACMU(DPDA1      ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA       )                          
*FOX  EJF1=E0F*(ONE+DPDA) ;                                             *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),EJF1       )                          
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;                                     *FOX
      CALL DAMUL(EJF1       ,EJF1       ,ISCRDA(  1+IDAA))              
      RSCRRI(  2+IDAA) = PMA         * PMA                              
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DAFUN('SQRT',ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))              
      CALL DACOP(ISCRDA(  4+IDAA),EJ1        )                          
*FOX  Y(1)=EJF0/EJF1*Y(1) ;                                             *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((1)))                     
*FOX  Y(2)=EJF0/EJF1*Y(2) ;                                             *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((2)))                     
          endif
          if(kzz.eq.27) then
            ! JBG bypass this element if 4D/5D case
            if(iclo6.eq.0) then
!                write(*,*)'Bypassing RF mult 4D or 5D case'
                goto 440
            endif
          xs=xsi(i)
          zs=zsi(i)
*FOX  XL=(X(1)-XS)*TILTC(I)+(X(2)-ZS)*TILTS(I) ;                        *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*XS         ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  2+IDAA),ONE*ZS         ,ISCRDA(  4+IDAA))     
      RSCRRI(  5+IDAA) = TILTC      (I          )                       
      RSCRRI(  6+IDAA) = TILTS      (I          )                       
      CALL DACMU(ISCRDA(  3+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),XL         )                          
*FOX  ZL=-(X(1)-XS)*TILTS(I)+(X(2)-ZS)*TILTC(I) ;                       *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*XS         ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  2+IDAA),ONE*ZS         ,ISCRDA(  4+IDAA))     
      RSCRRI(  5+IDAA) = TILTS      (I          )                       
      RSCRRI(  6+IDAA) = TILTC      (I          )                       
      CALL DACMU(ISCRDA(  3+IDAA),ONE*(-ONE       ),ISCRDA(  7+IDAA))   
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  9+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  9+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),ZL         )                          
*FOX  CRKVE=XL ;                                                        *FOX
      CALL DACOP(XL         ,CRKVE      )                               
*FOX  CIKVE=ZL ;                                                        *FOX
      CALL DACOP(ZL         ,CIKVE      )                               
*FOX  CRABAMP3=ED(IX)/(ONE+DPDA) ;                                      *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      RSCRRI(  2+IDAA) = ED         (IX         )                       
      CALL DADIC(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  3+IDAA),CRABAMP3   )                          
             crabfreq=ek(ix)*c1e3
             crabpht3=crabph3(ix)
*FOX  Y(1)=Y(1) + 2*(1/2.)*CRABAMP3*((CRKVE*CRKVE)-                     *FOX
*FOX  (CIKVE*CIKVE))*C1M3*                                              *FOX
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT3);               *FOX
      CALL DAMUL(CRKVE      ,CRKVE      ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,CIKVE      ,ISCRDA(  2+IDAA))              
      RSCRRI(  3+IDAA) = (1          ) / (2.         )                  
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  4+IDAA))          
      CALL DACDI(ISCRDA(  4+IDAA),ONE*CLIGHT     ,ISCRDA(  5+IDAA))     
      CALL DACMU(ISCRDA(  5+IDAA),ONE*CRABFREQ   ,ISCRDA(  6+IDAA))     
      CALL DACMU(ISCRDA(  6+IDAA),ONE*(2D0        ),ISCRDA(  7+IDAA))   
      CALL DACMU(ISCRDA(  7+IDAA),ONE*PI         ,ISCRDA(  8+IDAA))     
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  9+IDAA))    
      CALL DACAD(ISCRDA(  8+IDAA),ONE*CRABPHT3   ,ISCRDA( 10+IDAA))     
      CALL DACOP(Y          ((1)),ISCRDA( 11+IDAA))                     
      CALL DAFUN('COS ',ISCRDA( 10+IDAA),ISCRDA( 12+IDAA))              
      RSCRRI( 13+IDAA) = (2          ) * RSCRRI(  3+IDAA)               
      CALL DACMU(CRABAMP3   ,ONE*RSCRRI( 13+IDAA),ISCRDA( 14+IDAA))     
      CALL DAMUL(ISCRDA( 14+IDAA),ISCRDA(  9+IDAA),ISCRDA( 15+IDAA))    
      CALL DACMU(ISCRDA( 15+IDAA),ONE*C1M3       ,ISCRDA( 16+IDAA))     
      CALL DAMUL(ISCRDA( 16+IDAA),ISCRDA( 12+IDAA),ISCRDA( 17+IDAA))    
      CALL DAADD(ISCRDA( 11+IDAA),ISCRDA( 17+IDAA),ISCRDA( 18+IDAA))    
      CALL DACOP(ISCRDA( 18+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2) - 2*CRABAMP3*(CRKVE*CIKVE)*C1M3*                        *FOX
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT3);               *FOX
      CALL DAMUL(CRKVE      ,CIKVE      ,ISCRDA(  1+IDAA))              
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  2+IDAA))          
      CALL DACDI(ISCRDA(  2+IDAA),ONE*CLIGHT     ,ISCRDA(  3+IDAA))     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*CRABFREQ   ,ISCRDA(  4+IDAA))     
      CALL DACMU(ISCRDA(  4+IDAA),ONE*(2D0        ),ISCRDA(  5+IDAA))   
      CALL DACMU(ISCRDA(  5+IDAA),ONE*PI         ,ISCRDA(  6+IDAA))     
      CALL DACAD(ISCRDA(  6+IDAA),ONE*CRABPHT3   ,ISCRDA(  7+IDAA))     
      CALL DACOP(Y          ((2)),ISCRDA(  8+IDAA))                     
      CALL DAFUN('COS ',ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))              
      CALL DACMU(CRABAMP3   ,ONE*(2          ),ISCRDA( 10+IDAA))        
      CALL DAMUL(ISCRDA( 10+IDAA),ISCRDA(  1+IDAA),ISCRDA( 11+IDAA))    
      CALL DACMU(ISCRDA( 11+IDAA),ONE*C1M3       ,ISCRDA( 12+IDAA))     
      CALL DAMUL(ISCRDA( 12+IDAA),ISCRDA(  9+IDAA),ISCRDA( 13+IDAA))    
      CALL DASUB(ISCRDA(  8+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DACOP(ISCRDA( 14+IDAA),Y          ((2)))                     
*FOX  DPDA1=DPDA1 - 2*(1/6.)*(CRABAMP3)*(CRKVE*CRKVE*CRKVE-             *FOX
*FOX  3*CIKVE*CIKVE*CRKVE)*(((CRABFREQ*2D0)*PI)/CLIGHT)*                *FOX
*FOX  C1M6*
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI+CRABPHT3) ;                *FOX
      RSCRRI(  1+IDAA) = CRABFREQ    * (2D0        )                    
      RSCRRI(  2+IDAA) = RSCRRI(  1+IDAA) * PI                          
      RSCRRI(  3+IDAA) = (1          ) / (6.         )                  
      CALL DAMUL(CRKVE      ,CRKVE      ,ISCRDA(  4+IDAA))              
      CALL DAMUL(ISCRDA(  4+IDAA),CRKVE      ,ISCRDA(  5+IDAA))         
      CALL DACMU(CIKVE      ,ONE*(3          ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),CIKVE      ,ISCRDA(  7+IDAA))         
      CALL DAMUL(ISCRDA(  7+IDAA),CRKVE      ,ISCRDA(  8+IDAA))         
      RSCRRI(  9+IDAA) = RSCRRI(  2+IDAA) / CLIGHT                      
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA( 10+IDAA))          
      CALL DACDI(ISCRDA( 10+IDAA),ONE*CLIGHT     ,ISCRDA( 11+IDAA))     
      CALL DACMU(ISCRDA( 11+IDAA),ONE*CRABFREQ   ,ISCRDA( 12+IDAA))     
      CALL DACMU(ISCRDA( 12+IDAA),ONE*(2D0        ),ISCRDA( 13+IDAA))   
      CALL DACMU(ISCRDA( 13+IDAA),ONE*PI         ,ISCRDA( 14+IDAA))     
      CALL DASUB(ISCRDA(  5+IDAA),ISCRDA(  8+IDAA),ISCRDA( 15+IDAA))    
      CALL DACAD(ISCRDA( 14+IDAA),ONE*CRABPHT3   ,ISCRDA( 16+IDAA))     
      CALL DAFUN('SIN ',ISCRDA( 16+IDAA),ISCRDA( 17+IDAA))              
      RSCRRI( 18+IDAA) = (2          ) * RSCRRI(  3+IDAA)               
      CALL DACMU(CRABAMP3   ,ONE*RSCRRI( 18+IDAA),ISCRDA( 19+IDAA))     
      CALL DAMUL(ISCRDA( 19+IDAA),ISCRDA( 15+IDAA),ISCRDA( 20+IDAA))    
      CALL DACMU(ISCRDA( 20+IDAA),ONE*RSCRRI(  9+IDAA),ISCRDA( 21+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 21+IDAA),ONE*C1M6       ,ISCRDA( 22+IDAA))     
      CALL DAMUL(ISCRDA( 22+IDAA),ISCRDA( 17+IDAA),ISCRDA( 23+IDAA))    
      CALL DASUB(DPDA1      ,ISCRDA( 23+IDAA),ISCRDA( 24+IDAA))         
      CALL DACOP(ISCRDA( 24+IDAA),DPDA1      )                          
*FOX  EJF0=EJF1 ;                                                       *FOX
      CALL DACOP(EJF1       ,EJF0       )                               
*FOX  DPDA=DPDA1*C1M3 ;                                                 *FOX
      CALL DACMU(DPDA1      ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA       )                          
*FOX  EJF1=E0F*(ONE+DPDA) ;                                             *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),EJF1       )                          
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;                                     *FOX
      CALL DAMUL(EJF1       ,EJF1       ,ISCRDA(  1+IDAA))              
      RSCRRI(  2+IDAA) = PMA         * PMA                              
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DAFUN('SQRT',ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))              
      CALL DACOP(ISCRDA(  4+IDAA),EJ1        )                          
*FOX  Y(1)=EJF0/EJF1*Y(1) ;                                             *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((1)))                     
*FOX  Y(2)=EJF0/EJF1*Y(2) ;                                             *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((2)))                     
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
*FOX  XL=(X(1)-XS)*TILTC(I)+(X(2)-ZS)*TILTS(I) ;                        *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*XS         ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  2+IDAA),ONE*ZS         ,ISCRDA(  4+IDAA))     
      RSCRRI(  5+IDAA) = TILTC      (I          )                       
      RSCRRI(  6+IDAA) = TILTS      (I          )                       
      CALL DACMU(ISCRDA(  3+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),XL         )                          
*FOX  ZL=-(X(1)-XS)*TILTS(I)+(X(2)-ZS)*TILTC(I) ;                       *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*XS         ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  2+IDAA),ONE*ZS         ,ISCRDA(  4+IDAA))     
      RSCRRI(  5+IDAA) = TILTS      (I          )                       
      RSCRRI(  6+IDAA) = TILTC      (I          )                       
      CALL DACMU(ISCRDA(  3+IDAA),ONE*(-ONE       ),ISCRDA(  7+IDAA))   
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  9+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  9+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),ZL         )                          
*FOX  CRKVE=XL ;                                                        *FOX
      CALL DACOP(XL         ,CRKVE      )                               
*FOX  CIKVE=ZL ;                                                        *FOX
      CALL DACOP(ZL         ,CIKVE      )                               
*FOX  CRABAMP3=ED(IX)/(ONE+DPDA) ;                                      *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      RSCRRI(  2+IDAA) = ED         (IX         )                       
      CALL DADIC(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  3+IDAA),CRABAMP3   )                          
             crabfreq=ek(ix)*c1e3
             crabpht3=crabph3(ix)
*FOX  Y(2)=Y(2) - 2*(1/2.)*CRABAMP3*((CIKVE*CIKVE)-                     *FOX
*FOX  (CRKVE*CRKVE))*C1M3*                                              *FOX
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT3);               *FOX
      CALL DAMUL(CIKVE      ,CIKVE      ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CRKVE      ,CRKVE      ,ISCRDA(  2+IDAA))              
      RSCRRI(  3+IDAA) = (1          ) / (2.         )                  
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  4+IDAA))          
      CALL DACDI(ISCRDA(  4+IDAA),ONE*CLIGHT     ,ISCRDA(  5+IDAA))     
      CALL DACMU(ISCRDA(  5+IDAA),ONE*CRABFREQ   ,ISCRDA(  6+IDAA))     
      CALL DACMU(ISCRDA(  6+IDAA),ONE*(2D0        ),ISCRDA(  7+IDAA))   
      CALL DACMU(ISCRDA(  7+IDAA),ONE*PI         ,ISCRDA(  8+IDAA))     
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  9+IDAA))    
      CALL DACAD(ISCRDA(  8+IDAA),ONE*CRABPHT3   ,ISCRDA( 10+IDAA))     
      CALL DACOP(Y          ((2)),ISCRDA( 11+IDAA))                     
      CALL DAFUN('COS ',ISCRDA( 10+IDAA),ISCRDA( 12+IDAA))              
      RSCRRI( 13+IDAA) = (2          ) * RSCRRI(  3+IDAA)               
      CALL DACMU(CRABAMP3   ,ONE*RSCRRI( 13+IDAA),ISCRDA( 14+IDAA))     
      CALL DAMUL(ISCRDA( 14+IDAA),ISCRDA(  9+IDAA),ISCRDA( 15+IDAA))    
      CALL DACMU(ISCRDA( 15+IDAA),ONE*C1M3       ,ISCRDA( 16+IDAA))     
      CALL DAMUL(ISCRDA( 16+IDAA),ISCRDA( 12+IDAA),ISCRDA( 17+IDAA))    
      CALL DASUB(ISCRDA( 11+IDAA),ISCRDA( 17+IDAA),ISCRDA( 18+IDAA))    
      CALL DACOP(ISCRDA( 18+IDAA),Y          ((2)))                     
*FOX  Y(1)=Y(1) + 2*CRABAMP3*(CRKVE*CIKVE)*C1M3*                        *FOX
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT3);               *FOX
      CALL DAMUL(CRKVE      ,CIKVE      ,ISCRDA(  1+IDAA))              
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  2+IDAA))          
      CALL DACDI(ISCRDA(  2+IDAA),ONE*CLIGHT     ,ISCRDA(  3+IDAA))     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*CRABFREQ   ,ISCRDA(  4+IDAA))     
      CALL DACMU(ISCRDA(  4+IDAA),ONE*(2D0        ),ISCRDA(  5+IDAA))   
      CALL DACMU(ISCRDA(  5+IDAA),ONE*PI         ,ISCRDA(  6+IDAA))     
      CALL DACAD(ISCRDA(  6+IDAA),ONE*CRABPHT3   ,ISCRDA(  7+IDAA))     
      CALL DACOP(Y          ((1)),ISCRDA(  8+IDAA))                     
      CALL DAFUN('COS ',ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))              
      CALL DACMU(CRABAMP3   ,ONE*(2          ),ISCRDA( 10+IDAA))        
      CALL DAMUL(ISCRDA( 10+IDAA),ISCRDA(  1+IDAA),ISCRDA( 11+IDAA))    
      CALL DACMU(ISCRDA( 11+IDAA),ONE*C1M3       ,ISCRDA( 12+IDAA))     
      CALL DAMUL(ISCRDA( 12+IDAA),ISCRDA(  9+IDAA),ISCRDA( 13+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DACOP(ISCRDA( 14+IDAA),Y          ((1)))                     
*FOX  DPDA1=DPDA1 + 2*(1/6.)*(CRABAMP3)*(CIKVE*CIKVE*CIKVE-             *FOX
*FOX  3*CIKVE*CRKVE*CRKVE)*(((CRABFREQ*2D0)*PI)/CLIGHT)*                *FOX
*FOX  C1M6*
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI+CRABPHT3) ;                *FOX
      RSCRRI(  1+IDAA) = CRABFREQ    * (2D0        )                    
      RSCRRI(  2+IDAA) = RSCRRI(  1+IDAA) * PI                          
      RSCRRI(  3+IDAA) = (1          ) / (6.         )                  
      CALL DAMUL(CIKVE      ,CIKVE      ,ISCRDA(  4+IDAA))              
      CALL DAMUL(ISCRDA(  4+IDAA),CIKVE      ,ISCRDA(  5+IDAA))         
      CALL DACMU(CIKVE      ,ONE*(3          ),ISCRDA(  6+IDAA))        
      CALL DAMUL(ISCRDA(  6+IDAA),CRKVE      ,ISCRDA(  7+IDAA))         
      CALL DAMUL(ISCRDA(  7+IDAA),CRKVE      ,ISCRDA(  8+IDAA))         
      RSCRRI(  9+IDAA) = RSCRRI(  2+IDAA) / CLIGHT                      
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA( 10+IDAA))          
      CALL DACDI(ISCRDA( 10+IDAA),ONE*CLIGHT     ,ISCRDA( 11+IDAA))     
      CALL DACMU(ISCRDA( 11+IDAA),ONE*CRABFREQ   ,ISCRDA( 12+IDAA))     
      CALL DACMU(ISCRDA( 12+IDAA),ONE*(2D0        ),ISCRDA( 13+IDAA))   
      CALL DACMU(ISCRDA( 13+IDAA),ONE*PI         ,ISCRDA( 14+IDAA))     
      CALL DASUB(ISCRDA(  5+IDAA),ISCRDA(  8+IDAA),ISCRDA( 15+IDAA))    
      CALL DACAD(ISCRDA( 14+IDAA),ONE*CRABPHT3   ,ISCRDA( 16+IDAA))     
      CALL DAFUN('SIN ',ISCRDA( 16+IDAA),ISCRDA( 17+IDAA))              
      RSCRRI( 18+IDAA) = (2          ) * RSCRRI(  3+IDAA)               
      CALL DACMU(CRABAMP3   ,ONE*RSCRRI( 18+IDAA),ISCRDA( 19+IDAA))     
      CALL DAMUL(ISCRDA( 19+IDAA),ISCRDA( 15+IDAA),ISCRDA( 20+IDAA))    
      CALL DACMU(ISCRDA( 20+IDAA),ONE*RSCRRI(  9+IDAA),ISCRDA( 21+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 21+IDAA),ONE*C1M6       ,ISCRDA( 22+IDAA))     
      CALL DAMUL(ISCRDA( 22+IDAA),ISCRDA( 17+IDAA),ISCRDA( 23+IDAA))    
      CALL DAADD(DPDA1      ,ISCRDA( 23+IDAA),ISCRDA( 24+IDAA))         
      CALL DACOP(ISCRDA( 24+IDAA),DPDA1      )                          
*FOX  EJF0=EJF1 ;                                                       *FOX
      CALL DACOP(EJF1       ,EJF0       )                               
*FOX  DPDA=DPDA1*C1M3 ;                                                 *FOX
      CALL DACMU(DPDA1      ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA       )                          
*FOX  EJF1=E0F*(ONE+DPDA) ;                                             *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),EJF1       )                          
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;                                     *FOX
      CALL DAMUL(EJF1       ,EJF1       ,ISCRDA(  1+IDAA))              
      RSCRRI(  2+IDAA) = PMA         * PMA                              
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DAFUN('SQRT',ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))              
      CALL DACOP(ISCRDA(  4+IDAA),EJ1        )                          
*FOX  Y(1)=EJF0/EJF1*Y(1) ;                                             *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((1)))                     
*FOX  Y(2)=EJF0/EJF1*Y(2) ;                                             *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((2)))                     
          endif
          if(kzz.eq.28) then
            ! JBG bypass this element if 4D/5D case
            if(iclo6.eq.0) then
!                write(*,*)'Bypassing RF mult 4D or 5D case'
                goto 440
            endif
          xs=xsi(i)
          zs=zsi(i)
*FOX  XL=(X(1)-XS)*TILTC(I)+(X(2)-ZS)*TILTS(I) ;                        *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*XS         ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  2+IDAA),ONE*ZS         ,ISCRDA(  4+IDAA))     
      RSCRRI(  5+IDAA) = TILTC      (I          )                       
      RSCRRI(  6+IDAA) = TILTS      (I          )                       
      CALL DACMU(ISCRDA(  3+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),XL         )                          
*FOX  ZL=-(X(1)-XS)*TILTS(I)+(X(2)-ZS)*TILTC(I) ;                       *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*XS         ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  2+IDAA),ONE*ZS         ,ISCRDA(  4+IDAA))     
      RSCRRI(  5+IDAA) = TILTS      (I          )                       
      RSCRRI(  6+IDAA) = TILTC      (I          )                       
      CALL DACMU(ISCRDA(  3+IDAA),ONE*(-ONE       ),ISCRDA(  7+IDAA))   
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  9+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  9+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),ZL         )                          
*FOX  CRKVE=XL ;                                                        *FOX
      CALL DACOP(XL         ,CRKVE      )                               
*FOX  CIKVE=ZL ;                                                        *FOX
      CALL DACOP(ZL         ,CIKVE      )                               
*FOX  CRABAMP4=ED(IX)/(ONE+DPDA) ;                                      *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      RSCRRI(  2+IDAA) = ED         (IX         )                       
      CALL DADIC(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  3+IDAA),CRABAMP4   )                          
             crabfreq=ek(ix)*c1e3
             crabpht4=crabph4(ix)
*FOX  Y(1)=Y(1) + 6*(1/6.)*(CRABAMP4)*                                  *FOX
*FOX  (CRKVE*CRKVE*CRKVE-(3*CRKVE*CIKVE*CIKVE))*C1M3*C1M3*              *FOX
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT4);               *FOX
      CALL DACMU(CRKVE      ,ONE*(3          ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),CIKVE      ,ISCRDA(  2+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),CIKVE      ,ISCRDA(  3+IDAA))         
      RSCRRI(  4+IDAA) = (1          ) / (6.         )                  
      CALL DAMUL(CRKVE      ,CRKVE      ,ISCRDA(  5+IDAA))              
      CALL DAMUL(ISCRDA(  5+IDAA),CRKVE      ,ISCRDA(  6+IDAA))         
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  7+IDAA))          
      CALL DACDI(ISCRDA(  7+IDAA),ONE*CLIGHT     ,ISCRDA(  8+IDAA))     
      CALL DACMU(ISCRDA(  8+IDAA),ONE*CRABFREQ   ,ISCRDA(  9+IDAA))     
      CALL DACMU(ISCRDA(  9+IDAA),ONE*(2D0        ),ISCRDA( 10+IDAA))   
      CALL DACMU(ISCRDA( 10+IDAA),ONE*PI         ,ISCRDA( 11+IDAA))     
      CALL DASUB(ISCRDA(  6+IDAA),ISCRDA(  3+IDAA),ISCRDA( 12+IDAA))    
      CALL DACAD(ISCRDA( 11+IDAA),ONE*CRABPHT4   ,ISCRDA( 13+IDAA))     
      CALL DACOP(Y          ((1)),ISCRDA( 14+IDAA))                     
      CALL DAFUN('COS ',ISCRDA( 13+IDAA),ISCRDA( 15+IDAA))              
      RSCRRI( 16+IDAA) = (6          ) * RSCRRI(  4+IDAA)               
      CALL DACMU(CRABAMP4   ,ONE*RSCRRI( 16+IDAA),ISCRDA( 17+IDAA))     
      CALL DAMUL(ISCRDA( 17+IDAA),ISCRDA( 12+IDAA),ISCRDA( 18+IDAA))    
      CALL DACMU(ISCRDA( 18+IDAA),ONE*C1M3       ,ISCRDA( 19+IDAA))     
      CALL DACMU(ISCRDA( 19+IDAA),ONE*C1M3       ,ISCRDA( 20+IDAA))     
      CALL DAMUL(ISCRDA( 20+IDAA),ISCRDA( 15+IDAA),ISCRDA( 21+IDAA))    
      CALL DAADD(ISCRDA( 14+IDAA),ISCRDA( 21+IDAA),ISCRDA( 22+IDAA))    
      CALL DACOP(ISCRDA( 22+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2) - 6*(1/6.)*(CRABAMP4)*                                  *FOX
*FOX  (3*CIKVE*CRKVE*CRKVE-CIKVE*CIKVE*CIKVE)*C1M3*C1M3*                *FOX
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT4);               *FOX
      RSCRRI(  1+IDAA) = (1          ) / (6.         )                  
      CALL DACMU(CIKVE      ,ONE*(3          ),ISCRDA(  2+IDAA))        
      CALL DAMUL(ISCRDA(  2+IDAA),CRKVE      ,ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  3+IDAA),CRKVE      ,ISCRDA(  4+IDAA))         
      CALL DAMUL(CIKVE      ,CIKVE      ,ISCRDA(  5+IDAA))              
      CALL DAMUL(ISCRDA(  5+IDAA),CIKVE      ,ISCRDA(  6+IDAA))         
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  7+IDAA))          
      CALL DACDI(ISCRDA(  7+IDAA),ONE*CLIGHT     ,ISCRDA(  8+IDAA))     
      CALL DACMU(ISCRDA(  8+IDAA),ONE*CRABFREQ   ,ISCRDA(  9+IDAA))     
      CALL DACMU(ISCRDA(  9+IDAA),ONE*(2D0        ),ISCRDA( 10+IDAA))   
      CALL DACMU(ISCRDA( 10+IDAA),ONE*PI         ,ISCRDA( 11+IDAA))     
      CALL DASUB(ISCRDA(  4+IDAA),ISCRDA(  6+IDAA),ISCRDA( 12+IDAA))    
      CALL DACAD(ISCRDA( 11+IDAA),ONE*CRABPHT4   ,ISCRDA( 13+IDAA))     
      CALL DACOP(Y          ((2)),ISCRDA( 14+IDAA))                     
      CALL DAFUN('COS ',ISCRDA( 13+IDAA),ISCRDA( 15+IDAA))              
      RSCRRI( 16+IDAA) = (6          ) * RSCRRI(  1+IDAA)               
      CALL DACMU(CRABAMP4   ,ONE*RSCRRI( 16+IDAA),ISCRDA( 17+IDAA))     
      CALL DAMUL(ISCRDA( 17+IDAA),ISCRDA( 12+IDAA),ISCRDA( 18+IDAA))    
      CALL DACMU(ISCRDA( 18+IDAA),ONE*C1M3       ,ISCRDA( 19+IDAA))     
      CALL DACMU(ISCRDA( 19+IDAA),ONE*C1M3       ,ISCRDA( 20+IDAA))     
      CALL DAMUL(ISCRDA( 20+IDAA),ISCRDA( 15+IDAA),ISCRDA( 21+IDAA))    
      CALL DASUB(ISCRDA( 14+IDAA),ISCRDA( 21+IDAA),ISCRDA( 22+IDAA))    
      CALL DACOP(ISCRDA( 22+IDAA),Y          ((2)))                     
*FOX  DPDA1=DPDA1-6*(1/24.)*(CRABAMP4)*(CRKVE*CRKVE*CRKVE*CRKVE-        *FOX
*FOX  6*CRKVE*CRKVE*CIKVE*CIKVE+CIKVE*CIKVE*CIKVE*CIKVE)*               *FOX
*FOX  C1M9*(((CRABFREQ*2D0)*PI)/CLIGHT)*
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI+CRABPHT4) ;                *FOX
      RSCRRI(  1+IDAA) = CRABFREQ    * (2D0        )                    
      RSCRRI(  2+IDAA) = RSCRRI(  1+IDAA) * PI                          
      RSCRRI(  3+IDAA) = (1          ) / (24.        )                  
      CALL DAMUL(CRKVE      ,CRKVE      ,ISCRDA(  4+IDAA))              
      CALL DAMUL(ISCRDA(  4+IDAA),CRKVE      ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  5+IDAA),CRKVE      ,ISCRDA(  6+IDAA))         
      CALL DACMU(CRKVE      ,ONE*(6          ),ISCRDA(  7+IDAA))        
      CALL DAMUL(ISCRDA(  7+IDAA),CRKVE      ,ISCRDA(  8+IDAA))         
      CALL DAMUL(ISCRDA(  8+IDAA),CIKVE      ,ISCRDA(  9+IDAA))         
      CALL DAMUL(ISCRDA(  9+IDAA),CIKVE      ,ISCRDA( 10+IDAA))         
      CALL DAMUL(CIKVE      ,CIKVE      ,ISCRDA( 11+IDAA))              
      CALL DAMUL(ISCRDA( 11+IDAA),CIKVE      ,ISCRDA( 12+IDAA))         
      CALL DAMUL(ISCRDA( 12+IDAA),CIKVE      ,ISCRDA( 13+IDAA))         
      RSCRRI( 14+IDAA) = RSCRRI(  2+IDAA) / CLIGHT                      
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA( 15+IDAA))          
      CALL DACDI(ISCRDA( 15+IDAA),ONE*CLIGHT     ,ISCRDA( 16+IDAA))     
      CALL DACMU(ISCRDA( 16+IDAA),ONE*CRABFREQ   ,ISCRDA( 17+IDAA))     
      CALL DACMU(ISCRDA( 17+IDAA),ONE*(2D0        ),ISCRDA( 18+IDAA))   
      CALL DACMU(ISCRDA( 18+IDAA),ONE*PI         ,ISCRDA( 19+IDAA))     
      CALL DASUB(ISCRDA(  6+IDAA),ISCRDA( 10+IDAA),ISCRDA( 20+IDAA))    
      CALL DAADD(ISCRDA( 20+IDAA),ISCRDA( 13+IDAA),ISCRDA( 21+IDAA))    
      CALL DACAD(ISCRDA( 19+IDAA),ONE*CRABPHT4   ,ISCRDA( 22+IDAA))     
      CALL DAFUN('SIN ',ISCRDA( 22+IDAA),ISCRDA( 23+IDAA))              
      RSCRRI( 24+IDAA) = (6          ) * RSCRRI(  3+IDAA)               
      CALL DACMU(CRABAMP4   ,ONE*RSCRRI( 24+IDAA),ISCRDA( 25+IDAA))     
      CALL DAMUL(ISCRDA( 25+IDAA),ISCRDA( 21+IDAA),ISCRDA( 26+IDAA))    
      CALL DACMU(ISCRDA( 26+IDAA),ONE*C1M9       ,ISCRDA( 27+IDAA))     
      CALL DACMU(ISCRDA( 27+IDAA),ONE*RSCRRI( 14+IDAA),ISCRDA( 28+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA( 28+IDAA),ISCRDA( 23+IDAA),ISCRDA( 29+IDAA))    
      CALL DASUB(DPDA1      ,ISCRDA( 29+IDAA),ISCRDA( 30+IDAA))         
      CALL DACOP(ISCRDA( 30+IDAA),DPDA1      )                          
*FOX  EJF0=EJF1 ;                                                       *FOX
      CALL DACOP(EJF1       ,EJF0       )                               
*FOX  DPDA=DPDA1*C1M3 ;                                                 *FOX
      CALL DACMU(DPDA1      ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA       )                          
*FOX  EJF1=E0F*(ONE+DPDA) ;                                             *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),EJF1       )                          
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;                                     *FOX
      CALL DAMUL(EJF1       ,EJF1       ,ISCRDA(  1+IDAA))              
      RSCRRI(  2+IDAA) = PMA         * PMA                              
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DAFUN('SQRT',ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))              
      CALL DACOP(ISCRDA(  4+IDAA),EJ1        )                          
*FOX  Y(1)=EJF0/EJF1*Y(1) ;                                             *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((1)))                     
*FOX  Y(2)=EJF0/EJF1*Y(2) ;                                             *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((2)))                     
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
*FOX  XL=(X(1)-XS)*TILTC(I)+(X(2)-ZS)*TILTS(I) ;                        *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*XS         ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  2+IDAA),ONE*ZS         ,ISCRDA(  4+IDAA))     
      RSCRRI(  5+IDAA) = TILTC      (I          )                       
      RSCRRI(  6+IDAA) = TILTS      (I          )                       
      CALL DACMU(ISCRDA(  3+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),XL         )                          
*FOX  ZL=-(X(1)-XS)*TILTS(I)+(X(2)-ZS)*TILTC(I) ;                       *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*XS         ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  2+IDAA),ONE*ZS         ,ISCRDA(  4+IDAA))     
      RSCRRI(  5+IDAA) = TILTS      (I          )                       
      RSCRRI(  6+IDAA) = TILTC      (I          )                       
      CALL DACMU(ISCRDA(  3+IDAA),ONE*(-ONE       ),ISCRDA(  7+IDAA))   
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  9+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  9+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),ZL         )                          
*FOX  CRKVE=XL ;                                                        *FOX
      CALL DACOP(XL         ,CRKVE      )                               
*FOX  CIKVE=ZL ;                                                        *FOX
      CALL DACOP(ZL         ,CIKVE      )                               
*FOX  CRABAMP4=ED(IX)/(ONE+DPDA) ;                                      *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      RSCRRI(  2+IDAA) = ED         (IX         )                       
      CALL DADIC(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  3+IDAA),CRABAMP4   )                          
             crabfreq=ek(ix)*c1e3
             crabpht4=crabph4(ix)
*FOX  Y(1)=Y(1) + 6*(1/6.)*(CRABAMP4)*                                  *FOX
*FOX  (CIKVE*CIKVE*CIKVE-(3*CIKVE*CRKVE*CRKVE))*C1M3*C1M3*              *FOX
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT4);               *FOX
      CALL DACMU(CIKVE      ,ONE*(3          ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),CRKVE      ,ISCRDA(  2+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),CRKVE      ,ISCRDA(  3+IDAA))         
      RSCRRI(  4+IDAA) = (1          ) / (6.         )                  
      CALL DAMUL(CIKVE      ,CIKVE      ,ISCRDA(  5+IDAA))              
      CALL DAMUL(ISCRDA(  5+IDAA),CIKVE      ,ISCRDA(  6+IDAA))         
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  7+IDAA))          
      CALL DACDI(ISCRDA(  7+IDAA),ONE*CLIGHT     ,ISCRDA(  8+IDAA))     
      CALL DACMU(ISCRDA(  8+IDAA),ONE*CRABFREQ   ,ISCRDA(  9+IDAA))     
      CALL DACMU(ISCRDA(  9+IDAA),ONE*(2D0        ),ISCRDA( 10+IDAA))   
      CALL DACMU(ISCRDA( 10+IDAA),ONE*PI         ,ISCRDA( 11+IDAA))     
      CALL DASUB(ISCRDA(  6+IDAA),ISCRDA(  3+IDAA),ISCRDA( 12+IDAA))    
      CALL DACAD(ISCRDA( 11+IDAA),ONE*CRABPHT4   ,ISCRDA( 13+IDAA))     
      CALL DACOP(Y          ((1)),ISCRDA( 14+IDAA))                     
      CALL DAFUN('COS ',ISCRDA( 13+IDAA),ISCRDA( 15+IDAA))              
      RSCRRI( 16+IDAA) = (6          ) * RSCRRI(  4+IDAA)               
      CALL DACMU(CRABAMP4   ,ONE*RSCRRI( 16+IDAA),ISCRDA( 17+IDAA))     
      CALL DAMUL(ISCRDA( 17+IDAA),ISCRDA( 12+IDAA),ISCRDA( 18+IDAA))    
      CALL DACMU(ISCRDA( 18+IDAA),ONE*C1M3       ,ISCRDA( 19+IDAA))     
      CALL DACMU(ISCRDA( 19+IDAA),ONE*C1M3       ,ISCRDA( 20+IDAA))     
      CALL DAMUL(ISCRDA( 20+IDAA),ISCRDA( 15+IDAA),ISCRDA( 21+IDAA))    
      CALL DAADD(ISCRDA( 14+IDAA),ISCRDA( 21+IDAA),ISCRDA( 22+IDAA))    
      CALL DACOP(ISCRDA( 22+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2) + 6*(1/6.)*(CRABAMP4)*                                  *FOX
*FOX  (3*CRKVE*CIKVE*CIKVE-CRKVE*CRKVE*CRKVE)*C1M3*C1M3*                *FOX
*FOX  COS(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI + CRABPHT4);               *FOX
      RSCRRI(  1+IDAA) = (1          ) / (6.         )                  
      CALL DACMU(CRKVE      ,ONE*(3          ),ISCRDA(  2+IDAA))        
      CALL DAMUL(ISCRDA(  2+IDAA),CIKVE      ,ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  3+IDAA),CIKVE      ,ISCRDA(  4+IDAA))         
      CALL DAMUL(CRKVE      ,CRKVE      ,ISCRDA(  5+IDAA))              
      CALL DAMUL(ISCRDA(  5+IDAA),CRKVE      ,ISCRDA(  6+IDAA))         
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA(  7+IDAA))          
      CALL DACDI(ISCRDA(  7+IDAA),ONE*CLIGHT     ,ISCRDA(  8+IDAA))     
      CALL DACMU(ISCRDA(  8+IDAA),ONE*CRABFREQ   ,ISCRDA(  9+IDAA))     
      CALL DACMU(ISCRDA(  9+IDAA),ONE*(2D0        ),ISCRDA( 10+IDAA))   
      CALL DACMU(ISCRDA( 10+IDAA),ONE*PI         ,ISCRDA( 11+IDAA))     
      CALL DASUB(ISCRDA(  4+IDAA),ISCRDA(  6+IDAA),ISCRDA( 12+IDAA))    
      CALL DACAD(ISCRDA( 11+IDAA),ONE*CRABPHT4   ,ISCRDA( 13+IDAA))     
      CALL DACOP(Y          ((2)),ISCRDA( 14+IDAA))                     
      CALL DAFUN('COS ',ISCRDA( 13+IDAA),ISCRDA( 15+IDAA))              
      RSCRRI( 16+IDAA) = (6          ) * RSCRRI(  1+IDAA)               
      CALL DACMU(CRABAMP4   ,ONE*RSCRRI( 16+IDAA),ISCRDA( 17+IDAA))     
      CALL DAMUL(ISCRDA( 17+IDAA),ISCRDA( 12+IDAA),ISCRDA( 18+IDAA))    
      CALL DACMU(ISCRDA( 18+IDAA),ONE*C1M3       ,ISCRDA( 19+IDAA))     
      CALL DACMU(ISCRDA( 19+IDAA),ONE*C1M3       ,ISCRDA( 20+IDAA))     
      CALL DAMUL(ISCRDA( 20+IDAA),ISCRDA( 15+IDAA),ISCRDA( 21+IDAA))    
      CALL DAADD(ISCRDA( 14+IDAA),ISCRDA( 21+IDAA),ISCRDA( 22+IDAA))    
      CALL DACOP(ISCRDA( 22+IDAA),Y          ((2)))                     
*FOX  DPDA1=DPDA1+6*(1/6.)*(CRABAMP4)*(CRKVE*CRKVE*CRKVE*CIKVE-         *FOX
*FOX  CIKVE*CIKVE*CIKVE*CRKVE)*                                         *FOX
*FOX  C1M9*(((CRABFREQ*2D0)*PI)/CLIGHT)*
*FOX  SIN(SIGMDA/C1E3/CLIGHT*CRABFREQ*2D0*PI+CRABPHT4) ;                *FOX
      RSCRRI(  1+IDAA) = CRABFREQ    * (2D0        )                    
      RSCRRI(  2+IDAA) = RSCRRI(  1+IDAA) * PI                          
      RSCRRI(  3+IDAA) = (1          ) / (6.         )                  
      CALL DAMUL(CRKVE      ,CRKVE      ,ISCRDA(  4+IDAA))              
      CALL DAMUL(ISCRDA(  4+IDAA),CRKVE      ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  5+IDAA),CIKVE      ,ISCRDA(  6+IDAA))         
      CALL DAMUL(CIKVE      ,CIKVE      ,ISCRDA(  7+IDAA))              
      CALL DAMUL(ISCRDA(  7+IDAA),CIKVE      ,ISCRDA(  8+IDAA))         
      CALL DAMUL(ISCRDA(  8+IDAA),CRKVE      ,ISCRDA(  9+IDAA))         
      RSCRRI( 10+IDAA) = RSCRRI(  2+IDAA) / CLIGHT                      
      CALL DACDI(SIGMDA     ,ONE*C1E3       ,ISCRDA( 11+IDAA))          
      CALL DACDI(ISCRDA( 11+IDAA),ONE*CLIGHT     ,ISCRDA( 12+IDAA))     
      CALL DACMU(ISCRDA( 12+IDAA),ONE*CRABFREQ   ,ISCRDA( 13+IDAA))     
      CALL DACMU(ISCRDA( 13+IDAA),ONE*(2D0        ),ISCRDA( 14+IDAA))   
      CALL DACMU(ISCRDA( 14+IDAA),ONE*PI         ,ISCRDA( 15+IDAA))     
      CALL DASUB(ISCRDA(  6+IDAA),ISCRDA(  9+IDAA),ISCRDA( 16+IDAA))    
      CALL DACAD(ISCRDA( 15+IDAA),ONE*CRABPHT4   ,ISCRDA( 17+IDAA))     
      CALL DAFUN('SIN ',ISCRDA( 17+IDAA),ISCRDA( 18+IDAA))              
      RSCRRI( 19+IDAA) = (6          ) * RSCRRI(  3+IDAA)               
      CALL DACMU(CRABAMP4   ,ONE*RSCRRI( 19+IDAA),ISCRDA( 20+IDAA))     
      CALL DAMUL(ISCRDA( 20+IDAA),ISCRDA( 16+IDAA),ISCRDA( 21+IDAA))    
      CALL DACMU(ISCRDA( 21+IDAA),ONE*C1M9       ,ISCRDA( 22+IDAA))     
      CALL DACMU(ISCRDA( 22+IDAA),ONE*RSCRRI( 10+IDAA),ISCRDA( 23+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA( 23+IDAA),ISCRDA( 18+IDAA),ISCRDA( 24+IDAA))    
      CALL DAADD(DPDA1      ,ISCRDA( 24+IDAA),ISCRDA( 25+IDAA))         
      CALL DACOP(ISCRDA( 25+IDAA),DPDA1      )                          
*FOX  EJF0=EJF1 ;                                                       *FOX
      CALL DACOP(EJF1       ,EJF0       )                               
*FOX  DPDA=DPDA1*C1M3 ;                                                 *FOX
      CALL DACMU(DPDA1      ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA       )                          
*FOX  EJF1=E0F*(ONE+DPDA) ;                                             *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),EJF1       )                          
*FOX  EJ1=SQRT(EJF1*EJF1+PMA*PMA) ;                                     *FOX
      CALL DAMUL(EJF1       ,EJF1       ,ISCRDA(  1+IDAA))              
      RSCRRI(  2+IDAA) = PMA         * PMA                              
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DAFUN('SQRT',ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))              
      CALL DACOP(ISCRDA(  4+IDAA),EJ1        )                          
*FOX  Y(1)=EJF0/EJF1*Y(1) ;                                             *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((1)))                     
*FOX  Y(2)=EJF0/EJF1*Y(2) ;                                             *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DADIV(EJF0       ,EJF1       ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((2)))                     
          endif
        if(kzz.eq.22) then
          irrtr=imtr(ix)
*FOX  SIGMDA=SIGMDA+COTR(IRRTR,5)+RRTR(IRRTR,5,1)*X(1)+                 *FOX
*FOX  RRTR(IRRTR,5,2)*Y(1)+RRTR(IRRTR,5,3)*X(2)+RRTR(IRRTR,5,4)*Y(2) ;  *FOX
      RSCRRI(  1+IDAA) = COTR       (IRRTR      ,(5))                   
      RSCRRI(  2+IDAA) = RRTR       (IRRTR      ,(5),(1))               
      CALL DACOP(X          ((1)),ISCRDA(  3+IDAA))                     
      RSCRRI(  4+IDAA) = RRTR       (IRRTR      ,(5),(2))               
      CALL DACOP(Y          ((1)),ISCRDA(  5+IDAA))                     
      RSCRRI(  6+IDAA) = RRTR       (IRRTR      ,(5),(3))               
      CALL DACOP(X          ((2)),ISCRDA(  7+IDAA))                     
      RSCRRI(  8+IDAA) = RRTR       (IRRTR      ,(5),(4))               
      CALL DACOP(Y          ((2)),ISCRDA(  9+IDAA))                     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA( 10+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  5+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA( 11+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA( 12+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  9+IDAA),ONE*RSCRRI(  8+IDAA),ISCRDA( 13+IDAA))
     *                                                                  
      CALL DACAD(SIGMDA     ,ONE*RSCRRI(  1+IDAA),ISCRDA( 14+IDAA))     
      CALL DAADD(ISCRDA( 14+IDAA),ISCRDA( 10+IDAA),ISCRDA( 15+IDAA))    
      CALL DAADD(ISCRDA( 15+IDAA),ISCRDA( 11+IDAA),ISCRDA( 16+IDAA))    
      CALL DAADD(ISCRDA( 16+IDAA),ISCRDA( 12+IDAA),ISCRDA( 17+IDAA))    
      CALL DAADD(ISCRDA( 17+IDAA),ISCRDA( 13+IDAA),ISCRDA( 18+IDAA))    
      CALL DACOP(ISCRDA( 18+IDAA),SIGMDA     )                          
*FOX  PUX=X(1) ;                                                        *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUX        )                          
*FOX  PUZ=Y(1) ;                                                        *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUZ        )                          
*FOX  X(1)=COTR(IRRTR,1)+RRTR(IRRTR,1,1)*PUX+RRTR(IRRTR,1,2)*PUZ+       *FOX
*FOX  RRTR(IRRTR,1,6)*DPDA1 ;                                           *FOX
      RSCRRI(  1+IDAA) = COTR       (IRRTR      ,(1))                   
      RSCRRI(  2+IDAA) = RRTR       (IRRTR      ,(1),(1))               
      RSCRRI(  3+IDAA) = RRTR       (IRRTR      ,(1),(2))               
      RSCRRI(  4+IDAA) = RRTR       (IRRTR      ,(1),(6))               
      CALL DACMU(PUX        ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DACMU(PUZ        ,ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))     
      CALL DACMU(DPDA1      ,ONE*RSCRRI(  4+IDAA),ISCRDA(  7+IDAA))     
      CALL DACAD(ISCRDA(  5+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DAADD(ISCRDA(  9+IDAA),ISCRDA(  7+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),X          ((1)))                     
*FOX  Y(1)=COTR(IRRTR,2)+RRTR(IRRTR,2,1)*PUX+RRTR(IRRTR,2,2)*PUZ+       *FOX
*FOX  RRTR(IRRTR,2,6)*DPDA1 ;                                           *FOX
      RSCRRI(  1+IDAA) = COTR       (IRRTR      ,(2))                   
      RSCRRI(  2+IDAA) = RRTR       (IRRTR      ,(2),(1))               
      RSCRRI(  3+IDAA) = RRTR       (IRRTR      ,(2),(2))               
      RSCRRI(  4+IDAA) = RRTR       (IRRTR      ,(2),(6))               
      CALL DACMU(PUX        ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DACMU(PUZ        ,ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))     
      CALL DACMU(DPDA1      ,ONE*RSCRRI(  4+IDAA),ISCRDA(  7+IDAA))     
      CALL DACAD(ISCRDA(  5+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DAADD(ISCRDA(  9+IDAA),ISCRDA(  7+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),Y          ((1)))                     
*FOX  PUX=X(2) ;                                                        *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUX        )                          
*FOX  PUZ=Y(2) ;                                                        *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),PUZ        )                          
*FOX  X(2)=COTR(IRRTR,3)+RRTR(IRRTR,3,3)*PUX+RRTR(IRRTR,3,4)*PUZ+       *FOX
*FOX  RRTR(IRRTR,3,6)*DPDA1 ;                                           *FOX
      RSCRRI(  1+IDAA) = COTR       (IRRTR      ,(3))                   
      RSCRRI(  2+IDAA) = RRTR       (IRRTR      ,(3),(3))               
      RSCRRI(  3+IDAA) = RRTR       (IRRTR      ,(3),(4))               
      RSCRRI(  4+IDAA) = RRTR       (IRRTR      ,(3),(6))               
      CALL DACMU(PUX        ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DACMU(PUZ        ,ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))     
      CALL DACMU(DPDA1      ,ONE*RSCRRI(  4+IDAA),ISCRDA(  7+IDAA))     
      CALL DACAD(ISCRDA(  5+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DAADD(ISCRDA(  9+IDAA),ISCRDA(  7+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),X          ((2)))                     
*FOX  Y(2)=COTR(IRRTR,4)+RRTR(IRRTR,4,3)*PUX+RRTR(IRRTR,4,4)*PUZ+       *FOX
*FOX  RRTR(IRRTR,4,6)*DPDA1 ;                                           *FOX
      RSCRRI(  1+IDAA) = COTR       (IRRTR      ,(4))                   
      RSCRRI(  2+IDAA) = RRTR       (IRRTR      ,(4),(3))               
      RSCRRI(  3+IDAA) = RRTR       (IRRTR      ,(4),(4))               
      RSCRRI(  4+IDAA) = RRTR       (IRRTR      ,(4),(6))               
      CALL DACMU(PUX        ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DACMU(PUZ        ,ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))     
      CALL DACMU(DPDA1      ,ONE*RSCRRI(  4+IDAA),ISCRDA(  7+IDAA))     
      CALL DACAD(ISCRDA(  5+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DAADD(ISCRDA(  9+IDAA),ISCRDA(  7+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),Y          ((2)))                     
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
*FOX  EKK=(SMIDA(IPCH)*RATIOE(IX)+SMIZF(I))/(ONE+DPDA) ;                *FOX
      CALL DACOP(SMIDA      (IPCH       ),ISCRDA(  1+IDAA))             
      RSCRRI(  2+IDAA) = RATIOE     (IX         )                       
      RSCRRI(  3+IDAA) = SMIZF      (I          )                       
      CALL DACMU(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      CALL DACAD(ISCRDA(  4+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  6+IDAA))          
      CALL DADIV(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  7+IDAA))    
      CALL DACOP(ISCRDA(  7+IDAA),EKK        )                          
        else
*FOX  EKK=SMI(I)/(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      RSCRRI(  2+IDAA) = SMI        (I          )                       
      CALL DADIC(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  3+IDAA),EKK        )                          
        endif
        xs=xsi(i)
        zs=zsi(i)
*FOX  XL=(X(1)-XS)*TILTC(I)+(X(2)-ZS)*TILTS(I) ;                        *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*XS         ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  2+IDAA),ONE*ZS         ,ISCRDA(  4+IDAA))     
      RSCRRI(  5+IDAA) = TILTC      (I          )                       
      RSCRRI(  6+IDAA) = TILTS      (I          )                       
      CALL DACMU(ISCRDA(  3+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),XL         )                          
*FOX  ZL=-(X(1)-XS)*TILTS(I)+(X(2)-ZS)*TILTC(I) ;                       *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*XS         ,ISCRDA(  3+IDAA))     
      CALL DACSU(ISCRDA(  2+IDAA),ONE*ZS         ,ISCRDA(  4+IDAA))     
      RSCRRI(  5+IDAA) = TILTS      (I          )                       
      RSCRRI(  6+IDAA) = TILTC      (I          )                       
      CALL DACMU(ISCRDA(  3+IDAA),ONE*(-ONE       ),ISCRDA(  7+IDAA))   
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  9+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  9+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),ZL         )                          
*FOX  CRKVE=XL ;                                                        *FOX
      CALL DACOP(XL         ,CRKVE      )                               
*FOX  CIKVE=ZL ;                                                        *FOX
      CALL DACOP(ZL         ,CIKVE      )                               
        if(kzz.lt.0) goto 320
        goto(90,100,110,120,130,140,150,160,170,180,190,440,440,440,    &
     &       440,440,440,440,440,440,440,440,440,185,186),kzz
        goto 440
!--HORIZONTAL DIPOLE
   90   continue
*FOX  EKK=EKK*C1E3 ;                                                    *FOX
      CALL DACMU(EKK        ,ONE*C1E3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  Y(1)=Y(1)+EKK*TILTC(I) ;                                          *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = TILTC      (I          )                       
      CALL DACMU(EKK        ,ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))     
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*TILTS(I) ;                                          *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(EKK        ,ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))     
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),Y          ((2)))                     
        goto 440
!--NORMAL QUADRUPOLE
  100   continue
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;                  *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      RSCRRI(  3+IDAA) = (-ONE       ) * RSCRRI(  1+IDAA)               
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DAADD(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  7+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))         
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((2)))                     
        goto 440
!---NORMAL SEXTUPOLE
  110   continue
*FOX  EKK=EKK*C1M3 ;                                                    *FOX
      CALL DACMU(EKK        ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;                  *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      RSCRRI(  3+IDAA) = (-ONE       ) * RSCRRI(  1+IDAA)               
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DAADD(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  7+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))         
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((2)))                     
        goto 440
!--NORMAL OCTUPOLE
  120   continue
*FOX  EKK=EKK*C1M6 ;                                                    *FOX
      CALL DACMU(EKK        ,ONE*C1M6       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;                  *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      RSCRRI(  3+IDAA) = (-ONE       ) * RSCRRI(  1+IDAA)               
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DAADD(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  7+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))         
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((2)))                     
        goto 440
!--NORMAL DECAPOLE
  130   continue
*FOX  EKK=EKK*C1M9 ;                                                    *FOX
      CALL DACMU(EKK        ,ONE*C1M9       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;                  *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      RSCRRI(  3+IDAA) = (-ONE       ) * RSCRRI(  1+IDAA)               
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DAADD(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  7+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))         
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((2)))                     
        goto 440
!---NORMAL DODECAPOL
  140   continue
*FOX  EKK=EKK*C1M12 ;                                                   *FOX
      CALL DACMU(EKK        ,ONE*C1M12      ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;                  *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      RSCRRI(  3+IDAA) = (-ONE       ) * RSCRRI(  1+IDAA)               
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DAADD(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  7+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))         
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((2)))                     
        goto 440
!---NORMAL 14-POL
  150   continue
*FOX  EKK=EKK*C1M15 ;                                                   *FOX
      CALL DACMU(EKK        ,ONE*C1M15      ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;                  *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      RSCRRI(  3+IDAA) = (-ONE       ) * RSCRRI(  1+IDAA)               
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DAADD(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  7+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))         
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((2)))                     
        goto 440
!---NORMAL 16-POL
  160   continue
*FOX  EKK=EKK*C1M18 ;                                                   *FOX
      CALL DACMU(EKK        ,ONE*C1M18      ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;                  *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      RSCRRI(  3+IDAA) = (-ONE       ) * RSCRRI(  1+IDAA)               
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DAADD(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  7+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))         
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((2)))                     
        goto 440
!---NORMAL 18-POL
  170   continue
*FOX  EKK=EKK*C1M21 ;                                                   *FOX
      CALL DACMU(EKK        ,ONE*C1M21      ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;                  *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      RSCRRI(  3+IDAA) = (-ONE       ) * RSCRRI(  1+IDAA)               
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DAADD(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  7+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))         
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((2)))                     
        goto 440
!---NORMAL 20-POL
  180   continue
*FOX  EKK=EKK*C1M24 ;                                                   *FOX
      CALL DACMU(EKK        ,ONE*C1M24      ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(-TILTC(I)*CIKVE+TILTS(I)*CRKVE) ;                  *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      RSCRRI(  3+IDAA) = (-ONE       ) * RSCRRI(  1+IDAA)               
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DAADD(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  7+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))         
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((2)))                     
        goto 440
!--DIPEDGE ELEMENT
  185   continue
*FOX  Y(1)=Y(1)+(ED(IX)*TILTC(I)*CRKVE-EK(IX)*TILTS(I)*CIKVE)/          *FOX
*FOX  (ONE+DPDA) ;                                                      *FOX
      RSCRRI(  1+IDAA) = ED         (IX         )                       
      RSCRRI(  2+IDAA) = TILTC      (I          )                       
      RSCRRI(  3+IDAA) = EK         (IX         )                       
      RSCRRI(  4+IDAA) = TILTS      (I          )                       
      RSCRRI(  5+IDAA) = RSCRRI(  1+IDAA) * RSCRRI(  2+IDAA)            
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  5+IDAA),ISCRDA(  6+IDAA))     
      RSCRRI(  7+IDAA) = RSCRRI(  3+IDAA) * RSCRRI(  4+IDAA)            
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  7+IDAA),ISCRDA(  8+IDAA))     
      CALL DASUB(ISCRDA(  6+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 10+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA( 11+IDAA))                     
      CALL DADIV(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 12+IDAA))    
      CALL DAADD(ISCRDA( 11+IDAA),ISCRDA( 12+IDAA),ISCRDA( 13+IDAA))    
      CALL DACOP(ISCRDA( 13+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+(EK(IX)*TILTC(I)*CIKVE+ED(IX)*TILTS(I)*CRKVE)/          *FOX
*FOX  (ONE+DPDA) ;                                                      *FOX
      RSCRRI(  1+IDAA) = EK         (IX         )                       
      RSCRRI(  2+IDAA) = TILTC      (I          )                       
      RSCRRI(  3+IDAA) = ED         (IX         )                       
      RSCRRI(  4+IDAA) = TILTS      (I          )                       
      RSCRRI(  5+IDAA) = RSCRRI(  1+IDAA) * RSCRRI(  2+IDAA)            
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  5+IDAA),ISCRDA(  6+IDAA))     
      RSCRRI(  7+IDAA) = RSCRRI(  3+IDAA) * RSCRRI(  4+IDAA)            
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  7+IDAA),ISCRDA(  8+IDAA))     
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 10+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA( 11+IDAA))                     
      CALL DADIV(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 12+IDAA))    
      CALL DAADD(ISCRDA( 11+IDAA),ISCRDA( 12+IDAA),ISCRDA( 13+IDAA))    
      CALL DACOP(ISCRDA( 13+IDAA),Y          ((2)))                     
        goto 440
!--solenoid
  186   continue
*FOX  Y(1)=Y(1)-X(2)*ED(IX) ;                                           *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      RSCRRI(  3+IDAA) = ED         (IX         )                       
      CALL DACMU(ISCRDA(  2+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(ISCRDA(  5+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+X(1)*ED(IX) ;                                           *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((1)),ISCRDA(  2+IDAA))                     
      RSCRRI(  3+IDAA) = ED         (IX         )                       
      CALL DACMU(ISCRDA(  2+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(ISCRDA(  5+IDAA),Y          ((2)))                     
*FOX  CRKVE=Y(1)-X(1)*ED(IX)*EK(IX) ;                                   *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((1)),ISCRDA(  2+IDAA))                     
      RSCRRI(  3+IDAA) = ED         (IX         )                       
      RSCRRI(  4+IDAA) = EK         (IX         )                       
      CALL DACMU(ISCRDA(  2+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  5+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA(  6+IDAA))
     *                                                                  
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  6+IDAA),ISCRDA(  7+IDAA))    
      CALL DACOP(ISCRDA(  7+IDAA),CRKVE      )                          
*FOX  CIKVE=Y(2)-X(2)*ED(IX)*EK(IX) ;                                   *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      RSCRRI(  3+IDAA) = ED         (IX         )                       
      RSCRRI(  4+IDAA) = EK         (IX         )                       
      CALL DACMU(ISCRDA(  2+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  5+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA(  6+IDAA))
     *                                                                  
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  6+IDAA),ISCRDA(  7+IDAA))    
      CALL DACOP(ISCRDA(  7+IDAA),CIKVE      )                          
*FOX  Y(1)=CRKVE*COS(EK(IX))+CIKVE*SIN(EK(IX)) ;                        *FOX
      RSCRRI(  1+IDAA) = EK         (IX         )                       
      RSCRRI(  2+IDAA) = EK         (IX         )                       
      RSCRRI(  3+IDAA) = COS (RSCRRI(  1+IDAA))                         
      RSCRRI(  4+IDAA) = SIN (RSCRRI(  2+IDAA))                         
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  3+IDAA),ISCRDA(  5+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  4+IDAA),ISCRDA(  6+IDAA))     
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  7+IDAA))    
      CALL DACOP(ISCRDA(  7+IDAA),Y          ((1)))                     
*FOX  Y(2)=-CRKVE*SIN(EK(IX))+CIKVE*COS(EK(IX)) ;                       *FOX
      RSCRRI(  1+IDAA) = EK         (IX         )                       
      RSCRRI(  2+IDAA) = EK         (IX         )                       
      RSCRRI(  3+IDAA) = SIN (RSCRRI(  1+IDAA))                         
      RSCRRI(  4+IDAA) = COS (RSCRRI(  2+IDAA))                         
      CALL DACMU(CRKVE      ,ONE*(-ONE       ),ISCRDA(  5+IDAA))        
      CALL DACMU(ISCRDA(  5+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))
     *                                                                  
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  4+IDAA),ISCRDA(  7+IDAA))     
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((2)))                     
*FOX  CRKVE=X(1)*COS(EK(IX))+X(2)*SIN(EK(IX)) ;                         *FOX
      RSCRRI(  1+IDAA) = EK         (IX         )                       
      RSCRRI(  2+IDAA) = EK         (IX         )                       
      CALL DACOP(X          ((1)),ISCRDA(  3+IDAA))                     
      RSCRRI(  4+IDAA) = COS (RSCRRI(  1+IDAA))                         
      CALL DACOP(X          ((2)),ISCRDA(  5+IDAA))                     
      RSCRRI(  6+IDAA) = SIN (RSCRRI(  2+IDAA))                         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  5+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),CRKVE      )                          
*FOX  CIKVE=-X(1)*SIN(EK(IX))+X(2)*COS(EK(IX)) ;                        *FOX
      RSCRRI(  1+IDAA) = EK         (IX         )                       
      RSCRRI(  2+IDAA) = EK         (IX         )                       
      CALL DACOP(X          ((1)),ISCRDA(  3+IDAA))                     
      RSCRRI(  4+IDAA) = SIN (RSCRRI(  1+IDAA))                         
      CALL DACOP(X          ((2)),ISCRDA(  5+IDAA))                     
      RSCRRI(  6+IDAA) = COS (RSCRRI(  2+IDAA))                         
      CALL DACMU(ISCRDA(  3+IDAA),ONE*(-ONE       ),ISCRDA(  7+IDAA))   
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  5+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  9+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  9+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),CIKVE      )                          
*FOX  X(1)=CRKVE ;                                                      *FOX
      CALL DACOP(CRKVE      ,X          ((1)))                          
*FOX  X(2)=CIKVE ;                                                      *FOX
      CALL DACOP(CIKVE      ,X          ((2)))                          
*FOX  Y(1)=Y(1)+X(2)*ED(IX) ;                                           *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((2)),ISCRDA(  2+IDAA))                     
      RSCRRI(  3+IDAA) = ED         (IX         )                       
      CALL DACMU(ISCRDA(  2+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(ISCRDA(  5+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)-X(1)*ED(IX) ;                                           *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(X          ((1)),ISCRDA(  2+IDAA))                     
      RSCRRI(  3+IDAA) = ED         (IX         )                       
      CALL DACMU(ISCRDA(  2+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(ISCRDA(  5+IDAA),Y          ((2)))                     
        goto 440
  190   r0=ek(ix)
        nmz=nmu(ix)
          if(abs(dki(ix,1)).gt.pieni) then
            if(abs(dki(ix,3)).gt.pieni) then
*FOX  DKIP=DKI(IX,1)/(ONE+DPDA) ;                                       *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      RSCRRI(  2+IDAA) = DKI        (IX         ,(1))                   
      CALL DADIC(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  3+IDAA),DKIP       )                          
*FOX  Y(1)=Y(1)-(DKI(IX,1)/DKI(IX,3)*XL+DPDA*C1E3)*                     *FOX
*FOX  TILTC(I)*DKIP                                                     *FOX
*FOX  +C1E3*DKI(IX,1)/(ONE+DPDA)*(ONE-TILTC(I)) ;                       *FOX
      RSCRRI(  1+IDAA) = DKI        (IX         ,(1))                   
      RSCRRI(  2+IDAA) = DKI        (IX         ,(3))                   
      RSCRRI(  3+IDAA) = TILTC      (I          )                       
      RSCRRI(  4+IDAA) = RSCRRI(  1+IDAA) / RSCRRI(  2+IDAA)            
      CALL DACMU(XL         ,ONE*RSCRRI(  4+IDAA),ISCRDA(  5+IDAA))     
      CALL DACMU(DPDA       ,ONE*C1E3       ,ISCRDA(  6+IDAA))          
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  7+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  8+IDAA))          
      RSCRRI(  9+IDAA) = ONE         - RSCRRI(  3+IDAA)                 
      CALL DACOP(Y          ((1)),ISCRDA( 10+IDAA))                     
      RSCRRI( 11+IDAA) = TILTC      (I          )                       
      RSCRRI( 12+IDAA) = DKI        (IX         ,(1))                   
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RSCRRI( 11+IDAA),ISCRDA( 13+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA( 13+IDAA),DKIP       ,ISCRDA( 14+IDAA))         
      RSCRRI( 15+IDAA) = C1E3        * RSCRRI( 12+IDAA)                 
      CALL DADIC(ISCRDA(  8+IDAA),ONE*RSCRRI( 15+IDAA),ISCRDA( 16+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 16+IDAA),ONE*RSCRRI(  9+IDAA),ISCRDA( 17+IDAA))
     *                                                                  
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA( 14+IDAA),ISCRDA( 18+IDAA))    
      CALL DAADD(ISCRDA( 18+IDAA),ISCRDA( 17+IDAA),ISCRDA( 19+IDAA))    
      CALL DACOP(ISCRDA( 19+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)-(DKI(IX,1)/DKI(IX,3)*XL+DPDA*C1E3)*                     *FOX
*FOX  TILTS(I)*DKIP                                                     *FOX
*FOX  +C1E3*DKI(IX,1)/(ONE+DPDA)*TILTS(I) ;                             *FOX
      RSCRRI(  1+IDAA) = DKI        (IX         ,(1))                   
      RSCRRI(  2+IDAA) = DKI        (IX         ,(3))                   
      RSCRRI(  3+IDAA) = RSCRRI(  1+IDAA) / RSCRRI(  2+IDAA)            
      CALL DACMU(XL         ,ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))     
      CALL DACMU(DPDA       ,ONE*C1E3       ,ISCRDA(  5+IDAA))          
      CALL DAADD(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  7+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  8+IDAA))                     
      RSCRRI(  9+IDAA) = TILTS      (I          )                       
      RSCRRI( 10+IDAA) = DKI        (IX         ,(1))                   
      RSCRRI( 11+IDAA) = TILTS      (I          )                       
      CALL DACMU(ISCRDA(  6+IDAA),ONE*RSCRRI(  9+IDAA),ISCRDA( 12+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA( 12+IDAA),DKIP       ,ISCRDA( 13+IDAA))         
      RSCRRI( 14+IDAA) = C1E3        * RSCRRI( 10+IDAA)                 
      CALL DADIC(ISCRDA(  7+IDAA),ONE*RSCRRI( 14+IDAA),ISCRDA( 15+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 15+IDAA),ONE*RSCRRI( 11+IDAA),ISCRDA( 16+IDAA))
     *                                                                  
      CALL DASUB(ISCRDA(  8+IDAA),ISCRDA( 13+IDAA),ISCRDA( 17+IDAA))    
      CALL DAADD(ISCRDA( 17+IDAA),ISCRDA( 16+IDAA),ISCRDA( 18+IDAA))    
      CALL DACOP(ISCRDA( 18+IDAA),Y          ((2)))                     
            else
*FOX  Y(1)=Y(1)-DKI(IX,1)*DPDA*C1E3/(ONE+DPDA)*TILTC(I)                 *FOX
*FOX  +C1E3*DKI(IX,1)/(ONE+DPDA)*(ONE-TILTC(I)) ;                       *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  3+IDAA))          
      RSCRRI(  4+IDAA) = ONE         - RSCRRI(  1+IDAA)                 
      CALL DACOP(Y          ((1)),ISCRDA(  5+IDAA))                     
      RSCRRI(  6+IDAA) = DKI        (IX         ,(1))                   
      RSCRRI(  7+IDAA) = TILTC      (I          )                       
      RSCRRI(  8+IDAA) = DKI        (IX         ,(1))                   
      CALL DACMU(DPDA       ,ONE*RSCRRI(  6+IDAA),ISCRDA(  9+IDAA))     
      CALL DACMU(ISCRDA(  9+IDAA),ONE*C1E3       ,ISCRDA( 10+IDAA))     
      CALL DADIV(ISCRDA( 10+IDAA),ISCRDA(  2+IDAA),ISCRDA( 11+IDAA))    
      CALL DACMU(ISCRDA( 11+IDAA),ONE*RSCRRI(  7+IDAA),ISCRDA( 12+IDAA))
     *                                                                  
      RSCRRI( 13+IDAA) = C1E3        * RSCRRI(  8+IDAA)                 
      CALL DADIC(ISCRDA(  3+IDAA),ONE*RSCRRI( 13+IDAA),ISCRDA( 14+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 14+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA( 15+IDAA))
     *                                                                  
      CALL DASUB(ISCRDA(  5+IDAA),ISCRDA( 12+IDAA),ISCRDA( 16+IDAA))    
      CALL DAADD(ISCRDA( 16+IDAA),ISCRDA( 15+IDAA),ISCRDA( 17+IDAA))    
      CALL DACOP(ISCRDA( 17+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)-DKI(IX,1)*DPDA*C1E3/(ONE+DPDA)*TILTS(I)                 *FOX
*FOX  +C1E3*DKI(IX,1)/(ONE+DPDA)*TILTS(I) ;                             *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  3+IDAA))                     
      RSCRRI(  4+IDAA) = DKI        (IX         ,(1))                   
      RSCRRI(  5+IDAA) = TILTS      (I          )                       
      RSCRRI(  6+IDAA) = DKI        (IX         ,(1))                   
      RSCRRI(  7+IDAA) = TILTS      (I          )                       
      CALL DACMU(DPDA       ,ONE*RSCRRI(  4+IDAA),ISCRDA(  8+IDAA))     
      CALL DACMU(ISCRDA(  8+IDAA),ONE*C1E3       ,ISCRDA(  9+IDAA))     
      CALL DADIV(ISCRDA(  9+IDAA),ISCRDA(  1+IDAA),ISCRDA( 10+IDAA))    
      CALL DACMU(ISCRDA( 10+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA( 11+IDAA))
     *                                                                  
      RSCRRI( 12+IDAA) = C1E3        * RSCRRI(  6+IDAA)                 
      CALL DADIC(ISCRDA(  2+IDAA),ONE*RSCRRI( 12+IDAA),ISCRDA( 13+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 13+IDAA),ONE*RSCRRI(  7+IDAA),ISCRDA( 14+IDAA))
     *                                                                  
      CALL DASUB(ISCRDA(  3+IDAA),ISCRDA( 11+IDAA),ISCRDA( 15+IDAA))    
      CALL DAADD(ISCRDA( 15+IDAA),ISCRDA( 14+IDAA),ISCRDA( 16+IDAA))    
      CALL DACOP(ISCRDA( 16+IDAA),Y          ((2)))                     
            endif
            if(idp.eq.1.and.iabs(ition).eq.1) then
*FOX  SIGMDA=SIGMDA+RV*DKI(IX,1)*XL ;                                   *FOX
      RSCRRI(  1+IDAA) = DKI        (IX         ,(1))                   
      CALL DACMU(RV         ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DAMUL(ISCRDA(  2+IDAA),XL         ,ISCRDA(  3+IDAA))         
      CALL DAADD(SIGMDA     ,ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))         
      CALL DACOP(ISCRDA(  4+IDAA),SIGMDA     )                          
            endif
          endif
          if(abs(dki(ix,2)).gt.pieni) then
            if(abs(dki(ix,3)).gt.pieni) then
*FOX  DKIP=DKI(IX,2)/(ONE+DPDA) ;                                       *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      RSCRRI(  2+IDAA) = DKI        (IX         ,(2))                   
      CALL DADIC(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  3+IDAA),DKIP       )                          
*FOX  Y(1)=Y(1)+(DKI(IX,2)/DKI(IX,3)*ZL-DPDA*C1E3)*                     *FOX
*FOX  TILTS(I)*DKIP                                                     *FOX
*FOX  +C1E3*DKI(IX,2)/(ONE+DPDA)*TILTS(I) ;                             *FOX
      RSCRRI(  1+IDAA) = DKI        (IX         ,(2))                   
      RSCRRI(  2+IDAA) = DKI        (IX         ,(3))                   
      RSCRRI(  3+IDAA) = RSCRRI(  1+IDAA) / RSCRRI(  2+IDAA)            
      CALL DACMU(ZL         ,ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))     
      CALL DACMU(DPDA       ,ONE*C1E3       ,ISCRDA(  5+IDAA))          
      CALL DASUB(ISCRDA(  4+IDAA),ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  7+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  8+IDAA))                     
      RSCRRI(  9+IDAA) = TILTS      (I          )                       
      RSCRRI( 10+IDAA) = DKI        (IX         ,(2))                   
      RSCRRI( 11+IDAA) = TILTS      (I          )                       
      CALL DACMU(ISCRDA(  6+IDAA),ONE*RSCRRI(  9+IDAA),ISCRDA( 12+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA( 12+IDAA),DKIP       ,ISCRDA( 13+IDAA))         
      RSCRRI( 14+IDAA) = C1E3        * RSCRRI( 10+IDAA)                 
      CALL DADIC(ISCRDA(  7+IDAA),ONE*RSCRRI( 14+IDAA),ISCRDA( 15+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 15+IDAA),ONE*RSCRRI( 11+IDAA),ISCRDA( 16+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA( 13+IDAA),ISCRDA( 17+IDAA))    
      CALL DAADD(ISCRDA( 17+IDAA),ISCRDA( 16+IDAA),ISCRDA( 18+IDAA))    
      CALL DACOP(ISCRDA( 18+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)-(DKI(IX,2)/DKI(IX,3)*ZL-DPDA*C1E3)*                     *FOX
*FOX  TILTC(I)*DKIP                                                     *FOX
*FOX  -C1E3*DKI(IX,2)/(ONE+DPDA)*(ONE-TILTC(I)) ;                       *FOX
      RSCRRI(  1+IDAA) = DKI        (IX         ,(2))                   
      RSCRRI(  2+IDAA) = DKI        (IX         ,(3))                   
      RSCRRI(  3+IDAA) = TILTC      (I          )                       
      RSCRRI(  4+IDAA) = RSCRRI(  1+IDAA) / RSCRRI(  2+IDAA)            
      CALL DACMU(ZL         ,ONE*RSCRRI(  4+IDAA),ISCRDA(  5+IDAA))     
      CALL DACMU(DPDA       ,ONE*C1E3       ,ISCRDA(  6+IDAA))          
      CALL DASUB(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  7+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  8+IDAA))          
      RSCRRI(  9+IDAA) = ONE         - RSCRRI(  3+IDAA)                 
      CALL DACOP(Y          ((2)),ISCRDA( 10+IDAA))                     
      RSCRRI( 11+IDAA) = TILTC      (I          )                       
      RSCRRI( 12+IDAA) = DKI        (IX         ,(2))                   
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RSCRRI( 11+IDAA),ISCRDA( 13+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA( 13+IDAA),DKIP       ,ISCRDA( 14+IDAA))         
      RSCRRI( 15+IDAA) = C1E3        * RSCRRI( 12+IDAA)                 
      CALL DADIC(ISCRDA(  8+IDAA),ONE*RSCRRI( 15+IDAA),ISCRDA( 16+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 16+IDAA),ONE*RSCRRI(  9+IDAA),ISCRDA( 17+IDAA))
     *                                                                  
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA( 14+IDAA),ISCRDA( 18+IDAA))    
      CALL DASUB(ISCRDA( 18+IDAA),ISCRDA( 17+IDAA),ISCRDA( 19+IDAA))    
      CALL DACOP(ISCRDA( 19+IDAA),Y          ((2)))                     
            else
*FOX  Y(1)=Y(1)-DKI(IX,2)*DPDA*C1E3/(ONE+DPDA)*TILTS(I)                 *FOX
*FOX  +C1E3*DKI(IX,2)/(ONE+DPDA)*TILTS(I) ;                             *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  3+IDAA))                     
      RSCRRI(  4+IDAA) = DKI        (IX         ,(2))                   
      RSCRRI(  5+IDAA) = TILTS      (I          )                       
      RSCRRI(  6+IDAA) = DKI        (IX         ,(2))                   
      RSCRRI(  7+IDAA) = TILTS      (I          )                       
      CALL DACMU(DPDA       ,ONE*RSCRRI(  4+IDAA),ISCRDA(  8+IDAA))     
      CALL DACMU(ISCRDA(  8+IDAA),ONE*C1E3       ,ISCRDA(  9+IDAA))     
      CALL DADIV(ISCRDA(  9+IDAA),ISCRDA(  1+IDAA),ISCRDA( 10+IDAA))    
      CALL DACMU(ISCRDA( 10+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA( 11+IDAA))
     *                                                                  
      RSCRRI( 12+IDAA) = C1E3        * RSCRRI(  6+IDAA)                 
      CALL DADIC(ISCRDA(  2+IDAA),ONE*RSCRRI( 12+IDAA),ISCRDA( 13+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 13+IDAA),ONE*RSCRRI(  7+IDAA),ISCRDA( 14+IDAA))
     *                                                                  
      CALL DASUB(ISCRDA(  3+IDAA),ISCRDA( 11+IDAA),ISCRDA( 15+IDAA))    
      CALL DAADD(ISCRDA( 15+IDAA),ISCRDA( 14+IDAA),ISCRDA( 16+IDAA))    
      CALL DACOP(ISCRDA( 16+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+DKI(IX,2)*DPDA*C1E3/(ONE+DPDA)*TILTC(I)                 *FOX
*FOX  -C1E3*DKI(IX,2)/(ONE+DPDA)*(ONE-TILTC(I)) ;                       *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  3+IDAA))          
      RSCRRI(  4+IDAA) = ONE         - RSCRRI(  1+IDAA)                 
      CALL DACOP(Y          ((2)),ISCRDA(  5+IDAA))                     
      RSCRRI(  6+IDAA) = DKI        (IX         ,(2))                   
      RSCRRI(  7+IDAA) = TILTC      (I          )                       
      RSCRRI(  8+IDAA) = DKI        (IX         ,(2))                   
      CALL DACMU(DPDA       ,ONE*RSCRRI(  6+IDAA),ISCRDA(  9+IDAA))     
      CALL DACMU(ISCRDA(  9+IDAA),ONE*C1E3       ,ISCRDA( 10+IDAA))     
      CALL DADIV(ISCRDA( 10+IDAA),ISCRDA(  2+IDAA),ISCRDA( 11+IDAA))    
      CALL DACMU(ISCRDA( 11+IDAA),ONE*RSCRRI(  7+IDAA),ISCRDA( 12+IDAA))
     *                                                                  
      RSCRRI( 13+IDAA) = C1E3        * RSCRRI(  8+IDAA)                 
      CALL DADIC(ISCRDA(  3+IDAA),ONE*RSCRRI( 13+IDAA),ISCRDA( 14+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA( 14+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA( 15+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA( 12+IDAA),ISCRDA( 16+IDAA))    
      CALL DASUB(ISCRDA( 16+IDAA),ISCRDA( 15+IDAA),ISCRDA( 17+IDAA))    
      CALL DACOP(ISCRDA( 17+IDAA),Y          ((2)))                     
            endif
            if(idp.eq.1.and.iabs(ition).eq.1) then
*FOX  SIGMDA=SIGMDA-RV*DKI(IX,2)*ZL ;                                   *FOX
      RSCRRI(  1+IDAA) = DKI        (IX         ,(2))                   
      CALL DACMU(RV         ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DAMUL(ISCRDA(  2+IDAA),ZL         ,ISCRDA(  3+IDAA))         
      CALL DASUB(SIGMDA     ,ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))         
      CALL DACOP(ISCRDA(  4+IDAA),SIGMDA     )                          
            endif
          endif
        if(abs(r0).le.pieni.or.nmz.eq.0) goto 440
        if(nmz.ge.2) then
*FOX  YV1J=BBI(I,1)+BBI(I,2)*XL+AAI(I,2)*ZL ;                           *FOX
      RSCRRI(  1+IDAA) = BBI        (I          ,(1))                   
      RSCRRI(  2+IDAA) = BBI        (I          ,(2))                   
      RSCRRI(  3+IDAA) = AAI        (I          ,(2))                   
      CALL DACMU(XL         ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DACMU(ZL         ,ONE*RSCRRI(  3+IDAA),ISCRDA(  5+IDAA))     
      CALL DACAD(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  6+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DACOP(ISCRDA(  7+IDAA),YV1J       )                          
*FOX  YV2J=AAI(I,1)-BBI(I,2)*ZL+AAI(I,2)*XL ;                           *FOX
      RSCRRI(  1+IDAA) = AAI        (I          ,(1))                   
      RSCRRI(  2+IDAA) = BBI        (I          ,(2))                   
      RSCRRI(  3+IDAA) = AAI        (I          ,(2))                   
      CALL DACMU(ZL         ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DACMU(XL         ,ONE*RSCRRI(  3+IDAA),ISCRDA(  5+IDAA))     
      CALL DASUC(ISCRDA(  4+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  6+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))    
      CALL DACOP(ISCRDA(  7+IDAA),YV2J       )                          
*FOX  CRKVE=XL ;                                                        *FOX
      CALL DACOP(XL         ,CRKVE      )                               
*FOX  CIKVE=ZL ;                                                        *FOX
      CALL DACOP(ZL         ,CIKVE      )                               
          do 200 k=3,nmz
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  YV1J=YV1J+BBI(I,K)*CRKVE+AAI(I,K)*CIKVE ;                         *FOX
      RSCRRI(  1+IDAA) = BBI        (I          ,K          )           
      RSCRRI(  2+IDAA) = AAI        (I          ,K          )           
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(YV1J       ,ISCRDA(  3+IDAA),ISCRDA(  5+IDAA))         
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  4+IDAA),ISCRDA(  6+IDAA))    
      CALL DACOP(ISCRDA(  6+IDAA),YV1J       )                          
*FOX  YV2J=YV2J-BBI(I,K)*CIKVE+AAI(I,K)*CRKVE ;                         *FOX
      RSCRRI(  1+IDAA) = BBI        (I          ,K          )           
      RSCRRI(  2+IDAA) = AAI        (I          ,K          )           
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DASUB(YV2J       ,ISCRDA(  3+IDAA),ISCRDA(  5+IDAA))         
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  4+IDAA),ISCRDA(  6+IDAA))    
      CALL DACOP(ISCRDA(  6+IDAA),YV2J       )                          
  200     continue
*FOX  Y(1)=Y(1)+(TILTC(I)*YV1J-TILTS(I)*YV2J)/(ONE+DPDA) ;              *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(YV1J       ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(YV2J       ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DASUB(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  6+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  7+IDAA))                     
      CALL DADIV(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+(TILTC(I)*YV2J+TILTS(I)*YV1J)/(ONE+DPDA) ;              *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(YV2J       ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(YV1J       ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  6+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  7+IDAA))                     
      CALL DADIV(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),Y          ((2)))                     
        else
*FOX  Y(1)=Y(1)+(TILTC(I)*BBI(I,1)-TILTS(I)*AAI(I,1))/(ONE+DPDA) ;      *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = BBI        (I          ,(1))                   
      RSCRRI(  3+IDAA) = TILTS      (I          )                       
      RSCRRI(  4+IDAA) = AAI        (I          ,(1))                   
      RSCRRI(  5+IDAA) = RSCRRI(  1+IDAA) * RSCRRI(  2+IDAA)            
      RSCRRI(  6+IDAA) = RSCRRI(  3+IDAA) * RSCRRI(  4+IDAA)            
      RSCRRI(  7+IDAA) = RSCRRI(  5+IDAA) - RSCRRI(  6+IDAA)            
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  8+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  9+IDAA))                     
      CALL DADIC(ISCRDA(  8+IDAA),ONE*RSCRRI(  7+IDAA),ISCRDA( 10+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(ISCRDA( 11+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+(TILTC(I)*AAI(I,1)+TILTS(I)*BBI(I,1))/(ONE+DPDA) ;      *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = AAI        (I          ,(1))                   
      RSCRRI(  3+IDAA) = TILTS      (I          )                       
      RSCRRI(  4+IDAA) = BBI        (I          ,(1))                   
      RSCRRI(  5+IDAA) = RSCRRI(  1+IDAA) * RSCRRI(  2+IDAA)            
      RSCRRI(  6+IDAA) = RSCRRI(  3+IDAA) * RSCRRI(  4+IDAA)            
      RSCRRI(  7+IDAA) = RSCRRI(  5+IDAA) + RSCRRI(  6+IDAA)            
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  8+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  9+IDAA))                     
      CALL DADIC(ISCRDA(  8+IDAA),ONE*RSCRRI(  7+IDAA),ISCRDA( 10+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(ISCRDA( 11+IDAA),Y          ((2)))                     
        endif
        goto 440
!--SKEW ELEMENTS
  320   kzz=-kzz
        goto(330,340,350,360,370,380,390,400,410,420),kzz
        goto 440
!---VERTICAL DIPOLE
  330   continue
*FOX  EKK=EKK*C1E3 ;                                                    *FOX
      CALL DACMU(EKK        ,ONE*C1E3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  Y(1)=Y(1)-EKK*TILTS(I) ;                                          *FOX
      CALL DACOP(Y          ((1)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(EKK        ,ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))     
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*TILTC(I) ;                                          *FOX
      CALL DACOP(Y          ((2)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = TILTC      (I          )                       
      CALL DACMU(EKK        ,ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))     
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),Y          ((2)))                     
        goto 440
!---SKEW QUADRUPOLE
  340   continue
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DASUB(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((2)))                     
        goto 440
!---SKEW SEXTUPOLE
  350   continue
*FOX  EKK=EKK*C1M3 ;                                                    *FOX
      CALL DACMU(EKK        ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DASUB(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((2)))                     
        goto 440
!---SKEW OCTUPOLE
  360   continue
*FOX  EKK=EKK*C1M6 ;                                                    *FOX
      CALL DACMU(EKK        ,ONE*C1M6       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DASUB(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((2)))                     
        goto 440
!---SKEW DECAPOLE
  370   continue
*FOX  EKK=EKK*C1M9 ;                                                    *FOX
      CALL DACMU(EKK        ,ONE*C1M9       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DASUB(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((2)))                     
        goto 440
!---SKEW DODECAPOL
  380   continue
*FOX  EKK=EKK*C1M12 ;                                                   *FOX
      CALL DACMU(EKK        ,ONE*C1M12      ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DASUB(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((2)))                     
        goto 440
!---SKEW 14-POL
  390   continue
*FOX  EKK=EKK*C1M15 ;                                                   *FOX
      CALL DACMU(EKK        ,ONE*C1M15      ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DASUB(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((2)))                     
        goto 440
!---SKEW 16-POL
  400   continue
*FOX  EKK=EKK*C1M18 ;                                                   *FOX
      CALL DACMU(EKK        ,ONE*C1M18      ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DASUB(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((2)))                     
        goto 440
!---SKEW 18-POL
  410   continue
*FOX  EKK=EKK*C1M21 ;                                                   *FOX
      CALL DACMU(EKK        ,ONE*C1M21      ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DASUB(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((2)))                     
        goto 440
!---SKEW 20-POL
  420   continue
*FOX  EKK=EKK*C1M24 ;                                                   *FOX
      CALL DACMU(EKK        ,ONE*C1M24      ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),EKK        )                          
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  CRKVEUK=CRKVE*XL-CIKVE*ZL ;                                       *FOX
      CALL DAMUL(CRKVE      ,XL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,ZL         ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CRKVEUK    )                          
*FOX  CIKVE=CRKVE*ZL+CIKVE*XL ;                                         *FOX
      CALL DAMUL(CRKVE      ,ZL         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(CIKVE      ,XL         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),CIKVE      )                          
*FOX  CRKVE=CRKVEUK ;                                                   *FOX
      CALL DACOP(CRKVEUK    ,CRKVE      )                               
*FOX  Y(1)=Y(1)+EKK*(TILTC(I)*CIKVE-TILTS(I)*CRKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DASUB(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((1)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((1)))                     
*FOX  Y(2)=Y(2)+EKK*(TILTC(I)*CRKVE+TILTS(I)*CIKVE) ;                   *FOX
      RSCRRI(  1+IDAA) = TILTC      (I          )                       
      RSCRRI(  2+IDAA) = TILTS      (I          )                       
      CALL DACMU(CRKVE      ,ONE*RSCRRI(  1+IDAA),ISCRDA(  3+IDAA))     
      CALL DACMU(CIKVE      ,ONE*RSCRRI(  2+IDAA),ISCRDA(  4+IDAA))     
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(Y          ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(EKK        ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Y          ((2)))                     
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
*FOX  YP(1)=Y(1)*(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((1)))                     
*FOX  YP(2)=Y(2)*(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((2)))                     
*FOX  DPDA1=DPDA*C1E3 ;                                                 *FOX
      CALL DACMU(DPDA       ,ONE*C1E3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA1      )                          
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
*FOX  YP(1)=Y(1)*(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((1)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((1)))                     
*FOX  YP(2)=Y(2)*(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(Y          ((2)),ISCRDA(  2+IDAA))                     
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),YP         ((2)))                     
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
*FOX  CORRAU1(1)=X(1) ;                                                 *FOX
      CALL DACOP(X          ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRAU1    ((1)))                     
*FOX  CORRAU1(2)=YP(1) ;                                                *FOX
      CALL DACOP(YP         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRAU1    ((2)))                     
*FOX  CORRAU1(3)=X(2) ;                                                 *FOX
      CALL DACOP(X          ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRAU1    ((3)))                     
*FOX  CORRAU1(4)=YP(2) ;                                                *FOX
      CALL DACOP(YP         ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRAU1    ((4)))                     
*FOX  CORRAU1(5)=SIGMDA ;                                               *FOX
      CALL DACOP(SIGMDA     ,CORRAU1    ((5)))                          
*FOX  CORRAU1(6)=DPDA1 ;                                                *FOX
      CALL DACOP(DPDA1      ,CORRAU1    ((6)))                          
        do 435 kkk=1,6
          dpdav2(kkk)=dare(corrau1(kkk))
*FOX  CORRAU1(KKK)=CORRAU1(KKK)-DPDAV2(KKK) ;                           *FOX
      CALL DACOP(CORRAU1    (KKK        ),ISCRDA(  1+IDAA))             
      RSCRRI(  2+IDAA) = DPDAV2     (KKK        )                       
      CALL DACSU(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  3+IDAA),CORRAU1    (KKK        ))             
  435   continue
        if(ivar.gt.ivar1) then
*FOX  CORRAU1(7)=SMIDA(1) ;                                             *FOX
      CALL DACOP(SMIDA      ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRAU1    ((7)))                     
*FOX  CORRAU1(8)=SMIDA(2) ;                                             *FOX
      CALL DACOP(SMIDA      ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),CORRAU1    ((8)))                     
          dpdav=dare(smida(1))
*FOX  CORRNEW(7)=SMIDA(1)-DPDAV ;                                       *FOX
      CALL DACOP(SMIDA      ((1)),ISCRDA(  1+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*DPDAV      ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),CORRNEW    ((7)))                     
          dpdav=dare(smida(2))
*FOX  CORRNEW(8)=SMIDA(2)-DPDAV ;                                       *FOX
      CALL DACOP(SMIDA      ((2)),ISCRDA(  1+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*DPDAV      ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),CORRNEW    ((8)))                     
        endif
        call dacct(corrau1,nvar,corrnew,nvar,corrau2,nvar)
        do 436 kkk=1,6
*FOX  CORRAU2(KKK)=CORRAU2(KKK)+DPDAV2(KKK) ;                           *FOX
      CALL DACOP(CORRAU2    (KKK        ),ISCRDA(  1+IDAA))             
      RSCRRI(  2+IDAA) = DPDAV2     (KKK        )                       
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  3+IDAA),CORRAU2    (KKK        ))             
  436   continue
*FOX  CORRAU1(2)=CORRAU2(2)/(ONE+CORRAU2(6)) ;                          *FOX
      CALL DACOP(CORRAU2    ((6)),ISCRDA(  1+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  2+IDAA))     
      CALL DACOP(CORRAU2    ((2)),ISCRDA(  3+IDAA))                     
      CALL DADIV(ISCRDA(  3+IDAA),ISCRDA(  2+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),CORRAU1    ((2)))                     
*FOX  CORRAU1(4)=CORRAU2(4)/(ONE+CORRAU2(6)) ;                          *FOX
      CALL DACOP(CORRAU2    ((6)),ISCRDA(  1+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  2+IDAA))     
      CALL DACOP(CORRAU2    ((4)),ISCRDA(  3+IDAA))                     
      CALL DADIV(ISCRDA(  3+IDAA),ISCRDA(  2+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),CORRAU1    ((4)))                     
*FOX  X(1)=CORRAU2(1) ;                                                 *FOX
      CALL DACOP(CORRAU2    ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),X          ((1)))                     
*FOX  YP(1)=CORRAU2(2) ;                                                *FOX
      CALL DACOP(CORRAU2    ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),YP         ((1)))                     
*FOX  X(2)=CORRAU2(3) ;                                                 *FOX
      CALL DACOP(CORRAU2    ((3)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),X          ((2)))                     
*FOX  YP(2)=CORRAU2(4) ;                                                *FOX
      CALL DACOP(CORRAU2    ((4)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),YP         ((2)))                     
*FOX  SIGMDA=CORRAU2(5) ;                                               *FOX
      CALL DACOP(CORRAU2    ((5)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),SIGMDA     )                          
*FOX  DPDA1=CORRAU2(6) ;                                                *FOX
      CALL DACOP(CORRAU2    ((6)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),DPDA1      )                          
*FOX  DPDA=DPDA1*C1M3 ;                                                 *FOX
      CALL DACMU(DPDA1      ,ONE*C1M3       ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),DPDA       )                          
*FOX  Y(1)=YP(1)/(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(YP         ((1)),ISCRDA(  2+IDAA))                     
      CALL DADIV(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((1)))                     
*FOX  Y(2)=YP(2)/(ONE+DPDA) ;                                           *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACOP(YP         ((2)),ISCRDA(  2+IDAA))                     
      CALL DADIV(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),Y          ((2)))                     
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
        CALL DADAL(PZ      ,1)                                                  
        CALL DADAL(CRABAMP4,1)                                                  
        CALL DADAL(CRABAMP3,1)                                                  
        CALL DADAL(CRABAMP2,1)                                                  
        CALL DADAL(CRABAMP ,1)                                                  
        CALL DADAL(WY      ,1)                                                  
        CALL DADAL(WX      ,1)                                                  
        CALL DADAL(CRZBF   ,1)                                                  
        CALL DADAL(CBXBF   ,1)                                                  
        CALL DADAL(CRXBF   ,1)                                                  
        CALL DADAL(ZBBF    ,1)                                                  
        CALL DADAL(XBBF    ,1)                                                  
        CALL DADAL(ZRBF    ,1)                                                  
        CALL DADAL(CCCC    ,1)                                                  
        CALL DADAL(XRBF    ,1)                                                  
        CALL DADAL(TKBF    ,1)                                                  
        CALL DADAL(RHO2BF  ,1)                                                  
        CALL DADAL(CIKVEBF ,1)                                                  
        CALL DADAL(CRKVEBF ,1)                                                  
        CALL DADAL(YV2J    ,1)                                                  
        CALL DADAL(YV1J    ,1)                                                  
        CALL DADAL(CBZBF   ,1)                                                  
        CALL DADAL(CRKVEUK ,1)                                                  
        CALL DADAL(CIKVE   ,1)                                                  
        CALL DADAL(CRKVE   ,1)                                                  
        CALL DADAL(ZL      ,1)                                                  
        CALL DADAL(XL      ,1)                                                  
        CALL DADAL(EKK     ,1)                                                  
        CALL DADAL(EJF0    ,1)                                                  
        CALL DADAL(PUZ     ,1)                                                  
        CALL DADAL(PUX     ,1)                                                  
        CALL DADAL(TRACKI  ,1*(6))                                              
        CALL DADAL(BB      ,1*(11))                                             
        CALL DADAL(AA      ,1*(11))                                             
        CALL DADAL(CORRAU2 ,1*(MCOP))                                           
        CALL DADAL(CORRAU1 ,1*(MCOP))                                           
        CALL DADAL(CORRNEW ,1*(MCOP))                                           
        CALL DADAL(CORROLD ,1*(MCOP))                                           
        CALL DADAL(DKIP    ,1)                                                  
        CALL DADAL(YP      ,1*(2))                                              
        CALL DADAL(Y       ,1*(2))                                              
        CALL DADAL(X       ,1*(2))                                              
        CALL DADAL(SMIDA   ,1*(MCOR))                                           
        CALL DADAL(ASDAQ   ,1*(2)*(6))                                          
        CALL DADAL(ALDAQ   ,1*(2)*(6))                                          
        CALL DADAL(ASDA    ,1*(2)*(6))                                          
        CALL DADAL(ALDA    ,1*(2)*(6))                                          
        CALL DADAL(EJF1    ,1)                                                  
        CALL DADAL(EJ1     ,1)                                                  
        CALL DADAL(YY      ,1*(2))                                              
        CALL DADAL(XX      ,1*(2))                                              
        CALL DADAL(RV      ,1)                                                  
        CALL DADAL(DPDA1   ,1)                                                  
        CALL DADAL(DPDA    ,1)                                                  
        CALL DADAL(SIGMDA  ,1)                                                  
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
10030 format('|',i6,'|',a8,'|',f12.5,'|','X','|',f12.7,'|',f12.6,'|',   *FOX
     &f13.7,'|',f11.6,'|',f11.7,'|',f11.7,'|',f11.7,'|',f11.7,'|')
10040 format('|',6x,'|',8x,'|',12x,'|','Y','|',12x,'|',f12.6,'|', f13.7,*FOX
     &'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10050 format('|',6x,'|',a8,'|',12x,'|','S','|',12x,'|',f12.6,'|', f13.  &
     &7,'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10055 format('|',6x,'|',a8,'|',12x,'|','Y','|',12x,'|',f12.6,'|', f13.  *FOX
     &7,'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10060 format('|',6x,'|',8x,'|',12x,'|',102('-'))
10070 format('|',6x,'|',8x,'|',12x,'|','Y','|',f12.7,'|',f12.6,'|', f13.*FOX
     &7,'|',f11.6,'|',f11.7,'|',f11.7,'|',f11.7,'|',f11.7,'|')
10080 format('|',6x,'|',8x,'|',12x,'|','X','|',12x,'|',f12.6,'|', f13.7,*FOX
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
      double precision al,as,chi0,chid,dp1,dps,exz,sigm
      common/syos/as(6,2,npart,nele),al(6,2,npart,nele),sigm(mpa),      &
     &dps(mpa),idz(2)
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
*FOX{
      INTEGER SIGMDA  
      INTEGER DPDA    
      INTEGER DPDA1   
      INTEGER RV      
      INTEGER XX      (2)
      INTEGER YY      (2)
      INTEGER EJ1     
      INTEGER EJF1    
      INTEGER ALDA    (2,6)
      INTEGER ASDA    (2,6)
      INTEGER ALDAQ   (2,6)
      INTEGER ASDAQ   (2,6)
      INTEGER SMIDA   (MCOR)
      INTEGER XI      
      INTEGER YI      
      INTEGER ISCRDA, ISCRRI,IDAO
      DOUBLE PRECISION RSCRRI
      COMMON/DASCR/ISCRDA(100),RSCRRI(100),ISCRRI(100),IDAO
      if(1.eq.1) then                                                 
         CALL DAKEY('FOX V2.1')
         CALL DAALL(XI      ,1,'XI        ',NORD,NVAR)
         CALL DAALL(YI      ,1,'YI        ',NORD,NVAR)
      ENDIF
      IDAA = IDAO
*FOX}
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
       
       
*FOX  XX(1)=XX(1)*C1M3;                                                 *FOX
      CALL DACOP(XX         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1M3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),XX         ((1)))                     
*FOX  XX(2)=XX(2)*C1M3;                                                 *FOX
      CALL DACOP(XX         ((2)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1M3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),XX         ((2)))                     
*FOX  YY(1)=YY(1)*C1M3;                                                 *FOX
      CALL DACOP(YY         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1M3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),YY         ((1)))                     
*FOX  YY(2)=YY(2)*C1M3;                                                 *FOX
      CALL DACOP(YY         ((2)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1M3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),YY         ((2)))                     
       
!      CALL DRIFT(-EMBL/2)
*FOX  XX(1)=XX(1)-EMBL/TWO*YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-            *FOX
*FOX  YY(1)*YY(1)-YY(2)*YY(2)) ;                                        *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(XX         ((1)),ISCRDA( 12+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA( 13+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 14+IDAA))              
      RSCRRI( 15+IDAA) = EMBL        / TWO                              
      CALL DACMU(ISCRDA( 13+IDAA),ONE*RSCRRI( 15+IDAA),ISCRDA( 16+IDAA))
     *                                                                  
      CALL DADIV(ISCRDA( 16+IDAA),ISCRDA( 14+IDAA),ISCRDA( 17+IDAA))    
      CALL DASUB(ISCRDA( 12+IDAA),ISCRDA( 17+IDAA),ISCRDA( 18+IDAA))    
      CALL DACOP(ISCRDA( 18+IDAA),XX         ((1)))                     
*FOX  XX(2)=XX(2)-EMBL/TWO*YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-            *FOX
*FOX  YY(1)*YY(1)-YY(2)*YY(2)) ;                                        *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(XX         ((2)),ISCRDA( 12+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA( 13+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 14+IDAA))              
      RSCRRI( 15+IDAA) = EMBL        / TWO                              
      CALL DACMU(ISCRDA( 13+IDAA),ONE*RSCRRI( 15+IDAA),ISCRDA( 16+IDAA))
     *                                                                  
      CALL DADIV(ISCRDA( 16+IDAA),ISCRDA( 14+IDAA),ISCRDA( 17+IDAA))    
      CALL DASUB(ISCRDA( 12+IDAA),ISCRDA( 17+IDAA),ISCRDA( 18+IDAA))    
      CALL DACOP(ISCRDA( 18+IDAA),XX         ((2)))                     
!      CALL TILT(TX,TY)
*FOX  XX(2)=XX(2)-XX(1)*SIN(TX)*YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-       *FOX
*FOX  YY(2)*YY(2))/COS(ATAN(YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-           *FOX
*FOX  YY(1)*YY(1)-YY(2)*YY(2)))-TX) ;                                   *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(YY         ((1)),ISCRDA( 12+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))              
      CALL DADIV(ISCRDA( 12+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 15+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 16+IDAA))          
      CALL DACOP(YY         ((2)),ISCRDA( 17+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA( 18+IDAA))                     
      CALL DAFUN('ATAN',ISCRDA( 14+IDAA),ISCRDA( 19+IDAA))              
      CALL DAMUL(ISCRDA( 15+IDAA),ISCRDA( 16+IDAA),ISCRDA( 20+IDAA))    
      CALL DAMUL(ISCRDA( 17+IDAA),ISCRDA( 18+IDAA),ISCRDA( 21+IDAA))    
      CALL DASUB(ISCRDA( 20+IDAA),ISCRDA( 21+IDAA),ISCRDA( 22+IDAA))    
      CALL DACSU(ISCRDA( 19+IDAA),ONE*TX         ,ISCRDA( 23+IDAA))     
      CALL DACOP(XX         ((2)),ISCRDA( 24+IDAA))                     
      CALL DACOP(XX         ((1)),ISCRDA( 25+IDAA))                     
      RSCRRI( 26+IDAA) = SIN (TX         )                              
      CALL DACOP(YY         ((2)),ISCRDA( 27+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 22+IDAA),ISCRDA( 28+IDAA))              
      CALL DAFUN('COS ',ISCRDA( 23+IDAA),ISCRDA( 29+IDAA))              
      CALL DACMU(ISCRDA( 25+IDAA),ONE*RSCRRI( 26+IDAA),ISCRDA( 30+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA( 30+IDAA),ISCRDA( 27+IDAA),ISCRDA( 31+IDAA))    
      CALL DADIV(ISCRDA( 31+IDAA),ISCRDA( 28+IDAA),ISCRDA( 32+IDAA))    
      CALL DADIV(ISCRDA( 32+IDAA),ISCRDA( 29+IDAA),ISCRDA( 33+IDAA))    
      CALL DASUB(ISCRDA( 24+IDAA),ISCRDA( 33+IDAA),ISCRDA( 34+IDAA))    
      CALL DACOP(ISCRDA( 34+IDAA),XX         ((2)))                     
*FOX  XX(1)=XX(1)*(COS(TX)-SIN(TX)*TAN(ATAN(YY(1)/SQRT((ONE+DPDA)*      *FOX
*FOX  (ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TX)) ;                       *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(YY         ((1)),ISCRDA( 12+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))              
      CALL DADIV(ISCRDA( 12+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DAFUN('ATAN',ISCRDA( 14+IDAA),ISCRDA( 15+IDAA))              
      CALL DACSU(ISCRDA( 15+IDAA),ONE*TX         ,ISCRDA( 16+IDAA))     
      RSCRRI( 17+IDAA) = COS (TX         )                              
      RSCRRI( 18+IDAA) = SIN (TX         )                              
      CALL DAFUN('TAN ',ISCRDA( 16+IDAA),ISCRDA( 19+IDAA))              
      CALL DACMU(ISCRDA( 19+IDAA),ONE*RSCRRI( 18+IDAA),ISCRDA( 20+IDAA))
     *                                                                  
      CALL DASUC(ISCRDA( 20+IDAA),ONE*RSCRRI( 17+IDAA),ISCRDA( 21+IDAA))
     *                                                                  
      CALL DACOP(XX         ((1)),ISCRDA( 22+IDAA))                     
      CALL DAMUL(ISCRDA( 22+IDAA),ISCRDA( 21+IDAA),ISCRDA( 23+IDAA))    
      CALL DACOP(ISCRDA( 23+IDAA),XX         ((1)))                     
*FOX  YY(1)=SQRT((ONE+DPDA)*(ONE+DPDA)-YY(2)*YY(2))*SIN(ATAN(YY(1)/     *FOX
*FOX  SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TX) ;        *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(YY         ((1)),ISCRDA( 12+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))              
      CALL DADIV(ISCRDA( 12+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 15+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 16+IDAA))          
      CALL DACOP(YY         ((2)),ISCRDA( 17+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA( 18+IDAA))                     
      CALL DAFUN('ATAN',ISCRDA( 14+IDAA),ISCRDA( 19+IDAA))              
      CALL DAMUL(ISCRDA( 15+IDAA),ISCRDA( 16+IDAA),ISCRDA( 20+IDAA))    
      CALL DAMUL(ISCRDA( 17+IDAA),ISCRDA( 18+IDAA),ISCRDA( 21+IDAA))    
      CALL DASUB(ISCRDA( 20+IDAA),ISCRDA( 21+IDAA),ISCRDA( 22+IDAA))    
      CALL DACSU(ISCRDA( 19+IDAA),ONE*TX         ,ISCRDA( 23+IDAA))     
      CALL DAFUN('SQRT',ISCRDA( 22+IDAA),ISCRDA( 24+IDAA))              
      CALL DAFUN('SIN ',ISCRDA( 23+IDAA),ISCRDA( 25+IDAA))              
      CALL DAMUL(ISCRDA( 24+IDAA),ISCRDA( 25+IDAA),ISCRDA( 26+IDAA))    
      CALL DACOP(ISCRDA( 26+IDAA),YY         ((1)))                     
*FOX  XX(1)=XX(1)-XX(2)*SIN(TY)*YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-       *FOX
*FOX  YY(1)*YY(1))/COS(ATAN(YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-           *FOX
*FOX  YY(1)*YY(1)-YY(2)*YY(2)))-TY) ;                                   *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(YY         ((2)),ISCRDA( 12+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))              
      CALL DADIV(ISCRDA( 12+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 15+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 16+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA( 17+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA( 18+IDAA))                     
      CALL DAFUN('ATAN',ISCRDA( 14+IDAA),ISCRDA( 19+IDAA))              
      CALL DAMUL(ISCRDA( 15+IDAA),ISCRDA( 16+IDAA),ISCRDA( 20+IDAA))    
      CALL DAMUL(ISCRDA( 17+IDAA),ISCRDA( 18+IDAA),ISCRDA( 21+IDAA))    
      CALL DASUB(ISCRDA( 20+IDAA),ISCRDA( 21+IDAA),ISCRDA( 22+IDAA))    
      CALL DACSU(ISCRDA( 19+IDAA),ONE*TY         ,ISCRDA( 23+IDAA))     
      CALL DACOP(XX         ((1)),ISCRDA( 24+IDAA))                     
      CALL DACOP(XX         ((2)),ISCRDA( 25+IDAA))                     
      RSCRRI( 26+IDAA) = SIN (TY         )                              
      CALL DACOP(YY         ((1)),ISCRDA( 27+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 22+IDAA),ISCRDA( 28+IDAA))              
      CALL DAFUN('COS ',ISCRDA( 23+IDAA),ISCRDA( 29+IDAA))              
      CALL DACMU(ISCRDA( 25+IDAA),ONE*RSCRRI( 26+IDAA),ISCRDA( 30+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA( 30+IDAA),ISCRDA( 27+IDAA),ISCRDA( 31+IDAA))    
      CALL DADIV(ISCRDA( 31+IDAA),ISCRDA( 28+IDAA),ISCRDA( 32+IDAA))    
      CALL DADIV(ISCRDA( 32+IDAA),ISCRDA( 29+IDAA),ISCRDA( 33+IDAA))    
      CALL DASUB(ISCRDA( 24+IDAA),ISCRDA( 33+IDAA),ISCRDA( 34+IDAA))    
      CALL DACOP(ISCRDA( 34+IDAA),XX         ((1)))                     
*FOX  XX(2)=XX(2)*(COS(TY)-SIN(TY)*TAN(ATAN(YY(2)/SQRT((ONE+DPDA)*      *FOX
*FOX  (ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TY)) ;                       *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(YY         ((2)),ISCRDA( 12+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))              
      CALL DADIV(ISCRDA( 12+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DAFUN('ATAN',ISCRDA( 14+IDAA),ISCRDA( 15+IDAA))              
      CALL DACSU(ISCRDA( 15+IDAA),ONE*TY         ,ISCRDA( 16+IDAA))     
      RSCRRI( 17+IDAA) = COS (TY         )                              
      RSCRRI( 18+IDAA) = SIN (TY         )                              
      CALL DAFUN('TAN ',ISCRDA( 16+IDAA),ISCRDA( 19+IDAA))              
      CALL DACMU(ISCRDA( 19+IDAA),ONE*RSCRRI( 18+IDAA),ISCRDA( 20+IDAA))
     *                                                                  
      CALL DASUC(ISCRDA( 20+IDAA),ONE*RSCRRI( 17+IDAA),ISCRDA( 21+IDAA))
     *                                                                  
      CALL DACOP(XX         ((2)),ISCRDA( 22+IDAA))                     
      CALL DAMUL(ISCRDA( 22+IDAA),ISCRDA( 21+IDAA),ISCRDA( 23+IDAA))    
      CALL DACOP(ISCRDA( 23+IDAA),XX         ((2)))                     
*FOX  YY(2)=SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1))*SIN(ATAN(YY(2)/     *FOX
*FOX  SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TY) ;        *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(YY         ((2)),ISCRDA( 12+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))              
      CALL DADIV(ISCRDA( 12+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 15+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 16+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA( 17+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA( 18+IDAA))                     
      CALL DAFUN('ATAN',ISCRDA( 14+IDAA),ISCRDA( 19+IDAA))              
      CALL DAMUL(ISCRDA( 15+IDAA),ISCRDA( 16+IDAA),ISCRDA( 20+IDAA))    
      CALL DAMUL(ISCRDA( 17+IDAA),ISCRDA( 18+IDAA),ISCRDA( 21+IDAA))    
      CALL DASUB(ISCRDA( 20+IDAA),ISCRDA( 21+IDAA),ISCRDA( 22+IDAA))    
      CALL DACSU(ISCRDA( 19+IDAA),ONE*TY         ,ISCRDA( 23+IDAA))     
      CALL DAFUN('SQRT',ISCRDA( 22+IDAA),ISCRDA( 24+IDAA))              
      CALL DAFUN('SIN ',ISCRDA( 23+IDAA),ISCRDA( 25+IDAA))              
      CALL DAMUL(ISCRDA( 24+IDAA),ISCRDA( 25+IDAA),ISCRDA( 26+IDAA))    
      CALL DACOP(ISCRDA( 26+IDAA),YY         ((2)))                     
!      CALL DRIFT(LIN)
*FOX  XX(1)=XX(1)+LIN*YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-     *FOX
*FOX  YY(2)*YY(2)) ;                                                    *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(XX         ((1)),ISCRDA( 12+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA( 13+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 14+IDAA))              
      CALL DACMU(ISCRDA( 13+IDAA),ONE*LIN        ,ISCRDA( 15+IDAA))     
      CALL DADIV(ISCRDA( 15+IDAA),ISCRDA( 14+IDAA),ISCRDA( 16+IDAA))    
      CALL DAADD(ISCRDA( 12+IDAA),ISCRDA( 16+IDAA),ISCRDA( 17+IDAA))    
      CALL DACOP(ISCRDA( 17+IDAA),XX         ((1)))                     
*FOX  XX(2)=XX(2)+LIN*YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-     *FOX
*FOX  YY(2)*YY(2)) ;                                                    *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(XX         ((2)),ISCRDA( 12+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA( 13+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 14+IDAA))              
      CALL DACMU(ISCRDA( 13+IDAA),ONE*LIN        ,ISCRDA( 15+IDAA))     
      CALL DADIV(ISCRDA( 15+IDAA),ISCRDA( 14+IDAA),ISCRDA( 16+IDAA))    
      CALL DAADD(ISCRDA( 12+IDAA),ISCRDA( 16+IDAA),ISCRDA( 17+IDAA))    
      CALL DACOP(ISCRDA( 17+IDAA),XX         ((2)))                     
!      CALL KICK(L,CUR,LIN,RX,RY)
*FOX  XI=XX(1)-RX ;                                                     *FOX
      CALL DACOP(XX         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*RX         ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),XI         )                          
*FOX  YI=XX(2)-RY ;                                                     *FOX
      CALL DACOP(XX         ((2)),ISCRDA(  1+IDAA))                     
      CALL DACSU(ISCRDA(  1+IDAA),ONE*RY         ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),YI         )                          
*FOX  YY(1)=YY(1)-C1M7*CUR/CHI*XI/(XI*XI+YI*YI)*(SQRT((LIN+L)*(LIN+L)+  *FOX
*FOX  XI*XI+YI*YI)-SQRT((LIN-L)*(LIN-L)+XI*XI+YI*YI)) ;                 *FOX
      RSCRRI(  1+IDAA) = LIN         + L                                
      RSCRRI(  2+IDAA) = LIN         + L                                
      RSCRRI(  3+IDAA) = LIN         - L                                
      RSCRRI(  4+IDAA) = LIN         - L                                
      RSCRRI(  5+IDAA) = RSCRRI(  1+IDAA) * RSCRRI(  2+IDAA)            
      CALL DAMUL(XI         ,XI         ,ISCRDA(  6+IDAA))              
      CALL DAMUL(YI         ,YI         ,ISCRDA(  7+IDAA))              
      RSCRRI(  8+IDAA) = RSCRRI(  3+IDAA) * RSCRRI(  4+IDAA)            
      CALL DAMUL(XI         ,XI         ,ISCRDA(  9+IDAA))              
      CALL DAMUL(YI         ,YI         ,ISCRDA( 10+IDAA))              
      CALL DACAD(ISCRDA(  6+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA( 11+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA( 11+IDAA),ISCRDA(  7+IDAA),ISCRDA( 12+IDAA))    
      CALL DACAD(ISCRDA(  9+IDAA),ONE*RSCRRI(  8+IDAA),ISCRDA( 13+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA( 13+IDAA),ISCRDA( 10+IDAA),ISCRDA( 14+IDAA))    
      CALL DAFUN('SQRT',ISCRDA( 12+IDAA),ISCRDA( 15+IDAA))              
      CALL DAFUN('SQRT',ISCRDA( 14+IDAA),ISCRDA( 16+IDAA))              
      CALL DAMUL(XI         ,XI         ,ISCRDA( 17+IDAA))              
      CALL DAMUL(YI         ,YI         ,ISCRDA( 18+IDAA))              
      CALL DAADD(ISCRDA( 17+IDAA),ISCRDA( 18+IDAA),ISCRDA( 19+IDAA))    
      CALL DASUB(ISCRDA( 15+IDAA),ISCRDA( 16+IDAA),ISCRDA( 20+IDAA))    
      CALL DACOP(YY         ((1)),ISCRDA( 21+IDAA))                     
      RSCRRI( 22+IDAA) = C1M7        * CUR                              
      RSCRRI( 23+IDAA) = RSCRRI( 22+IDAA) / CHI                         
      CALL DACMU(XI         ,ONE*RSCRRI( 23+IDAA),ISCRDA( 24+IDAA))     
      CALL DADIV(ISCRDA( 24+IDAA),ISCRDA( 19+IDAA),ISCRDA( 25+IDAA))    
      CALL DAMUL(ISCRDA( 25+IDAA),ISCRDA( 20+IDAA),ISCRDA( 26+IDAA))    
      CALL DASUB(ISCRDA( 21+IDAA),ISCRDA( 26+IDAA),ISCRDA( 27+IDAA))    
      CALL DACOP(ISCRDA( 27+IDAA),YY         ((1)))                     
*FOX  YY(2)=YY(2)-C1M7*CUR/CHI*YI/(XI*XI+YI*YI)*(SQRT((LIN+L)*(LIN+L)+  *FOX
*FOX  XI*XI+YI*YI)-SQRT((LIN-L)*(LIN-L)+XI*XI+YI*YI)) ;                 *FOX
      RSCRRI(  1+IDAA) = LIN         + L                                
      RSCRRI(  2+IDAA) = LIN         + L                                
      RSCRRI(  3+IDAA) = LIN         - L                                
      RSCRRI(  4+IDAA) = LIN         - L                                
      RSCRRI(  5+IDAA) = RSCRRI(  1+IDAA) * RSCRRI(  2+IDAA)            
      CALL DAMUL(XI         ,XI         ,ISCRDA(  6+IDAA))              
      CALL DAMUL(YI         ,YI         ,ISCRDA(  7+IDAA))              
      RSCRRI(  8+IDAA) = RSCRRI(  3+IDAA) * RSCRRI(  4+IDAA)            
      CALL DAMUL(XI         ,XI         ,ISCRDA(  9+IDAA))              
      CALL DAMUL(YI         ,YI         ,ISCRDA( 10+IDAA))              
      CALL DACAD(ISCRDA(  6+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA( 11+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA( 11+IDAA),ISCRDA(  7+IDAA),ISCRDA( 12+IDAA))    
      CALL DACAD(ISCRDA(  9+IDAA),ONE*RSCRRI(  8+IDAA),ISCRDA( 13+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA( 13+IDAA),ISCRDA( 10+IDAA),ISCRDA( 14+IDAA))    
      CALL DAFUN('SQRT',ISCRDA( 12+IDAA),ISCRDA( 15+IDAA))              
      CALL DAFUN('SQRT',ISCRDA( 14+IDAA),ISCRDA( 16+IDAA))              
      CALL DAMUL(XI         ,XI         ,ISCRDA( 17+IDAA))              
      CALL DAMUL(YI         ,YI         ,ISCRDA( 18+IDAA))              
      CALL DAADD(ISCRDA( 17+IDAA),ISCRDA( 18+IDAA),ISCRDA( 19+IDAA))    
      CALL DASUB(ISCRDA( 15+IDAA),ISCRDA( 16+IDAA),ISCRDA( 20+IDAA))    
      CALL DACOP(YY         ((2)),ISCRDA( 21+IDAA))                     
      RSCRRI( 22+IDAA) = C1M7        * CUR                              
      RSCRRI( 23+IDAA) = RSCRRI( 22+IDAA) / CHI                         
      CALL DACMU(YI         ,ONE*RSCRRI( 23+IDAA),ISCRDA( 24+IDAA))     
      CALL DADIV(ISCRDA( 24+IDAA),ISCRDA( 19+IDAA),ISCRDA( 25+IDAA))    
      CALL DAMUL(ISCRDA( 25+IDAA),ISCRDA( 20+IDAA),ISCRDA( 26+IDAA))    
      CALL DASUB(ISCRDA( 21+IDAA),ISCRDA( 26+IDAA),ISCRDA( 27+IDAA))    
      CALL DACOP(ISCRDA( 27+IDAA),YY         ((2)))                     
!     A = A-1E-7*2*L*CUR/CHI*XI/(XI*XI+YI*YI)                           *FOX
!     B = B-1E-7*2*L*CUR/CHI*YI/(XI*XI+YI*YI)                           *FOX
!      
!     THIS HAPPENS WHEN THE EMBEDDING DRIFT LENGTH IS PRACTICALLY INFINITE.
!     CALL DRIFT(LEFF-LIN)
*FOX  XX(1)=XX(1)+(LEFF-LIN)*YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-          *FOX
*FOX  YY(1)*YY(1)-YY(2)*YY(2)) ;                                        *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      RSCRRI( 10+IDAA) = LEFF        - LIN                              
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 11+IDAA))    
      CALL DASUB(ISCRDA( 11+IDAA),ISCRDA(  9+IDAA),ISCRDA( 12+IDAA))    
      CALL DACOP(XX         ((1)),ISCRDA( 13+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA( 14+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 12+IDAA),ISCRDA( 15+IDAA))              
      CALL DACMU(ISCRDA( 14+IDAA),ONE*RSCRRI( 10+IDAA),ISCRDA( 16+IDAA))
     *                                                                  
      CALL DADIV(ISCRDA( 16+IDAA),ISCRDA( 15+IDAA),ISCRDA( 17+IDAA))    
      CALL DAADD(ISCRDA( 13+IDAA),ISCRDA( 17+IDAA),ISCRDA( 18+IDAA))    
      CALL DACOP(ISCRDA( 18+IDAA),XX         ((1)))                     
*FOX  XX(2)=XX(2)+(LEFF-LIN)*YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-          *FOX
*FOX  YY(1)*YY(1)-YY(2)*YY(2)) ;                                        *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      RSCRRI( 10+IDAA) = LEFF        - LIN                              
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 11+IDAA))    
      CALL DASUB(ISCRDA( 11+IDAA),ISCRDA(  9+IDAA),ISCRDA( 12+IDAA))    
      CALL DACOP(XX         ((2)),ISCRDA( 13+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA( 14+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 12+IDAA),ISCRDA( 15+IDAA))              
      CALL DACMU(ISCRDA( 14+IDAA),ONE*RSCRRI( 10+IDAA),ISCRDA( 16+IDAA))
     *                                                                  
      CALL DADIV(ISCRDA( 16+IDAA),ISCRDA( 15+IDAA),ISCRDA( 17+IDAA))    
      CALL DAADD(ISCRDA( 13+IDAA),ISCRDA( 17+IDAA),ISCRDA( 18+IDAA))    
      CALL DACOP(ISCRDA( 18+IDAA),XX         ((2)))                     
!      CALL INVTILT(TX,TY)
*FOX  XX(1)=XX(1)+XX(2)*SIN(TY)*YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-       *FOX
*FOX  YY(1)*YY(1))/COS(ATAN(YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-           *FOX
*FOX  YY(1)*YY(1)-YY(2)*YY(2)))+TY) ;                                   *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(YY         ((2)),ISCRDA( 12+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))              
      CALL DADIV(ISCRDA( 12+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 15+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 16+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA( 17+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA( 18+IDAA))                     
      CALL DAFUN('ATAN',ISCRDA( 14+IDAA),ISCRDA( 19+IDAA))              
      CALL DAMUL(ISCRDA( 15+IDAA),ISCRDA( 16+IDAA),ISCRDA( 20+IDAA))    
      CALL DAMUL(ISCRDA( 17+IDAA),ISCRDA( 18+IDAA),ISCRDA( 21+IDAA))    
      CALL DASUB(ISCRDA( 20+IDAA),ISCRDA( 21+IDAA),ISCRDA( 22+IDAA))    
      CALL DACAD(ISCRDA( 19+IDAA),ONE*TY         ,ISCRDA( 23+IDAA))     
      CALL DACOP(XX         ((1)),ISCRDA( 24+IDAA))                     
      CALL DACOP(XX         ((2)),ISCRDA( 25+IDAA))                     
      RSCRRI( 26+IDAA) = SIN (TY         )                              
      CALL DACOP(YY         ((1)),ISCRDA( 27+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 22+IDAA),ISCRDA( 28+IDAA))              
      CALL DAFUN('COS ',ISCRDA( 23+IDAA),ISCRDA( 29+IDAA))              
      CALL DACMU(ISCRDA( 25+IDAA),ONE*RSCRRI( 26+IDAA),ISCRDA( 30+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA( 30+IDAA),ISCRDA( 27+IDAA),ISCRDA( 31+IDAA))    
      CALL DADIV(ISCRDA( 31+IDAA),ISCRDA( 28+IDAA),ISCRDA( 32+IDAA))    
      CALL DADIV(ISCRDA( 32+IDAA),ISCRDA( 29+IDAA),ISCRDA( 33+IDAA))    
      CALL DAADD(ISCRDA( 24+IDAA),ISCRDA( 33+IDAA),ISCRDA( 34+IDAA))    
      CALL DACOP(ISCRDA( 34+IDAA),XX         ((1)))                     
*FOX  XX(2)=XX(2)*(COS(TY)+SIN(TY)*TAN(ATAN(YY(2)/SQRT((ONE+DPDA)*      *FOX
*FOX  (ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))+TY)) ;                       *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(YY         ((2)),ISCRDA( 12+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))              
      CALL DADIV(ISCRDA( 12+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DAFUN('ATAN',ISCRDA( 14+IDAA),ISCRDA( 15+IDAA))              
      CALL DACAD(ISCRDA( 15+IDAA),ONE*TY         ,ISCRDA( 16+IDAA))     
      RSCRRI( 17+IDAA) = COS (TY         )                              
      RSCRRI( 18+IDAA) = SIN (TY         )                              
      CALL DAFUN('TAN ',ISCRDA( 16+IDAA),ISCRDA( 19+IDAA))              
      CALL DACMU(ISCRDA( 19+IDAA),ONE*RSCRRI( 18+IDAA),ISCRDA( 20+IDAA))
     *                                                                  
      CALL DACAD(ISCRDA( 20+IDAA),ONE*RSCRRI( 17+IDAA),ISCRDA( 21+IDAA))
     *                                                                  
      CALL DACOP(XX         ((2)),ISCRDA( 22+IDAA))                     
      CALL DAMUL(ISCRDA( 22+IDAA),ISCRDA( 21+IDAA),ISCRDA( 23+IDAA))    
      CALL DACOP(ISCRDA( 23+IDAA),XX         ((2)))                     
*FOX  YY(2)=SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1))*SIN(ATAN(YY(2)/     *FOX
*FOX  SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))+TY) ;        *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(YY         ((2)),ISCRDA( 12+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))              
      CALL DADIV(ISCRDA( 12+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 15+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 16+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA( 17+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA( 18+IDAA))                     
      CALL DAFUN('ATAN',ISCRDA( 14+IDAA),ISCRDA( 19+IDAA))              
      CALL DAMUL(ISCRDA( 15+IDAA),ISCRDA( 16+IDAA),ISCRDA( 20+IDAA))    
      CALL DAMUL(ISCRDA( 17+IDAA),ISCRDA( 18+IDAA),ISCRDA( 21+IDAA))    
      CALL DASUB(ISCRDA( 20+IDAA),ISCRDA( 21+IDAA),ISCRDA( 22+IDAA))    
      CALL DACAD(ISCRDA( 19+IDAA),ONE*TY         ,ISCRDA( 23+IDAA))     
      CALL DAFUN('SQRT',ISCRDA( 22+IDAA),ISCRDA( 24+IDAA))              
      CALL DAFUN('SIN ',ISCRDA( 23+IDAA),ISCRDA( 25+IDAA))              
      CALL DAMUL(ISCRDA( 24+IDAA),ISCRDA( 25+IDAA),ISCRDA( 26+IDAA))    
      CALL DACOP(ISCRDA( 26+IDAA),YY         ((2)))                     
*FOX  XX(2)=XX(2)+XX(1)*SIN(TX)*YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-       *FOX
*FOX  YY(2)*YY(2))/COS(ATAN(YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-           *FOX
*FOX  YY(1)*YY(1)-YY(2)*YY(2)))+TX) ;                                   *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(YY         ((1)),ISCRDA( 12+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))              
      CALL DADIV(ISCRDA( 12+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 15+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 16+IDAA))          
      CALL DACOP(YY         ((2)),ISCRDA( 17+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA( 18+IDAA))                     
      CALL DAFUN('ATAN',ISCRDA( 14+IDAA),ISCRDA( 19+IDAA))              
      CALL DAMUL(ISCRDA( 15+IDAA),ISCRDA( 16+IDAA),ISCRDA( 20+IDAA))    
      CALL DAMUL(ISCRDA( 17+IDAA),ISCRDA( 18+IDAA),ISCRDA( 21+IDAA))    
      CALL DASUB(ISCRDA( 20+IDAA),ISCRDA( 21+IDAA),ISCRDA( 22+IDAA))    
      CALL DACAD(ISCRDA( 19+IDAA),ONE*TX         ,ISCRDA( 23+IDAA))     
      CALL DACOP(XX         ((2)),ISCRDA( 24+IDAA))                     
      CALL DACOP(XX         ((1)),ISCRDA( 25+IDAA))                     
      RSCRRI( 26+IDAA) = SIN (TX         )                              
      CALL DACOP(YY         ((2)),ISCRDA( 27+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 22+IDAA),ISCRDA( 28+IDAA))              
      CALL DAFUN('COS ',ISCRDA( 23+IDAA),ISCRDA( 29+IDAA))              
      CALL DACMU(ISCRDA( 25+IDAA),ONE*RSCRRI( 26+IDAA),ISCRDA( 30+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA( 30+IDAA),ISCRDA( 27+IDAA),ISCRDA( 31+IDAA))    
      CALL DADIV(ISCRDA( 31+IDAA),ISCRDA( 28+IDAA),ISCRDA( 32+IDAA))    
      CALL DADIV(ISCRDA( 32+IDAA),ISCRDA( 29+IDAA),ISCRDA( 33+IDAA))    
      CALL DAADD(ISCRDA( 24+IDAA),ISCRDA( 33+IDAA),ISCRDA( 34+IDAA))    
      CALL DACOP(ISCRDA( 34+IDAA),XX         ((2)))                     
*FOX  XX(1)=XX(1)*(COS(TX)+SIN(TX)*TAN(ATAN(YY(1)/SQRT((ONE+DPDA)*      *FOX
*FOX  (ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))+TX)) ;                       *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(YY         ((1)),ISCRDA( 12+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))              
      CALL DADIV(ISCRDA( 12+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DAFUN('ATAN',ISCRDA( 14+IDAA),ISCRDA( 15+IDAA))              
      CALL DACAD(ISCRDA( 15+IDAA),ONE*TX         ,ISCRDA( 16+IDAA))     
      RSCRRI( 17+IDAA) = COS (TX         )                              
      RSCRRI( 18+IDAA) = SIN (TX         )                              
      CALL DAFUN('TAN ',ISCRDA( 16+IDAA),ISCRDA( 19+IDAA))              
      CALL DACMU(ISCRDA( 19+IDAA),ONE*RSCRRI( 18+IDAA),ISCRDA( 20+IDAA))
     *                                                                  
      CALL DACAD(ISCRDA( 20+IDAA),ONE*RSCRRI( 17+IDAA),ISCRDA( 21+IDAA))
     *                                                                  
      CALL DACOP(XX         ((1)),ISCRDA( 22+IDAA))                     
      CALL DAMUL(ISCRDA( 22+IDAA),ISCRDA( 21+IDAA),ISCRDA( 23+IDAA))    
      CALL DACOP(ISCRDA( 23+IDAA),XX         ((1)))                     
*FOX  YY(1)=SQRT((ONE+DPDA)*(ONE+DPDA)-YY(2)*YY(2))*SIN(ATAN(YY(1)/     *FOX
*FOX  SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))+TX) ;        *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(YY         ((1)),ISCRDA( 12+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))              
      CALL DADIV(ISCRDA( 12+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 15+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA( 16+IDAA))          
      CALL DACOP(YY         ((2)),ISCRDA( 17+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA( 18+IDAA))                     
      CALL DAFUN('ATAN',ISCRDA( 14+IDAA),ISCRDA( 19+IDAA))              
      CALL DAMUL(ISCRDA( 15+IDAA),ISCRDA( 16+IDAA),ISCRDA( 20+IDAA))    
      CALL DAMUL(ISCRDA( 17+IDAA),ISCRDA( 18+IDAA),ISCRDA( 21+IDAA))    
      CALL DASUB(ISCRDA( 20+IDAA),ISCRDA( 21+IDAA),ISCRDA( 22+IDAA))    
      CALL DACAD(ISCRDA( 19+IDAA),ONE*TX         ,ISCRDA( 23+IDAA))     
      CALL DAFUN('SQRT',ISCRDA( 22+IDAA),ISCRDA( 24+IDAA))              
      CALL DAFUN('SIN ',ISCRDA( 23+IDAA),ISCRDA( 25+IDAA))              
      CALL DAMUL(ISCRDA( 24+IDAA),ISCRDA( 25+IDAA),ISCRDA( 26+IDAA))    
      CALL DACOP(ISCRDA( 26+IDAA),YY         ((1)))                     
!     CALL SHIFT(-EMBL*TAN(TX),-EMBL*TAN(TY)/COS(TX))
*FOX  XX(1)=XX(1)+EMBL*TAN(TX) ;                                        *FOX
      CALL DACOP(XX         ((1)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = TAN (TX         )                              
      RSCRRI(  3+IDAA) = EMBL        * RSCRRI(  2+IDAA)                 
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  4+IDAA),XX         ((1)))                     
*FOX  XX(2)=XX(2)+EMBL*TAN(TY)/COS(TX) ;                                *FOX
      CALL DACOP(XX         ((2)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = TAN (TY         )                              
      RSCRRI(  3+IDAA) = COS (TX         )                              
      RSCRRI(  4+IDAA) = EMBL        * RSCRRI(  2+IDAA)                 
      RSCRRI(  5+IDAA) = RSCRRI(  4+IDAA) / RSCRRI(  3+IDAA)            
      CALL DACAD(ISCRDA(  1+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA(  6+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  6+IDAA),XX         ((2)))                     
!     CALL DRIFT(-EMBL/2)
*FOX  XX(1)=XX(1)-EMBL/TWO*YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-            *FOX
*FOX  YY(1)*YY(1)-YY(2)*YY(2)) ;                                        *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(XX         ((1)),ISCRDA( 12+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA( 13+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 14+IDAA))              
      RSCRRI( 15+IDAA) = EMBL        / TWO                              
      CALL DACMU(ISCRDA( 13+IDAA),ONE*RSCRRI( 15+IDAA),ISCRDA( 16+IDAA))
     *                                                                  
      CALL DADIV(ISCRDA( 16+IDAA),ISCRDA( 14+IDAA),ISCRDA( 17+IDAA))    
      CALL DASUB(ISCRDA( 12+IDAA),ISCRDA( 17+IDAA),ISCRDA( 18+IDAA))    
      CALL DACOP(ISCRDA( 18+IDAA),XX         ((1)))                     
*FOX  XX(2)=XX(2)-EMBL/TWO*YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-            *FOX
*FOX  YY(1)*YY(1)-YY(2)*YY(2)) ;                                        *FOX
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACAD(DPDA       ,ONE*ONE        ,ISCRDA(  2+IDAA))          
      CALL DACOP(YY         ((1)),ISCRDA(  3+IDAA))                     
      CALL DACOP(YY         ((1)),ISCRDA(  4+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA(  6+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  7+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(XX         ((2)),ISCRDA( 12+IDAA))                     
      CALL DACOP(YY         ((2)),ISCRDA( 13+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 11+IDAA),ISCRDA( 14+IDAA))              
      RSCRRI( 15+IDAA) = EMBL        / TWO                              
      CALL DACMU(ISCRDA( 13+IDAA),ONE*RSCRRI( 15+IDAA),ISCRDA( 16+IDAA))
     *                                                                  
      CALL DADIV(ISCRDA( 16+IDAA),ISCRDA( 14+IDAA),ISCRDA( 17+IDAA))    
      CALL DASUB(ISCRDA( 12+IDAA),ISCRDA( 17+IDAA),ISCRDA( 18+IDAA))    
      CALL DACOP(ISCRDA( 18+IDAA),XX         ((2)))                     
       
*FOX  XX(1)=XX(1)*C1E3;                                                 *FOX
      CALL DACOP(XX         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1E3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),XX         ((1)))                     
*FOX  XX(2)=XX(2)*C1E3;                                                 *FOX
      CALL DACOP(XX         ((2)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1E3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),XX         ((2)))                     
*FOX  YY(1)=YY(1)*C1E3;                                                 *FOX
      CALL DACOP(YY         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1E3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),YY         ((1)))                     
*FOX  YY(2)=YY(2)*C1E3;                                                 *FOX
      CALL DACOP(YY         ((2)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*C1E3       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),YY         ((2)))                     
       
        CALL DADAL(YI      ,1)                                                  
        CALL DADAL(XI      ,1)                                                  
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
      double precision al,as,chi0,chid,dp1,dps,exz,sigm
      common/syos/as(6,2,npart,nele),al(6,2,npart,nele),sigm(mpa),      &
     &dps(mpa),idz(2)
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
*FOX{
      INTEGER SIGMDA  
      INTEGER DPDA    
      INTEGER DPDA1   
      INTEGER RV      
      INTEGER XX      (2)
      INTEGER YY      (2)
      INTEGER EJ1     
      INTEGER EJF1    
      INTEGER ALDA    (2,6)
      INTEGER ASDA    (2,6)
      INTEGER ALDAQ   (2,6)
      INTEGER ASDAQ   (2,6)
      INTEGER SMIDA   (MCOR)
      INTEGER ISCRDA, ISCRRI,IDAO
      DOUBLE PRECISION RSCRRI
      COMMON/DASCR/ISCRDA(100),RSCRRI(100),ISCRRI(100),IDAO
      if(1.eq.1) then                                                 
         CALL DAKEY('FOX V2.1')
      ENDIF
      IDAA = IDAO
*FOX}
!-----------------------------------------------------------------------
      ix=ixcav
      if(kz(ix).eq.12) then
*FOX  EJ1=EJ1+ED(IX)*SIN(HSYC(IX)*SIGMDA/C1E3*                          *FOX
*FOX  ITIONC(IX)+PHASC(IX)) ;
      RSCRRI(  1+IDAA) = HSYC       (IX         )                       
      ISCRRI(  2+IDAA) = ITIONC     (IX         )                       
      RSCRRI(  3+IDAA) = PHASC      (IX         )                       
      CALL DACMU(SIGMDA     ,ONE*RSCRRI(  1+IDAA),ISCRDA(  4+IDAA))     
      CALL DACDI(ISCRDA(  4+IDAA),ONE*C1E3       ,ISCRDA(  5+IDAA))     
      CALL DACMU(ISCRDA(  5+IDAA),ONE*ISCRRI(  2+IDAA),ISCRDA(  6+IDAA))
     *                                                                  
      CALL DACAD(ISCRDA(  6+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      RSCRRI(  8+IDAA) = ED         (IX         )                       
      CALL DAFUN('SIN ',ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))              
      CALL DACMU(ISCRDA(  9+IDAA),ONE*RSCRRI(  8+IDAA),ISCRDA( 10+IDAA))
     *                                                                  
      CALL DAADD(EJ1        ,ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))         
      CALL DACOP(ISCRDA( 11+IDAA),EJ1        )                          
      else
*FOX  EJ1=EJ1+HSY(1)*SIN(HSY(3)*SIGMDA/C1E3*ITION+PHAS) ;               *FOX
      RSCRRI(  1+IDAA) = HSY        ((3))                               
      CALL DACMU(SIGMDA     ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACDI(ISCRDA(  2+IDAA),ONE*C1E3       ,ISCRDA(  3+IDAA))     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*ITION      ,ISCRDA(  4+IDAA))     
      CALL DACAD(ISCRDA(  4+IDAA),ONE*PHAS       ,ISCRDA(  5+IDAA))     
      RSCRRI(  6+IDAA) = HSY        ((1))                               
      CALL DAFUN('SIN ',ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))              
      CALL DACMU(ISCRDA(  7+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAADD(EJ1        ,ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))         
      CALL DACOP(ISCRDA(  9+IDAA),EJ1        )                          
      endif
*FOX  EJF1=SQRT(EJ1*EJ1-PMA*PMA) ;                                      *FOX
      CALL DAMUL(EJ1        ,EJ1        ,ISCRDA(  1+IDAA))              
      RSCRRI(  2+IDAA) = PMA         * PMA                              
      CALL DACSU(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DAFUN('SQRT',ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))              
      CALL DACOP(ISCRDA(  4+IDAA),EJF1       )                          
*FOX  DPDA1=(EJF1-E0F)/E0F*C1E3 ;                                       *FOX
      CALL DACSU(EJF1       ,ONE*E0F        ,ISCRDA(  1+IDAA))          
      CALL DACDI(ISCRDA(  1+IDAA),ONE*E0F        ,ISCRDA(  2+IDAA))     
      CALL DACMU(ISCRDA(  2+IDAA),ONE*C1E3       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),DPDA1      )                          
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
!   XX, YY    (REAL)    ARGUMENT TO CERF.                              **FOX
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
      double precision al,as,chi0,chid,dp1,dps,exz,sigm
      common/syos/as(6,2,npart,nele),al(6,2,npart,nele),sigm(mpa),      &
     &dps(mpa),idz(2)
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
*FOX{
      INTEGER XX      
      INTEGER YY      
      INTEGER WX      
      INTEGER WY      
      INTEGER X       
      INTEGER Y       
      INTEGER Q       
      INTEGER H       
      INTEGER XH      
      INTEGER YH      
      INTEGER RX      (33)
      INTEGER RY      (33)
      INTEGER TX      
      INTEGER TN      
      INTEGER TY      
      INTEGER SAUX    
      INTEGER SX      
      INTEGER SY      
      INTEGER XL      
      INTEGER ISCRDA, ISCRRI,IDAO
      DOUBLE PRECISION RSCRRI
      COMMON/DASCR/ISCRDA(100),RSCRRI(100),ISCRRI(100),IDAO
      if(1.eq.1) then                                                 
         CALL DAKEY('FOX V2.1')
         CALL DAALL(X       ,1,'X         ',NORD,NVAR)
         CALL DAALL(Y       ,1,'Y         ',NORD,NVAR)
         CALL DAALL(Q       ,1,'Q         ',NORD,NVAR)
         CALL DAALL(H       ,1,'H         ',NORD,NVAR)
         CALL DAALL(XH      ,1,'XH        ',NORD,NVAR)
         CALL DAALL(YH      ,1,'YH        ',NORD,NVAR)
         CALL DAALL(RX      ,1*(33),'RX        ',NORD,NVAR)
         CALL DAALL(RY      ,1*(33),'RY        ',NORD,NVAR)
         CALL DAALL(TX      ,1,'TX        ',NORD,NVAR)
         CALL DAALL(TN      ,1,'TN        ',NORD,NVAR)
         CALL DAALL(TY      ,1,'TY        ',NORD,NVAR)
         CALL DAALL(SAUX    ,1,'SAUX      ',NORD,NVAR)
         CALL DAALL(SX      ,1,'SX        ',NORD,NVAR)
         CALL DAALL(SY      ,1,'SY        ',NORD,NVAR)
         CALL DAALL(XL      ,1,'XL        ',NORD,NVAR)
      ENDIF
      IDAA = IDAO
*FOX}
!-----------------------------------------------------------------------
*FOX  X=XX ;                                                            *FOX
      CALL DACOP(XX         ,X          )                               
*FOX  Y=YY ;                                                            *FOX
      CALL DACOP(YY         ,Y          )                               
      if(dare(x).lt.zero) then
        write(*,*) ' Problem in DA complex error function: dare(x) < 0'
*FOX    X=-X ;                                                          *FOX
      CALL DACMU(X          ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),X          )                          
      endif
      if(dare(y).lt.zero) then
        write(*,*) ' Problem in DA complex error function: dare(y) < 0'
*FOX    Y=-Y ;                                                          *FOX
      CALL DACMU(Y          ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),Y          )                          
      endif
      if(dare(y).lt.ylim.and.dare(x).lt.xlim) then
*FOX    Q=(ONE-Y/YLIM)*SQRT(ONE-X*X/XLIM/XLIM) ;                        *FOX
      CALL DACDI(Y          ,ONE*YLIM       ,ISCRDA(  1+IDAA))          
      CALL DAMUL(X          ,X          ,ISCRDA(  2+IDAA))              
      CALL DACDI(ISCRDA(  2+IDAA),ONE*XLIM       ,ISCRDA(  3+IDAA))     
      CALL DACDI(ISCRDA(  3+IDAA),ONE*XLIM       ,ISCRDA(  4+IDAA))     
      CALL DASUC(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  5+IDAA))     
      CALL DASUC(ISCRDA(  4+IDAA),ONE*ONE        ,ISCRDA(  6+IDAA))     
      CALL DAFUN('SQRT',ISCRDA(  6+IDAA),ISCRDA(  7+IDAA))              
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),Q          )                          
*FOX    DUM=3.2D0 ;
      DUM         = (3.2D0      )                                       
*FOX    H=ONE/(DUM*Q) ;                                                 *FOX
      CALL DACMU(Q          ,ONE*DUM        ,ISCRDA(  1+IDAA))          
      CALL DADIC(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),H          )                          
        nc=7+int(23.0d0*dare(q))
*FOX    XL=EXP((1-NC)*LOG(H)) ;                                         *FOX
      ISCRRI(  1+IDAA) = (1          ) - NC                             
      CALL DAFUN('LOG ',H          ,ISCRDA(  2+IDAA))                   
      CALL DACMU(ISCRDA(  2+IDAA),ONE*ISCRRI(  1+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DAFUN('EXP ',ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))              
      CALL DACOP(ISCRDA(  4+IDAA),XL         )                          
*FOX    XH=Y+HALF/H ;                                                   *FOX
      CALL DADIC(H          ,ONE*HALF       ,ISCRDA(  1+IDAA))          
      CALL DAADD(Y          ,ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),XH         )                          
*FOX    YH=X ;                                                          *FOX
      CALL DACOP(X          ,YH         )                               
        nuu=10+int(21.0d0*dare(q))
        nuu1=nuu+1
*FOX    RX(NUU1)=ZERO ;                                                 *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(RX         (NUU1       ),RSCRRI(100))                  
*FOX    RY(NUU1)=ZERO ;                                                 *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(RY         (NUU1       ),RSCRRI(100))                  
        do 10 n=nuu,1,-1
          n1=n+1
*FOX      TX=XH+N*RX(N1) ;                                              *FOX
      CALL DACOP(RX         (N1         ),ISCRDA(  1+IDAA))             
      CALL DACMU(ISCRDA(  1+IDAA),ONE*N          ,ISCRDA(  2+IDAA))     
      CALL DAADD(XH         ,ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),TX         )                          
*FOX      TY=YH-N*RY(N1) ;                                              *FOX
      CALL DACOP(RY         (N1         ),ISCRDA(  1+IDAA))             
      CALL DACMU(ISCRDA(  1+IDAA),ONE*N          ,ISCRDA(  2+IDAA))     
      CALL DASUB(YH         ,ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),TY         )                          
*FOX      TN=TX*TX+TY*TY ;                                              *FOX
      CALL DAMUL(TX         ,TX         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(TY         ,TY         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),TN         )                          
*FOX      RX(N)=HALF*TX/TN ;                                            *FOX
      CALL DACMU(TX         ,ONE*HALF       ,ISCRDA(  1+IDAA))          
      CALL DADIV(ISCRDA(  1+IDAA),TN         ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),RX         (N          ))             
*FOX      RY(N)=HALF*TY/TN ;                                            *FOX
      CALL DACMU(TY         ,ONE*HALF       ,ISCRDA(  1+IDAA))          
      CALL DADIV(ISCRDA(  1+IDAA),TN         ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),RY         (N          ))             
   10   continue
*FOX    SX=ZERO ;                                                       *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(SX         ,RSCRRI(100))                               
*FOX    SY=ZERO ;                                                       *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(SY         ,RSCRRI(100))                               
        do 20 n=nc,1,-1
*FOX      SAUX=SX+XL ;                                                  *FOX
      CALL DAADD(SX         ,XL         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),SAUX       )                          
*FOX      SX=RX(N)*SAUX-RY(N)*SY ;                                      *FOX
      CALL DACOP(RX         (N          ),ISCRDA(  1+IDAA))             
      CALL DACOP(RY         (N          ),ISCRDA(  2+IDAA))             
      CALL DAMUL(ISCRDA(  1+IDAA),SAUX       ,ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),SY         ,ISCRDA(  4+IDAA))         
      CALL DASUB(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(ISCRDA(  5+IDAA),SX         )                          
*FOX      SY=RX(N)*SY+RY(N)*SAUX ;                                      *FOX
      CALL DACOP(RX         (N          ),ISCRDA(  1+IDAA))             
      CALL DACOP(RY         (N          ),ISCRDA(  2+IDAA))             
      CALL DAMUL(ISCRDA(  1+IDAA),SY         ,ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),SAUX       ,ISCRDA(  4+IDAA))         
      CALL DAADD(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(ISCRDA(  5+IDAA),SY         )                          
*FOX      XL=H*XL ;                                                     *FOX
      CALL DAMUL(H          ,XL         ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),XL         )                          
   20   continue
*FOX    WX=CC*SX ;                                                      *FOX
      CALL DACMU(SX         ,ONE*CC         ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),WX         )                          
*FOX    WY=CC*SY ;                                                      *FOX
      CALL DACMU(SY         ,ONE*CC         ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),WY         )                          
      else
*FOX    XH=Y ;                                                          *FOX
      CALL DACOP(Y          ,XH         )                               
*FOX    YH=X ;                                                          *FOX
      CALL DACOP(X          ,YH         )                               
*FOX    RX(1)=ZERO ;                                                    *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(RX         ((1)),RSCRRI(100))                          
*FOX    RY(1)=ZERO ;                                                    *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(RY         ((1)),RSCRRI(100))                          
        do 30 n=9,1,-1
*FOX      TX=XH+N*RX(1) ;                                               *FOX
      CALL DACOP(RX         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*N          ,ISCRDA(  2+IDAA))     
      CALL DAADD(XH         ,ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),TX         )                          
*FOX      TY=YH-N*RY(1) ;                                               *FOX
      CALL DACOP(RY         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*N          ,ISCRDA(  2+IDAA))     
      CALL DASUB(YH         ,ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),TY         )                          
*FOX      TN=TX*TX+TY*TY ;                                              *FOX
      CALL DAMUL(TX         ,TX         ,ISCRDA(  1+IDAA))              
      CALL DAMUL(TY         ,TY         ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),TN         )                          
*FOX      RX(1)=HALF*TX/TN ;                                            *FOX
      CALL DACMU(TX         ,ONE*HALF       ,ISCRDA(  1+IDAA))          
      CALL DADIV(ISCRDA(  1+IDAA),TN         ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),RX         ((1)))                     
*FOX      RY(1)=HALF*TY/TN ;                                            *FOX
      CALL DACMU(TY         ,ONE*HALF       ,ISCRDA(  1+IDAA))          
      CALL DADIV(ISCRDA(  1+IDAA),TN         ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),RY         ((1)))                     
   30   continue
*FOX    WX=CC*RX(1) ;                                                   *FOX
      CALL DACOP(RX         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*CC         ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),WX         )                          
*FOX    WY=CC*RY(1) ;                                                   *FOX
      CALL DACOP(RY         ((1)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*CC         ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),WY         )                          
      endif
!      if(dare(y).eq.0.) then
!*FOX    WX=EXP(-X*X) ;                                                 *FOX
!      endif
!hr05 if(dare(yy).lt.0.) then
      if(dare(yy).lt.0.d0) then                                          !hr05
*FOX    WX=TWO*EXP(Y*Y-X*X)*COS(TWO*X*Y)-WX ;                           *FOX
      CALL DAMUL(Y          ,Y          ,ISCRDA(  1+IDAA))              
      CALL DAMUL(X          ,X          ,ISCRDA(  2+IDAA))              
      CALL DACMU(X          ,ONE*TWO        ,ISCRDA(  3+IDAA))          
      CALL DAMUL(ISCRDA(  3+IDAA),Y          ,ISCRDA(  4+IDAA))         
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DAFUN('EXP ',ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))              
      CALL DAFUN('COS ',ISCRDA(  4+IDAA),ISCRDA(  7+IDAA))              
      CALL DACMU(ISCRDA(  6+IDAA),ONE*TWO        ,ISCRDA(  8+IDAA))     
      CALL DAMUL(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DASUB(ISCRDA(  9+IDAA),WX         ,ISCRDA( 10+IDAA))         
      CALL DACOP(ISCRDA( 10+IDAA),WX         )                          
*FOX    WY=-TWO*EXP(Y*Y-X*X)*SIN(TWO*X*Y)-WY ;                          *FOX
      CALL DAMUL(Y          ,Y          ,ISCRDA(  1+IDAA))              
      CALL DAMUL(X          ,X          ,ISCRDA(  2+IDAA))              
      CALL DACMU(X          ,ONE*TWO        ,ISCRDA(  3+IDAA))          
      CALL DAMUL(ISCRDA(  3+IDAA),Y          ,ISCRDA(  4+IDAA))         
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DAFUN('EXP ',ISCRDA(  5+IDAA),ISCRDA(  6+IDAA))              
      CALL DAFUN('SIN ',ISCRDA(  4+IDAA),ISCRDA(  7+IDAA))              
      RSCRRI(  8+IDAA) = (-ONE       ) * TWO                            
      CALL DACMU(ISCRDA(  6+IDAA),ONE*RSCRRI(  8+IDAA),ISCRDA(  9+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA(  9+IDAA),ISCRDA(  7+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),WY         ,ISCRDA( 11+IDAA))         
      CALL DACOP(ISCRDA( 11+IDAA),WY         )                          
!hr05   if(dare(xx).gt.0.) then
        if(dare(xx).gt.0.d0) then                                        !hr05
*FOX      WY=-WY ;                                                      *FOX
      CALL DACMU(WY         ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),WY         )                          
        endif
      else
!hr05   if(dare(xx).lt.0.) then
        if(dare(xx).lt.0.d0) then                                        !hr05
*FOX      WY=-WY ;                                                      *FOX
      CALL DACMU(WY         ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),WY         )                          
        endif
      endif
        CALL DADAL(XL      ,1)                                                  
        CALL DADAL(SY      ,1)                                                  
        CALL DADAL(SX      ,1)                                                  
        CALL DADAL(SAUX    ,1)                                                  
        CALL DADAL(TY      ,1)                                                  
        CALL DADAL(TN      ,1)                                                  
        CALL DADAL(TX      ,1)                                                  
        CALL DADAL(RY      ,1*(33))                                             
        CALL DADAL(RX      ,1*(33))                                             
        CALL DADAL(YH      ,1)                                                  
        CALL DADAL(XH      ,1)                                                  
        CALL DADAL(H       ,1)                                                  
        CALL DADAL(Q       ,1)                                                  
        CALL DADAL(Y       ,1)                                                  
        CALL DADAL(X       ,1)                                                  
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
*FOX{
      INTEGER TRACK   (6)
      INTEGER ISCRDA, ISCRRI,IDAO
      DOUBLE PRECISION RSCRRI
      COMMON/DASCR/ISCRDA(100),RSCRRI(100),ISCRRI(100),IDAO
      if(1.eq.1) then                                                 
         CALL DAKEY('FOX V2.1')
      ENDIF
      IDAA = IDAO
*FOX}
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
*FOX{
      INTEGER TRACK   (6)
      INTEGER A       
      INTEGER H       
      INTEGER SQR1A   
      INTEGER A1      
      INTEGER HD1     
      INTEGER H1X     
      INTEGER H1Y     
      INTEGER H1Z     
      INTEGER X1      
      INTEGER Y1      
      INTEGER DET     
      INTEGER H1      
      INTEGER ISCRDA, ISCRRI,IDAO
      DOUBLE PRECISION RSCRRI
      COMMON/DASCR/ISCRDA(100),RSCRRI(100),ISCRRI(100),IDAO
      if(1.eq.1) then                                                 
         CALL DAKEY('FOX V2.1')
         CALL DAALL(A       ,1,'A         ',NORD,NVAR)
         CALL DAALL(H       ,1,'H         ',NORD,NVAR)
         CALL DAALL(SQR1A   ,1,'SQR1A     ',NORD,NVAR)
         CALL DAALL(A1      ,1,'A1        ',NORD,NVAR)
         CALL DAALL(HD1     ,1,'HD1       ',NORD,NVAR)
         CALL DAALL(H1X     ,1,'H1X       ',NORD,NVAR)
         CALL DAALL(H1Y     ,1,'H1Y       ',NORD,NVAR)
         CALL DAALL(H1Z     ,1,'H1Z       ',NORD,NVAR)
         CALL DAALL(X1      ,1,'X1        ',NORD,NVAR)
         CALL DAALL(Y1      ,1,'Y1        ',NORD,NVAR)
         CALL DAALL(DET     ,1,'DET       ',NORD,NVAR)
         CALL DAALL(H1      ,1,'H1        ',NORD,NVAR)
      ENDIF
      IDAA = IDAO
*FOX}
!-----------------------------------------------------------------------
*FOX    H=TRACK(6)+ONE-SQRT((ONE+TRACK(6))*(ONE+TRACK(6))-              *FOX
*FOX    TRACK(2)*TRACK(2)-TRACK(4)*TRACK(4)) ;                          *FOX
      CALL DACOP(TRACK      ((6)),ISCRDA(  1+IDAA))                     
      CALL DACOP(TRACK      ((6)),ISCRDA(  2+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  3+IDAA))     
      CALL DACAD(ISCRDA(  2+IDAA),ONE*ONE        ,ISCRDA(  4+IDAA))     
      CALL DACOP(TRACK      ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(TRACK      ((2)),ISCRDA(  6+IDAA))                     
      CALL DACOP(TRACK      ((4)),ISCRDA(  7+IDAA))                     
      CALL DACOP(TRACK      ((4)),ISCRDA(  8+IDAA))                     
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  9+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA( 10+IDAA))    
      CALL DAMUL(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 11+IDAA))    
      CALL DASUB(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 12+IDAA))    
      CALL DASUB(ISCRDA( 12+IDAA),ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))    
      CALL DACOP(TRACK      ((6)),ISCRDA( 14+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 13+IDAA),ISCRDA( 15+IDAA))              
      CALL DACAD(ISCRDA( 14+IDAA),ONE*ONE        ,ISCRDA( 16+IDAA))     
      CALL DASUB(ISCRDA( 16+IDAA),ISCRDA( 15+IDAA),ISCRDA( 17+IDAA))    
      CALL DACOP(ISCRDA( 17+IDAA),H          )                          
*FOX    TRACK(6)=TRACK(6)-CALPHA*TPHI*TRACK(2)                          *FOX
*FOX              -TRACK(4)*SALPHA*TPHI+H*TPHI*TPHI ;                   *FOX
      CALL DACOP(TRACK      ((6)),ISCRDA(  1+IDAA))                     
      CALL DACOP(TRACK      ((2)),ISCRDA(  2+IDAA))                     
      CALL DACOP(TRACK      ((4)),ISCRDA(  3+IDAA))                     
      RSCRRI(  4+IDAA) = CALPHA      * TPHI                             
      CALL DACMU(ISCRDA(  2+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  3+IDAA),ONE*SALPHA     ,ISCRDA(  6+IDAA))     
      CALL DACMU(ISCRDA(  6+IDAA),ONE*TPHI       ,ISCRDA(  7+IDAA))     
      CALL DACMU(H          ,ONE*TPHI       ,ISCRDA(  8+IDAA))          
      CALL DACMU(ISCRDA(  8+IDAA),ONE*TPHI       ,ISCRDA(  9+IDAA))     
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  5+IDAA),ISCRDA( 10+IDAA))    
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA(  7+IDAA),ISCRDA( 11+IDAA))    
      CALL DAADD(ISCRDA( 11+IDAA),ISCRDA(  9+IDAA),ISCRDA( 12+IDAA))    
      CALL DACOP(ISCRDA( 12+IDAA),TRACK      ((6)))                     
*FOX    TRACK(2)=(TRACK(2)-TPHI*H*CALPHA)/CPHI ;                        *FOX
      CALL DACOP(TRACK      ((2)),ISCRDA(  1+IDAA))                     
      CALL DACMU(H          ,ONE*TPHI       ,ISCRDA(  2+IDAA))          
      CALL DACMU(ISCRDA(  2+IDAA),ONE*CALPHA     ,ISCRDA(  3+IDAA))     
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACDI(ISCRDA(  4+IDAA),ONE*CPHI       ,ISCRDA(  5+IDAA))     
      CALL DACOP(ISCRDA(  5+IDAA),TRACK      ((2)))                     
*FOX    TRACK(4)=(TRACK(4)-TPHI*H*SALPHA)/CPHI ;                        *FOX
      CALL DACOP(TRACK      ((4)),ISCRDA(  1+IDAA))                     
      CALL DACMU(H          ,ONE*TPHI       ,ISCRDA(  2+IDAA))          
      CALL DACMU(ISCRDA(  2+IDAA),ONE*SALPHA     ,ISCRDA(  3+IDAA))     
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACDI(ISCRDA(  4+IDAA),ONE*CPHI       ,ISCRDA(  5+IDAA))     
      CALL DACOP(ISCRDA(  5+IDAA),TRACK      ((4)))                     
*FOX    HD1=SQRT((ONE+TRACK(6))*(ONE+TRACK(6))-TRACK(2)*TRACK(2)-       *FOX
*FOX    TRACK(4)*TRACK(4)) ;                                            *FOX
      CALL DACOP(TRACK      ((6)),ISCRDA(  1+IDAA))                     
      CALL DACOP(TRACK      ((6)),ISCRDA(  2+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  3+IDAA))     
      CALL DACAD(ISCRDA(  2+IDAA),ONE*ONE        ,ISCRDA(  4+IDAA))     
      CALL DACOP(TRACK      ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(TRACK      ((2)),ISCRDA(  6+IDAA))                     
      CALL DACOP(TRACK      ((4)),ISCRDA(  7+IDAA))                     
      CALL DACOP(TRACK      ((4)),ISCRDA(  8+IDAA))                     
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  9+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA( 10+IDAA))    
      CALL DAMUL(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 11+IDAA))    
      CALL DASUB(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 12+IDAA))    
      CALL DASUB(ISCRDA( 12+IDAA),ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))    
      CALL DAFUN('SQRT',ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))              
      CALL DACOP(ISCRDA( 14+IDAA),HD1        )                          
*FOX    H1X=TRACK(2)/HD1 ;                                              *FOX
      CALL DACOP(TRACK      ((2)),ISCRDA(  1+IDAA))                     
      CALL DADIV(ISCRDA(  1+IDAA),HD1        ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),H1X        )                          
*FOX    H1Y=TRACK(4)/HD1 ;                                              *FOX
      CALL DACOP(TRACK      ((4)),ISCRDA(  1+IDAA))                     
      CALL DADIV(ISCRDA(  1+IDAA),HD1        ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),H1Y        )                          
*FOX    H1Z=ONE-(ONE+TRACK(6))/HD1 ;                                    *FOX
      CALL DACOP(TRACK      ((6)),ISCRDA(  1+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  2+IDAA))     
      CALL DADIV(ISCRDA(  2+IDAA),HD1        ,ISCRDA(  3+IDAA))         
      CALL DASUC(ISCRDA(  3+IDAA),ONE*ONE        ,ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),H1Z        )                          
*FOX    X1=CALPHA*TPHI*TRACK(5)+(ONE+CALPHA*SPHI*H1X)*TRACK(1)          *FOX
*FOX       +TRACK(3)*SALPHA*SPHI*H1X ;                                  *FOX
      RSCRRI(  1+IDAA) = CALPHA      * SPHI                             
      CALL DACMU(H1X        ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACAD(ISCRDA(  2+IDAA),ONE*ONE        ,ISCRDA(  3+IDAA))     
      CALL DACOP(TRACK      ((5)),ISCRDA(  4+IDAA))                     
      CALL DACOP(TRACK      ((1)),ISCRDA(  5+IDAA))                     
      CALL DACOP(TRACK      ((3)),ISCRDA(  6+IDAA))                     
      RSCRRI(  7+IDAA) = CALPHA      * TPHI                             
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  7+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  5+IDAA),ISCRDA(  9+IDAA))    
      CALL DACMU(ISCRDA(  6+IDAA),ONE*SALPHA     ,ISCRDA( 10+IDAA))     
      CALL DACMU(ISCRDA( 10+IDAA),ONE*SPHI       ,ISCRDA( 11+IDAA))     
      CALL DAMUL(ISCRDA( 11+IDAA),H1X        ,ISCRDA( 12+IDAA))         
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  9+IDAA),ISCRDA( 13+IDAA))    
      CALL DAADD(ISCRDA( 13+IDAA),ISCRDA( 12+IDAA),ISCRDA( 14+IDAA))    
      CALL DACOP(ISCRDA( 14+IDAA),X1         )                          
*FOX    Y1=SALPHA*TPHI*TRACK(5)+(ONE+SALPHA*SPHI*H1Y)*TRACK(3)          *FOX
*FOX       +TRACK(1)*CALPHA*SPHI*H1Y ;                                  *FOX
      RSCRRI(  1+IDAA) = SALPHA      * SPHI                             
      CALL DACMU(H1Y        ,ONE*RSCRRI(  1+IDAA),ISCRDA(  2+IDAA))     
      CALL DACAD(ISCRDA(  2+IDAA),ONE*ONE        ,ISCRDA(  3+IDAA))     
      CALL DACOP(TRACK      ((5)),ISCRDA(  4+IDAA))                     
      CALL DACOP(TRACK      ((3)),ISCRDA(  5+IDAA))                     
      CALL DACOP(TRACK      ((1)),ISCRDA(  6+IDAA))                     
      RSCRRI(  7+IDAA) = SALPHA      * TPHI                             
      CALL DACMU(ISCRDA(  4+IDAA),ONE*RSCRRI(  7+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  5+IDAA),ISCRDA(  9+IDAA))    
      CALL DACMU(ISCRDA(  6+IDAA),ONE*CALPHA     ,ISCRDA( 10+IDAA))     
      CALL DACMU(ISCRDA( 10+IDAA),ONE*SPHI       ,ISCRDA( 11+IDAA))     
      CALL DAMUL(ISCRDA( 11+IDAA),H1Y        ,ISCRDA( 12+IDAA))         
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  9+IDAA),ISCRDA( 13+IDAA))    
      CALL DAADD(ISCRDA( 13+IDAA),ISCRDA( 12+IDAA),ISCRDA( 14+IDAA))    
      CALL DACOP(ISCRDA( 14+IDAA),Y1         )                          
*FOX    TRACK(5)=TRACK(5)/CPHI+H1Z*(SPHI*CALPHA*TRACK(1)                *FOX
*FOX       +SPHI*SALPHA*TRACK(3)) ;                                     *FOX
      CALL DACOP(TRACK      ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(TRACK      ((3)),ISCRDA(  2+IDAA))                     
      RSCRRI(  3+IDAA) = SPHI        * CALPHA                           
      CALL DACMU(ISCRDA(  1+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  4+IDAA))
     *                                                                  
      RSCRRI(  5+IDAA) = SPHI        * SALPHA                           
      CALL DACMU(ISCRDA(  2+IDAA),ONE*RSCRRI(  5+IDAA),ISCRDA(  6+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  4+IDAA),ISCRDA(  6+IDAA),ISCRDA(  7+IDAA))    
      CALL DACOP(TRACK      ((5)),ISCRDA(  8+IDAA))                     
      CALL DACDI(ISCRDA(  8+IDAA),ONE*CPHI       ,ISCRDA(  9+IDAA))     
      CALL DAMUL(H1Z        ,ISCRDA(  7+IDAA),ISCRDA( 10+IDAA))         
      CALL DAADD(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(ISCRDA( 11+IDAA),TRACK      ((5)))                     
*FOX    TRACK(1)=X1 ;                                                   *FOX
      CALL DACOP(X1         ,TRACK      ((1)))                          
*FOX    TRACK(3)=Y1 ;                                                   *FOX
      CALL DACOP(Y1         ,TRACK      ((3)))                          
        CALL DADAL(H1      ,1)                                                  
        CALL DADAL(DET     ,1)                                                  
        CALL DADAL(Y1      ,1)                                                  
        CALL DADAL(X1      ,1)                                                  
        CALL DADAL(H1Z     ,1)                                                  
        CALL DADAL(H1Y     ,1)                                                  
        CALL DADAL(H1X     ,1)                                                  
        CALL DADAL(HD1     ,1)                                                  
        CALL DADAL(A1      ,1)                                                  
        CALL DADAL(SQR1A   ,1)                                                  
        CALL DADAL(H       ,1)                                                  
        CALL DADAL(A       ,1)                                                  
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
*FOX{
      INTEGER TRACK   (6)
      INTEGER S       
      INTEGER SP      
      INTEGER DD      
      INTEGER DUM     (13)
      INTEGER SX      
      INTEGER SY      
      INTEGER SEPX    
      INTEGER SEPY    
      INTEGER SEPX0   
      INTEGER SEPY0   
      INTEGER COSTH   
      INTEGER SINTH   
      INTEGER COSTHP  
      INTEGER SINTHP  
      INTEGER BBFX    
      INTEGER BBFY    
      INTEGER BBF0    
      INTEGER BBGX    
      INTEGER BBGY    
      INTEGER ISCRDA, ISCRRI,IDAO
      DOUBLE PRECISION RSCRRI
      COMMON/DASCR/ISCRDA(100),RSCRRI(100),ISCRRI(100),IDAO
      if(1.eq.1) then                                                 
         CALL DAKEY('FOX V2.1')
         CALL DAALL(S       ,1,'S         ',NORD,NVAR)
         CALL DAALL(SP      ,1,'SP        ',NORD,NVAR)
         CALL DAALL(DD      ,1,'DD        ',NORD,NVAR)
         CALL DAALL(DUM     ,1*(13),'DUM       ',NORD,NVAR)
         CALL DAALL(SX      ,1,'SX        ',NORD,NVAR)
         CALL DAALL(SY      ,1,'SY        ',NORD,NVAR)
         CALL DAALL(SEPX    ,1,'SEPX      ',NORD,NVAR)
         CALL DAALL(SEPY    ,1,'SEPY      ',NORD,NVAR)
         CALL DAALL(SEPX0   ,1,'SEPX0     ',NORD,NVAR)
         CALL DAALL(SEPY0   ,1,'SEPY0     ',NORD,NVAR)
         CALL DAALL(COSTH   ,1,'COSTH     ',NORD,NVAR)
         CALL DAALL(SINTH   ,1,'SINTH     ',NORD,NVAR)
         CALL DAALL(COSTHP  ,1,'COSTHP    ',NORD,NVAR)
         CALL DAALL(SINTHP  ,1,'SINTHP    ',NORD,NVAR)
         CALL DAALL(BBFX    ,1,'BBFX      ',NORD,NVAR)
         CALL DAALL(BBFY    ,1,'BBFY      ',NORD,NVAR)
         CALL DAALL(BBF0    ,1,'BBF0      ',NORD,NVAR)
         CALL DAALL(BBGX    ,1,'BBGX      ',NORD,NVAR)
         CALL DAALL(BBGY    ,1,'BBGY      ',NORD,NVAR)
      ENDIF
      IDAA = IDAO
*FOX}
!-----------------------------------------------------------------------
      do 2000 jsli=1,nsli
*FOX    S=(TRACK(5)-STAR(3,JSLI))*HALF ;                                *FOX
      CALL DACOP(TRACK      ((5)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = STAR       ((3),JSLI       )                   
      CALL DACSU(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  3+IDAA),ONE*HALF       ,ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),S          )                          
*FOX    SP=S/CPHI2 ;                                                    *FOX
      CALL DACDI(S          ,ONE*CPHI2      ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),SP         )                          
*FOX    DUM(1)=BCU(IBB,1)+TWO*BCU(IBB,4)*SP+BCU(IBB,6)*SP*SP ;          *FOX
      RSCRRI(  1+IDAA) = BCU        (IBB        ,(1))                   
      RSCRRI(  2+IDAA) = BCU        (IBB        ,(4))                   
      RSCRRI(  3+IDAA) = BCU        (IBB        ,(6))                   
      RSCRRI(  4+IDAA) = TWO         * RSCRRI(  2+IDAA)                 
      CALL DACMU(SP         ,ONE*RSCRRI(  4+IDAA),ISCRDA(  5+IDAA))     
      CALL DACMU(SP         ,ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))     
      CALL DAMUL(ISCRDA(  6+IDAA),SP         ,ISCRDA(  7+IDAA))         
      CALL DACAD(ISCRDA(  5+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),DUM        ((1)))                     
*FOX    DUM(2)=BCU(IBB,2)+TWO*BCU(IBB,9)*SP+BCU(IBB,10)*SP*SP ;         *FOX
      RSCRRI(  1+IDAA) = BCU        (IBB        ,(2))                   
      RSCRRI(  2+IDAA) = BCU        (IBB        ,(9))                   
      RSCRRI(  3+IDAA) = BCU        (IBB        ,(10))                  
      RSCRRI(  4+IDAA) = TWO         * RSCRRI(  2+IDAA)                 
      CALL DACMU(SP         ,ONE*RSCRRI(  4+IDAA),ISCRDA(  5+IDAA))     
      CALL DACMU(SP         ,ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))     
      CALL DAMUL(ISCRDA(  6+IDAA),SP         ,ISCRDA(  7+IDAA))         
      CALL DACAD(ISCRDA(  5+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),DUM        ((2)))                     
*FOX    DUM(3)=BCU(IBB,3)+(BCU(IBB,5)+BCU(IBB,7))*SP+                   *FOX
*FOX    BCU(IBB,8)*SP*SP ;                                              *FOX
      RSCRRI(  1+IDAA) = BCU        (IBB        ,(5))                   
      RSCRRI(  2+IDAA) = BCU        (IBB        ,(7))                   
      RSCRRI(  3+IDAA) = RSCRRI(  1+IDAA) + RSCRRI(  2+IDAA)            
      RSCRRI(  4+IDAA) = BCU        (IBB        ,(3))                   
      RSCRRI(  5+IDAA) = BCU        (IBB        ,(8))                   
      CALL DACMU(SP         ,ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))     
      CALL DACMU(SP         ,ONE*RSCRRI(  5+IDAA),ISCRDA(  7+IDAA))     
      CALL DAMUL(ISCRDA(  7+IDAA),SP         ,ISCRDA(  8+IDAA))         
      CALL DACAD(ISCRDA(  6+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA(  9+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  9+IDAA),ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))    
      CALL DACOP(ISCRDA( 10+IDAA),DUM        ((3)))                     
*FOX    DUM(4)=DUM(1)-DUM(2) ;                                          *FOX
      CALL DACOP(DUM        ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(DUM        ((2)),ISCRDA(  2+IDAA))                     
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),DUM        ((4)))                     
*FOX    DUM(5)=DUM(4)*DUM(4)+FOUR*DUM(3)*DUM(3) ;                       *FOX
      CALL DACOP(DUM        ((4)),ISCRDA(  1+IDAA))                     
      CALL DACOP(DUM        ((4)),ISCRDA(  2+IDAA))                     
      CALL DACOP(DUM        ((3)),ISCRDA(  3+IDAA))                     
      CALL DACOP(DUM        ((3)),ISCRDA(  4+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DACMU(ISCRDA(  3+IDAA),ONE*FOUR       ,ISCRDA(  6+IDAA))     
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  4+IDAA),ISCRDA(  7+IDAA))    
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),DUM        ((5)))                     
        if(ibbc.eq.1.and.(abs(dare(dum(4))).gt.pieni.and.               &
     &abs(dare(dum(5))).gt.pieni)) then
          ibbc1=1
*FOX    DUM(5)=SQRT(DUM(5)) ;                                           *FOX
      CALL DACOP(DUM        ((5)),ISCRDA(  1+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))              
      CALL DACOP(ISCRDA(  2+IDAA),DUM        ((5)))                     
        else
          ibbc1=0
        endif
*FOX    SEPX0=TRACK(1)+TRACK(2)*S-STAR(1,JSLI) ;                        *FOX
      CALL DACOP(TRACK      ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(TRACK      ((2)),ISCRDA(  2+IDAA))                     
      RSCRRI(  3+IDAA) = STAR       ((1),JSLI       )                   
      CALL DAMUL(ISCRDA(  2+IDAA),S          ,ISCRDA(  4+IDAA))         
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACSU(ISCRDA(  5+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  6+IDAA),SEPX0      )                          
*FOX    SEPY0=TRACK(3)+TRACK(4)*S-STAR(2,JSLI) ;                        *FOX
      CALL DACOP(TRACK      ((3)),ISCRDA(  1+IDAA))                     
      CALL DACOP(TRACK      ((4)),ISCRDA(  2+IDAA))                     
      RSCRRI(  3+IDAA) = STAR       ((2),JSLI       )                   
      CALL DAMUL(ISCRDA(  2+IDAA),S          ,ISCRDA(  4+IDAA))         
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACSU(ISCRDA(  5+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  6+IDAA),SEPY0      )                          
        if(ibbc1.eq.1) then
          sfac=one
!hr05     if(dare(dum(4)).lt.zero) sfac=-one
          if(dare(dum(4)).lt.zero) sfac=(-1d0*one)                       !hr05
*FOX    DUM(6)=SFAC*DUM(4)/DUM(5) ;                                     *FOX
      CALL DACOP(DUM        ((4)),ISCRDA(  1+IDAA))                     
      CALL DACOP(DUM        ((5)),ISCRDA(  2+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*SFAC       ,ISCRDA(  3+IDAA))     
      CALL DADIV(ISCRDA(  3+IDAA),ISCRDA(  2+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),DUM        ((6)))                     
*FOX    DUM(7)=DUM(1)+DUM(2) ;                                          *FOX
      CALL DACOP(DUM        ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(DUM        ((2)),ISCRDA(  2+IDAA))                     
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),DUM        ((7)))                     
*FOX    COSTH=HALF*(ONE+DUM(6)) ;                                       *FOX
      CALL DACOP(DUM        ((6)),ISCRDA(  1+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  2+IDAA))     
      CALL DACMU(ISCRDA(  2+IDAA),ONE*HALF       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),COSTH      )                          
          if(abs(dare(costh)).gt.pieni) then
*FOX    COSTH=SQRT(COSTH) ;                                             *FOX
      CALL DAFUN('SQRT',COSTH      ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),COSTH      )                          
          else
*FOX    COSTH=ZERO ;                                                    *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(COSTH      ,RSCRRI(100))                               
          endif
*FOX    SINTH=HALF*(ONE-DUM(6)) ;                                       *FOX
      CALL DACOP(DUM        ((6)),ISCRDA(  1+IDAA))                     
      CALL DASUC(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  2+IDAA))     
      CALL DACMU(ISCRDA(  2+IDAA),ONE*HALF       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),SINTH      )                          
          if(abs(dare(sinth)).gt.pieni) then
*FOX    SINTH=-SFAC*SQRT(SINTH) ;                                       *FOX
      CALL DAFUN('SQRT',SINTH      ,ISCRDA(  1+IDAA))                   
      RSCRRI(  2+IDAA) = (-ONE       ) * SFAC                           
      CALL DACMU(ISCRDA(  1+IDAA),ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  3+IDAA),SINTH      )                          
          else
*FOX    SINTH=ZERO ;                                                    *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(SINTH      ,RSCRRI(100))                               
          endif
          if(dare(dum(3)).lt.zero) then
*FOX    SINTH=-SINTH ;                                                  *FOX
      CALL DACMU(SINTH      ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),SINTH      )                          
          endif
*FOX    SY=SFAC*DUM(5) ;                                                *FOX
      CALL DACOP(DUM        ((5)),ISCRDA(  1+IDAA))                     
      CALL DACMU(ISCRDA(  1+IDAA),ONE*SFAC       ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),SY         )                          
*FOX    SX=(DUM(7)+SY)*HALF ;                                           *FOX
      CALL DACOP(DUM        ((7)),ISCRDA(  1+IDAA))                     
      CALL DAADD(ISCRDA(  1+IDAA),SY         ,ISCRDA(  2+IDAA))         
      CALL DACMU(ISCRDA(  2+IDAA),ONE*HALF       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),SX         )                          
*FOX    SY=(DUM(7)-SY)*HALF ;                                           *FOX
      CALL DACOP(DUM        ((7)),ISCRDA(  1+IDAA))                     
      CALL DASUB(ISCRDA(  1+IDAA),SY         ,ISCRDA(  2+IDAA))         
      CALL DACMU(ISCRDA(  2+IDAA),ONE*HALF       ,ISCRDA(  3+IDAA))     
      CALL DACOP(ISCRDA(  3+IDAA),SY         )                          
*FOX    SEPX=SEPX0*COSTH+SEPY0*SINTH ;                                  *FOX
      CALL DAMUL(SEPX0      ,COSTH      ,ISCRDA(  1+IDAA))              
      CALL DAMUL(SEPY0      ,SINTH      ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),SEPX       )                          
*FOX    SEPY=-SEPX0*SINTH+SEPY0*COSTH ;                                 *FOX
      CALL DACMU(SEPX0      ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),SINTH      ,ISCRDA(  2+IDAA))         
      CALL DAMUL(SEPY0      ,COSTH      ,ISCRDA(  3+IDAA))              
      CALL DAADD(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),SEPY       )                          
        else
*FOX    SX=DUM(1) ;                                                     *FOX
      CALL DACOP(DUM        ((1)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),SX         )                          
*FOX    SY=DUM(2) ;                                                     *FOX
      CALL DACOP(DUM        ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(ISCRDA(  1+IDAA),SY         )                          
*FOX    SEPX=SEPX0 ;                                                    *FOX
      CALL DACOP(SEPX0      ,SEPX       )                               
*FOX    SEPY=SEPY0 ;                                                    *FOX
      CALL DACOP(SEPY0      ,SEPY       )                               
        endif
        if(dare(sx).gt.dare(sy)) then
          call bbff(sepx,sepy,sx,sy,bbfx,bbfy,bbgx,bbgy)
        else
          call bbff(sepy,sepx,sy,sx,bbfy,bbfx,bbgy,bbgx)
        endif
*FOX    BBFX=F*BBFX ;                                                   *FOX
      CALL DACMU(BBFX       ,ONE*F          ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),BBFX       )                          
*FOX    BBFY=F*BBFY ;                                                   *FOX
      CALL DACMU(BBFY       ,ONE*F          ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),BBFY       )                          
*FOX    BBGX=F*BBGX ;                                                   *FOX
      CALL DACMU(BBGX       ,ONE*F          ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),BBGX       )                          
*FOX    BBGY=F*BBGY ;                                                   *FOX
      CALL DACMU(BBGY       ,ONE*F          ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),BBGY       )                          
        if(ibbc1.eq.1) then
*FOX    DUM(8)=TWO*(BCU(IBB,4)-BCU(IBB,9)+(BCU(IBB,6)-BCU(IBB,10))*SP) ;*FOX
      RSCRRI(  1+IDAA) = BCU        (IBB        ,(6))                   
      RSCRRI(  2+IDAA) = BCU        (IBB        ,(10))                  
      RSCRRI(  3+IDAA) = RSCRRI(  1+IDAA) - RSCRRI(  2+IDAA)            
      RSCRRI(  4+IDAA) = BCU        (IBB        ,(4))                   
      RSCRRI(  5+IDAA) = BCU        (IBB        ,(9))                   
      CALL DACMU(SP         ,ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))     
      RSCRRI(  7+IDAA) = RSCRRI(  4+IDAA) - RSCRRI(  5+IDAA)            
      CALL DACAD(ISCRDA(  6+IDAA),ONE*RSCRRI(  7+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DACMU(ISCRDA(  8+IDAA),ONE*TWO        ,ISCRDA(  9+IDAA))     
      CALL DACOP(ISCRDA(  9+IDAA),DUM        ((8)))                     
*FOX    DUM(9)=BCU(IBB,5)+BCU(IBB,7)+TWO*BCU(IBB,8)*SP ;                *FOX
      RSCRRI(  1+IDAA) = BCU        (IBB        ,(5))                   
      RSCRRI(  2+IDAA) = BCU        (IBB        ,(7))                   
      RSCRRI(  3+IDAA) = BCU        (IBB        ,(8))                   
      RSCRRI(  4+IDAA) = TWO         * RSCRRI(  3+IDAA)                 
      CALL DACMU(SP         ,ONE*RSCRRI(  4+IDAA),ISCRDA(  5+IDAA))     
      RSCRRI(  6+IDAA) = RSCRRI(  1+IDAA) + RSCRRI(  2+IDAA)            
      CALL DACAD(ISCRDA(  5+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  7+IDAA),DUM        ((9)))                     
*FOX    DUM(10)=(DUM(4)*DUM(8)+FOUR*DUM(3)*DUM(9))/                     *FOX
*FOX    DUM(5)/DUM(5)/DUM(5) ;                                          *FOX
      CALL DACOP(DUM        ((4)),ISCRDA(  1+IDAA))                     
      CALL DACOP(DUM        ((8)),ISCRDA(  2+IDAA))                     
      CALL DACOP(DUM        ((3)),ISCRDA(  3+IDAA))                     
      CALL DACOP(DUM        ((9)),ISCRDA(  4+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DACMU(ISCRDA(  3+IDAA),ONE*FOUR       ,ISCRDA(  6+IDAA))     
      CALL DAMUL(ISCRDA(  6+IDAA),ISCRDA(  4+IDAA),ISCRDA(  7+IDAA))    
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  7+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(DUM        ((5)),ISCRDA(  9+IDAA))                     
      CALL DACOP(DUM        ((5)),ISCRDA( 10+IDAA))                     
      CALL DACOP(DUM        ((5)),ISCRDA( 11+IDAA))                     
      CALL DADIV(ISCRDA(  8+IDAA),ISCRDA(  9+IDAA),ISCRDA( 12+IDAA))    
      CALL DADIV(ISCRDA( 12+IDAA),ISCRDA( 10+IDAA),ISCRDA( 13+IDAA))    
      CALL DADIV(ISCRDA( 13+IDAA),ISCRDA( 11+IDAA),ISCRDA( 14+IDAA))    
      CALL DACOP(ISCRDA( 14+IDAA),DUM        ((10)))                    
*FOX    DUM(11)=SFAC*(DUM(8)/DUM(5)-DUM(4)*DUM(10)) ;                   *FOX
      CALL DACOP(DUM        ((8)),ISCRDA(  1+IDAA))                     
      CALL DACOP(DUM        ((5)),ISCRDA(  2+IDAA))                     
      CALL DACOP(DUM        ((4)),ISCRDA(  3+IDAA))                     
      CALL DACOP(DUM        ((10)),ISCRDA(  4+IDAA))                    
      CALL DADIV(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  6+IDAA))    
      CALL DASUB(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA(  7+IDAA))    
      CALL DACMU(ISCRDA(  7+IDAA),ONE*SFAC       ,ISCRDA(  8+IDAA))     
      CALL DACOP(ISCRDA(  8+IDAA),DUM        ((11)))                    
*FOX    DUM(12)=BCU(IBB,4)+BCU(IBB,9)+(BCU(IBB,6)+BCU(IBB,10))*SP ;     *FOX
      RSCRRI(  1+IDAA) = BCU        (IBB        ,(6))                   
      RSCRRI(  2+IDAA) = BCU        (IBB        ,(10))                  
      RSCRRI(  3+IDAA) = RSCRRI(  1+IDAA) + RSCRRI(  2+IDAA)            
      RSCRRI(  4+IDAA) = BCU        (IBB        ,(4))                   
      RSCRRI(  5+IDAA) = BCU        (IBB        ,(9))                   
      CALL DACMU(SP         ,ONE*RSCRRI(  3+IDAA),ISCRDA(  6+IDAA))     
      RSCRRI(  7+IDAA) = RSCRRI(  4+IDAA) + RSCRRI(  5+IDAA)            
      CALL DACAD(ISCRDA(  6+IDAA),ONE*RSCRRI(  7+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  8+IDAA),DUM        ((12)))                    
*FOX    DUM(13)=SFAC*(DUM(4)*DUM(8)*HALF+TWO*DUM(3)*DUM(9))/DUM(5) ;    *FOX
      CALL DACOP(DUM        ((4)),ISCRDA(  1+IDAA))                     
      CALL DACOP(DUM        ((8)),ISCRDA(  2+IDAA))                     
      CALL DACOP(DUM        ((3)),ISCRDA(  3+IDAA))                     
      CALL DACOP(DUM        ((9)),ISCRDA(  4+IDAA))                     
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  5+IDAA))    
      CALL DACMU(ISCRDA(  5+IDAA),ONE*HALF       ,ISCRDA(  6+IDAA))     
      CALL DACMU(ISCRDA(  3+IDAA),ONE*TWO        ,ISCRDA(  7+IDAA))     
      CALL DAMUL(ISCRDA(  7+IDAA),ISCRDA(  4+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  6+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(DUM        ((5)),ISCRDA( 10+IDAA))                     
      CALL DACMU(ISCRDA(  9+IDAA),ONE*SFAC       ,ISCRDA( 11+IDAA))     
      CALL DADIV(ISCRDA( 11+IDAA),ISCRDA( 10+IDAA),ISCRDA( 12+IDAA))    
      CALL DACOP(ISCRDA( 12+IDAA),DUM        ((13)))                    
          if(abs(dare(costh)).gt.pieni) then
*FOX    COSTHP=DUM(11)/FOUR/COSTH ;                                     *FOX
      CALL DACOP(DUM        ((11)),ISCRDA(  1+IDAA))                    
      CALL DACDI(ISCRDA(  1+IDAA),ONE*FOUR       ,ISCRDA(  2+IDAA))     
      CALL DADIV(ISCRDA(  2+IDAA),COSTH      ,ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),COSTHP     )                          
          else
*FOX    COSTHP=ZERO ;                                                   *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(COSTHP     ,RSCRRI(100))                               
          endif
          if(abs(dare(sinth)).gt.pieni) then
*FOX    SINTHP=-DUM(11)/FOUR/SINTH ;                                    *FOX
      CALL DACOP(DUM        ((11)),ISCRDA(  1+IDAA))                    
      CALL DACMU(ISCRDA(  1+IDAA),ONE*(-ONE       ),ISCRDA(  2+IDAA))   
      CALL DACDI(ISCRDA(  2+IDAA),ONE*FOUR       ,ISCRDA(  3+IDAA))     
      CALL DADIV(ISCRDA(  3+IDAA),SINTH      ,ISCRDA(  4+IDAA))         
      CALL DACOP(ISCRDA(  4+IDAA),SINTHP     )                          
          else
*FOX    SINTHP=ZERO ;                                                   *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(SINTHP     ,RSCRRI(100))                               
          endif
*FOX    TRACK(6)=TRACK(6)-                                              *FOX
*FOX    (BBFX*(COSTHP*SEPX0+SINTHP*SEPY0)+                              *FOX
*FOX    BBFY*(-SINTHP*SEPX0+COSTHP*SEPY0)+                              *FOX
*FOX    BBGX*(DUM(12)+DUM(13))+BBGY*(DUM(12)-DUM(13)))/                 *FOX
*FOX    CPHI*HALF ;
      CALL DACOP(DUM        ((12)),ISCRDA(  1+IDAA))                    
      CALL DACOP(DUM        ((13)),ISCRDA(  2+IDAA))                    
      CALL DACOP(DUM        ((12)),ISCRDA(  3+IDAA))                    
      CALL DACOP(DUM        ((13)),ISCRDA(  4+IDAA))                    
      CALL DAMUL(COSTHP     ,SEPX0      ,ISCRDA(  5+IDAA))              
      CALL DAMUL(SINTHP     ,SEPY0      ,ISCRDA(  6+IDAA))              
      CALL DACMU(SINTHP     ,ONE*(-ONE       ),ISCRDA(  7+IDAA))        
      CALL DAMUL(ISCRDA(  7+IDAA),SEPX0      ,ISCRDA(  8+IDAA))         
      CALL DAMUL(COSTHP     ,SEPY0      ,ISCRDA(  9+IDAA))              
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA( 10+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  9+IDAA),ISCRDA( 11+IDAA))    
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA( 12+IDAA))    
      CALL DASUB(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA( 13+IDAA))    
      CALL DAMUL(BBFX       ,ISCRDA( 10+IDAA),ISCRDA( 14+IDAA))         
      CALL DAMUL(BBFY       ,ISCRDA( 11+IDAA),ISCRDA( 15+IDAA))         
      CALL DAMUL(BBGX       ,ISCRDA( 12+IDAA),ISCRDA( 16+IDAA))         
      CALL DAMUL(BBGY       ,ISCRDA( 13+IDAA),ISCRDA( 17+IDAA))         
      CALL DAADD(ISCRDA( 14+IDAA),ISCRDA( 15+IDAA),ISCRDA( 18+IDAA))    
      CALL DAADD(ISCRDA( 18+IDAA),ISCRDA( 16+IDAA),ISCRDA( 19+IDAA))    
      CALL DAADD(ISCRDA( 19+IDAA),ISCRDA( 17+IDAA),ISCRDA( 20+IDAA))    
      CALL DACOP(TRACK      ((6)),ISCRDA( 21+IDAA))                     
      CALL DACDI(ISCRDA( 20+IDAA),ONE*CPHI       ,ISCRDA( 22+IDAA))     
      CALL DACMU(ISCRDA( 22+IDAA),ONE*HALF       ,ISCRDA( 23+IDAA))     
      CALL DASUB(ISCRDA( 21+IDAA),ISCRDA( 23+IDAA),ISCRDA( 24+IDAA))    
      CALL DACOP(ISCRDA( 24+IDAA),TRACK      ((6)))                     
*FOX    BBF0=BBFX ;                                                     *FOX
      CALL DACOP(BBFX       ,BBF0       )                               
*FOX    BBFX=BBF0*COSTH-BBFY*SINTH ;                                    *FOX
      CALL DAMUL(BBF0       ,COSTH      ,ISCRDA(  1+IDAA))              
      CALL DAMUL(BBFY       ,SINTH      ,ISCRDA(  2+IDAA))              
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),BBFX       )                          
*FOX    BBFY=BBF0*SINTH+BBFY*COSTH ;                                    *FOX
      CALL DAMUL(BBF0       ,SINTH      ,ISCRDA(  1+IDAA))              
      CALL DAMUL(BBFY       ,COSTH      ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),BBFY       )                          
        else
*FOX    TRACK(6)=TRACK(6)-                                              *FOX
*FOX    (BBGX*(BCU(IBB,4)+BCU(IBB,6)*SP)+                               *FOX
*FOX    BBGY*(BCU(IBB,9)+BCU(IBB,10)*SP))/CPHI ;                        *FOX
      RSCRRI(  1+IDAA) = BCU        (IBB        ,(4))                   
      RSCRRI(  2+IDAA) = BCU        (IBB        ,(6))                   
      RSCRRI(  3+IDAA) = BCU        (IBB        ,(9))                   
      RSCRRI(  4+IDAA) = BCU        (IBB        ,(10))                  
      CALL DACMU(SP         ,ONE*RSCRRI(  2+IDAA),ISCRDA(  5+IDAA))     
      CALL DACMU(SP         ,ONE*RSCRRI(  4+IDAA),ISCRDA(  6+IDAA))     
      CALL DACAD(ISCRDA(  5+IDAA),ONE*RSCRRI(  1+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DACAD(ISCRDA(  6+IDAA),ONE*RSCRRI(  3+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DAMUL(BBGX       ,ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))         
      CALL DAMUL(BBGY       ,ISCRDA(  8+IDAA),ISCRDA( 10+IDAA))         
      CALL DAADD(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 11+IDAA))    
      CALL DACOP(TRACK      ((6)),ISCRDA( 12+IDAA))                     
      CALL DACDI(ISCRDA( 11+IDAA),ONE*CPHI       ,ISCRDA( 13+IDAA))     
      CALL DASUB(ISCRDA( 12+IDAA),ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))    
      CALL DACOP(ISCRDA( 14+IDAA),TRACK      ((6)))                     
        endif
*FOX    TRACK(6)=TRACK(6)-(BBFX*(TRACK(2)-BBFX*HALF)+                   *FOX
*FOX    BBFY*(TRACK(4)-BBFY*HALF))*HALF ;                               *FOX
      CALL DACOP(TRACK      ((2)),ISCRDA(  1+IDAA))                     
      CALL DACOP(TRACK      ((4)),ISCRDA(  2+IDAA))                     
      CALL DACMU(BBFX       ,ONE*HALF       ,ISCRDA(  3+IDAA))          
      CALL DACMU(BBFY       ,ONE*HALF       ,ISCRDA(  4+IDAA))          
      CALL DASUB(ISCRDA(  1+IDAA),ISCRDA(  3+IDAA),ISCRDA(  5+IDAA))    
      CALL DASUB(ISCRDA(  2+IDAA),ISCRDA(  4+IDAA),ISCRDA(  6+IDAA))    
      CALL DAMUL(BBFX       ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DAMUL(BBFY       ,ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))         
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(TRACK      ((6)),ISCRDA( 10+IDAA))                     
      CALL DACMU(ISCRDA(  9+IDAA),ONE*HALF       ,ISCRDA( 11+IDAA))     
      CALL DASUB(ISCRDA( 10+IDAA),ISCRDA( 11+IDAA),ISCRDA( 12+IDAA))    
      CALL DACOP(ISCRDA( 12+IDAA),TRACK      ((6)))                     
*FOX    TRACK(1)=TRACK(1)+S*BBFX ;                                      *FOX
      CALL DACOP(TRACK      ((1)),ISCRDA(  1+IDAA))                     
      CALL DAMUL(S          ,BBFX       ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),TRACK      ((1)))                     
*FOX    TRACK(2)=TRACK(2)-BBFX ;                                        *FOX
      CALL DACOP(TRACK      ((2)),ISCRDA(  1+IDAA))                     
      CALL DASUB(ISCRDA(  1+IDAA),BBFX       ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),TRACK      ((2)))                     
*FOX    TRACK(3)=TRACK(3)+S*BBFY ;                                      *FOX
      CALL DACOP(TRACK      ((3)),ISCRDA(  1+IDAA))                     
      CALL DAMUL(S          ,BBFY       ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),TRACK      ((3)))                     
*FOX    TRACK(4)=TRACK(4)-BBFY ;                                        *FOX
      CALL DACOP(TRACK      ((4)),ISCRDA(  1+IDAA))                     
      CALL DASUB(ISCRDA(  1+IDAA),BBFY       ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),TRACK      ((4)))                     
 2000 continue
        CALL DADAL(BBGY    ,1)                                                  
        CALL DADAL(BBGX    ,1)                                                  
        CALL DADAL(BBF0    ,1)                                                  
        CALL DADAL(BBFY    ,1)                                                  
        CALL DADAL(BBFX    ,1)                                                  
        CALL DADAL(SINTHP  ,1)                                                  
        CALL DADAL(COSTHP  ,1)                                                  
        CALL DADAL(SINTH   ,1)                                                  
        CALL DADAL(COSTH   ,1)                                                  
        CALL DADAL(SEPY0   ,1)                                                  
        CALL DADAL(SEPX0   ,1)                                                  
        CALL DADAL(SEPY    ,1)                                                  
        CALL DADAL(SEPX    ,1)                                                  
        CALL DADAL(SY      ,1)                                                  
        CALL DADAL(SX      ,1)                                                  
        CALL DADAL(DUM     ,1*(13))                                             
        CALL DADAL(DD      ,1)                                                  
        CALL DADAL(SP      ,1)                                                  
        CALL DADAL(S       ,1)                                                  
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
*FOX{
      INTEGER TRACK   (6)
      INTEGER A1      
      INTEGER H1      
      INTEGER SQR1A   
      INTEGER H1D     
      INTEGER H1X     
      INTEGER H1Y     
      INTEGER H1Z     
      INTEGER DET     
      INTEGER X1      
      INTEGER Y1      
      INTEGER Z1      
      INTEGER ISCRDA, ISCRRI,IDAO
      DOUBLE PRECISION RSCRRI
      COMMON/DASCR/ISCRDA(100),RSCRRI(100),ISCRRI(100),IDAO
      if(1.eq.1) then                                                 
         CALL DAKEY('FOX V2.1')
         CALL DAALL(A1      ,1,'A1        ',NORD,NVAR)
         CALL DAALL(H1      ,1,'H1        ',NORD,NVAR)
         CALL DAALL(SQR1A   ,1,'SQR1A     ',NORD,NVAR)
         CALL DAALL(H1D     ,1,'H1D       ',NORD,NVAR)
         CALL DAALL(H1X     ,1,'H1X       ',NORD,NVAR)
         CALL DAALL(H1Y     ,1,'H1Y       ',NORD,NVAR)
         CALL DAALL(H1Z     ,1,'H1Z       ',NORD,NVAR)
         CALL DAALL(DET     ,1,'DET       ',NORD,NVAR)
         CALL DAALL(X1      ,1,'X1        ',NORD,NVAR)
         CALL DAALL(Y1      ,1,'Y1        ',NORD,NVAR)
         CALL DAALL(Z1      ,1,'Z1        ',NORD,NVAR)
      ENDIF
      IDAA = IDAO
*FOX}
!-----------------------------------------------------------------------
*FOX    H1D=SQRT((ONE+TRACK(6))*(ONE+TRACK(6))-                         *FOX
*FOX    TRACK(2)*TRACK(2)-TRACK(4)*TRACK(4)) ;                          *FOX
      CALL DACOP(TRACK      ((6)),ISCRDA(  1+IDAA))                     
      CALL DACOP(TRACK      ((6)),ISCRDA(  2+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  3+IDAA))     
      CALL DACAD(ISCRDA(  2+IDAA),ONE*ONE        ,ISCRDA(  4+IDAA))     
      CALL DACOP(TRACK      ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(TRACK      ((2)),ISCRDA(  6+IDAA))                     
      CALL DACOP(TRACK      ((4)),ISCRDA(  7+IDAA))                     
      CALL DACOP(TRACK      ((4)),ISCRDA(  8+IDAA))                     
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  9+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA( 10+IDAA))    
      CALL DAMUL(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 11+IDAA))    
      CALL DASUB(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 12+IDAA))    
      CALL DASUB(ISCRDA( 12+IDAA),ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))    
      CALL DAFUN('SQRT',ISCRDA( 13+IDAA),ISCRDA( 14+IDAA))              
      CALL DACOP(ISCRDA( 14+IDAA),H1D        )                          
*FOX    H1X=TRACK(2)/H1D ;                                              *FOX
      CALL DACOP(TRACK      ((2)),ISCRDA(  1+IDAA))                     
      CALL DADIV(ISCRDA(  1+IDAA),H1D        ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),H1X        )                          
*FOX    H1Y=TRACK(4)/H1D ;                                              *FOX
      CALL DACOP(TRACK      ((4)),ISCRDA(  1+IDAA))                     
      CALL DADIV(ISCRDA(  1+IDAA),H1D        ,ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),H1Y        )                          
*FOX    H1Z=ONE-(ONE+TRACK(6))/H1D ;                                    *FOX
      CALL DACOP(TRACK      ((6)),ISCRDA(  1+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  2+IDAA))     
      CALL DADIV(ISCRDA(  2+IDAA),H1D        ,ISCRDA(  3+IDAA))         
      CALL DASUC(ISCRDA(  3+IDAA),ONE*ONE        ,ISCRDA(  4+IDAA))     
      CALL DACOP(ISCRDA(  4+IDAA),H1Z        )                          
*FOX    H1=(TRACK(6)+ONE-SQRT((ONE+TRACK(6))*(ONE+TRACK(6))-            *FOX
*FOX    TRACK(2)*TRACK(2)-TRACK(4)*TRACK(4)))*CPHI*CPHI ;               *FOX
      CALL DACOP(TRACK      ((6)),ISCRDA(  1+IDAA))                     
      CALL DACOP(TRACK      ((6)),ISCRDA(  2+IDAA))                     
      CALL DACAD(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  3+IDAA))     
      CALL DACAD(ISCRDA(  2+IDAA),ONE*ONE        ,ISCRDA(  4+IDAA))     
      CALL DACOP(TRACK      ((2)),ISCRDA(  5+IDAA))                     
      CALL DACOP(TRACK      ((2)),ISCRDA(  6+IDAA))                     
      CALL DACOP(TRACK      ((4)),ISCRDA(  7+IDAA))                     
      CALL DACOP(TRACK      ((4)),ISCRDA(  8+IDAA))                     
      CALL DAMUL(ISCRDA(  3+IDAA),ISCRDA(  4+IDAA),ISCRDA(  9+IDAA))    
      CALL DAMUL(ISCRDA(  5+IDAA),ISCRDA(  6+IDAA),ISCRDA( 10+IDAA))    
      CALL DAMUL(ISCRDA(  7+IDAA),ISCRDA(  8+IDAA),ISCRDA( 11+IDAA))    
      CALL DASUB(ISCRDA(  9+IDAA),ISCRDA( 10+IDAA),ISCRDA( 12+IDAA))    
      CALL DASUB(ISCRDA( 12+IDAA),ISCRDA( 11+IDAA),ISCRDA( 13+IDAA))    
      CALL DACOP(TRACK      ((6)),ISCRDA( 14+IDAA))                     
      CALL DAFUN('SQRT',ISCRDA( 13+IDAA),ISCRDA( 15+IDAA))              
      CALL DACAD(ISCRDA( 14+IDAA),ONE*ONE        ,ISCRDA( 16+IDAA))     
      CALL DASUB(ISCRDA( 16+IDAA),ISCRDA( 15+IDAA),ISCRDA( 17+IDAA))    
      CALL DACMU(ISCRDA( 17+IDAA),ONE*CPHI       ,ISCRDA( 18+IDAA))     
      CALL DACMU(ISCRDA( 18+IDAA),ONE*CPHI       ,ISCRDA( 19+IDAA))     
      CALL DACOP(ISCRDA( 19+IDAA),H1         )                          
*FOX    DET=ONE/CPHI+TPHI*(H1X*CALPHA+H1Y*SALPHA-H1Z*SPHI) ;            *FOX
      CALL DACMU(H1X        ,ONE*CALPHA     ,ISCRDA(  1+IDAA))          
      CALL DACMU(H1Y        ,ONE*SALPHA     ,ISCRDA(  2+IDAA))          
      CALL DACMU(H1Z        ,ONE*SPHI       ,ISCRDA(  3+IDAA))          
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  4+IDAA))    
      CALL DASUB(ISCRDA(  4+IDAA),ISCRDA(  3+IDAA),ISCRDA(  5+IDAA))    
      RSCRRI(  6+IDAA) = ONE         / CPHI                             
      CALL DACMU(ISCRDA(  5+IDAA),ONE*TPHI       ,ISCRDA(  7+IDAA))     
      CALL DACAD(ISCRDA(  7+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  8+IDAA))
     *                                                                  
      CALL DACOP(ISCRDA(  8+IDAA),DET        )                          
*FOX    X1= TRACK(1)*(ONE/CPHI+SALPHA*(H1Y-H1Z*SALPHA*SPHI)*TPHI)       *FOX
*FOX       +TRACK(3)*SALPHA*TPHI*(-H1X+H1Z*CALPHA*SPHI)                 *FOX
*FOX       -TRACK(5)*(CALPHA+H1Y*CALPHA*SALPHA*SPHI                     *FOX
*FOX       -H1X*SALPHA*SALPHA*SPHI)*TPHI ;                              *FOX
      CALL DACMU(H1Z        ,ONE*SALPHA     ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*SPHI       ,ISCRDA(  2+IDAA))     
      CALL DASUB(H1Y        ,ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))         
      RSCRRI(  4+IDAA) = ONE         / CPHI                             
      CALL DACMU(ISCRDA(  3+IDAA),ONE*SALPHA     ,ISCRDA(  5+IDAA))     
      CALL DACMU(ISCRDA(  5+IDAA),ONE*TPHI       ,ISCRDA(  6+IDAA))     
      CALL DACMU(H1X        ,ONE*(-ONE       ),ISCRDA(  7+IDAA))        
      CALL DACMU(H1Z        ,ONE*CALPHA     ,ISCRDA(  8+IDAA))          
      CALL DACMU(ISCRDA(  8+IDAA),ONE*SPHI       ,ISCRDA(  9+IDAA))     
      CALL DACMU(H1Y        ,ONE*CALPHA     ,ISCRDA( 10+IDAA))          
      CALL DACMU(ISCRDA( 10+IDAA),ONE*SALPHA     ,ISCRDA( 11+IDAA))     
      CALL DACMU(ISCRDA( 11+IDAA),ONE*SPHI       ,ISCRDA( 12+IDAA))     
      CALL DACMU(H1X        ,ONE*SALPHA     ,ISCRDA( 13+IDAA))          
      CALL DACMU(ISCRDA( 13+IDAA),ONE*SALPHA     ,ISCRDA( 14+IDAA))     
      CALL DACMU(ISCRDA( 14+IDAA),ONE*SPHI       ,ISCRDA( 15+IDAA))     
      CALL DACAD(ISCRDA(  6+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA( 16+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  7+IDAA),ISCRDA(  9+IDAA),ISCRDA( 17+IDAA))    
      CALL DACAD(ISCRDA( 12+IDAA),ONE*CALPHA     ,ISCRDA( 18+IDAA))     
      CALL DASUB(ISCRDA( 18+IDAA),ISCRDA( 15+IDAA),ISCRDA( 19+IDAA))    
      CALL DACOP(TRACK      ((1)),ISCRDA( 20+IDAA))                     
      CALL DACOP(TRACK      ((3)),ISCRDA( 21+IDAA))                     
      CALL DACOP(TRACK      ((5)),ISCRDA( 22+IDAA))                     
      CALL DAMUL(ISCRDA( 20+IDAA),ISCRDA( 16+IDAA),ISCRDA( 23+IDAA))    
      CALL DACMU(ISCRDA( 21+IDAA),ONE*SALPHA     ,ISCRDA( 24+IDAA))     
      CALL DACMU(ISCRDA( 24+IDAA),ONE*TPHI       ,ISCRDA( 25+IDAA))     
      CALL DAMUL(ISCRDA( 25+IDAA),ISCRDA( 17+IDAA),ISCRDA( 26+IDAA))    
      CALL DAMUL(ISCRDA( 22+IDAA),ISCRDA( 19+IDAA),ISCRDA( 27+IDAA))    
      CALL DACMU(ISCRDA( 27+IDAA),ONE*TPHI       ,ISCRDA( 28+IDAA))     
      CALL DAADD(ISCRDA( 23+IDAA),ISCRDA( 26+IDAA),ISCRDA( 29+IDAA))    
      CALL DASUB(ISCRDA( 29+IDAA),ISCRDA( 28+IDAA),ISCRDA( 30+IDAA))    
      CALL DACOP(ISCRDA( 30+IDAA),X1         )                          
*FOX    Y1= TRACK(1)*CALPHA*TPHI*(-H1Y+H1Z*SALPHA*SPHI)                 *FOX
*FOX       +TRACK(3)*(ONE/CPHI+CALPHA*(H1X-H1Z*CALPHA*SPHI)*TPHI)       *FOX
*FOX       -TRACK(5)*(SALPHA-H1Y*CALPHA*CALPHA*SPHI                     *FOX
*FOX       +H1X*CALPHA*SALPHA*SPHI)*TPHI ;                              *FOX
      CALL DACMU(H1Z        ,ONE*CALPHA     ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*SPHI       ,ISCRDA(  2+IDAA))     
      CALL DASUB(H1X        ,ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))         
      CALL DACMU(H1Y        ,ONE*(-ONE       ),ISCRDA(  4+IDAA))        
      CALL DACMU(H1Z        ,ONE*SALPHA     ,ISCRDA(  5+IDAA))          
      CALL DACMU(ISCRDA(  5+IDAA),ONE*SPHI       ,ISCRDA(  6+IDAA))     
      RSCRRI(  7+IDAA) = ONE         / CPHI                             
      CALL DACMU(ISCRDA(  3+IDAA),ONE*CALPHA     ,ISCRDA(  8+IDAA))     
      CALL DACMU(ISCRDA(  8+IDAA),ONE*TPHI       ,ISCRDA(  9+IDAA))     
      CALL DACMU(H1Y        ,ONE*CALPHA     ,ISCRDA( 10+IDAA))          
      CALL DACMU(ISCRDA( 10+IDAA),ONE*CALPHA     ,ISCRDA( 11+IDAA))     
      CALL DACMU(ISCRDA( 11+IDAA),ONE*SPHI       ,ISCRDA( 12+IDAA))     
      CALL DACMU(H1X        ,ONE*CALPHA     ,ISCRDA( 13+IDAA))          
      CALL DACMU(ISCRDA( 13+IDAA),ONE*SALPHA     ,ISCRDA( 14+IDAA))     
      CALL DACMU(ISCRDA( 14+IDAA),ONE*SPHI       ,ISCRDA( 15+IDAA))     
      CALL DAADD(ISCRDA(  4+IDAA),ISCRDA(  6+IDAA),ISCRDA( 16+IDAA))    
      CALL DACAD(ISCRDA(  9+IDAA),ONE*RSCRRI(  7+IDAA),ISCRDA( 17+IDAA))
     *                                                                  
      CALL DASUC(ISCRDA( 12+IDAA),ONE*SALPHA     ,ISCRDA( 18+IDAA))     
      CALL DAADD(ISCRDA( 18+IDAA),ISCRDA( 15+IDAA),ISCRDA( 19+IDAA))    
      CALL DACOP(TRACK      ((1)),ISCRDA( 20+IDAA))                     
      CALL DACOP(TRACK      ((3)),ISCRDA( 21+IDAA))                     
      CALL DACOP(TRACK      ((5)),ISCRDA( 22+IDAA))                     
      CALL DACMU(ISCRDA( 20+IDAA),ONE*CALPHA     ,ISCRDA( 23+IDAA))     
      CALL DACMU(ISCRDA( 23+IDAA),ONE*TPHI       ,ISCRDA( 24+IDAA))     
      CALL DAMUL(ISCRDA( 24+IDAA),ISCRDA( 16+IDAA),ISCRDA( 25+IDAA))    
      CALL DAMUL(ISCRDA( 21+IDAA),ISCRDA( 17+IDAA),ISCRDA( 26+IDAA))    
      CALL DAMUL(ISCRDA( 22+IDAA),ISCRDA( 19+IDAA),ISCRDA( 27+IDAA))    
      CALL DACMU(ISCRDA( 27+IDAA),ONE*TPHI       ,ISCRDA( 28+IDAA))     
      CALL DAADD(ISCRDA( 25+IDAA),ISCRDA( 26+IDAA),ISCRDA( 29+IDAA))    
      CALL DASUB(ISCRDA( 29+IDAA),ISCRDA( 28+IDAA),ISCRDA( 30+IDAA))    
      CALL DACOP(ISCRDA( 30+IDAA),Y1         )                          
*FOX    Z1=-TRACK(1)*H1Z*CALPHA*SPHI                                    *FOX
*FOX       -TRACK(3)*H1Z*SALPHA*SPHI                                    *FOX
*FOX       +TRACK(5)*(ONE+H1X*CALPHA*SPHI+H1Y*SALPHA*SPHI) ;            *FOX
      CALL DACMU(H1X        ,ONE*CALPHA     ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*SPHI       ,ISCRDA(  2+IDAA))     
      CALL DACMU(H1Y        ,ONE*SALPHA     ,ISCRDA(  3+IDAA))          
      CALL DACMU(ISCRDA(  3+IDAA),ONE*SPHI       ,ISCRDA(  4+IDAA))     
      CALL DACAD(ISCRDA(  2+IDAA),ONE*ONE        ,ISCRDA(  5+IDAA))     
      CALL DAADD(ISCRDA(  5+IDAA),ISCRDA(  4+IDAA),ISCRDA(  6+IDAA))    
      CALL DACOP(TRACK      ((1)),ISCRDA(  7+IDAA))                     
      CALL DACOP(TRACK      ((3)),ISCRDA(  8+IDAA))                     
      CALL DACOP(TRACK      ((5)),ISCRDA(  9+IDAA))                     
      CALL DACMU(ISCRDA(  7+IDAA),ONE*(-ONE       ),ISCRDA( 10+IDAA))   
      CALL DAMUL(ISCRDA( 10+IDAA),H1Z        ,ISCRDA( 11+IDAA))         
      CALL DACMU(ISCRDA( 11+IDAA),ONE*CALPHA     ,ISCRDA( 12+IDAA))     
      CALL DACMU(ISCRDA( 12+IDAA),ONE*SPHI       ,ISCRDA( 13+IDAA))     
      CALL DAMUL(ISCRDA(  8+IDAA),H1Z        ,ISCRDA( 14+IDAA))         
      CALL DACMU(ISCRDA( 14+IDAA),ONE*SALPHA     ,ISCRDA( 15+IDAA))     
      CALL DACMU(ISCRDA( 15+IDAA),ONE*SPHI       ,ISCRDA( 16+IDAA))     
      CALL DAMUL(ISCRDA(  9+IDAA),ISCRDA(  6+IDAA),ISCRDA( 17+IDAA))    
      CALL DASUB(ISCRDA( 13+IDAA),ISCRDA( 16+IDAA),ISCRDA( 18+IDAA))    
      CALL DAADD(ISCRDA( 18+IDAA),ISCRDA( 17+IDAA),ISCRDA( 19+IDAA))    
      CALL DACOP(ISCRDA( 19+IDAA),Z1         )                          
*FOX    TRACK(1)=X1/DET ;                                               *FOX
      CALL DADIV(X1         ,DET        ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),TRACK      ((1)))                     
*FOX    TRACK(3)=Y1/DET ;                                               *FOX
      CALL DADIV(Y1         ,DET        ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),TRACK      ((3)))                     
*FOX    TRACK(5)=Z1/DET ;                                               *FOX
      CALL DADIV(Z1         ,DET        ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),TRACK      ((5)))                     
*FOX    TRACK(6)=TRACK(6)+CALPHA*SPHI*TRACK(2)                          *FOX
*FOX            +SALPHA*SPHI*TRACK(4) ;                                 *FOX
      CALL DACOP(TRACK      ((6)),ISCRDA(  1+IDAA))                     
      CALL DACOP(TRACK      ((2)),ISCRDA(  2+IDAA))                     
      CALL DACOP(TRACK      ((4)),ISCRDA(  3+IDAA))                     
      RSCRRI(  4+IDAA) = CALPHA      * SPHI                             
      CALL DACMU(ISCRDA(  2+IDAA),ONE*RSCRRI(  4+IDAA),ISCRDA(  5+IDAA))
     *                                                                  
      RSCRRI(  6+IDAA) = SALPHA      * SPHI                             
      CALL DACMU(ISCRDA(  3+IDAA),ONE*RSCRRI(  6+IDAA),ISCRDA(  7+IDAA))
     *                                                                  
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  5+IDAA),ISCRDA(  8+IDAA))    
      CALL DAADD(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),TRACK      ((6)))                     
*FOX    TRACK(2)=(TRACK(2)+CALPHA*SPHI*H1)*CPHI ;                       *FOX
      CALL DACOP(TRACK      ((2)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = CALPHA      * SPHI                             
      CALL DACMU(H1         ,ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))     
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACMU(ISCRDA(  4+IDAA),ONE*CPHI       ,ISCRDA(  5+IDAA))     
      CALL DACOP(ISCRDA(  5+IDAA),TRACK      ((2)))                     
*FOX    TRACK(4)=(TRACK(4)+SALPHA*SPHI*H1)*CPHI ;                       *FOX
      CALL DACOP(TRACK      ((4)),ISCRDA(  1+IDAA))                     
      RSCRRI(  2+IDAA) = SALPHA      * SPHI                             
      CALL DACMU(H1         ,ONE*RSCRRI(  2+IDAA),ISCRDA(  3+IDAA))     
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACMU(ISCRDA(  4+IDAA),ONE*CPHI       ,ISCRDA(  5+IDAA))     
      CALL DACOP(ISCRDA(  5+IDAA),TRACK      ((4)))                     
        CALL DADAL(Z1      ,1)                                                  
        CALL DADAL(Y1      ,1)                                                  
        CALL DADAL(X1      ,1)                                                  
        CALL DADAL(DET     ,1)                                                  
        CALL DADAL(H1Z     ,1)                                                  
        CALL DADAL(H1Y     ,1)                                                  
        CALL DADAL(H1X     ,1)                                                  
        CALL DADAL(H1D     ,1)                                                  
        CALL DADAL(SQR1A   ,1)                                                  
        CALL DADAL(H1      ,1)                                                  
        CALL DADAL(A1      ,1)                                                  
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
*FOX{
      INTEGER SEPX    
      INTEGER SEPY    
      INTEGER SIGXX   
      INTEGER SIGYY   
      INTEGER BBFX    
      INTEGER BBFY    
      INTEGER BBGX    
      INTEGER BBGY    
      INTEGER SIGXY   
      INTEGER EXPFAC  
      INTEGER X       
      INTEGER FAC     
      INTEGER FAC2    
      INTEGER CONST   
      INTEGER ARG1X   
      INTEGER ARG1Y   
      INTEGER ARG2X   
      INTEGER ARG2Y   
      INTEGER WX1     
      INTEGER WY1     
      INTEGER WX2     
      INTEGER WY2     
      INTEGER COMFAC  
      INTEGER COMFAC2 
      INTEGER ISCRDA, ISCRRI,IDAO
      DOUBLE PRECISION RSCRRI
      COMMON/DASCR/ISCRDA(100),RSCRRI(100),ISCRRI(100),IDAO
      if(1.eq.1) then                                                 
         CALL DAKEY('FOX V2.1')
         CALL DAALL(SIGXY   ,1,'SIGXY     ',NORD,NVAR)
         CALL DAALL(EXPFAC  ,1,'EXPFAC    ',NORD,NVAR)
         CALL DAALL(X       ,1,'X         ',NORD,NVAR)
         CALL DAALL(FAC     ,1,'FAC       ',NORD,NVAR)
         CALL DAALL(FAC2    ,1,'FAC2      ',NORD,NVAR)
         CALL DAALL(CONST   ,1,'CONST     ',NORD,NVAR)
         CALL DAALL(ARG1X   ,1,'ARG1X     ',NORD,NVAR)
         CALL DAALL(ARG1Y   ,1,'ARG1Y     ',NORD,NVAR)
         CALL DAALL(ARG2X   ,1,'ARG2X     ',NORD,NVAR)
         CALL DAALL(ARG2Y   ,1,'ARG2Y     ',NORD,NVAR)
         CALL DAALL(WX1     ,1,'WX1       ',NORD,NVAR)
         CALL DAALL(WY1     ,1,'WY1       ',NORD,NVAR)
         CALL DAALL(WX2     ,1,'WX2       ',NORD,NVAR)
         CALL DAALL(WY2     ,1,'WY2       ',NORD,NVAR)
         CALL DAALL(COMFAC  ,1,'COMFAC    ',NORD,NVAR)
         CALL DAALL(COMFAC2 ,1,'COMFAC2   ',NORD,NVAR)
      ENDIF
      IDAA = IDAO
*FOX}
!-----------------------------------------------------------------------
      if(dare(sigxx).eq.dare(sigyy)) then
*FOX    X=SEPX*SEPX+SEPY*SEPY ;                                         *FOX
      CALL DAMUL(SEPX       ,SEPX       ,ISCRDA(  1+IDAA))              
      CALL DAMUL(SEPY       ,SEPY       ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),X          )                          
        if(abs(dare(sigxx)+dare(sigyy)).gt.pieni) then
*FOX      CONST=X/(SIGXX+SIGYY) ;                                       *FOX
      CALL DAADD(SIGXX      ,SIGYY      ,ISCRDA(  1+IDAA))              
      CALL DADIV(X          ,ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),CONST      )                          
        else
*FOX      CONST=ZERO ;                                                  *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(CONST      ,RSCRRI(100))                               
        endif
*FOX    EXPFAC=EXP(-CONST) ;                                            *FOX
      CALL DACMU(CONST      ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAFUN('EXP ',ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))              
      CALL DACOP(ISCRDA(  2+IDAA),EXPFAC     )                          
        if(abs(dare(x)).gt.pieni) then
*FOX      BBFX=TWO*SEPX*(ONE-EXPFAC)/X ;                                *FOX
      CALL DASUC(EXPFAC     ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(SEPX       ,ONE*TWO        ,ISCRDA(  2+IDAA))          
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DADIV(ISCRDA(  3+IDAA),X          ,ISCRDA(  4+IDAA))         
      CALL DACOP(ISCRDA(  4+IDAA),BBFX       )                          
*FOX      BBFY=TWO*SEPY*(ONE-EXPFAC)/X ;                                *FOX
      CALL DASUC(EXPFAC     ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(SEPY       ,ONE*TWO        ,ISCRDA(  2+IDAA))          
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DADIV(ISCRDA(  3+IDAA),X          ,ISCRDA(  4+IDAA))         
      CALL DACOP(ISCRDA(  4+IDAA),BBFY       )                          
*FOX      COMFAC=-SEPX*BBFX+SEPY*BBFY ;                                 *FOX
      CALL DACMU(SEPX       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),BBFX       ,ISCRDA(  2+IDAA))         
      CALL DAMUL(SEPY       ,BBFY       ,ISCRDA(  3+IDAA))              
      CALL DAADD(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),COMFAC     )                          
          if(dare(sigxx).lt.zero) then
*FOX        SIGXX=-SIGXX ;                                              *FOX
      CALL DACMU(SIGXX      ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),SIGXX      )                          
          endif
          if(dare(sigyy).lt.zero) then
*FOX        SIGYY=-SIGYY ;                                              *FOX
      CALL DACMU(SIGYY      ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),SIGYY      )                          
          endif
*FOX      COMFAC2=(SIGXX+SIGYY)*(SIGXX+SIGYY) ;                         *FOX
      CALL DAADD(SIGXX      ,SIGYY      ,ISCRDA(  1+IDAA))              
      CALL DAADD(SIGXX      ,SIGYY      ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),COMFAC2    )                          
*FOX      BBGX=(COMFAC+FOUR*SEPX*SEPX*CONST/X*EXPFAC)/(TWO*X) ;         *FOX
      CALL DACMU(SEPX       ,ONE*FOUR       ,ISCRDA(  1+IDAA))          
      CALL DAMUL(ISCRDA(  1+IDAA),SEPX       ,ISCRDA(  2+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),CONST      ,ISCRDA(  3+IDAA))         
      CALL DADIV(ISCRDA(  3+IDAA),X          ,ISCRDA(  4+IDAA))         
      CALL DAMUL(ISCRDA(  4+IDAA),EXPFAC     ,ISCRDA(  5+IDAA))         
      CALL DACMU(X          ,ONE*TWO        ,ISCRDA(  6+IDAA))          
      CALL DAADD(COMFAC     ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DADIV(ISCRDA(  7+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),BBGX       )                          
*FOX      BBGY=(-COMFAC+FOUR*SEPY*SEPY*CONST/X*EXPFAC)/(TWO*X) ;        *FOX
      CALL DACMU(COMFAC     ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACMU(SEPY       ,ONE*FOUR       ,ISCRDA(  2+IDAA))          
      CALL DAMUL(ISCRDA(  2+IDAA),SEPY       ,ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  3+IDAA),CONST      ,ISCRDA(  4+IDAA))         
      CALL DADIV(ISCRDA(  4+IDAA),X          ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  5+IDAA),EXPFAC     ,ISCRDA(  6+IDAA))         
      CALL DACMU(X          ,ONE*TWO        ,ISCRDA(  7+IDAA))          
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DADIV(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),BBGY       )                          
        else
*FOX      BBFX=ZERO ;                                                   *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(BBFX       ,RSCRRI(100))                               
*FOX      BBFY=ZERO ;                                                   *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(BBFY       ,RSCRRI(100))                               
*FOX      BBGX=ZERO ;                                                   *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(BBGX       ,RSCRRI(100))                               
*FOX      BBGY=ZERO ;                                                   *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(BBGY       ,RSCRRI(100))                               
        endif
      else
*FOX    X=SEPX*SEPX/SIGXX+SEPY*SEPY/SIGYY ;                             *FOX
      CALL DAMUL(SEPX       ,SEPX       ,ISCRDA(  1+IDAA))              
      CALL DADIV(ISCRDA(  1+IDAA),SIGXX      ,ISCRDA(  2+IDAA))         
      CALL DAMUL(SEPY       ,SEPY       ,ISCRDA(  3+IDAA))              
      CALL DADIV(ISCRDA(  3+IDAA),SIGYY      ,ISCRDA(  4+IDAA))         
      CALL DAADD(ISCRDA(  2+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(ISCRDA(  5+IDAA),X          )                          
*FOX    FAC2=TWO*(SIGXX-SIGYY) ;                                        *FOX
      CALL DASUB(SIGXX      ,SIGYY      ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*TWO        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),FAC2       )                          
        if(dare(sigxx).lt.dare(sigyy)) then
*FOX      FAC2=TWO*(SIGYY-SIGXX) ;                                      *FOX
      CALL DASUB(SIGYY      ,SIGXX      ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*TWO        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),FAC2       )                          
        endif
*FOX    FAC=SQRT(FAC2) ;                                                *FOX
      CALL DAFUN('SQRT',FAC2       ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),FAC        )                          
*FOX    CONST=SQRPI2/FAC ;                                              *FOX
      CALL DADIC(FAC        ,ONE*SQRPI2     ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),CONST      )                          
*FOX    SIGXY=SQRT(SIGXX/SIGYY) ;                                       *FOX
      CALL DADIV(SIGXX      ,SIGYY      ,ISCRDA(  1+IDAA))              
      CALL DAFUN('SQRT',ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))              
      CALL DACOP(ISCRDA(  2+IDAA),SIGXY      )                          
*FOX    ARG1X=(SEPX/FAC) ;                                              *FOX
      CALL DADIV(SEPX       ,FAC        ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ARG1X      )                          
        if(dare(sepx).lt.zero) then
*FOX      ARG1X=-(SEPX/FAC) ;                                           *FOX
      CALL DADIV(SEPX       ,FAC        ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*(-ONE       ),ISCRDA(  2+IDAA))   
      CALL DACOP(ISCRDA(  2+IDAA),ARG1X      )                          
        endif
*FOX    ARG1Y=(SEPY/FAC) ;                                              *FOX
      CALL DADIV(SEPY       ,FAC        ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ARG1Y      )                          
        if(dare(sepy).lt.zero) then
*FOX      ARG1Y=-(SEPY/FAC) ;                                           *FOX
      CALL DADIV(SEPY       ,FAC        ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*(-ONE       ),ISCRDA(  2+IDAA))   
      CALL DACOP(ISCRDA(  2+IDAA),ARG1Y      )                          
        endif
        call errff(arg1x,arg1y,wy1,wx1)
        if(dare(x).lt.hundred) then
*FOX      EXPFAC=EXP(-X*HALF) ;                                         *FOX
      CALL DACMU(X          ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACMU(ISCRDA(  1+IDAA),ONE*HALF       ,ISCRDA(  2+IDAA))     
      CALL DAFUN('EXP ',ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))              
      CALL DACOP(ISCRDA(  3+IDAA),EXPFAC     )                          
*FOX      ARG2X=ARG1X/SIGXY ;                                           *FOX
      CALL DADIV(ARG1X      ,SIGXY      ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ARG2X      )                          
*FOX      ARG2Y=ARG1Y*SIGXY ;                                           *FOX
      CALL DAMUL(ARG1Y      ,SIGXY      ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ARG2Y      )                          
          call errff(arg2x,arg2y,wy2,wx2)
*FOX      BBFX=CONST*(WX1-EXPFAC*WX2) ;                                 *FOX
      CALL DAMUL(EXPFAC     ,WX2        ,ISCRDA(  1+IDAA))              
      CALL DASUB(WX1        ,ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))         
      CALL DAMUL(CONST      ,ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),BBFX       )                          
*FOX      BBFY=CONST*(WY1-EXPFAC*WY2) ;                                 *FOX
      CALL DAMUL(EXPFAC     ,WY2        ,ISCRDA(  1+IDAA))              
      CALL DASUB(WY1        ,ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))         
      CALL DAMUL(CONST      ,ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),BBFY       )                          
          if(dare(sepx).lt.zero) then
*FOX        BBFX=-BBFX ;                                                *FOX
      CALL DACMU(BBFX       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),BBFX       )                          
          endif
          if(dare(sepy).lt.zero) then
*FOX        BBFY=-BBFY ;                                                *FOX
      CALL DACMU(BBFY       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),BBFY       )                          
          endif
*FOX      COMFAC=SEPX*BBFX+SEPY*BBFY ;                                  *FOX
      CALL DAMUL(SEPX       ,BBFX       ,ISCRDA(  1+IDAA))              
      CALL DAMUL(SEPY       ,BBFY       ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),COMFAC     )                          
*FOX      BBGX=-(COMFAC+TWO*(EXPFAC/SIGXY-ONE))/FAC2 ;                  *FOX
      CALL DADIV(EXPFAC     ,SIGXY      ,ISCRDA(  1+IDAA))              
      CALL DACSU(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  2+IDAA))     
      CALL DACMU(ISCRDA(  2+IDAA),ONE*TWO        ,ISCRDA(  3+IDAA))     
      CALL DAADD(COMFAC     ,ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))         
      CALL DACMU(ISCRDA(  4+IDAA),ONE*(-ONE       ),ISCRDA(  5+IDAA))   
      CALL DADIV(ISCRDA(  5+IDAA),FAC2       ,ISCRDA(  6+IDAA))         
      CALL DACOP(ISCRDA(  6+IDAA),BBGX       )                          
*FOX      BBGY= (COMFAC+TWO*(EXPFAC*SIGXY-ONE))/FAC2 ;                  *FOX
      CALL DAMUL(EXPFAC     ,SIGXY      ,ISCRDA(  1+IDAA))              
      CALL DACSU(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  2+IDAA))     
      CALL DACMU(ISCRDA(  2+IDAA),ONE*TWO        ,ISCRDA(  3+IDAA))     
      CALL DAADD(COMFAC     ,ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))         
      CALL DADIV(ISCRDA(  4+IDAA),FAC2       ,ISCRDA(  5+IDAA))         
      CALL DACOP(ISCRDA(  5+IDAA),BBGY       )                          
        else
*FOX      BBFX=CONST*WX1 ;                                              *FOX
      CALL DAMUL(CONST      ,WX1        ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),BBFX       )                          
*FOX      BBFY=CONST*WY1 ;                                              *FOX
      CALL DAMUL(CONST      ,WY1        ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),BBFY       )                          
          if(dare(sepx).lt.zero) then
*FOX        BBFX=-BBFX ;                                                *FOX
      CALL DACMU(BBFX       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),BBFX       )                          
          endif
          if(dare(sepy).lt.zero) then
*FOX        BBFY=-BBFY ;                                                *FOX
      CALL DACMU(BBFY       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),BBFY       )                          
          endif
*FOX      COMFAC=SEPX*BBFX+SEPY*BBFY ;                                  *FOX
      CALL DAMUL(SEPX       ,BBFX       ,ISCRDA(  1+IDAA))              
      CALL DAMUL(SEPY       ,BBFY       ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),COMFAC     )                          
*FOX      BBGX=-(COMFAC-TWO)/FAC2 ;                                     *FOX
      CALL DACSU(COMFAC     ,ONE*TWO        ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*(-ONE       ),ISCRDA(  2+IDAA))   
      CALL DADIV(ISCRDA(  2+IDAA),FAC2       ,ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),BBGX       )                          
*FOX      BBGY= -BBGX ;                                                 *FOX
      CALL DACMU(BBGX       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),BBGY       )                          
        endif
      endif
        CALL DADAL(COMFAC2 ,1)                                                  
        CALL DADAL(COMFAC  ,1)                                                  
        CALL DADAL(WY2     ,1)                                                  
        CALL DADAL(WX2     ,1)                                                  
        CALL DADAL(WY1     ,1)                                                  
        CALL DADAL(WX1     ,1)                                                  
        CALL DADAL(ARG2Y   ,1)                                                  
        CALL DADAL(ARG2X   ,1)                                                  
        CALL DADAL(ARG1Y   ,1)                                                  
        CALL DADAL(ARG1X   ,1)                                                  
        CALL DADAL(CONST   ,1)                                                  
        CALL DADAL(FAC2    ,1)                                                  
        CALL DADAL(FAC     ,1)                                                  
        CALL DADAL(X       ,1)                                                  
        CALL DADAL(EXPFAC  ,1)                                                  
        CALL DADAL(SIGXY   ,1)                                                  
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
      double precision al,as,chi0,chid,dp1,dps,exz,sigm
      common/syos/as(6,2,npart,nele),al(6,2,npart,nele),sigm(mpa),      &
     &dps(mpa),idz(2)
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
      chd(1)='  D-X '                                                   *FOX
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
