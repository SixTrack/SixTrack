+dk nwrtcoll
! General routines:
! collimate_init()
! collimate_start_sample()
! collimate_start_turn()
! collimate_start_element()
! collimate_start_collimator()
! collimate_do_collimator()
! collimate_end_collimator()
! collimate_end_element()
! collimate_end_turn()
! collimate_end_sample()
! collimate_exit()
!
! To stop a messy future, each of these should contain calls to
! implementation specific functions: e.g. collimate_init_k2(), etc.
! These should contain the "real" code.

! In addition, these files contain:
! 1: The RNG used in collimation.
! 2: A bunch distribution generator

!>
!! collimate_init()
!! This routine is called once at the start of the simulation and
!! can be used to do any initial configuration and/or file loading.
!<
      subroutine collimate_init()
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
      integer i,ix,j,jb,jj,jx,kpz,kzz,napx0,nbeaux,nmz,nthinerr
      double precision benkcc,cbxb,cbzb,cikveb,crkveb,crxb,crzb,r0,r000,&
     &r0a,r2b,rb,rho2b,rkb,tkb,xbb,xrb,zbb,zrb
      logical lopen
+ca parpro
+ca parnum
+ca common
+ca commons
+ca commont1
+ca commondl
+ca commonxz
+ca commonta
+ca commonmn
+ca commonm1
+ca commontr
+ca beamdim
      dimension nbeaux(nbb)
+if collimat
+ca collpara
+ca dbtrthin
+ca database
+ca dbcommon
+ca dblinopt
+ca dbpencil
+ca info
+ca dbcolcom
+ca dbdaten
+ei
+if bnlelens
+ca rhicelens
+ei

+if g4collimat
      double precision g4_ecut
      integer g4_physics
+ei

      open(unit=outlun, file='colltrack.out')

      write(lout,*) '         -------------------------------'
      write(lout,*)
      write(lout,*) '          Program      C O L L T R A C K '
      write(lout,*)
      write(lout,*) '            R. Assmann           -    AB/ABP'
      write(lout,*) '            C. Bracco            -    AB/ABP'
      write(lout,*) '            V. Previtali         -    AB/ABP'
      write(lout,*) '            S. Redaelli          -    AB/OP'
      write(lout,*) '            G. Robert-Demolaize  -    BNL'
      write(lout,*) '            A. Rossi             -    AB/ABP'
      write(lout,*) '            T. Weiler            -    IEKP'
      write(lout,*) '                 CERN 2001 - 2009'
      write(lout,*)
      write(lout,*) '         -------------------------------'
      write(lout,*) 'Collimation version of Sixtrack running... 08/2009'

      write(outlun,*)
      write(outlun,*)
      write(outlun,*) '         -------------------------------'
      write(outlun,*)
      write(outlun,*) '         Program      C O L L T R A C K '
      write(outlun,*)
      write(outlun,*) '            R. Assmann       -    AB/ABP'
      write(outlun,*) '             C.Bracco        -    AB/ABP'
      write(outlun,*) '           V. Previtali      -    AB/ABP'
      write(outlun,*) '           S. Redaelli       -    AB/OP'
      write(outlun,*) '      G. Robert-Demolaize    -    BNL'
      write(outlun,*) '             A. Rossi        -    AB/ABP'
      write(outlun,*) '             T. Weiler       -    IEKP'
      write(outlun,*)
      write(outlun,*) '                 CERN 2001 - 2009'
      write(outlun,*)
      write(outlun,*) '         -------------------------------'
      write(outlun,*)
      write(outlun,*)
!
      write(lout,*) '                     R. Assmann, F. Schmidt, CERN'
      write(lout,*) '                           C. Bracco,        CERN'
      write(lout,*) '                           V. Previtali,     CERN'
      write(lout,*) '                           S. Redaelli,      CERN'
      write(lout,*) '                       G. Robert-Demolaize,  BNL'
      write(lout,*) '                           A. Rossi,         CERN'
      write(lout,*) '                           T. Weiler         IEKP'

      write(lout,*)
      write(lout,*) 'Generating particle distribution at FIRST element!'
      write(lout,*) 'Optical functions obtained from Sixtrack internal!'
      write(lout,*) 'Emittance and energy obtained from Sixtrack input!'
      write(lout,*)
      write(lout,*)
      
      write(lout,*) 'Info: Betax0   [m]    ', tbetax(1)
      write(lout,*) 'Info: Betay0   [m]    ', tbetay(1)
      write(lout,*) 'Info: Alphax0         ', talphax(1)
      write(lout,*) 'Info: Alphay0         ', talphay(1)
      write(lout,*) 'Info: Orbitx0  [mm]   ', torbx(1)
      write(lout,*) 'Info: Orbitxp0 [mrad] ', torbxp(1)
      write(lout,*) 'Info: Orbity0  [mm]   ', torby(1)
      write(lout,*) 'Info: Orbitpy0 [mrad] ', torbyp(1)
      write(lout,*) 'Info: Emitx0_dist [um]', remitx_dist
      write(lout,*) 'Info: Emity0_dist [um]', remity_dist
      write(lout,*) 'Info: Emitx0_collgap [um]', remitx_collgap
      write(lout,*) 'Info: Emity0_collgap [um]', remity_collgap
      write(lout,*) 'Info: E0       [MeV]  ', e0
      write(lout,*)
      write(lout,*)
!
      myemitx0_dist = remitx_dist*1d-6
      myemity0_dist = remity_dist*1d-6
      myemitx0_collgap = remitx_collgap*1d-6
      myemity0_collgap = remity_collgap*1d-6

      myalphax = talphax(1)
      myalphay = talphay(1)
      mybetax  = tbetax(1)
      mybetay  = tbetay(1)

!07-2006      myenom   = e0
!      MYENOM   = 1.001*E0
!
      if (myemitx0_dist.le.0.d0 .or. myemity0_dist.le.0.d0
     &.or. myemitx0_collgap.le.0.d0 .or. myemity0_collgap.le.0.d0) then
        write(lout,*)
     &       'ERR> EMITTANCES NOT DEFINED! CHECK COLLIMAT BLOCK!'
        write(lout,*)"ERR> EXPECTED FORMAT OF LINE 9 IN COLLIMAT BLOCK:"
        write(lout,*)
     & "emitnx0_dist  emitny0_dist  emitnx0_collgap  emitny0_collgap"

        write(lout,*) "ERR> ALL EMITTANCES SHOULD BE NORMALIZED.",
     & "FIRST PUT EMITTANCE FOR DISTRIBTION GENERATION, ",
     & "THEN FOR COLLIMATOR POSITION ETC. UNITS IN [MM*MRAD]."
        write(lout,*) "ERR> EXAMPLE:"
        write(lout,*) "2.5 2.5 3.5 3.5"
        call prror(-1)
      endif
!
!++  Calculate the gammas
!
      mygammax = (1d0+myalphax**2)/mybetax
      mygammay = (1d0+myalphay**2)/mybetay
!
!++  Number of points and generate distribution
!
!GRD SEMI-AUTOMATIC INPUT
!      NLOOP=10
!      MYNEX=6.003
!      MYDEX=0.0015
!      MYNEY=6.003
!      MYDEY=0.0015
!      DO_COLL=1
!      NSIG_PRIM=5.
!      NSIG_SEC=6.
      rselect=64
!
      write(lout,*) 'INFO>  NLOOP     = ', nloop
      write(lout,*) 'INFO>  DO_THISDIS     = ', do_thisdis
      write(lout,*) 'INFO>  MYNEX     = ', mynex
      write(lout,*) 'INFO>  MYDEX     = ', mdex
      write(lout,*) 'INFO>  MYNEY     = ', myney
      write(lout,*) 'INFO>  MYDEY     = ', mdey
      write(lout,*) 'INFO>  FILENAME_DIS     = ', filename_dis
      write(lout,*) 'INFO>  ENERROR     = ', enerror
      write(lout,*) 'INFO>  BUNCHLENGTH     = ', bunchlength
      write(lout,*) 'INFO>  RSELECT   = ', int(rselect)
      write(lout,*) 'INFO>  DO_COLL   = ', do_coll
!APRIL2005
!+if cr
!      write(lout,*) 'INFO>  NSIG_PRIM = ', nsig_prim
!+ei
!+if .not.cr
!      write(*,*) 'INFO>  NSIG_PRIM = ', nsig_prim
!+ei
!+if cr
!      write(lout,*) 'INFO>  NSIG_SEC  = ', nsig_sec
!+ei
!+if .not.cr
!      write(*,*) 'INFO>  NSIG_SEC  = ', nsig_sec
!+ei
      write(lout,*) 'INFO>  DO_NSIG   = ', do_nsig
      write(lout,*) 'INFO>  NSIG_TCP3    = ', nsig_tcp3
      write(lout,*) 'INFO>  NSIG_TCSG3   = ', nsig_tcsg3
      write(lout,*) 'INFO>  NSIG_TCSM3   = ', nsig_tcsm3
      write(lout,*) 'INFO>  NSIG_TCLA3   = ', nsig_tcla3
      write(lout,*) 'INFO>  NSIG_TCP7    = ', nsig_tcp7
      write(lout,*) 'INFO>  NSIG_TCSG7   = ', nsig_tcsg7
      write(lout,*) 'INFO>  NSIG_TCSM7   = ', nsig_tcsm7
      write(lout,*) 'INFO>  NSIG_TCLA7   = ', nsig_tcla7
      write(lout,*) 'INFO>  NSIG_TCLP    = ', nsig_tclp
      write(lout,*) 'INFO>  NSIG_TCLI    = ', nsig_tcli
!      write(lout,*) 'INFO>  NSIG_TCTH    = ', nsig_tcth
!      write(lout,*) 'INFO>  NSIG_TCTV    = ', nsig_tctv
      write(lout,*) 'INFO>  NSIG_TCTH1   = ', nsig_tcth1
      write(lout,*) 'INFO>  NSIG_TCTV1   = ', nsig_tctv1
      write(lout,*) 'INFO>  NSIG_TCTH2   = ', nsig_tcth2
      write(lout,*) 'INFO>  NSIG_TCTV2   = ', nsig_tctv2
      write(lout,*) 'INFO>  NSIG_TCTH5   = ', nsig_tcth5
      write(lout,*) 'INFO>  NSIG_TCTV5   = ', nsig_tctv5
      write(lout,*) 'INFO>  NSIG_TCTH8   = ', nsig_tcth8
      write(lout,*) 'INFO>  NSIG_TCTV8   = ', nsig_tctv8
!
      write(lout,*) 'INFO>  NSIG_TCDQ    = ', nsig_tcdq
      write(lout,*) 'INFO>  NSIG_TCSTCDQ = ', nsig_tcstcdq
      write(lout,*) 'INFO>  NSIG_TDI     = ', nsig_tdi
      write(lout,*) 'INFO>  NSIG_TCXRP   = ', nsig_tcxrp
      write(lout,*) 'INFO>  NSIG_TCRYP   = ', nsig_tcryo
!APRIL2005
!SEPT2005
      write(lout,*)
      write(lout,*) 'INFO> INPUT PARAMETERS FOR THE SLICING:'
      write(lout,*)
      write(lout,*) 'INFO>  N_SLICES    = ', n_slices
      write(lout,*) 'INFO>  SMIN_SLICES = ',smin_slices
      write(lout,*) 'INFO>  SMAX_SLICES = ',smax_slices
      write(lout,*) 'INFO>  RECENTER1   = ',recenter1
      write(lout,*) 'INFO>  RECENTER2   = ',recenter2
      write(lout,*)
      write(lout,*) 'INFO>  FIT1_1   = ',fit1_1
      write(lout,*) 'INFO>  FIT1_2   = ',fit1_2
      write(lout,*) 'INFO>  FIT1_3   = ',fit1_3
      write(lout,*) 'INFO>  FIT1_4   = ',fit1_4
      write(lout,*) 'INFO>  FIT1_5   = ',fit1_5
      write(lout,*) 'INFO>  FIT1_6   = ',fit1_6
      write(lout,*) 'INFO>  SCALING1 = ',ssf1
      write(lout,*)
      write(lout,*) 'INFO>  FIT2_1   = ',fit2_1
      write(lout,*) 'INFO>  FIT2_2   = ',fit2_2
      write(lout,*) 'INFO>  FIT2_3   = ',fit2_3
      write(lout,*) 'INFO>  FIT2_4   = ',fit2_4
      write(lout,*) 'INFO>  FIT2_5   = ',fit2_5
      write(lout,*) 'INFO>  FIT2_6   = ',fit2_6
      write(lout,*) 'INFO>  SCALING2 = ',ssf2
      write(lout,*)

!SEPT2005
!
! HERE WE CHECK IF THE NEW INPUT IS READ CORRECTLY
!
      write(lout,*) 'INFO>  EMITXN0_DIST      = ', emitnx0_dist
      write(lout,*) 'INFO>  EMITYN0_DIST      = ', emitny0_dist
      write(lout,*) 'INFO>  EMITXN0_COLLGAP   = ', emitnx0_collgap
      write(lout,*) 'INFO>  EMITYN0_COLLGAP   = ', emitny0_collgap
      write(lout,*)
      write(lout,*) 'INFO>  DO_SELECT         = ', do_select
      write(lout,*) 'INFO>  DO_NOMINAL        = ', do_nominal
      write(lout,*) 'INFO>  RND_SEED          = ', rnd_seed
      write(lout,*) 'INFO>  DOWRITE_DIST      = ', dowrite_dist
      write(lout,*) 'INFO>  NAME_SEL          = ', name_sel
      write(lout,*) 'INFO>  DO_ONESIDE        = ', do_oneside
      write(lout,*) 'INFO>  DOWRITE_IMPACT    = ', dowrite_impact
      write(lout,*) 'INFO>  DOWRITE_SECONDARY = ', dowrite_secondary
      write(lout,*) 'INFO>  DOWRITE_AMPLITUDE = ', dowrite_amplitude
      write(lout,*)
      write(lout,*) 'INFO>  XBEAT             = ', xbeat
      write(lout,*) 'INFO>  XBEATPHASE        = ', xbeatphase
      write(lout,*) 'INFO>  YBEAT             = ', ybeat
      write(lout,*) 'INFO>  YBEATPHASE        = ', ybeatphase
      write(lout,*)
      write(lout,*) 'INFO>  C_RMSTILT_PRIM     = ', c_rmstilt_prim
      write(lout,*) 'INFO>  C_RMSTILT_SEC      = ', c_rmstilt_sec
      write(lout,*) 'INFO>  C_SYSTILT_PRIM     = ', c_systilt_prim
      write(lout,*) 'INFO>  C_SYSTILT_SEC      = ', c_systilt_sec
      write(lout,*) 'INFO>  C_RMSOFFSET_PRIM   = ', c_rmsoffset_prim
      write(lout,*) 'INFO>  C_SYSOFFSET_PRIM   = ', c_sysoffset_prim
      write(lout,*) 'INFO>  C_RMSOFFSET_SEC    = ', c_rmsoffset_sec
      write(lout,*) 'INFO>  C_SYSOFFSET_SEC    = ', c_sysoffset_sec
      write(lout,*) 'INFO>  C_OFFSETTITLT_SEED = ', c_offsettilt_seed
      write(lout,*) 'INFO>  C_RMSERROR_GAP     = ', c_rmserror_gap
      write(lout,*) 'INFO>  DO_MINGAP          = ', do_mingap
      write(lout,*)
      write(lout,*) 'INFO>  RADIAL            = ', radial
      write(lout,*) 'INFO>  NR                = ', nr
      write(lout,*) 'INFO>  NDR               = ', ndr
      write(lout,*)
      write(lout,*) 'INFO>  DRIFTSX           = ', driftsx
      write(lout,*) 'INFO>  DRIFTSY           = ', driftsy
      write(lout,*) 'INFO>  CUT_INPUT         = ', cut_input
      write(lout,*) 'INFO>  SYSTILT_ANTISYMM  = ', systilt_antisymm
      write(lout,*)
      write(lout,*) 'INFO>  IPENCIL           = ', ipencil
      write(lout,*) 'INFO>  PENCIL_OFFSET     = ', pencil_offset
      write(lout,*) 'INFO>  PENCIL_RMSX       = ', pencil_rmsx
      write(lout,*) 'INFO>  PENCIL_RMSY       = ', pencil_rmsy
      write(lout,*) 'INFO>  PENCIL_DISTR      = ', pencil_distr
      write(lout,*)
      write(lout,*) 'INFO>  COLL_DB           = ', coll_db
      write(lout,*) 'INFO>  IBEAM             = ', ibeam
      write(lout,*)
      write(lout,*) 'INFO>  DOWRITETRACKS     = ', dowritetracks
      write(lout,*)
      write(lout,*) 'INFO>  CERN              = ', cern
      write(lout,*)
      write(lout,*) 'INFO>  CASTORDIR     = ', castordir
      write(lout,*)
      write(lout,*) 'INFO>  JOBNUMBER     = ', jobnumber
      write(lout,*)
      write(lout,*) 'INFO>  CUTS     = ', sigsecut2, sigsecut3
      write(lout,*)
!
      mynp = nloop*napx
!
      napx00 = napx
!
      write(lout,*) 'INFO>  NAPX     = ', napx, mynp
      write(lout,*) 'INFO>  Sigma_x0 = ', sqrt(mybetax*myemitx0_dist)
      write(lout,*) 'INFO>  Sigma_y0 = ', sqrt(mybetay*myemity0_dist)
!
! HERE WE SET THE MARKER FOR INITIALIZATION:
!
      firstrun = .true.
!
! ...and here is implemented colltrack's beam distribution:
!
!
!++  Initialize random number generator
!
!      IF (FIRSTRUN) THEN
        if (rnd_seed.eq.0) rnd_seed = mclock_liar()
        if (rnd_seed.lt.0) rnd_seed = abs(rnd_seed)
        rnd_lux = 3
        rnd_k1  = 0
        rnd_k2  = 0
        call rluxgo(rnd_lux, rnd_seed, rnd_k1, rnd_k2)
        write(lout,*)
        write(outlun,*) 'INFO>  rnd_seed: ', rnd_seed
!      ENDIF
!GRD-SR, 09-02-2006
!Call distribution routines only if collimation block is in fort.3, otherwise
!the standard sixtrack would be prevented by the 'stop' command
      if(do_coll) then
!GRD-SR
      if (radial) then
         call    makedis_radial(mynp, myalphax, myalphay, mybetax,
     &        mybetay, myemitx0_dist, myemity0_dist, myenom, nr, ndr,
     &        myx, myxp, myy, myyp, myp, mys)
      else
         if (do_thisdis.eq.1) then
            call makedis(mynp, myalphax, myalphay, mybetax, mybetay,
     &           myemitx0_dist, myemity0_dist,
     &           myenom, mynex, mdex, myney, mdey,
     &           myx, myxp, myy, myyp, myp, mys)
         elseif(do_thisdis.eq.2) then
            call makedis_st(mynp, myalphax, myalphay, mybetax, mybetay,
     &           myemitx0_dist, myemity0_dist,
     &           myenom, mynex, mdex, myney, mdey,
     &           myx, myxp, myy, myyp, myp, mys)
         elseif(do_thisdis.eq.3) then
            call makedis_de(mynp, myalphax, myalphay, mybetax, mybetay,
     &           myemitx0_dist, myemity0_dist,
     &           myenom, mynex, mdex, myney, mdey,
     &           myx, myxp, myy, myyp, myp, mys,enerror,bunchlength)
         elseif(do_thisdis.eq.4) then
            call readdis(filename_dis,
     &           mynp, myx, myxp, myy, myyp, myp, mys)
         elseif(do_thisdis.eq.5) then
            call makedis_ga(mynp, myalphax, myalphay, mybetax,
     &           mybetay, myemitx0_dist, myemity0_dist,
     &           myenom, mynex, mdex, myney, mdey,
     &           myx, myxp, myy, myyp, myp, mys,
     &           enerror, bunchlength )
         elseif(do_thisdis.eq.6) then
            call readdis_norm(filename_dis, 
     &           mynp, myalphax, myalphay, mybetax, mybetay,
     &           myemitx0_dist, myemity0_dist, myenom, 
     &           myx, myxp, myy, myyp, myp, mys, 
     &           enerror, bunchlength)

         else
            write(lout,*) 'INFO> review your distribution parameters !!'
            call prror(-1)
         endif
!
       endif
!
!GRD-SR,09-02-2006
      endif
!GRD-SR
!++  Reset distribution for pencil beam
!
       if (ipencil.gt.0) then
         write(lout,*) 'WARN>  Distributions reset to pencil beam!'
         write(lout,*)
         write(outlun,*) 'WARN>  Distributions reset to pencil beam!'
         do j = 1, mynp
            myx(j)  = 0d0
            myxp(j) = 0d0
            myy(j)  = 0d0
            myyp(j) = 0d0
         end do
       endif

!++  Optionally write the generated particle distribution
      open(unit=52,file='dist0.dat')
       if (dowrite_dist) then
        do j = 1, mynp
          write(52,'(6(1X,E15.7))') myx(j), myxp(j), myy(j), myyp(j),   &
!     SR, 11-08-2005
     &          mys(j), myp(j)
        end do
       endif
      close(52)
!
!++  Initialize efficiency array
!
      do i=1,iu
      sum_ax(i)   = 0d0
      sqsum_ax(i) = 0d0
      sum_ay(i)   = 0d0
      sqsum_ay(i) = 0d0
      nampl(i)    = 0d0
      sampl(i)    = 0d0
      end do
!
      nspx = 0d0
      nspy = 0d0
!
      np0  = mynp
!
      ax0  = myalphax
      bx0  = mybetax
      mux0 = mux(1)
      ay0  = myalphay
      by0  = mybetay
      muy0 = muy(1)
      iturn = 1
      ie    = 1
      n_tot_absorbed = 0
      
      if (int(mynp/napx00) .eq. 0) then
         write (lout,*) ""
         write (lout,*) "********************************************"
         write (lout,*) "Error in setting up collimation tracking:"
         write (lout,*) "Number of samples is zero!"
         write (lout,*) "Did you forget the COLL block in fort.3?"
         write (lout,*) "If you want to do standard (not collimation)"//
     &                  " tracking, please use the standard SixTrack."
         write (lout,*) "Value of do_coll = ", do_coll
         write (lout,*) "Value of mynp    = ", mynp
         write (lout,*) "Value of napx00  = ", napx00
         write (lout,*) "********************************************"
         call prror(-1)
      endif

!++  Read collimator database
      call readcollimator

!Then do any implementation specific initial loading
+if collimate_k2
      call collimate_init_k2
+ei
+if merlinscatter
      call collimate_init_merlin
+ei
+if g4collimat
!! This function lives in the G4Interface.cpp file in the g4collimat folder
!! Accessed by linking libg4collimat.a
!! Set the energy cut at 70% - i.e. 30% energy loss

      g4_ecut = 0.7

!! Select the physics engine to use
!! 0 = FTFP_BERT
!! 1 = QGSP_BERT

      g4_physics = 0
      call g4_collimation_init(e0, rnd_seed, g4_ecut, g4_physics)
+ei
      end

!>
!! collimate_start_sample()
!! This routine is called from trauthin before each sample
!! is injected into thin 6d
!<
      subroutine collimate_start_sample(nsample)
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
      integer i,ix,j,jj,jx,kpz,kzz,napx0,nbeaux,nmz,nthinerr,nsample
      double precision benkcc,cbxb,cbzb,cikveb,crkveb,crxb,crzb,r0,r000,&
     &r0a,r2b,rb,rho2b,rkb,tkb,xbb,xrb,zbb,zrb
      logical lopen
+ca parpro
+ca parnum
+ca common
+ca commons
+ca commont1
+ca commondl
+ca commonxz
+ca commonta
+ca commonmn
+ca commonm1
+ca commontr
+ca beamdim
      dimension nbeaux(nbb)
+if collimat
+ca collpara
+ca dbtrthin
+ca database
+ca dbcommon
+ca dblinopt
+ca dbpencil
+ca info
+ca dbthin6d
+ca dbcolcom
+ei
+if bnlelens
+ca rhicelens
+ei
+if cr
+ca crco
+ei
!GRD
      j = nsample
      samplenumber=j

! HERE WE OPEN ALL THE NEEDED OUTPUT FILES

      open(unit=99,file='betatron.dat')

      open(unit=42, file='beta_beat.dat')
      write(42,*)                                                       &
     &'# 1=s 2=bx/bx0 3=by/by0 4=sigx0 5=sigy0 6=crot 7=acalc'

      open(unit=44, file='survival.dat') ! RB, DM: 2014 bug fix
      write(44,*)                                                       &
     &'# 1=turn 2=n_particle'

      open(unit=43, file='collgaps.dat')
      if(firstrun) write(43,*)                                          &
     &'# ID name  angle[rad]  betax[m]  betay[m] ',                     &
     &'halfgap[m]  Material  Length[m]  sigx[m]  sigy[m] ',             &
     &'tilt1[rad] tilt2[rad] nsig'

!      if (dowrite_impact) then
!        open(unit=46, file='coll_impact.dat')
!        write(46,*)                                                     &
!     &'# 1=sample 2=iturn 3=icoll 4=nimp 5=nabs 6=imp_av 7=imp_sig'
!      endif
!
      open(unit=40, file='collimator-temp.db')
!
!      open(unit=47, file='tertiary.dat')
!      write(47,*)                                                       &
!     &'# 1=x 2=xp 3=y 4=yp 5=p 6=Ax 7=Axd 8=Ay 9=Ar 10=Ard'
!
!      if (dowrite_secondary) then
!        open(unit=48, file='secondary.dat')
!        write(48,'(2a)')                                                &
!     &'# 1=x 2=xp 3=y 4=yp 5=p 6=Ax 7=Axd 8=Ay 9=Ar 10=Ard'
!      endif
!
! TW06/08 added ouputfile for real collimator settings (incluing slicing, ...)
      open(unit=55, file='collsettings.dat')
      if(firstrun) write(55,*)                                          &
     &'# name  slicenumber  halfgap[m]  gap_offset[m] ',                &
     &'tilt jaw1[rad]  tilt jaw2[rad] length[m] material'               &
! TW06/08
      if (dowrite_impact) then
        open(unit=49,file='impact.dat')
        write(49,*)                                                     &
     &'# 1=impact 2=divergence'
      endif


      if (dowritetracks) then
!GRD SPECIAL FILE FOR SECONDARY HALO
        if (cern) then
        open(unit=41,file='stuff')
        write(41,*) samplenumber
        close(41)
        open(unit=41,file='stuff')
        read(41,*) smpl
        close(41)
        pfile(1:8) = 'tracks2.'
        if(samplenumber.le.9) then
           pfile(9:9) = smpl
           pfile(10:13) = '.dat'
        elseif(samplenumber.gt.9.and.samplenumber.le.99) then
           pfile(9:10) = smpl
           pfile(11:14) = '.dat'
        elseif(samplenumber.gt.99.and.                                  &
     &samplenumber.le.int(mynp/napx00)) then
           pfile(9:11) = smpl
           pfile(12:15) = '.dat'
        endif

        if(samplenumber.le.9)                                           &
     &open(unit=38,file=pfile(1:13))

        if(samplenumber.gt.9.and.samplenumber.le.99)                    &
     &open(unit=38,file=pfile(1:14))

        if(samplenumber.gt.99.and.                                      &
     &samplenumber.le.int(mynp/napx00))                                 &
     &open(unit=38,file=pfile(1:15))
        else
          open(unit=38,file='tracks2.dat')
        endif !end if (cern)

        if(firstrun) write(38,*)                                        &
     &'# 1=name 2=turn 3=s 4=x 5=xp 6=y 7=yp 8=DE/E 9=type'

!AUGUST2006:write pencul sheet beam coordiantes to file ---- TW
        open(unit=9997, file='pencilbeam_distr.dat')
        if(firstrun) write(9997,*) 'x    xp    y      yp'

      endif !end if (dowritetracks) then


!GRD-SR,09-02-2006 => new series of output controlled by the 'dowrite_impact flag
      if(do_select) then
        open(unit=45, file='coll_ellipse.dat')
        if (firstrun) then
           write(45,*)                                                  &
     &          '#  1=name 2=x 3=y 4=xp 5=yp 6=E 7=s 8=turn 9=halo ',   &
     & '10=nabs_type'
        endif
      endif

      if(dowrite_impact) then
        open(unit=46, file='all_impacts.dat')
        open(unit=47, file='all_absorptions.dat')
        open(unit=48, file='FLUKA_impacts.dat')
! RB: adding output files FLUKA_impacts_all.dat and Coll_Scatter.dat
        open(unit=4801, file='FLUKA_impacts_all.dat')
        open(unit=3998, file='Coll_Scatter.dat')
        open(unit=39, file='FirstImpacts.dat')

        if (firstrun) then
          write(46,'(a)') '# 1=name 2=turn 3=s'
          write(47,'(a)') '# 1=name 2=turn 3=s'
          write(48,'(a)')                                               &
     &'# 1=icoll 2=c_rotation 3=s 4=x 5=xp 6=y 7=yp 8=nabs 9=np 10=turn'
          write(39,*)                                                   &
     &     '%1=name,2=iturn, 3=icoll, 4=nabs, 5=s_imp[m], 6=s_out[m], ',&
     &     '7=x_in(b!)[m], 8=xp_in, 9=y_in, 10=yp_in, ',                &
     &     '11=x_out [m], 12=xp_out, 13=y_out, 14=yp_out'

! RB: write headers in new output files
          write(4801,'(a)')                                             &
     &'# 1=icoll 2=c_rotation 3=s 4=x 5=xp 6=y 7=yp 8=nabs 9=np 10=turn'
          write(3998,*)
     &     "#1=icoll, 2=iturn, 3=np, 4=nabs (1:Nuclear-Inelastic,2:Nucle
     &ar-Elastic,3:pp-Elastic,4:Single-Diffractive,5:Coulomb), 5=dp, 6=d
     &x', 7=dy'"
        endif
      endif

      if(name_sel(1:3).eq.'COL') then
        open(unit=555, file='RHIClosses.dat')
        if(firstrun) write(555,'(a)')                                   &
     &'# 1=name 2=turn 3=s 4=x 5=xp 6=y 7=yp 8=dp/p 9=type'
      endif


!++  Copy new particles to tracking arrays. Also add the orbit offset at
!++  start of ring!

            do i = 1, napx00
              xv(1,i)  = 1d3*myx(i+(j-1)*napx00)  +torbx(1)              !hr08
              yv(1,i)  = 1d3*myxp(i+(j-1)*napx00) +torbxp(1)             !hr08
              xv(2,i)  = 1d3*myy(i+(j-1)*napx00)  +torby(1)              !hr08
              yv(2,i)  = 1d3*myyp(i+(j-1)*napx00) +torbyp(1)             !hr08

!JULY2005 assignation of the proper bunch length
              sigmv(i) = mys(i+(j-1)*napx00)
              ejv(i)   = myp(i+(j-1)*napx00)

!GRD FOR NOT FAST TRACKING ONLY
              ejfv(i)=sqrt(ejv(i)**2-pma**2)                             !hr08
              rvv(i)=(ejv(i)*e0f)/(e0*ejfv(i))
              dpsv(i)=(ejfv(i)-e0f)/e0f
              oidpsv(i)=one/(one+dpsv(i))
              dpsv1(i)=(dpsv(i)*c1e3)*oidpsv(i)                          !hr08

              nlostp(i)=i
              do ieff =1, numeff
                counted_r(i,ieff) = 0
                counted_x(i,ieff) = 0
                counted_y(i,ieff) = 0

                do ieffdpop =1, numeffdpop
                  counted2d(i,ieff,ieffdpop) = 0
                end do
              end do

              do ieffdpop =1, numeffdpop
                counteddpop(i,ieffdpop) = 0
              end do

            end do

!!!!!!!!!!!!!!!!!!!!!!START THIN6D CUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!++  Some initialization
      do i = 1, numeff
        rsig(i) = (dble(i)/2d0 - 0.5d0) + 5d0                           !hr08
      enddo

      dpopbins(1)= 1d-4

      do i = 2, numeffdpop
        dpopbins(i)= dble(i-1)*4d-4
      enddo

      firstcoll = .true.


!GRD HERE WE NEED TO INITIALIZE SOME COLLIMATION PARAMETERS
      napx = napx00
      do j = 1, napx
         part_hit_pos(j)    = 0
         part_hit_turn(j)   = 0
         part_abs_pos(j)    = 0
         part_abs_turn(j)   = 0
         part_select(j) = 1
         part_indiv(j)  = -1e-6
         part_linteract(j) = 0d0
         part_impact(j) = 0
!      enddo
!++ Moved initialization to the start of EACH set, RA/GRD 14/6/04
!      do j=1,napx
        tertiary(j)=0
        secondary(j)=0
        other(j)=0
        nabs_type(j) = 0
!GRD HERE WE INITIALIZE THE VALUES OF IPART(j)
        ipart(j) = j
        flukaname(j) = 0
      end do
!GRD



!++  This we only do once, for the first call to this routine. Numbers
!++  are saved in memory to use exactly the same info for each sample.
!++  COMMON block to decide for first usage and to save coll info.
      if (firstrun) then
      !Reading of collimation database moved to subroutine collimate_init
+if beamgas
!YIL call beam gas initiation routine
      call beamGasInit(myenom)
+ei

      write(lout,*) 'number of collimators', db_ncoll
      do icoll = 1, db_ncoll
         write(lout,*) 'COLLIMATOR', icoll, ' ', db_name1(icoll)
         write(lout,*) 'collimator', icoll, ' ', db_name2(icoll)
      end do

!******write settings for alignment error in colltrack.out file
      write(outlun,*) ' '
      write(outlun,*) 'Alignment errors settings (tilt, offset,...)'
      write(outlun,*) ' '
      write(outlun,*) 'SETTING> c_rmstilt_prim   : ', c_rmstilt_prim
      write(outlun,*) 'SETTING> c_rmstilt_sec    : ', c_rmstilt_sec
      write(outlun,*) 'SETTING> c_systilt_prim   : ', c_systilt_prim
      write(outlun,*) 'SETTING> c_systilt_sec    : ', c_systilt_sec
      write(outlun,*) 'SETTING> c_rmsoffset_prim : ', c_rmsoffset_prim
      write(outlun,*) 'SETTING> c_rmsoffset_sec  : ', c_rmsoffset_sec
      write(outlun,*) 'SETTING> c_sysoffset_prim : ', c_sysoffset_prim
      write(outlun,*) 'SETTING> c_sysoffset_sec  : ', c_sysoffset_sec
      write(outlun,*) 'SETTING> c_offsettilt seed: ', c_offsettilt_seed
      write(outlun,*) 'SETTING> c_rmserror_gap   : ', c_rmserror_gap
      write(outlun,*) 'SETTING> do_mingap        : ', do_mingap
      write(outlun,*) ' '

!     TW - 01/2007
!     added offset and random_seed for tilt and offset
!*****intialize random generator with offset_seed
      c_offsettilt_seed = abs(c_offsettilt_seed)
      rnd_lux = 3
      rnd_k1  = 0
      rnd_k2  = 0
      call rluxgo(rnd_lux, c_offsettilt_seed, rnd_k1, rnd_k2)
!      write(outlun,*) 'INFO>  c_offsettilt seed: ', c_offsettilt_seed

! reset counter to assure starting at the same position in case of
! using rndm5 somewhere else in the code before
      zbv = rndm5(1)

!++  Generate random tilts (Gaussian distribution plus systematic)
!++  Do this only for the first call of this routine (first sample)
!++  Keep all collimator database info and errors in memeory (COMMON
!++  block) in order to re-use exactly the same information for every
!++  sample.
         if (c_rmstilt_prim.gt.0.d0 .or. c_rmstilt_sec.gt.0.d0 .or.     &!hr08
     &        c_systilt_prim.ne.0.d0 .or. c_systilt_sec.ne.0.d0) then    !hr08
            do icoll = 1, db_ncoll
               if (db_name1(icoll)(1:3).eq.'TCP') then
                  c_rmstilt = c_rmstilt_prim
                  c_systilt = c_systilt_prim
               else
                  c_rmstilt = c_rmstilt_sec
                  c_systilt = c_systilt_sec
               endif
               db_tilt(icoll,1) = c_systilt+c_rmstilt*myran_gauss(3d0)
               if (systilt_antisymm) then
                  db_tilt(icoll,2) =                                    &
     &                 -1d0*c_systilt+c_rmstilt*myran_gauss(3d0)
!    &                 c_rmstilt*myran_gauss(3d0)-c_systilt              !hr01
               else
                  db_tilt(icoll,2) =                                    &
     &                 c_systilt+c_rmstilt*myran_gauss(3d0)
               endif
               write(outlun,*) 'INFO>  Collimator ', db_name1(icoll),   &
     &              ' jaw 1 has tilt [rad]: ', db_tilt(icoll,1)
               write(outlun,*) 'INFO>  Collimator ', db_name1(icoll),   &
     &              ' jaw 2 has tilt [rad]: ', db_tilt(icoll,2)
            end do
         endif

!++  Generate random offsets (Gaussian distribution plus systematic)
!++  Do this only for the first call of this routine (first sample)
!++  Keep all collimator database info and errors in memeory (COMMON
!++  block) in order to re-use exactly the same information for every
!++  sample and throughout a all run.
         if (c_sysoffset_prim.ne.0.d0 .or. c_sysoffset_sec.ne.0.d0 .or. &!hr08
     &        c_rmsoffset_prim.gt.0.d0.or.c_rmsoffset_sec.gt.0.d0) then  !hr08
            do icoll = 1, db_ncoll 
               if (db_name1(icoll)(1:3).eq.'TCP') then
                  db_offset(icoll) = c_sysoffset_prim +                 &
     &                 c_rmsoffset_prim*myran_gauss(3d0)
               else
                  db_offset(icoll) = c_sysoffset_sec +                  &
     &                 c_rmsoffset_sec*myran_gauss(3d0)
               endif
               write(outlun,*) 'INFO>  offset: ', db_name1(icoll),      &
     &              db_offset(icoll)
            end do
         endif

!++  Generate random offsets (Gaussian distribution)
!++  Do this only for the first call of this routine (first sample)
!++  Keep all collimator database info and errors in memeory (COMMON
!++  block) in order to re-use exactly the same information for every
!++  sample and throughout a all run.
!         if (c_rmserror_gap.gt.0.) then
!            write(outlun,*) 'INFO> c_rmserror_gap = ',c_rmserror_gap
            do icoll = 1, db_ncoll 
               gap_rms_error(icoll) = c_rmserror_gap * myran_gauss(3d0)
               write(outlun,*) 'INFO>  gap_rms_error: ',                &
     &              db_name1(icoll),gap_rms_error(icoll)
            end do
!---- creating a file with beta-functions at TCP/TCS 
         open(unit=10000, file='twisslike.out')
         open(unit=10001, file='sigmasettings.out')
         mingap = 20

         do j=1,iu
! this transformation gives the right marker/name to the corresponding 
! beta-dunctions or vice versa ;)
            if(ic(j).le.nblo) then
               do jb=1,mel(ic(j))
                  myix=mtyp(ic(j),jb)
               enddo
            else
               myix=ic(j)-nblo
            endif

! Using same code-block as below to evalute the collimator opening
! for each collimator, this is needed to get the smallest collimator gap
! in principal only looking for primary and secondary should be enough
! JULY 2008 added changes (V6.503) for names in TCTV -> TCTVA and TCTVB 
! both namings before and after V6.503 can be used 
            if ( bez(myix)(1:2).eq.'TC'                                 &
     &           .or. bez(myix)(1:2).eq.'tc'                            &
     &           .or. bez(myix)(1:2).eq.'TD'                            &
     &           .or. bez(myix)(1:2).eq.'td'                            &
     &           .or. bez(myix)(1:3).eq.'COL'                           &
     &           .or. bez(myix)(1:3).eq.'col') then
               if(bez(myix)(1:3).eq.'TCP' .or.                          &
     &              bez(myix)(1:3).eq.'tcp') then
                  if(bez(myix)(7:9).eq.'3.B' .or.                       &
     &                 bez(myix)(7:9).eq.'3.b') then
                     nsig = nsig_tcp3
                  else
                     nsig = nsig_tcp7
                  endif
               elseif(bez(myix)(1:4).eq.'TCSG' .or.                     &
     &                 bez(myix)(1:4).eq.'tcsg') then
                  if(bez(myix)(8:10).eq.'3.B' .or.                      &
     &                 bez(myix)(8:10).eq.'3.b' .or.                    &
     &                 bez(myix)(9:11).eq.'3.B' .or.                    &
     &                 bez(myix)(9:11).eq.'3.b') then
                     nsig = nsig_tcsg3
                  else
                     nsig = nsig_tcsg7
                  endif
                  if((bez(myix)(5:6).eq.'.4'.and.bez(myix)(8:9).eq.'6.')&
     &                 ) then
                     nsig = nsig_tcstcdq
                  endif
               elseif(bez(myix)(1:4).eq.'TCSP' .or.                        &
     &                 bez(myix)(1:4).eq.'tcsp') then
                  if(bez(myix)(9:11).eq.'6.B'.or.
     &                 bez(myix)(9:11).eq.'6.b') then
                     nsig = nsig_tcstcdq
                  endif
               elseif(bez(myix)(1:4).eq.'TCSM' .or.                     &
     &                 bez(myix)(1:4).eq.'tcsm') then
                  if(bez(myix)(8:10).eq.'3.B' .or.                      &
     &                 bez(myix)(8:10).eq.'3.b' .or.                    &
     &                 bez(myix)(9:11).eq.'3.B' .or.                    &
     &                 bez(myix)(9:11).eq.'3.b') then
                     nsig = nsig_tcsm3
                  else
                     nsig = nsig_tcsm7
                  endif
               elseif(bez(myix)(1:4).eq.'TCLA' .or.                     &
     &                 bez(myix)(1:4).eq.'tcla') then
                  if(bez(myix)(9:11).eq.'7.B' .or.                      &
     &                 bez(myix)(9:11).eq.'7.b') then
                     nsig = nsig_tcla7
                  else
                     nsig = nsig_tcla3
                  endif
               elseif(bez(myix)(1:4).eq.'TCDQ' .or.                     &
     &                 bez(myix)(1:4).eq.'tcdq') then
                  nsig = nsig_tcdq
! YIL11: Checking only the IR value for TCT's..
               elseif(bez(myix)(1:4).eq.'TCTH' .or.                     &
     &                bez(myix)(1:4).eq.'tcth' .or.
     &                bez(myix)(1:5).eq.'TCTPH' .or.                    & 
     &                bez(myix)(1:5).eq.'tctph') then                   &
                  if(bez(myix)(8:8).eq.'1' .or.                         &                                                                                                                                       
     &                 bez(myix)(9:9).eq.'1' ) then
                     nsig = nsig_tcth1
                  elseif(bez(myix)(8:8).eq.'2' .or.                     &                                                                                                                                           
     &                 bez(myix)(9:9).eq.'2' ) then
                     nsig = nsig_tcth2
                  elseif(bez(myix)(8:8).eq.'5'.or.                      &
     &                 bez(myix)(9:9).eq.'5' ) then
                     nsig = nsig_tcth5
                  elseif(bez(myix)(8:8).eq.'8' .or.                     &                                                                                                                                           
     &                 bez(myix)(9:9).eq.'8' ) then
                     nsig = nsig_tcth8
                  endif
               elseif(bez(myix)(1:4).eq.'TCTV' .or.                     &
     &                bez(myix)(1:4).eq.'tctv'.or.
     &                bez(myix)(1:5).eq.'TCTPV' .or.                    &
     &                bez(myix)(1:5).eq.'tctpv' ) then
                  if(bez(myix)(8:8).eq.'1' .or.                         &
     &                 bez(myix)(9:9).eq.'1' ) then
                     nsig = nsig_tctv1
                  elseif(bez(myix)(8:8).eq.'2' .or.                     &
     &                 bez(myix)(9:9).eq.'2' ) then
                     nsig = nsig_tctv2
                  elseif(bez(myix)(8:8).eq.'5' .or.                     &
     &                 bez(myix)(9:9).eq.'5' ) then
                     nsig = nsig_tctv5
                  elseif(bez(myix)(8:8).eq.'8' .or.                     &
     &                 bez(myix)(9:9).eq.'8' ) then
                     nsig = nsig_tctv8
                  endif
               elseif(bez(myix)(1:3).eq.'TDI' .or.                      &
     &                 bez(myix)(1:3).eq.'tdi') then
                  nsig = nsig_tdi
               elseif(bez(myix)(1:4).eq.'TCLP' .or.                     &
     &                 bez(myix)(1:4).eq.'tclp' .or.                    &
     &                 bez(myix)(1:4).eq.'TCL.' .or.                    &
     &                 bez(myix)(1:4).eq.'tcl.'.or.                     &                                                                                                  
     &                 bez(myix)(1:4).eq.'TCLX' .or.                    &                                                                                                                   
     &                 bez(myix)(1:4).eq.'tclx') then
                  nsig = nsig_tclp
               elseif(bez(myix)(1:4).eq.'TCLI' .or.                     &
     &                 bez(myix)(1:4).eq.'tcli') then
                  nsig = nsig_tcli
               elseif(bez(myix)(1:4).eq.'TCXR' .or.                     &
     &                 bez(myix)(1:4).eq.'tcxr') then
                  nsig = nsig_tcxrp
!     TW 04/2008 ---- start adding TCRYO
               elseif(bez(myix)(1:5).eq.'TCRYO' .or.                    &
     &                 bez(myix)(1:5).eq.'tcryo'.or.                    &
     &                 bez(myix)(1:5).eq.'TCLD.' .or.                   &
     &                 bez(myix)(1:5).eq.'tcld.') then
                  nsig = nsig_tcryo
!     TW 04/2008 ---- end adding TCRYO
               elseif(bez(myix)(1:3).eq.'COL' .or.                      &
     &                 bez(myix)(1:3).eq.'col') then
                  if(bez(myix)(1:4).eq.'COLM' .or.                      &
     &                 bez(myix)(1:4).eq.'colm' .or.                    &
     &                 bez(myix)(1:5).eq.'COLH0' .or.                   &
     &                 bez(myix)(1:5).eq.'colh0') then
                     nsig = nsig_tcth1
                  elseif(bez(myix)(1:5).eq.'COLV0' .or.                 &
     &                    bez(myix)(1:5).eq.'colv0') then
                     nsig = nsig_tcth2
                  elseif(bez(myix)(1:5).eq.'COLH1' .or.                 &
     &                    bez(myix)(1:5).eq.'colh1') then
!     JUNE2005   HERE WE USE NSIG_TCTH2 AS THE OPENING IN THE VERTICAL
!     JUNE2005   PLANE FOR THE PRIMARY COLLIMATOR OF RHIC; NSIG_TCTH5 STANDS
!     JUNE2005   FOR THE OPENING OF THE FIRST SECONDARY COLLIMATOR OF RHIC
                     nsig = nsig_tcth5
                  elseif(bez(myix)(1:5).eq.'COLV1' .or.                 &
     &                    bez(myix)(1:5).eq.'colv1') then
                     nsig = nsig_tcth8
                  elseif(bez(myix)(1:5).eq.'COLH2' .or.                 &
     &                    bez(myix)(1:5).eq.'colh2') then
                     nsig = nsig_tctv1
                  endif
!     JUNE2005   END OF DEDICATED TREATMENT OF RHIC OPENINGS
               else
                  write(lout,*) "WARNING: Problem detected while "//
     &                 "writing twisslike.out' and " //
     &                 "'sigmasettings.out': Collimator name '" //
     &                 bez(myix) // "' was not recognized!"
                  write(lout,*) " ->Setting nsig=1000.0."
                  nsig = 1000.0
               endif

               do i = 1, db_ncoll
! start searching minimum gap 
                  if ((db_name1(i)(1:11).eq.bez(myix)(1:11)) .or.       &
     &                 (db_name2(i)(1:11).eq.bez(myix)(1:11))) then
                     if ( db_length(i) .gt. 0d0 ) then
                        nsig_err = nsig + gap_rms_error(i)
! jaw 1 on positive side x-axis
                        gap_h1 = nsig_err - sin(db_tilt(i,1))*          &
     &                       db_length(i)/2
                        gap_h2 = nsig_err + sin(db_tilt(i,1))*          &
     &                       db_length(i)/2
! jaw 2 on negative side of x-axis (see change of sign comapred 
! to above code lines, alos have a look to setting of tilt angle)
                        gap_h3 = nsig_err + sin(db_tilt(i,2))*          &
     &                       db_length(i)/2
                        gap_h4 = nsig_err - sin(db_tilt(i,2))*          &
     &                       db_length(i)/2
! find minumum halfgap
! --- searching for smallest halfgap
!! ---scaling for beta beat needed? 
!                        if (do_nominal) then                            
!                           bx_dist = db_bx(icoll) * scale_bx / scale_bx0
!                           by_dist = db_by(icoll) * scale_by / scale_by0
!                        else
!                           bx_dist = tbetax(j) * scale_bx / scale_bx0
!                           by_dist = tbetay(j) * scale_by / scale_by0
!                        endif
                        if (do_nominal) then                            
                           bx_dist = db_bx(icoll) 
                           by_dist = db_by(icoll)
                        else
                           bx_dist = tbetax(j)
                           by_dist = tbetay(j)
                        endif

                        sig_offset = db_offset(i) /                     &
     &                       (sqrt(bx_dist**2 * cos(db_rotation(i))**2  &
     &                       + by_dist**2 * sin(db_rotation(i))**2 ))
                        write(10000,*) bez(myix),tbetax(j),tbetay(j),   &
     &                       torbx(j),torby(j), nsig, gap_rms_error(i)
                        write(10001,*) bez(myix), gap_h1, gap_h2,       & 
     &                       gap_h3, gap_h4, sig_offset, db_offset(i),  &
     &                       nsig, gap_rms_error(i)
                        if ((gap_h1 + sig_offset) .le. mingap) then
                           mingap = gap_h1 + sig_offset
                           coll_mingap_id = i
                           coll_mingap1 = db_name1(i)
                           coll_mingap2 = db_name2(i) 
                        elseif ((gap_h2 + sig_offset) .le. mingap) then
                           mingap = gap_h2 + sig_offset
                           coll_mingap_id = i
                           coll_mingap1 = db_name1(i)
                           coll_mingap2 = db_name2(i) 
                        elseif ((gap_h3 - sig_offset) .le. mingap) then
                           mingap = gap_h3 - sig_offset
                           coll_mingap_id = i
                           coll_mingap1 = db_name1(i)
                           coll_mingap2 = db_name2(i) 
                        elseif ((gap_h4 - sig_offset) .le. mingap) then
                           mingap = gap_h4 - sig_offset
                           coll_mingap_id = i
                           coll_mingap1 = db_name1(i)
                           coll_mingap2 = db_name2(i)
                        endif
                     endif
                  endif
               enddo !do i = 1, db_ncoll

! could be done more elegant the above code to search the minimum gap
! and should also consider the jaw tilt
            endif
         enddo !do j=1,iu
         write(10000,*) coll_mingap_id,coll_mingap1,coll_mingap2,       &
     &        mingap
         write(10000,*) 'INFO> IPENCIL initial ',ipencil

! if pencil beam is used and on collimator with smallest gap the
! distribution should be generated, set ipencil to coll_mingap_id    
         if (ipencil.gt.0 .and. do_mingap) then
            ipencil = coll_mingap_id
         endif

         write(10000,*) 'INFO> IPENCIL new (if do_mingap) ',ipencil
         write(10001,*) coll_mingap_id,coll_mingap1,coll_mingap2,       &
     &        mingap

! if pencil beam is used and on collimator with smallest gap the
! distribution should be generated, set ipencil to coll_mingap_id    
         write(10001,*) 'INFO> IPENCIL new (if do_mingap) ',ipencil
         write(10001,*) 'INFO> rnd_seed is (before reinit)',rnd_seed

         close(10000)
         close(10001)

!****** re-intialize random generator with rnd_seed 
!       reinit with initial value used in first call  
         rnd_lux = 3
         rnd_k1  = 0
         rnd_k2  = 0
         call rluxgo(rnd_lux, rnd_seed, rnd_k1, rnd_k2)
! TW - 01/2007

!GRD INITIALIZE LOCAL ADDITIVE PARAMETERS, ie THE ONE WE DON'T WANT
!GRD TO KEEP OVER EACH LOOP
         do j=1,napx
           tertiary(j)=0
           secondary(j)=0
           other(j)=0
           nabs_type(j) = 0
         end do

         do k = 1, numeff
           neff(k)  = 0d0
           neffx(k) = 0d0
           neffy(k) = 0d0

             do j = 1, numeffdpop
              neff2d(k,j) = 0d0
             enddo
         enddo

         do k = 1, numeffdpop
           neffdpop(k)  = 0d0
           npartdpop(k) = 0
         enddo

         do j=1,max_ncoll
           cn_impact(j) = 0
           cn_absorbed(j) = 0
           csum(j) = 0d0
           csqsum(j) = 0d0
         enddo

!++ End of first call stuff (end of first run)
      endif

!GRD NOW WE CAN BEGIN THE LOOPS
      end

!>
!! collimate_start_collimator()
!! This routine is called each time we hit a collimator
!<
      subroutine collimate_start_collimator(stracki)
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
      integer i,ix,j,jj,jx,kpz,kzz,napx0,nbeaux,nmz,nthinerr
      double precision benkcc,cbxb,cbzb,cikveb,crkveb,crxb,crzb,r0,r000,&
     &r0a,r2b,rb,rho2b,rkb,tkb,xbb,xrb,zbb,zrb
      logical lopen
+ca parpro
+ca parnum
+ca common
+ca commons
+ca commont1
+ca commondl
+ca commonxz
+ca commonta
+ca commonmn
+ca commonm1
+ca commontr
+ca beamdim
+ca commonex
      dimension nbeaux(nbb)
+if collimat
+ca collpara
+ca dbtrthin
+ca database
+ca dbcommon
+ca dblinopt
+ca dbpencil
+ca info
+ca dbcolcom
+ca dbthin6d

+ei

+if bnlelens
+ca rhicelens
+ei

      double precision c5m4,stracki

+if fast
      c5m4=5.0d-4
+ei

            if(bez(myix)(1:3).eq.'TCP' .or.                             &
     &           bez(myix)(1:3).eq.'tcp') then

              if(bez(myix)(7:9).eq.'3.B' .or.                           &
     &             bez(myix)(7:9).eq.'3.b') then
                nsig = nsig_tcp3
              else
                nsig = nsig_tcp7
              endif

            elseif(bez(myix)(1:4).eq.'TCSG' .or.                        &
     &             bez(myix)(1:4).eq.'tcsg') then

              if(bez(myix)(8:10).eq.'3.B' .or.                          &
     &             bez(myix)(8:10).eq.'3.b' .or.                        &
     &             bez(myix)(9:11).eq.'3.B' .or.                        &
     &             bez(myix)(9:11).eq.'3.b') then
                nsig = nsig_tcsg3
              else
                nsig = nsig_tcsg7
              endif

              if((bez(myix)(5:6).eq.'.4'.and.bez(myix)(8:9).eq.'6.')    &
     &             ) then
                nsig = nsig_tcstcdq
              endif

            elseif(bez(myix)(1:4).eq.'TCSP' .or.                        &
     &             bez(myix)(1:4).eq.'tcsp') then

             if(bez(myix)(9:11).eq.'6.B'.or.
     &          bez(myix)(9:11).eq.'6.b') then
                nsig = nsig_tcstcdq
              endif

            elseif(bez(myix)(1:4).eq.'TCSM' .or.                        &
     &             bez(myix)(1:4).eq.'tcsm') then

              if(bez(myix)(8:10).eq.'3.B' .or.                          &
     &             bez(myix)(8:10).eq.'3.b' .or.                        &
     &             bez(myix)(9:11).eq.'3.B' .or.                        &
     &             bez(myix)(9:11).eq.'3.b') then
                nsig = nsig_tcsm3
              else
                nsig = nsig_tcsm7
              endif

            elseif(bez(myix)(1:4).eq.'TCLA' .or.                        &
     &             bez(myix)(1:4).eq.'tcla') then

              if(bez(myix)(9:11).eq.'7.B' .or.                          &
     &             bez(myix)(9:11).eq.'7.b') then
                nsig = nsig_tcla7
              else
                nsig = nsig_tcla3
              endif

            elseif(bez(myix)(1:4).eq.'TCDQ' .or.                        &
     &             bez(myix)(1:4).eq.'tcdq') then
              nsig = nsig_tcdq
! YIL11: Checking only the IR value for TCT's..
            elseif(bez(myix)(1:4).eq.'TCTH' .or.                        &
     &                bez(myix)(1:4).eq.'tcth' .or.
     &                bez(myix)(1:5).eq.'TCTPH' .or.                    & 
     &                bez(myix)(1:5).eq.'tctph') then                   &

                  if(bez(myix)(8:8).eq.'1' .or.                         &
     &                 bez(myix)(9:9).eq.'1' ) then
                     nsig = nsig_tcth1
                  elseif(bez(myix)(8:8).eq.'2' .or.                     &
     &                 bez(myix)(9:9).eq.'2' ) then
                     nsig = nsig_tcth2
                  elseif(bez(myix)(8:8).eq.'5'.or.                      &
     &                 bez(myix)(9:9).eq.'5' ) then
                     nsig = nsig_tcth5
                  elseif(bez(myix)(8:8).eq.'8' .or.                     &
     &                 bez(myix)(9:9).eq.'8' ) then
                     nsig = nsig_tcth8
                  endif

            elseif(bez(myix)(1:4).eq.'TCTV' .or.                        &
     &                bez(myix)(1:4).eq.'tctv'.or.
     &                bez(myix)(1:5).eq.'TCTPV' .or.                    &
     &                bez(myix)(1:5).eq.'tctpv' ) then

                  if(bez(myix)(8:8).eq.'1' .or.                         &
     &                 bez(myix)(9:9).eq.'1' ) then
                     nsig = nsig_tctv1
                  elseif(bez(myix)(8:8).eq.'2' .or.                     &
     &                 bez(myix)(9:9).eq.'2' ) then
                     nsig = nsig_tctv2
                  elseif(bez(myix)(8:8).eq.'5' .or.                     &
     &                 bez(myix)(9:9).eq.'5' ) then
                     nsig = nsig_tctv5
                  elseif(bez(myix)(8:8).eq.'8' .or.                     &
     &                 bez(myix)(9:9).eq.'8' ) then
                     nsig = nsig_tctv8
                  endif

            elseif(bez(myix)(1:3).eq.'TDI' .or.                         &
     &             bez(myix)(1:3).eq.'tdi') then
              nsig = nsig_tdi
            elseif(bez(myix)(1:4).eq.'TCLP' .or.                        &
     &             bez(myix)(1:4).eq.'tclp' .or.                        &
     &             bez(myix)(1:4).eq.'TCL.' .or.                        &
     &             bez(myix)(1:4).eq.'tcl.'.or.                         &
     &             bez(myix)(1:4).eq.'TCLX' .or.                        &
     &             bez(myix)(1:4).eq.'tclx') then
              nsig = nsig_tclp
            elseif(bez(myix)(1:4).eq.'TCLI' .or.                        &
     &             bez(myix)(1:4).eq.'tcli') then
              nsig = nsig_tcli
            elseif(bez(myix)(1:4).eq.'TCXR' .or.                        &
     &             bez(myix)(1:4).eq.'tcxr') then
              nsig = nsig_tcxrp
            elseif(bez(myix)(1:5).eq.'TCRYO' .or.                       &
     &             bez(myix)(1:5).eq.'tcryo'.or.
     &             bez(myix)(1:5).eq.'TCLD.' .or.                       &
     &             bez(myix)(1:5).eq.'tcld.') then
              nsig = nsig_tcryo
            elseif(bez(myix)(1:3).eq.'COL' .or.                         &
     &             bez(myix)(1:3).eq.'col') then

              if(bez(myix)(1:4).eq.'COLM' .or.                          &
     &             bez(myix)(1:4).eq.'colm' .or.                        &
     &             bez(myix)(1:5).eq.'COLH0' .or.                       &
     &             bez(myix)(1:5).eq.'colh0') then
                nsig = nsig_tcth1
              elseif(bez(myix)(1:5).eq.'COLV0' .or.                     &
     &               bez(myix)(1:5).eq.'colv0') then
                nsig = nsig_tcth2
              elseif(bez(myix)(1:5).eq.'COLH1' .or.                     &
     &               bez(myix)(1:5).eq.'colh1') then
!     JUNE2005   HERE WE USE NSIG_TCTH2 AS THE OPENING IN THE VERTICAL
!     JUNE2005   PLANE FOR THE PRIMARY COLLIMATOR OF RHIC; NSIG_TCTH5 STANDS
!     JUNE2005   FOR THE OPENING OF THE FIRST SECONDARY COLLIMATOR OF RHIC
                nsig = nsig_tcth5
              elseif(bez(myix)(1:5).eq.'COLV1' .or.                     &
     &               bez(myix)(1:5).eq.'colv1') then
                nsig = nsig_tcth8
              elseif(bez(myix)(1:5).eq.'COLH2' .or.                     &
     &               bez(myix)(1:5).eq.'colh2') then
                nsig = nsig_tctv1
              endif

            else
              if(firstrun.and.iturn.eq.1) then
                 write(lout,*) "WARNING: When setting opening for the"//
     &                " collimator named '" // bez(myix) //
     &                "' from fort.3, the name was not recognized."
                 write(lout,*) " -> Setting nsig=1000.0."
              endif
              nsig=1000.0
!JUNE2005   END OF DEDICATED TREATMENT OF RHIC OPENINGS
            endif

!++  Write trajectory for any selected particle
        c_length = 0d0

!     SR, 23-11-2005: To avoid binary entries in 'amplitude.dat'
        if ( firstrun ) then
          if (rselect.gt.0 .and. rselect.lt.65) then
            do j = 1, napx
              xj     = (xv(1,j)-torbx(ie))/1d3
              xpj    = (yv(1,j)-torbxp(ie))/1d3
              yj     = (xv(2,j)-torby(ie))/1d3
              ypj    = (yv(2,j)-torbyp(ie))/1d3
              pj     = ejv(j)/1d3

              if (iturn.eq.1.and.j.eq.1) then
                sum_ax(ie)=0d0
                sum_ay(ie)=0d0
              endif

!-- DRIFT PART
              if (stracki.eq.0.) then
                if(iexact.eq.0) then
                  xj  = xj + 0.5d0*c_length*xpj
                  yj  = yj + 0.5d0*c_length*ypj
                else
                  zpj=sqrt(1d0-xpj**2-ypj**2)
                  xj = xj + 0.5d0*c_length*(xpj/zpj)
                  yj = yj + 0.5d0*c_length*(ypj/zpj)
                endif
              endif
!
              gammax = (1d0 + talphax(ie)**2)/tbetax(ie)
              gammay = (1d0 + talphay(ie)**2)/tbetay(ie)
!
              if (part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
          nspx    = sqrt(                                               &
     &abs( gammax*(xj)**2 +                                             &
     &2d0*talphax(ie)*xj*xpj +                                          &
     &tbetax(ie)*xpj**2 )/myemitx0_collgap
     &)
                nspy    = sqrt(                                         &
     &abs( gammay*(yj)**2 +                                             &
     &2d0*talphay(ie)*yj*ypj +                                          &
     &tbetay(ie)*ypj**2 )/myemity0_collgap
     &)
                sum_ax(ie)   = sum_ax(ie) + nspx
                sqsum_ax(ie) = sqsum_ax(ie) + nspx**2
                sum_ay(ie)   = sum_ay(ie) + nspy
                sqsum_ay(ie) = sqsum_ay(ie) + nspy**2
                nampl(ie)    = nampl(ie) + 1
              else
                nspx = 0d0
                nspy = 0d0
              endif
                sampl(ie)    = totals
                ename(ie)    = bez(myix)(1:16)
            end do
          endif
      endif

!GRD HERE WE LOOK FOR ADEQUATE DATABASE INFORMATION
          found = .false.
!     SR, 01-09-2005: to set found = .TRUE., add the condition L>0!!
          do j = 1, db_ncoll
            if ((db_name1(j)(1:11).eq.bez(myix)(1:11)) .or.             &
     &          (db_name2(j)(1:11).eq.bez(myix)(1:11))) then
               if ( db_length(j) .gt. 0d0 ) then
                 found = .true.
                 icoll = j
               endif
            endif
          end do
          if (.not. found .and. firstrun .and. iturn.eq.1) then
            write(lout,*)
     &           'ERR>  Collimator not found in colldb: ', bez(myix)
      endif

      end

!>
!! collimate_do_collimator()
!! This routine is calls the actual scattering functions
!<
      subroutine collimate_do_collimator(stracki)
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
      integer i,ix,j,jj,jx,kpz,kzz,napx0,nbeaux,nmz,nthinerr
      double precision benkcc,cbxb,cbzb,cikveb,crkveb,crxb,crzb,r0,r000,&
     &r0a,r2b,rb,rho2b,rkb,tkb,xbb,xrb,zbb,zrb
      logical lopen
+ca parpro
+ca parnum
+ca common
+ca commons
+ca commont1
+ca commondl
+ca commonxz
+ca commonta
+ca commonmn
+ca commonm1
+ca commontr
+ca beamdim
+ca commonex
      dimension nbeaux(nbb)
+if collimat
+ca collpara
+ca dbtrthin
+ca database
+ca dbcommon
+ca dblinopt
+ca dbpencil
+ca info
+ca dbcolcom
+ca dbthin6d

+ei

+if bnlelens
+ca rhicelens
+ei

      double precision c5m4,stracki

+if g4collimat
      integer g4_lostc
      double precision x_tmp,y_tmp,xp_tmp,yp_tmp
+ei

+if fast
      c5m4=5.0d-4
+ei

!-----------------------------------------------------------------------
!GRD NEW COLLIMATION PARAMETERS
!-----------------------------------------------------------------------
!++  Get the aperture from the beta functions and emittance
!++  A simple estimate of beta beating can be included that
!++  has twice the betatron phase advance
         if(.not. do_nsig) nsig = db_nsig(icoll)
+if crlibm
          scale_bx = (1d0 + xbeat*sin_rn(4*pi*mux(ie)+
+ei
+if .not.crlibm
          scale_bx = (1d0 + xbeat*sin(4*pi*mux(ie)+                     &
+ei
     &xbeatphase)  )
+if crlibm
          scale_by = (1d0 + ybeat*sin_rn(4*pi*muy(ie)+
+ei
+if .not.crlibm
          scale_by = (1d0 + ybeat*sin(4*pi*muy(ie)+                     &
+ei
     &ybeatphase)  )
!
          if (firstcoll) then
            scale_bx0 = scale_bx
            scale_by0 = scale_by
            firstcoll = .false.
          endif
!-------------------------------------------------------------------
!++  Assign nominal OR design beta functions for later
          if (do_nominal) then
            bx_dist = db_bx(icoll) * scale_bx / scale_bx0
            by_dist = db_by(icoll) * scale_by / scale_by0
          else
            bx_dist = tbetax(ie) * scale_bx / scale_bx0
            by_dist = tbetay(ie) * scale_by / scale_by0
          endif

!++  Write beam ellipse at selected collimator
! ---- changed name_sel(1:11) name_sel(1:12) to be checked if feasible!!
          if (                                                          &
     &         ((db_name1(icoll).eq.name_sel(1:12))                     &
     &         .or.(db_name2(icoll).eq.name_sel(1:12)))                 &
     &         .and. dowrite_dist) then
!          if (firstrun .and.                                            &
!     &         ((db_name1(icoll).eq.name_sel(1:11))                     &
!     &         .or.(db_name2(icoll).eq.name_sel(1:11)))                 &
!     &         .and. dowrite_dist) then
! --- get halo on each turn
!     &.and. iturn.eq.1 .and. dowrite_dist) then
! --- put open and close at the pso. where it is done for the 
! --- other files belonging to dowrite_impact flag !(may not a good loc.)
!            open(unit=45, file='coll_ellipse.dat')
!            write(45,'(a)')                                             &
!     &'#  1=x 2=y 3=xp 4=yp 5=E 6=s'
            do j = 1, napx
            write(45,'(1X,I8,6(1X,E15.7),3(1X,I4,1X,I4))')              &
     &ipart(j)+100*samplenumber,xv(1,j), xv(2,j), yv(1,j), yv(2,j),     &
     &ejv(j), mys(j),iturn,secondary(j)+tertiary(j)+other(j),           &
     &nabs_type(j)
            end do
!            close(45)
          endif

!-------------------------------------------------------------------
!++  Output to temporary database and screen
          if (iturn.eq.1.and.firstrun) then
            write(40,*) '# '
            write(40,*) db_name1(icoll)(1:11)
            write(40,*) db_material(icoll)
            write(40,*) db_length(icoll)
            write(40,*) db_rotation(icoll)
            write(40,*) db_offset(icoll)
            write(40,*) tbetax(ie)
            write(40,*) tbetay(ie)
!
            write(outlun,*) ' '
            write(outlun,*)   'Collimator information: '
            write(outlun,*) ' '
            write(outlun,*) 'Name:                '                     &
     &, db_name1(icoll)(1:11)
            write(outlun,*) 'Material:            '                     &
     &, db_material(icoll)
            write(outlun,*) 'Length [m]:          '                     &
     &, db_length(icoll)
            write(outlun,*) 'Rotation [rad]:      '                     &
     &, db_rotation(icoll)
            write(outlun,*) 'Offset [m]:          '                     &
     &,db_offset(icoll)
            write(outlun,*) 'Design beta x [m]:   '                     &
     &,db_bx(icoll)
            write(outlun,*) 'Design beta y [m]:   '                     &
     &,db_by(icoll)
            write(outlun,*) 'Optics beta x [m]:   '                     &
     &,tbetax(ie)
            write(outlun,*) 'Optics beta y [m]:   '                     &
     &,tbetay(ie)
!          else
!            write(lout,*) 'C WRITE', iturn,firstrun
          endif


!-------------------------------------------------------------------
!++  Calculate aperture of collimator
!JUNE2005   HERE ONE HAS TO HAVE PARTICULAR TREATMENT OF THE OPENING OF
!JUNE2005   THE PRIMARY COLLIMATOR OF RHIC
         if(db_name1(icoll)(1:4).ne.'COLM') then
          nsig = nsig + gap_rms_error(icoll)
          xmax = nsig*sqrt(bx_dist*myemitx0_collgap)
          ymax = nsig*sqrt(by_dist*myemity0_collgap)
          xmax_pencil = (nsig+pencil_offset)*                           &
     &sqrt(bx_dist*myemitx0_collgap)
          ymax_pencil = (nsig+pencil_offset)*                           &
     &sqrt(by_dist*myemity0_collgap)
          xmax_nom = db_nsig(icoll)*sqrt(db_bx(icoll)*myemitx0_collgap)
          ymax_nom = db_nsig(icoll)*sqrt(db_by(icoll)*myemity0_collgap)
          c_rotation = db_rotation(icoll)
          c_length   = db_length(icoll)
          c_material = db_material(icoll)
          c_offset   = db_offset(icoll)
          c_tilt(1)  = db_tilt(icoll,1)
          c_tilt(2)  = db_tilt(icoll,2)

+if crlibm
          calc_aperture = sqrt( xmax**2 * cos_rn(c_rotation)**2         &
+ei
+if .not.crlibm
          calc_aperture = sqrt( xmax**2 * cos(c_rotation)**2            &
+ei
+if crlibm
     &                    + ymax**2 * sin_rn(c_rotation)**2 )
+ei
+if .not.crlibm
     &                    + ymax**2 * sin(c_rotation)**2 )
+ei

+if crlibm
          nom_aperture = sqrt( xmax_nom**2 * cos_rn(c_rotation)**2      &
     &                   + ymax_nom**2 * sin_rn(c_rotation)**2 )
+ei
+if .not.crlibm
          nom_aperture = sqrt( xmax_nom**2 * cos(c_rotation)**2         &
     &                   + ymax_nom**2 * sin(c_rotation)**2 )
+ei
!
            pencil_aperture =                                           &
+if crlibm
     &                    sqrt( xmax_pencil**2 * cos_rn(c_rotation)**2  &
     &                    + ymax_pencil**2 * sin_rn(c_rotation)**2 )
+ei
+if .not.crlibm
     &                    sqrt( xmax_pencil**2 * cos(c_rotation)**2     &
     &                    + ymax_pencil**2 * sin(c_rotation)**2 )
+ei

!++  Get x and y offsets at collimator center point
+if crlibm
            x_pencil(icoll) = xmax_pencil * (cos_rn(c_rotation))
            y_pencil(icoll) = ymax_pencil * (sin_rn(c_rotation))
+ei
+if .not.crlibm
            x_pencil(icoll) = xmax_pencil * (cos(c_rotation))
            y_pencil(icoll) = ymax_pencil * (sin(c_rotation))
+ei

!++  Get corresponding beam angles (uses xp_max)
          xp_pencil(icoll) =                                            &
     &              -1d0 * sqrt(myemitx0_collgap/tbetax(ie))*talphax(ie)
     &                   * xmax / sqrt(myemitx0_collgap*tbetax(ie))
     
          yp_pencil(icoll) =                                            &
     &              -1d0 * sqrt(myemity0_collgap/tbetay(ie))*talphay(ie)
     &                   * ymax / sqrt(myemity0_collgap*tbetay(ie))

! that the way xp is calculated for makedis subroutines !!!!
!        if (rndm4().gt.0.5) then
!          myxp(j)  = sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-        &
!     &myalphax*myx(j)/mybetax
!        else
!          myxp(j)  = -1*sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-     &
!     &myalphax*myx(j)/mybetax
!        endif
!            xp_pencil(icoll) =                                          &
!     &           sqrt(sqrt((myemitx0/tbetax(ie)                         &
!     &           -x_pencil(icoll)**2/tbetax(ie)**2)**2))                &
!     &           -talphax(ie)*x_pencil(icoll)/tbetax(ie)
!            write(*,*) " ************************************ "
!            write(*,*) myemitx0/tbetax(ie)                              &
!     &           -x_pencil(icoll)**2/tbetax(ie)**2
!            write(*,*)sqrt(sqrt((myemitx0/tbetax(ie)                    &
!     &           -x_pencil(icoll)**2/tbetax(ie)**2)**2))
!            write(*,*) -talphax(ie)*x_pencil(icoll)/tbetax(ie)
!            write(*,*) sqrt(myemitx0/tbetax(ie))*talphax(ie)            &
!     &                   * x_pencil(icoll) / sqrt(myemitx0*tbetax(ie))
!            write(*,*)  sqrt(sqrt((myemitx0/tbetax(ie)                  &
!     &           -x_pencil(icoll)**2/tbetax(ie)**2)**2))                &
!     &           -talphax(ie)*x_pencil(icoll)/tbetax(ie)
!            write(*,*) xp_pencil(icoll)
!            write(*,*) " ************************************ "
!
!            yp_pencil(icoll) =                                          &
!     &           sqrt(sqrt((myemity0/tbetay(ie)                         &
!     &           -y_pencil(icoll)**2/tbetay(ie)**2)**2))                &
!     &           -talphay(ie)*y_pencil(icoll)/tbetay(ie)
!!
            xp_pencil0 = xp_pencil(icoll)
            yp_pencil0 = yp_pencil(icoll)

            pencil_dx(icoll)  =                                         &
+if crlibm
     &                     sqrt( xmax_pencil**2 * cos_rn(c_rotation)**2 &
     &                     + ymax_pencil**2 * sin_rn(c_rotation)**2 )   &
+ei
+if .not.crlibm
     &                     sqrt( xmax_pencil**2 * cos(c_rotation)**2    &
     &                     + ymax_pencil**2 * sin(c_rotation)**2 )      &
+ei
     &                     - calc_aperture
!++ TW -- tilt for of jaw for pencil beam
!++ as in Ralphs orig routine, but not in collimate subroutine itself
!            nprim = 3
!            if ( (icoll.eq.ipencil) &
!     &           icoll.le.nprim .and. (j.ge.(icoll-1)*nev/nprim)        &
!     &           .and. (j.le.(icoll)*nev/nprim))) then
! this is done for every bunch (64 particle bucket)
! important: Sixtrack calculates in "mm" and collimate2 in "m"
! therefore 1E-3 is used to  
            if ((icoll.eq.ipencil).and.(iturn.eq.1).and.
     &           (pencil_distr.ne.3)) then ! RB: added condition that pencil_distr.ne.3 in order to do the tilt

!!               write(*,*) " ************************************** "
!!               write(*,*) " * INFO> seting tilt for pencil beam  * "
!!               write(*,*) " ************************************** "
!     c_tilt(1) =  (xp_pencil0*cos(c_rotation)                  &
+if crlibm
!adriana
               c_tilt(1) = c_tilt(1) + (xp_pencil0*cos_rn(c_rotation)    &
     &                     + sin_rn(c_rotation)*yp_pencil0)
+ei
+if .not.crlibm
               c_tilt(1) = c_tilt(1) + (xp_pencil0*cos(c_rotation)       &
     &                     + sin(c_rotation)*yp_pencil0)
+ei
               write(lout,*) "INFO> Changed tilt1  ICOLL  to  ANGLE  ",  &
     &              icoll, c_tilt(1)
!
!! respects if the tilt symmetric or not, for systilt_antiymm c_tilt is 
!! -systilt + rmstilt otherwise +systilt + rmstilt
!!               if (systilt_antisymm) then
!! to align the jaw/pencil to the beam always use the minus regardless which 
!! orientation of the jaws was used (symmetric/antisymmetric) 
!                c_tilt(2) =  -1.*(xp_pencil0*cos(c_rotation)             &
+if crlibm
!adriana
                c_tilt(2) = c_tilt(2) -1.*(xp_pencil0*cos_rn(c_rotation)  &
     &                 + sin_rn(c_rotation)*yp_pencil0)
+ei
+if .not.crlibm
                c_tilt(2) = c_tilt(2) -1.*(xp_pencil0*cos(c_rotation)     &
     &                 + sin(c_rotation)*yp_pencil0)
+ei
!!               else
!!                  c_tilt(2) = c_tilt(2) + (xp_pencil0*cos(c_rotation)   &
!!     &                 + sin(c_rotation)*yp_pencil0)
!!               endif
               write(lout,*) "INFO> Changed tilt2  ICOLL  to  ANGLE  ",   &
     &              icoll, c_tilt(2)
            endif

!++ TW -- tilt angle changed (added to genetated on if spec. in fort.3) 

!JUNE2005   HERE IS THE SPECIAL TREATMENT...
         elseif(db_name1(icoll)(1:4).eq.'COLM') then

            xmax = nsig_tcth1*sqrt(bx_dist*myemitx0_collgap)
            ymax = nsig_tcth2*sqrt(by_dist*myemity0_collgap)

            c_rotation = db_rotation(icoll)
            c_length   = db_length(icoll)
            c_material = db_material(icoll)
            c_offset   = db_offset(icoll)
            c_tilt(1)  = db_tilt(icoll,1)
            c_tilt(2)  = db_tilt(icoll,2)
            calc_aperture = xmax
            nom_aperture = ymax
         endif





!-------------------------------------------------------------------
!++  Further output
        if(firstrun) then
          if (iturn.eq.1) then
            write(outlun,*) xp_pencil(icoll), yp_pencil(icoll),         &
     &pencil_dx(icoll)
            write(outlun,'(a,i4)') 'Collimator number:   '              &
     &,icoll
            write(outlun,*) 'Beam size x [m]:     '                     &
     &,sqrt(tbetax(ie)*myemitx0_collgap), "(from collgap emittance)"
            write(outlun,*) 'Beam size y [m]:     '                     &
     &,sqrt(tbetay(ie)*myemity0_collgap), "(from collgap emittance)"
            write(outlun,*) 'Divergence x [urad]:     '                 &
     &,1d6*xp_pencil(icoll)
            write(outlun,*) 'Divergence y [urad]:     '                 &
     &,1d6*yp_pencil(icoll)
            write(outlun,*) 'Aperture (nom) [m]:  '                     &
     &,nom_aperture
            write(outlun,*) 'Aperture (cal) [m]:  '                     &
     &,calc_aperture
            write(outlun,*) 'Collimator halfgap [sigma]:  '             &
     &,nsig
            write(outlun,*) 'RMS error on halfgap [sigma]:  '           &
     &,gap_rms_error(icoll)
            write(outlun,*) ' '

            write(43,'(i10,1x,a,4(1x,e19.10),1x,a,6(1x,e13.5))')
     &icoll,db_name1(icoll)(1:12),                                      &
     &db_rotation(icoll),                                               &
     &tbetax(ie), tbetay(ie), calc_aperture,                            &
     &db_material(icoll),                                               &
     &db_length(icoll),                                                 &
     &sqrt(tbetax(ie)*myemitx0_collgap),                                &
     &sqrt(tbetay(ie)*myemity0_collgap),                                &
     &db_tilt(icoll,1),                                                 &
     &db_tilt(icoll,2),                                                 &
     &nsig

! coll settings file
            if(n_slices.le.1) then
            write(55,'(a,1x,i10,5(1x,e13.5),1x,a)')                     &
     &db_name1(icoll)(1:12),                                            &
     &n_slices,calc_aperture,                                           &
     &db_offset(icoll),                                                 &
     &db_tilt(icoll,1),                                                 &
     &db_tilt(icoll,2),                                                 &
     &db_length(icoll),                                                 &
     &db_material(icoll)
         endif
          endif
        endif

!++  Assign aperture which we define as the FULL width (factor 2)!!!
!JUNE2005 AGAIN, SOME SPECIFIC STUFF FOR RHIC
         if(db_name1(icoll)(1:4).eq.'COLM') then
            c_aperture = 2d0*calc_aperture
            nom_aperture = 2d0*nom_aperture
         elseif(db_name1(icoll)(1:4).ne.'COLM') then
            c_aperture = 2d0*calc_aperture
         endif
!JUNE2005
          c_aperture = 2d0*calc_aperture
!          IF(IPENCIL.GT.zero) THEN
!          C_APERTURE = 2.*pencil_aperture
!GRD-------------------------------------------------------------------
      if(firstrun.and.iturn.eq.1.and.icoll.eq.7) then
        open(unit=99,file='distsec')
        do j=1,napx
          write(99,'(4(1X,E15.7))') xv(1,j),yv(1,j),xv(2,j),yv(2,j)
        enddo
        close(99)
      endif
!GRD-------------------------------------------------------------------


! RB: addition matched halo sampled directly on the TCP using pencil beam flag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if ((iturn.eq.1).and.(ipencil.eq.icoll).and.
     &         (pencil_distr.eq.3)) then

!     create distribution where the normalized distance between jaw and beam is the smallest - this is where particles will first impact:
!     without imperfections, it is:
!              -- at the face of the collimator for the case of beta'<0 (POSITIVE alpha - beam converging) and 
!              -- at the exit of the collimator for the case of beta'>0 (NEGATIVE alpha beam diverging)

!     with imperfections: include errors on gap, tilt and offset. We have to calculate the normalized distance to each corner separately!

!     First: calculate optical parameters at start and end of collimator (half a collimator length upstream and downstream of present s-position)
!     Assuming a purely vertical or horizontal halo - need to add more conditions for other cases!
             
!     Using standard twiss transfer matrix for a drift : ( new_halo_model_checks.nb )
!     at start of collimator:
             ldrift = -c_length / 2.d0 !Assign the drift length over which the optics functions are propagated
             betax1 = tbetax(ie) - 2*ldrift*talphax(ie) + 
     &            (ldrift**2 * (1+talphax(ie)**2))/tbetax(ie) 
             betay1 = tbetay(ie) - 2*ldrift*talphay(ie) + 
     &            (ldrift**2 * (1+talphay(ie)**2))/tbetay(ie)

             alphax1 = talphax(ie) - 
     &            (ldrift*(1+talphax(ie)**2))/tbetax(ie)
             alphay1 = talphay(ie) - 
     &            (ldrift*(1+talphay(ie)**2))/tbetay(ie)

!     at end of collimator:
             ldrift = c_length / 2.d0
             betax2 = tbetax(ie) - 2*ldrift*talphax(ie) + 
     &            (ldrift**2 * (1+talphax(ie)**2))/tbetax(ie) 
             betay2 = tbetay(ie) - 2*ldrift*talphay(ie) + 
     &            (ldrift**2 * (1+talphay(ie)**2))/tbetay(ie)

             alphax2 = talphax(ie) - 
     &            (ldrift*(1+talphax(ie)**2))/tbetax(ie)
             alphay2 = talphay(ie) - 
     &            (ldrift*(1+talphay(ie)**2))/tbetay(ie)

!     calculate beam size at start and end of collimator. account for collimation plane
             if((mynex.gt.0).and.(myney.eq.0.0)) then  ! horizontal halo 
                beamsize1 = sqrt(betax1 * myemitx0_collgap)
                beamsize2 = sqrt(betax2 * myemitx0_collgap)
             elseif((mynex.eq.0).and.(myney.gt.0.0)) then   ! vertical halo
                beamsize1 = sqrt(betay1 * myemity0_collgap)
                beamsize2 = sqrt(betay2 * myemity0_collgap)
             else
                write(lout,*)
     &               "attempting to use a halo not purely in the "//
     &               "horizontal or vertical plane with pencil_dist=3"//
     &               " - abort."
               call prror(-1)
             endif
             
!     calculate offset from tilt of positive and negative jaws, at start and end
!     remember: tilt angle is defined such that one corner stays at nominal position, the other corner is more open

!     jaw in positive x (or y):
             if (c_tilt(1).ge.0) then
                tiltOffsPos1 = 0.d0
                tiltOffsPos2 = abs(sin(c_tilt(1))) * c_length
             else
                tiltOffsPos1 = abs(sin(c_tilt(1))) * c_length
                tiltOffsPos2 = 0.d0
             endif

!     jaw in negative x (or y):
             if (c_tilt(2).ge.0) then
                tiltOffsNeg1 = abs(sin(c_tilt(2))) * c_length
                tiltOffsNeg2 = 0.d0
             else
                tiltOffsNeg1 = 0.d0
                tiltOffsNeg2 = abs(sin(c_tilt(2))) * c_length
             endif

!     calculate half distance from jaws to beam center (in units of beam sigma) at the beginning of the collimator, positive and neg jaws. 
            Nap1pos=(c_aperture/2d0 + c_offset + tiltOffsPos1)/beamsize1
            Nap2pos=(c_aperture/2d0 + c_offset + tiltOffsPos2)/beamsize2
            Nap1neg=(c_aperture/2d0 - c_offset + tiltOffsNeg1)/beamsize1
            Nap2neg=(c_aperture/2d0 - c_offset + tiltOffsNeg2)/beamsize2

! debugging output - can be removed when not needed
!            write(7878,*) c_tilt(1),c_tilt(2),c_offset
!       write(7878,*) tiltOffsPos1,tiltOffsPos2,tiltOffsNeg1,tiltOffsNeg2
!            write(7878,*) Nap1pos,Nap2pos,Nap1neg,Nap2neg
!            write(7878,*) min(Nap1pos,Nap2pos,Nap1neg,Nap2neg)
!            write(7878,*) mynex * sqrt(tbetax(ie)/betax1)

!     Minimum normalized distance from jaw to beam center - this is the n_sigma at which the halo should be generated
            minAmpl = min(Nap1pos,Nap2pos,Nap1neg,Nap2neg) 

!     Assign amplitudes in x and y for the halo generation function
            if((mynex.gt.0).and.(myney.eq.0.0)) then ! horizontal halo 
               mynex2 = minAmpl 
            elseif((mynex.eq.0).and.(myney.gt.0.0)) then ! vertical halo
               myney2 = minAmpl
            endif               ! other cases taken care of above - in these cases, program has already stopped            

!     assign optics parameters to use for the generation of the starting halo - at start or end of collimator
             if((minAmpl.eq.Nap1pos).or.(minAmpl.eq.Nap1neg)) then ! min normalized distance occurs at start of collimator
                mybetax=betax1
                mybetay=betay1
                myalphax=alphax1
                myalphay=alphay1
                ldrift = -c_length / 2.d0
             else               ! min normalized distance occurs at end of collimator
                mybetax=betax2
                mybetay=betay2
                myalphax=alphax2
                myalphay=alphay2
                ldrift = c_length / 2.d0
             endif

             write(7878,*) napx,myalphax,myalphay, mybetax, mybetay,
     &            myemitx0_collgap, myemity0_collgap,
     &            myenom, mynex2, mdex, myney2,mdey

!     create new pencil beam distribution with spread at start or end of collimator at the minAmpl
!     note: if imperfections are active, equal amounts of particles are still generated on the two jaws.
!     but it might be then that only one jaw is hit on the first turn, thus only by half of the particles
!     the particle generated on the other side will then hit the same jaw several turns later, possibly smearing the impact parameter
!     This could possibly be improved in the future.
             call makedis_coll(napx,myalphax,myalphay, mybetax, mybetay,
     &            myemitx0_collgap, myemity0_collgap,
     &            myenom, mynex2, mdex, myney2,mdey,
     &            myx, myxp, myy, myyp, myp, mys)
             
             do j = 1, napx
                xv(1,j)  = 1d3*myx(j)  + torbx(ie) 
                yv(1,j)  = 1d3*myxp(j) + torbxp(ie)             
                xv(2,j)  = 1d3*myy(j)  + torby(ie)              
                yv(2,j)  = 1d3*myyp(j) + torbyp(ie)             
                sigmv(j) = mys(j)
                ejv(j)   = myp(j)

!     as main routine will track particles back half a collimator length (to start of jaw), 
!     track them now forward (if generated at face) or backward (if generated at end) 
!     1/2 collimator length to center of collimator (ldrift pos or neg)
                xv(1,j)  = xv(1,j) - ldrift*yv(1,j)
                xv(2,j)  = xv(2,j) - ldrift*yv(2,j)

!     write out distribution - generated either at the BEGINNING or END of the collimator
                write(4997,'(6(1X,E15.7))') myx(j), myxp(j), myy(j), 
     &               myyp(j), mys(j), myp(j)
             enddo
          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! end RB addition

!++  Copy particle data to 1-dim array and go back to meters
            do j = 1, napx
              rcx(j)  = (xv(1,j)-torbx(ie))/1d3
              rcxp(j) = (yv(1,j)-torbxp(ie))/1d3
              rcy(j)  = (xv(2,j)-torby(ie))/1d3
              rcyp(j) = (yv(2,j)-torbyp(ie))/1d3
              rcp(j)  = ejv(j)/1d3
              rcs(j)  = 0d0
              part_hit_before_turn(j) = part_hit_turn(j)
              part_hit_before_pos(j)  = part_hit_pos(j)
              rcx0(j)  = rcx(j)
              rcxp0(j) = rcxp(j)
              rcy0(j)  = rcy(j)
              rcyp0(j) = rcyp(j)
              rcp0(j)  = rcp(j)
              ejf0v(j) = ejfv(j)

!++  For zero length element track back half collimator length
!  DRIFT PART
              if (stracki.eq.0.) then
                if(iexact.eq.0) then
                  rcx(j)  = rcx(j) - 0.5d0*c_length*rcxp(j)
                  rcy(j)  = rcy(j) - 0.5d0*c_length*rcyp(j)
                else
                  zpj=sqrt(1d0-rcxp(j)**2-rcyp(j)**2)
                  rcx(j) = rcx(j) - 0.5d0*c_length*(rcxp(j)/zpj)
                  rcy(j) = rcy(j) - 0.5d0*c_length*(rcyp(j)/zpj)
                endif
              else
                write(lout,*) "ERROR: Non-zero length collimator!"
                call prror(-1)
              endif

              flukaname(j) = ipart(j)+100*samplenumber
            end do

!++  Do the collimation tracking
             enom_gev = myenom*1d-3

!++  Allow primaries to be one-sided, if requested
          if ((db_name1(icoll)(1:3).eq.'TCP' .or.                       &
     &db_name1(icoll)(1:3).eq.'COL')                                    &
     &.and. do_oneside) then
            onesided = .true.
          else
            onesided = .false.
          endif


!GRD HERE IS THE MAJOR CHANGE TO THE CODE: IN ORDER TO TRACK PROPERLY THE
!GRD SPECIAL RHIC PRIMARY COLLIMATOR, IMPLEMENTATION OF A DEDICATED ROUTINE
          if (found) then
            if(db_name1(icoll)(1:4).eq.'COLM') then
               call collimaterhic(c_material,                           &
     &              c_length, c_rotation,                               &
     &              c_aperture, nom_aperture,                           &
     &              c_offset, c_tilt,                                   &
     &              rcx, rcxp, rcy, rcyp,                               &
     &              rcp, rcs, napx, enom_gev,                           &
     &              part_hit_pos,part_hit_turn,                         &
     &              part_abs_pos,part_abs_turn,                         &
     &              part_impact, part_indiv, part_linteract,            &
     &              onesided,                                           &
!GRD let's also add the FLUKA possibility
     &              flukaname)
            else

!GRD-SR, 09-02-2006
!Force the treatment of the TCDQ equipment as a onsided collimator.
!Both for Beam 1 and Beam 2, the TCDQ is at positive x side.
!              if(db_name1(icoll)(1:4).eq.'TCDQ' ) onesided = .true.
! to treat all collimators onesided 
! -> only for worst case TCDQ studies
               if(db_name1(icoll)(1:4).eq.'TCDQ') onesided = .true.
               if(db_name1(icoll)(1:5).eq.'TCXRP') onesided = .true.
!GRD-SR

!==> SLICE here is possible
!
!     SR, 29-08-2005: Slice the collimator jaws in 'n_slices' pieces
!     using two 4th-order polynomial fits. For each slices, the new
!     gaps and centre are calculates
!     It is assumed that the jaw point closer to the beam defines the
!     nominal aperture.
!
!     SR, 01-09-2005: new official version - input assigned through
!     the 'fort.3' file.
!               if (n_slices.gt.1d0 .and.                                &
!     &              totals.gt.smin_slices .and.                         &
!     &              totals.lt.smax_slices .and.                         &
!     &              db_name1(icoll)(1:4).eq.'TCSG' ) then
!                  if (firstrun) then
!                  write(*,*) 'INFOslice - Collimator ',
!     &              db_name1(icoll), ' sliced in ',n_slices,
!     &              ' pieces!'
!                  endif
!CB
               if (n_slices.gt.1d0 .and.                                &
     &              totals.gt.smin_slices .and.                         &
     &              totals.lt.smax_slices .and.                         &
     &             (db_name1(icoll)(1:4).eq.'TCSG'                      &
     &             .or. db_name1(icoll)(1:3).eq.'TCP'                   &
     &             .or. db_name1(icoll)(1:4).eq.'TCLA'                  &
     &             .or. db_name1(icoll)(1:3).eq.'TCT'                   &
     &             .or. db_name1(icoll)(1:4).eq.'TCLI'                  &
     &             .or. db_name1(icoll)(1:4).eq.'TCL.'
!     RB: added slicing of TCRYO as well    
     &             .or. db_name1(icoll)(1:5).eq.'TCRYO')) then
                      
                  if (firstrun) then
                     write(lout,*) 'INFO> slice - Collimator ',         &
     &                    db_name1(icoll), ' sliced in ',n_slices,      &
     &                    ' pieces !'
                  endif
!
!!     In this preliminary try, all secondary collimators are sliced.
!!     Slice only collimators with finite length!!
!               if (db_name1(icoll)(1:4).eq.'TCSG' .and.
!     &              c_length.gt.0d0 ) then
!!     Slice the primaries, to have more statistics faster!
!!               if (db_name1(icoll)(1:3).eq.'TCP' .and.
!!     +              c_length.gt.0d0 ) then
!!
!!
!!     Calculate longitudinal positions of slices and corresponding heights
!!     and angles from the fit parameters.
!!     -> MY NOTATION: y1_sl: jaw at x > 0; y2_sl: jaw at x < 0;
!!     Note: here, take (n_slices+1) points in order to calculate the
!!           tilt angle of the last slice!!
!                  do jjj=1,n_slices+1
!                     x_sl(jjj) = (jjj-1) * c_length / dble(n_slices)
!                     y1_sl(jjj) =  fit1_1 +                             &
!     &                    fit1_2*x_sl(jjj) +                            &
!     &                    fit1_3*(x_sl(jjj)**2) +                       &
!     &                    fit1_4*(x_sl(jjj)**3) +                       &
!     &                    fit1_5*(x_sl(jjj)**4) +                       &
!     &                    fit1_6*(x_sl(jjj)**5)
!
!                     y2_sl(jjj) = -1d0 * (fit2_1 +                      &
!     &                    fit2_2*x_sl(jjj) +                            &
!     &                    fit2_3*(x_sl(jjj)**2) +                       &
!     &                    fit2_4*(x_sl(jjj)**3) +                       &
!     &                    fit2_5*(x_sl(jjj)**4) +                       &
!     &                    fit2_6*(x_sl(jjj)**5))
!                  enddo

!     CB:10-2007 deformation of the jaws scaled with length
               do jjj=1,n_slices+1
                  x_sl(jjj) = (jjj-1) * c_length / dble(n_slices)
                  y1_sl(jjj) =  fit1_1 +                                &
     &                 fit1_2*x_sl(jjj) +                               &
     &                 fit1_3/c_length*(x_sl(jjj)**2) +                 &
     &                 fit1_4*(x_sl(jjj)**3) +                          &
     &                 fit1_5*(x_sl(jjj)**4) +                          &
     &                 fit1_6*(x_sl(jjj)**5)
!     
                  y2_sl(jjj) = -1d0 * (fit2_1 +                         &
     &                 fit2_2*x_sl(jjj) +                               &
     &                 fit2_3/c_length*(x_sl(jjj)**2) +                 &
     &                 fit2_4*(x_sl(jjj)**3) +                          &
     &                 fit2_5*(x_sl(jjj)**4) +                          &
     &                 fit2_6*(x_sl(jjj)**5))
               enddo

!     Apply the slicing scaling factors (ssf's):
!     
!                  do jjj=1,n_slices+1
!                     y1_sl(jjj) = ssf1 * y1_sl(jjj)
!                     y2_sl(jjj) = ssf2 * y2_sl(jjj)
!                  enddo

!     CB:10-2007 coordinates rotated of the tilt 
                  do jjj=1,n_slices+1
                     y1_sl(jjj) = ssf1 * y1_sl(jjj)
                     y2_sl(jjj) = ssf2 * y2_sl(jjj)
! CB code
                     x1_sl(jjj)=x_sl(jjj)*cos(db_tilt(icoll,1))-        &
     &                    y1_sl(jjj)*sin(db_tilt(icoll,1))
                     x2_sl(jjj)=x_sl(jjj)*cos(db_tilt(icoll,2))-        &
     &                    y2_sl(jjj)*sin(db_tilt(icoll,2))
                     y1_sl(jjj) = y1_sl(jjj)*cos(db_tilt(icoll,1))+     &
     &                    x_sl(jjj)*sin(db_tilt(icoll,1))
                     y2_sl(jjj) = y2_sl(jjj)*cos(db_tilt(icoll,2))+     &
     &                    x_sl(jjj)*sin(db_tilt(icoll,2))
                  enddo

!     Sign of the angle defined differently for the two jaws!
                  do jjj=1,n_slices
                     angle1(jjj) = (( y1_sl(jjj+1) - y1_sl(jjj) ) /     &
     &                    ( x1_sl(jjj+1)-x1_sl(jjj) ))
                     angle2(jjj) =(( y2_sl(jjj+1) - y2_sl(jjj) ) /      &
     &                    ( x2_sl(jjj+1)-x2_sl(jjj) ))
                  enddo
!
!     Sign of the angle defined differently for the two jaws!
!                  do jjj=1,n_slices
!                     angle1(jjj) = ( y1_sl(jjj+1) - y1_sl(jjj) ) /     &
!     &                    (c_length / dble(n_slices) )
!                     angle2(jjj) = ( y2_sl(jjj+1) - y2_sl(jjj) ) /     &
!     &                    (c_length / dble(n_slices) )
!                  enddo
!     For both jaws, look for the 'deepest' point (closest point to beam)
!     Then, shift the vectors such that this closest point defines
!     the nominal aperture
!     Index here must go up to (n_slices+1) in case the last point is the
!     closest (and also for the later calculation of 'a_tmp1' and 'a_tmp2')

!     SR, 01-09-2005: add the recentring flag, as given in 'fort.3' to
!     choose whether recentre the deepest point or not
                  max_tmp = 1e6
                  do jjj=1, n_slices+1
                     if ( y1_sl(jjj).lt.max_tmp ) then
                        max_tmp = y1_sl(jjj)
                     endif
                  enddo

                  do jjj=1, n_slices+1
                     y1_sl(jjj) = y1_sl(jjj) - max_tmp * recenter1      &
     &                    + 0.5 *c_aperture
                  enddo
                  max_tmp = -1e6

                  do jjj=1, n_slices+1
                     if ( y2_sl(jjj).gt.max_tmp ) then
                        max_tmp = y2_sl(jjj)
                     endif
                  enddo

                  do jjj=1, n_slices+1
                     y2_sl(jjj) = y2_sl(jjj) - max_tmp * recenter2      &
     &                    - 0.5 *c_aperture
                  enddo

!!     Check the collimator jaw surfaces (beam frame, before taking into
!!     account the azimuthal angle of the collimator)
                  if (firstrun) then
                    write(lout,*) 'Slicing collimator ',db_name1(icoll)
                     do jjj=1,n_slices
                       write(lout,*) x_sl(jjj), y1_sl(jjj), y2_sl(jjj), &
     &                   angle1(jjj), angle2(jjj), db_tilt(icoll,1),    &
     &                   db_tilt(icoll,2)
                     enddo
                  endif
!
!!     Check the calculation of slice gap and centre
!                  if (firstrun) then
!                     write(*,*) 'Verify centre and gap!'
!                     do jjj=1,n_slices
!                        if ( angle1(jjj).gt.0d0 ) then
!                           a_tmp1 = y1_sl(jjj)
!                        else
!                           a_tmp1 = y1_sl(jjj+1)
!                        endif
!                        if ( angle2(jjj).lt.0d0 ) then
!                           a_tmp2 = y2_sl(jjj)
!                        else
!                           a_tmp2 = y2_sl(jjj+1)
!                        endif
!                        write(*,*) a_tmp1 - a_tmp2,
!     +                       0.5 * ( a_tmp1 + a_tmp2 )
!                     enddo
!                  endif
!
!     Now, loop over the number of slices and call collimate2 each time!
!     For each slice, the corresponding offset and angle are to be used.
                  do jjj=1,n_slices

!     First calculate aperture and centre of the slice
!     Note that:
!     (1)due to our notation for the angle sign,
!     the rotation point of the slice (index j or j+1)
!     DEPENDS on the angle value!!
!     (2) New version of 'collimate2' is required: one must pass
!     the slice number in order the calculate correctly the 's'
!     coordinate in the impact files.

!     Here, 'a_tmp1' and 'a_tmp2' are, for each slice, the closest
!     corners to the beam
                        if ( angle1(jjj).gt.0d0 ) then
                           a_tmp1 = y1_sl(jjj)
                        else
                           a_tmp1 = y1_sl(jjj+1)
                        endif
                        if ( angle2(jjj).lt.0d0 ) then
                           a_tmp2 = y2_sl(jjj)
                        else
                           a_tmp2 = y2_sl(jjj+1)
                        endif
!!     Write down the information on slice centre and offset
!                     if (firstrun) then
!                        write(*,*) 'Processing slice number ',jjj,
!     &                       ' of ',n_slices,' for the collimator ',
!     &                       db_name1(icoll)
!                        write(*,*) 'Aperture [m]= ',
!     &                       a_tmp1 - a_tmp2
!                        write(*,*) 'Offset [m]  = ',
!     &                       0.5 * ( a_tmp1 + a_tmp2 )
!                     endif
!!
!     Be careful! the initial tilt must be added!
!     We leave it like this for the moment (no initial tilt)
!                     c_tilt(1) = c_tilt(1) + angle1(jjj)
!                     c_tilt(2) = c_tilt(2) + angle2(jjj)
                     c_tilt(1) = angle1(jjj)
                     c_tilt(2) = angle2(jjj)
!     New version of 'collimate2' is required: one must pass the
!     slice number in order the calculate correctly the 's'
!     coordinate in the impact files.
!     +                    a_tmp1 - a_tmp2,
!     +                    0.5 * ( a_tmp1 + a_tmp2 ),
! -- TW SEP07 added compatility for tilt, gap and ofset errors to slicing
! -- TW gaprms error is already included in the c_aperture used above  
! -- TW tilt error is added to y1_sl and y2_sl therfore included in 
! -- TW angle1 and angle2 no additinal changes needed 
! -- TW offset error directly added to call of collimate2

! --- TW JUNE08 
                     if (firstrun) then
                        write(55,'(a,1x,i10,5(1x,e13.5),1x,a)')         &
     &                       db_name1(icoll)(1:12),                     &
     &                       jjj,                                       &
     &                       (a_tmp1 - a_tmp2)/2d0,                     &
     &                       0.5 * (a_tmp1 + a_tmp2) + c_offset,        &
     &                       c_tilt(1),                                 &
     &                       c_tilt(2),                                 &
     &                       c_length / dble(n_slices),                 & 
     &                       db_material(icoll)
                     endif
! --- TW JUNE08 
                     call collimate2(c_material,                        &
     &                    c_length / dble(n_slices),                    &
     &                    c_rotation,                                   &
     &                    a_tmp1 - a_tmp2,                              &
     &                    0.5 * ( a_tmp1 + a_tmp2 ) + c_offset,         &
     &                    c_tilt,                                       &
     &                    rcx, rcxp, rcy, rcyp,                         &
     &                    rcp, rcs, napx, enom_gev,                     &
     &                    part_hit_pos, part_hit_turn,                  &
     &                    part_abs_pos, part_abs_turn,                  &
     &                    part_impact, part_indiv,                      &
     &                    part_linteract, onesided, flukaname,          &
     &                    secondary,                                    &
     &                    jjj, nabs_type)
                  enddo
               else
!     Treatment of non-sliced collimators

+if g4collimat
!! Add the geant4 geometry
        if(firstrun.and.iturn.eq.1) then
          call g4_add_collimator(db_name1(icoll), c_material, c_length,
     & c_aperture, c_rotation, c_offset)
        endif

!! Here we do the real collimation
!! First set the correct collimator
        call g4_set_collimator(db_name1(icoll))
        call FLUSH()

!! Loop over all our particles
        g4_lostc = 0
        do j = 1, napx
          if (part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
!! Rotate particles in the frame of the collimator
!! There is more precision if we do it here rather
!! than in the g4 geometry
            x_tmp = rcx(j)
            y_tmp = rcy(j)
            xp_tmp = rcxp(j)
            yp_tmp = rcyp(j)
          rcx(j) = x_tmp*cos(c_rotation) +sin(c_rotation)*y_tmp
          rcy(j) = y_tmp*cos(c_rotation) -sin(c_rotation)*x_tmp
          rcxp(j) = xp_tmp*cos(c_rotation)+sin(c_rotation)*yp_tmp
          rcyp(j) = yp_tmp*cos(c_rotation)-sin(c_rotation)*xp_tmp

!! Call the geant4 collimation function
          call g4_collimate(rcx(j), rcy(j), rcxp(j), rcyp(j), rcp(j))

!! Get the particle back + information
          call g4_collimate_return(rcx(j), rcy(j), rcxp(j), rcyp(j),
     & rcp(j),part_hit(j), part_abs(j), part_impact(j), part_indiv(j), !TODO - fixme, part_hit/part_abs!!
     & part_linteract(j))

!! Rotate back into the accelerator frame
      x_tmp = rcx(j)
      y_tmp = rcy(j)
      xp_tmp = rcxp(j)
      yp_tmp = rcyp(j)
      rcx(j)=x_tmp*cos(-1d0*c_rotation) +sin(-1d0*c_rotation)*y_tmp
      rcy(j)=y_tmp*cos(-1d0*c_rotation) -sin(-1d0*c_rotation)*x_tmp
      rcxp(j)=xp_tmp*cos(-1d0*c_rotation)+sin(-1d0*c_rotation)*yp_tmp
      rcyp(j)=yp_tmp*cos(-1d0*c_rotation)-sin(-1d0*c_rotation)*xp_tmp

          if(part_hit_pos(j).ne.0 .and. part_hit_turn(j).ne.0) then
             part_hit_pos = ie
             part_hit_turn = iturn
          endif

          if(part_abs_pos(j).ne.0 .and. part_abs_turn(j).ne.0) then
            if(dowrite_impact) then
!! FLUKA_impacts.dat
      write(48,'(i4,(1x,f6.3),(1x,f8.6),4(1x,e19.10),i2,2(1x,i7))')     &
     &icoll,c_rotation,                                                 &
     &0.0,                                                              &!hr09
     &0.0,0.0,0.0,0.0,                                                  &
     &part_abs(j),flukaname(j),iturn  !!! TODO -  It should not be part_abs here.
              endif

            part_abs_pos(j)  = ie
            part_abs_turn(j) = iturn
            rcx(j) = 99.99d-3
            rcy(j) = 99.99d-3
            g4_lostc = g4_lostc + 1
          endif

          call FLUSH()

          endif !part_abs_pos(j) .ne. 0 .and. part_abs_turn(j) .ne. 0
        enddo
!      write(lout,*) 'COLLIMATOR LOSSES ', db_name1(icoll), g4_lostc
+ei
+if .not.g4collimat
                  call collimate2(c_material, c_length, c_rotation,     &
     &                 c_aperture, c_offset, c_tilt,                    &
     &                 rcx, rcxp, rcy, rcyp,                            &
     &                 rcp, rcs, napx, enom_gev,                        &
     &                 part_hit_pos,part_hit_turn,                      &
     &                 part_abs_pos, part_abs_turn,                     &
     &                 part_impact, part_indiv, part_linteract,         &
     &                 onesided, flukaname, secondary, 1, nabs_type)    &
+ei
               endif !if (n_slices.gt.1d0 .and.

               endif !if(db_name1(icoll)(1:4).eq.'COLM') then
          endif !if (found) then
      end

            
!>
!! collimate_end_collimator()
!! This routine is called at the exit of a collimator
!<
      subroutine collimate_end_collimator()
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
      integer i,ix,j,jj,jx,kpz,kzz,napx0,nbeaux,nmz,nthinerr
      double precision benkcc,cbxb,cbzb,cikveb,crkveb,crxb,crzb,r0,r000,&
     &r0a,r2b,rb,rho2b,rkb,tkb,xbb,xrb,zbb,zrb
      logical lopen
+ca parpro
+ca parnum
+ca common
+ca commons
+ca commont1
+ca commondl
+ca commonxz
+ca commonta
+ca commonmn
+ca commonm1
+ca commontr
+ca beamdim
+ca commonex
      dimension nbeaux(nbb)
+if collimat
+ca collpara
+ca dbtrthin
+ca database
+ca dbcommon
+ca dblinopt
+ca dbpencil
+ca info
+ca dbcolcom
+ca dbthin6d

+ei

+if bnlelens
+ca rhicelens
+ei

      double precision c5m4,stracki

+if fast
      c5m4=5.0d-4
+ei



!++  Output information:
!++
!++  PART_HIT_POS (MAX_NPART)  Hit flag for last hit
!++  PART_HIT_TURN(MAX_NPART)  Hit flag for last hit
!++  PART_ABS_POS (MAX_NPART)  Abs flag
!++  PART_ABS_TURN(MAX_NPART)  Abs flag
!++  PART_IMPACT  (MAX_NPART)  Impact parameter (0 for inner face)
!++  PART_INDIV   (MAX_NPART)  Divergence of impacting particles
!------------------------------------------------------------------------------
!++  Calculate average impact parameter and save info for all
!++  collimators. Copy information back and do negative drift.
          n_impact = 0
          n_absorbed = 0
          sum      = 0d0
          sqsum    = 0d0

!++  Copy particle data back and do path length stuff; check for absorption
!++  Add orbit offset back.
          do j = 1, napx

!APRIL2005 IN ORDER TO GET RID OF NUMERICAL ERRORS, JUST DO THE TREATMENT FOR
!APRIL2005 IMPACTING PARTICLES...
             if (part_hit_pos(j) .eq.ie .and.
     &           part_hit_turn(j).eq.iturn    ) then
!++  For zero length element track back half collimator length
! DRIFT PART
!       write(lout,*) j, ' hit ', part_hit(j)
              if (stracki.eq.0.) then
!       write(lout,*) j, ' backtrack '
                if(iexact.eq.0) then
                  rcx(j)  = rcx(j) - 0.5d0*c_length*rcxp(j)
                  rcy(j)  = rcy(j) - 0.5d0*c_length*rcyp(j)
                else
                  zpj=sqrt(1d0-rcxp(j)**2-rcyp(j)**2)
                  rcx(j) = rcx(j) - 0.5d0*c_length*(rcxp(j)/zpj)
                  rcy(j) = rcy(j) - 0.5d0*c_length*(rcyp(j)/zpj)
                endif
              endif

!++  Now copy data back to original verctor
              xv(1,j) = rcx(j)*1d3  +torbx(ie)
              yv(1,j) = rcxp(j)*1d3 +torbxp(ie)
              xv(2,j) = rcy(j)*1d3  +torby(ie)
              yv(2,j) = rcyp(j)*1d3 +torbyp(ie)
              ejv(j) = rcp(j)*1d3

!++  Energy update, as recommended by Frank
              ejfv(j)=sqrt(ejv(j)*ejv(j)-pma*pma)
              rvv(j)=(ejv(j)*e0f)/(e0*ejfv(j))
              dpsv(j)=(ejfv(j)-e0f)/e0f
              oidpsv(j)=one/(one+dpsv(j))
              dpsv1(j)=dpsv(j)*c1e3*oidpsv(j)
              yv(1,j)=ejf0v(j)/ejfv(j)*yv(1,j)
              yv(2,j)=ejf0v(j)/ejfv(j)*yv(2,j)

!APRIL2005 ...OTHERWISE JUST GET BACK FORMER COORDINATES
            else
              xv(1,j) = rcx0(j)*1d3+torbx(ie)
              yv(1,j) = rcxp0(j)*1d3+torbxp(ie)
              xv(2,j) = rcy0(j)*1d3+torby(ie)
              yv(2,j) = rcyp0(j)*1d3+torbyp(ie)
              ejv(j) = rcp0(j)*1d3
            endif
!APRIL2005
! 
!TW for roman pot checking
!            if(icoll.eq.73) then
!               do j = 1,napx 
!                  write(9998,*)flukaname(j),rcx0(j),rcy0(j),rcx(j),     &
!     &rcy(j),rcxp0(j),rcyp0(j),rcxp(j),rcyp(j)
!               enddo
!            elseif(icoll.eq.74) then
!               do j = 1,napx 
!                  write(9999,*)flukaname(j),rcx0(j),rcy0(j),rcx(j),     &
!     &rcy(j),rcxp0(j),rcyp0(j),rcxp(j),rcyp(j)
!               enddo
!            endif
!
!++  Write trajectory for any selected particle
!
!!            if (firstrun) then
!!              if (rselect.gt.0 .and. rselect.lt.65) then
!            DO j = 1, NAPX
!
!!              xj     = (xv(1,j)-torbx(ie))/1d3
!!              xpj    = (yv(1,j)-torbxp(ie))/1d3
!!              yj     = (xv(2,j)-torby(ie))/1d3
!!              ypj    = (yv(2,j)-torbyp(ie))/1d3
!!              pj     = ejv(j)/1d3
!GRD
!07-2006 TEST
!!              if (iturn.eq.1.and.j.eq.1) then
!!              sum_ax(ie)=0d0
!!              sum_ay(ie)=0d0
!!              endif
!GRD
!
!!              gammax = (1d0 + talphax(ie)**2)/tbetax(ie)
!!              gammay = (1d0 + talphay(ie)**2)/tbetay(ie)
!
!!             if (part_abs(j).eq.0) then
!!          nspx    = sqrt(                                               &
!!     &abs( gammax*(xj)**2 +                                             &
!!     &2d0*talphax(ie)*xj*xpj +                                          &
!!     &tbetax(ie)*xpj**2 )/myemitx0                                      &
!!     &)
!!                nspy    = sqrt(                                         &
!!     &abs( gammay*(yj)**2 +                                             &
!!     &2d0*talphay(ie)*yj*ypj +                                          &
!!     &tbetay(ie)*ypj**2 )/myemity0                                      &
!!     &)

!++  First check for particle interaction at this collimator and this turn
            if (part_hit_pos (j).eq.ie .and.
     &          part_hit_turn(j).eq.iturn    ) then

!++  Fill the change in particle angle into histogram
              if(dowrite_impact) then
                write(46,'(i8,1x,i4,1x,f8.2)')                          &
     &               ipart(j)+100*samplenumber,iturn,sampl(ie)
              endif

              ! Particle has impacted
              if(part_abs_pos(j) .ne.0 .and.
     &           part_abs_turn(j).ne.0      ) then
                if(dowrite_impact) then
                  write(47,'(i8,1x,i4,1x,f8.2)')                        &
     &ipart(j)+100*samplenumber,iturn,sampl(ie)
                endif
+if hdf5
       hdfpid=ipart(j)+100*samplenumber
       hdfturn=iturn
       hdfs=sampl(ie)-0.5*c_length
       hdfx=(rcx0(j)*1d3+torbx(ie))-0.5*c_length*(rcxp0(j)*1d3+         &
     &      torbxp(ie))
       hdfxp=rcxp0(j)*1d3+torbxp(ie)
       hdfy=(rcy0(j)*1d3+torby(ie))-0.5*c_length*(rcyp0(j)*1d3+         &
     &      torbyp(ie))
       hdfyp=rcyp0(j)*1d3+torbyp(ie)
       hdfdee=(ejv(j)-myenom)/myenom
       hdftyp=secondary(j)+tertiary(j)+other(j)
       CALL APPENDREADING(hdfpid,hdfturn,hdfs,hdfx,hdfxp,hdfy,hdfyp,    &
     &                    hdfdee,hdftyp)
+ei
+if .not.hdf5
       write(38,'(1x,i8,1x,i4,1x,f10.2,4(1x,e11.5),1x,e11.3,1x,i4)')
     &ipart(j)+100*samplenumber,iturn,sampl(ie)-0.5*c_length,           &
     &(rcx0(j)*1d3+torbx(ie))-0.5*c_length*(rcxp0(j)*1d3+torbxp(ie)),   &
     &rcxp0(j)*1d3+torbxp(ie),                                          &
     &(rcy0(j)*1d3+torby(ie))-0.5*c_length*(rcyp0(j)*1d3+torbyp(ie)),   &
     &rcyp0(j)*1d3+torbyp(ie),                                          &
     &(ejv(j)-myenom)/myenom,secondary(j)+tertiary(j)+other(j)
+ei

              !Here we've found a newly hit particle
              elseif (part_abs_pos (j).eq.0 .and.
     &                part_abs_turn(j).eq.0       ) then
                 xkick = rcxp(j) - rcxp0(j)
                 ykick = rcyp(j) - rcyp0(j)
                 
                 if (db_name1(icoll)(1:3).eq.'TCP'.or.                  &
     &                db_name1(icoll)(1:4).eq.'COLM'.or.                &
     &                db_name1(icoll)(1:5).eq.'COLH0'.or.               &
     &                db_name1(icoll)(1:5).eq.'COLV0') then
                    secondary(j) = 1
                 elseif (db_name1(icoll)(1:3).eq.'TCS'.or.              &
     &                   db_name1(icoll)(1:4).eq.'COLH1'.or.            &
     &                   db_name1(icoll)(1:4).eq.'COLV1'.or.            &
     &                   db_name1(icoll)(1:4).eq.'COLH2') then
                    tertiary(j)  = 2
                 elseif ((db_name1(icoll)(1:3).eq.'TCL').or.            &
     &                   (db_name1(icoll)(1:3).eq.'TCT').or.            &
     &                   (db_name1(icoll)(1:3).eq.'TCD').or.            &
     &                   (db_name1(icoll)(1:3).eq.'TDI')) then
                    other(j)     = 4
                 endif
              else
                 write(lout,*) "Error in collimate_end_collimator"
                 write(lout,*) "Particle cannot be both absorbed"//
     &                " and not absorbed."
                 write(lout,*) part_abs_pos (j),  part_abs_turn(j)
                 call prror(-1)
              endif

!GRD THIS LOOP MUST NOT BE WRITTEN INTO THE "IF(FIRSTRUN)" LOOP !!!!!
      if (dowritetracks) then
        if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
          if ((secondary(j).eq.1.or.tertiary(j).eq.2.or.other(j).eq.4)  &
     & .and.(xv(1,j).lt.99d0 .and. xv(2,j).lt.99d0) .and.               &
!GRD HERE WE APPLY THE SAME KIND OF CUT THAN THE SIGSECUT PARAMETER
     &(                                                                 &
     &((                                                                &
     &(xv(1,j)*1d-3)**2                                                 &
     &/                                                                 &
     &(tbetax(ie)*myemitx0_collgap)
     &).ge.dble(sigsecut2)).or.                                         &
     &((                                                                &
     &(xv(2,j)*1d-3)**2                                                 &
     &/                                                                 &
     &(tbetay(ie)*myemity0_collgap)
     &).ge.dble(sigsecut2)).or.                                         &
     &(((xv(1,j)*1d-3)**2/(tbetax(ie)*myemitx0_collgap))+
     &((xv(2,j)*1d-3)**2/(tbetay(ie)*myemity0_collgap))
     &.ge.sigsecut3)                                                    &
     &) ) then

          xj     = (xv(1,j)-torbx(ie))/1d3
          xpj    = (yv(1,j)-torbxp(ie))/1d3
          yj     = (xv(2,j)-torby(ie))/1d3
          ypj    = (yv(2,j)-torbyp(ie))/1d3

+if hdf5
!       We write trajectories before and after element in this case.
       hdfpid=ipart(j)+100*samplenumber
       hdfturn=iturn
       hdfs=sampl(ie)-0.5*c_length
       hdfx=  ! xv(1,j)-0.5*c_length*yv(1,j)
     &    (rcx0(j)*1d3+torbx(ie))-0.5*c_length*(rcxp0(j)*1d3+torbxp(ie))
       hdfxp= ! yv(1,j)
     &    rcxp0(j)*1d3+torbxp(ie)
       hdfy=  ! xv(2,j)-0.5*c_length*yv(2,j)
     &    (rcy0(j)*1d3+torby(ie))-0.5*c_length*(rcyp0(j)*1d3+torbyp(ie))
       hdfyp= ! yv(2,j)
     &    rcyp0(j)*1d3+torbyp(ie)
       hdfdee=(ejv(j)-myenom)/myenom
       hdftyp=secondary(j)+tertiary(j)+other(j)
       call APPENDREADING(hdfpid,hdfturn,hdfs,hdfx,hdfxp,hdfy,hdfyp,    &
     &                    hdfdee,hdftyp)
       hdfs=sampl(ie)+0.5*c_length
       hdfx=xv(1,j)+0.5*c_length*yv(1,j)
       hdfxp=yv(1,j)
       hdfy=xv(2,j)+0.5*c_length*yv(2,j)
       hdfyp=yv(2,j)
       call APPENDREADING(hdfpid,hdfturn,hdfs,hdfx,hdfxp,hdfy,hdfyp,    &
     &                    hdfdee,hdftyp)
     
+ei
+if .not.hdf5
       write(38,'(1x,i8,1x,i4,1x,f10.2,4(1x,e11.5),1x,e11.3,1x,i4)')
     &ipart(j)+100*samplenumber,iturn,sampl(ie)-0.5*c_length,           &
     &(rcx0(j)*1d3+torbx(ie))-0.5*c_length*(rcxp0(j)*1d3+torbxp(ie)),   &
     &rcxp0(j)*1d3+torbxp(ie),                                          &
     &(rcy0(j)*1d3+torby(ie))-0.5*c_length*(rcyp0(j)*1d3+torbyp(ie)),   &
     &rcyp0(j)*1d3+torbyp(ie),                                          &
     &(ejv(j)-myenom)/myenom,secondary(j)+tertiary(j)+other(j)
       
       write(38,'(1x,i8,1x,i4,1x,f10.2,4(1x,e11.5),1x,e11.3,1x,i4)')
     &ipart(j)+100*samplenumber,iturn,sampl(ie)+0.5*c_length,           &
     &xv(1,j)+0.5*c_length*yv(1,j),yv(1,j),                             &
     &xv(2,j)+0.5*c_length*yv(2,j),yv(2,j),(ejv(j)-myenom)/myenom,      &
     &secondary(j)+tertiary(j)+other(j)
+ei
          endif
        endif
      endif

!++  Calculate impact observables, fill histograms, save collimator info, ...
              n_impact = n_impact + 1
              sum = sum + part_impact(j)
              sqsum = sqsum + part_impact(j)**2
              cn_impact(icoll) = cn_impact(icoll) + 1
              csum(icoll) = csum(icoll) + part_impact(j)
              csqsum(icoll) = csqsum(icoll) + part_impact(j)**2

!++  If the interacting particle was lost, add-up counters for absorption
!++  Note: a particle with x/y >= 99. never hits anything any more in
!++        the logic of this program. Be careful to always fulfill this!
              if (part_abs_pos(j).ne.0 .and. part_abs_turn(j).ne.0) then
                 n_absorbed = n_absorbed + 1
                 cn_absorbed(icoll) = cn_absorbed(icoll) + 1
                 n_tot_absorbed = n_tot_absorbed + 1
                 iturn_last_hit = part_hit_before_turn(j)
                 iturn_absorbed = part_hit_turn(j)
                 if (iturn_last_hit.eq.0) then
                    iturn_last_hit = iturn_absorbed
                    iturn_survive  = iturn_absorbed - iturn_last_hit
                 endif
              endif
                 
!++  End of check for hit this turn and element
           endif
        end do ! end do j = 1, napx

!++  Calculate statistical observables and save into files...
          if (n_impact.gt.0) then
            average = sum/n_impact

            if (sqsum/n_impact.ge.average**2) then
              sigma = sqrt(sqsum/n_impact - average**2)
            else
              sigma = 0d0
            endif

          else
            average = 0d0
            sigma   = 0d0
          endif

          if (cn_impact(icoll).gt.0) then
            caverage(icoll) = csum(icoll)/cn_impact(icoll)

            if ((caverage(icoll)**2).gt.                                &
     &(csqsum(icoll)/cn_impact(icoll))) then
               csigma(icoll) = 0
            else
              csigma(icoll) = sqrt(csqsum(icoll)/                       &
     &cn_impact(icoll) - caverage(icoll)**2)
            endif

          endif

!-----------------------------------------------------------------
!++  For a  S E L E C T E D  collimator only consider particles that
!++  were scattered on this selected collimator at the first turn. All
!++  other particles are discarded.
!++  - This is switched on with the DO_SELECT flag in the input file.
!++  - Note that the part_select(j) flag defaults to 1 for all particles.

! should name_sel(1:11) extended to allow longer names as done for 
! coll the coll_ellipse.dat file !!!!!!!!
           if (((db_name1(icoll).eq.name_sel(1:11))                     &
     &.or.(db_name2(icoll).eq.name_sel(1:11)))                          &
     &.and. iturn.eq.1  ) then
            num_selhit = 0
            num_surhit = 0
            num_selabs = 0

            do j = 1, napx
               if( part_hit_pos (j).eq.ie .and.
     &             part_hit_turn(j).eq.iturn    ) then
               
                num_selhit = num_selhit+1
                if (part_abs_pos(j) .eq.0 .and.
     &              part_abs_turn(j).eq.0       ) then
                  num_surhit = num_surhit+1
                else
                  num_selabs = num_selabs + 1
                endif

!++  If we want to select only partciles interacting at the specified
!++  collimator then remove all other particles and reset the number
!++  of the absorbed particles to the selected collimator.
              elseif (do_select.and.firstrun) then
                part_select(j) = 0
                n_tot_absorbed = num_selabs
              endif
            end do

!++  Calculate average impact parameter and save distribution into file
!++  only for selected collimator
            n_impact = 0
            sum      = 0d0
            sqsum    = 0d0

            do j = 1, napx
               if ( part_hit_pos (j).eq.ie .and.
     &              part_hit_turn(j).eq.iturn    ) then
                  if (part_impact(j).lt.-0.5d0) then
                     write(lout,*)
     &                    'ERR>  Found invalid impact parameter!',
     &                    part_impact(j)
                     write(outlun,*)
     &                    'ERR>  Invalid impact parameter!',
     &                    part_impact(j)
                     call prror(-1)
                  endif
                  n_impact = n_impact + 1
                  sum = sum + part_impact(j)
                  sqsum = sqsum + part_impact(j)**2
                  if (part_hit_pos (j).ne.0 .and.
     &                part_hit_turn(j).ne.0 .and.
     &                dowrite_impact              ) then
                     write(49,*) part_impact(j), part_indiv(j)
                  endif
               endif
            end do
            if (n_impact.gt.0) then
              average = sum/n_impact
              if(sqsum/n_impact.ge.average**2) then
           sigma = sqrt(sqsum/n_impact - average**2)
              else
                 sigma = 0d0
              endif
            endif
!
!++  Some information
            write(lout,*)
     &'INFO>  Selected collimator had N hits. N: ',                     &
     &num_selhit
            write(lout,*)
     &'INFO>  Number of impacts                : ',                     &
     &n_impact
            write(lout,*)
     &'INFO>  Number of escaped protons        : ',                     &
     &num_surhit
            write(lout,*)
     &'INFO>  Average impact parameter [m]     : ',                     &
     &average
            write(lout,*)
     &'INFO>  Sigma impact parameter [m]       : ',                     &
     &sigma

            if (dowrite_impact) close(49)

!++  End of    S E L E C T E D   collimator
          endif

      end

!>
!! collimate_end_sample()
!! This routine is called from trauthin after each sample
!! has been tracked by thin6d
!<
      subroutine collimate_end_sample(j)
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
      integer i,ix,j,jj,jx,kpz,kzz,napx0,nbeaux,nmz,nthinerr
      double precision benkcc,cbxb,cbzb,cikveb,crkveb,crxb,crzb,r0,r000,&
     &r0a,r2b,rb,rho2b,rkb,tkb,xbb,xrb,zbb,zrb
      logical lopen
+ca parpro
+ca parnum
+ca common
+ca commons
+ca commont1
+ca commondl
+ca commonxz
+ca commonta
+ca commonmn
+ca commonm1
+ca commontr
+ca beamdim
      dimension nbeaux(nbb)

+if collimat
+ca collpara
+ca dbtrthin
+ca database
+ca dbcommon
+ca dblinopt
+ca dbpencil
+ca info
+ca dbcolcom
+ei

+if bnlelens
+ca rhicelens
+ei

!++  Save particle offsets to a file
      close(42)
      close(44)

      if(dowrite_impact) close(49)

      if(dowritetracks) then
        if(cern) close(38)
      endif

!------------------------------------------------------------------------
!++  Write the number of absorbed particles
      write(outlun,*) 'INFO>  Number of impacts             : ',        &
     &n_tot_absorbed+nsurvive_end
      write(outlun,*) 'INFO>  Number of impacts at selected : ',        &
     &num_selhit
      write(outlun,*) 'INFO>  Number of surviving particles : ',        &
     &nsurvive_end
      write(outlun,*) 'INFO>  Number of absorbed particles  : ',        &
     &n_tot_absorbed

      write(outlun,*)
!GRD UPGRADE JANUARY 2005
      if(n_tot_absorbed.ne.0) then                                       !hr08
      write(outlun,*) ' INFO>  Eff_r @  8 sigma    [e-4] : ',           &
     &(neff(5)/dble(n_tot_absorbed))/1d-4                                !hr08
      write(outlun,*) ' INFO>  Eff_r @ 10 sigma    [e-4] : ',           &
     &(neff(9)/dble(n_tot_absorbed))/1d-4                                !hr08
      write(outlun,*) ' INFO>  Eff_r @ 10-20 sigma [e-4] : ',           &
     &((neff(9)-neff(19))/(dble(n_tot_absorbed)))/1d-4                   !hr08
      write(outlun,*)
      write(outlun,*) neff(5)/dble(n_tot_absorbed),                     &
     &neff(9)/dble(n_tot_absorbed),                                     &
     &(neff(9)-neff(19))/(dble(n_tot_absorbed)), ' !eff'
      write(outlun,*)

!UPGRADE JANUARY 2005
      else
          write(lout,*) 'NO PARTICLE ABSORBED'
      endif

      write(lout,*)
      write(lout,*) 'INFO>  Number of impacts             : ',          &
     &n_tot_absorbed+nsurvive_end
      write(lout,*) 'INFO>  Number of impacts at selected : ',
     &num_selhit
      write(lout,*) 'INFO>  Number of surviving particles : ',          &
     &nsurvive_end
      write(lout,*) 'INFO>  Number of absorbed particles  : ',
     &n_tot_absorbed
      write(lout,*)
!GRD UPGRADE JANUARY 2005
      if(n_tot_absorbed.ne.0) then                                       !hr08
      write(lout,*) ' INFO>  Eff_r @  8 sigma    [e-4] : ',
     &(neff(5)/dble(n_tot_absorbed))/1d-4                               !hr08
      write(lout,*) ' INFO>  Eff_r @ 10 sigma    [e-4] : ',
     &(neff(9)/dble(n_tot_absorbed))/1d-4                                !hr08
      write(lout,*) ' INFO>  Eff_r @ 10-20 sigma [e-4] : ',
     &((neff(9)-neff(19))/dble(n_tot_absorbed))/1d-4                     !hr08
      write(lout,*)
!UPGRADE JANUARY 2005
      else
          write(lout,*) 'NO PARTICLE ABSORBED'
      endif
!
!------------------------------------------------------------------------
!++  Write efficiency file
!
      inquire( unit=1991, opened=lopen)
      if (lopen) then
         write(lout,*)
     &        "ERROR in efficiency.dat: FILE 1991 already taken"
        call prror(-1)
      endif
      open(unit=1991, file='efficiency.dat')
!UPGRADE JANUARY 2005
      if(n_tot_absorbed.ne.0) then
      write(1991,*)                                                     &
     &'# 1=rad_sigma 2=frac_x 3=frac_y 4=frac_r'
      do k=1,numeff
        write(1991,'(7(1x,e15.7),1x,I5)') rsig(k),                      &
     &neffx(k)/dble(n_tot_absorbed),                                    &
     &neffy(k)/dble(n_tot_absorbed),                                    &
     &neff(k)/dble(n_tot_absorbed),                                     &
     &neffx(k),                                                         &
     &neffy(k),                                                         &
     &neff(k), n_tot_absorbed
      end do
      else
          write(lout,*) 'NO PARTICLE ABSORBED'
      endif
!END OF UPGRADE
      close(1991)
!!------------------------------------------------------------------------
!++  Write efficiency vs dp/p file

      inquire( unit=1992, opened=lopen )
      if (lopen) then
           write(lout,*)
     &          "ERROR in efficiency_dpop.dat: FILE 1992 already taken"
        call prror(-1)
      endif

      open(unit=1992, file='efficiency_dpop.dat')
!UPGRADE 4/11/2014
      if(n_tot_absorbed.ne.0) then
      write(1992,*)                                                       &
     &'# 1=dp/p 2=n_dpop/tot_nabs 3=n_dpop 4=tot_nabs 5=npart' 

      do k=1,numeffdpop
        write(1992,'(3(1x,e15.7),2(1x,I5))') dpopbins(k),               &
     &neffdpop(k)/dble(n_tot_absorbed),                                 &
     &neffdpop(k), n_tot_absorbed, npartdpop(k)
      end do

      else
          write(lout,*) 'NO PARTICLE ABSORBED'
      endif
!END OF UPGRADE
      close(1992)
!!------------------------------------------------------------------------
!++  Write 2D efficiency file (eff vs. A_r and dp/p)
      inquire( unit=1993, opened=lopen )
      if (lopen) then
         write(lout,*)
     &        "ERROR in efficiency_2d.dat:FILE 1993 already taken"
        call prror(-1)
      endif

      open(unit=1993, file='efficiency_2d.dat')
      if(n_tot_absorbed.ne.0) then
      write(1993,*)                                                       &
     &'# 1=rad_sigma 2=dp/p 3=n/tot_nabs 4=n 5=tot_nabs' 
      do i=1,numeff
        do k=1,numeffdpop
          write(1993,'(4(1x,e15.7),1(1x,I5))') rsig(i),  dpopbins(k),    &
     &neff2d(i,k)/dble(n_tot_absorbed),                                 &
     &neff2d(i,k), n_tot_absorbed
        end do
      end do
      else
          write(lout,*) 'NO PARTICLE ABSORBED'
      endif
!END OF UPGRADE
      close(1993)
!!------------------------------------------------------------------------
!------------------------------------------------------------------------
!++  Write collimation summary file
!
      open(unit=50, file='coll_summary.dat')
      write(50,*)                                                       &
     &'# 1=icoll 2=collname 3=nimp 4=nabs 5=imp_av 6=imp_sig 7=length'

      do icoll = 1, db_ncoll
        if(db_length(icoll).gt.0d0) then
          write(50,'(i4,1x,a,2(1x,i5),2(1x,e15.7),3x,f3.1)')            &
     &icoll, db_name1(icoll),cn_impact(icoll), cn_absorbed(icoll),      &
     &caverage(icoll), csigma(icoll),db_length(icoll)
        endif
      end do

      close(50)
      end

!>
!! collimate_exit()
!! This routine is called once at the end of the simulation and
!! can be used to do any final postrocessing and/or file saving.
!<
      subroutine collimate_exit()
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
      integer i,ix,j,jb,jj,jx,kpz,kzz,napx0,nbeaux,nmz,nthinerr
      double precision benkcc,cbxb,cbzb,cikveb,crkveb,crxb,crzb,r0,r000,&
     &r0a,r2b,rb,rho2b,rkb,tkb,xbb,xrb,zbb,zrb
      logical lopen
+ca parpro
+ca parnum
+ca common
+ca commons
+ca commont1
+ca commondl
+ca commonxz
+ca commonta
+ca commonmn
+ca commonm1
+ca commontr
+ca beamdim
      dimension nbeaux(nbb)
+if collimat
+ca collpara
+ca dbtrthin
+ca database
+ca dbcommon
+ca dblinopt
+ca dbpencil
+ca info
+ca dbcolcom
+ei
+if bnlelens
+ca rhicelens
+ei

      close(outlun)
      close(43)

      if(dowritetracks) then
        if(.not. cern) close(38)
        if(name_sel(1:3).eq.'COL') close(555)
      endif

      if(do_select) then
         close(45)
      endif

      if(dowrite_impact) then
        close(46)
        close(47)
        close(48)
        close(4801)
        close(3998)
        close(39)
      endif

      open(unit=56, file='amplitude.dat')
      open(unit=51, file='amplitude2.dat')
!UPGRADE JANUARY 2005
      open(unit=57, file='betafunctions.dat')
      if(dowrite_amplitude) then
      write(56,*)                                                       &
     &'# 1=ielem 2=name 3=s 4=AX_AV 5=AX_RMS 6=AY_AV 7=AY_RMS',         &
     &'8=alphax 9=alphay 10=betax 11=betay 12=orbitx',                  &
     &'13=orbity 14=tdispx 15=tdispy',                                  &
     &'16=xbob 17=ybob 18=xpbob 19=ypbob'
      do i=1,iu
        write(56,'(i4, (1x,a16), 17(1x,e20.13))')                       &!hr08
     &i, ename(i), sampl(i),                                            &!hr08
     &sum_ax(i)/dble(max(nampl(i),1)),                                  &!hr08
     &sqrt(abs((sqsum_ax(i)/dble(max(nampl(i),1)))-                     &!hr08
     &(sum_ax(i)/dble(max(nampl(i),1)))**2)),                           &!hr08
     &sum_ay(i)/dble(max(nampl(i),1)),                                  &!hr08
     &sqrt(abs((sqsum_ay(i)/dble(max(nampl(i),1)))-                     &!hr08
     &(sum_ay(i)/dble(max(nampl(i),1)))**2)),                           &!hr08
     &talphax(i), talphay(i),                                           &!hr08
     &tbetax(i), tbetay(i), torbx(i), torby(i),                         &!hr08
     &tdispx(i), tdispy(i),                                             &!hr08
     &xbob(i),ybob(i),xpbob(i),ypbob(i)                                  !hr08
      end do
!GRD
      write(51,*)                                                       &
     &'# 1=ielem 2=name 3=s 4=ORBITX',                                  &
     &'5=orbity 6=tdispx 7=tdispy',                                     &
     &'8=xbob 9=ybob 10=xpbob 11=ypbob'
      do i=1,iu
        write(51,'(i4, (1x,a16), 9(1x,e15.7))')                         &
     &i, ename(i), sampl(i),                                            &
     &torbx(i), torby(i),                                               &
     &tdispx(i), tdispy(i),                                             &
     &xbob(i),ybob(i),xpbob(i),ypbob(i)
      end do
!GRD UPGRADE
      write(57,*)                                                       &
     &'# 1=ielem 2=name       3=s             4=TBETAX(m)     5=TBETAY(m
     &)     6=TORBX(mm)    7=TORBY(mm)     8=TORBXP(mrad)   9=TORBYP(mrad
     &)  10=TDISPX(m)  11=MUX()    12=MUY()'


      do i=1,iu
        write(57,'(i4, (1x,a16), 10(1x,e15.7))')                         &
     &      i, ename(i), sampl(i),                                            &
     &      tbetax(i), tbetay(i), 
     &      torbx(i), torby(i), torbxp(i), torbyp(i), tdispx(i), mux(i), 
     &      muy(i)            ! RB: added printout of closed orbit and angle

      end do
      endif

      close(56)
      close(51)
      close(57)
!GRD END OF UPGRADE
!GRD
!      DO J=1,iu
!        DO I=1,numl
!        xaveragesumoverturns(j)  = xaverage(j,i)
!     &                             + xaverage(j,MAX((i-1),1))
!        yaveragesumoverturns(j)  = yaverage(j,i)
!     &                             + yaverage(j,MAX((i-1),1))
!        xpaveragesumoverturns(j) = xpaverage(j,i)
!     &                             + xpaverage(j,MAX((i-1),1))
!        ypaveragesumoverturns(j) = ypaverage(j,i)
!     &                             + ypaverage(j,MAX((i-1),1))
!        END DO
!        xclosedorbitcheck(j)=(xaveragesumoverturns(j)
!     &                        +xaverage(j,numl))/(2*numl)
!        yclosedorbitcheck(j)=(yaveragesumoverturns(j)
!     &                        +yaverage(j,numl))/(2*numl)
!        xpclosedorbitcheck(j)=(xpaveragesumoverturns(j)
!     &                        +xpaverage(j,numl))/(2*numl)
!        ypclosedorbitcheck(j)=(ypaveragesumoverturns(j)
!     &                        +ypaverage(j,numl))/(2*numl)
!      END DO
!
!      OPEN(unit=99, file='xchecking.dat')
!      WRITE(99,*) '# 1=s 2=x 3=xp 4=y 5=yp'
!      DO J=1,iu
!      WRITE(99,'(i, 5(1x,e15.7))')
!     &     j, SAMPL(j),
!     &     xclosedorbitcheck(j), xpclosedorbitcheck(j),
!     &     yclosedorbitcheck(j), ypclosedorbitcheck(j)
!      END DO
!      CLOSE(99)
!GRD
!GRD WE CAN ALSO MAKE AN ORBIT CHECKING
!GRD
      open(unit=99, file='orbitchecking.dat')
      write(99,*) '# 1=s 2=torbitx 3=torbity'
      do j=1,iu
      write(99,'(i4, 3(1x,e15.7))')                                     &
     &j, sampl(j),torbx(j), torby(j)
      end do
      close(99)

+if g4collimat
      call g4_terminate()
+ei
      end

!>
!! This routine is called at the start of each tracking turn
!<
      subroutine collimate_start_turn(n)
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
      integer i,ix,j,jb,jj,jx,kpz,kzz,napx0,nbeaux,nmz,nthinerr
      double precision benkcc,cbxb,cbzb,cikveb,crkveb,crxb,crzb,r0,r000,&
     &r0a,r2b,rb,rho2b,rkb,tkb,xbb,xrb,zbb,zrb
      logical lopen
+ca parpro
+ca parnum
+ca common
+ca commons
+ca commont1
+ca commondl
+ca commonxz
+ca commonta
+ca commonmn
+ca commonm1
+ca commontr
+ca beamdim
      dimension nbeaux(nbb)
+if collimat
+ca collpara
+ca dbtrthin
+ca database
+ca dbcommon
+ca dblinopt
+ca dbpencil
+ca info
+ca dbcolcom
+ei
      integer n

      iturn=n
      totals=0d0 !This keeps track of the s position of the current element,
                   ! which is also done by cadcum
      end

!>
!! This routine is called at the start of every element
!<
      subroutine collimate_start_element(i)
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
      integer i,ix,j,jj,jx,kpz,kzz,napx0,nbeaux,nmz,nthinerr
      double precision benkcc,cbxb,cbzb,cikveb,crkveb,crxb,crzb,r0,r000,&
     &r0a,r2b,rb,rho2b,rkb,tkb,xbb,xrb,zbb,zrb
      logical lopen
+ca parpro
+ca parnum
+ca common
+ca commons
+ca commont1
+ca commondl
+ca commonxz
+ca commonta
+ca commonmn
+ca commonm1
+ca commontr
+ca beamdim
      dimension nbeaux(nbb)
+if collimat
+ca collpara
+ca dbtrthin
+ca database
+ca dbcommon
+ca dblinopt
+ca dbpencil
+ca info
+ca dbcolcom
+ca dbthin6d
+ei
+if bnlelens
+ca rhicelens
+ei

        ie=i
!++  For absorbed particles set all coordinates to zero. Also
!++  include very large offsets, let's say above 100mm or
!++  100mrad.
          do j = 1, napx
            if ( (part_abs_pos(j).ne.0 .and. part_abs_turn(j).ne.0) .or.&
     &xv(1,j).gt.100d0 .or.                                             &
     &yv(1,j).gt.100d0 .or.                                             &
     &xv(2,j).gt.100d0 .or.                                             &
     &yv(2,j).gt.100d0) then
              xv(1,j) = 0d0
              yv(1,j) = 0d0
              xv(2,j) = 0d0
              yv(2,j) = 0d0
              ejv(j)  = myenom
              sigmv(j)= 0d0
              part_abs_pos(j)=ie
              part_abs_turn(j)=iturn
              secondary(j) = 0
              tertiary(j)  = 0
              other(j) = 0
              nabs_type(j) = 0
            endif
          end do

!GRD SAVE COORDINATES OF PARTICLE 1 TO CHECK ORBIT
          if(firstrun) then
            xbob(ie)=xv(1,1)
            ybob(ie)=xv(2,1)
            xpbob(ie)=yv(1,1)
            ypbob(ie)=yv(2,1)
          endif

!++  Here comes sixtrack stuff
          if(ic(i).le.nblo) then
            do jb=1,mel(ic(i))
               myix=mtyp(ic(i),jb)
            enddo
          else
              myix=ic(i)-nblo
          endif

!++  Make sure we go into collimation routine for any definition
!++  of collimator element, relying on element name instead.
          if (                                                          &
!GRD HERE ARE SOME CHANGES TO MAKE RHIC TRAKING AVAILABLE
!APRIL2005
     &(bez(myix)(1:3).eq.'TCP'.or.bez(myix)(1:3).eq.'tcp') .or.         &
     &(bez(myix)(1:3).eq.'TCS'.or.bez(myix)(1:3).eq.'tcs') .or.         &
!UPGRADE January 2005
     &(bez(myix)(1:3).eq.'TCL'.or.bez(myix)(1:3).eq.'tcl') .or.         &
     &(bez(myix)(1:3).eq.'TCT'.or.bez(myix)(1:3).eq.'tct') .or.         &
     &(bez(myix)(1:3).eq.'TCD'.or.bez(myix)(1:3).eq.'tcd') .or.         &
     &(bez(myix)(1:3).eq.'TDI'.or.bez(myix)(1:3).eq.'tdi') .or.         &
! UPGRADE MAI 2006 -> TOTEM
     &(bez(myix)(1:3).eq.'TCX'.or.bez(myix)(1:3).eq.'tcx') .or.         &
! TW 04/2008 adding TCRYO 
     &(bez(myix)(1:3).eq.'TCR'.or.bez(myix)(1:3).eq.'tcr') .or.         &
!RHIC
     &(bez(myix)(1:3).eq.'COL'.or.bez(myix)(1:3).eq.'col') ) then
!GRD     write(*,*) bez(myix),'found!!'
!APRIL2005
         myktrack = 1
          else
            myktrack = ktrack(i)
          endif
!          write(*,*) 'ralph>  Element name: ', bez(myix), ktrack(i),
!     &                myktrack

      end

!>
!! collimate_end_element()
!! This routine is called at the end of every element
!<
      subroutine collimate_end_element
      implicit none

+ca crcoall
+if crlibm
+ca crlibco
+ei
      integer ix,j,jj,jx,kpz,kzz,napx0,nbeaux,nmz,nthinerr
      double precision benkcc,cbxb,cbzb,cikveb,crkveb,crxb,crzb,r0,r000,&
     &r0a,r2b,rb,rho2b,rkb,tkb,xbb,xrb,zbb,zrb
      logical lopen
+ca parpro
+ca parnum
+ca common
+ca commons
+ca commont1
+ca commondl
+ca commonxz
+ca commonta
+ca commonmn
+ca commonm1
+ca commontr
+ca beamdim
      dimension nbeaux(nbb)
+if collimat
+ca collpara
+ca dbtrthin
+ca database
+ca dbcommon
+ca dblinopt
+ca dbpencil
+ca info
+ca dbcolcom
+ca dbthin6d
+ei
+if bnlelens
+ca rhicelens
+ei

      if (firstrun) then
        if (rselect.gt.0 .and. rselect.lt.65) then
          do j = 1, napx
            xj     = (xv(1,j)-torbx(ie))/1d3
            xpj    = (yv(1,j)-torbxp(ie))/1d3
            yj     = (xv(2,j)-torby(ie))/1d3
            ypj    = (yv(2,j)-torbyp(ie))/1d3
            pj     = ejv(j)/1d3

            if (iturn.eq.1.and.j.eq.1) then
              sum_ax(ie)=0d0
              sum_ay(ie)=0d0
            endif

            if (tbetax(ie).gt.0.) then
              gammax = (1d0 + talphax(ie)**2)/tbetax(ie)
              gammay = (1d0 + talphay(ie)**2)/tbetay(ie)
            else
              gammax = (1d0 + talphax(ie-1)**2)/tbetax(ie-1)
              gammay = (1d0 + talphay(ie-1)**2)/tbetay(ie-1)
            endif

            if (part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
              if(tbetax(ie).gt.0.) then
                nspx    = sqrt(                                         &
     &               abs( gammax*(xj)**2 +                              &
     &               2d0*talphax(ie)*xj*xpj +                           &
     &               tbetax(ie)*xpj**2 )/myemitx0_collgap
     &               )
                nspy    = sqrt(                                         &
     &               abs( gammay*(yj)**2 +                              &
     &               2d0*talphay(ie)*yj*ypj +                           &
     &               tbetay(ie)*ypj**2 )/myemity0_collgap
     &               )
              else
                nspx    = sqrt(                                         &
     &               abs( gammax*(xj)**2 +                              &
     &               2d0*talphax(ie-1)*xj*xpj +                         &
     &               tbetax(ie-1)*xpj**2 )/myemitx0_collgap
     &               )
                nspy    = sqrt(                                         &
     &               abs( gammay*(yj)**2 +                              &
     &               2d0*talphay(ie-1)*yj*ypj +                         &
     &               tbetay(ie-1)*ypj**2 )/myemity0_collgap
     &               )
              endif

              sum_ax(ie)   = sum_ax(ie) + nspx
              sqsum_ax(ie) = sqsum_ax(ie) + nspx**2
              sum_ay(ie)   = sum_ay(ie) + nspy
              sqsum_ay(ie) = sqsum_ay(ie) + nspy**2
              nampl(ie)    = nampl(ie) + 1
            else
              nspx = 0d0
              nspy = 0d0
            endif
              sampl(ie)    = totals
              ename(ie)    = bez(myix)(1:16)
          end do
        endif
      endif

!GRD THIS LOOP MUST NOT BE WRITTEN INTO THE "IF(FIRSTRUN)" LOOP !!!!
          if (dowritetracks) then
            do j = 1, napx
              xj     = (xv(1,j)-torbx(ie))/1d3
              xpj    = (yv(1,j)-torbxp(ie))/1d3
              yj     = (xv(2,j)-torby(ie))/1d3
              ypj    = (yv(2,j)-torbyp(ie))/1d3

              arcdx = 2.5d0
              arcbetax = 180d0

                if (xj.le.0.) then
                  xdisp = xj + (pj-myenom)/myenom * arcdx               &
     &* sqrt(tbetax(ie)/arcbetax)
                else
                  xdisp = xj - (pj-myenom)/myenom * arcdx               &
     &* sqrt(tbetax(ie)/arcbetax)
                endif
                xndisp = xj
                nspxd   = sqrt(                                         &
     &abs(gammax*xdisp**2 + 2d0*talphax(ie)*xdisp*xpj                   &
     &+ tbetax(ie)*xpj**2)/myemitx0_collgap
     &)
                nspx    = sqrt(                                         &
     &abs( gammax*xndisp**2 + 2d0*talphax(ie)*xndisp*                   &
     &xpj + tbetax(ie)*xpj**2 )/myemitx0_collgap
     &)
                nspy    = sqrt(                                         &
     &abs( gammay*yj**2 + 2d0*talphay(ie)*yj                            &
     &*ypj + tbetay(ie)*ypj**2 )/myemity0_collgap
     &)

         if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
         if ((secondary(j).eq.1.or.tertiary(j).eq.2.or.other(j).eq.4)
     & .and.(xv(1,j).lt.99d0 .and. xv(2,j).lt.99d0) .and.
!GRD HERE WE APPLY THE SAME KIND OF CUT THAN THE SIGSECUT PARAMETER
     &(
     &((
     &(xv(1,j)*1d-3)**2
     &/
     &(tbetax(ie)*myemitx0_collgap)
     &).ge.dble(sigsecut2)).or.
     &((
     &(xv(2,j)*1d-3)**2
     &/
     &(tbetay(ie)*myemity0_collgap)
     &).ge.dble(sigsecut2)).or.
     &(((xv(1,j)*1d-3)**2/(tbetax(ie)*myemitx0_collgap))+
     &((xv(2,j)*1d-3)**2/(tbetay(ie)*myemity0_collgap))
     &.ge.sigsecut3)
     &) ) then
                xj     = (xv(1,j)-torbx(ie))/1d3
                xpj    = (yv(1,j)-torbxp(ie))/1d3
                yj     = (xv(2,j)-torby(ie))/1d3
                ypj    = (yv(2,j)-torbyp(ie))/1d3
+if hdf5
       hdfpid=ipart(j)+100*samplenumber
       hdfturn=iturn
       hdfs=sampl(ie)
       hdfx=xv(1,j)
       hdfxp=yv(1,j)
       hdfy=xv(2,j)
       hdfyp=yv(2,j)
       hdfdee=(ejv(j)-myenom)/myenom
       hdftyp=secondary(j)+tertiary(j)+other(j)
       call APPENDREADING(hdfpid,hdfturn,hdfs,hdfx,hdfxp,hdfy,hdfyp,    &
     &                    hdfdee,hdftyp)
+ei
+if .not.hdf5
       write(38,'(1x,i8,1x,i4,1x,f10.2,4(1x,e11.5),1x,e11.3,1x,i4)')
     &ipart(j)+100*samplenumber,iturn,sampl(ie),                        &
     &xv(1,j),yv(1,j),                                                  &
     &xv(2,j),yv(2,j),(ejv(j)-myenom)/myenom,                           &
     &secondary(j)+tertiary(j)+other(j)
+ei
              endif
         endif
            end do

!!JUNE2005 here I close the "if(dowritetracks)" outside of the firstrun flag
      endif
!JUNE2005
      end


!>
!! collimate_end_turn()
!! This routine is called at the end of every turn
!<
      subroutine collimate_end_turn
      implicit none

+ca crcoall
+if crlibm
+ca crlibco
+ei
      integer i,ix,j,jj,jx,kpz,kzz,napx0,nbeaux,nmz,nthinerr
      double precision benkcc,cbxb,cbzb,cikveb,crkveb,crxb,crzb,r0,r000,&
     &r0a,r2b,rb,rho2b,rkb,tkb,xbb,xrb,zbb,zrb
      logical lopen
+ca parpro
+ca parnum
+ca common
+ca commons
+ca commont1
+ca commondl
+ca commonxz
+ca commonta
+ca commonmn
+ca commonm1
+ca commontr
+ca beamdim
      dimension nbeaux(nbb)
+if collimat
+ca collpara
+ca dbtrthin
+ca database
+ca dbcommon
+ca dblinopt
+ca dbpencil
+ca info
+ca dbcolcom
+ca dbthin6d
+ei
+if bnlelens
+ca rhicelens
+ei

!__________________________________________________________________
!++  Now do analysis at selected elements...

!++  Save twiss functions of present element
        ax0  = talphax(ie)
        bx0  = tbetax(ie)
        mux0 = mux(ie)
        ay0  = talphay(ie)
        by0  = tbetay(ie)
        muy0 = muy(ie)

!GRD GET THE COORDINATES OF THE PARTICLES AT THE IEth ELEMENT:
        do j = 1,napx
              xgrd(j)  = xv(1,j)
              xpgrd(j) = yv(1,j)
              ygrd(j)  = xv(2,j)
              ypgrd(j) = yv(2,j)

              xineff(j)  = xv(1,j) - torbx (ie)
              xpineff(j) = yv(1,j) - torbxp(ie)
              yineff(j)  = xv(2,j) - torby (ie)
              ypineff(j) = yv(2,j) - torbyp(ie)

              pgrd(j)  = ejv(j)
!APRIL2005
              ejfvgrd(j) = ejfv(j)
              sigmvgrd(j) = sigmv(j)
              rvvgrd(j) = rvv(j)
              dpsvgrd(j) = dpsv(j)
              oidpsvgrd(j) = oidpsv(j)
              dpsv1grd(j) = dpsv1(j)
!APRIL2005
!GRD IMPORTANT: ALL PARTICLES ABSORBED ARE CONSIDERED TO BE LOST,
!GRD SO WE GIVE THEM A LARGE OFFSET
           if (part_abs_pos(j).ne.0 .and. part_abs_turn(j).ne.0) then
              xgrd(j)  = 99.5d0
              ygrd(j)  = 99.5d0
           endif
        end do

!++  For LAST ELEMENT in the ring calculate the number of surviving
!++  particles and save into file versus turn number
        if (ie.eq.iu) then
          nsurvive = 0
          do j = 1, napx
            if (xgrd(j).lt.99d0 .and. ygrd(j).lt.99d0) then
              nsurvive = nsurvive + 1
            endif
          end do

          write(44,'(2i7)') iturn, nsurvive
          if (iturn.eq.numl) then
            nsurvive_end = nsurvive_end + nsurvive
          endif
        endif

!=======================================================================
!++  Do collimation analysis at element 20 ("zero" turn) or LAST
!++  ring element.

!++  If selecting, look at number of scattered particles at selected
!++  collimator. For the "zero" turn consider the information at element
!++  20 (before collimation), otherwise take information at last ring
!++  element.
        if (do_coll .and.                                               &
     &(  (iturn.eq.1 .and. ie.eq.20) .or.                               &
     &(ie.eq.iu)  )    ) then

!++  Calculate gammas
!------------------------------------------------------------------------
          gammax = (1 + talphax(ie)**2)/tbetax(ie)
          gammay = (1 + talphay(ie)**2)/tbetay(ie)

!________________________________________________________________________
!++  Loop over all particles.
          do j = 1, napx
!
!------------------------------------------------------------------------
!++  Save initial distribution of particles that were scattered on
!++  the first turn at the selected primary collimator
!
!            IF (DOWRITE_DIST .AND. DO_SELECT .AND. ITURN.EQ.1 .AND.
!     &          PART_SELECT(j).EQ.1) THEN
!              WRITE(987,'(4(1X,E15.7))') X00(J), XP00(J),
!     &                                        Y00(J), YP00(J)
!            ENDIF
!------------------------------------------------------------------------
!++  Do the binning in amplitude, only considering particles that were
!++  not absorbed before.

            if (xgrd(j).lt.99d0 .and. ygrd(j) .lt.99d0 .and.            &
     &(part_select(j).eq.1 .or. ie.eq.20)) then

!++  Normalized amplitudes are calculated

!++  Allow to apply some dispersive offset. Take arc dispersion (2m) and
!++  normalize with arc beta_x function (180m).
              arcdx    = 2.5d0
              arcbetax = 180d0
              xdisp = abs(xgrd(j)*1d-3) +                               &
     &abs((pgrd(j)-myenom)/myenom)*arcdx                                &
     &* sqrt(tbetax(ie)/arcbetax)
              nspx    = sqrt(                                           &
     &abs(gammax*xdisp**2 +                                             &
     &2d0*talphax(ie)*xdisp*(xpgrd(j)*1d-3)+                            &
     &tbetax(ie)*(xpgrd(j)*1d-3)**2 )/myemitx0_collgap
     &)
              nspy    = sqrt(                                           &
     &abs( gammay*(ygrd(j)*1d-3)**2 +                                   &
     &2d0*talphay(ie)*(ygrd(j)*1d-3*ypgrd(j)*1d-3)                      &
     &+ tbetay(ie)*(ypgrd(j)*1d-3)**2 )/myemity0_collgap
     &)

!++  Populate the efficiency arrays at the end of each turn...
! Modified by M.Fiascaris, July 2016
              if (ie.eq.iu) then
                do ieff = 1, numeff
                  if (counted_r(j,ieff).eq.0 .and.                      &
     &sqrt(                                                             &
     &((xineff(j)*1d-3)**2 +                                            &
     & (talphax(ie)*xineff(j)*1d-3 + tbetax(ie)*xpineff(j)*1d-3)**2)    &
     &/                                                                 &
     &(tbetax(ie)*myemitx0_collgap)                                     &
     &+                                                                 &
     &((yineff(j)*1d-3)**2 +                                            &
     & (talphay(ie)*yineff(j)*1d-3 + tbetay(ie)*ypineff(j)*1d-3)**2)    &
     &/                                                                 &
     &(tbetay(ie)*myemity0_collgap)                                     &
     &).ge.rsig(ieff)) then
                    neff(ieff) = neff(ieff)+1d0
                    counted_r(j,ieff)=1
                  endif

!++ 2D eff
        do ieffdpop =1, numeffdpop
          if (counted2d(j,ieff,ieffdpop).eq.0 .and.
     &abs((ejv(j)-myenom)/myenom).ge.dpopbins(ieffdpop)) then     
            neff2d(ieff,ieffdpop) = neff2d(ieff,ieffdpop)+1d0
            counted2d(j,ieff,ieffdpop)=1
          endif
        end do

        if (counted_x(j,ieff).eq.0 .and.                                &
     &sqrt(                                                             &
     &((xineff(j)*1d-3)**2 +                                            &
     & (talphax(ie)*xineff(j)*1d-3 + tbetax(ie)*xpineff(j)*1d-3)**2)    &
     &/                                                                 &
     &(tbetax(ie)*myemitx0_collgap)                                     &
     &).ge.rsig(ieff)) then
                    neffx(ieff) = neffx(ieff) + 1d0
                    counted_x(j,ieff)=1
        endif

        if (counted_y(j,ieff).eq.0 .and.
     &sqrt(                                                             &
     &((yineff(j)*1d-3)**2 +                                            &
     & (talphay(ie)*yineff(j)*1d-3 + tbetay(ie)*ypineff(j)*1d-3)**2)    &
     &/                                                                 &
     &(tbetay(ie)*myemity0_collgap)                                     &
     &).ge.rsig(ieff)) then
          neffy(ieff) = neffy(ieff) + 1d0
          counted_y(j,ieff)=1
        endif
      end do !do ieff = 1, numeff

            do ieffdpop = 1, numeffdpop
              if (counteddpop(j,ieffdpop).eq.0) then
                dpopmin = 0d0
                mydpop = abs((ejv(j)-myenom)/myenom)
                if (ieffdpop.gt.1) dpopmin = dpopbins(ieffdpop-1)
                  dpopmax = dpopbins(ieffdpop)
                if (mydpop.ge.dpopmin .and. mydpop.lt.mydpop) then
                   npartdpop(ieffdpop)=npartdpop(ieffdpop)+1
                endif
              endif

             if (counteddpop(j,ieffdpop).eq.0 .and.
     &abs((ejv(j)-myenom)/myenom).ge.dpopbins(ieffdpop)) then
               neffdpop(ieffdpop) = neffdpop(ieffdpop)+1d0
               counteddpop(j,ieffdpop)=1
             endif
            end do
          endif

!++  Do an emittance drift
              driftx = driftsx*sqrt(tbetax(ie)*myemitx0_collgap)
              drifty = driftsy*sqrt(tbetay(ie)*myemity0_collgap)
              if (ie.eq.iu) then
                dnormx  = driftx / sqrt(tbetax(ie)*myemitx0_collgap)
                dnormy  = drifty / sqrt(tbetay(ie)*myemity0_collgap)
                xnorm  = (xgrd(j)*1d-3) /
     &                   sqrt(tbetax(ie)*myemitx0_collgap)
                xpnorm = (talphax(ie)*(xgrd(j)*1d-3)+                   &
     &tbetax(ie)*(xpgrd(j)*1d-3)) /                                     &
     &sqrt(tbetax(ie)*myemitx0_collgap)
+if crlibm
                xangle = atan2_rn(xnorm,xpnorm)
+ei
+if .not.crlibm
                xangle = atan2(xnorm,xpnorm)
+ei
+if crlibm
                xnorm  = xnorm  + dnormx*sin_rn(xangle)
+ei
+if .not.crlibm
                xnorm  = xnorm  + dnormx*sin(xangle)
+ei
+if crlibm
                xpnorm = xpnorm + dnormx*cos_rn(xangle)
+ei
+if .not.crlibm
                xpnorm = xpnorm + dnormx*cos(xangle)
+ei
                xgrd(j)   = 1000d0 *
     &                     (xnorm * sqrt(tbetax(ie)*myemitx0_collgap) )
                xpgrd(j)  = 1000d0 *
     &                     ( (xpnorm*sqrt(tbetax(ie)*myemitx0_collgap)
     &-talphax(ie)*xgrd(j)*1d-3)/tbetax(ie))

                ynorm  = (ygrd(j)*1d-3)
     &               / sqrt(tbetay(ie)*myemity0_collgap)
                ypnorm = (talphay(ie)*(ygrd(j)*1d-3)+                   &
     &tbetay(ie)*(ypgrd(j)*1d-3)) /                                     &
     &sqrt(tbetay(ie)*myemity0_collgap)
+if crlibm
                yangle = atan2_rn(ynorm,ypnorm)
+ei
+if .not.crlibm
                yangle = atan2(ynorm,ypnorm)
+ei
+if crlibm
                ynorm  = ynorm  + dnormy*sin_rn(yangle)
+ei
+if .not.crlibm
                ynorm  = ynorm  + dnormy*sin(yangle)
+ei
+if crlibm
                ypnorm = ypnorm + dnormy*cos_rn(yangle)
+ei
+if .not.crlibm
                ypnorm = ypnorm + dnormy*cos(yangle)
+ei
                ygrd(j)   = 1000d0 *
     &                     (ynorm * sqrt(tbetay(ie)*myemity0_collgap) )
                ypgrd(j)  = 1000d0 * 
     &                     ( (ypnorm*sqrt(tbetay(ie)*myemity0_collgap)
     &-talphay(ie)*ygrd(j)*1d-3)/tbetay(ie))

                endif

!------------------------------------------------------------------------
!++  End of check for selection flag and absorption
            endif

!++  End of do loop over particles
          end do

!_________________________________________________________________
!
!++  End of collimation efficiency analysis for selected particles
        end if
!------------------------------------------------------------------
!++  For LAST ELEMENT in the ring compact the arrays by moving all
!++  lost particles to the end of the array.
        if (ie.eq.iu) then
          imov = 0
          do j = 1, napx
            if (xgrd(j).lt.99d0 .and. ygrd(j).lt.99d0) then
              imov = imov + 1
              xgrd(imov)           = xgrd(j)
              ygrd(imov)           = ygrd(j)
              xpgrd(imov)          = xpgrd(j)
              ypgrd(imov)          = ypgrd(j)
              pgrd(imov)           = pgrd(j)
              ejfvgrd(imov)        = ejfvgrd(j)
              sigmvgrd(imov)       = sigmvgrd(j)
              rvvgrd(imov)         = rvvgrd(j)
              dpsvgrd(imov)        = dpsvgrd(j)
              oidpsvgrd(imov)      = oidpsvgrd(j)
              dpsv1grd(imov)       = dpsv1grd(j)
              part_hit_pos(imov)   = part_hit_pos(j)
              part_hit_turn(imov)  = part_hit_turn(j)
              part_abs_pos(imov)   = part_abs_pos(j)
              part_abs_turn(imov)  = part_abs_turn(j)
              part_select(imov) = part_select(j)
              part_impact(imov) = part_impact(j)
              part_indiv(imov)  = part_indiv(j)
              part_linteract(imov)  = part_linteract(j)
              part_hit_before_pos(imov)  = part_hit_before_pos(j)
              part_hit_before_turn(imov) = part_hit_before_turn(j)
              secondary(imov) = secondary(j)
              tertiary(imov) = tertiary(j)
              other(imov) = other(j)
              nabs_type(imov) = nabs_type(j)
!GRD HERE WE ADD A MARKER FOR THE PARTICLE FORMER NAME
              ipart(imov) = ipart(j)
              flukaname(imov) = flukaname(j)
!KNS: Also compact nlostp (used for standard LOST calculations + output)
              nlostp(imov) = nlostp(j)
              do ieff = 1, numeff
                counted_r(imov,ieff) = counted_r(j,ieff)
                counted_x(imov,ieff) = counted_x(j,ieff)
                counted_y(imov,ieff) = counted_y(j,ieff)
              end do
            endif
          end do
          write(lout,*) 'INFO>  Compacted the particle distributions: ',
     &napx, ' -->  ', imov, ", turn =",iturn
          napx = imov
        endif

!------------------------------------------------------------------------
!++  Write final distribution
      if (dowrite_dist.and.(ie.eq.iu).and.(iturn.eq.numl)) then
        open(unit=9998, file='distn.dat')
        write(9998,*)
     &'# 1=x 2=xp 3=y 4=yp 5=z 6 =E'

        do j = 1, napx
          write(9998,'(6(1X,E15.7))') (xgrd(j)-torbx(1))/1d3,           &
     &(xpgrd(j)-torbxp(1))/1d3, (ygrd(j)-torby(1))/1d3,                 &
     &(ypgrd(j)-torbyp(1))/1d3,sigmvgrd(j),ejfvgrd(j)
        end do

        close(9998)
      endif

!GRD NOW ONE HAS TO COPY BACK THE NEW DISTRIBUTION TO ITS "ORIGINAL NAME"
!GRD AT THE END OF EACH TURN
      if (ie.eq.iu) then
         do j = 1,napx
            xv(1,j) = xgrd(j)
            yv(1,j) = xpgrd(j)
            xv(2,j) = ygrd(j)
            yv(2,j) = ypgrd(j)
            ejv(j)  = pgrd(j)
            ejfv(j)   = ejfvgrd(j)
            sigmv(j)  = sigmvgrd(j)
            rvv(j)    = rvvgrd(j)
            dpsv(j)   = dpsvgrd(j)
            oidpsv(j) = oidpsvgrd(j)
            dpsv1(j)  = dpsv1grd(j)
         end do
      endif

         if (firstrun) then
       if (rselect.gt.0 .and. rselect.lt.65) then
            do j = 1, napx
              xj     = (xv(1,j)-torbx(ie))/1d3
              xpj    = (yv(1,j)-torbxp(ie))/1d3
              yj     = (xv(2,j)-torby(ie))/1d3
              ypj    = (yv(2,j)-torbyp(ie))/1d3
              pj     = ejv(j)/1d3
              if (iturn.eq.1.and.j.eq.1) then
                sum_ax(ie)=0d0
                sum_ay(ie)=0d0
              endif

              if (tbetax(ie).gt.0.) then
          gammax = (1d0 + talphax(ie)**2)/tbetax(ie)
                gammay = (1d0 + talphay(ie)**2)/tbetay(ie)
              else
          gammax = (1d0 + talphax(ie-1)**2)/tbetax(ie-1)
          gammay = (1d0 + talphay(ie-1)**2)/tbetay(ie-1)
              endif
!
              if (part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
                if(tbetax(ie).gt.0.) then
          nspx    = sqrt(                                               &
     &abs( gammax*(xj)**2 +                                             &
     &2d0*talphax(ie)*xj*xpj +                                          &
     &tbetax(ie)*xpj**2 )/myemitx0_collgap
     &)
                nspy    = sqrt(                                         &
     &abs( gammay*(yj)**2 +                                             &
     &2d0*talphay(ie)*yj*ypj +                                          &
     &tbetay(ie)*ypj**2 )/myemity0_collgap
     &)
                else
          nspx    = sqrt(                                               &
     &abs( gammax*(xj)**2 +                                             &
     &2d0*talphax(ie-1)*xj*xpj +                                        &
     &tbetax(ie-1)*xpj**2 )/myemitx0_collgap
     &)
                nspy    = sqrt(                                         &
     &abs( gammay*(yj)**2 +                                             &
     &2d0*talphay(ie-1)*yj*ypj +                                        &
     &tbetay(ie-1)*ypj**2 )/myemity0_collgap
     &)
                endif
                sum_ax(ie)   = sum_ax(ie) + nspx
                sqsum_ax(ie) = sqsum_ax(ie) + nspx**2
                sum_ay(ie)   = sum_ay(ie) + nspy
                sqsum_ay(ie) = sqsum_ay(ie) + nspy**2
                nampl(ie)    = nampl(ie) + 1
                sampl(ie)    = totals
                ename(ie)    = bez(myix)(1:16)
              else
                nspx = 0d0
                nspy = 0d0
              endif
            end do
          endif
         endif

!GRD THIS LOOP MUST NOT BE WRITTEN INTO THE "IF(FIRSTRUN)" LOOP !!!!
          if (dowritetracks) then
            do j = 1, napx
              xj     = (xv(1,j)-torbx(ie))/1d3
              xpj    = (yv(1,j)-torbxp(ie))/1d3
              yj     = (xv(2,j)-torby(ie))/1d3
              ypj    = (yv(2,j)-torbyp(ie))/1d3
              arcdx = 2.5d0
              arcbetax = 180d0

                if (xj.le.0.) then
                  xdisp = xj + (pj-myenom)/myenom * arcdx               &
     &* sqrt(tbetax(ie)/arcbetax)
                else
                  xdisp = xj - (pj-myenom)/myenom * arcdx               &
     &* sqrt(tbetax(ie)/arcbetax)
                endif

                xndisp = xj
                nspxd   = sqrt(                                         &
     &abs(gammax*xdisp**2 + 2d0*talphax(ie)*xdisp*xpj                   &
     &+ tbetax(ie)*xpj**2)/myemitx0_collgap
     &)
                nspx    = sqrt(                                         &
     &abs( gammax*xndisp**2 + 2d0*talphax(ie)*xndisp*                   &
     &xpj + tbetax(ie)*xpj**2 )/myemitx0_collgap
     &)
                nspy    = sqrt(                                         &
     &abs( gammay*yj**2 + 2d0*talphay(ie)*yj                            &
     &*ypj + tbetay(ie)*ypj**2 )/myemity0_collgap
     &)

         if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
        if ((secondary(j).eq.1.or.tertiary(j).eq.2.or.other(j).eq.4)
     &.and.(xv(1,j).lt.99d0 .and. xv(2,j).lt.99d0) .and.
!GRD HERE WE APPLY THE SAME KIND OF CUT THAN THE SIGSECUT PARAMETER
     &(
     &((
     &(xv(1,j)*1d-3)**2
     &/
     &(tbetax(ie)*myemitx0_collgap)
     &).ge.dble(sigsecut2)).or.
     &((
     &(xv(2,j)*1d-3)**2
     &/
     &(tbetay(ie)*myemity0_collgap)
     &).ge.dble(sigsecut2)).or.
     &(((xv(1,j)*1d-3)**2/(tbetax(ie)*myemitx0_collgap))+
     &((xv(2,j)*1d-3)**2/(tbetay(ie)*myemity0_collgap))
     &.ge.sigsecut3)
     &) ) then
                xj     = (xv(1,j)-torbx(ie))/1d3
                xpj    = (yv(1,j)-torbxp(ie))/1d3
                yj     = (xv(2,j)-torby(ie))/1d3
                ypj    = (yv(2,j)-torbyp(ie))/1d3
+if hdf5
       hdfpid=ipart(j)+100*samplenumber
       hdfturn=iturn
       hdfs=sampl(ie)
       hdfx=xv(1,j)
       hdfxp=yv(1,j)
       hdfy=xv(2,j)
       hdfyp=yv(2,j)
       hdfdee=(ejv(j)-myenom)/myenom
       hdftyp=secondary(j)+tertiary(j)+other(j)
       call APPENDREADING(hdfpid,hdfturn,hdfs,hdfx,hdfxp,hdfy,hdfyp,    &
     &                    hdfdee,hdftyp)
+ei
+if .not.hdf5
       write(38,'(1x,i8,1x,i4,1x,f10.2,4(1x,e11.5),1x,e11.3,1x,i4)')
     &ipart(j)+100*samplenumber,iturn,sampl(ie),                        &
     &xv(1,j),yv(1,j),                                                  &
     &xv(2,j),yv(2,j),(ejv(j)-myenom)/myenom,                           &
     &secondary(j)+tertiary(j)+other(j)
+ei
                endif
              endif
            end do
          endif
!=======================================================================
      end




!>
!! "Merlin" scattering collimation configuration
!! This routine pre-calcuates some varibles for
!! the nuclear properties
!<
      subroutine collimate_init_merlin()
      implicit none
      integer i
      double precision CalcElectronDensity,CalcPlasmaEnergy
+ca interac
      ! compute the electron densnity and plasma energy for each material
      do i=1, nmat
         edens(i) = CalcElectronDensity(zatom(i),rho(i),anuc(i))
         pleng(i) = CalcPlasmaEnergy(edens(i))
      end do
      end

!>
!! K2 scattering collimation configuration
!<
      subroutine collimate_init_k2()
!nothing currently
      end

!>
!! subroutine collimate2(c_material, c_length, c_rotation,           &
!!-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----
!!----                                                                    -----
!!-----  NEW ROUTINES PROVIDED FOR THE COLLIMATION STUDIES VIA SIXTRACK   -----
!!-----                                                                   -----
!!-----          G. ROBERT-DEMOLAIZE, November 1st, 2004                  -----
!!-----                                                                   -----
!!-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----
!!++  Based on routines by JBJ. Changed by RA 2001.
!!GRD
!!GRD MODIFIED VERSION FOR COLLIMATION SYSTEM: G. ROBERT-DEMOLAIZE
!!GRD
!!
!!++  - Deleted all HBOOK stuff.
!!++  - Deleted optics routine and all parser routines.
!!++  - Replaced RANMAR call by RANLUX call
!!++  - Included RANLUX code from CERNLIB into source
!!++  - Changed dimensions from CGen(100,nmat) to CGen(200,nmat)
!!++  - Replaced FUNPRE with FUNLXP
!!++  - Replaced FUNRAN with FUNLUX
!!++  - Included all CERNLIB code into source: RANLUX, FUNLXP, FUNLUX,
!!++                                         FUNPCT, FUNLZ, RADAPT,
!!++                                           RGS56P
!!++    with additional entries:             RLUXIN, RLUXUT, RLUXAT,
!!++                                           RLUXGO
!!++
!!++  - Changed program so that Nev is total number of particles
!!++    (scattered and not-scattered)
!!++  - Added debug comments
!!++  - Put real dp/dx
!<
      subroutine collimate2(c_material, c_length, c_rotation,           &
     &c_aperture, c_offset, c_tilt,x_in, xp_in, y_in,yp_in,p_in, s_in,  &
     &     np, enom,                                                    &
     &     lhit_pos, lhit_turn,                                         &
     &     part_abs_pos_local, part_abs_turn_local,                     &
     &     impact, indiv, lint, onesided,name,                          &
     &     flagsec, j_slices, nabs_type)
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei

+ca commonex
+ca parpro
+ca collpara
+ca dbpencil
+ca info
+ca dbcollim
+ca flukavars

+ca database
+ca collMatNum

      double precision x_flk,xp_flk,y_flk,yp_flk,zpj

      double precision s_impact
      integer flagsec(npart)

!     SR, 18-08-2005: add temporary variable to write in FirstImpacts
!     the initial distribution of the impacting particles in the
!     collimator frame.
      double precision xinn,xpinn,yinn,ypinn

!     SR, 29-08-2005: add the slice number to calculate the impact
!     location within the collimator.
!     j_slices = 1 for the a non sliced collimator!
      integer j_slices

      save
!      write(lout,*) 'In col2 ', c_material, c_length, c_aperture,       &
!     &c_offset, c_tilt, x_in, xp_in, y_in,p_in, np, enom
!=======================================================================
! Be=1 Al=2 Cu=3 W=4 Pb=5
! LHC uses:    Al, 0.2 m
!              Cu, 1.0 m

      if (c_material.eq.'BE') then
         mat = 1
      elseif (c_material.eq.'AL') then
         mat = 2
      elseif (c_material.eq.'CU') then
         mat = 3
      elseif (c_material.eq.'W') then
         mat = 4
      elseif (c_material.eq.'PB') then
         mat = 5
      elseif (c_material.eq.'C') then
         mat = 6
      elseif (c_material.eq.'C2') then
         mat = 7
      elseif (c_material.eq.'MoGR') then
         mat = 8
      elseif (c_material.eq.'CuCD') then
         mat = 9
      elseif (c_material.eq.'Mo') then
         mat = 10
      elseif (c_material.eq.'Glid') then
         mat = 11
      elseif (c_material.eq.'Iner') then
         mat = 12
!02/2008 TW added vacuum and black absorber (was missing)
      elseif (c_material.eq.'VA') then
         mat = nmat-1
      elseif (c_material.eq.'BL') then
         mat = nmat
      else
         write(lout,*)
         write(lout,*) 'ERR>  In subroutine collimate2:'
         write(lout,*) 'ERR>  Material "', c_material, '" not found.'
         write(lout,*) 'ERR>  Check your CollDB! Stopping now.'
         call prror(-1)
      endif
!
        length  = c_length
        nev = np
        p0  = enom

!++  Initialize scattering processes
      call scatin(p0)

! EVENT LOOP,  initial distribution is here a flat distribution with
! xmin=x-, xmax=x+, etc. from the input file

      nhit    = 0
      fracab  = 0d0
      mirror  = 1d0

!==> SLICE here

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do j = 1, nev
!
! SR-GRD (04-08-2005):
!        Don't do scattering process for particles already absorbed
         if (       part_abs_pos_local(j)  .ne. 0
     &        .and. part_abs_turn_local(j) .ne. 0)
     &        goto 777

        impact(j) = -1d0
        lint(j)   = -1d0
        indiv(j)  = -1d0

        x   = x_in(j)
        xp  = xp_in(j)
        z   = y_in(j)
        zp  = yp_in(j)
        p   = p_in(j)
!        sp  = s_in(J)
        sp   = 0d0
        dpop = (p - p0)/p0
!APRIL2005
        x_flk  = 0d0
        y_flk  = 0d0
        xp_flk = 0d0
        yp_flk = 0d0
!APRIL2005
!
!++  Transform particle coordinates to get into collimator coordinate
!++  system
!
!++  First check whether particle was lost before
!
!        if (x.lt.99d-3 .and. z.lt.99d-3) then
!
!++  First do rotation into collimator frame
!
+if crlibm
          x  = x_in(j)*cos_rn(c_rotation) +sin_rn(c_rotation)*y_in(j)
+ei
+if .not.crlibm
          x  = x_in(j)*cos(c_rotation) +sin(c_rotation)*y_in(j)
+ei
+if crlibm
          z  = y_in(j)*cos_rn(c_rotation) -sin_rn(c_rotation)*x_in(j)
+ei
+if .not.crlibm
          z  = y_in(j)*cos(c_rotation) -sin(c_rotation)*x_in(j)
+ei
+if crlibm
          xp = xp_in(j)*cos_rn(c_rotation)+sin_rn(c_rotation)*yp_in(j)
+ei
+if .not.crlibm
          xp = xp_in(j)*cos(c_rotation)+sin(c_rotation)*yp_in(j)
+ei
+if crlibm
          zp = yp_in(j)*cos_rn(c_rotation)-sin_rn(c_rotation)*xp_in(j)
+ei
+if .not.crlibm
          zp = yp_in(j)*cos(c_rotation)-sin(c_rotation)*xp_in(j)
+ei
!
!++  For one-sided collimators consider only positive X. For negative
!++  X jump to the next particle
!
          if ((onesided .and. x.lt.0d0).and.
     &         ((icoll.ne.ipencil) .or. (iturn.ne.1))) goto 777 ! RB: adding exception from goto if it's 
!
!++  Now mirror at the horizontal axis for negative X offset
!
          if (x.lt.0d0) then                                             !hr09
            mirror = -1d0
            tiltangle = -1d0*c_tilt(2)
          endif
          if (x.ge.0d0) then                                             !hr09
            mirror = 1d0
            tiltangle = c_tilt(1)
          endif
          x  = mirror * x
          xp = mirror * xp
!
!          if (j.eq.1) then
!             write(*,*) 'INFOtilt',
!     &            icoll, j_slices, c_tilt(1), c_tilt(2),
!     &            mirror, tiltangle, c_offset, c_aperture/2
!          endif
!
!++  Shift with opening and offset
!
          x  = (x - c_aperture/2d0) - mirror*c_offset                    !hr09
!
!++  Include collimator tilt
!
          if (tiltangle.gt.0.d0) then                                    !hr09
            xp = xp - tiltangle
          endif
          if (tiltangle.lt.0.) then
+if crlibm
            x  = x + sin_rn(tiltangle) * c_length
+ei
+if .not.crlibm
            x  = x + sin(tiltangle) * c_length
+ei
            xp = xp - tiltangle
          endif
!
!++  For selected collimator, first turn reset particle distribution
!++  to simple pencil beam
!
! -- TW why did I set this to 0, seems to be needed for getting 
!       right amplitude => no "tilt" of jaw for the first turn !!!!
!          c_tilt(1) = 0d0
!          c_tilt(2) = 0d0
!
          nprim = 3

          if (( (icoll.eq.ipencil .and. iturn.eq.1) .or. (iturn.eq.1     &
     &.and. ipencil.eq.999 .and. icoll.le.nprim .and.                   &
     &(j.ge.(icoll-1)*nev/nprim) .and. (j.le.(icoll)*nev/nprim))).and.
     &         (pencil_distr.ne.3)) then   ! RB addition : don't go in this if-statement if pencil_distr=3. This distribution is generated in main loop instead
! -- TW why did I set this to 0, seems to be needed for getting 
!       right amplitude => no "tilt" of jaw for the first turn !!!!
          c_tilt(1) = 0d0
          c_tilt(2) = 0d0
!
!
!AUGUST2006: Standard pencil beam as implemented by GRD ------- TW
!
             if (pencil_rmsx.eq.0.d0 .and. pencil_rmsy.eq.0.d0) then     !hr09
                x    = pencil_dx(icoll)
                xp   = 0.
                z    = 0.
                zp   = 0.
!                dpop = 0.
             endif
!
!AUGUST2006: Rectangular (pencil-beam) sheet-beam with  ------ TW
!            pencil_offset is the rectangulars center
!            pencil_rmsx defines spread of impact parameter
!            pencil_rmsy defines spread parallel to jaw surface
! 
            if (pencil_distr.eq.0 .and.(pencil_rmsx.ne.0.               &
     &.or.pencil_rmsy.ne.0.)) then
! how to assure that all generated particles are on the jaw ?!
                x    = pencil_dx(icoll)                                 &
     &                 + pencil_rmsx*(dble(rndm4())-0.5d0)               !hr09
                xp   = 0.
                z    = pencil_rmsy*(dble(rndm4())-0.5d0)                 !hr09
                zp   = 0.
!                dpop = 0.
             endif
!
!AUGUST2006: Gaussian (pencil-beam) sheet-beam with ------- TW
!            pencil_offset is the mean  gaussian distribution
!            pencil_rmsx defines spread of impact parameter
!            pencil_rmsy defines spread parallel to jaw surface
! 
            if (pencil_distr.eq.1 .and.(pencil_rmsx.ne.0.d0             &!hr09
     &.or.pencil_rmsy.ne.0.d0 )) then                                    !hr09
                x    = pencil_dx(icoll) + pencil_rmsx*ran_gauss(2d0)
! all generated particles are on the jaw now
                x    = sqrt(x**2)
                xp   = 0.
                z    = pencil_rmsy*ran_gauss(2d0)
                zp   = 0.
!                dpop = 0.                
             endif
!
!AUGUST2006: Gaussian (pencil-beam) sheet-beam with ------- TW
!            pencil_offset is the mean  gaussian distribution
!            pencil_rmsx defines spread of impact parameter
!                        here pencil_rmsx is not gaussian!!!
!            pencil_rmsy defines spread parallel to jaw surface
! 
            if (pencil_distr.eq.2 .and.(pencil_rmsx.ne.0.d0             &!hr09
     &.or.pencil_rmsy.ne.0.d0 )) then                                    !hr09
                x    = pencil_dx(icoll)                                 &!hr09
     &              + pencil_rmsx*(dble(rndm4())-0.5d0)                  !hr09
! all generated particles are on the jaw now
                x    = sqrt(x**2)
                xp   = 0.
                z    = pencil_rmsy*ran_gauss(2d0)
                zp   = 0.
!                dpop = 0.                
             endif

!
!JULY2007: Selection of pos./neg. jaw  implemented by GRD ---- TW
!
! ensure that for onesided only particles on pos. jaw are created
             if (onesided) then
                mirror = 1d0
             else
!     if(rndm4().lt.0.5) mirror = -1d0
!     if(rndm4().ge.0.5) mirror = 1d0  => using two different random
                if(rndm4().lt.0.5) then 
                   mirror = -1d0
                else 
                   mirror = 1d0
                endif
             endif
!    
! -- TW SEP07 if c_tilt is set to zero before entering pencil beam 
!             section the assigning of the tilt will result in 
!             assigning zeros 
             if (mirror.lt.0d0) then                                     !hr09
!!     tiltangle = -1d0*c_tilt(2)
                tiltangle = c_tilt(2)
             else 
                tiltangle = c_tilt(1)
             endif
!!!!--- commented this out since particle is tilted after leaving 
!!!!--- collimator -> remove this  code fragment in final verion
!!             x  = mirror * x
!!             xp = mirror * xp
!
!++  Include collimator tilt
! this is propably not correct
!
!             xp =  (xp_pencil0*cos(c_rotation)+                         &
!     &            sin(c_rotation)*yp_pencil0)
!             if (tiltangle.gt.0.) then
!                xp = xp - tiltangle
!!             endif
!!             elseif (tiltangle.lt.0.) then
!             else
!               x  = x + sin(tiltangle) * c_length
!               xp = xp - tiltangle
!             endif
!
       write(9997,'(f10.8,(2x,f10.8),(2x,f10.8),(2x,f10.8)(2x,f10.8))')     
     &            x, xp, z, zp, tiltangle

!
      endif


!          if(rndm4().lt.0.5) mirror = -abs(mirror)
!          if(rndm4().ge.0.5) mirror = abs(mirror)
!        endif
!
!     SR, 18-08-2005: after finishing the coordinate transformation,
!     or the coordinate manipulations in case of pencil beams,
!     write down the initial coordinates of the impacting particles
          xinn  = x
          xpinn = xp
          yinn  = z
          ypinn = zp
!
!++  Possibility to slice here (RA,SR: 29-08-2005)
!
!
!++  particle passing above the jaw are discarded => take new event
!++  entering by the face, shorten the length (zlm) and keep track of
!++  entrance longitudinal coordinate (keeps) for histograms
!
!++  The definition is that the collimator jaw is at x>=0.
!
!++  1) Check whether particle hits the collimator
!
          hit     =  .false.
          s       =  0.d0                                                !hr09
          keeps   =  0.d0                                                !hr09
          zlm     =  -1d0 * length
!
          if (x.ge.0.d0) then                                            !hr09
!
!++  Particle hits collimator and we assume interaction length ZLM equal
!++  to collimator length (what if it would leave collimator after
!++  small length due to angle???)
!
            zlm = length
            impact(j) = x
            indiv(j) = xp
          else if (xp.le.0.d0) then                                      !hr09
!
!++  Particle does not hit collimator. Interaction length ZLM is zero.
!
            zlm = 0d0
          else
!
!++  Calculate s-coordinate of interaction point
!
            s = (-1d0*x) / xp
            if (s.le.0) then
              write(lout,*) 'S.LE.0 -> This should not happen'
              call prror(-1)
            endif
!
            if (s .lt. length) then
              zlm = length - s
              impact(j) = 0d0
              indiv(j) = xp
            else
              zlm = 0d0
            endif
!
          endif
!
!++  First do the drift part
! DRIFT PART
          drift_length = length - zlm
          if (drift_length.gt.0.d0) then                                 !hr09
            if(iexact.eq.0) then
              x  = x + xp* drift_length
              z  = z + zp * drift_length
              sp = sp + drift_length
            else
              zpj = sqrt(1d0-xp**2-zp**2)
              x = x + drift_length*(xp/zpj)
              z = z + drift_length*(zp/zpj)
              sp = sp + drift_length
            endif
          endif
!
!++  Now do the scattering part
!
          if (zlm.gt.0.) then
!JUNE2005
            s_impact = sp
!JUNE2005
            nhit = nhit + 1
!            WRITE(*,*) J,X,XP,Z,ZP,SP,DPOP
!     RB: add new input arguments to jaw icoll,iturn,ipart for writeout
            call jaw(s, nabs, icoll,iturn,name(j),dowrite_impact)

            nabs_type(j) = nabs
!JUNE2005
!JUNE2005 SR+GRD: CREATE A FILE TO CHECK THE VALUES OF IMPACT PARAMETERS
!JUNE2005
!     SR, 29-08-2005: Add to the longitudinal coordinates the position
!     of the slice beginning

            if(dowrite_impact) then
              if(flagsec(j).eq.0) then
               write(39,'(i5,1x,i7,1x,i2,1x,i1,2(1x,f5.3),8(1x,e17.9))')&
     &               name(j),iturn,icoll,nabs,                          &
     &               s_impact + (dble(j_slices)-1d0) * c_length,        &!hr09
     &               s+sp + (dble(j_slices)-1d0) * c_length,            &!hr09
     &               xinn,xpinn,yinn,ypinn,                             &
     &               x,xp,z,zp
              endif
            endif
!!     SR, 18-08-2005: add also the initial coordinates of the
!!                     impacting particles!
!            if(flagsec(j).eq.0) then
!              write(333,'(i5,1x,i7,1x,i2,1x,i1,2(1x,f5.3),8(1x,e17.9))')&
!     +              name(j),iturn,icoll,nabs,s_impact,s+sp,
!     +              xinn,xpinn,yinn,ypinn,
!     +              x,xp,z,zp
!            endif
!     !Old format...
!            if(flagsec(j).eq.0) then
!              write(333,'(i5,1x,i4,1x,i2,1x,i1,2(1x,f5.3),2(1x,e16.7))')
!     &name(j),iturn,icoll,nabs,s_impact,s+sp,impact(j),x
!            endif
!JUNE2005
            
            lhit_pos(j) = ie
            lhit_turn(j) = iturn
            
!-- September2006  TW added from Ralphs code
!--------------------------------------------------------------
!++ Change tilt for pencil beam impact
!
!            if ( (icoll.eq.ipencil                                      &
!     &           .and. iturn.eq.1)   .or.                               &
!     &           (iturn.eq.1 .and. ipencil.eq.999 .and.                 &
!     &                             icoll.le.nprim .and.                 &
!     &            (j.ge.(icoll-1)*nev/nprim) .and.                      &
!     &            (j.le.(icoll)*nev/nprim)                              &
!     &           )  ) then
!
!               if (.not. changed_tilt1(icoll) .and. mirror.gt.0.) then 
! ----- Maybe a warning would be nice that c_tilt is overwritten !!!!!
! changed xp_pencil0(icoll) to xp_pencil0 due to definition mismatch
! this has to be solved if necassary and understood 
!                 c_tilt(1) = xp_pencil0(icoll)*cos(c_rotation)+         &
!     &                       sin(c_rotation)*yp_pencil0(icoll)
!                 c_tilt(1) = xp_pencil0*cos(c_rotation)+                &
!     &                       sin(c_rotation)*yp_pencil0
!                 write(*,*) "INFO> Changed tilt1  ICOLL  to  ANGLE  ",  &
!     &                   icoll, c_tilt(1), j
!                 changed_tilt1(icoll) = .true.
!               elseif (.not. changed_tilt2(icoll)                       &
!     &                                   .and. mirror.lt.0.) then
! changed xp_pencil0(icoll) to xp_pencil0 due to definition mismatch
! this has to be solved if necassary and understood 
!                 c_tilt(2) = -1.*(xp_pencil0(icoll)*cos(c_rotation)+    &
!     &                       sin(c_rotation)*yp_pencil0(icoll))
!                 c_tilt(2) = -1.*(xp_pencil0*cos(c_rotation)+           &
!     &                       sin(c_rotation)*yp_pencil0)
!                 write(*,*) "INFO> Changed tilt2  ICOLL  to  ANGLE  ",  &
!     &                   icoll, c_tilt(2), j
!                 changed_tilt2(icoll) = .true.
!               endif
!            endif
!
!----------------------------------------------------------------
!-- September 2006
!
!++  If particle is absorbed then set x and y to 99.99 mm
!     SR: before assigning new (x,y) for nabs=1, write the
!     inelastic impact file .

!     RB: writeout should be done for both inelastic and single diffractive. doing all transformations in x_flk and making the set to 99.99 mm conditional for nabs=1
!!! /* start RB fix */

! transform back to lab system for writeout. 
! keep x,y,xp,yp unchanged for continued tracking, store lab system variables in x_flk etc

            x_flk = xInt
            xp_flk = xpInt


+if crlibm
            if (tiltangle.gt.0.) then
               x_flk  = x_flk  + tiltangle*(sInt+sp)
               xp_flk = xp_flk + tiltangle
            elseif (tiltangle.lt.0.d0) then !hr09
               xp_flk = xp_flk + tiltangle
               x_flk  = x_flk - sin_rn(tiltangle) * ( length-(sInt+sp) )
            endif
            x_flk = (x_flk + c_aperture/2d0) + mirror*c_offset !hr09
            x_flk    = mirror * x_flk
            xp_flk   = mirror * xp_flk
            y_flk  = yInt  *cos_rn(-1d0*c_rotation) -                         &          
     &           x_flk  *sin_rn(-1d0*c_rotation)
            yp_flk = ypInt *cos_rn(-1d0*c_rotation) -                         &
     &           xp_flk *sin_rn(-1d0*c_rotation)
            x_flk  = x_flk  *cos_rn(-1d0*c_rotation) +                         &
     &           yInt  *sin_rn(-1d0*c_rotation)
            xp_flk = xp_flk *cos_rn(-1d0*c_rotation) +                         &
     &           ypInt *sin_rn(-1d0*c_rotation)

+ei
+if .not.crlibm
            if (tiltangle.gt.0.) then
               x_flk  = x_flk  + tiltangle*(sInt+sp)
               xp_flk = xp_flk + tiltangle
            elseif (tiltangle.lt.0.d0) then !hr09
               xp_flk = xp_flk + tiltangle
               x_flk  = x_flk - sin(tiltangle) * ( length -(sInt+sp) )
            endif
            x_flk = (x_flk + c_aperture/2d0) + mirror*c_offset !hr09
            x_flk    = mirror * x_flk
            xp_flk   = mirror * xp_flk
            y_flk  = yInt  *cos(-1d0*c_rotation) -                         &          
     &           x_flk  *sin(-1d0*c_rotation)
            yp_flk = ypInt *cos(-1d0*c_rotation) -                         &
     &           xp_flk *sin(-1d0*c_rotation)
            x_flk  = x_flk  *cos(-1d0*c_rotation) +                         &
     &           yInt  *sin(-1d0*c_rotation)
            xp_flk = xp_flk *cos(-1d0*c_rotation) +                         &
     &           ypInt *sin(-1d0*c_rotation)

+ei

! write out all impacts to all_impacts.dat
            if(dowrite_impact) then
         write(4801,'(i4,(1x,f6.3),(1x,f8.6),4(1x,e19.10),i2,2(1x,i7))')     &  
     &              icoll,c_rotation,                                        &  
     &              sInt + sp + (dble(j_slices)-1d0) * c_length,                         & !hr09
     &              x_flk*1d3, xp_flk*1d3, y_flk*1d3, yp_flk*1d3,                     &
     &              nabs,name(j),iturn
            endif

! standard FLUKA_impacts writeout of inelastic and single diffractive
            if ((nabs.eq.1).OR.(nabs.eq.4)) then                 

!     SR, 29-08-2005: Include the slice numer!
              if(dowrite_impact) then
      write(48,'(i4,(1x,f6.3),(1x,f8.6),4(1x,e19.10),i2,2(1x,i7))')     &
     &icoll,c_rotation,                                                 &
     &sInt + sp + (dble(j_slices)-1d0) * c_length,                         &!hr09
     &x_flk*1d3, xp_flk*1d3, y_flk*1d3, yp_flk*1d3,                     &
     &nabs,name(j),iturn
              endif
!
!     Finally, the actual coordinate change to 99 mm
              if (nabs.eq.1) then
                 fracab = fracab + 1
                 x = 99.99d-3
                 z = 99.99d-3
                 part_abs_pos_local(j) = ie
                 part_abs_turn_local(j) = iturn
                 lint(j) = zlm
              endif
            endif
          endif
!!! /* end RB fix */

!
!++  Do the rest drift, if particle left collimator early
!  DRIFT PART
          if (nabs.ne.1 .and. zlm.gt.0.) then
             drift_length = (length-(s+sp))
             if (drift_length.gt.1d-15) then
               if(iexact.eq.0) then
                 x  = x + xp * drift_length
                 z  = z + zp * drift_length
                 sp = sp + drift_length
               else
                 zpj = sqrt(1d0-xp**2-zp**2)
                 x = x + drift_length*(xp/zpj)
                 z = z + drift_length*(zp/zpj)
                 sp = sp + drift_length
               endif
             endif
             lint(j) = zlm - drift_length
          endif
!
!++  Transform back to particle coordinates with opening and offset
!
          if (x.lt.99.0d-3) then
!
!++  Include collimator tilt
!
             if (tiltangle.gt.0.d0) then                                 !hr09
                x  = x  + tiltangle*c_length
                xp = xp + tiltangle
             elseif (tiltangle.lt.0.d0) then                             !hr09
                x  = x + tiltangle*c_length
                xp = xp + tiltangle
!
+if crlibm
                x  = x - sin_rn(tiltangle) * c_length
+ei
+if .not.crlibm
                x  = x - sin(tiltangle) * c_length
+ei
            endif
!
!++  Transform back to particle coordinates with opening and offset
!
            z00 = z
            x00 = x + mirror*c_offset
            x = (x + c_aperture/2d0) + mirror*c_offset                   !hr09
!
!++  Now mirror at the horizontal axis for negative X offset
!
            x    = mirror * x
            xp   = mirror * xp
!
!++  Last do rotation into collimator frame
!
+if crlibm
            x_in(j)  = x  *cos_rn(-1d0*c_rotation) +                    &
+ei
+if .not.crlibm
            x_in(j)  = x  *cos(-1d0*c_rotation) +                       &
+ei
+if crlibm
     &z  *sin_rn(-1d0*c_rotation)
+ei
+if .not.crlibm
     &z  *sin(-1d0*c_rotation)
+ei
+if crlibm
            y_in(j)  = z  *cos_rn(-1d0*c_rotation) -                    &
+ei
+if .not.crlibm
            y_in(j)  = z  *cos(-1d0*c_rotation) -                       &
+ei
+if crlibm
     &x  *sin_rn(-1d0*c_rotation)
+ei
+if .not.crlibm
     &x  *sin(-1d0*c_rotation)
+ei
+if crlibm
            xp_in(j) = xp *cos_rn(-1d0*c_rotation) +                    &
+ei
+if .not.crlibm
            xp_in(j) = xp *cos(-1d0*c_rotation) +                       &
+ei
+if crlibm
     &zp *sin_rn(-1d0*c_rotation)
+ei
+if .not.crlibm
     &zp *sin(-1d0*c_rotation)
+ei
+if crlibm
            yp_in(j) = zp *cos_rn(-1d0*c_rotation) -                    &
+ei
+if .not.crlibm
            yp_in(j) = zp *cos(-1d0*c_rotation) -                       &
+ei
+if crlibm
     &xp *sin_rn(-1d0*c_rotation)
+ei
+if .not.crlibm
     &xp *sin(-1d0*c_rotation)
+ei
!
            if (( (icoll.eq.ipencil                                      &
     &.and. iturn.eq.1)   .or.                                          &
     &(iturn.eq.1 .and. ipencil.eq.999 .and.                            &
     &icoll.le.nprim .and.                                              &
     &(j.ge.(icoll-1)*nev/nprim) .and.                                  &
     &(j.le.(icoll)*nev/nprim)                                          &
     &)  ).and.(pencil_distr.ne.3)) then    ! RB: adding condition that this shouldn't be done if pencil_distr=3
!
               x00  = mirror * x00
+if crlibm
               x_in(j)  = x00  *cos_rn(-1d0*c_rotation) +
+ei
+if .not.crlibm
               x_in(j)  = x00  *cos(-1d0*c_rotation) +
+ei
+if crlibm
     &z00  *sin_rn(-1d0*c_rotation)
+ei
+if .not.crlibm
     &z00  *sin(-1d0*c_rotation)
+ei
+if crlibm
               y_in(j)  = z00  *cos_rn(-1d0*c_rotation) -               &
+ei
+if .not.crlibm
               y_in(j)  = z00  *cos(-1d0*c_rotation) -                  &
+ei
+if crlibm
     &x00  *sin_rn(-1d0*c_rotation)
+ei
+if .not.crlibm
     &x00  *sin(-1d0*c_rotation)
+ei
!
               xp_in(j) = xp_in(j) + mirror*xp_pencil0
               yp_in(j) = yp_in(j) + mirror*yp_pencil0
               x_in(j) = x_in(j) + mirror*x_pencil(icoll)
               y_in(j) = y_in(j) + mirror*y_pencil(icoll)
            endif
!
            p_in(j) = (1d0 + dpop) * p0
!     SR, 30-08-2005: add the initial position of the slice
            s_in(j) = sp + (dble(j_slices)-1d0) * c_length               !hr09
!            s_in(j) = s_in(j) + sp
!
          else
            x_in(j)  = x
            y_in(j)  = z
          endif
!
! output for comparing the particle in accelerator frame 
!
c$$$          if(dowrite_impact) then
c$$$             write(9996,'(i5,1x,i7,1x,i2,1x,i1,2(1x,f5.3),8(1x,e17.9))')  &
c$$$     &            name(j),iturn,icoll,nabs,                             &
c$$$     &            s_in(j),                                              &
c$$$     &            s+sp + (dble(j_slices)-1d0) * c_length,               &!hr09
c$$$     &            x_in(j),xp_in(j),y_in(j),yp_in(j),                    &
c$$$     &            x,xp,z,zp
c$$$          endif
!
!++  End of check for particles not being lost before
!
!        endif
!
!        IF (X.GT.99.00) WRITE(*,*) 'After : ', X, X_IN(J)
!
!++  End of loop over all particles
!
 777  continue
      end do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!      WRITE(*,*) 'Number of particles:            ', Nev
!      WRITE(*,*) 'Number of particle hits:        ', Nhit
!      WRITE(*,*) 'Number of absorped particles:   ', fracab
!      WRITE(*,*) 'Number of escaped particles:    ', Nhit-fracab
!      WRITE(*,*) 'Fraction of absorped particles: ', 100.*fracab/Nhit
!
      end

!>
!! collimaterhic()
!! ???
!<
      subroutine collimaterhic(c_material, c_length, c_rotation,        &
     &     c_aperture, n_aperture,                                      &
     &     c_offset, c_tilt,                                            &
     &     x_in, xp_in, y_in,                                           &
     &     yp_in, p_in, s_in, np, enom,                                 &
     &     lhit_pos,lhit_turn,                                          &
     &     part_abs_pos_local, part_abs_turn_local,                     &
     &     impact, indiv, lint, onesided,                               &
     &     name)
!
!++  Based on routines by JBJ. Changed by RA 2001.
!
!++  - Deleted all HBOOK stuff.
!++  - Deleted optics routine and all parser routines.
!++  - Replaced RANMAR call by RANLUX call
!++  - Included RANLUX code from CERNLIB into source
!++  - Changed dimensions from CGen(100,nmat) to CGen(200,nmat)
!++  - Replaced FUNPRE with FUNLXP
!++  - Replaced FUNRAN with FUNLUX
!++  - Included all CERNLIB code into source: RANLUX, FUNLXP, FUNLUX,
!++                                         FUNPCT, FUNLZ, RADAPT,
!++                                           RGS56P
!++    with additional entries:             RLUXIN, RLUXUT, RLUXAT,
!++                                           RLUXGO
!++
!++  - Changed program so that Nev is total number of particles
!++    (scattered and not-scattered)
!++  - Added debug comments
!++  - Put real dp/dx
!
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
!
      double precision sx, sz
!
+ca parpro
+ca collpara
+ca dbpencil
+ca info
+ca dbcollim
!
+ca database
!
      double precision x_flk,xp_flk,y_flk,yp_flk
!JUNE2005
      double precision n_aperture  !aperture in m for the vertical plane
!JUNE2005
!DEBUG
      integer event
!DEBUG
      save
!=======================================================================
! Be=1 Al=2 Cu=3 W=4 Pb=5
!
! LHC uses:    Al, 0.2 m
!              Cu, 1.0 m
!
      if (c_material.eq.'BE') then
         mat = 1
      elseif (c_material.eq.'AL') then
         mat = 2
      elseif (c_material.eq.'CU') then
         mat = 3
      elseif (c_material.eq.'W') then
         mat = 4
      elseif (c_material.eq.'PB') then
         mat = 5
      elseif (c_material.eq.'C') then
         mat = 6
      elseif (c_material.eq.'C2') then
         mat = 7
      elseif (c_material.eq.'MoGR') then
         mat = 8
      elseif (c_material.eq.'CuCD') then
         mat = 9
      elseif (c_material.eq.'Mo') then
         mat = 10
      elseif (c_material.eq.'Glid') then
         mat = 11
      elseif (c_material.eq.'Iner') then
         mat = 12
      else
         write(lout,*) 'ERR>  Material not found. STOP', c_material
         call prror(-1)
      endif
!
        length  = c_length
        nev = np
        p0  = enom
!
!++  Initialize scattering processes
!
      call scatin(p0)

! EVENT LOOP,  initial distribution is here a flat distribution with
! xmin=x-, xmax=x+, etc. from the input file
!
      nhit    = 0
      fracab  = 0.d0                                                     !hr09
      mirror  = 1.d0                                                     !hr09
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do j = 1, nev
!
        impact(j) = -1.d0                                                !hr09
        lint(j)   = -1.d0                                                !hr09
        indiv(j)  = -1.d0                                                !hr09
!
        x   = x_in(j)
        xp  = xp_in(j)
        z   = y_in(j)
        zp  = yp_in(j)
        p   = p_in(j)
!        sp  = s_in(J)
        sp   = 0.
        dpop = (p - p0)/p0
!
!++  Transform particle coordinates to get into collimator coordinate
!++  system
!
!++  First check whether particle was lost before
!
!        if (x.lt.99.0*1e-3 .and. z.lt.99.0*1e-3) then
        if (x.lt.99.0*1d-3 .and. z.lt.99.0*1d-3) then
!
!++  First do rotation into collimator frame
!
!JUNE2005
!JUNE2005 CHANGE TO MAKE THE RHIC TREATMENT EASIER...
!JUNE2005
!+if crlibm
!          x  = x_in(j)*cos_rn(c_rotation) +sin_rn(c_rotation)*y_in(j)
!+ei
!+if .not.crlibm
!          x  = x_in(j)*cos(c_rotation) +sin(c_rotation)*y_in(j)
!+ei
!+if crlibm
!          z  = y_in(j)*cos_rn(c_rotation) -sin_rn(c_rotation)*x_in(j)
!+ei
!+if .not.crlibm
!          z  = y_in(j)*cos(c_rotation) -sin(c_rotation)*x_in(j)
!+ei
!+if crlibm
!          xp = xp_in(j)*cos_rn(c_rotation)+sin_rn(c_rotation)*yp_in(j)
!+ei
!+if .not.crlibm
!          xp = xp_in(j)*cos(c_rotation)+sin(c_rotation)*yp_in(j)
!+ei
!+if crlibm
!          zp = yp_in(j)*cos_rn(c_rotation)-sin_rn(c_rotation)*xp_in(j)
!+ei
!+if .not.crlibm
!          zp = yp_in(j)*cos(c_rotation)-sin(c_rotation)*xp_in(j)
!+ei
          x  = -1d0*x_in(j)
          z  = -1d0*y_in(j)
          xp = -1d0*xp_in(j)
          zp = -1d0*yp_in(j)
!JUNE2005
!
!++  For one-sided collimators consider only positive X. For negative
!++  X jump to the next particle
!
!GRD          IF (ONESIDED .AND. X.LT.0) GOTO 777
!JUNE2005          if (onesided .and. x.lt.0d0 .or. z.gt.0d0) goto 777
          if (onesided .and. (x.lt.0d0 .and. z.gt.0d0)) goto 777
!
!++  Now mirror at the horizontal axis for negative X offset
!
!GRD
!GRD THIS WE HAVE TO COMMENT OUT IN CASE OF RHIC BECAUSE THERE ARE
!GRD ONLY ONE-SIDED COLLIMATORS
!GRD
!          IF (X.LT.0) THEN
!            MIRROR = -1.
!            tiltangle = -1.*C_TILT(2)
!          ELSE
!            MIRROR = 1.
            tiltangle = c_tilt(1)
!          ENDIF
!          X  = MIRROR * X
!          XP = MIRROR * XP
!GRD
!
!++  Shift with opening and offset
!
          x  = (x - c_aperture/2d0) - mirror*c_offset                    !hr09
!GRD
!GRD SPECIAL FEATURE TO TAKE INTO ACCOUNT THE PARTICULAR SHAPE OF RHIC PRIMARY COLLIMATORS
!GRD
!JUNE2005  HERE WE ADD THE ABILITY TO HAVE 2 DIFFERENT OPENINGS FOR THE TWO PLANES
!JUNE2005  OF THE PRIMARY COLLIMATOR OF RHIC
!JUNE2005
!          z  = z + c_aperture/2 + mirror*c_offset
          z  = (z + n_aperture/2d0) + mirror*c_offset                    !hr09
!JUNE2005
!          if(iturn.eq.1)                                                &
!     &write(*,*) 'check ',x,xp,z,zp,c_aperture,n_aperture
!JUNE2005
!
!++  Include collimator tilt
!
          if (tiltangle.gt.0.) then
            xp = xp - tiltangle
          elseif (tiltangle.lt.0.) then
+if crlibm
            x  = x + sin_rn(tiltangle) * c_length
+ei
+if .not.crlibm
            x  = x + sin(tiltangle) * c_length
+ei
            xp = xp - tiltangle
          endif
!
!++  For selected collimator, first turn reset particle distribution
!++  to simple pencil beam
!
            nprim = 3
            if ( (icoll.eq.ipencil                                      &
     &.and. iturn.eq.1) .or.                                            &
     &(iturn.eq.1 .and. ipencil.eq.999 .and.                            &
     &icoll.le.nprim .and.                                              &
     &(j.ge.(icoll-1)*nev/nprim) .and.                                  &
     &(j.le.(icoll)*nev/nprim)                                          &
     &)  ) then
              x    = pencil_dx(icoll)
              xp   = 0.d0                                                !hr09
              z    = 0.d0                                                !hr09
              zp   = 0.d0                                                !hr09
              dpop = 0.d0                                                !hr09
              if(rndm4().lt.0.5) mirror = -1d0*abs(mirror)               !hr09
              if(rndm4().ge.0.5) mirror = abs(mirror)
            endif
!
!++  particle passing above the jaw are discarded => take new event
!++  entering by the face, shorten the length (zlm) and keep track of
!++  entrance longitudinal coordinate (keeps) for histograms
!
!++  The definition is that the collimator jaw is at x>=0.
!
!++  1) Check whether particle hits the collimator
!
          hit     =  .false.
          s       =  0.d0                                                !hr09
          keeps   =  0.d0                                                !hr09
          zlm     =  -1.0d0 * length
!
!GRD
!JUNE2005          if (x.ge.0d0 .and. z.le.0d0) then
          if (x.ge.0d0 .and. z.le.0d0) then
             goto 10
!
!++  Particle hits collimator and we assume interaction length ZLM equal
!++  to collimator length (what if it would leave collimator after
!++  small length due to angle???)
!
!JUNE2005
!            zlm = length
!            impact(j) = max(x,(-1d0*z))
!            if(impact(j).eq.x) then
!               indiv(j) = xp
!            else
!               indiv(j) = zp
!            endif
!          endif
!JUNE2005
!GRD
!JUNE2005          if(x.lt.0d0.and.z.gt.0d0.and.xp.le.0d0.and.zp.ge.0d0) then
          elseif(x.lt.0d0.and.z.gt.0d0.and.xp.le.0d0                    &
     &.and.zp.ge.0d0) then
             goto 20
!GRD
!JUNE2005          if(x.lt.0d0.and.z.gt.0d0.and.xp.le.0d0.and.zp.ge.0d0) then
!
!++  Particle does not hit collimator. Interaction length ZLM is zero.
!
!JUNE2005            zlm = 0.
!JUNE2005          endif
!GRD
!JUNE2005          if (x.lt.0d0.and.z.gt.0d0.and.xp.gt.0d0.and.zp.ge.0d0) then
!JUNE2005
!            zlm = 0.
!          endif
!JUNE2005
!
!JUNE2005
!JUNE2005 THAT WAS PIECE OF CAKE; NOW COMES THE TRICKY PART...
!JUNE2005
!JUNE2005 THE IDEA WOULD BE TO FIRST LIST ALL THE IMPACT
!JUNE2005 POSSIBILITIES, THEN SEND VIA GOTO TO THE CORRECT
!JUNE2005 TREATMENT
!JUNE2005
          elseif((x.lt.0d0).and.(z.le.0d0)) then
             goto 100
          elseif((x.ge.0d0).and.(z.gt.0d0)) then
             goto 200
          elseif((x.lt.0d0).and.(xp.gt.0d0)) then
             goto 300
          elseif((z.gt.0d0).and.(zp.lt.0d0)) then
             goto 400
          endif
!GRD
 10         continue
            event = 10
            zlm = length
            impact(j) = max(x,(-1d0*z))
            if(impact(j).eq.x) then
               indiv(j) = xp
            else
               indiv(j) = zp
            endif
            goto 999
!GRD
 20         continue
            event = 20
            zlm = 0.d0                                                   !hr09
            goto 999
!GRD
 100        continue
            event = 100
            zlm = length
            impact(j) = -1d0*z
            indiv(j) = zp
            goto 999
!GRD
 200        continue
            event = 200
            zlm = length
            impact(j) = x
            indiv(j) = xp
            goto 999
!GRD
!JUNE2005
!JUNE2005 HERE ONE HAS FIRST TO CHECK IF THERE'S NOT A HIT IN THE
!JUNE2005 OTHER PLANE AT THE SAME TIME
!JUNE2005
 300        continue
            event = 300
            if(z.gt.0d0.and.zp.lt.0d0) goto 500
!
!++  Calculate s-coordinate of interaction point
!
            s = (-1.0d0*x) / xp
            if (s.le.0d0) then
              write(lout,*) 'S.LE.0 -> This should not happen (1)'
              call prror(-1)
            endif
!
            if (s .lt. length) then
              zlm = length - s
              impact(j) = 0.d0                                           !hr09
              indiv(j) = xp
            else
              zlm = 0.d0                                                 !hr09
            endif
            goto 999
!GRD
 400        continue
            event = 400
!JUNE2005          if (x.lt.0d0.and.z.gt.0d0.and.xp.le.0d0.and.zp.lt.0d0) then
!
!++  Calculate s-coordinate of interaction point
!
            s = (-1.0d0*z) / zp
            if (s.le.0) then
              write(lout,*) 'S.LE.0 -> This should not happen (2)'
              call prror(-1)
            endif
!
            if (s .lt. length) then
              zlm = length - s
              impact(j) = 0.d0                                           !hr09
              indiv(j) = zp
            else
              zlm = 0.d0                                                 !hr09
            endif
!JUNE2005          endif
!GRD
            goto 999
!GRD
!GRD
!JUNE2005          if (x.lt.0d0.and.z.gt.0d0.and.xp.gt.0d0.and.zp.lt.0d0) then
 500        continue
            event = 500
!
!++  Calculate s-coordinate of interaction point
!
            sx = (-1.0d0*x) / xp
            sz = (-1.0d0*z) / zp
!
            if(sx.lt.sz) s=sx
            if(sx.ge.sz) s=sz
!
            if (s.le.0d0) then
              write(lout,*) 'S.LE.0 -> This should not happen (3)'
              call prror(-1)
            endif
!
            if (s .lt. length) then
              zlm = length - s
              impact(j) = 0.d0                                           !hr09
              if(s.eq.sx) then
                indiv(j) = xp
              else
                indiv(j) = zp
              endif
            else
              zlm = 0.d0                                                 !hr09
            endif
!
!JUNE2005          endif
!GRD
!GRD
 999      continue
!JUNE2005
!          write(*,*) 'event ',event,x,xp,z,zp
!          if(impact(j).lt.0d0) then
!             if(impact(j).ne.-1d0)                                      &
!     &write(*,*) 'argh! ',impact(j),x,xp,z,zp,s,event
!          endif
!          if(impact(j).ge.0d0) then
!      write(*,*) 'impact! ',impact(j),x,xp,z,zp,s,event
!          endif
!JUNE2005
!
!++  First do the drift part
!
          drift_length = length - zlm
          if (drift_length.gt.0.d0) then                                 !hr09
            x  = x + xp* drift_length
            z  = z + zp * drift_length
            sp = sp + drift_length
          endif
!
!++  Now do the scattering part
!
          if (zlm.gt.0.d0) then                                          !hr09
            nhit = nhit + 1
!            WRITE(*,*) J,X,XP,Z,ZP,SP,DPOP
!DEBUG
!            write(*,*) 'abs?',s,zlm
!DEBUG
!JUNE2005
!JUNE2005 IN ORDER TO HAVE A PROPER TREATMENT IN THE CASE OF THE VERTICAL
!JUNE2005 PLANE, CHANGE AGAIN THE FRAME FOR THE SCATTERING SUBROUTINES...
!JUNE2005
            if(event.eq.100.or.event.eq.400) then
!GRD first go back into normal frame...
               x = (x + c_aperture/2d0) + mirror*c_offset                !hr09
               z = (z - n_aperture/2d0) - mirror*c_offset                !hr09
               x = -1d0*x
               xp = -1d0*xp
               z = -1d0*z
               zp = -1d0*zp
!GRD ...then do as for a vertical collimator
               x = z
               xp = zp
               z = -1d0*x
               zp = -1d0*x
               x  = (x - n_aperture/2d0) - mirror*c_offset               !hr09
               z  = (z + c_aperture/2d0) + mirror*c_offset               !hr09
            endif
!JUNE2005
!     RB: add new input arguments to jaw icoll,iturn,ipart for writeout
            call jaw(s, nabs, icoll, iturn, name(j), dowrite_impact)

!DEBUG
!            write(*,*) 'abs?',nabs
!DEBUG
!JUNE2005
!JUNE2005 ...WITHOUT FORGETTING TO GO BACK TO THE "ORIGINAL" FRAME AFTER THE
!JUNE2005 ROUTINES, SO AS TO AVOID RIDICULOUS VALUES FOR KICKS IN EITHER PLANE
            if(event.eq.100.or.event.eq.400) then
!GRD first go back into normal frame...
               x = (x + n_aperture/2d0) + mirror*c_offset                !hr09
               z = (z - c_aperture/2d0) - mirror*c_offset                !hr09
               x = -1d0*z
               xp = -1d0*zp
               z = x
               zp = xp
!GRD ...then go back to face the horizontal jaw at 180 degrees
               x = -1d0*x
               xp = -1d0*xp
               z = -1d0*z
               zp = -1d0*zp
               x  = (x - c_aperture/2d0) - mirror*c_offset               !hr09
               z  = (z + n_aperture/2d0) + mirror*c_offset               !hr09
            endif
!JUNE2005
            lhit_pos(j)  = ie
            lhit_turn(j) = iturn
!
!++  If particle is absorbed then set x and y to 99.99 mm
!
            if (nabs.eq.1) then
!APRIL2005
!TO WRITE FLUKA INPUT CORRECTLY, WE HAVE TO GO BACK IN THE MACHINE FRAME
            if (tiltangle.gt.0.d0) then                                  !hr09
              x  = x  + tiltangle*c_length
              xp = xp + tiltangle
            elseif (tiltangle.lt.0.d0) then                              !hr09
              x  = x + tiltangle*c_length
              xp = xp + tiltangle
!
              x  = x - sin(tiltangle) * c_length
            endif
!
!++  Transform back to particle coordinates with opening and offset
!
            x = (x + c_aperture/2d0) + mirror*c_offset                   !hr09
!GRD
!JUNE2005  OF COURSE WE ADAPT ALSO THE PREVIOUS CHANGE WHEN SHIFTING BACK
!JUNE2005  TO  THE ACCELERATOR FRAME...
!            z = z - c_aperture/2 - mirror*c_offset
            z = (z - n_aperture/2d0) - mirror*c_offset                   !hr09
!JUNE2005
!
!++   Last do rotation into collimator frame
!
                  x_flk  = -1d0*x
                  y_flk  = -1d0*z
                  xp_flk = -1d0*xp
                  yp_flk = -1d0*zp
!NOW WE CAN WRITE THE COORDINATES OF THE LOST PARTICLES
              if(dowrite_impact) then
      write(48,'(i4,(2x,f5.3),(2x,f8.6),4(1x,e16.7),2x,i2,2x,i5)')      &
     &icoll,c_rotation,s+sp,                                            &
     &x_flk*1d3, xp_flk*1d3, y_flk*1d3, yp_flk*1d3,                     &
     &nabs,name(j)
              endif
!APRIL2005
              fracab = fracab + 1
!              x = 99.99*1e-3
!              z = 99.99*1e-3
              x = 99.99*1.0d-3
              z = 99.99*1.0d-3
              part_abs_pos_local(j) = ie
              part_abs_turn_local(j) = iturn
              lint(j) = zlm
            endif
          endif
!
!++  Do the rest drift, if particle left collimator early
!
          if (nabs.ne.1 .and. zlm.gt.0.d0) then                          !hr09
            drift_length = (length-(s+sp))
!            if (drift_length.gt.1.e-15) then
            if (drift_length.gt.1.0d-15) then
!              WRITE(*,*) J, DRIFT_LENGTH
              x  = x + xp * drift_length
              z  = z + zp * drift_length
              sp = sp + drift_length
            endif
            lint(j) = zlm - drift_length
          endif
!
!++  Transform back to particle coordinates with opening and offset
!
!          if (x.lt.99.0*1e-3 .and. z.lt.99.0*1e-3) then
          if (x.lt.99.0*1d-3 .and. z.lt.99.0*1d-3) then
!
!++  Include collimator tilt
!
            if (tiltangle.gt.0.d0) then                                  !hr09
              x  = x  + tiltangle*c_length
              xp = xp + tiltangle
            elseif (tiltangle.lt.0.d0) then                              !hr09
              x  = x + tiltangle*c_length
              xp = xp + tiltangle
!
+if crlibm
              x  = x - sin_rn(tiltangle) * c_length
+ei
+if .not.crlibm
              x  = x - sin(tiltangle) * c_length
+ei
            endif
!
!++  Transform back to particle coordinates with opening and offset
!
            z00 = z
            x00 = x + mirror*c_offset
            x = (x + c_aperture/2d0) + mirror*c_offset                   !hr09
!GRD
!JUNE2005  OF COURSE WE ADAPT ALSO THE PREVIOUS CHANGE WHEN SHIFTING BACK
!JUNE2005  TO  THE ACCELERATOR FRAME...
!            z = z - c_aperture/2 - mirror*c_offset
            z = (z - n_aperture/2d0) - mirror*c_offset                   !hr09
!JUNE2005
!
!++  Now mirror at the horizontal axis for negative X offset
!
            x    = mirror * x
            xp   = mirror * xp
!
!++  Last do rotation into collimator frame
!
!JUNE2005
!+if crlibm
!            x_in(j)  = x  *cos_rn(-1.*c_rotation) +                     &
!+ei
!+if .not.crlibm
!            x_in(j)  = x  *cos(-1.*c_rotation) +                        &
!+ei
!+if crlibm
!     &z  *sin_rn(-1.*c_rotation)
!+ei
!+if .not.crlibm
!     &z  *sin(-1.*c_rotation)
!+ei
!+if crlibm
!            y_in(j)  = z  *cos_rn(-1.*c_rotation) -                     &
!+ei
!+if .not.crlibm
!            y_in(j)  = z  *cos(-1.*c_rotation) -                        &
!+ei
!+if crlibm
!     &x  *sin_rn(-1.*c_rotation)
!+ei
!+if .not.crlibm
!     &x  *sin(-1.*c_rotation)
!+ei
!+if crlibm
!            xp_in(j) = xp *cos_rn(-1.*c_rotation) +                     &
!+ei
!+if .not.crlibm
!            xp_in(j) = xp *cos(-1.*c_rotation) +                        &
!+ei
!+if crlibm
!     &zp *sin_rn(-1.*c_rotation)
!+ei
!+if .not.crlibm
!     &zp *sin(-1.*c_rotation)
!+ei
!+if crlibm
!            yp_in(j) = zp *cos_rn(-1.*c_rotation) -                     &
!+ei
!+if .not.crlibm
!            yp_in(j) = zp *cos(-1.*c_rotation) -                        &
!+ei
!+if crlibm
!     &xp *sin_rn(-1.*c_rotation)
!+ei
!+if .not.crlibm
!     &xp *sin(-1.*c_rotation)
!+ei
            x_in(j) = -1d0*x
            y_in(j) = -1d0*z
            xp_in(j) = -1d0*xp
            yp_in(j) = -1d0*zp
!JUNE2005
!
            if ( (icoll.eq.ipencil                                      &
     &.and. iturn.eq.1)   .or.                                          &
     &(iturn.eq.1 .and. ipencil.eq.999 .and.                            &
     &icoll.le.nprim .and.                                              &
     &(j.ge.(icoll-1)*nev/nprim) .and.                                  &
     &(j.le.(icoll)*nev/nprim)                                          &
     &)  ) then
!
               x00  = mirror * x00
+if crlibm
               x_in(j)  = x00  *cos_rn(-1.d0*c_rotation) +              &!hr09
+ei
+if .not.crlibm
               x_in(j)  = x00  *cos(-1.d0*c_rotation) +                 &!hr09
+ei
+if crlibm
     &z00  *sin_rn(-1.d0*c_rotation)                                     !hr09
+ei
+if .not.crlibm
     &z00  *sin(-1.d0*c_rotation)                                        !hr09
+ei
+if crlibm
               y_in(j)  = z00  *cos_rn(-1.d0*c_rotation) -              &!hr09
+ei
+if .not.crlibm
               y_in(j)  = z00  *cos(-1.d0*c_rotation) -                 &!hr09
+ei
+if crlibm
     &x00  *sin_rn(-1.d0*c_rotation)                                     !hr09
+ei
+if .not.crlibm
     &x00  *sin(-1.d0*c_rotation)                                        !hr09
+ei
!
               xp_in(j) = xp_in(j) + mirror*xp_pencil0
               yp_in(j) = yp_in(j) + mirror*yp_pencil0
               x_in(j) = x_in(j) + mirror*x_pencil(icoll)
               y_in(j) = y_in(j) + mirror*y_pencil(icoll)
            endif
!
            p_in(j) = (1d0 + dpop) * p0                                  !hr09
            s_in(j) = s_in(j) + sp
!
          else
            x_in(j)  = x
            y_in(j)  = z
          endif
!
!++  End of check for particles not being lost before
!
        endif
!
!        IF (X.GT.99.00) WRITE(*,*) 'After : ', X, X_IN(J)
!
!++  End of loop over all particles
!
 777  continue
      end do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!      WRITE(*,*) 'Number of particles:            ', Nev
!      WRITE(*,*) 'Number of particle hits:        ', Nhit
!      WRITE(*,*) 'Number of absorped particles:   ', fracab
!      WRITE(*,*) 'Number of escaped particles:    ', Nhit-fracab
!      WRITE(*,*) 'Fraction of absorped particles: ', 100.*fracab/Nhit
!
      end
!
!-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----
!! END collimaterhic()

!>
!! ichoix(ma)
!! Select a scattering type (elastic, sd, inelastic, ...)
!<
      function ichoix(ma)
      implicit none
+if crlibm
+ca crlibco
+ei
+ca interac
      integer ma,i,ichoix
      double precision aran
      real rndm4
      aran=dble(rndm4())
      i=1
  10  if ( aran.gt.cprob(i,ma) ) then
          i=i+1
          goto 10
      endif
      ichoix=i
      return
      end

!>
!! gettran(inter,xmat,p)
!! This function determines: GETTRAN - rms transverse momentum transfer
!! Note: For single-diffractive scattering the vector p of momentum
!! is modified (energy loss is applied)
!<
      function gettran(inter,xmat,p)

      implicit none
+if crlibm
+ca crlibco
+ei
+ca interac
      integer inter,length,xmat
      double precision p,gettran,t,xm2,bsd
      real rndm4,truth,xran(1)

! inter=2: Nuclear Elastic, 3: pp Elastic, 4: Single Diffractive, 5:Coulomb
+if .not.merlinscatter
      if        ( inter.eq.2 ) then
+if crlibm
           gettran = (-1d0*log_rn(dble(rndm4())))/bn(xmat)               !hr09
+ei
+if .not.crlibm
           gettran = (-1d0*log(dble(rndm4())))/bn(xmat)                  !hr09
+ei
!
         elseif ( inter .eq. 3 ) then
+if crlibm
           gettran = (-1d0*log_rn(dble(rndm4())))/bpp                    !hr09
+ei
+if .not.crlibm
           gettran = (-1d0*log(dble(rndm4())))/bpp                       !hr09
+ei
!
         elseif ( inter .eq. 4 ) then
+if crlibm
           xm2 = exp_rn( dble(rndm4()) * xln15s )
+ei
+if .not.crlibm
           xm2 = exp( dble(rndm4()) * xln15s )
+ei
           p = p  *(1.d0 - xm2/ecmsq)
           if ( xm2 .lt. 2.d0 ) then
                bsd = 2.d0 * bpp
              elseif (( xm2 .ge. 2.d0 ).and. ( xm2 .le. 5.d0 )) then
                bsd = ((106.d0-17.d0*xm2) *  bpp )/ 36.d0                !hr09
              elseif ( xm2 .gt. 5.d0 ) then
                bsd = (7.d0 * bpp) / 12.d0                               !hr09
           endif
+if crlibm
           gettran = (-1d0*log_rn(dble(rndm4())))/bsd                    !hr09
+ei
+if .not.crlibm
           gettran = (-1d0*log(dble(rndm4())))/bsd                       !hr09
+ei

         elseif ( inter.eq.5 ) then
           length=1
           call funlux( cgen(1,mat) , xran, length)
           truth=xran(1)
           t=dble(truth)                                                 !hr09
           gettran = t
      endif
+ei
+if merlinscatter

      if ( inter.eq.2 ) then
           gettran = (-1d0*log(dble(rndm4())))/bn(xmat)                  !hr09

      elseif ( inter .eq. 3 ) then
           call merlinscatter_get_elastic_t(gettran)

      elseif ( inter .eq. 4 ) then
           call merlinscatter_get_sd_xi(xm2)
           call merlinscatter_get_sd_t(gettran)
           p = p  * (1.d0 - (xm2/ecmsq))

      elseif ( inter.eq.5 ) then
           length=1
           call funlux( cgen(1,mat) , xran, length)
           truth=xran(1)
           t=dble(truth)                                                 !hr09
           gettran = t
      endif

+ei
      return
      end

!>
!! tetat(t,p,tx,tz)
!! ???
!!
!<
      subroutine tetat(t,p,tx,tz)
      implicit none
+if crlibm
+ca crlibco
+ei
      double precision t,p,tx,tz,va,vb,va2,vb2,r2,teta
      real rndm4
      teta = sqrt(t)/p
! Generate sine and cosine of an angle uniform in [0,2pi](see RPP)
   10 va  =(2d0*dble(rndm4()))-1d0                                       !hr09
      vb = dble(rndm4())
      va2 = va**2
      vb2 = vb**2
      r2 = va2 + vb2
      if ( r2.gt.1.d0) go to 10
      tx = teta * ((2.d0*va)*vb) / r2                                    !hr09
      tz = teta * (va2 - vb2) / r2
      return
      end

!>
!! ruth(t)
!! Calculate the rutherford scattering cross section
!<
      function ruth(t)
      implicit none
+if crlibm
+ca crlibco
+ei
+ca interac
      real ruth,t
      double precision cnorm,cnform
      parameter(cnorm=2.607d-5,cnform=0.8561d3) ! DM: changed 2.607d-4 to 2.607d-5 to fix Rutherford bug

+if crlibm
      ruth=real((cnorm*exp_rn(((-1d0*dble(t))*cnform)*emr(mcurr)**2))*  &!hr09
     &(zatom(mcurr)/dble(t))**2)
+ei
+if .not.crlibm
      ruth=((cnorm*exp(((-1d0*dble(t))*cnform)*emr(mcurr)**2))*         &!hr09
     &(zatom(mcurr)/dble(t))**2)                                         !hr09
+ei
      end

!>
!! block data scdata
!! Cross section inputs and material property database
!! GRD CHANGED ON 2/2003 TO INCLUDE CODE FOR C, C2 from JBJ (rwa)
!<
      block data scdata
      implicit none
+if crlibm
+ca crlibco
+ei
+ca interac
      integer i
! Total number of materials are defined in nmat
! Number of real materials are defined in nrmat
! The last materials in nmat are 'vacuum' and 'black',see in sub. SCATIN

! Reference data at pRef=450Gev
      data (mname(i),i=1,nrmat)/ 'Be','Al','Cu','W','Pb','C','C2',      &
     & 'MoGR','CuCD', 'Mo', 'Glid', 'Iner'/

      data mname(nmat-1), mname(nmat)/'vacu','blac'/

!GRD IMPLEMENT CHANGES FROM JBJ, 2/2003 RWA
      data (anuc(i),i=1,5)/ 9.01d0,26.98d0,63.55d0,183.85d0,207.19d0/
      data (anuc(i),i=6,7)/12.01d0,12.01d0/
      data (anuc(i),i=8,nrmat)/13.53d0,25.24d0,95.96d0,63.15d0,166.7d0/

      data (zatom(i),i=1,5)/ 4d0, 13d0, 29d0, 74d0, 82d0/
      data (zatom(i),i=6,7)/ 6d0, 6d0/
      data (zatom(i),i=8,nrmat)/ 6.65d0, 11.9d0, 42d0, 28.8d0, 67.7d0/

      data (rho(i),i=1,5)/ 1.848d0, 2.70d0, 8.96d0, 19.3d0, 11.35d0/
      data (rho(i),i=6,7)/ 1.67d0, 4.52d0/
      data (rho(i),i=8,nrmat)/ 2.5d0, 5.4d0, 10.22d0, 8.93d0, 18d0/

      data (radl(i),i=1,5)/ 0.353d0,0.089d0,0.0143d0,0.0035d0,0.0056d0/
      data (radl(i),i=6,7)/ 0.2557d0, 0.094d0/
      data (radl(i),i=8,nrmat)/ 0.1193d0, 0.0316d0, 0.0096d0, 0.0144d0, &
     & 0.00385d0/
      data radl(nmat-1),radl(nmat)/ 1.d12, 1.d12 /

!MAY06-GRD value for Tungsten (W) not stated
      data (emr(i),i=1,5)/  0.22d0, 0.302d0, 0.366d0, 0.520d0, 0.542d0/
      data (emr(i),i=6,7)/  0.25d0, 0.25d0/
      data (emr(i),i=8,nrmat)/ 0.25d0, 0.308d0, 0.481d0, 0.418d0,       &
     & 0.578d0/

      data tlcut / 0.0009982d0/
      data (hcut(i),i=1,5)/0.02d0, 0.02d0, 3*0.01d0/
      data (hcut(i),i=6,7)/0.02d0, 0.02d0/
      data (hcut(i),i=8,nrmat)/0.02d0, 0.02d0, 0.02d0, 0.02d0, 0.02d0/
!

      data (dpodx(i),i=1,5)/ .55d0, .81d0, 2.69d0, 5.79d0, 3.4d0 /
      data (dpodx(i),i=6,7)/ .75d0, 1.5d0 /

!October 2013
!Mean excitation energy (GeV) values added by Claudia for Bethe-Bloch implementation:
      data (exenergy(i),i=1,5)/ 63.7e-9,166e-9, 322e-9, 727e-9, 823e-9 /
      data (exenergy(i),i=6,7)/ 78e-9, 78.0e-9 /
      data (exenergy(i),i=8,nrmat)/ 87.1e-9, 152.9e-9, 424e-9, 320.8e-9,&
     & 682.2e-9/
 
! All cross-sections are in barns,nuclear values from RPP at 20geV
! Coulomb is integerated above t=tLcut[Gev2] (+-1% out Gauss mcs)

! in Cs and CsRef,1st index: Cross-sections for processes
! 0:Total, 1:absorption, 2:nuclear elastic, 3:pp or pn elastic
! 4:Single Diffractive pp or pn, 5:Coulomb for t above mcs

! Claudia 2013: updated cross section values. Unit: Barn. New 2013:
      data csref(0,1),csref(1,1),csref(5,1)/0.271d0, 0.192d0, 0.0035d-2/
      data csref(0,2),csref(1,2),csref(5,2)/0.643d0, 0.418d0, 0.034d-2/
      data csref(0,3),csref(1,3),csref(5,3)/1.253d0, 0.769d0, 0.153d-2/
      data csref(0,4),csref(1,4),csref(5,4)/2.765d0, 1.591d0, 0.768d-2/
      data csref(0,5),csref(1,5),csref(5,5)/3.016d0, 1.724d0, 0.907d-2/
      data csref(0,6),csref(1,6),csref(5,6)/0.337d0, 0.232d0, 0.0076d-2/
      data csref(0,7),csref(1,7),csref(5,7)/0.337d0, 0.232d0, 0.0076d-2/
      data csref(0,8),csref(1,8),csref(5,8)/0.362d0, 0.247d0, 0.0094d-2/
      data csref(0,9),csref(1,9),csref(5,9)/0.572d0, 0.370d0, 0.0279d-2/
      data csref(0,10),csref(1,10),csref(5,10)/1.713d0,1.023d0,0.265d-2/
      data csref(0,11),csref(1,11),csref(5,11)/1.246d0,0.765d0,0.139d-2/
      data csref(0,12),csref(1,12),csref(5,12)/2.548d0,1.473d0,0.574d-2/

! pp cross-sections and parameters for energy dependence
      data pptref,pperef,sdcoe,pref/0.04d0,0.007d0,0.00068d0,450.0d0/
      data pptco,ppeco,freeco/0.05788d0,0.04792d0,1.618d0/

! Nuclear elastic slope from Schiz et al.,PRD 21(3010)1980
!MAY06-GRD value for Tungsten (W) not stated
      data (bnref(i),i=1,5)/74.7d0,120.3d0,217.8d0,440.3d0,455.3d0/
      data (bnref(i),i=6,7)/70.d0, 70.d0/
      data (bnref(i),i=8,nrmat)/ 76.7d0, 115.0d0, 273.9d0, 208.7d0,      &
     & 392.1d0/
!GRD LAST 2 ONES INTERPOLATED

! Cprob to choose an interaction in iChoix
      data (cprob(0,i),i=1,nmat)/nmat*0.0d0/
      data (cprob(5,i),i=1,nmat)/nmat*1.0d0/
      end

!>
!! scatin(plab)
!! Configure the K2 scattering routine cross sections
!!
!<
      subroutine scatin(plab)
      implicit none
+if merlinscatter
+ca database
+ei
+if crlibm
+ca crlibco
+ei
+ca interac
      integer ma,i
      double precision plab
      real ruth,tlow,thigh
      external ruth
!
      ecmsq = (2d0 * 0.93828d0) * plab                                   !hr09
+if .not.merlinscatter
+if crlibm
      xln15s=log_rn(0.15d0*ecmsq)                                        !hr09
+ei
+if .not.crlibm
      xln15s=log(0.15d0*ecmsq)                                           !hr09
+ei

+if crlibm
!Claudia Fit from COMPETE collaboration points "arXiv:hep-ph/0206172v1 19Jun2002"
      pptot=0.041084d0-0.0023302d0*log_rn(ecmsq)+0.00031514d0*
     &  log_rn(ecmsq)**2

!Claudia used the fit from TOTEM for ppel (in barn)
      ppel=(11.7d0-1.59d0*log_rn(ecmsq)+0.134d0*log_rn(ecmsq)**2)/1000

!Claudia updated SD cross that cointains renormalized pomeron flux (in barn)
      ppsd=(4.3d0+0.3d0*log_rn(ecmsq))/1000

!Claudia new fit for the slope parameter with new data at sqrt(s)=7 TeV from TOTEM
      bpp=7.156d0+1.439d0*log_rn(sqrt(ecmsq))
+ei

+if .not.crlibm
!Claudia Fit from COMPETE collaboration points "arXiv:hep-ph/0206172v1 19Jun2002"
      pptot=0.041084d0-0.0023302d0*log(ecmsq)+0.00031514d0*log(ecmsq)**2

!Claudia used the fit from TOTEM for ppel (in barn)
      ppel=(11.7d0-1.59d0*log(ecmsq)+0.134d0*log(ecmsq)**2)/1000

!Claudia updated SD cross that cointains renormalized pomeron flux (in barn)
      ppsd=(4.3d0+0.3d0*log(ecmsq))/1000
+ei
+ei
+if merlinscatter !No crlibm...
      call merlinscatter_setup(plab,rnd_seed)
      call merlinscatter_setdata(pptot,ppel,ppsd)
+ei
!Claudia new fit for the slope parameter with new data at sqrt(s)=7 TeV from TOTEM
      bpp=7.156d0+1.439d0*log(sqrt(ecmsq))
      
! unmeasured tungsten data,computed with lead data and power laws
      bnref(4) = bnref(5)*(anuc(4) / anuc(5))**(2d0/3d0)
      emr(4) = emr(5) * (anuc(4)/anuc(5))**(1d0/3d0)
   10 format(/' ppRef TOT El     ',4f12.6//)
   11 format(/' pp    TOT El Sd b',4f12.6//)

! Compute cross-sections (CS) and probabilities + Interaction length
! Last two material treated below statement number 100

      tlow=real(tlcut)                                                   !hr09
      do 100 ma=1,nrmat
        mcurr=ma
! prepare for Rutherford differential distribution
        thigh=real(hcut(ma))                                             !hr09
        call funlxp ( ruth , cgen(1,ma) ,tlow, thigh )

! freep: number of nucleons involved in single scattering
        freep(ma) = freeco * anuc(ma)**(1d0/3d0)

! compute pp and pn el+single diff contributions to cross-section
! (both added : quasi-elastic or qel later)
        cs(3,ma) = freep(ma) * ppel
        cs(4,ma) = freep(ma) * ppsd

! correct TOT-CSec for energy dependence of qel
! TOT CS is here without a Coulomb contribution
        cs(0,ma) = csref(0,ma) + freep(ma) * (pptot - pptref)
        bn(ma) = (bnref(ma) * cs(0,ma)) / csref(0,ma)                    !hr09
! also correct inel-CS
        cs(1,ma) = (csref(1,ma) * cs(0,ma)) / csref(0,ma)                !hr09
!
! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
        cs(2,ma) = ((cs(0,ma) - cs(1,ma)) - cs(3,ma)) - cs(4,ma)         !hr09
        cs(5,ma) = csref(5,ma)
! Now add Coulomb
        cs(0,ma) = cs(0,ma) + cs(5,ma)
! Interaction length in meter
      xintl(ma) = (0.01d0*anuc(ma))/(((fnavo * rho(ma))*cs(0,ma))*1d-24) !hr09

   20   format(/1x,a4,' Int.Len. ',f10.6,' CsTot',2f12.4/)

   21   format('  bN freep',2 f12.6,'   emR ',f7.4/)

! Filling CProb with cumulated normalised Cross-sections
        do 50 i=1,4
          cprob(i,ma)=cprob(i-1,ma)+cs(i,ma)/cs(0,ma)

 50     continue

   22   format(i4,' prob CS CsRref',3(f12.5,2x))
  100 continue

! Last two materials for 'vaccum' (nmat-1) and 'full black' (nmat)
      cprob(1,nmat-1)=1d0
      cprob(1,nmat)=1d0
      xintl(nmat-1)=1d12
      xintl(nmat)=0.0d0
  120 format(/1x,a4,' Int.Len. ',e10.3/)
      return
      end

!>
!! jaw(s,nabs,icoll,iturn,ipart,dowrite_impact)
!! ???
!!     RB: adding as input arguments to jaw variables icoll,iturn,ipart
!!         these are only used for the writeout of particle histories
!!
!!++  Input:   ZLM is interaction length
!!++           MAT is choice of material
!!
!!++  Output:  nabs = 1   Particle is absorped
!!++           nabs = 4   Single-diffractive scattering
!!++           dpop       Adjusted for momentum loss (dE/dx)
!!++           s          Exit longitudinal position
!!
!!++  Physics:  If monte carlo interaction length greater than input
!!++            interaction length, then use input interaction length
!!++            Is that justified???
!!
!!     nabs=1....absorption
!!
!<
      subroutine jaw(s,nabs,icoll,iturn,ipart,dowrite_impact)

      implicit none
+if crlibm
+ca crlibco
+ei
!
+ca interac
+ca flukavars
      integer nabs,inter,ichoix,iturn,icoll,ipart,nabs_tmp ! RB: added variables icoll,iturn,ipart for writeout
      logical dowrite_impact
      double precision  m_dpodx, mc_int_l,s_in,I,c_material     !CT, RB, DM
      double precision p,rlen,s,t,gettran,dxp,dzp,p1,zpBef,xpBef,pBef
      real get_dpodx
      real rndm4
!...cne=1/(sqrt(b))
!...dpodx=dE/(dx*c)

!++  Note that the input parameter is dpop. Here the momentum p is
!++  constructed out of this input.

      p=p0*(1.d0+dpop)
      nabs=0
      nabs_tmp=nabs

      if(mat.eq.nmat) then
!++  Collimator treated as black absorber
        nabs=1
        nabs_tmp=nabs
        s=0d0

        if(dowrite_impact) then
! write Coll_Scatter.dat for complete scattering histories
           write(3998,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e14.6))')           &  
     &     icoll,iturn,ipart,nabs_tmp,-1.d0,0.d0,0.d0
        end if
        return
      else if(mat.eq.nmat-1) then
!++  Collimator treated as drift
        s=zlm
        x=x+s*xp
        z=z+s*zp

        if(dowrite_impact) then
! write Coll_Scatter.dat for complete scattering histories
           write(3998,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e14.6))')           &  
     &     icoll,iturn,ipart,nabs_tmp,-1.d0,0.d0,0.d0
        end if

        return
      end if

!++  Initialize the interaction length to input interaction length
      rlen=zlm

!++  Do a step for a point-like interaction. This is a loop with
!++  label 10!!!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!++  Get monte-carlo interaction length.

+if crlibm
10    zlm1=(-1d0*xintl(mat))*log_rn(dble(rndm4()))                       !hr09
+ei
+if .not.crlibm
10    zlm1=(-1d0*xintl(mat))*log(dble(rndm4()))                          !hr09
+ei
      nabs_tmp=0 !! type of interaction reset before following scattering process 

      xpBef=xp ! save angles and momentum before scattering
      zpBef=zp
      pBef=p

      if(zlm1.gt.rlen) then

!++  If the monte-carlo interaction length is longer than the
!++  remaining collimator length, then put it to the remaining
!++  length, do multiple coulomb scattering and return.
!++  LAST STEP IN ITERATION LOOP


       zlm1=rlen
       call mcs(s)
       s=(zlm-rlen)+s                                                    !hr09
!       m_dpodx=get_dpodx(p,mat) ! Claudia 2013
+if .not.merlinscatter
       call calc_ion_loss(mat,p,rlen,m_dpodx)  ! DM routine to include tail
       p=p-m_dpodx*s
+ei
+if merlinscatter
!void calc_ion_loss_merlin_(double* p, double* ElectronDensity, double* PlasmaEnergy, double* MeanIonisationEnergy, double* result)
      call merlinscatter_calc_ion_loss(p,edens(mat),                     &
     & pleng(mat),exenergy(mat),s,m_dpodx)
       p=p-m_dpodx
+ei

       dpop=(p-p0)/p0
       if(dowrite_impact) then
! write Coll_Scatter.dat for complete scattering histories
          write(3998,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e18.10))')           &
     &    icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef
       end if
       return
      end if

!++  Otherwise do multi-coulomb scattering.
!++  REGULAR STEP IN ITERATION LOOP

      call mcs(s)

!++  Check if particle is outside of collimator (X.LT.0) after
!++  MCS. If yes, calculate output longitudinal position (s),
!++  reduce momentum (output as dpop) and return.
!++  PARTICLE LEFT COLLIMATOR BEFORE ITS END.

      if(x.le.0d0) then
       s=(zlm-rlen)+s                                                    !hr09

+if .not.merlinscatter
       call calc_ion_loss(mat,p,rlen,m_dpodx)
       p=p-m_dpodx*s
+ei
+if merlinscatter
       call merlinscatter_calc_ion_loss(p,edens(mat),                    &
     & pleng(mat),exenergy(mat),s,m_dpodx)
       p=p-m_dpodx
+ei
       dpop=(p-p0)/p0

       if(dowrite_impact) then
! write Coll_Scatter.dat for complete scattering histories
          write(3998,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e18.10))')           &
     &    icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef
       end if

       return
      end if

!++  Check whether particle is absorbed. If yes, calculate output
!++  longitudinal position (s), reduce momentum (output as dpop)
!++  and return.
!++  PARTICLE WAS ABSORPED INSIDE COLLIMATOR DURING MCS.

      inter=ichoix(mat)

      nabs=inter
      nabs_tmp=nabs

!     RB, DM: save coordinates before interaction for writeout to FLUKA_impacts.dat
      xInt=x
      xpInt=xp
      yInt=z
      ypInt=zp
      sInt=(zlm-rlen)+zlm1                                                 !hr09

      if(inter.eq.1) then
       s=(zlm-rlen)+zlm1                                                 !hr09

+if .not.merlinscatter
       call calc_ion_loss(mat,p,rlen,m_dpodx)
       p=p-m_dpodx*s
+ei
+if merlinscatter
       call merlinscatter_calc_ion_loss(p,edens(mat),                    &
     & pleng(mat),exenergy(mat),s,m_dpodx)
       p=p-m_dpodx
+ei

       dpop=(p-p0)/p0

! write Coll_Scatter.dat for complete scattering histories
           write(3998,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e14.6))')           &  
     &     icoll,iturn,ipart,nabs_tmp,-1.d0,0.d0,0.d0
       return
      end if

!++  Now treat the other types of interaction, as determined by ICHOIX:

!++      Nuclear-Elastic:          inter = 2
!++      pp Elastic:               inter = 3
!++      Single-Diffractive:       inter = 4    (changes momentum p)
!++      Coulomb:                  inter = 5

!++  As the single-diffractive interaction changes the momentum, save
!++  input momentum in p1.
      p1 = p

!++  Gettran returns some monte carlo number, that, as I believe, gives
!++  the rms transverse momentum transfer.
      t = gettran(inter,mat,p)

!++  Tetat calculates from the rms transverse momentum transfer in
!++  monte-carlo fashion the angle changes for x and z planes. The
!++  angle change is proportional to SQRT(t) and 1/p, as expected.
      call tetat(t,p,dxp,dzp)

!++  Apply angle changes
      xp=xp+dxp
      zp=zp+dzp

!++  Treat single-diffractive scattering.
      if(inter.eq.4) then

!++ added update for s
        s=(zlm-rlen)+zlm1                                                !hr09
        xpsd=dxp
        zpsd=dzp
        psd=p1
!
!++  Add this code to get the momentum transfer also in the calling
!++  routine...
        dpop=(p-p0)/p0

      end if

!++  Calculate the remaining interaction length and close the iteration
!++  loop.
      if(dowrite_impact) then
! write Coll_Scatter.dat for complete scattering histories. 
! Includes changes in angle from both
         write(3998,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e18.10))')            & 
     &   icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef
      end if
      rlen=rlen-zlm1
      goto 10

      end

!>
!! jaw0(s,nabs)
!! ???
!!
!!++  Input:   ZLM is interaction length
!!++           MAT is choice of material
!!
!!++  Output:  nabs = 1   Particle is absorped
!!++           nabs = 4   Single-diffractive scattering
!!++           dpop       Adjusted for momentum loss (dE/dx)
!!++           s          Exit longitudinal position
!!
!!++  Physics:  If monte carlo interaction length greater than input
!!++            interaction length, then use input interaction length
!!++            Is that justified???
!!
!!     nabs=1....absorption
!!
!<
      subroutine jaw0(s,nabs)
      implicit none
+if crlibm
+ca crlibco
+ei
!
+ca interac
      integer nabs,inter,ichoix,icoll,iturn,ipart
      double precision p,rlen,s,t,gettran,dxp,dzp,p1
      real rndm4
!...cne=1/(sqrt(b))
!...dpodx=dE/(dx*c)
      p=p0/(1.d0-dpop)
      nabs=0
      if(mat.eq.nmat) then
!
!++  Collimator treated as black absorber
!
        nabs=1
        s=0d0
        return
      else if(mat.eq.nmat-1) then
!
!++  Collimator treated as drift
!
        s=zlm
        x=x+s*xp
        z=z+s*zp
        return
      end if
!
!++  Initialize the interaction length to input interaction length
!
      rlen=zlm
!
!++  Do a step for a point-like interaction. This is a loop with
!++  label 10!!!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!++  Get monte-carlo interaction length.
!
+if crlibm
10    zlm1=(-1d0*xintl(mat))*log_rn(dble(rndm4()))                       !hr09
+ei
+if .not.crlibm
10    zlm1=(-1d0*xintl(mat))*log(dble(rndm4()))
+ei
!
      if(zlm1.gt.rlen) then
!
!++  If the monte-carlo interaction length is shorter than the
!++  remaining collimator length, then put it to the remaining
!++  length, do multiple coulomb scattering and return.
!++  LAST STEP IN ITERATION LOOP
!
       zlm1=rlen
       call mcs(s)
       s=(zlm-rlen)+s                                                    !hr09
       p=p-dpodx(mat)*s
       dpop=1.d0-p0/p
       return
      end if
!
!++  Otherwise do multi-coulomb scattering.
!++  REGULAR STEP IN ITERATION LOOP
!
      call mcs(s)
!
!++  Check if particle is outside of collimator (X.LT.0) after
!++  MCS. If yes, calculate output longitudinal position (s),
!++  reduce momentum (output as dpop) and return.
!++  PARTICLE LEFT COLLIMATOR BEFORE ITS END.
!
      if(x.le.0.d0) then
       s=(zlm-rlen)+s                                                    !hr09
       p=p-dpodx(mat)*s
       dpop=1.d0-p0/p
       return
      end if
!
!++  Check whether particle is absorbed. If yes, calculate output
!++  longitudinal position (s), reduce momentum (output as dpop)
!++  and return.
!++  PARTICLE WAS ABSORPED INSIDE COLLIMATOR DURING MCS.
!
      inter=ichoix(mat)
      if(inter.eq.1) then
       nabs=1
       s=(zlm-rlen)+zlm1                                                 !hr09
       p=p-dpodx(mat)*s
       dpop=1.d0-p0/p
! write Coll_Scatter.dat for complete scattering histories
       write(3998,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e14.6))')               &  
     & icoll,iturn,ipart,nabs,-1.d0,0.d0,0.d0
       return
      end if
!
!++  Now treat the other types of interaction, as determined by ICHOIX:
!
!++      Nuclear-Elastic:          inter = 2
!++      pp Elastic:               inter = 3
!++      Single-Diffractive:       inter = 4    (changes momentum p)
!++      Coulomb:                  inter = 5
!
!++  As the single-diffractive interaction changes the momentum, save
!++  input momentum in p1.
!
      p1 = p
!
!++  Gettran returns some monte carlo number, that, as I believe, gives
!++  the rms transverse momentum transfer.
!
      t = gettran(inter,mat,p)
!
!++  Tetat calculates from the rms transverse momentum transfer in
!++  monte-carlo fashion the angle changes for x and z planes. The
!++  angle change is proportional to SQRT(t) and 1/p, as expected.
!
      call tetat(t,p,dxp,dzp)
!
!++  Apply angle changes
!
      xp=xp+dxp
      zp=zp+dzp
!
!++  Treat single-diffractive scattering.
!
      if(inter.eq.4) then
        nabs=4
        xpsd=dxp
        zpsd=dzp
        psd=p1
      end if
!
!++  Calculate the remaining interaction length and close the iteration
!++  loop.
!
      rlen=rlen-zlm1
      goto 10
!
      end

!>
!! mcs(s)
!!++  Input:   zlm1   Monte-carlo interaction length
!!
!!++  Output:  s      Longitudinal position
!!++           p0     Reference momentum
!!++           dpop   Relative momentum offset
!!
!!     collimator: x>0 and y<zlm1
!<
      subroutine mcs(s)
      implicit none
+if crlibm
+ca crlibco
+ei
!      save h,dh,bn
+ca interac
      double precision h,dh,theta,rlen0,rlen,ae,be,bn0,s
      double precision radl_mat,rad_len ! Claudia 2013 added variables


!   bn=sqrt(3)/(number of sigmas for s-determination(=4))
      data h/.001d0/dh/.0001d0/bn0/.4330127019d0/

      radl_mat=radl(mat)
      theta=13.6d-3/(p0*(1.d0+dpop))      !Claudia added log part
      rad_len=radl(mat)                    !Claudia

      x=(x/theta)/radl(mat)                                              !hr09
      xp=xp/theta
      z=(z/theta)/radl(mat)                                              !hr09
      zp=zp/theta
      rlen0=zlm1/radl(mat)
      rlen=rlen0
10    ae=bn0*x
      be=bn0*xp
      call soln3(ae,be,dh,rlen,s)
      if(s.lt.h) s=h
      call scamcs(x,xp,s,radl_mat)
      if(x.le.0.d0) then
       s=(rlen0-rlen)+s                                                  !hr09
       goto 20
      end if
      if(s+dh.ge.rlen) then
       s=rlen0
       goto 20
      end if
      rlen=rlen-s
      goto 10
20    call scamcs(z,zp,s,radl_mat)
      s=s*radl(mat)
      x=(x*theta)*radl(mat)                                              !hr09
      xp=xp*theta
      z=(z*theta)*radl(mat)                                              !hr09
      zp=zp*theta
      end

!>
!! scamcs(xx,xxp,s,radl_mat)
!! ???
!<
      subroutine scamcs(xx,xxp,s,radl_mat)
      implicit none
+if crlibm
+ca crlibco
+ei
      double precision v1,v2,r2,a,z1,z2,ss,s,xx,xxp,x0,xp0
      double precision radl_mat
      real rndm4
      x0=xx
      xp0=xxp
5     v1=2d0*dble(rndm4())-1d0
      v2=2d0*dble(rndm4())-1d0
      r2=v1**2+v2**2                                                     !hr09
      if(r2.ge.1.d0) goto 5
+if crlibm
      a=sqrt((-2.d0*log_rn(r2))/r2)                                      !hr09
+ei
+if .not.crlibm
      a=sqrt((-2.d0*log(r2))/r2)                                         !hr09
+ei
      z1=v1*a
      z2=v2*a
      ss=sqrt(s)    
+if crlibm
      xx=x0+s*(xp0+(.5d0*ss)*(1+0.038*log_rn(s))*(z2+z1*.577350269d0)) !Claudia: added logarithmic part in mcs formula                                                     !hr09
      xxp=xp0+ss*z2*(1+0.038*log_rn(s))  

+ei
+if .not.crlibm
      xx=x0+s*(xp0+(.5d0*ss)*(1+0.038*log(s))*(z2+z1*.577350269d0)) !Claudia: added logarithmic part in mcs formula                                                     !hr09
      xxp=xp0+ss*z2*(1+0.038*log(s))  
+ei
      end

!>
!! soln3(a,b,dh,smax,s)
!! ???
!<
      subroutine soln3(a,b,dh,smax,s)
      implicit none
+if crlibm
+ca crlibco
+ei
      double precision b,a,s,smax,c,dh
      if(b.eq.0.d0) then
       s=a**0.6666666666666667d0
!      s=a**(2.d0/3.d0)
       if(s.gt.smax) s=smax
       return
      end if
      if(a.eq.0.d0) then
       if(b.gt.0.d0) then
         s=b**2
       else
         s=0.d0
       end if
       if(s.gt.smax) s=smax
       return
      end if
      if(b.gt.0.d0) then
       if(smax**3.le.(a+b*smax)**2) then
        s=smax
        return
       else
        s=smax*.5d0
        call iterat(a,b,dh,s)
       end if
      else
       c=(-1d0*a)/b
       if(smax.lt.c) then
        if(smax**3.le.(a+b*smax)**2) then
         s=smax
         return
        else
         s=smax*.5d0
         call iterat(a,b,dh,s)
        end if
       else
        s=c*.5d0
        call iterat(a,b,dh,s)
       end if
      end if
      end


      subroutine iterat(a,b,dh,s)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      double precision ds,s,a,b,dh

      ds=s
10    ds=ds*.5d0
      if(s**3.lt.(a+b*s)**2) then
        s=s+ds
      else
        s=s-ds
      end if
      if(ds.lt.dh) then
        return
      else
        goto 10
      end if
      end

!>
!! get_dpodx(p,mat_i)
!! calculate mean ionization energy loss according to Bethe-Bloch
!<
      function get_dpodx(p,mat_i)          !Claudia
      implicit none
      integer mat
+ca collMatNum
      common/materia/mat
      double precision anuc,zatom,rho,emr,exenergy
      double precision PE,me,mp,K,gamma_p
      common/mater/anuc(nmat),zatom(nmat),rho(nmat),emr(nmat)
      common/meanexen/exenergy(nmat)
      double precision beta_p,gamma_s,beta_s,me2,mp2,T,part_1,part_2,   &
     &I_s,delta
      parameter(me=0.510998910e-3,mp=938.272013e-3,K=0.307075)
      double precision p
      integer mat_i
      double precision dpodx,get_dpodx   
+if crlibm
+ca crlibco
+ei

      mp2=mp**2
      me2=me**2
      beta_p=1.
      gamma_p=p/mp
      beta_s=beta_p**2
      gamma_s=gamma_p**2
      T=(2*me*beta_s*gamma_s)/(1+(2*gamma_p*me/mp)+me2/mp2)
      PE=dsqrt(rho(mat_i)*zatom(mat_i)/anuc(mat_i))*28.816e-9
      I_s=exenergy(mat_i)**2
      part_1=K*zatom(mat_i)/(anuc(mat_i)*beta_s)
+if crlibm 
      delta=log_rn(PE/exenergy(mat_i))+log_rn(beta_p*gamma_p)-0.5
      part_2=0.5*log_rn((2*me*beta_s*gamma_s*T)/I_s)
+ei                   
+if .not.crlibm 
      delta=log(PE/exenergy(mat_i))+log(beta_p*gamma_p)-0.5
      part_2=0.5*log((2*me*beta_s*gamma_s*T)/I_s)
+ei
      get_dpodx = part_1*(part_2-beta_s-delta)*rho(mat_i)*1.e-1
      return
      end

!>
!! CalcElectronDensity(AtomicNumber, Density, AtomicMass)
!! Function to calculate the electron density in a material
!! Should give the number per cubic meter
!<
      function CalcElectronDensity(AtomicNumber, Density, AtomicMass)
      implicit none
      double precision AtomicNumber, Density, AtomicMass
      double precision Avogadro
      double precision CalcElectronDensity
      double precision PartA, PartB
      parameter (Avogadro = 6.022140857e23)
      PartA = AtomicNumber * Avogadro * Density
      !1e-6 factor converts to n/m^-3
      PartB = AtomicMass * 1e-6
      CalcElectronDensity = PartA/PartB
      return
      end

!>
!! CalcPlasmaEnergy(ElectronDensity)
!! Function to calculate the plasma energy in a material
!! CalculatePlasmaEnergy = (PlanckConstantBar * sqrt((ElectronDensity *(ElectronCharge**2)) / (ElectronMass * FreeSpacePermittivity)))/ElectronCharge*eV;
!<
      function CalcPlasmaEnergy(ElectronDensity)
      implicit none
      double precision ElectronDensity
      double precision CalcPlasmaEnergy
      double precision sqrtAB,PartA,PartB,FSPC2

      !Values from the 2016 PDG
      double precision PlanckConstantBar,ElectronCharge,ElectronMass
      double precision ElectronCharge2
      double precision FreeSpacePermittivity,FreeSpacePermeability
      double precision SpeedOfLight,SpeedOfLight2

      parameter (PlanckConstantBar = 1.054571800e-34)
      parameter (ElectronCharge = 1.6021766208e-19)
      parameter (ElectronCharge2 = ElectronCharge*ElectronCharge)
      parameter (ElectronMass = 9.10938356e-31)
      parameter (SpeedOfLight = 299792458)
      parameter (SpeedOfLight2 = SpeedOfLight*SpeedOfLight)

      parameter (FreeSpacePermeability = 16.0e-7*atan(1.0)) ! Henry per meter
      parameter (FSPC2 = FreeSpacePermeability*SpeedOfLight2)
      parameter (FreeSpacePermittivity = 1.0/FSPC2)
      parameter (PartB = ElectronMass * FreeSpacePermittivity)

      PartA = ElectronDensity * ElectronCharge2

      sqrtAB = sqrt(PartA/PartB)
      CalcPlasmaEnergy=PlanckConstantBar*sqrtAB/ElectronCharge*1e-9
      return
      end

!>
!! CALC_ION_LOSS(IS,PC,DZ,EnLo)
!! subroutine for the calculazion of the energy loss by ionization
!! Either mean energy loss from Bethe-Bloch, or higher energy loss, according to finite probability from cross section
!! written by DM for crystals, introduced in main code by RB
!<
      SUBROUTINE CALC_ION_LOSS(IS,PC,DZ,EnLo)

! IS material ID
! PC momentum in GeV
! DZ length traversed in material (meters)
! EnLo energy loss in GeV/meter

      IMPLICIT none
      integer IS
+ca collMatNum
      double precision PC,DZ,EnLo,exenergy,exEn
      double precision k,re,me,mp !Daniele: parameters for dE/dX calculation (const,electron radius,el. mass, prot.mass)
      double precision enr,mom,betar,gammar,bgr !Daniele: energy,momentum,beta relativistic, gamma relativistic
      double precision Tmax,plen !Daniele: maximum energy tranfer in single collision, plasma energy (see pdg)
      double precision thl,Tt,cs_tail,prob_tail
      double precision ranc
      REAL*4 RNDM4
      double precision anuc,zatom,rho,emr


      common/meanexen/exenergy(nmat)

      common/mater/anuc(nmat),zatom(nmat),rho(nmat),emr(nmat)
!      common/betheBl/enr,mom,gammar,betar,bgr,exEn,Tmax,plen

      data k/0.307075/      !constant in front bethe-bloch [MeV g^-1 cm^2]
      data re/2.818d-15/  !electron radius [m]
      data me/0.510998910/ !electron mass [MeV/c^2]
      data mp/938.272013/ !proton mass [MeV/c^2]

+if crlibm
+ca crlibco
+ei

      mom=PC*1.0d3              ! [GeV/c] -> [MeV/c]
      enr=(mom*mom+mp*mp)**0.5  ! [MeV]
      gammar=enr/mp
      betar=mom/enr
      bgr=betar*gammar
      
      ! mean excitation energy - convert to MeV
      exEn=exenergy(IS)*1.0d3

      ! Tmax is max energy loss from kinematics
      Tmax=(2.0d0*me*bgr**2)/(1.0d0+2*gammar*me/mp+(me/mp)**2) ![MeV]

      ! plasma energy - see PDG 2010 table 27.1
      plen=((rho(IS)*zatom(IS)/anuc(IS))**0.5)*28.816d-6 ![MeV]
      
      ! calculate threshold energy
      ! Above this threshold, the cross section for high energy loss is calculated and then a random number is generated to determine if tail energy loss should be applied, or only mean from Bethe-Bloch
      ! below threshold, only the standard bethe-bloch is used (all particles get average energy loss)

      ! thl is 2* width of landau distribution (as in fig 27.7 in PDG 2010). See Alfredo's presentation for derivation
      thl= 4.0d0*k*zatom(IS)*DZ*100.0d0*rho(IS)/(anuc(IS)*betar**2) ![MeV]
!     write(3456,*) thl     ! should typically be >0.06MeV for approximations to be valid - check!

+if crlibm
      ! Bethe Bloch mean energy loss
      EnLo=((k*zatom(IS))/(anuc(IS)*betar**2))*
     +     (0.5*log_rn((2.0d0*me*bgr*bgr*Tmax)/(exEn*exEn))
     +     -betar**2.0-log_rn(plen/exEn)-log_rn(bgr)+0.5);

      EnLo=EnLo*rho(IS)*0.1*DZ  ![GeV]

      ! threshold Tt is bethe bloch + 2*width of Landau distribution
      Tt=EnLo*1000.0d0+thl      ![MeV]

       ! cross section - see Alfredo's presentation for derivation
       cs_tail=((k*zatom(IS))/(anuc(IS)*betar**2))*
     + ((0.5*((1.0d0/Tt)-(1.0d0/Tmax)))-
     + (log_rn(Tmax/Tt)*(betar**2)/(2.0d0*Tmax))+
     + ((Tmax-Tt)/(4.0d0*(gammar**2)*(mp**2))))

       ! probability of being in tail: cross section * density * path length
       prob_tail=cs_tail*rho(IS)*DZ*100.0d0;

       ranc=dble(rndm4())

       ! determine based on random number if tail energy loss occurs.
       if(ranc.lt.prob_tail)then
         EnLo=((k*zatom(IS))/(anuc(IS)*betar**2))*
     +   (0.5*log_rn((2.0d0*me*bgr*bgr*Tmax)/(exEn*exEn))
     +   -betar**2.0-log_rn(plen/exEn)-log_rn(bgr)+0.5+
     +   (TMax**2)/(8.0d0*(gammar**2)*(mp**2)));


+ei                   
+if .not.crlibm 
      ! Bethe Bloch mean energy loss
      EnLo=((k*zatom(IS))/(anuc(IS)*betar**2))*
     +     (0.5*log((2.0d0*me*bgr*bgr*Tmax)/(exEn*exEn))
     +     -betar**2.0-log(plen/exEn)-log(bgr)+0.5);

      EnLo=EnLo*rho(IS)*0.1*DZ  ![GeV]

      ! threshold Tt is bethe bloch + 2*width of Landau distribution
      Tt=EnLo*1000.0d0+thl      ![MeV]

       ! cross section - see Alfredo's presentation for derivation
       cs_tail=((k*zatom(IS))/(anuc(IS)*betar**2))*
     + ((0.5*((1.0d0/Tt)-(1.0d0/Tmax)))-
     + (log(Tmax/Tt)*(betar**2)/(2.0d0*Tmax))+
     + ((Tmax-Tt)/(4.0d0*(gammar**2)*(mp**2))))

       ! probability of being in tail: cross section * density * path length
       prob_tail=cs_tail*rho(IS)*DZ*100.0d0;

       ranc=dble(rndm4())

       ! determine based on random number if tail energy loss occurs.
       if(ranc.lt.prob_tail)then
         EnLo=((k*zatom(IS))/(anuc(IS)*betar**2))*
     +   (0.5*log((2.0d0*me*bgr*bgr*Tmax)/(exEn*exEn))
     +   -betar**2.0-log(plen/exEn)-log(bgr)+0.5+
     +   (TMax**2)/(8.0d0*(gammar**2)*(mp**2)));

+ei


         EnLo=EnLo*rho(IS)*0.1 ![GeV/m]

       else 
          ! if tial energy loss does not occur, just use the standard Bethe Bloch
         EnLo=EnLo/DZ  ![GeV/m]
       endif
     
c      write(*,*)cs_tail,prob_tail,ranc,EnLo*DZ

      RETURN
      END

      subroutine makedis(mynp, myalphax, myalphay, mybetax, mybetay,    &
     &myemitx0, myemity0, myenom, mynex, mdex, myney, mdey,             &
     &myx, myxp, myy, myyp, myp, mys)
!
!  Generate distribution
!
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
!
+ca collpara
+ca dbmkdist
      double precision pi
!
      save
!-----------------------------------------------------------------------
!++  Generate particle distribution
!
!
!++  Generate random distribution, assuming optical parameters at IP1
!
!
!++  Calculate the gammas
!
+if crlibm
      pi=4d0*atan_rn(1d0)
+ei
+if .not.crlibm
      pi=4d0*atan(1d0)
+ei
      mygammax = (1d0+myalphax**2)/mybetax
      mygammay = (1d0+myalphay**2)/mybetay
!++TW 11/07 reset j, helps if subroutine is called twice 
! was done during try to reset distribution, still needed
! will this subroutine ever called twice? 
      j = 0
!
!
!++  Number of points and generate distribution
!
      write(lout,*)
      write(lout,*) 'Generation of particle distribution Version 1'
      write(lout,*)
      write(lout,*) 'This routine generates particles in phase space'
      write(lout,*) 'X/XP and Y/YP ellipses, as defined in the input'
      write(lout,*) 'parameters. Distribution is flat in the band.'
      write(lout,*) 'X and Y are fully uncorrelated.'
      write(lout,*)
!
      write(outlun,*)
      write(outlun,*) 'Generation of particle distribution Version 1'
      write(outlun,*)
      write(outlun,*) 'This routine generates particles in phase space'
      write(outlun,*) 'X/XP and Y/YP ellipses, as defined in the input'
      write(outlun,*) 'parameters. Distribution is flat in the band.'
      write(outlun,*) 'X and Y are fully uncorrelated.'
      write(outlun,*)
      write(outlun,*) 'INFO>  Number of particles   = ', mynp
      write(outlun,*) 'INFO>  Av number of x sigmas = ', mynex
      write(outlun,*) 'INFO>  +- spread in x sigmas = ', mdex
      write(outlun,*) 'INFO>  Av number of y sigmas = ', myney
      write(outlun,*) 'INFO>  +- spread in y sigmas = ', mdey
      write(outlun,*) 'INFO>  Nominal beam energy   = ', myenom
      write(outlun,*) 'INFO>  Sigma_x0 = ', sqrt(mybetax*myemitx0)
      write(outlun,*) 'INFO>  Sigma_y0 = ', sqrt(mybetay*myemity0)
      write(outlun,*) 'INFO>  Beta x   = ', mybetax
      write(outlun,*) 'INFO>  Beta y   = ', mybetay
      write(outlun,*) 'INFO>  Alpha x  = ', myalphax
      write(outlun,*) 'INFO>  Alpha y  = ', myalphay
      write(outlun,*)
!
      do while (j.lt.mynp)
!
        j = j + 1
        myemitx = myemitx0*(mynex + ((2d0*dble(rndm4()-0.5))*mdex) )**2  !hr09
        xsigmax = sqrt(mybetax*myemitx)
+if crlibm
        myx(j)   = xsigmax * sin_rn((2d0*pi)*dble(rndm4()))              !hr09
+ei
+if .not.crlibm
        myx(j)   = xsigmax * sin((2d0*pi)*dble(rndm4()))                 !hr09
+ei
        if (rndm4().gt.0.5) then
          myxp(j)  = sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-        &
     &(myalphax*myx(j))/mybetax                                          !hr09
        else
          myxp(j)  = -1d0*sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-   &!hr09
     &(myalphax*myx(j))/mybetax                                          !hr09
        endif
!
        myemity = myemity0*(myney + ((2d0*dble(rndm4()-0.5))*mdey) )**2  !hr09
        ysigmay = sqrt(mybetay*myemity)
+if crlibm
        myy(j)   = ysigmay * sin_rn((2d0*pi)*dble(rndm4()))              !hr09
+ei
+if .not.crlibm
        myy(j)   = ysigmay * sin((2d0*pi)*dble(rndm4()))                 !hr09
+ei
        if (rndm4().gt.0.5) then
          myyp(j)  = sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-        &
     &(myalphay*myy(j))/mybetay                                          !hr09
        else
          myyp(j)  = -1d0*sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-   &
     &(myalphay*myy(j))/mybetay                                          !hr09
        endif
!
!APRIL2005 TEST FOR FATS FLAG
        myp(j)   = myenom
!        if(j.eq.1) then
!          myp(j)   = myenom*(1-0.01)
!!       do j=2,mynp
!        else
!          myp(j) = myp(1) + (j-1)*2d0*0.01*myenom/(mynp-1)
!        endif
!APRIL2005 END OF TEST SECTION
        mys(j)   = 0d0
!
!++  Dangerous stuff, just for the moment
!
        if (cut_input) then
          if ( (.not. (myy(j).lt.-.008d-3 .and. myyp(j).lt.0.1d-3 .and. &
     &myyp(j).gt.0d0) ) .and.                                           &
     &(.not. (myy(j).gt..008d-3 .and. myyp(j).gt.-0.1d-3 .and.          &
     &myyp(j).lt.0d0) ) ) then
            j = j - 1
          endif
        endif
!
      end do
!
      return
      end
!
!========================================================================
!
! SR, 08-05-2005: Add the finite beam size in the othe dimension
      subroutine makedis_st(mynp, myalphax, myalphay, mybetax, mybetay, &
     &     myemitx0, myemity0, myenom, mynex, mdex, myney, mdey,        &
     &     myx, myxp, myy, myyp, myp, mys)

!     Uses the old routine 'MAKEDIS' for the halo plane and adds the
!     transverse beam size in the other plane (matched distrubutions
!     are generated starting from thetwiss functions).
!     If 'mynex' and 'myney' are BOTH set to zero, nominal bunches
!     centred in the aperture centre are generated. (SR, 08-05-2005)
!
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
!
+ca collpara
+ca dbmkdist
      double precision pi
!
      double precision iix, iiy, phix, phiy
!
      save
!
!-----------------------------------------------------------------------
!++  Generate particle distribution
!
!
!++  Generate random distribution, assuming optical parameters at IP1
!
!++  Calculate the gammas
!
      write(lout,*) '  New routine to add the finite beam size in the'
      write(lout,*) '  other dimension (SR, 08-06-2005).'

      pi=4d0*atan(1d0)
!
      mygammax = (1d0+myalphax**2)/mybetax
      mygammay = (1d0+myalphay**2)/mybetay
!
      do j=1, mynp
         if ((mynex.gt.0d0).and.(myney.eq.0d0)) then
            myemitx = myemitx0*(mynex+((2d0*dble(rndm4()-0.5))*mdex))**2 !hr09
            xsigmax = sqrt(mybetax*myemitx)
+if crlibm
            myx(j)   = xsigmax * sin_rn((2d0*pi)*dble(rndm4()))          !hr09
+ei
+if .not.crlibm
            myx(j)   = xsigmax * sin((2d0*pi)*dble(rndm4()))             !hr09
+ei
            if (rndm4().gt.0.5) then
              myxp(j) = sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-     &
     &              (myalphax*myx(j))/mybetax                            !hr09
            else
              myxp(j) = -1d0*sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-&
     &              (myalphax*myx(j))/mybetax                            !hr09
            endif
!
            phiy = (2d0*pi)*dble(rndm4())                                !hr09
!
+if crlibm
            iiy = (-1d0*myemity0) * log_rn( dble(rndm4()) )              !hr09
!
            myy(j) = sqrt((2d0*iiy)*mybetay) * cos_rn(phiy)              !hr09
            myyp(j) = (-1d0*sqrt((2d0*iiy)/mybetay)) * (sin_rn(phiy) +  &!hr09
     &           myalphay * cos_rn(phiy))                                !hr09
+ei
+if .not.crlibm
            iiy = (-1d0*myemity0) * log( dble(rndm4()) )                 !hr09
!
            myy(j) = sqrt((2d0*iiy)*mybetay) * cos(phiy)                 !hr09
            myyp(j) = (-1d0*sqrt((2d0*iiy)/mybetay)) * (sin(phiy) +     &!hr09
     &           myalphay * cos(phiy))                                   !hr09
+ei
         elseif ( mynex.eq.0d0.and.myney.gt.0d0 ) then                   !hr09
            myemity = myemity0*(myney+((2d0*dble(rndm4()-0.5))*mdey))**2 !hr09
            ysigmay = sqrt(mybetay*myemity)
+if crlibm
            myy(j)   = ysigmay * sin_rn((2d0*pi)*dble(rndm4()))          !hr09
+ei
+if .not.crlibm
            myy(j)   = ysigmay * sin((2d0*pi)*dble(rndm4()))             !hr09
+ei
            if (rndm4().gt.0.5) then
              myyp(j) = sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-     &!hr09
     &              (myalphay*myy(j))/mybetay                            !hr09
            else
              myyp(j) = -1d0*sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-&!hr09
     &              (myalphay*myy(j))/mybetay                            !hr09
            endif
!
            phix = (2d0*pi)*dble(rndm4())                                !hr09
+if crlibm
            iix = (-1d0* myemitx0) * log_rn( dble(rndm4()) )             !hr09
!
            myx(j) = sqrt((2d0*iix)*mybetax) * cos_rn(phix)              !hr09
            myxp(j) = (-1d0*sqrt((2d0*iix)/mybetax)) * (sin_rn(phix) +  &!hr09
     &           myalphax * cos_rn(phix))                                !hr09
+ei
+if .not.crlibm
            iix = (-1d0* myemitx0) * log( dble(rndm4()) )                !hr09
!
            myx(j) = sqrt((2d0*iix)*mybetax) * cos(phix)                 !hr09
            myxp(j) = (-1d0*sqrt((2d0*iix)/mybetax)) * (sin(phix) +     &!hr09
     &           myalphax * cos(phix))                                   !hr09
+ei
         elseif ( mynex.eq.0d0.and.myney.eq.0d0 ) then                   !hr09
            phix = (2d0*pi)*dble(rndm4())                                !hr09
+if crlibm
            iix = (-1d0* myemitx0) * log_rn( dble(rndm4()) )             !hr09
!
            myx(j) = sqrt((2d0*iix)*mybetax) * cos_rn(phix)              !hr09
            myxp(j) = (-1d0*sqrt((2d0*iix)/mybetax)) * (sin_rn(phix) +  &!hr09
     &           myalphax * cos_rn(phix))                                !hr09
+ei
+if .not.crlibm
            iix = (-1d0*myemitx0) * log( dble(rndm4()) )                 !hr09
!
            myx(j) = sqrt((2d0*iix)*mybetax) * cos(phix)                 !hr09
            myxp(j) = (-1d0*sqrt((2d0*iix)/mybetax)) * (sin(phix) +     &!hr09
     &           myalphax * cos(phix))                                   !hr09
+ei
            phiy = (2d0*pi)*dble(rndm4())                                !hr09
+if crlibm
            iiy = (-1d0*myemity0) * log_rn( dble(rndm4()) )              !hr09
            myy(j) = sqrt((2d0*iiy)*mybetay) * cos_rn(phiy)              !hr09
            myyp(j) = (-1d0*sqrt((2d0*iiy)/mybetay)) * (sin_rn(phiy) +  &!hr09
     &           myalphay * cos_rn(phiy))                                !hr09
+ei
+if .not.crlibm
            iiy = (-1d0*myemity0) * log( dble(rndm4()) )                 !hr09
!
            myy(j) = sqrt((2d0*iiy)*mybetay) * cos(phiy)                 !hr09
            myyp(j) = (-1d0*sqrt((2d0*iiy)/mybetay)) * (sin(phiy) +     &!hr09
     &           myalphay * cos(phiy))                                   !hr09
+ei
         else
            write(lout,*) "Error - beam parameters not correctly set!"
         endif
!
         myp(j)   = myenom
         mys(j)   = 0d0
!
      end do
!
      return
      end

!========================================================================
!
!     RB: new routine to sample part of matched phase ellipse which is outside 
!     the cut of the jaws
!     Assuming cut of the jaw at mynex for hor plane.
!     largest amplitude outside of jaw is mynex + mdex.  Analog for vertical plane.

!     same routine as makedis_st, but rejection sampling to get
!     only particles hitting the collimator on the same turn. 

!     Treat as a pencil beam in main routine. 

      subroutine makedis_coll(mynp,myalphax, myalphay, mybetax, mybetay, &
     &     myemitx0, myemity0, myenom, mynex, mdex, myney, mdey,        &
     &     myx, myxp, myy, myyp, myp, mys)
 
      implicit none
+ca crcoall
+ca collpara
+ca dbmkdist

      double precision pi, iix, iiy, phix,phiy,cutoff
      
      save
!
!-----------------------------------------------------------------------
!++  Generate particle distribution
!
!++  Calculate the gammas
!
      write(lout,*) '  RB 2013: new pencil beam routine'
      pi=4d0*atan(1d0)
!
      mygammax = (1d0+myalphax**2)/mybetax
      mygammay = (1d0+myalphay**2)/mybetay

! calcualte cutoff in x or y from the collimator jaws. 
      if ((mynex.gt.0d0).and.(myney.eq.0d0)) then
         cutoff=mynex*sqrt(mybetax*myemitx0)
      else
         cutoff=myney*sqrt(mybetay*myemity0)
      endif

!
      do j=1, mynp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if ((mynex.gt.0d0).and.(myney.eq.0d0)) then  ! halo in x
 887        continue
            myemitx = myemitx0*(mynex+(dble(rndm4())*mdex))**2  
            xsigmax = sqrt(mybetax*myemitx)
            myx(j)   = xsigmax * sin((2d0*pi)*dble(rndm4()))       
            if (abs(myx(j)).lt.cutoff) goto 887
            if (rndm4().gt.0.5) then
              myxp(j) = sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-     &
     &              (myalphax*myx(j))/mybetax                             
            else
              myxp(j) = -1d0*sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-&
     &              (myalphax*myx(j))/mybetax                             
            endif
            phiy = (2d0*pi)*dble(rndm4())                                 
            iiy = (-1d0*myemity0) * log( dble(rndm4()) )                  
            myy(j) = sqrt((2d0*iiy)*mybetay) * cos(phiy)                  
            myyp(j) = (-1d0*sqrt((2d0*iiy)/mybetay)) * (sin(phiy) +     & 
     &           myalphay * cos(phiy))                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         elseif ( mynex.eq.0d0.and.myney.gt.0d0 ) then  ! halo in y
 886        continue
            myemity = myemity0*(myney+(dble(rndm4())*mdey))**2  
            ysigmay = sqrt(mybetay*myemity)
            myy(j)   = ysigmay * sin((2d0*pi)*dble(rndm4()))              
            if (abs(myy(j)).lt.cutoff) goto 886
            if (rndm4().gt.0.5) then
              myyp(j) = sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-     & 
     &              (myalphay*myy(j))/mybetay                             
            else
              myyp(j) = -1d0*sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-& 
     &              (myalphay*myy(j))/mybetay                             
            endif
            phix = (2d0*pi)*dble(rndm4())                                 
            iix = (-1d0* myemitx0) * log( dble(rndm4()) )                 
            myx(j) = sqrt((2d0*iix)*mybetax) * cos(phix)                  
            myxp(j) = (-1d0*sqrt((2d0*iix)/mybetax)) * (sin(phix) +     & 
     &           myalphax * cos(phix))                             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
         elseif ( mynex.eq.0d0.and.myney.eq.0d0 ) then  ! nominal bunches centered in the aperture - can't apply rejection sampling. return with error
            write(lout,*) "Stop in makedis_coll. attempting to use halo type 
     &3 with Gaussian dist. "
            call prror(-1)
c$$$            phix = (2d0*pi)*dble(rndm4())                                 
c$$$            iix = (-1d0*myemitx0) * log( dble(rndm4()) )                  
c$$$            myx(j) = sqrt((2d0*iix)*mybetax) * cos(phix)                  
c$$$            myxp(j) = (-1d0*sqrt((2d0*iix)/mybetax)) * (sin(phix) +     & 
c$$$     &           myalphax * cos(phix))                                    
c$$$            phiy = (2d0*pi)*dble(rndm4())                                 
c$$$            iiy = (-1d0*myemity0) * log( dble(rndm4()) )                  
c$$$
c$$$            myy(j) = sqrt((2d0*iiy)*mybetay) * cos(phiy)                  
c$$$
c$$$            myyp(j) = (-1d0*sqrt((2d0*iiy)/mybetay)) * (sin(phiy) +     & 
c$$$     &           myalphay * cos(phiy))                                    

         else
            write(lout,*) "Error - beam parameters not correctly set!"
         endif
!
         myp(j)   = myenom
         mys(j)   = 0d0
!
      end do
!
      return
      end
!


!
!========================================================================
!
! SR, 09-05-2005: Add the energy spread and the finite bunch length.
!                 Gaussian distributions assumed
      subroutine makedis_de(mynp, myalphax, myalphay, mybetax, mybetay, &
     &     myemitx0, myemity0, myenom, mynex, mdex, myney, mdey,        &
     &     myx, myxp, myy, myyp, myp, mys,                              &
     &     enerror,bunchlength)

!     Uses the old routine 'MAKEDIS' for the halo plane and adds the
!     transverse beam size in the other plane (matched distrubutions
!     are generated starting from thetwiss functions).
!     If 'mynex' and 'myney' are BOTH set to zero, nominal bunches
!     centred in the aperture centre are generated. (SR, 08-05-2005)
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
!
+ca collpara
+ca dbmkdist
      double precision pi
!
      double precision ran_gauss
      double precision iix, iiy, phix, phiy
      double precision enerror, bunchlength
      double precision en_error, bunch_length
!
      double precision long_cut
      double precision a_st, b_st
!
      save
!-----------------------------------------------------------------------
!++  Generate particle distribution
!
!
!++  Generate random distribution, assuming optical parameters at IP1
!
!++  Calculate the gammas
      pi=4d0*atan(1d0)
!
      mygammax = (1d0+myalphax**2)/mybetax
      mygammay = (1d0+myalphay**2)/mybetay

!     Assign bunch length and dp/p depending on the energy
!     Check if the units in metres are correct!
!GRD      if ( myenom.eq.7e6 ) then
!GRD         en_error     = 1.129e-4
!GRD         bunch_length = 7.55e-2
!GRD      elseif ( myenom.eq.4.5e5 ) then
!GRD         en_error     = 3.06e-4
!GRD         bunch_length = 11.24e-2
!GRD      else
      en_error = enerror
      bunch_length = bunchlength
      
!GRD         write(lout,*)"Warning-Energy different from LHC inj or top!"
!GRD         write(lout,*)"  => 7TeV values of dp/p and bunch length used!"
!GRD      endif
!GRD
      write(lout,*) "Generation of bunch with dp/p and length:"
      write(lout,*) "  RMS bunch length  = ", bunch_length
      write(lout,*) "  RMS energy spread = ", en_error

      do j=1, mynp
         if ((mynex.gt.0d0).and.(myney.eq.0d0)) then
            myemitx = myemitx0*(mynex+((2d0*dble(rndm4()-0.5))*mdex))**2 !hr09
            xsigmax = sqrt(mybetax*myemitx)
+if crlibm
            myx(j)   = xsigmax * sin_rn((2d0*pi)*dble(rndm4()))          !hr09
+ei
+if .not.crlibm
            myx(j)   = xsigmax * sin((2d0*pi)*dble(rndm4()))             !hr09
+ei
            if (rndm4().gt.0.5) then
              myxp(j) = sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-     &!hr09
     &              (myalphax*myx(j))/mybetax                            !hr09
            else
              myxp(j) = -1d0*sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-&!hr09
     &              (myalphax*myx(j))/mybetax                            !hr09
            endif
!
            phiy = (2d0*pi)*dble(rndm4())                                !hr09
!
+if crlibm
            iiy = (-1d0*myemity0) * log_rn( dble(rndm4()) )              !hr09
!
            myy(j) = sqrt((2d0*iiy)*mybetay) * cos_rn(phiy)              !hr09
            myyp(j) = (-1d0*sqrt((2d0*iiy)/mybetay)) * (sin_rn(phiy) +  &!hr09
     &           myalphay * cos_rn(phiy))                                !hr09
+ei
+if .not.crlibm
            iiy = (-1d0*myemity0) * log( dble(rndm4()) )                 !hr09
!
            myy(j) = sqrt((2d0*iiy)*mybetay) * cos(phiy)                 !hr09
            myyp(j) = (-1d0*sqrt((2d0*iiy)/mybetay)) * (sin(phiy) +     &!hr09
     &           myalphay * cos(phiy))                                   !hr09
+ei
         elseif ( mynex.eq.0d0.and.myney.gt.0d0 ) then
            myemity = myemity0*(myney+((2d0*dble(rndm4()-0.5))*mdey))**2 !hr09
            ysigmay = sqrt(mybetay*myemity)
+if crlibm
            myy(j)   = ysigmay * sin_rn((2d0*pi)*dble(rndm4()))          !hr09
+ei
+if .not.crlibm
            myy(j)   = ysigmay * sin((2d0*pi)*dble(rndm4()))             !hr09
+ei
            if (rndm4().gt.0.5) then
              myyp(j) = sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-     &!hr09
     &              (myalphay*myy(j))/mybetay                            !hr09
            else
              myyp(j) = -1d0*sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-&!hr09
     &              (myalphay*myy(j))/mybetay                            !hr09
            endif
!
            phix = (2d0*pi)*dble(rndm4())                                !hr09
+if crlibm
            iix = (-1d0*myemitx0) * log_rn( dble(rndm4()) )              !hr09
!
            myx(j) = sqrt((2d0*iix)*mybetax) * cos_rn(phix)              !hr09
            myxp(j) = (-1d0*sqrt((2d0*iix)/mybetax)) * (sin_rn(phix) +  &!hr09
     &           myalphax * cos_rn(phix))                                !hr09
+ei
+if .not.crlibm
            iix = (-1d0*myemitx0) * log( dble(rndm4()) )                 !hr09
!
            myx(j) = sqrt((2d0*iix)*mybetax) * cos(phix)                 !hr09
            myxp(j) = (-1d0*sqrt((2d0*iix)/mybetax)) * (sin(phix) +     &!hr09
     &           myalphax * cos(phix))                                   !hr09
+ei
         elseif ( mynex.eq.0d0.and.myney.eq.0d0 ) then                   !hr09
            phix = (2d0*pi)*dble(rndm4())                                !hr09
+if crlibm
            iix = (-1d0*myemitx0) * log_rn( dble(rndm4()) )              !hr09
!
            myx(j) = sqrt((2d0*iix)*mybetax) * cos_rn(phix)              !hr09
            myxp(j) = (-1d0*sqrt((2d0*iix)/mybetax)) * (sin_rn(phix) +  &!hr09
     &           myalphax * cos_rn(phix))                                !hr09
+ei
+if .not.crlibm
            iix = (-1d0*myemitx0) * log( dble(rndm4()) )                 !hr09
!
            myx(j) = sqrt((2d0*iix)*mybetax) * cos(phix)                 !hr09
            myxp(j) = (-1d0*sqrt((2d0*iix)/mybetax)) * (sin(phix) +     &!hr09
     &           myalphax * cos(phix))                                   !hr09
+ei
            phiy = (2d0*pi)*dble(rndm4())                                !hr09
+if crlibm
            iiy = (-1d0*myemity0) * log_rn( dble(rndm4()) )              !hr09
            myy(j) = sqrt((2d0*iiy)*mybetay) * cos_rn(phiy)              !hr09
            myyp(j) = (-1d0*sqrt((2d0*iiy)/mybetay)) * (sin_rn(phiy) +  &!hr09
     &           myalphay * cos_rn(phiy))                                !hr09
+ei
+if .not.crlibm
            iiy = (-1d0*myemity0) * log( dble(rndm4()) )                 !hr09
!
            myy(j) = sqrt((2d0*iiy)*mybetay) * cos(phiy)                 !hr09
            myyp(j) = (-1d0*sqrt((2d0*iiy)/mybetay)) * (sin(phiy) +     &!hr09
     &           myalphay * cos(phiy))                                   !hr09
+ei
         else
            write(lout,*) "Error - beam parameters not correctly set!"
         endif
!
      end do
! SR, 11-08-2005 For longitudinal phase-space, add a cut at 2 sigma
!++   1st: generate mynpnumbers within the chose cut
      long_cut = 2
      j = 1
      do while (j.le.mynp)
         a_st = ran_gauss(5d0)
         b_st = ran_gauss(5d0)
         do while ((a_st**2+b_st**2).gt.long_cut**2)
            a_st = ran_gauss(5d0)
            b_st = ran_gauss(5d0)
         enddo
         mys(j) = a_st
         myp(j) = b_st
         j = j + 1
      enddo
!++   2nd: give the correct values
      do j=1,mynp
         myp(j) = myenom * (1d0 + myp(j) * en_error)
         mys(j) = bunch_length * mys(j)
      enddo
!
      return
      end
!
!========================================================================
!
      subroutine readdis(filename_dis,mynp,myx,myxp,myy,myyp,myp,mys)
!
!     SR, 09-08-2005
!     Format for the input file:
!               x, y   -> [ m ]
!               xp, yp -> [ rad ]
!               s      -> [ mm ]
!               DE     -> [ MeV ]
!
      implicit none

+ca crcoall
+if crlibm
+ca crlibco
+ei
+ca collpara
+ca dbmkdist

      character*80   filename_dis
      
      logical lopen
      integer stat
      
      save

      write(lout,*) "Reading input bunch from file ", filename_dis

      inquire( unit=53, opened=lopen )
      if (lopen) then
         write(lout,*) "ERROR in subroutine readdis: "//
     &        "FORTRAN Unit 53 was already open!"
         goto 20
      endif
      open(unit=53, file=filename_dis, iostat=stat,
     &     status="OLD",action="read")
      if (stat.ne.0)then
         write(lout,*) "Error in subroutine readdis: "//
     &        "Could not open the file."
         write(lout,*) "Got iostat=",stat
         goto 20
      endif

      do j=1,mynp
         read(53,*,end=10,err=20) myx(j), myxp(j), myy(j), myyp(j),     &
     &        mys(j), myp(j)
      enddo
      
 10   mynp = j - 1
      write(lout,*) "Number of particles read from the file = ",mynp

      close(53)

      return
      
 20   continue

      write(lout,*) "I/O Error on Unit 53 in subroutine readdis"
      call prror(-1)
      
      end

!
!========================================================================
!
      subroutine readdis_norm(filename_dis, 
     &           mynp, myalphax, myalphay, mybetax, mybetay,
     &           myemitx, myemity, myenom, 
     &           myx, myxp, myy, myyp, myp, mys, enerror,
     &           bunchlength)
!     Format for the input file:
!               x, y   -> [ sigma ]
!               xp, yp -> [ sigma ]
!               s      -> [ sigma ]
!               DE     -> [ sigma ]
!
      implicit none

+ca crcoall
+if crlibm
+ca crlibco
+ei
+ca collpara
+ca dbmkdist

+ca parpro
+ca commonmn
+ca common

      character*80   filename_dis
      double precision enerror, bunchlength
      
      logical lopen
      integer stat
      
      double precision normx, normy, normxp, normyp, normp, norms
      double precision myemitz

      write(lout,*) "Reading input bunch from file ", filename_dis

      if (iclo6.eq.0) then
         write(lout,*) "ERROR DETECTED: Incompatible flag           "
         write(lout,*) "in line 2 of the TRACKING block             "
         write(lout,*) "of fort.3 for calculating the closed orbit  "
         write(lout,*) "(iclo6 must not be =0). When using an input "
         write(lout,*) "distribution in normalized coordinates for  "
         write(lout,*) "collimation the closed orbit is needed for a"
         write(lout,*) "correct TAS matrix for coordinate transform."
         call prror(-1)
      endif

      inquire( unit=53, opened=lopen )
      if (lopen) then
         write(lout,*) "ERROR in subroutine readdis: "//
     &        "FORTRAN Unit 53 was already open!"
         goto 20
      endif
      open(unit=53, file=filename_dis, iostat=stat,
     &     status="OLD",action="read")
      if (stat.ne.0)then
         write(lout,*) "Error in subroutine readdis: "//
     &        "Could not open the file."
         write(lout,*) "Got iostat=",stat
         goto 20
      endif

      do j=1,mynp
         read(53,*,end=10,err=20) normx, normxp, normy,
     &     normyp, norms, normp
! A normalized distribution with x,xp,y,yp,z,zp is read and 
! transformed with the TAS matrix T , which is the transformation matrix
! from normalized to physical coordinates it is scaled with the geometric
! emittances in diag matrix S. x = T*S*normx
! units of TAS matrix # m,rad,m,rad,m,1
! The collimation coordinates/units are
! x[m], x'[rad], y[m], y'[rad]$, sig[mm], dE [MeV].

!         write(lout,*) " myenom [MeV]= ",myenom
!         write(lout,*) " myemitx [m]= ",myemitx
!         write(lout,*) " myemity [m]= ",myemity
!         write(lout,*) " bunchlength [mm]= ",bunchlength
!         write(lout,*) " enerror = ",enerror
                  
         !convert bunchlength from [mm] to [m]
         ! enerror is the energy spread
         myemitz  = bunchlength * 0.001d0 * enerror


! scaling the TAS matrix entries of the longitudinal coordinate. tas(ia,j,k)  ia=the particle for which the tas was written

         myx(j)   = 
     &     normx  * sqrt(myemitx)*tas(1,1,1) + 
     &     normxp * sqrt(myemitx)*tas(1,1,2) +
     &     normy  * sqrt(myemity)*tas(1,1,3) +
     &     normyp * sqrt(myemity)*tas(1,1,4) +
     &     norms  * sqrt(myemitz)*tas(1,1,5) +
     &     normp  * sqrt(myemitz)*0.001d0*tas(1,1,6)
         myxp(j)  = 
     &     normx  * sqrt(myemitx)*tas(1,2,1) + 
     &     normxp * sqrt(myemitx)*tas(1,2,2) +
     &     normy  * sqrt(myemity)*tas(1,2,3) +
     &     normyp * sqrt(myemity)*tas(1,2,4) +
     &     norms  * sqrt(myemitz)*tas(1,2,5) +
     &     normp  * sqrt(myemitz)*0.001d0*tas(1,2,6)
         myy(j)   = 
     &     normx  * sqrt(myemitx)*tas(1,3,1) + 
     &     normxp * sqrt(myemitx)*tas(1,3,2) +
     &     normy  * sqrt(myemity)*tas(1,3,3) +
     &     normyp * sqrt(myemity)*tas(1,3,4) +
     &     norms  * sqrt(myemitz)*tas(1,3,5) +
     &     normp  * sqrt(myemitz)*0.001d0*tas(1,3,6)
         myyp(j)  = 
     &     normx  * sqrt(myemitx)*tas(1,4,1) + 
     &     normxp * sqrt(myemitx)*tas(1,4,2) +
     &     normy  * sqrt(myemity)*tas(1,4,3) +
     &     normyp * sqrt(myemity)*tas(1,4,4) +
     &     norms  * sqrt(myemitz)*tas(1,4,5) +
     &     normp  * sqrt(myemitz)*0.001d0*tas(1,4,6)
         mys(j)   = 
     &     normx  * sqrt(myemitx)*tas(1,5,1) + 
     &     normxp * sqrt(myemitx)*tas(1,5,2) +
     &     normy  * sqrt(myemity)*tas(1,5,3) +
     &     normyp * sqrt(myemity)*tas(1,5,4) +
     &     norms  * sqrt(myemitz)*tas(1,5,5) +
     &     normp  * sqrt(myemitz)*0.001d0*tas(1,5,6)
         myp(j)   = 
     &     normx  * sqrt(myemitx)*1000.d0*tas(1,6,1) + 
     &     normxp * sqrt(myemitx)*1000.d0*tas(1,6,2) +
     &     normy  * sqrt(myemity)*1000.d0*tas(1,6,3) +
     &     normyp * sqrt(myemity)*1000.d0*tas(1,6,4) +
     &     norms  * sqrt(myemitz)*1000.d0*tas(1,6,5) +
     &     normp  * sqrt(myemitz)*tas(1,6,6)

! add the momentum
! convert to canonical variables
! dE/E with unit [1] from the closed orbit is added 
!For the 4D coordinates the closed orbit
! will be added by SixTrack itself later on.
         myxp(j)  = myxp(j)*(1.d0+myp(j)+clop6v(3,1))
         myyp(j)  = myyp(j)*(1.d0+myp(j)+clop6v(3,1))
! unit conversion for collimation [m] to [mm]
         mys(j)   = mys(j)*1000.d0
         myp(j)   = myenom*(1.d0+myp(j))

      enddo
      
 10   mynp = j - 1
      write(lout,*) "Number of particles read from the file = ",mynp

      close(53)

      return
      
 20   continue
      write(lout,*) "I/O Error on Unit 53 in subroutine readdis"
      call prror(-1)
      
      end

!
!========================================================================
!
      subroutine makedis_radial(mynp, myalphax, myalphay, mybetax,      &
     &mybetay, myemitx0, myemity0, myenom, nr, ndr,myx, myxp, myy,      &
     &myyp, myp, mys)
!
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
!
+ca collpara
+ca dbmkdist
      double precision pi
!
      save
!-----------------------------------------------------------------------
!++  Generate particle distribution
!
!
!++  Generate random distribution, assuming optical parameters at IP1
!
!++  Calculate the gammas
!
+if crlibm
      pi=4d0*atan_rn(1d0)
+ei
+if .not.crlibm
      pi=4d0*atan(1d0)
+ei
      mygammax = (1d0+myalphax**2)/mybetax
      mygammay = (1d0+myalphay**2)/mybetay
!
!++  Number of points and generate distribution
!
      mynex = nr/sqrt(2d0)
      mdex = ndr/sqrt(2d0)
      myney = nr/sqrt(2d0)
      mdey = ndr/sqrt(2d0)
!
      write(lout,*)
      write(lout,*) 'Generation of particle distribution Version 2'
      write(lout,*)
      write(lout,*) 'This routine generates particles in that are fully'
      write(lout,*) 'correlated between X and Y.'
      write(lout,*)
!
      write(outlun,*)
      write(outlun,*) 'Generation of particle distribution Version 2'
      write(outlun,*)
      write(outlun,*)                                                   &
     &'This routine generates particles in that are fully'
      write(outlun,*) 'correlated between X and Y.'
      write(outlun,*)
      write(outlun,*)
      write(outlun,*) 'INFO>  Number of particles   = ', mynp
      write(outlun,*) 'INFO>  Av number of x sigmas = ', mynex
      write(outlun,*) 'INFO>  +- spread in x sigmas = ', mdex
      write(outlun,*) 'INFO>  Av number of y sigmas = ', myney
      write(outlun,*) 'INFO>  +- spread in y sigmas = ', mdey
      write(outlun,*) 'INFO>  Nominal beam energy   = ', myenom
      write(outlun,*) 'INFO>  Sigma_x0 = ', sqrt(mybetax*myemitx0)
      write(outlun,*) 'INFO>  Sigma_y0 = ', sqrt(mybetay*myemity0)
      write(outlun,*)
!
      do while (j.lt.mynp)
!
        j = j + 1
        myemitx = myemitx0*(mynex + ((2d0*dble(rndm4()-0.5))*mdex) )**2  !hr09
        xsigmax = sqrt(mybetax*myemitx)
+if crlibm
        myx(j)   = xsigmax * sin_rn((2d0*pi)*dble(rndm4()))              !hr09
+ei
+if .not.crlibm
        myx(j)   = xsigmax * sin((2d0*pi)*dble(rndm4()))                 !hr09
+ei
        if (rndm4().gt.0.5) then
          myxp(j)  = sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-        &!hr09
     &(myalphax*myx(j))/mybetax                                          !hr09
        else
          myxp(j)  = -1d0*sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-   &!hr09
     &(myalphax*myx(j))/mybetax                                          !hr09
        endif
!
        myemity = myemity0*(myney + ((2d0*dble(rndm4()-0.5))*mdey) )**2  !hr09
        ysigmay = sqrt(mybetay*myemity)
+if crlibm
        myy(j)   = ysigmay * sin_rn((2d0*pi)*dble(rndm4()))              !hr09
+ei
+if .not.crlibm
        myy(j)   = ysigmay * sin((2d0*pi)*dble(rndm4()))                 !hr09
+ei
        if (rndm4().gt.0.5) then
          myyp(j)  = sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-        &!hr09
     &(myalphay*myy(j))/mybetay                                          !hr09
        else
          myyp(j)  = -1d0*sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-   &!hr09
     &(myalphay*myy(j))/mybetay                                          !hr09
        endif
!
!APRIL2005
        myp(j)   = myenom
!        if(j.eq.1) then
!          myp(j)   = myenom*(1-0.05)
!!       do j=2,mynp
!        else
!          myp(j) = myp(1) + (j-1)*2d0*0.05*myenom/(mynp-1)
!        endif
!APRIL2005
        mys(j)   = 0d0
!
!++  Dangerous stuff, just for the moment
!
!        IF ( (.NOT. (Y(j).LT.-.008e-3 .AND. YP(j).LT.0.1e-3 .AND.
!     1               YP(j).GT.0.0) ) .AND.
!     2       (.NOT. (Y(j).GT..008e-3 .AND. YP(j).GT.-0.1e-3 .AND.
!     3               YP(j).LT.0.0) ) ) THEN
!          J = J - 1
!        ENDIF
!
      end do
!
      return
      end

!>
!! \brief The routine makes an initial Gaussian distribution
!! 
!!     Uses the old routine 'MAKEDIS' for the halo plane and adds the\n
!!     transverse beam size in the other plane (matched distrubutions\n
!!     are generated starting from the twiss functions).\n
!!     If 'mynex' and 'myney' are BOTH set to zero, nominal bunches\n
!!     centred in the aperture centre are generated. (SR, 08-05-2005)
!!     
!!     YIL EDIT 2010: particle 0 is always on orbit...
!! 
!! @author Javier Barranco <jbarranc@cern.ch>
!! @param mynp
!! @param myalphax
!! @param myalphay
!! @param mybetax
!! @param mybetay
!! @param myemitx0
!! @param myemity0
!! @param myenom
!! @param mynex
!! @param mdex
!! @param myney
!! @param mdey
!! @param myx
!! @param myxp
!! @param myy
!! @param myyp
!! @param myp
!! @param mys
!! @param enerror
!! @param bunchlength
!!
!! @date Last modified: 06. August 2009
!! @see ran_gauss
!! 
!<
      subroutine makedis_ga( mynp, myalphax, myalphay, mybetax,
     & mybetay, myemitx0, myemity0, myenom, mynex, mdex, myney, mdey,
     &     myx, myxp, myy, myyp, myp, mys,
     &     enerror, bunchlength )
!
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
!
+ca collpara
+ca dbmkdist

! !YIL debug july 2010
+ca parpro
+ca commont1

      double precision pi
!YIL march2010 edit: was missing enerror, bunchlength etc... 
! no common block for these parameters?
!
      double precision ran_gauss, gauss_rand
      double precision iix, iiy, phix, phiy
      double precision enerror, bunchlength
      double precision en_error, bunch_length
!
      double precision long_cut
      double precision a_st, b_st
      integer startpar
!
      save

!-----------------------------------------------------------------------
!++  Generate particle distribution
!
!
!++  Generate random distribution, assuming optical parameters at IP1
!
!++  Calculate the gammas
      pi=4d0*atan(1d0)
!
      mygammax = (1d0+myalphax**2)/mybetax
      mygammay = (1d0+myalphay**2)/mybetay
      en_error = enerror
      bunch_length = bunchlength

      write (lout,*) "Generation of bunch with dp/p and length:"
      write (lout,*) "  RMS bunch length  = ", bunch_length
      write (lout,*) "  RMS energy spread = ", en_error
! JBG August 2007
      write (lout,*)
      write (lout,*) "   ***STEP 1 for Gaussian Beam***"
      write (lout,*)
      write (lout,*) "   Beam generated with 5 sigma cut"
      write (lout,*)
      write (lout,*) "  Parameters used for Distribution Generation"
      write (lout,*) "  BetaX =", mybetax    
      write (lout,*) "  BetaY =", mybetay
      write (lout,*) "  EmittanceX =", myemitx0
      write (lout,*) "  EmittanceY =", myemity0
      write (lout,*)
      
      startpar=1
+if beamgas
      ! YIL July 2010 first particle on orbit
      !  initial xangle (if any) is not
      !  yet applied at this point...
      !  so we can set all to 0.
      startpar=2
      myx(1)=0.0d0                                                       !hr13
      myy(1)=0.0d0                                                       !hr13
      myxp(1)=0.0d0                                                      !hr13
      myyp(1)=0.0d0                                                      !hr13
      myp(1) = myenom
      mys(1) = 0d0
      !YIL end edit July 2010
+ei
      do j=startpar, mynp
! JBG July 2007    
! Option added for septum studies
!
            myemitx=myemitx0
            xsigmax = sqrt(mybetax*myemitx)
            myx(j)  = xsigmax * ran_gauss(mynex)
            myxp(j) = ran_gauss(mynex)*sqrt(myemitx/mybetax)-((myalphax*&!hr13
     &myx(j))/mybetax)                                                   !hr13
!    
!            if (rndm4().gt.0.5) then
!              myxp(j)  = sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-       
!     &              myalphax*myx(j)/mybetax
!              write(*,*)'Xp pos: ',myxp(j)
!            else
!              myxp(j)  = -1*sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-    
!     &              myalphax*myx(j)/mybetax
!              write(*,*)'Xp neg: ',myxp(j)
!            endif
!
           myemity=myemity0
           ysigmay = sqrt(mybetay*myemity)
!        write(*,*)'Sigma Y: ',ysigmay
            myy(j)   = ysigmay * ran_gauss(myney)
      myyp(j) = ran_gauss(myney)*sqrt(myemity/mybetay)-                 &!hr13
     &((myalphay*myy(j))/mybetay)                                        !hr13
    
!            myy(j)   = ysigmay * sin(2d0*pi*rndm4())
!            if (rndm4().gt.0.5) then
!              myyp(j)  = sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-        &
!     &              myalphay*myy(j)/mybetay
!            else
!              myyp(j)  = -1*sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-     &
!     &              myalphay*myy(j)/mybetay
!            endif
!
      end do
! SR, 11-08-2005 For longitudinal phase-space, add a cut at 2 sigma
!
!++   1st: generate mynpnumbers within the chosen cut
!
      long_cut = 2
      j = startpar
      do while (j.le.mynp)
         a_st = ran_gauss(5d0)
         b_st = ran_gauss(5d0)
         do while ((a_st*a_st+b_st*b_st).gt.long_cut*long_cut)
            a_st = ran_gauss(5d0)
            b_st = ran_gauss(5d0)
         enddo
         mys(j) = a_st
         myp(j) = b_st
         j = j + 1
      enddo
!++   2nd: give the correct values
      do j=startpar,mynp
         myp(j) = myenom * (1d0 + myp(j) * en_error)
         mys(j) = bunch_length * mys(j)
      enddo
!
      return
      end subroutine
! end of subroutine makedis_ga




      function rndm4()
      implicit none
+if crlibm
+ca crlibco
+ei
      integer len, in
      real rndm4, a
      save IN,a
      parameter ( len =  30000 )
      dimension a(len)
      data in/1/
!
      if ( in.eq.1 ) then
         call ranlux(a,len)
         rndm4=a(1)
         in=2
!        write(6,'('' LEN: '',i5)')LEN
      else
         rndm4=a(in)
         in=in+1
         if(in.eq.len+1)in=1
      endif
      return
      end
!
!
!ccccccccccccccccccccccccccccccccccccccc
!-TW-01/2007
! function rndm5(irnd) , irnd = 1 will reset 
! inn counter => enables reproducible set of 
! random unmbers
!cccccccccccccccccccccccccccccccccc
!
      function rndm5(irnd)
      implicit none
      integer len, inn, irnd
      real rndm5, a
      save
      parameter ( len =  30000 )
      dimension a(len)
      data inn/1/
!
! reset inn to 1 enable reproducible random numbers
      if ( irnd .eq. 1) inn = 1
      if ( inn.eq.1 ) then
         call ranlux(a,len)
         rndm5=a(1)
         inn=2
      else
         rndm5=a(inn)
         inn=inn+1
         if(inn.eq.len+1)inn=1
      endif
      return
      end
!
!ccccccccccccccccccccccccccccccccccccccc 
!
!
      double precision function myran_gauss(cut)
!*********************************************************************
!
! myran_gauss - will generate a normal distribution from a uniform
!     distribution between [0,1].
!     See "Communications of the ACM", V. 15 (1972), p. 873.
!
!     cut - double precision - cut for distribution in units of sigma
!     the cut must be greater than 0.5
!
!     changed rndm4 to rndm5(irnd) and defined flag as true 
! 
!*********************************************************************
      implicit none
      
+if crlibm
+ca crlibco
+ei
      logical flag
      real rndm5
      double precision x, u1, u2, twopi, r,cut
      save
      
      flag = .true. !Does this initialize only once, or is it executed every pass?
                    !See ran_gauss(cut)

+if crlibm
      twopi=8d0*atan_rn(1d0)                                             !hr09
+ei
+if .not.crlibm
      twopi=8d0*atan(1d0)
+ei
 1    if (flag) then
         r = dble(rndm5(0))
         r = max(r, 0.5d0**32)
         r = min(r, 1d0-0.5d0**32)
+if crlibm
         u1 = sqrt(-2d0*log_rn( r ))
+ei
+if .not.crlibm
         u1 = sqrt(-2d0*log( r ))
+ei
         u2 = dble(rndm5(0))
+if crlibm
         x = u1 * cos_rn(twopi*u2)
      else
         x = u1 * sin_rn(twopi*u2)
+ei
+if .not.crlibm
         x = u1 * cos(twopi*u2)
      else
         x = u1 * sin(twopi*u2)
+ei
      endif
      
      flag = .not. flag
      
!     cut the distribution if cut > 0.5
      if (cut .gt. 0.5d0 .and. abs(x) .gt. cut) goto 1
      
      myran_gauss = x
      return
      end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine ranlux(rvec,lenv)
!         Subtract-and-borrow random number generator proposed by
!         Marsaglia and Zaman, implemented by F. James with the name
!         RCARRY in 1991, and later improved by Martin Luescher
!         in 1993 to produce "Luxury Pseudorandom Numbers".
!     Fortran 77 coded by F. James, 1993
!
!   LUXURY LEVELS.
!   ------ ------      The available luxury levels are:
!
!  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
!           and Zaman, very long period, but fails many tests.
!  level 1  (p=48): considerable improvement in quality over level 0,
!           now passes the gap test, but still fails spectral test.
!  level 2  (p=97): passes all known tests, but theoretically still
!           defective.
!  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
!           correlations have very small chance of being observed.
!  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
!
!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!  Calling sequences for RANLUX:                                  ++
!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
!!!                   32-bit random floating point numbers between  ++
!!!                   zero (not included) and one (also not incl.). ++
!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
!!!               which is integer between zero and MAXLEV, or if   ++
!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
!!!               should be set to zero unless restarting at a break++
!!!               point given by output of RLUXAT (see RLUXAT).     ++
!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
!!!               which can be used to restart the RANLUX generator ++
!!!               at the current point by calling RLUXGO.  K1 and K2++
!!!               specify how many numbers were generated since the ++
!!!               initialization with LUX and INT.  The restarting  ++
!!!               skips over  K1+K2*E9   numbers, so it can be long.++
!!!   A more efficient but less convenient way of restarting is by: ++
!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
!!!                 32-bit integer seeds, to be used for restarting ++
!!!      ISVEC must be dimensioned 25 in the calling program        ++
!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
      integer lenv,isdext,iseeds,maxlev,ndskip,itwo24,next,j24,i24,     &
     &inseed,mkount,kount,in24,nskip,lxdflt,jsdflt,jseed,lp,i,k,icons,  &
     &inner,izip,izip2,ivec,isk,igiga,isd,k2,k1,inout,lout2,ins,lux,ilx,&
     &iouter
      real rvec,seeds,twop12,twom12,twom24,carry,uni
      dimension rvec(lenv)
      dimension seeds(24), iseeds(24), isdext(25)
      parameter (maxlev=4, lxdflt=3)
      dimension ndskip(0:maxlev)
      dimension next(24)
      parameter (twop12=4096., igiga=1000000000,jsdflt=314159265)
      parameter (itwo24=2**24, icons=2147483563)
      save notyet, i24, j24, carry, seeds, twom24, twom12, luxlev
      save nskip, ndskip, in24, next, kount, mkount, inseed
      integer luxlev
      logical notyet
      data notyet, luxlev, in24, kount, mkount /.true., lxdflt, 0,0,0/
      data i24,j24,carry/24,10,0./
!                               default
!  Luxury Level   0     1     2   *3*    4
      data ndskip/0,   24,   73,  199,  365 /
!Corresponds to p=24    48    97   223   389
!     time factor 1     2     3     6    10   on slow workstation
!                 1    1.5    2     3     5   on fast mainframe
!
!  NOTYET is .TRUE. if no initialization has been performed yet.
!              Default Initialization by Multiplicative Congruential
      if (notyet) then
         notyet = .false.
         jseed = jsdflt
         inseed = jseed
         write(lout,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',jseed
         luxlev = lxdflt
         nskip = ndskip(luxlev)
         lp = nskip + 24
         in24 = 0
         kount = 0
         mkount = 0
!         WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
!     &        LUXLEV,'      p =',LP
            twom24 = 1.
         do 25 i= 1, 24
            twom24 = twom24 * 0.5
         k = jseed/53668
         jseed = 40014*(jseed-k*53668) -k*12211
         if (jseed .lt. 0)  jseed = jseed+icons
         iseeds(i) = mod(jseed,itwo24)
   25    continue
         twom12 = twom24 * 4096.
         do 50 i= 1,24
         seeds(i) = real(iseeds(i))*twom24
         next(i) = i-1
   50    continue
         next(1) = 24
         i24 = 24
         j24 = 10
         carry = 0.
         if (seeds(24) .eq. 0.) carry = twom24
      endif
!
!          The Generator proper: "Subtract-with-borrow",
!          as proposed by Marsaglia and Zaman,
!          Florida State University, March, 1989
!
      do 100 ivec= 1, lenv
      uni = seeds(j24) - seeds(i24) - carry
      if (uni .lt. 0.)  then
         uni = uni + 1.
         carry = twom24
      else
         carry = 0.
      endif
      seeds(i24) = uni
      i24 = next(i24)
      j24 = next(j24)
      rvec(ivec) = uni
!  small numbers (with less than 12 "significant" bits) are "padded".
      if (uni .lt. twom12)  then
         rvec(ivec) = rvec(ivec) + twom24*seeds(j24)
!        and zero is forbidden in case someone takes a logarithm
         if (rvec(ivec) .eq. 0.)  rvec(ivec) = twom24*twom24
      endif
!        Skipping to luxury.  As proposed by Martin Luscher.
      in24 = in24 + 1
      if (in24 .eq. 24)  then
         in24 = 0
         kount = kount + nskip
         do 90 isk= 1, nskip
         uni = seeds(j24) - seeds(i24) - carry
         if (uni .lt. 0.)  then
            uni = uni + 1.
            carry = twom24
         else
            carry = 0.
         endif
         seeds(i24) = uni
         i24 = next(i24)
         j24 = next(j24)
   90    continue
      endif
  100 continue
      kount = kount + lenv
      if (kount .ge. igiga)  then
         mkount = mkount + 1
         kount = kount - igiga
      endif
      return
!
!           Entry to input and float integer seeds from previous run
      entry rluxin(isdext)
         notyet = .false.
         twom24 = 1.
         do 195 i= 1, 24
         next(i) = i-1
  195    twom24 = twom24 * 0.5
         next(1) = 24
         twom12 = twom24 * 4096.
      write(lout,*) ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
      write(lout,'(5X,5I12)') isdext
      do 200 i= 1, 24
      seeds(i) = real(isdext(i))*twom24
  200 continue
      carry = 0.
      if (isdext(25) .lt. 0)  carry = twom24
      isd = iabs(isdext(25))
      i24 = mod(isd,100)
      isd = isd/100
      j24 = mod(isd,100)
      isd = isd/100
      in24 = mod(isd,100)
      isd = isd/100
      luxlev = isd
        if (luxlev .le. maxlev) then
          nskip = ndskip(luxlev)
          write(lout,'(A,I2)')' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
     &luxlev
        else  if (luxlev .ge. 24) then
          nskip = luxlev - 24
          write(lout,'(A,I5)')' RANLUX P-VALUE SET BY RLUXIN TO:',luxlev
        else
          nskip = ndskip(maxlev)
          write(lout,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',luxlev
          luxlev = maxlev
        endif
      inseed = -1
      return
!
!                    Entry to ouput seeds as integers
      entry rluxut(isdext)
      do 300 i= 1, 24
         isdext(i) = int(seeds(i)*twop12*twop12)
  300 continue
      isdext(25) = i24 + 100*j24 + 10000*in24 + 1000000*luxlev
      if (carry .gt. 0.)  isdext(25) = -isdext(25)
      return
!
!     Entry to output the "convenient" restart point
!     Note: The first argument was originall called "lout";
!     however this conflicts with the variable name used for selecting output unit.
!     It was therefore renamed to "lout2".
      entry rluxat(lout2,inout,k1,k2)
      lout2 = luxlev
      inout = inseed
      k1 = kount
      k2 = mkount
      return
!
!                    Entry to initialize from one or three integers
      entry rluxgo(lux,ins,k1,k2)
         if (lux .lt. 0) then
            luxlev = lxdflt
         else if (lux .le. maxlev) then
            luxlev = lux
         else if (lux .lt. 24 .or. lux .gt. 2000) then
            luxlev = maxlev
            write(lout,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',lux
         else
            luxlev = lux
            do 310 ilx= 0, maxlev
              if (lux .eq. ndskip(ilx)+24)  luxlev = ilx
  310       continue
         endif
      if (luxlev .le. maxlev)  then
         nskip = ndskip(luxlev)
         write(lout,'(A,I2,A,I4)')                                      &
     &' RANLUX LUXURY LEVEL SET BY RLUXGO :',
     &luxlev,'     P=', nskip+24
      else
          nskip = luxlev - 24
          write(lout,'(A,I5)')                                          &
     &' RANLUX P-VALUE SET BY RLUXGO TO:',luxlev
      endif
      in24 = 0
      if (ins .lt. 0)  write(lout,*)
     &' Illegal initialization by RLUXGO, negative input seed'
      if (ins .gt. 0)  then
        jseed = ins
        write(lout,'(A,3I12)')                                          &
     &' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
     &jseed, k1,k2
      else
        jseed = jsdflt
        write(lout,*)' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
      endif
      inseed = jseed
      notyet = .false.
      twom24 = 1.
         do 325 i= 1, 24
           twom24 = twom24 * 0.5
         k = jseed/53668
         jseed = 40014*(jseed-k*53668) -k*12211
         if (jseed .lt. 0)  jseed = jseed+icons
         iseeds(i) = mod(jseed,itwo24)
  325    continue
      twom12 = twom24 * 4096.
         do 350 i= 1,24
         seeds(i) = real(iseeds(i))*twom24
         next(i) = i-1
  350    continue
      next(1) = 24
      i24 = 24
      j24 = 10
      carry = 0.
      if (seeds(24) .eq. 0.) carry = twom24
!        If restarting at a break point, skip K1 + IGIGA*K2
!        Note that this is the number of numbers delivered to
!        the user PLUS the number skipped (if luxury .GT. 0).
      kount = k1
      mkount = k2
      if (k1+k2 .ne. 0)  then
        do 500 iouter= 1, k2+1
          inner = igiga
          if (iouter .eq. k2+1)  inner = k1
          do 450 isk= 1, inner
            uni = seeds(j24) - seeds(i24) - carry
            if (uni .lt. 0.)  then
               uni = uni + 1.
               carry = twom24
            else
               carry = 0.
            endif
            seeds(i24) = uni
            i24 = next(i24)
            j24 = next(j24)
  450     continue
  500   continue
!         Get the right value of IN24 by direct calculation
        in24 = mod(kount, nskip+24)
        if (mkount .gt. 0)  then
           izip = mod(igiga, nskip+24)
           izip2 = mkount*izip + in24
           in24 = mod(izip2, nskip+24)
        endif
!       Now IN24 had better be between zero and 23 inclusive
        if (in24 .gt. 23) then
           write(lout,'(A/A,3I11,A,I5)')
     &'  Error in RESTARTING with RLUXGO:','  The values', ins,         &
     &k1, k2, ' cannot occur at luxury level', luxlev
           in24 = 0
        endif
      endif
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine funlxp (func,xfcum,x2low,x2high)
!         F. JAMES,   Sept, 1994
!
!         Prepares the user function FUNC for FUNLUX
!         Inspired by and mostly copied from FUNPRE and FUNRAN
!         except that
!    1. FUNLUX uses RANLUX underneath,
!    2. FUNLXP expands the first and last bins to cater for
!              functions with long tails on left and/or right,
!    3. FUNLXP calls FUNPCT to do the actual finding of percentiles.
!    4. both FUNLXP and FUNPCT use RADAPT for Gaussian integration.
!
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
      external func
      integer ifunc,ierr
      real x2high,x2low,xfcum,rteps,xhigh,xlow,xrange,uncert,x2,tftot1, &
     &x3,tftot2,func
+ca funint
      dimension xfcum(200)
      parameter (rteps=0.0002)
      save ifunc
      data ifunc/0/
      ifunc = ifunc + 1
!         FIND RANGE WHERE FUNCTION IS NON-ZERO.
      call funlz(func,x2low,x2high,xlow,xhigh)
      xrange = xhigh-xlow
      if(xrange .le. 0.)  then
        write(lout,'(A,2G15.5)') ' FUNLXP finds function range .LE.0',
     &xlow,xhigh
        go to 900
      endif
      call radapt(func,xlow,xhigh,1,rteps,0.,tftot ,uncert)
!      WRITE(6,1003) IFUNC,XLOW,XHIGH,TFTOT
 1003 format(' FUNLXP: integral of USER FUNCTION',                      &
     &i3,' from ',e12.5,' to ',e12.5,' is ',e14.6)
!
!      WRITE (6,'(A,A)') ' FUNLXP preparing ',
!     + 'first the whole range, then left tail, then right tail.'
      call funpct(func,ifunc,xlow,xhigh,xfcum,1,99,tftot,ierr)
      if (ierr .gt. 0)  go to 900
      x2 = xfcum(3)
      call radapt(func,xlow,x2,1,rteps,0.,tftot1 ,uncert)
      call funpct(func,ifunc,xlow,x2 ,xfcum,101,49,tftot1,ierr)
      if (ierr .gt. 0)  go to 900
      x3 = xfcum(98)
      call radapt(func,x3,xhigh,1,rteps,0.,tftot2 ,uncert)
      call funpct(func,ifunc,x3,xhigh,xfcum,151,49,tftot2,ierr)
      if (ierr .gt. 0)  go to 900
!      WRITE(6,1001) IFUNC,XLOW,XHIGH
 1001 format(' FUNLXP has prepared USER FUNCTION',i3,                   &
     &' between',g12.3,' and',g12.3,' for FUNLUX')
      return
  900 continue
      write(lout,*) ' Fatal error in FUNLXP. FUNLUX will not work.'
      end
!
      subroutine funpct(func,ifunc,xlow,xhigh,xfcum,nlo,nbins,tftot,    &
     &ierr)
!        Array XFCUM is filled from NLO to NLO+NBINS, which makes
!        the number of values NBINS+1, or the number of bins NBINS
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
      external func
      integer ierr,nbins,nlo,ifunc,nz,ibin,maxz,iz,nitmax,ihome
      real tftot,xhigh,xlow,func,xfcum,rteps,tpctil,tz,tzmax,x,f,tcum,  &
     &x1,f1,dxmax,fmin,fminz,xincr,tincr,xbest,dtbest,tpart,x2,precis,  &
     &refx,uncert,tpart2,dtpar2,dtabs,aberr
      dimension xfcum(*)
      parameter (rteps=0.005, nz=10, maxz=20, nitmax=6,precis=1e-6)
!      DOUBLE PRECISION TPCTIL, TZ, TCUM, XINCR, DTABS,
!     &  TINCR, TZMAX, XBEST, DTBEST, DTPAR2
!
      ierr = 0
      if (tftot .le. 0.) go to 900
      tpctil = tftot/real(nbins)                                         !hr09
      tz = tpctil/real(nz)
      tzmax = tz * 2.
      xfcum(nlo) = xlow
      xfcum(nlo+nbins) = xhigh
      x = xlow
      f = func(x)
      if (f .lt. 0.) go to 900
!         Loop over percentile bins
      do 600 ibin = nlo, nlo+nbins-2
      tcum = 0.
      x1 = x
      f1 = f
      dxmax = (xhigh -x) / nz
      fmin = tz/dxmax
      fminz = fmin
!         Loop over trapezoids within a supposed percentil
      do 500 iz= 1, maxz
      xincr = tz/max(f1,fmin,fminz)
  350 x = x1 + xincr
      f = func(x)
      if (f .lt. 0.) go to 900
      tincr = ((x-x1) * 0.5) * (f+f1)                                    !hr09
      if (tincr .lt. tzmax) go to 370
      xincr = xincr * 0.5
      go to 350
  370 continue
      tcum = tcum + tincr
      if (tcum .ge. tpctil*0.99) go to 520
      fminz = (tz*f)/ (tpctil-tcum)                                      !hr09
      f1 = f
      x1 = x
  500 continue
      write(lout,*) ' FUNLUX:  WARNING. FUNPCT fails trapezoid.'
!         END OF TRAPEZOID LOOP
!         Adjust interval using Gaussian integration with
!             Newton corrections since F is the derivative
  520 continue
      x1 = xfcum(ibin)
      xbest = x
      dtbest = tpctil
      tpart = tpctil
!         Allow for maximum NITMAX more iterations on RADAPT
      do 550 ihome= 1, nitmax
  535 xincr = (tpctil-tpart) / max(f,fmin)
      x = xbest + xincr
      x2 = x
        if (ihome .gt. 1 .and. x2 .eq. xbest) then
        write(lout,'(A,G12.3)')
     &' FUNLUX: WARNING from FUNPCT: insufficient precision at X=',x
        go to 580
        endif
      refx = abs(x)+precis
      call radapt(func,x1,x2,1,rteps,0.,tpart2,uncert)
      dtpar2 = tpart2-tpctil
      dtabs = abs(dtpar2)
      if(abs(xincr)/refx .lt. precis) goto 545
      if(dtabs .lt. dtbest) goto 545
      xincr = xincr * 0.5
      goto 535
  545 dtbest = dtabs
      xbest = x
      tpart = tpart2
      f = func(x)
      if(f .lt. 0.) goto 900
      if(dtabs .lt. rteps*tpctil) goto 580
  550 continue
      write(lout,'(A,I4)')
     &' FUNLUX: WARNING from FUNPCT: cannot converge, bin',ibin
!
  580 continue
      xincr = (tpctil-tpart) / max(f,fmin)
      x = xbest + xincr
      xfcum(ibin+1) = x
      f = func(x)
      if(f .lt. 0.) goto 900
  600 continue
!         END OF LOOP OVER BINS
      x1 = xfcum((nlo+nbins)-1)                                          !hr09
      x2 = xhigh
      call radapt(func,x1,x2,1,rteps,0.,tpart ,uncert)
      aberr = abs(tpart-tpctil)/tftot
!      WRITE(6,1001) IFUNC,XLOW,XHIGH
      if(aberr .gt. rteps)  write(lout,1002) aberr
      return
  900 write(lout,1000) x,f
      ierr = 1
      return
 1000 format(/' FUNLUX fatal error in FUNPCT: function negative:'/      &
     &,' at X=',e15.6,', F=',e15.6/)
! 1001 FORMAT(' FUNPCT has prepared USER FUNCTION',I3,
!     + ' between',G12.3,' and',G12.3,' for FUNLUX.')
 1002 format(' WARNING: Relative error in cumulative distribution',     &
     &' may be as big as',f10.7)
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine funlux(array,xran,len)
!         Generation of LEN random numbers in any given distribution,
!         by 4-point interpolation in the inverse cumulative distr.
!         which was previously generated by FUNLXP
      implicit none
+if crlibm
+ca crlibco
+ei
+ca funint
      integer len,ibuf,j,j1
      real array,xran,gap,gapinv,tleft,bright,gaps,gapins,x,p,a,b
      dimension array(200)
      dimension xran(len)
!  Bin width for main sequence, and its inverse
      parameter (gap= 1./99.,  gapinv=99.)
!  Top of left tail, bottom of right tail (each tail replaces 2 bins)
      parameter (tleft= 2./99.,bright=97./99.)
!  Bin width for minor sequences (tails), and its inverse
      parameter (gaps=tleft/49.,  gapins=1./gaps)
!
!   The array ARRAY is assumed to have the following structure:
!        ARRAY(1-100) contains the 99 bins of the inverse cumulative
!                     distribution of the entire function.
!        ARRAY(101-150) contains the 49-bin blowup of main bins
!                       1 and 2 (left tail of distribution)
!        ARRAY(151-200) contains the 49-bin blowup of main bins
!                       98 and 99 (right tail of distribution)
!
      call ranlux(xran,len)

      do 500 ibuf= 1, len
      x = xran(ibuf)
      j = int(  x    *gapinv) + 1
      if (j .lt. 3)  then
         j1 = int( x *gapins)
             j = j1 + 101
             j = max(j,102)
             j = min(j,148)
         p = (   x -gaps*real(j1-1)) * gapins                            !hr09
         a = (p+1.0) * array(j+2) - (p-2.0)*array(j-1)
         b = (p-1.0) * array(j) - p * array(j+1)
      xran(ibuf) = ((a*p)*(p-1.0))*0.16666667 + ((b*(p+1.))*(p-2.))*0.5  !hr09
      else if (j .gt. 97)  then
         j1 = int((x-bright)*gapins)
             j = j1 + 151
             j = max(j,152)
             j = min(j,198)
         p = ((x -bright) -gaps*(j1-1)) * gapins                         !hr09
         a = (p+1.0) * array(j+2) - (p-2.0)*array(j-1)
         b = (p-1.0) * array(j) - p * array(j+1)
      xran(ibuf) = ((a*p)*(p-1.0))*0.16666667 + ((b*(p+1.))*(p-2.))*0.5  !hr09
      else
!      J = MAX(J,2)
!      J = MIN(J,98)
         p = (   x -gap*real(j-1)) * gapinv                              !hr09
         a = (p+1.) * array(j+2) - (p-2.)*array(j-1)
         b = (p-1.) * array(j) - p * array(j+1)
      xran(ibuf) = ((a*p)*(p-1.))*0.16666667 + ((b*(p+1.))*(p-2.))*0.5   !hr09
      endif
  500 continue
      tftot = x
      return
      end
      subroutine funlz(func,x2low,x2high,xlow,xhigh)
!         FIND RANGE WHERE FUNC IS NON-ZERO.
!         WRITTEN 1980, F. JAMES
!         MODIFIED, NOV. 1985, TO FIX BUG AND GENERALIZE
!         TO FIND SIMPLY-CONNECTED NON-ZERO REGION (XLOW,XHIGH)
!         ANYWHERE WITHIN THE GIVEN REGION (X2LOW,H2HIGH).
!            WHERE 'ANYWHERE' MEANS EITHER AT THE LOWER OR UPPER
!            EDGE OF THE GIVEN REGION, OR, IF IN THE MIDDLE,
!            COVERING AT LEAST 1% OF THE GIVEN REGION.
!         OTHERWISE IT IS NOT GUARANTEED TO FIND THE NON-ZERO REGION.
!         IF FUNCTION EVERYWHERE ZERO, FUNLZ SETS XLOW=XHIGH=0.
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
      external func
      integer logn,nslice,i,k
      real xhigh,xlow,x2high,x2low,func,xmid,xh,xl,xnew
      xlow = x2low
      xhigh = x2high
!         FIND OUT IF FUNCTION IS ZERO AT ONE END OR BOTH
      xmid = xlow
      if (func(xlow) .gt. 0.) go to 120
      xmid = xhigh
      if (func(xhigh) .gt. 0.)  go to 50
!         FUNCTION IS ZERO AT BOTH ENDS,
!         LOOK FOR PLACE WHERE IT IS NON-ZERO.
      do 30 logn= 1, 7
      nslice = 2**logn
      do 20 i= 1, nslice, 2
      xmid = xlow + (real(i) * (xhigh-xlow)) / real(nslice)              !hr09
      if (func(xmid) .gt. 0.)  go to 50
   20 continue
   30 continue
!         FALLING THROUGH LOOP MEANS CANNOT FIND NON-ZERO VALUE
      write(lout,554)
      write(lout,555) xlow, xhigh
      xlow = 0.
      xhigh = 0.
      go to 220
!
   50 continue
!         DELETE 'LEADING' ZERO RANGE
      xh = xmid
      xl = xlow
      do 70 k= 1, 20
      xnew = 0.5*(xh+xl)
      if (func(xnew) .eq. 0.) go to 68
      xh = xnew
      go to 70
   68 xl = xnew
   70 continue
      xlow = xl
      write(lout,555) x2low,xlow
  120 continue
      if (func(xhigh) .gt. 0.) go to 220
!         DELETE 'TRAILING' RANGE OF ZEROES
      xl = xmid
      xh = xhigh
      do 170 k= 1, 20
      xnew = 0.5*(xh+xl)
      if (func(xnew) .eq. 0.) go to 168
      xl = xnew
      go to 170
  168 xh = xnew
  170 continue
      xhigh = xh
      write(lout,555) xhigh, x2high
!
  220 continue
      return
  554 format('0CANNOT FIND NON-ZERO FUNCTION VALUE')
  555 format(' FUNCTION IS ZERO FROM X=',e12.5,' TO ',e12.5)
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine radapt(f,a,b,nseg,reltol,abstol,res,err)

!     RES = Estimated Integral of F from A to B,
!     ERR = Estimated absolute error on RES.
!     NSEG  specifies how the adaptation is to be done:
!        =0   means use previous binning,
!        =1   means fully automatic, adapt until tolerance attained.
!        =n>1 means first split interval into n equal segments,
!             then adapt as necessary to attain tolerance.
!     The specified tolerances are:
!            relative: RELTOL ;  absolute: ABSTOL.
!        It stop s when one OR the other is satisfied, or number of
!        segments exceeds NDIM.  Either TOLA or TOLR (but not both!)
!        can be set to zero, in which case only the other is used.

      implicit none
+if crlibm
+ca crlibco
+ei
      external f
      integer nseg,ndim,nter,nsegd,i,iter,ibig
      real err,res,abstol,reltol,b,a,xlo,xhi,tval,ters,te,root,xhib,    &
     &bin,xlob,bige,hf,xnew,r1,f
      double precision tvals,terss

      parameter (ndim=100)
      parameter (r1 = 1., hf = r1/2.)

      dimension xlo(ndim),xhi(ndim),tval(ndim),ters(ndim)
      save xlo,xhi,tval,ters,nter
      data nter /0/

      if(nseg .le. 0)  then
       if(nter .eq. 0) then
        nsegd=1
        go to 2
       endif
       tvals=0d0
       terss=0d0
       do 1 i = 1,nter
       call rgs56p(f,xlo(i),xhi(i),tval(i),te)
       ters(i)=te**2
       tvals=tvals+dble(tval(i))                                         !hr09
       terss=terss+dble(ters(i))
    1  continue
       root= real(sqrt(2.d0*terss))                                      !hr09
       go to 9
      endif
      nsegd=min(nseg,ndim)
    2 xhib=a
      bin=(b-a)/real(nsegd)                                              !hr09
      do 3 i = 1,nsegd
      xlo(i)=xhib
      xlob=xlo(i)
      xhi(i)=xhib+bin
      if(i .eq. nsegd) xhi(i)=b
      xhib=xhi(i)
      call rgs56p(f,xlob,xhib,tval(i),te)
      ters(i)=te**2
    3 continue
      nter=nsegd
      do 4 iter = 1,ndim
      tvals=dble(tval(1))                                                !hr09
      terss=dble(ters(1))                                                !hr09
      do 5 i = 2,nter
      tvals=tvals+dble(tval(i))                                          !hr09
      terss=terss+dble(ters(i))                                          !hr09
    5 continue
      root= real(sqrt(2.d0*terss))                                       !hr09
      if(root .le. abstol .or. root .le. real(dble(reltol)*abs(tvals))) &!hr09
     &go to 9                                                            !hr09
      if(nter .eq. ndim) go to 9
      bige=ters(1)
      ibig=1
      do 6 i = 2,nter
      if(ters(i) .gt. bige) then
       bige=ters(i)
       ibig=i
      endif
    6 continue
      nter=nter+1
      xhi(nter)=xhi(ibig)
      xnew=hf*(xlo(ibig)+xhi(ibig))
      xhi(ibig)=xnew
      xlo(nter)=xnew
      call rgs56p(f,xlo(ibig),xhi(ibig),tval(ibig),te)
      ters(ibig)=te**2
      call rgs56p(f,xlo(nter),xhi(nter),tval(nter),te)
      ters(nter)=te**2
    4 continue
    9 res=real(tvals)                                                    !hr09
      err=root
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rgs56p(f,a,b,res,err)
      implicit none
+if crlibm
+ca crlibco
+ei
      integer i
      real err,res,b,a,f,w6,x6,w5,x5,rang,r1,hf
      double precision e5,e6

      parameter (r1 = 1., hf = r1/2.)
      dimension x5(5),w5(5),x6(6),w6(6)

      data (x5(i),w5(i),i=1,5)                                          &
     &/4.6910077030668004e-02, 1.1846344252809454e-01,                  &
     &2.3076534494715846e-01, 2.3931433524968324e-01,                   &
     &5.0000000000000000e-01, 2.8444444444444444e-01,                   &
     &7.6923465505284154e-01, 2.3931433524968324e-01,                   &
     &9.5308992296933200e-01, 1.1846344252809454e-01/

      data (x6(i),w6(i),i=1,6)                                          &
     &/3.3765242898423989e-02, 8.5662246189585178e-02,                  &
     &1.6939530676686775e-01, 1.8038078652406930e-01,                   &
     &3.8069040695840155e-01, 2.3395696728634552e-01,                   &
     &6.1930959304159845e-01, 2.3395696728634552e-01,                   &
     &8.3060469323313225e-01, 1.8038078652406930e-01,                   &
     &9.6623475710157601e-01, 8.5662246189585178e-02/

      rang=b-a
      e5=0d0
      e6=0d0
      do 1 i = 1,5
      e5=e5+dble(w5(i)*f(a+rang*x5(i)))                                  !hr09
      e6=e6+dble(w6(i)*f(a+rang*x6(i)))                                  !hr09
    1 continue
      e6=e6+dble(w6(6)*f(a+rang*x6(6)))
      res=real((dble(hf)*(e6+e5))*dble(rang))                            !hr09
      err=real(abs((e6-e5)*dble(rang)))                                  !hr09
      return
      end
!GRD
!
!*********************************************************************
!
! Define INTEGER function MCLOCK that can differ from system to system
!
!*********************************************************************
!
      integer function mclock_liar( )
!
      implicit none
+ca crcoall
+if crlibm
+ca crlibco
+ei
      save
!
      integer    mclock
      integer    count_rate, count_max
      logical    clock_ok
!
!        MCLOCK_LIAR = MCLOCK()
!
      clock_ok = .true.
!
      if (clock_ok) then
!
         call system_clock( mclock, count_rate, count_max )
         if ( count_max .eq. 0 ) then
            clock_ok = .false.
            write(lout,*)'INFO>  System Clock not present or not',
     &' Responding'
            write(lout,*)'INFO>  R.N.G. Reseed operation disabled.'
         endif
!
      endif
!
      mclock_liar = mclock
!
      return
      end
      double precision function ran_gauss(cut)
!*********************************************************************
!
! RAN_GAUSS - will generate a normal distribution from a uniform
!   distribution between [0,1].
!   See "Communications of the ACM", V. 15 (1972), p. 873.
!
! cut - double precision - cut for distribution in units of sigma
!                the cut must be greater than 0.5
!
!*********************************************************************
      implicit none
+if crlibm
+ca crlibco
+ei

      logical flag
      DATA flag/.TRUE./
      real rndm4
      double precision x, u1, u2, twopi, r,cut
      save
      
+if crlibm
            twopi=8d0*atan_rn(1d0) !Why not 2*pi, where pi is in block "common"?
+ei
+if .not.crlibm
            twopi=8d0*atan(1d0)
+ei
    1       if (flag) then
              r = dble(rndm4( ))
              r = max(r, 0.5d0**32)
              r = min(r, 1d0-0.5d0**32)
+if crlibm
              u1 = sqrt(-2d0*log_rn( r ))
+ei
+if .not.crlibm
              u1 = sqrt(-2d0*log( r ))
+ei
              u2 = dble(rndm4( ))
+if crlibm
              x = u1 * cos_rn(twopi*u2)
+ei
+if .not.crlibm
              x = u1 * cos(twopi*u2)
+ei
            else
+if crlibm
              x = u1 * sin_rn(twopi*u2)
+ei
+if .not.crlibm
              x = u1 * sin(twopi*u2)
+ei
            endif

          flag = .not. flag

!  cut the distribution if cut > 0.5
          if (cut .gt. 0.5d0 .and. abs(x) .gt. cut) goto 1

          ran_gauss = x
        return
      end

!>
!! readcollimator()
!! This routine is called once at the start of the simulation and
!! is used to read the collimator settings input file
!<
        subroutine readcollimator
!
        integer I,J,K
+ca crcoall
+ca parpro
+ca collpara
+ca database
+ca dbcommon
!

      logical lopen

      save
!
!--------------------------------------------------------------------
!++  Read collimator database
!
!      write(*,*) 'reading collimator database'
      inquire( unit=53, opened=lopen )
      if (lopen) then
         write(lout,*) "ERROR in subroutine readcollimator: "//
     &        "FORTRAN Unit 53 was already open!"
         call prror(-1)
      endif

      open(unit=53,file=coll_db, iostat=ios,
     &     status="OLD",action="read")
      if (ios.ne.0)then
         write(lout,*) "Error in subroutine readcollimator: "//
     &        "Could not open the file ",coll_db
         write(lout,*) "Got iostat=",ios
         call prror(-1)
      endif
!
!      write(*,*) 'inside collimator database'
!      I = 0
      read(53,*)
      read(53,*,iostat=ios) db_ncoll
      write(lout,*) 'number of collimators = ',db_ncoll
!     write(*,*) 'ios = ',ios
      if (ios.ne.0) then
        write(outlun,*) 'ERR>  Problem reading collimator DB ',ios
        call prror(-1)
      endif
      if (db_ncoll.gt.max_ncoll) then
         write(lout,*) 'ERR> db_ncoll > max_ncoll '
         call prror(-1)
      endif
!
      do j=1,db_ncoll
      read(53,*)
!GRD
!GRD ALLOW TO RECOGNIZE BOTH CAPITAL AND NORMAL LETTERS
!GRD
        read(53,*,iostat=ios) db_name1(j)
!        write(*,*) 'ios = ',ios
        if (ios.ne.0) then
          write(outlun,*) 'ERR>  Problem reading collimator DB ', j,ios
          call prror(-1)
        endif
!
        read(53,*,iostat=ios) db_name2(j)
!        write(*,*) 'ios = ',ios
        if (ios.ne.0) then
          write(outlun,*) 'ERR>  Problem reading collimator DB ', j,ios
          call prror(-1)
        endif
!
        read(53,*,iostat=ios) db_nsig(j)
!        write(*,*) 'ios = ',ios
        if (ios.ne.0) then
          write(outlun,*) 'ERR>  Problem reading collimator DB ', j,ios
          call prror(-1)
        endif
!GRD
        read(53,*,iostat=ios) db_material(j)
!        write(*,*) 'ios = ',ios
        if (ios.ne.0) then
          write(outlun,*) 'ERR>  Problem reading collimator DB ', j,ios
          call prror(-1)
        endif
        read(53,*,iostat=ios) db_length(j)
!        write(*,*) 'ios = ',ios
        if (ios.ne.0) then
          write(outlun,*) 'ERR>  Problem reading collimator DB ', j,ios
          call prror(-1)
        endif
        read(53,*,iostat=ios) db_rotation(j)
!        write(*,*) 'ios = ',ios
        if (ios.ne.0) then
          write(outlun,*) 'ERR>  Problem reading collimator DB ', j,ios
          call prror(-1)
        endif
        read(53,*,iostat=ios) db_offset(j)
!        write(*,*) 'ios = ',ios
        if (ios.ne.0) then
          write(outlun,*) 'ERR>  Problem reading collimator DB ', j,ios
          call prror(-1)
        endif
        read(53,*,iostat=ios) db_bx(j)
!        write(*,*) 'ios = ',ios
        if (ios.ne.0) then
          write(outlun,*) 'ERR>  Problem reading collimator DB ', j,ios
          call prror(-1)
        endif
        read(53,*,iostat=ios) db_by(j)
!        write(*,*) 'ios = ',ios
        if (ios.ne.0) then
          write(outlun,*) 'ERR>  Problem reading collimator DB ', j,ios
          call prror(-1)
        endif
      enddo
!
      close(53)
!
      end

