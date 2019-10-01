! ============================================================================ !
!
!  Crystal Collimation Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ============================================================================ !

module coll_crystal

  use floatPrecision

  implicit none

contains

end module coll_crystal

!===============================================================
! MDA: added icoll and MAX_NPART parameters
!
      SUBROUTINE collimate_cry(icoll, iturn, ie, name_coll,C_MATERIAL, C_LENGTH, &
     &                   C_ROTATION, &
     &                   C_APERTURE, C_OFFSET, C_TILT, &
     &                   X_IN, XP_IN, Y_IN,   &
     &                   YP_IN, P_IN, S_IN, NP, ENOM, LHIT, LHIT_TURN, &
     &                   PART_ABS, PART_ABS_TURN, IMPACT, INDIV, LINT, &
     &                   BX,BY,AX, &
     &                   AY,EMITX0,EMITY0, &
     &                   name,flagsec,dowrite_impact,MAX_NPART, Cry_tilt, Cry_length)!

!
!++  Based on routines by JBJ-R.assmann... Re-written for the crystal
! case by V.previtali in september 2008
!
!++  - Deleted all HBOOK stuff.
!++  - Deleted optics routine and all parser routines.
!++  - Replaced RANMAR call by RANLUX call
!++  - Included RANLUX code from CERNLIB into source
!++  - Changed dimensions from CGen(100,nmat) to CGen(200,nmat)
!++  - Replaced FUNPRE with FUNLXP
!++  - Replaced FUNRAN with FUNLUX
!++  - Included all CERNLIB code into source: RANLUX, FUNLXP, FUNLUX,
!++	                                      FUNPCT, FUNLZ, RADAPT,
!++                                           RGS56P
!++	 with additional entries:             RLUXIN, RLUXUT, RLUXAT,
!++                                           RLUXGO
!++
!++  - Changed program so that Nev is total number of particles
!++    (scattered and not-scattered)
!++  - Added debug comments
!++  - Put real dp/dx
!
!
!   !!!!!!!!   KNOWN ISSUES:  !!!!!!!
!      1- the pencil beam is not working in the 6track version (I did not change
!      the colltrack version yet)
!
      use coll_db
      use mod_ranlux
      use mod_funlux
      use coll_k2

      IMPLICIT NONE
!
      integer     MAX_NCOLL
      PARAMETER     (MAX_NCOLL = 99)
!
         integer         icoll      ! MDA: collimator ID to read parameters from new CollDB format
         integer         mat
         integer         Nev
         integer         j
         integer         nabs
         integer         NHIT
         integer         MAX_NPART
         integer         NP
!        PARAMETER         (MAX_NPART=1500000)
!
         double precision            p0
         double precision            xp_pencil0(MAX_NCOLL)
         double precision            yp_pencil0(MAX_NCOLL)
         double precision            x_pencil(MAX_NCOLL)
         double precision            y_pencil(MAX_NCOLL)
         double precision            zlm
         double precision            x,xp
         double precision            shift
         double precision            x_shift, xp_shift,s_shift !coordinates after shift/rotation
         double precision            x_rot, xp_rot,s_rot
         double precision            x_temp, xp_temp,s_temp !!all the _temp variables are used when you hit the cry from below
         double precision            tilt_int, x_int,xp_int,s_int       !all the _int variables are used when you hit the cry from below (int=interaction point)
         double precision            x00
         double precision            z
         double precision            z00
         double precision            zp
         double precision            p
         double precision            dpop
         double precision            s
         double precision            a_eq,b_eq,c_eq,Delta
         double precision            ENOM
!
         double precision            x_PRINT,xp_PRINT,y_PRINT,yp_PRINT
         double precision            AX,BX,AY,BY
         double precision            X_NORM,XP_NORM,Y_NORM,YP_NORM
         double precision            EMITX0
         double precision            EMITY0
!
         integer         LHIT(MAX_NPART)
         integer         LHIT_TURN(MAX_NPART)
         integer         PART_ABS(MAX_NPART)
         integer         PART_ABS_TURN(MAX_NPART)
!
!
         double precision            x_in(MAX_NPART)
         double precision            xp_in(MAX_NPART)
         double precision            y_in(MAX_NPART)
         double precision            yp_in(MAX_NPART)
         double precision            p_in(MAX_NPART)    !be careful: [Gev]
         double precision            s_in(MAX_NPART)
!        adding variables for the pencil beam. Variables in the absolute reference frame.
         double precision            x_in0(MAX_NPART)
         double precision            xp_in0(MAX_NPART)
         double precision            y_in0(MAX_NPART)
         double precision            yp_in0(MAX_NPART)
         double precision            p_in0(MAX_NPART)    !be careful: [Gev]
         double precision            s_in0(MAX_NPART)
         double precision            s_impact
         integer                     flagsec(MAX_NPART)
         logical                     dowrite_impact

!
         double precision            IMPACT(MAX_NPART)
         double precision            INDIV(MAX_NPART)
         double precision            LINT(MAX_NPART)
         integer                     name(MAX_NPART)
!
         double precision            x_out(MAX_NPART)
         double precision            xp_out(MAX_NPART)
         double precision            y_out(MAX_NPART)
         double precision            yp_out(MAX_NPART)
         double precision            p_out(MAX_NPART)
         double precision            s_out(MAX_NPART)
!
         double precision            fracab
         double precision            drift_length
         double precision            mirror
         double precision            tiltangle
         double precision            tiltangle2
!
         CHARACTER(LEN=4) C_MATERIAL     !Material
         double precision      C_LENGTH       !Length in m
         double precision      C_ROTATION     !Rotation angle vs vertical in radian
         double precision      C_APERTURE     !Aperture in m
         double precision      C_OFFSET       !Offset in m
         double precision      C_TILT(2)      !Tilt in radian
         double precision      C_TILT0(2)      !Tilt in radian
         double precision      PENCIL_DX(MAX_NCOLL)
         double precision      PENCIL_SPREAD(MAX_NCOLL)
         double precision      cry_bend

!
         LOGICAL           CHANGED_TILT1(MAX_NCOLL)
         LOGICAL           CHANGED_TILT2(MAX_NCOLL)
!
         COMMON /tilt/ CHANGED_TILT1, CHANGED_TILT2
!
!         REAL      rndm4
!         REAL      RAN_GAUSS
!
         common/materia/mat
         common/nommom/p0
!
        integer ie,iturn,nabs_total
!        COMMON  /INFO/ IE, ITURN, nabs_total
!
        integer   IPENCIL
!        integer   ICOLL
!        common  /icoll/  icoll

        common  /pencil/  xp_pencil0,yp_pencil0,pencil_dx,ipencil
!
        COMMON  /PENCIL2/ X_PENCIL, Y_PENCIL
!
!
        double precision        AMPLZ
!
      double precision XP_tangent
!
      double precision Rcurv,C_xmax,C_ymax           !crystal geometrical parameters  - be careful! are in [m]
      double precision Alayer                           !amorphous layer [mm]
      integer C_orient                           !crystal orientation [0-1]
      double precision Cry_tilt_part                    !crystal tilt [rad]
      double precision Cry_tilt                         !crystal tilt [rad]
      double precision Cry_tilt0                        !tilt of the crystal for having channeling (determined by divergence of the beam) [rad]
      double precision Cry_length                       !original length (from the db) [m]
      double precision miscut
!                                           !instead of this parameter, I use (mat-7)
!

      integer bool_proc(1000)
      integer bool_proc_old(1000)
      integer n_chan
      integer n_VR
      integer n_amorphous
      character(len=*) name_coll             ! MDA: length to be defined by input
      CHARACTER(len=50) PROC                     !string that contains the physical process
!                                            !=' ' if the particle does not pass by crystal, ='*' if there is interaction
      logical  bool_create
      common /miscut/ miscut
!      common /Par_Cry1/ Cry_length, Rcurv,C_xmax,C_ymax,Alayer,C_orient
      common /Par_Cry1/ Rcurv,C_xmax,C_ymax,Alayer,C_orient

!      common /Par_Cry2/ Cry_tilt,Cry_tilt0
      common /Process/ bool_proc,bool_create
      common /Process_old/ bool_proc_old
      common/Proc2/PROC
      logical write_c_out             !if set to 1, coordinates of the
      common /outputs/ write_c_out

!------daniele----------
! adding new variables for debugged studies on process experienced

      integer idx_proc                !daniele
      integer   samplenumber                !daniele
      character(len=4) smpl                !daniele
      character(len=80) pfile                !daniele
      common /samplenumber/ pfile,smpl,samplenumber                !daniele

!----------------------

! MDA: parameters assignment from new CollDB format

      Rcurv = cdb_cryBend(icoll)
!      write(*,*) "Rcurv", Rcurv
      Alayer = cdb_cryThick(icoll)
      C_xmax = cdb_cryXDim(icoll)
      C_ymax = cdb_cryYDim(icoll)
      C_orient = cdb_cryOrient(icoll)
!      Cry_tilt = cdb_cryTilt(icoll)
!      Cry_length = cdb_cLength(icoll)
      miscut = cdb_cryMisCut(icoll)

!      write(*,*) Rcurv, Alayer, C_xmax, C_ymax, C_orient, Cry_tilt, Cry_length, miscut
!      stop

!
  666 FORMAT(A,1x,e20.10,1x,A,1x,e20.10,A,1x,e20.10,1x,A,1x,e20.10)
!=======================================================================
!
!      write(*,*) 'enter collimate_cry routine, rotation ',C_ROTATION

!      open(unit=9999,file='debug.dat')
      IF (C_MATERIAL(1:4).eq.'Si')THEN           ! MDA: changed the label for silicon from CRY-Si to Si to fit new format
           mat = 8
      ELSEIF (C_MATERIAL.eq.'CRY-W')THEN
           mat = 9
      ELSEIF (C_MATERIAL.eq.'CRY-C')THEN
           mat = 10
      ELSEIF (C_MATERIAL.eq.'CRY-Ge')THEN
           mat = 11
      ELSE
           WRITE(*,*) 'ERR>', C_MATERIAL, ' Material not found. STOP'
           STOP
      ENDIF
!

!          write(*,*) "debug - length  bent" , C_LENGTH
!          write(*,*) "debug - length  unbent" , CRY_LENGTH
!      write(*,*) "Cry_length", Cry_length
      Cry_bend =  Cry_length/Rcurv !cry_length longitudinal estension of the straight crystal
!          write(*,*) "debug - bend angle" , CRY_BEND
!          write(*,*) "debug - C_xmax" , C_xmax
!          write(*,*) "debug - miscut angle" , miscut
      if (c_length .gt. 0.) then
       NEV = NP
       P0  = ENOM
       C_TILT0(1) = C_TILT(1)
       C_TILT0(2) = C_TILT(2)
       tiltangle=C_TILT0(1)
!
!++  Initialize scattering processes
!
       call k2coll_scatin(p0)

! EVENT LOOP,  initial distribution is here a flat distribution with
! xmin=x-, xmax=x+, etc. from the input file
!
      nhit    = 0
      fracab  = 0.
      n_chan  = 0          !valentina :initialize to zero the counters for crystal effects
      n_VR    = 0          !
      n_amorphous = 0      !
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!        WRITE(*,*) ITURN,ICOLL
      do j = 1, nev
!
        impact(j) = -1.
        lint(j)   = -1.
        indiv(j)  = -1.

        idx_proc = name(j)-100*samplenumber

        if (ITURN .eq. 1) then                                                          !daniele
                bool_proc_old(idx_proc)=-1                                                     !daniele
        else                                                                            !daniele
               if (bool_proc(idx_proc).ne.-1)     &
     &                       bool_proc_old(idx_proc)=bool_proc(idx_proc)                  !daniele
        endif                                                                           !daniele
        PROC='out' !the default process is 'out'                                        !daniele
        bool_proc(idx_proc)=-1                                                                 !daniele


        nabs = 0
!
!        xp_in(j) = 1.4291566891454061E-005_fPrec
!        x_in(j) = 2.3678668975404154E-003_fPrec
!        y_in(j) = 9.7959100975990727E-005_fPrec
!        yp_in(j) = -3.0041857465904482E-006_fPrec
!        p_in(j) = 6499.9999993204920_fPrec
!        s_in(j) = 0.0000000000000000_fPrec

        write(*,*) "Coordinates after entering"
        write(*,*) x_in(j), xp_in(j), y_in(j), yp_in(j)
        write(*,*) p_in(j), s_in(j)

        s   = 0
        x   = x_in(j)
        xp  = xp_in(j)
        z   = y_in(j)
        zp  = yp_in(j)
        p   = p_in(j)
!
        x_temp=0.
        x_int=0.
        x_rot=0
        x_shift=0
        s_temp=0.
        s_int=0.
        s_rot=0
        s_shift=0
        xp_int=0.
        xp_temp=0.
        xp_rot=0
        xp_shift=0
        shift=0.
        tilt_int=0.
        dpop = (p -p0)/p0
!        Cry_tilt_part=Cry_tilt                         !crystal tilt [rad]
!             write(*,*)'collimator',mat,'particle',j,'x',x,
!     1       'in sigma', x/sqrt(bx)/sqrt(EMITX0),'sig_sk',
!     2       C_APERTURE/2
!




!---------------DANIELE-------------------------
! corrected position for variable association for FirstImpact.dat
! also for a vertical crystal case.
! Only x_in0 (i.e. b) have to be assigned after the change of reference frame
!-----------------------------------------------

        s_in0(j)   = s_in(j)                         !daniele
!        x_in0(j)   = x                         !daniele
        xp_in0(j)  = xp                         !daniele
        y_in0(j)   = z                         !daniele
        yp_in0(j)  = zp                         !daniele
        p_in0(j)   = p                         !daniele

!/*----------------DANIELE------------------*/




!++  transform particle coordinates to get into collimator coordinate
!++  system
!
!++  first check whether particle was lost before
!
        if (x.lt.99.0*1d-3 .and. z.lt.99.0*1d-3) then
! /
!++  first do rotation into collimator frame
!

!          write(*,*) "debug - rotazione" ,c_rotation
          x  = x_in(j)*cos(c_rotation) +sin(c_rotation)*y_in(j)
          z  = y_in(j)*cos(c_rotation) -sin(c_rotation)*x_in(j)
          xp = xp_in(j)*cos(c_rotation)+sin(c_rotation)*yp_in(j)
          zp = yp_in(j)*cos(c_rotation)-sin(c_rotation)*xp_in(j)

!          write(*,*) x, xp, z, zp

!
!++  for one-sided collimators consider only positive x. for negative
!++  x jump to the next particle

!          if ((name_coll(1:10) .eq. "CRY.SPSEXP")) then
!          mirror=-1
!!             write(*,*) "valentina crsitallo SPS riconosciuto"
!             if (x .gt. 0) then
!                goto 777
!             endif
!          else
!             mirror=1
!             if  (x.lt.0 ) then
!                 goto 777
!             endif
!          endif
!          x  = mirror * x
!          xp = mirror * xp



!++Shift with opening and offset
!
          X  = X - C_APERTURE/2 - C_OFFSET
!
!++  Include collimator tilt
!
!          write(*,*) "debug - collimator tilt" ,tiltangle
          IF (tiltangle.GT.0.) THEN
            XP = XP - tiltangle
          ELSEIF (tiltangle.LT.0.) THEN
            X  = X + SIN(tiltangle) * C_LENGTH
            XP = XP - tiltangle
          ENDIF
!
!++  For selected collimator, first turn reset particle distribution
!++  to simple pencil beam
!
          IF (ipencil .eq. 0) bool_create=.true.
!
!! pencil beam stuffffffffffffffffffffffffffffffffff  (~80 lines)
!c ----------------------------pencil beam generatio at the crystal--------------------------------------
!          IF ( (ICOLL.EQ.IPENCIL
!     1           .AND. ITURN.EQ.1)
!     2           ) THEN
!            X    = pencil_dx(ICOLL)
!            XP = cry_tilt0 !valentina (if the beam is generated @ crystal, I want it to have the natural beam divergence)
!            bool_create=.true.
!c
!            AMPLZ = RAN_GAUSS(3.)!*pencil_spread(ICOLL)
!            Z     = AMPLZ
!            ZP   = 0.
!            dpop = 0.
!c
!            if (write_c_out) then
!c++             I want to write the original coordinates in a file... I have to transform back all coordinates....Include collimator tilt
!c
!              IF (tiltangle.GT.0.) THEN
!                x_print  = X  + tiltangle*C_LENGTH
!                XP_print = XP + tiltangle
!              ELSEIF (tiltangle.LT.0.) THEN
!                x_print  = X + tiltangle*C_LENGTH
!                XP_print = XP + tiltangle
!                x_print  = X - SIN(tiltangle) * C_LENGTH
!              ELSE
!                x_print = x
!                xp_print = xp
!              ENDIF
!
!c++  Transform back to particle coordinates with opening and offset
!c
!              x_print = x_print + C_APERTURE/2 + MIRROR*C_OFFSET
!c++  Now mirror at the horizontal axis for negative X offset
!c
!              x_print    = MIRROR * x_print
!              XP_print   = MIRROR * xp_print
!
!c++  Last do rotation into collimator frame
!c
!              Y_print  = Z  *COS(-1.*C_ROTATION) -
!     1                     X_print  *SIN(-1.*C_ROTATION)
!              x_print  = x_print  *COS(-1.*C_ROTATION) +
!     1                     Z  *SIN(-1.*C_ROTATION)
!              YP_print = ZP *COS(-1.*C_ROTATION) -
!     1                     XP_print *SIN(-1.*C_ROTATION)
!              XP_print = XP_print *COS(-1.*C_ROTATION) +
!     1                     ZP *SIN(-1.*C_ROTATION)
!c              WRITE(*,'(i4,2x,i4,2x,a,2x,5(f15.8,2x))')
!c     2          ITURN,ICOLL,C_MATERIAL,x_print,XP_print,Y_print,YP_print
!c     3          ,P
!C
!c++ Then drift forward of half lenght
!c
!              x_print= x_print+xp_print*C_length/2        !valentina drift removed
!              y_print= y_print+yp_print*C_length/2
!              WRITE(881,'(i4,2x,i4,2x,a,2x,5(f15.8,2x))')
!     2        ITURN,ICOLL,C_MATERIAL,x_print,XP_print,Y_print,YP_print
!     3        ,P
!
!              X_NORM=x_print/ SQRT(bx)/sqrt(EMITX0)
!              XP_NORM=(x_print*ax+xp_print*bx)/SQRT(bx)/sqrt(EMITX0)
!              Y_NORM=y_print/ SQRT(by)/sqrt(EMITY0)
!              YP_NORM=(y_print*AY+yp_print*by)/SQRT(by)/sqrt(EMITY0)
!c
!              WRITE(883,'(i4,2x,i4,2x,a,2x,7(f15.8,2x))')
!     2          ITURN,ICOLL,C_MATERIAL,X_NORM,XP_NORM,Y_NORM,YP_NORM,
!     3          SQRT(X_NORM**2+XP_NORM**2),
!     4          SQRT(Y_NORM**2+YP_NORM**2),P
!c
!c
!            endif
!          ENDIF
!c.--------------------end of pencil beam stuff-------------------------------


!-valentina for the first impact file
!        s_in0(j)   = s_in(j)                         !daniele SEE COMMENTS ABOVE
        x_in0(j)   = x                         !daniele
!        xp_in0(j)  = xp                         !daniele
!        y_in0(j)   = z                         !daniele
!        yp_in0(j)  = zp                         !daniele
!        p_in0(j)   = p                         !daniele


!        write(*,*) "debug - coll RFS" , s_in(j),x,xp,z,zp




!-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o


!
! crystal!!! transform in the crystal rerference system
!    1st transformation: shift of the center of my reference frame
!          write(*,*) "debug - cry tilt" , Cry_tilt
!          write(*,*) "debug - cry tilt 0" , Cry_tilt0

         if (Cry_tilt .lt. 0) then
           S_shift=S
           write(*,*) j,'- s=',s
           shift=Rcurv*(1-cos(Cry_tilt))
           if (Cry_tilt .lt. (-Cry_bend) ) then
                shift= ( Rcurv * &
     &          ( cos ( - Cry_tilt) &
     &          - cos( Cry_bend - Cry_tilt) ) )
           endif
           X_shift=X-shift
         else
           S_shift=S
           X_shift=X
         endif
          write(*,*) "debug - S shift" ,  S_shift
          write(*,*) "debug - X shift" ,  X_shift
!
!    2nd transformation: rotation
         S_rot =X_shift*sin(Cry_tilt)+S_shift*cos(Cry_tilt)
         X_rot = X_shift*cos(Cry_tilt)-S_shift*sin(Cry_tilt)
         XP_rot= XP - Cry_tilt
          write(*,*) "debug - cry tilt", Cry_tilt
          write(*,*) "debug - S rot" ,  S_rot
          write(*,*) "debug - X rot" ,  X_rot
          write(*,*) "debug - XP rot" ,  XP_rot
!    3rd transformation: drift to the new coordinate s=0
         XP=XP_rot
         X= X_rot - XP_rot*S_rot
         Z= Z - ZP*S_rot
         S=0
          write(*,*) "debug - S cryRF" ,  S_rot
          write(*,*) "debug - X cryRF" ,  X_rot
          write(*,*) "debug - XP cryRF" ,  XP_rot
!
!  NOW CHECK IF THE PARTICLE HIT the crystal
!

!          write(*,*) "debug - checking if crystal hit"

         if (x .ge. 0 .and. x.lt.C_xmax) then !Daniele: check that part. hit cry
           s_impact=s_in0(j) !(for the first impact)
!           write(*,*) "HIT", x, C_xmax
          write(*,*) "HIT"
!           stop
!           write(*,*)'hit the cry entrance face', x, C_xmax
!           write(*,*)'impact at s,x = ', s_impact,x_in0(j)
!           write(*,*)'with angle xp = ',xp
!           write(*,*)'s before', s
           CALL CRYST(mat-7,X,XP,Z,ZP,p,cry_length)
           s=Rcurv*sin(cry_bend)
           zlm=Rcurv*sin(cry_bend)
!           write(*,*) 'process:',PROC
!           write(*,*)'s after', s
!           write(*,*)'hit the crystal'
           if (PROC(1:3).ne.'out')then
             NHIT = NHIT + 1
             LHIT(j) = ie
             LHIT_TURN(j) = ITURN
             IMPACT(j) = X_in0(j)
             INDIV(j) = XP_in0(j)
           endif
         else
           XP_tangent=sqrt((-2*X*Rcurv+X**2)/(Rcurv**2))
!           write(*,*) "NOT hit the cry entrance face", x, C_xmax
!           stop
!           write(*,*)j,'-','tangent',xp_tangent,'angle',xp
!           write(*,*)'s tan',Rcurv*sin(XP_tangent)
!           write(*,*) 's tot', c_length,Rcurv*sin(cry_bend)
           if ( XP .ge. XP_tangent  ) then

! if it hits the crystal, calculate in which point and apply the
! transformation and drift to that point
             a_eq=(1.+xp**2)
             b_eq=2.*xp*(x-Rcurv)
             c_eq=-2.*x*Rcurv+x**2
             Delta=b_eq**2-4.*(a_eq*c_eq)
             S_int=(-b_eq-sqrt(Delta))/(2.*a_eq)
!             write(*,*)'s int',S_int
             if (S_int .lt. Rcurv*sin(cry_bend)) then
!  transform to a new ref system:shift and rotate
               X_int=XP*S_int+X
               XP_int=XP
               Z=Z+ZP*S_int
               X=0
               S=0
!               tilt_int=2*X_int/S_int
               tilt_int=S_int/Rcurv
               XP=XP-tilt_int
!               write(*,*)'hit the cry from below!!!'
!               write(*,*)'tilt int',tilt_int,'over',cry_bend
!               write(*,*)'tilt bending',Cry_length/Rcurv,
!     &         'total tilt', cry_tilt-cry_tilt0,
!     &         'int. tilt', tilt_int
!               s_impact=Rcurv*(sin(Cry_length/Rcurv)
!     &           -sin(Cry_length/Rcurv-tilt_int))!(for the first impact)
!               x_in0(j)=Rcurv*(1-cos(Cry_length/Rcurv-tilt_int))
!               write(*,*)'impact at s,x = ', s_impact,x_in0(j)
!               write(*,*)'with angle xp = ',xp
!               write(*,*) "debug - S minicry" ,  S
!               write(*,*) "debug - X minicry" ,  X
!               write(*,*) "debug - XP minicry" ,  XP
! call cry routine
              CALL CRYST(mat-7,X,XP,Z,ZP,p,(cry_length-(tilt_int*Rcurv)))
              s=Rcurv*sin(cry_bend-tilt_int)
              zlm=Rcurv*sin(cry_bend-tilt_int)
!              write(*,*) 'process:',PROC
!              write(*,*) "debug - S minicry 2" ,  S
!              write(*,*) "debug - X minicry 2" ,  X
!              write(*,*) "debug - XP minicry 2" ,  XP
              if (PROC(1:3).ne.'out')then
                 X_rot=X_int
                 S_rot=S_int
                 XP_rot=XP_int
                 S_shift=S_rot*cos(-Cry_tilt)+X_rot*sin(-Cry_tilt)
                 X_shift=-S_rot*sin(-Cry_tilt)+X_rot*cos(-Cry_tilt)
                 XP_shift=XP_rot + Cry_tilt
                 if (Cry_tilt .lt. 0) then
                   S_impact=S_shift
                   X_in0(j)=X_shift+shift
                   XP_in0(j)=XP_shift
                 else
                   X_in0(j)=X_shift
                   S_impact=S_shift
                   XP_in0(j)=XP_shift
                 endif
                 NHIT = NHIT + 1
                 LHIT(j) = ie
                 LHIT_TURN(j) = ITURN
                 IMPACT(j) = X_in0(j)
                 INDIV(j) = XP_in0(j)
               endif
!           write(*,*)'s after', s
! un-rotate
               X_temp=X
               S_temp=S
               XP_temp=XP
               S=S_temp*cos(-tilt_int)+X_temp*sin(-tilt_int)
               X=-S_temp*sin(-tilt_int)+X_temp*cos(-tilt_int)
               XP=XP_temp + tilt_int
!     2nd: shift back the 2 axis
               X=X+X_int
               S=S+S_int
!               write(*,*)'s after', s
             else
!               write(*,*)'treat the drift'
               S=Rcurv*sin(cry_length/Rcurv)
               X=X+S*XP
               Z=Z+S*ZP
             endif
           else
!              write(*,*) 'just the drift'
             S=Rcurv*sin(cry_length/Rcurv)
             X=X+S*XP
             Z=Z+S*ZP
           endif
         endif
!               WRITE(*,*)'X1_cry',X,'Z1_Cry',Z,'XP1_Cry',XP,'ZP1_Cry',ZP
!     1         ,'s',s, Cry_tilt
!




! trasform back from the crystal to the collimator reference system
!    1st: un-rotate the coordinates
               X_rot=X
               S_rot=S
               XP_rot=XP
!               write(*,*) "debug - S cryRF 2" ,  S_rot
!               write(*,*) "debug - X cryRF 2" ,  X_rot
!               write(*,*) "debug - XP cryRF 2" ,  XP_rot
               S_shift=S_rot*cos(-Cry_tilt)+X_rot*sin(-Cry_tilt)
               X_shift=-S_rot*sin(-Cry_tilt)+X_rot*cos(-Cry_tilt)
               XP_shift=XP_rot + Cry_tilt
!     2nd: shift back the reference frame
               if (Cry_tilt .lt. 0) then
                 S=S_shift
                 X=X_shift+shift
                 XP=XP_shift
               else
                 X=X_shift
                 S=S_shift
                 XP=XP_shift
               endif
!     3rd: shift to new S=Length position
               X=XP*(c_length-S)+X
               Z=ZP*(c_length-S)+Z
               S=c_length


!               write(*,*) "debug - S Coll RF 2" ,  S_rot
!               write(*,*) "debug - X Coll RF 2" ,  X_rot
!               write(*,*) "debug - XP Coll RF 2" ,  XP_rot
!
!          WRITE(*,*)'X1_coll',X,'Z1_coll',Z,'XP1_coll',XP,'ZP1_coll',ZP
!     1         ,'s' ,s

               NABS=0

!               if (PROC(1:2).eq.'AM')then
!                 bool_proc(j)=1
!                 n_amorphous = n_amorphous + 1
!               elseif (PROC(1:2).eq.'VR') then
!                 bool_proc(j)=2
!                 n_VR = n_VR + 1
!               elseif (PROC(1:2).eq.'CH')then
!                 bool_proc(j) = 3
!                 n_chan = n_Chan + 1
!               elseif (PROC(1:2).eq.'VC') then
!                 bool_proc(j)=3
!                 n_chan = n_Chan + 1
!               elseif (PROC(1:3).eq.'out')then
!                 bool_proc(j)=-1
!               elseif (PROC(1:8).eq.'absorbed') then
!                 bool_proc(j)=5
!                 NABS=1
!               elseif (PROC(1:2).eq.'DC')then
!                 bool_proc(j)=6
!               elseif (PROC(1:3).eq.'mcs')then        !daniele
!                 bool_proc(j)=7
!               elseif (PROC(1:4).eq.'diff')then       !daniele
!                 bool_proc(j)=8
!               else
!                write(*,*)'???????????????????',PROC(1:2)
!                stop
!               endif



!  ----------------------DANIELE----------
! debugged assignation of the process experienced in the crystal
! ----------------

!               if (PROC(1:3).ne.'out')then
!                 write(*,*) PROC
!              endif

               if (PROC(1:2).eq.'AM')then
                 bool_proc(idx_proc)=1
                 n_amorphous = n_amorphous + 1
               elseif (PROC(1:2).eq.'VR') then
                 bool_proc(idx_proc)=2
                 n_VR = n_VR + 1
               elseif (PROC(1:2).eq.'CH')then
                 bool_proc(idx_proc) = 3
                 n_chan = n_Chan + 1
               elseif (PROC(1:2).eq.'VC') then
                 bool_proc(idx_proc)=4
                 !n_chan = n_Chan + 1
               elseif (PROC(1:3).eq.'out')then
                 bool_proc(idx_proc)=-1
               elseif (PROC(1:8).eq.'absorbed') then
                 bool_proc(idx_proc)=5
                 NABS=1
               elseif (PROC(1:2).eq.'DC')then
                 bool_proc(idx_proc)=6
               elseif (PROC(1:3).eq.'pne')then
                 bool_proc(idx_proc)=7
               elseif (PROC(1:3).eq.'ppe')then
                 bool_proc(idx_proc)=8
               elseif (PROC(1:4).eq.'diff')then
                 bool_proc(idx_proc)=9
               elseif (PROC(1:4).eq.'ruth')then
                 bool_proc(idx_proc)=10
               elseif (PROC(1:11).eq.'ch_absorbed') then
                 bool_proc(idx_proc)=15
                 NABS=1
               elseif (PROC(1:6).eq.'ch_pne')then
                 bool_proc(idx_proc)=17
               elseif (PROC(1:6).eq.'ch_ppe')then
                 bool_proc(idx_proc)=18
               elseif (PROC(1:7).eq.'ch_diff')then
                 bool_proc(idx_proc)=19
               elseif (PROC(1:7).eq.'ch_ruth')then
                 bool_proc(idx_proc)=20
               elseif (PROC(1:4).eq.'TRVR')then
                 bool_proc(idx_proc)=100
               elseif (PROC(1:4).eq.'TRAM')then
                 bool_proc(idx_proc)=101
               else
                write(*,*)'???????????????????',PROC(1:2)
                stop

               endif

!               if (PROC(1:3).ne.'out')then
!                 write(*,*) bool_proc(idx_proc)
!               endif

! ---------------------END DANIELE--------------


!===========================
!++  Transform back to particle coordinates with opening and offset
!
       IF (PART_ABS(j).eq.0) THEN
!
!++  Include collimator tilt
!
         IF (tiltangle.GT.0.) THEN
           X  = X  + tiltangle*C_LENGTH
           XP = XP + tiltangle
         ELSEIF (tiltangle.LT.0.) THEN
           X  = X + tiltangle*C_LENGTH
           XP = XP + tiltangle
!
           X  = X - SIN(tiltangle) * C_LENGTH
         ENDIF
!
!++  Transform back to particle coordinates with opening and offset
!
         Z00 = Z
         X00 = X + MIRROR*C_OFFSET
         X = X + C_APERTURE/2 + MIRROR*C_OFFSET
!
!++  Now mirror at the horizontal axis for negative X offset
!
         X    = MIRROR * X
         XP   = MIRROR * XP


!
!++  Last do rotation into collimator frame
!
         X_IN(J)  = X  *COS(-1.*C_ROTATION) + Z  *SIN(-1.*C_ROTATION)
         Y_IN(J)  = Z  *COS(-1.*C_ROTATION) - X  *SIN(-1.*C_ROTATION)
         XP_IN(J) = XP *COS(-1.*C_ROTATION) + ZP *SIN(-1.*C_ROTATION)
         YP_IN(J) = ZP *COS(-1.*C_ROTATION) - XP *SIN(-1.*C_ROTATION)

!----- other pencil beam stuff-------
         IF ( ICOLL.EQ.IPENCIL) then
           X00  = MIRROR * X00
           X_IN(J)  = X00  *COS(-1.*C_ROTATION) + Z00  *SIN(-1.*C_ROTATION)
           Y_IN(J)  = Z00  *COS(-1.*C_ROTATION) - X00  *SIN(-1.*C_ROTATION)
!
           XP_IN(J) = XP_IN(J) + MIRROR*XP_PENCIL0(ICOLL)
           YP_IN(J) = YP_IN(J) + MIRROR*YP_PENCIL0(ICOLL)
           X_IN(J) = X_IN(J) + MIRROR*X_PENCIL(ICOLL)
           Y_IN(J) = Y_IN(J) + MIRROR*Y_PENCIL(ICOLL)

           IF (.NOT. CHANGED_TILT1(ICOLL) .AND. MIRROR.GT.0.) THEN
                   WRITE (*,*) 'NEVER!!!'
                   C_TILT(1) = XP_PENCIL0(ICOLL)*COS(C_ROTATION)+SIN(C_ROTATION)*YP_PENCIL0(ICOLL)
                   WRITE(*,*) 'INFO> Changed tilt1  ICOLL  to  ANGLE  ',ICOLL, C_TILT(1)
!
                   CHANGED_TILT1(ICOLL) = .true.
           ELSEIF (.NOT. CHANGED_TILT2(ICOLL) .AND. MIRROR.LT.0.) THEN
                   C_TILT(2) = -1.*(XP_PENCIL0(ICOLL)*COS(C_ROTATION)+SIN(C_ROTATION)*YP_PENCIL0(ICOLL))
                   WRITE(*,*) 'INFO> Changed tilt2  ICOLL  to  ANGLE  ',ICOLL, C_TILT(2)
!
                   CHANGED_TILT2(ICOLL) = .true.
           ENDIF
         ELSE
           C_TILT(1) = 0.
           C_TILT(2) = 0.
           CHANGED_TILT1(ICOLL) = .true.
           CHANGED_TILT2(ICOLL) = .true.
         ENDIF
!------------------- end pencil beam stuff-----------------

!-----------daniele------------------
! debugged assignation of p after passage in the crystal
!-------------

!         p_in(J) = (1 + dpop) * p0      !daniele
         p_in(J) = P                   !daniele
         s_in(J) = s_in(J) + S
!
          if (nabs.eq.1) then
             fracab = fracab + 1
             part_abs(j) = ie
             part_abs_turn(j) = iturn
             lint(j) = zlm
             if (dowrite_impact) then
               write(48,'(i4,(1x,f6.3),(1x,f8.6),4(1x,e19.10),i2, &
     &         2(1x,i7))')                                         &
     &         icoll,c_rotation,                                    &
     &         s,                 &
     &         x_in(j)*1d3, xp_in(j)*1d3, y_in(j)*1d3, yp_in(j)*1d3,  &
     &         nabs,name(j),iturn
             endif
             x = 99.99d-3
             z = 99.99d-3
          endif
       ENDIF



!-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o-~-o
! valentina First impact file
!
!                write(9999,*) "dowrite impact value", dowrite_impact
             if(flagsec(j).eq.0 .and. PROC(1:3) .ne. 'out') then
               flagsec(j)=1
               if (dowrite_impact) then
               write(39,'(i5,1x,i7,1x,i2,1x,i1,2(1x,f7.6),8(1x,e17.9))') &
     &               name(j),iturn,icoll,nabs,                          &
     &               s_impact,          &
     &               s, &
     &               x_in0(j),xp_in0(j),y_in0(j),yp_in0(j), &
     &               x_in(j),xp_in(j),y_in(j),yp_in(j)
               endif
             endif
!
!++  End of check for particles not being lost before   (see @330)
!
        ENDIF
!++  End of loop over all particles
!
 777  END DO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!        if (nhit .gt. 0.) then
!          WRITE(*,*) 'Collimator:                    ',ICOLL
!          WRITE(*,*) 'Number of particles:           ', Nev
!          WRITE(*,*) 'Number of particle hits:       ', Nhit
!          WRITE(*,*) 'Number of absorped particles:   ', fracab
!          WRITE(*,*) 'Number of escaped particles:    ', Nhit-fracab
!          WRITE(*,*) 'Fraction of absorbed particles:', 100.*fracab/Nhit
!            WRITE(*,*)'Fraction of channeled particles:',100*n_chan/Nhit
!            WRITE(*,*)'Fraction of VR particles:       ',100*n_VR/Nhit
!            WRITE(*,*)'Fraction of amorphous process:  ',100*n_amorphous
!     1/Nhit
!        endif
!
      endif  !collimator with length = 0
      return
      end


!.**************************************************************************
!     SUBROUTINE FOR THE MOVEMENTS OF THE PARTICLES IN THE CRYSTAL
!.**************************************************************************
      SUBROUTINE CRYST(IS,x,xp,y,yp,PC,Length)

      use mod_ranlux
      use mod_funlux
!
!     Simple tranport protons in crystal 2
!-----------------------------------------------------------C
!      J -   number of element                              C
!      S - longitudinal coordinate
!      IS -   number of substance 1-4: Si,W,C,Ge(110)       C
!      x,xp,y,yp - coordinates at input of crystal          C
!      PC -   momentum of particle*c [GeV]                  C
!      W  -   weigth of particle                            C
!-----------------------------------------------------------C
!
!
      IMPLICIT none
!
!
      double precision Rcurv,Length,C_xmax,C_ymax !crystal geometrical parameters
                                                  ! [m],[m],[m],[m],[rad]
      double precision ymax,ymin       !crystal geometrical parameters
      double precision s_length             !element length along s
      double precision Alayer               !amorphous layer [m]
      integer C_orient                      !crystal orientation
      integer IS                            !index of the material
!      integer counter
      double precision  DLRI(4),DLYI(4),AI(4),DES(4)!cry parameters:see line~270
      double precision DESt                  ! Daniele: changed energy loss by ionization now calculated and not tabulated
      integer NAM,ZN                        !switch on/off the nuclear
                                            !interaction (NAM) and the MCS (ZN)
      double precision x,xp,y,yp,PC         !coordinates of the particle
                                            ![m],[rad],[m],[rad],[GeV]
      double precision x0,y0                !coordinates of the particle [m]
      double precision s                    !long coordinates of the particle [m]
      double precision a_eq,b_eq,c_eq,Delta !second order equation param.
      double precision Ang_rms, Ang_avr     !Volume reflection mean angle [rad]
      double precision c_v1                 !fitting coefficient
      double precision c_v2                 !fitting coefficient
      double precision Dechan               !probability for dechanneling
      double precision Lrefl, Srefl         !distance of the reflection point [m]
      double precision Vcapt                !volume capture probability
      double precision Chann                !channeling probability
      double precision N_atom               !probability for entering
                                            !channeling near atomic planes
      double precision Dxp                  !variation in angle
      double precision xpcrit               !critical angle for curved crystal[rad]
      double precision xpcrit0              !critical angle for str. crystal [rad]
      double precision Rcrit                !critical curvature radius [m]
      double precision ratio                !x=Rcurv/Rcrit
      double precision Cry_length
      double precision TLdech2              !tipical dechanneling length(1) [m]
      double precision TLdech1              !tipical dechanneling length(2) [m]
      double precision tdech, Ldech,Sdech   !angle, lenght, and S coordinate
                                            !of dechanneling point
      double precision Rlength, Red_S       !reduced length/s coordinate
                                            !(in case of dechanneling)
      double precision Am_length            !Amorphous length
      double precision Length_xs, Length_ys !Amorphous length
      double precision miscut               !miscut angle in rad
      double precision L_chan, tchan
      double precision xp_rel               !xp-miscut angle in mrad
!      REAL RNDM                           !random numbers
!      REAL      rndm4
!      REAL      RAN_GAUSS
      double precision eUm(4)                !maximum potential
      CHARACTER(LEN=50) PROC        !string that contains the physical process

      double precision rho(4),z(4),Ime(4) !Daniele: material parameters for dE/dX calculation
      double precision k,re,me,mp !Daniele: parameters for dE/dX calculation (const,electron radius,el. mass, prot.mass)
      double precision enr,mom,betar,gammar,bgr !Daniele: energy,momentum,beta relativistic, gamma relativistic
      double precision Tmax,plen !Daniele: maximum energy tranfer in single collision, plasma energy (see pdg)
      double precision anuc_cry2(4),aTF,dP,const_dech,u1,xpin,ypin
      double precision alpha !Daniele: par for new chann prob
      double precision Pvr !Daniele: prob for VR->AM transition

      double precision emr_cry(4)

!      common /Par_Cry1/ Cry_length, Rcurv,C_xmax,C_ymax,Alayer,C_orient
      common /Par_Cry1/ Rcurv,C_xmax,C_ymax,Alayer,C_orient
      common /miscut/ miscut
      common /Proc2/PROC
      COMMON/NPC/     NAM,ZN
      COMMON/CRYS/    DLRI,DLYI,AI,DES
      COMMON/eUc/     eUm  !
      common/ion/rho,z,Ime,k,re,me,mp,anuc_cry2,emr_cry
      common/ion2/enr,mom,gammar,betar,bgr,Tmax,plen
      common/dech/aTF,dP,u1

!      common/utils/ counter
!
!
      NAM=1 !switch on/off the nuclear interaction (NAM) and the MCS (ZN)
      ZN=1

!-------------Daniele: dE/dX and dechanneling length calculation--------------------


      mom=PC*1.0d3 ! [GeV/c] -> [MeV/c]
      enr=(mom*mom+mp*mp)**0.5 ! [MeV]
      gammar=enr/mp
      betar=mom/enr
      bgr=betar*gammar

      Tmax=(2.0d0*me*bgr**2)/(1.0d0+2*gammar*me/mp+(me/mp)**2) ![MeV]

      plen=((rho(IS)*z(IS)/anuc_cry2(IS))**0.5)*28.816d-6 ![MeV]

!      DESt=((k*z(IS))/(anuc_cry2(IS)*betar**2))*
!     + (0.5*log((2.0d0*me*bgr*bgr*Tmax)/(Ime(IS)*Ime(IS)))
!     + -betar**2.0-log(plen/Ime(IS))-log(bgr)+0.5);

!      DESt=DESt*rho(IS)*0.1 ![GeV/m]

      const_dech=(256.0/(9.0*(4.D0*DATAN(1.D0))**2))* &
     & (1.0/(log(2.0*me*gammar/Ime(IS))-1.0))*((aTF*dP)/(re*me))  ![m/MeV]

      const_dech=const_dech*1.0d3    ![m/GeV]

!      write(*,*)DESt, const_dech

!----------------------------------------------------------


!      miscut=0.001000
!
!      write(*,*)"last miscut angle =",miscut
!      write(*,*) 'enter crystal subroutine'
!      write(*,*) 'particle energy Gev :', PC
!      write(*,*) 'x_initial :', x
!      write(*,*) 'Length [m]:', Length
!      write(*,*) 'Random:', rndm4()
!      write(*,*)'xp',xp,'x',x , 's', s
      s=0
!      write(*,*) "s_length", Rcurv, length
      s_length=Rcurv*(sin(length/Rcurv)) !
      L_chan=length

      if ( miscut .lt. 0 &
     &     .and. x .gt. 0 & !should be useless
     &     .and. x .lt. -length*tan(miscut)) then
            L_chan=-x/sin(miscut)
      endif
      tchan=L_chan/Rcurv
      xp_rel=xp-miscut
!  FIRST CASE: p don't interact with crystal
      ymin = - C_ymax/2
      ymax =  C_ymax/2
      IF (y.LT.ymin .or. y.GT.ymax .or. x.gt.C_xmax) THEN
        x = x+xp*s_length
        y = y+yp*s_length
        PROC='out'
        GOTO 111
! SECOND CASE: p hits the amorphous layer
      ELSEIF ( (x.LT.Alayer) .or.  ((y-ymin).LT.Alayer) .or. ((ymax-y).lt.Alayer)  ) THEN
        x0=x
        y0=y
        a_eq=(1+(xp)**2)
        b_eq=(2*(x)*(xp)-2*(xp)*Rcurv)
        c_eq=(x)**2-2*(x)*Rcurv
        Delta=b_eq**2-4*a_eq*c_eq
        s=((-b_eq+sqrt(Delta))/(2*a_eq))
        if (s .ge. s_length) s=s_length
        x=(xp)*s+x0
        Length_xs=sqrt((x-x0)**2+s**2)
        if ( (yp .ge.0 .and. (y+yp*s).le.ymax)) then
          Length_ys = yp*Length_xs
        elseif (yp.lt.0 .and. (y+yp*s).ge. ymin) then
          Length_ys = yp*Length_xs
        else
          s=(ymax-y)/yp
          Length_ys = sqrt((ymax-y)**2+s**2)
          x=x0+xp*s
          Length_xs=sqrt((x-x0)**2+s**2)
        endif
        Am_length   = sqrt(Length_xs**2+Length_ys**2)
        s=s/2
        x=x0+xp*s
        y=y0+yp*s
        PROC='AM'
!        CALL MOVE_AM_(IS,NAM,Am_Length,DES(IS),DLYi(IS),DLRi(IS),xp,yp
!     + ,PC)
!        write(*,*) "Am_length", Am_Length
        CALL CALC_ION_LOSS_CRY(IS,PC,AM_Length,DESt)
        CALL MOVE_AM_(IS,NAM,Am_Length,DESt,DLYi(IS),DLRi(IS),xp,yp,PC)
        x=x+xp*(s_length-s)
        y=y+yp*(s_length-s)
        GOTO 111
      ELSEIF ((x.GT.(C_xmax-Alayer)) .and. x.LT.(C_xmax)  ) THEN
        PROC='AM'
!        CALL MOVE_AM_(IS,NAM,s_length,DES(IS),DLYi(IS),DLRi(IS), xp,yp,PC)
!        write(*,*) "s_length 1137", s_length
        CALL CALC_ION_LOSS_CRY(IS,PC,s_length,DESt)
        CALL MOVE_AM_(IS,NAM,s_length,DESt,DLYi(IS),DLRi(IS), xp,yp,PC)
        WRITE(*,*)'Fix here!'
        GOTO 111
      ELSE
        DESt = 0   ! temporary solution in case DESt is not initialized
      END IF
!
! THIRD CASE: the p interacts with the crystal.
!. Define typical angles/probabilities for orientation 110
!
      xpcrit0 = (2.e-9*eUm(IS)/PC)**0.5       ! critical angle (rad) for
                                              ! straight crystals
      Rcrit  = PC/(2.e-6*eUm(IS))*AI(IS)      ! critical curvature radius [m]
                                              ! if R>Rcritical=>no channeling is
                                              ! possible (ratio<1)
      ratio = Rcurv/Rcrit                     ! parameter Rcry/Rcritical
!      write(*,*) "Critical Radius: ",Rcrit
      xpcrit = xpcrit0*(Rcurv-Rcrit)/Rcurv    ! critical angle for curved crystal
!----------------valentina approx-----------
!      xpcrit = xpcrit0*(1-(Rcrit/Rcurv))**0.5
!      write(*,*)(Rcurv-Rcrit)/Rcurv,(1-(Rcrit/Rcurv))**0.5

                                              ! NB: if ratio<1 => xpcrit<0
      c_v1 = 1.7                              ! fitting coefficient ??!
      c_v2 = -1.5                             ! fitting coefficient ???
      if (ratio .le. 1.) then                 ! case 1:no possibile channeling
        Ang_rms = c_v1*0.42*xpcrit0*sin(1.4*ratio)  ! rms scattering
        Ang_avr = c_v2*xpcrit0*0.05*ratio           ! average angle reflection
        Vcapt = 0.0                                 ! probability of VC
        elseif (ratio .le. 3) then              ! case 2: strongly bent xstal
          Ang_rms = c_v1*0.42*xpcrit0*sin(1.571*0.3*ratio+0.85)! rms scattering
          Ang_avr = c_v2*xpcrit0*(0.1972*ratio-0.1472)  ! avg angle reflection
!          Vcapt   = 0.01*(ratio-0.7)/(PC**2)
          Vcapt   = 0.0007*(ratio-0.7)/PC**0.2 !correction by sasha drozdin/armen
          !K=0.00070 is taken based on simulations using CATCH.f (V.Biryukov)
        else                                       ! case 3: Rcry >> Rcrit
          Ang_rms = c_v1*xpcrit0*(1./ratio)        !
          Ang_avr = c_v2*xpcrit0*(1.-1.6667/ratio) ! average angle for VR
!          Vcapt = 0.01*(ratio-0.7)/(PC**2)        ! probability for VC
          Vcapt = 0.0007*(ratio-0.7)/PC**0.2  !correction by sasha drozdin/armen
          ! K=0.0007 is taken based on simulations using CATCH.f (V.Biryukov)
      endif
!c----------------valentina approx-----------
!      Ang_avr=-(xpcrit+xpcrit0)
!      Ang_rms=(xpcrit0-xpcrit)/2
!c-----------end valentina approx--------------

!      write(*,*) "Rcrit" , Rcrit,"Rcurv",Rcurv,
!     c "Ratio: ",ratio,"average VR angle:", ang_avr*1e6,"+-",
!     c ang_rms*1e6, "ang crit:", xpcrit0*1e6,xpcrit*1e6
!
      if(C_orient .eq. 2) then
        Ang_avr = Ang_avr * 0.93                     ! for (111)
        Ang_rms = Ang_rms * 1.05
        xpcrit  = xpcrit * 0.98
      endif
!
!. case 3-1: channeling
!      IF (abs(xp_rel) .lt. xpcrit) THEN              ! if R' < R'c (ok CH) (1)
!        Chann  = (xpcrit**2-xp_rel**2)**0.5/xpcrit0  ! probability of CH/VC  OCCHIO DANIELE LINEA COMMENTATA PER PROVARE LA SUCCESSIVA PROB.
!        N_atom = 0.1                                ! probability of entering OCCHIO DANIELE, prob cambiata in accordo con nuova chann
!--------------DAN CHAN prob------
!         alpha = xp_rel/xpcrit
!         Chann = ((0.64-(1/ratio)*(1/ratio)*0.64)**0.5)*     !DANIELE saturation at 80%
!     &           (1-alpha*alpha)**0.5

!         Chann = ((0.8-(1/ratio)*(1/ratio)*0.8)**0.5)*     !DANIELE saturation at 90%
!     &           (1-alpha*alpha)**0.5

!         N_atom=0.14                                      !DANIELE for sat. at 80%
!         N_atom=0.1                                       !DANIELE for sat. at 90%

!      IF (abs(xp_rel) .lt. xpcrit) THEN
!         Chann = ((0.9-(1/ratio)*(1/ratio)*0.9)**0.5)*     !DANIELE saturation at 95%
!     &           (1-alpha*alpha)**0.5
!        N_atom = 0.1

      IF (abs(xp_rel) .lt. xpcrit) THEN
         alpha = xp_rel/xpcrit
         Chann = ((0.9-alpha*alpha*0.9)**0.5)*(1-(1/ratio))**0.5         !DANIELE saturation at 95%
         N_atom = 0.1

!      IF (abs(xp_rel) .lt. xpcrit) THEN
!         alpha = xp_rel/xpcrit
!         Chann = ((0.8-alpha*alpha*0.8)**0.5)*(1-(1/ratio))**0.5         !DANIELE saturation at 90%
!         N_atom = 0.1

!--------------end DAN CHAN prob------




                                                     ! close to atomic planes
        IF (rndm4() .le. Chann) then      ! if they can channel: 2 options
                                          ! option 1:channeling
!          TLdech1= 0.00054*PC*(1.-1./ratio)**2 ! calculate dechanneling length
!          TLdech1= 0.0005*PC*(1.-1./ratio)**2 !calculate tipical dech. length(m)
          TLdech1= const_dech*PC*(1.-1./ratio)**2 !Daniele: updated calculate tipical dech. length(m)
          IF (rndm4() .le. N_atom) then
!            TLdech1= 0.000004*PC*(1.-1./ratio)**2! calculate tipical dechanneling
                                                  !length near atomic planes(m)
           !next line new from sasha
!            TLdech1= 2.0e-6*PC*(1.-1./ratio)**2  ! dechanneling length (m)
            TLdech1= (const_dech/200.d0)*PC*(1.-1./ratio)**2  ! Daniele: updated dechanneling length (m)

                               !for short crystal for high amplitude particles
          ENDIF

!          TLdech1=TLdech1/100 !!!!CHECK

          Dechan = -log(rndm4())                 ! probability of dechanneling
          Ldech  = TLdech1*Dechan                ! actual dechan. length
                     ! careful: the dechanneling lentgh is along the trajectory
                     ! of the particle -not along the longitudinal coordinate...
          if(Ldech .LT. L_chan) THEN
           PROC='DC'                  !
            Dxp= Ldech/Rcurv             ! change angle from channeling [mrad]
            Sdech=Ldech*cos(miscut+0.5*Dxp)

            x  = x+ Ldech*(sin(0.5*Dxp+miscut))   ! trajectory at channeling exit
            xp = xp + Dxp + 2.0*(rndm4()-0.5)*xpcrit
            y= y + yp * Sdech
            write(*,*) "Ldech", Ldech
            write(*,*) "DESt", DESt
!            CALL CALC_ION_LOSS_CRY(IS,PC,Ldech,DESt)
            PC = PC - 0.5*DESt*Ldech          ! energy loss to ionization while in CH [GeV]

            x = x + 0.5*(s_length-Sdech)*xp
            y = y + 0.5*(s_length-Sdech)*yp

            write(*,*) "s_length-Sdech", s_length-Sdech
!            CALL CALC_ION_LOSS_CRY(IS,PC,s_length-Sdech,DESt)
            CALL MOVE_AM_(IS,NAM,s_length-Sdech,DESt,DLYi(IS),DLRi(IS),xp,yp,PC)
           !next line new from sasha
!            PC = PC - 0.5*DES(IS)*y          ! energy loss to ionization [GeV]
            x = x + 0.5*(s_length-Sdech)*xp
            y = y + 0.5*(s_length-Sdech)*yp
          else
            PROC='CH'
            xpin=XP
            ypin=YP

!            write(*,*) PROC
!            write(*,*) "MOVE CH", IS, NAM, L_chan, X, XP, YP, PC, Rcurv, Rcrit
            CALL MOVE_CH_(IS,NAM,L_chan,X,XP,YP,PC,Rcurv,Rcrit)  !daniele:check if a nuclear interaction happen while in CH

!            write(*,*) PROC

            if(PROC(1:2).ne.'CH')then             !daniele: if an nuclear interaction happened, move until the middle with initial xp,yp
            x = x + 0.5 * L_chan * xpin           !then propagate until the "crystal exit" with the new xp,yp
            y = y + 0.5 * L_chan * ypin           !accordingly with the rest of the code in "thin lens approx"
            x = x + 0.5 * L_chan * XP
            y = y + 0.5 * L_chan * YP
!            write(*,*) "Length", Length
            CALL CALC_ION_LOSS_CRY(IS,PC,Length,DESt)
            PC = PC - DESt*Length       ! energy loss to ionization [GeV]
            else
            Dxp= L_chan/Rcurv + 0.5*RAN_GAUSS(1.0d0)*xpcrit ! change angle[rad]
            xp = Dxp
            !next line new from sasha
            x  = x+ L_chan*(sin(0.5*Dxp+miscut)) ! trajectory at channeling exit
!            xp = xp + Dxp + 2.0*(rndm4()-0.5)*xpcrit
            y = y + s_length * yp
!            !next line new from sasha
!            PC = PC - 0.5*DES(IS)*Length       ! energy loss to ionization [GeV]
!            write(*,*) "Length", Length
            CALL CALC_ION_LOSS_CRY(IS,PC,Length,DESt)
            PC = PC - 0.5*DESt*Length       ! energy loss to ionization [GeV]
            endif
          endif
        ELSE                                   !option 2: VR
                                               ! good for channeling
                                               ! but don't channel         (1-2)
          PROC='VR'                            !volume reflection at the surface
!          Dxp=0.5*(xp_rel/xpcrit+1)*Ang_avr
!          xp=xp+Dxp+Ang_rms*RAN_GAUSS(1.)
            !next line new from sasha
          xp=xp+0.45*(xp/xpcrit+1)*Ang_avr
          x = x + 0.5*s_length * xp
          y = y + 0.5*s_length * yp
!          CALL MOVE_AM_(IS,NAM,s_length,DES(IS),DLYi(IS),DLRi(IS),  xp ,yp,PC)
!          write(*,*) "s_length 1319", s_length
          CALL CALC_ION_LOSS_CRY(IS,PC,s_length,DESt)
          CALL MOVE_AM_(IS,NAM,s_length,DESt,DLYi(IS),DLRi(IS),xp ,yp,PC)
          x = x + 0.5*s_length * xp
          y = y + 0.5*s_length * yp
        ENDIF                                    !
! case 3-2: no good for channeling. check if the  can VR
      ELSE
        Lrefl =  (xp_rel)*Rcurv                  ! distance of refl. point [m]
!        Srefl = sin(xp) * Lrefl
        Srefl = sin(xp_rel/2+miscut) * Lrefl
        if(Lrefl .gt. 0. .and. Lrefl .lt. Length) then
                ! VR point inside
                !2 options: volume capture and volume reflection
          IF (rndm4() .gt. Vcapt .or. ZN .eq. 0.) THEN   !opt. 1: VR
           PROC='VR'
            x = x + xp * Srefl
            y = y + yp * Srefl
            Dxp= Ang_avr
            xp = xp + Dxp + Ang_rms*RAN_GAUSS(1.0d0)
            x = x + 0.5* xp * (s_length - Srefl)
            y = y + 0.5* yp * (s_length - Srefl)     !
!            CALL MOVE_AM_(IS,NAM,s_length-Srefl,DES(IS),DLYi(IS),DLRi(IS),xp ,yp,PC)
!            write(*,*) "s_length-Srefl", s_length-Srefl
            CALL CALC_ION_LOSS_CRY(IS,PC,s_length-Srefl,DESt)
            CALL MOVE_AM_(IS,NAM,s_length-Srefl,DESt,DLYi(IS),DLRi(IS),xp ,yp,PC)
            x = x + 0.5 * xp * (s_length - Srefl)
            y = y + 0.5 * yp * (s_length - Srefl)
          ELSE                                      !opt 2: VC
            x = x + xp * Srefl
            y = y + yp * Srefl
!            TLdech2= 0.00011*PC**0.25*(1.-1./ratio)**2 ! dechanneling length(m)
!            Dechan = log(1.-rndm4())
!            Ldech  = -TLdech2*Dechan
           !next 2 lines new from sasha - different dechanneling
           !probability
 !           TLdech2= 0.01*PC*(1.-1./ratio)**2   ! typical dechanneling length(m)
 !           Ldech  = 0.005*TLdech2*(sqrt(0.01-log(rndm4())) -0.1)**2 ! DC length
            TLdech2= (const_dech/10.0d0)*PC*(1.-1./ratio)**2   ! Daniele: updated typical dechanneling length(m)
            Ldech  = TLdech2*(sqrt(0.01-log(rndm4())) -0.1)**2 ! daniele: updated DC length
            tdech=Ldech/Rcurv
            Sdech=Ldech*cos(xp+0.5*tdech)
            IF(Ldech .LT. (Length-Lrefl)) then
              PROC='DC'
              Dxp= Ldech/Rcurv + 0.5*ran_gauss(1.0d0)*xpcrit
              x  = x+ Ldech*(sin(0.5*Dxp+xp))   ! trajectory at channeling exit
              y = y + Sdech * yp
              xp =  Dxp
              Red_S = s_length-Srefl -Sdech
              x = x + 0.5 * xp * Red_S
              y = y + 0.5 * yp * Red_S
!              write(*,*) "Srefl", Srefl
              CALL CALC_ION_LOSS_CRY(IS,PC,Srefl,DESt)
              PC=PC - DESt * Srefl !Daniele: "added" energy loss before capture
!              write(*,*) "Sdech", Sdech
              CALL CALC_ION_LOSS_CRY(IS,PC,Sdech,DESt)
              PC=PC - 0.5 * DESt * Sdech !Daniele: "added" energy loss while captured
!              CALL MOVE_AM_(IS,NAM,Red_S,DES(IS),DLYi(IS),DLRi(IS),xp,yp,PC)
!              write(*,*) "Red_S", Red_S
              CALL CALC_ION_LOSS_CRY(IS,PC,Red_S,DESt)
              CALL MOVE_AM_(IS,NAM,Red_S,DESt,DLYi(IS),DLRi(IS),xp,yp,PC)
              x = x + 0.5 * xp * Red_S
              y = y + 0.5 * yp * Red_S
            else
              PROC='VC'
              Rlength = Length-Lrefl
              tchan = Rlength / Rcurv
              Red_S=Rlength*cos(xp+0.5*tchan)
              write(*,*) "Lrefl", Lrefl
              CALL CALC_ION_LOSS_CRY(IS,PC,Lrefl,DESt)
              PC=PC - DESt*Lrefl  !Daniele: "added" energy loss before capture
              xpin=XP
              ypin=YP
              CALL MOVE_CH_(IS,NAM,Rlength,X,XP,YP,PC,Rcurv,Rcrit)  !daniele:check if a nuclear interaction happen while in CH

              if(PROC(1:2).ne.'VC')then             !daniele: if an nuclear interaction happened, move until the middle with initial xp,yp
              x = x + 0.5 * Rlength * xpin           !then propagate until the "crystal exit" with the new xp,yp
              y = y + 0.5 * Rlength * ypin           !accordingly with the rest of the code in "thin lens approx"
              x = x + 0.5 * Rlength * XP
              y = y + 0.5 * Rlength * YP
!              write(*,*) "Rlength", Rlength
              CALL CALC_ION_LOSS_CRY(IS,PC,Rlength,DESt)
              PC=PC - DESt*Rlength
              else
              Dxp = (Length-Lrefl)/Rcurv
              x  = x+ sin(0.5*Dxp+xp)*Rlength     ! trajectory at channeling exit
              y = y + red_S * yp
              xp =  Length/Rcurv + 0.5*ran_gauss(1.0d0)*xpcrit ! [mrad]
!              write(*,*) "Rlength", Rlength
              CALL CALC_ION_LOSS_CRY(IS,PC,Rlength,DESt)
              PC=PC - 0.5*DESt*Rlength  !Daniele: "added" energy loss once captured
              endif
            endif
          ENDIF
!.  case 3-3: move in amorphous substance (big input angles)---------------  MODIFIED FOR TRANSITION VRAM DANIELE
        else
           if(xp_rel .GT. L_chan/Rcurv+2.0*xpcrit .or. xp_rel .lt. -xpcrit) then
             PROC='AM'
             x = x + 0.5 * s_length * xp
             y = y + 0.5 * s_length * yp
            if(ZN .gt. 0) then
!           CALL MOVE_AM_(IS,NAM,s_length,DES(IS),DLYi(IS),DLRi(IS), xp,yp,PC)
!              write(*,*) "s_length 1421", s_length
             CALL CALC_ION_LOSS_CRY(IS,PC,s_length,DESt)
             CALL MOVE_AM_(IS,NAM,s_length,DESt,DLYi(IS),DLRi(IS), xp,yp,PC)
            endif
            x = x + 0.5 * s_length * xp
            y = y + 0.5 * s_length * yp
          else
!            Pvr=0.5*erf((xp_rel-(L_chan/Rcurv)-xpcrit)/(2.0*4.0*xpcrit*xpcrit)**0.5)+0.5
!            Pvr=0.2+(0.6/(2.0*xpcrit))*(xp_rel-(L_chan/Rcurv))
            Pvr=((xp_rel-(L_chan/Rcurv))/(2.0*xpcrit))
            if(rndm4() .gt. Pvr) then
            PROC='TRVR'
            x = x + xp * Srefl
            y = y + yp * Srefl
!            Dxp=(2.0*Ang_avr-Ang_rms)/(2.0*L_chan/Rcurv+2.0*xpcrit)*
!     +         xp_rel+L_chan/Rcurv*(2.0*Ang_avr-Ang_rms)/
!     +         (2.0*L_chan/Rcurv+2.0*xpcrit)-Ang_avr
            Dxp=-3.0*Ang_rms*xp_rel/(2.0*xpcrit)+Ang_avr+(3.0*Ang_rms*(L_chan/Rcurv)/(2.0*xpcrit))
            write(*,*) xp_rel, Dxp, Ang_avr, Ang_rms
!            xp=xp+Dxp+Ang_rms*RAN_GAUSS(1.)
            xp=xp+Dxp
            x=x+0.5*xp*(s_length-Srefl)
            y=y+0.5*yp*(s_length-Srefl)
!            write(*,*) "s_length-Srefl", s_length-Srefl
            CALL CALC_ION_LOSS_CRY(IS,PC,s_length-Srefl,DESt)
            CALL MOVE_AM_(IS,NAM,s_length-Srefl,DESt,DLYi(IS),DLRi(IS),xp ,yp,PC)
            x = x + 0.5 * xp * (s_length - Srefl)
            y = y + 0.5 * yp * (s_length - Srefl)
            else
            PROC='TRAM'
            x = x + xp * Srefl
            y = y + yp * Srefl
!            Dxp=-1.0*(13.6/PC)*SQRT(s_length/DLRi(IS))*xp_rel/
!     +      (2.0*xpcrit)+(13.6/PC)*SQRT(s_length/DLRi(IS))*
!     +      (1.0+(L_chan/Rcurv)/(2.0*xpcrit))
            Dxp=-1.0*(13.6/PC)*SQRT(s_length/DLRi(IS))*0.001*xp_rel/    &
     &          (2.0*xpcrit)+(13.6/PC)*SQRT(s_length/DLRi(IS))*0.001*   &
     &          (1.0+(L_chan/Rcurv)/(2.0*xpcrit))
!            write(*,*) xp_rel, Dxp
!            xp=xp+Dxp+(13.6/PC)*SQRT(s_length/DLRi(IS))*RAN_GAUSS(1.)
            xp=xp+Dxp
            x=x+0.5*xp*(s_length-Srefl)
            y=y+0.5*yp*(s_length-Srefl)
!            write(*,*) "s_length-Srefl", s_length-Srefl
            CALL CALC_ION_LOSS_CRY(IS,PC,s_length-Srefl,DESt)
            CALL MOVE_AM_(IS,NAM,s_length-Srefl,DESt,DLYi(IS),DLRi(IS),xp ,yp,PC)
            x = x + 0.5 * xp * (s_length - Srefl)
            y = y + 0.5 * yp * (s_length - Srefl)
            endif
          endif
        endif
       ENDIF

!      if (counter .eq. 0) then
111   write(833,*)'crystal parameters:\n Length:',Length, '\n Rcurv:'   &
     & , Rcurv ,'\n Critical Radius:', Rcrit, 'ratio',ratio             &
     &, '\n Critical angle for straight:',                              &
     & xpcrit0,'\n critical angle for curved crystal:', xpcrit,' \n Leng&
     &th:', Length, '\n xmax:', C_xmax, ' ymax:', ymax, '  C_orient: '  &
     &, C_orient                                                        &
     &, '\n Avg angle reflection:', Ang_avr, '\n full channeling angle: &
     &',(Length/Rcurv)
!      counter=1
!      endif
!
!      WRITE(*,*)'Crystal process: ',PROC
!      write(*,*)'xp_final :', xp
!      WRITE(*,*)'Crystal process: ',PROC,'Chann Angle',Ch_angle/1000,   +
!     1'Critical angle: ', xpcrit/1000
!     2 DLRI(IS),DLYI(IS),AI(IS),DES(IS),eUm(IS),IS,ZN,NAM,C_orient     !,W
      END

!.**************************************************************************
!     subroutine for the calculazion of the energy loss by ionization
!.**************************************************************************
      SUBROUTINE CALC_ION_LOSS_CRY(IS,PC,DZ,EnLo)

      use mod_ranlux
      use mod_funlux

      IMPLICIT none
      integer IS
      double precision PC,DZ,EnLo
      double precision rho(4),z(4),Ime(4) !Daniele: material parameters for dE/dX calculation
      double precision k,re,me,mp !Daniele: parameters for dE/dX calculation (const,electron radius,el. mass, prot.mass)
      double precision enr,mom,betar,gammar,bgr !Daniele: energy,momentum,beta relativistic, gamma relativistic
      double precision Tmax,plen !Daniele: maximum energy tranfer in single collision, plasma energy (see pdg)
      double precision anuc_cry2(4)
      double precision emr_cry(4)
      double precision thl,Tt,cs_tail,prob_tail
      double precision ranc
!      REAL RNDM4

      common/ion/rho,z,Ime,k,re,me,mp,anuc_cry2,emr_cry
      common/ion2/enr,mom,gammar,betar,bgr,Tmax,plen

      write(*,*) "Entering CALC_ION_LOSS_CRY"
!      write(*,*) "thl", DZ
       thl= 4.0d0*k*z(IS)*DZ*100.0d0*rho(IS)/(anuc_cry2(IS)*betar**2) ![MeV]

!      write(*,*) "EnLo", me, bgr, plen, Ime(IS)
       EnLo=((k*z(IS))/(anuc_cry2(IS)*betar**2))*          &
     & (0.5*log((2.0d0*me*bgr*bgr*Tmax)/(Ime(IS)*Ime(IS))) &
     & -betar**2.0-log(plen/Ime(IS))-log(bgr)+0.5);

!      write(*,*) "EnLo"
       EnLo=EnLo*rho(IS)*0.1*DZ ![GeV]

!      write(*,*) "Tt", EnLo, thl
       Tt=EnLo*1000.0d0+thl  ![MeV]

!      write(*,*) "cs_tail", k, z(IS), anuc_cry2(IS), betar, Tt, Tmax, gammar, mp
       cs_tail=((k*z(IS))/(anuc_cry2(IS)*betar**2))* &
     & ((0.5*((1.0d0/Tt)-(1.0d0/Tmax)))-             &
     & (log(Tmax/Tt)*(betar**2)/(2.0d0*Tmax))+       &
     & ((Tmax-Tt)/(4.0d0*(gammar**2)*(mp**2))))

!      write(*,*) "prob_tail"
       prob_tail=cs_tail*rho(IS)*DZ*100.0d0;

!      write(*,*) "ranc"
       ranc=dble(rndm4())

       if(ranc.lt.prob_tail)then
         EnLo=((k*z(IS))/(anuc_cry2(IS)*betar**2))*          &
     &   (0.5*log((2.0d0*me*bgr*bgr*Tmax)/(Ime(IS)*Ime(IS))) &
     &   -betar**2.0-log(plen/Ime(IS))-log(bgr)+0.5+         &
     &   (TMax**2)/(8.0d0*(gammar**2)*(mp**2)));

         EnLo=EnLo*rho(IS)*0.1 ![GeV/m]

       else
         EnLo=EnLo/DZ  ![GeV/m]
       endif

!      write(*,*)cs_tail,prob_tail,ranc,EnLo*DZ
       write(*,*) "CALC_ION_LOSS_CRY done"

      RETURN
      END


!.**************************************************************************
!     subroutine for the movement in the amorphous
!.**************************************************************************
      SUBROUTINE MOVE_AM_(IS,NAM,DZ,DEI,DLY,DLr, XP,YP,PC)

      use mod_ranlux
      use mod_funlux

!. Moving in amorphous substance...........................
      IMPLICIT none
      integer IS,NAM
      double precision DZ,DEI,DLY,DLr, XP,YP,PC
      double precision DLAI(4),SAI(4),DES(4)
      double precision DLRI(4),DLYI(4),AI(4)
      double precision AM(30),QP(30),NPAI
      double precision Dc(4),eUm(4)
      double precision DYA,W_p
!      REAL RNDM4
!         REAL      RAN_GAUSS
      CHARACTER(LEN=50) PROC              !string that contains the physical process
      COMMON /ALAST/DLAI,SAI
      COMMON/CRYS/ DLRI,DLYI,AI,DES
      common /Proc2/PROC


!------- adding variables --------


      double precision cprob_cry(0:5,1:4)
      double precision cs(0:5,1:4)
      double precision csref_cry(0:5,1:4)
      double precision freep(4)
      double precision anuc_cry(4)
      double precision collnt_cry(4)
      double precision bn(4)
      double precision bnref_cry(4)

      double precision freeco_cry,ppel,ppsd,pptot,pptref_cry,pptco_cry
      double precision pperef_cry,sdcoe_cry,pref_cry,ppeco_cry

      double precision ecmsq,xln15s,bpp,xm2,bsd,t,teta,va,vb,va2,vb2
      double precision tx,tz,r2,zlm,f

      integer i,ichoix
      double precision aran

      common/scat_cry/cprob_cry,csref_cry,collnt_cry,bnref_cry,anuc_cry
      common/scat_cry/pptref_cry,pperef_cry,sdcoe_cry,pref_cry
      common/scat_cry/pptco_cry,ppeco_cry,freeco_cry

      double precision tlcut_cry,hcut_cry(4)
      double precision cgen_cry(200,4),tlow,thigh,ruth_cry
      external ruth_cry

      common/ruth_scat_cry/tlcut_cry,hcut_cry
      common/ruth_scat_cry/cgen_cry,mcurr_cry

      integer length_cry,mcurr_cry
      double precision xran_cry(1)

!-------daniele--------------
! adding new variables to study the energy loss when the routine is called
!---------------------------

      double precision PC_in_dan     !daniele
      integer PROC_dan    !daniele
      double precision xp_in,yp_in,kxmcs,kymcs

      write(*,*) "Entering MOVE_AM"


      PC_in_dan=PC                   !daniele

      xp_in=XP
      yp_in=YP


!---------------daniele------------
! NEW TREATMENT OF SCATTERING ROUTINE BASED ON STANDARD SIXTRACK ROUTINE
!----------------------------------


!------- useful calculations for cross-section and event topology calculation --------------

      ecmsq = 2 * 0.93828d0 * PC
      write(*,*) ecmsq, PC
      xln15s=log(0.15*ecmsq)
! pp(pn) data
!      pptot = pptref_cry *(PC / pref_cry)** pptco_cry
!      ppel = pperef_cry *(PC / pref_cry)** ppeco_cry
!      ppsd = sdcoe_cry * log(0.15d0 * ecmsq)
!      bpp = 8.5d0 + 1.086d0 * log(sqrt(ecmsq))

!new models, see Claudia's thesis

      pptot=0.041084d0-0.0023302d0*log(ecmsq)+0.00031514d0*log(ecmsq)**2
      ppel=(11.7d0-1.59d0*log(ecmsq)+0.134d0*log(ecmsq)**2)/1000
      ppsd=(4.3d0+0.3d0*log(ecmsq))/1000
      bpp=7.156d0+1.439d0*log(sqrt(ecmsq))

!------ distribution for Ruth. scatt.---------

      tlow=tlcut_cry
      mcurr_cry=IS
      thigh=hcut_cry(IS)
!      write(*,*)tlow,thigh
      call funlxp(ruth_cry,cgen_cry(1,IS),tlow,thigh)



!---------- cross-section calculation -----------

!
! freep: number of nucleons involved in single scattering
        freep(IS) = freeco_cry * anuc_cry(IS)**(1d0/3d0)
! compute pp and pn el+single diff contributions to cross-section
! (both added : quasi-elastic or qel later)
        cs(3,IS) = freep(IS) * ppel
        cs(4,IS) = freep(IS) * ppsd
!
! correct TOT-CSec for energy dependence of qel
! TOT CS is here without a Coulomb contribution
        cs(0,IS) = csref_cry(0,IS) + freep(IS) * (pptot - pptref_cry)
        bn(IS) = bnref_cry(IS) * cs(0,IS) / csref_cry(0,IS)
! also correct inel-CS
        cs(1,IS) = csref_cry(1,IS) * cs(0,IS) / csref_cry(0,IS)
!
! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
        cs(2,IS) = cs(0,IS) - cs(1,IS) - cs(3,IS) - cs(4,IS)
        cs(5,IS) = csref_cry(5,IS)
! Now add Coulomb
        cs(0,IS) = cs(0,IS) + cs(5,IS)
! Interaction length in meter
!        xintl(ma) = 0.01d0*anuc(ma)/(fnavo * rho(ma)*cs(0,ma)*1d-24)  !don't need at the moment, take it from pdg


! Calculate cumulative probability

        do i=1,4
          cprob_cry(i,IS)=cprob_cry(i-1,IS)+cs(i,IS)/cs(0,IS)
        end do


!--------- Multiple Coulomb Scattering ---------

      xp=xp*1000
      yp=yp*1000
!      write(*,*)'xp initial:', xp, 'yp initial', yp
!. DEI - dE/dx stoping energy
      PC    = PC - DEI*DZ    ! energy lost because of ionization process[GeV]

!. DYA - rms of coloumb scattering
      DYA = (13.6/PC)*SQRT(DZ/DLr)             !MCS (mrad)
!      write(*,*)'dya=',dya

      kxmcs = DYA*RAN_GAUSS(1.0d0)
      kymcs = DYA*RAN_GAUSS(1.0d0)

      XP = xp+kxmcs
      yp = yp+kymcs

!c     XP    = XP+DYA*RAN_GAUSS(1.)
!c      YP    = YP+DYA*RAN_GAUSS(1.)

      if(NAM .eq. 0) return                    !turn on/off nuclear interactions


!--------- Can nuclear interaction happen? -----


      zlm=-collnt_cry(IS)*log(dble(rndm4()))

!      zlm=0.0

      if(zlm.lt.DZ) then

!--------- Choose nuclear interaction --------


      aran=dble(rndm4())
      i=1
  10  if ( aran.gt.cprob_cry(i,IS) ) then
          i=i+1
          goto 10
      endif
      ichoix=i


!-------- Do the interaction ----------

!      ichoix=2

      if (ichoix.eq.1) then

        PROC = 'absorbed'            !deep inelastic, impinging p disappeared

      endif

      if (ichoix.eq.2) then          ! p-n elastic
        PROC = 'pne'
        t = -log(dble(rndm4()))/bn(IS)
      endif

      if ( ichoix .eq. 3 ) then   !p-p elastic
        PROC='ppe'
        t = -log(dble(rndm4()))/bpp
      endif

      if ( ichoix .eq. 4 ) then   !single diffractive
        PROC='diff'
        xm2 = exp( dble(rndm4()) * xln15s )
        PC = PC  *(1.d0 - xm2/ecmsq)
          if ( xm2 .lt. 2.d0 ) then
             bsd = 2.d0 * bpp
           elseif (( xm2 .ge. 2.d0 ).and. ( xm2 .le. 5.d0 )) then
                bsd = (106.d0-17.d0*xm2) *  bpp / 36.d0
           elseif ( xm2 .gt. 5.d0 ) then
                bsd = 7.d0 * bpp / 12.d0
          endif
        t = -log(dble(rndm4()))/bsd
      endif

      if ( ichoix .eq. 5 ) then
        PROC='ruth'
           length_cry=1
           call funlux( cgen_cry(1,IS) ,xran_cry,length_cry)
           t=xran_cry(1)
      endif

!---------- calculate the related kick -----------


      if ( ichoix .eq. 4) then
        teta = sqrt(t)/PC_in_dan                !DIFF has changed PC!!!
      else
        teta = sqrt(t)/PC
      endif


! 100  va  =2d0*rndm4()-1d0
!      vb = dble(rndm4())
!      va2 = va*va
!      vb2 = vb*vb
!      r2 = va2 + vb2
!      if ( r2.gt.1.d0) go to 100
!      tx = teta * (2.d0*va*vb) / r2
!      tz = teta * (va2 - vb2) / r2

!c 100  va  =2d0*rndm4()-1d0
!c      vb  =2d0*rndm4()-1d0
!c      va2 = va*va
!c      vb2 = vb*vb
!c      r2 = va2 + vb2
!c      if ( r2.gt.1.d0) go to 100
!c      f = SQRT(-2d0*LOG(r2)/r2)
!c      tx = teta*f*va
!c      tz = teta*f*vb


      tx= teta * RAN_GAUSS(1.0d0)
      tz= teta * RAN_GAUSS(1.0d0)

      tx = tx * 1000
      tz = tz * 1000

!---------- change p angle --------

      XP = XP + tx
      YP = YP + tz

      end if                !close the if(zlm.lt.DZ)

      xp=xp/1000
      yp=yp/1000


!-----------------------------------------------
! print out the energy loss and process experienced
!-------------------

               if (PROC(1:2).eq.'AM')then
                 PROC_dan=1
               elseif (PROC(1:2).eq.'VR') then
                 PROC_dan=2
               elseif (PROC(1:2).eq.'CH')then
                 PROC_dan=3
               elseif (PROC(1:2).eq.'VC') then
                 PROC_dan=4
               elseif (PROC(1:3).eq.'out')then
                 PROC_dan=-1
               elseif (PROC(1:8).eq.'absorbed') then
                 PROC_dan=5
               elseif (PROC(1:2).eq.'DC')then
                 PROC_dan=6
               elseif (PROC(1:3).eq.'pne')then
                 PROC_dan=7
               elseif (PROC(1:3).eq.'ppe')then
                 PROC_dan=8
               elseif (PROC(1:4).eq.'diff')then
                 PROC_dan=9
               elseif (PROC(1:4).eq.'ruth')then
                 PROC_dan=10
               elseif (PROC(1:4).eq.'TRVR')then
                 PROC_dan=100
               elseif (PROC(1:4).eq.'TRAM')then
                 PROC_dan=101
               endif

!      write(889,*) PC_in_dan, PC, PROC_dan, tx,tz,
!     & xp-xp_in-(kxmcs/1000), yp-yp_in-(kymcs/1000)


               write(*,*) "MOVE_AM done"

      RETURN
      END


!--------- Multiple Coulomb Scattering ---------

!c      xp=xp*1000
!c      yp=yp*1000
!      write(*,*)'xp initial:', xp, 'yp initial', yp
!. DEI - dE/dx stoping energy
!c      PC    = PC - DEI*DZ    ! energy lost beacause of ionization process[GeV]

!. DYA - rms of coloumb scattering
!c      DYA = (13.6/PC)*SQRT(DZ/DLr)             !MCS (mrad)
!      write(*,*)'dya=',dya

!c      kxmcs = DYA*RAN_GAUSS(1.)
!c      kymcs = DYA*RAN_GAUSS(1.)

!c      XP = xp+kxmcs
!c      yp = yp+kymcs

!      XP    = XP+DYA*RAN_GAUSS(1.)
!      YP    = YP+DYA*RAN_GAUSS(1.)

!c      xp = xp/1000
!c      yp = yp/1000

!--------END DANIELE------------------------------------


!.**************************************************************************
!     DANIELE subroutine for check if a nuclear interaction happen while in channeling
!.**************************************************************************
      SUBROUTINE MOVE_CH_(IS,NAM,DZ,X,XP,YP,PC,R,Rc)

      use mod_ranlux
      use mod_funlux

      IMPLICIT none
      integer IS,NAM
      double precision DZ,X,XP,YP,PC,R,Rc
      double precision DLAI(4),SAI(4),DES(4)
      double precision DLRI(4),DLYI(4),AI(4)
      double precision AM(30),QP(30),NPAI
      double precision Dc(4),eUm(4)
      double precision DYA,W_p
!      REAL RNDM4
!         REAL      RAN_GAUSS
      CHARACTER(LEN=50) PROC              !string that contains the physical process
      COMMON /ALAST/DLAI,SAI
      COMMON/CRYS/ DLRI,DLYI,AI,DES
      common /Proc2/PROC
      COMMON/eUc/eUm

!------- adding variables --------


      double precision cprob_cry(0:5,1:4)
      double precision cs(0:5,1:4)
      double precision csref_cry(0:5,1:4)
      double precision freep(4)
      double precision anuc_cry(4)
      double precision collnt_cry(4)
      double precision bn(4)
      double precision bnref_cry(4)

      double precision freeco_cry,ppel,ppsd,pptot,pptref_cry,pptco_cry
      double precision pperef_cry,sdcoe_cry,pref_cry,ppeco_cry

      double precision ecmsq,xln15s,bpp,xm2,bsd,t,teta,va,vb,va2,vb2
      double precision tx,tz,r2,zlm,f

      integer i,ichoix
      double precision aran

      common/scat_cry/cprob_cry,csref_cry,collnt_cry,bnref_cry,anuc_cry
      common/scat_cry/pptref_cry,pperef_cry,sdcoe_cry,pref_cry
      common/scat_cry/pptco_cry,ppeco_cry,freeco_cry

      double precision tlcut_cry,hcut_cry(4)
      double precision cgen_cry(200,4),tlow,thigh,ruth_cry
      external ruth_cry

      common/ruth_scat_cry/tlcut_cry,hcut_cry
      common/ruth_scat_cry/cgen_cry,mcurr_cry

      integer length_cry,mcurr_cry
      double precision xran_cry(1)



!-------daniele--------------
! adding new variables for calculation of average density seen
!---------------------------

      integer Np
      double precision aTF,dP,N_am,rho_min,rho_Max,u1,avrrho
      double precision x_i,Ueff,pv,Et,Ec,x_min,x_Max,nuc_cl_l
      double precision xminU,Umin

      common/dech/aTF,dP,u1

      double precision rho(4),z(4),Ime(4)
      double precision k,re,me,mp
      double precision anuc_cry2(4)
      double precision emr_cry(4)

      common/ion/rho,z,Ime,k,re,me,mp,anuc_cry2,emr_cry

      double precision csref_tot_rsc,csref_inel_rsc

!-------daniele--------------
! adding new variables to study the energy loss when the routine is called
!---------------------------


      double precision PC_in_dan     !daniele
      integer PROC_dan    !daniele
      double precision xp_in,yp_in,kxmcs,kymcs

      write(*,*) "Entering MOVE_CH"


      PC_in_dan=PC                   !daniele

      xp_in=XP
      yp_in=YP


!---------------daniele------------
! NEW TREATMENT OF SCATTERING ROUTINE BASED ON STANDARD SIXTRACK ROUTINE
!----------------------------------


!------- useful calculations for cross-section and event topology calculation --------------

      write(*,*) "ecmsq"
      ecmsq = 2 * 0.93828d0 * PC
      write(*,*) "xln15s"
      xln15s=log(0.15*ecmsq)
! pp(pn) data
!      pptot = pptref_cry *(PC / pref_cry)** pptco_cry
!      ppel = pperef_cry *(PC / pref_cry)** ppeco_cry
!      ppsd = sdcoe_cry * log(0.15d0 * ecmsq)
!      bpp = 8.5d0 + 1.086d0 * log(sqrt(ecmsq))

!new models, see Claudia's thesis

      write(*,*) "models for cross sections"
      pptot=0.041084d0-0.0023302d0*log(ecmsq)+0.00031514d0*log(ecmsq)**2
      ppel=(11.7d0-1.59d0*log(ecmsq)+0.134d0*log(ecmsq)**2)/1000
      ppsd=(4.3d0+0.3d0*log(ecmsq))/1000
      bpp=7.156d0+1.439d0*log(sqrt(ecmsq))

!------ distribution for Ruth. scatt.---------

      write(*,*) "tlow"
      tlow=tlcut_cry
      write(*,*) "mcurr"
      mcurr_cry=IS
      write(*,*) "thigh"
      thigh=hcut_cry(IS)
!      write(*,*)tlow,thigh
      write(*,*) "funlxp"
      call funlxp(ruth_cry,cgen_cry(1,IS),tlow,thigh)

!---------- rescale the total and inelastic cross-section accordigly to the average density seen

      write(*,*) "rescaling calculation"
      x_i=X

      Np=INT(x_i/dP)       !calculate in which cristalline plane the particle enters
      x_i=x_i-Np*dP        !rescale the incoming x at the left crystalline plane
      x_i=x_i-(dP/2.d0)    !rescale the incoming x in the middle of crystalline planes

      pv=PC*1.d9*PC*1.d9/sqrt(PC*1.d9*PC*1.d9+93828.d6*93828.d6)        !calculate pv=P/E

      Ueff=eUm(IS)*(2.d0*x_i/dP)*(2.d0*x_i/dP)+pv*x_i/R   !calculate effective potential

      Et=pv*XP*XP/2.+Ueff                                !calculate transverse energy

      Ec=eUm(IS)*(1.d0-Rc/R)*(1.d0-Rc/R)                 !calculate critical energy in bent crystals

!---- to avoid negative Et

      xminU=-dP*dP*PC*1e9/(8.d0*eUm(IS)*R)
      Umin=abs(eUm(IS)*(2.d0*xminU/dP)*(2.d0*xminU/dP)+pv*xminU/R)
      write(*,*) "Umin", Umin
      Et=Et+Umin
      Ec=Ec+Umin

!-----

!      write(*,*)xminU,Umin,Et,Ec,pv,XP,Ueff

      x_min=-(dP/2.d0)*Rc/R-(dP/2.d0)*sqrt(Et/Ec);          !calculate min e max of the trajectory between crystalline planes
      x_Max=-(dP/2.d0)*Rc/R+(dP/2.d0)*sqrt(Et/Ec);

      x_min=x_min-dP/2.d0;                                    !change ref. frame and go back with 0 on the crystalline plane on the left
      x_Max=x_Max-dP/2.d0;
      write(*,*) "x_min", x_min
      write(*,*) "x_Max", x_Max

      N_am=rho(IS)*6.022d23*1.d6/anuc_cry2(IS)           !calculate the "normal density" in m^-3

      write(*,*) N_am, dP, u1
      rho_Max=N_am*dP/2.d0*(erf(x_Max/sqrt(2.d0*u1*u1)) &      !calculate atomic density at min and max of the trajectory oscillation
     & -erf((dP-x_Max)/sqrt(2*u1*u1)))

      rho_min=N_am*dP/2.d0*(erf(x_min/sqrt(2.d0*u1*u1)) &
     & -erf((dP-x_min)/sqrt(2*u1*u1)));

     write(*,*) "rho_Max", rho_Max
     write(*,*) "rho_min", rho_min
      avrrho=(rho_Max-rho_min)/(x_Max-x_min)                        !"zero-approximation" of average nuclear density seen along the trajectory
      write(*,*) "avrrho", avrrho

      write(*,*) "N_am", N_am
      avrrho=2.d0*avrrho/N_am
      write(*,*) "avrrho", avrrho

      write(*,*) "csref_cry", csref_cry(0,IS)
      csref_tot_rsc=csref_cry(0,IS)*avrrho        !rescaled total ref cs
      write(*,*) "csref_tot_rsc", csref_tot_rsc
      csref_inel_rsc=csref_cry(1,IS)*avrrho        !rescaled inelastic ref cs

!      write(889,*) x_i,pv,Ueff,Et,Ec,N_am,avrrho,
!     + csref_tot_rsc,csref_inel_rsc

!---------- cross-section calculation -----------

      write(*,*) "cross sections calculation"
!
! freep: number of nucleons involved in single scattering
      write(*,*) "freep"
        freep(IS) = freeco_cry * anuc_cry(IS)**(1d0/3d0)
! compute pp and pn el+single diff contributions to cross-section
! (both added : quasi-elastic or qel later)
        write(*,*) "cs(3)"
        cs(3,IS) = freep(IS) * ppel
        write(*,*) "cs(4)"
        cs(4,IS) = freep(IS) * ppsd
!
! correct TOT-CSec for energy dependence of qel
! TOT CS is here without a Coulomb contribution
        write(*,*) "cs(0)"
        cs(0,IS) = csref_tot_rsc + freep(IS) * (pptot - pptref_cry)
        write(*,*) "bnref_cry", bnref_cry(IS)
        write(*,*) "cs(0)", cs(0,IS)
        write(*,*) "csref_tot_rsc", csref_tot_rsc
!        write(*,*) "bn(IS)"
!        bn(IS) = bnref_cry(IS) * cs(0,IS) / csref_tot_rsc
! also correct inel-CS
        write(*,*) "cs(1)"
        if (csref_tot_rsc == 0) then
          cs(1,IS) = 0
        else
          cs(1,IS) = csref_inel_rsc * cs(0,IS) / csref_tot_rsc
        end if
!
! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
        write(*,*) "cs(2)"
        cs(2,IS) = cs(0,IS) - cs(1,IS) - cs(3,IS) - cs(4,IS)
        write(*,*) "cs(5)"
        cs(5,IS) = csref_cry(5,IS)
! Now add Coulomb
        write(*,*) "cs(0)"
        cs(0,IS) = cs(0,IS) + cs(5,IS)
! Interaction length in meter
!        xintl(ma) = 0.01d0*anuc(ma)/(fnavo * rho(ma)*cs(0,ma)*1d-24)  !don't need at the moment, take it from pdg


! Calculate cumulative probability
        write(*,*) "cumulative probability"

        if (cs(0,IS) == 0) then
          do i=1,4
            cprob_cry(i,IS)=cprob_cry(i-1,IS)
          end do
        else
          do i=1,4
            cprob_cry(i,IS)=cprob_cry(i-1,IS)+cs(i,IS)/cs(0,IS)
          end do
        end if


!--------- Multiple Coulomb Scattering ---------

      xp=xp*1000
      yp=yp*1000

!-------- not needed, energy loss by ionization taken into account in the main body


!      write(*,*)'xp initial:', xp, 'yp initial', yp
!. DEI - dE/dx stoping energy
!      PC    = PC - DEI*DZ    ! energy lost beacause of ionization process[GeV]

!. DYA - rms of coloumb scattering
!      DYA = (13.6/PC)*SQRT(DZ/DLr)             !MCS (mrad)
!      write(*,*)'dya=',dya

!      kxmcs = DYA*RAN_GAUSS(1.)
!      kymcs = DYA*RAN_GAUSS(1.)

!      XP = xp+kxmcs
!      yp = yp+kymcs

!c     XP    = XP+DYA*RAN_GAUSS(1.)
!c      YP    = YP+DYA*RAN_GAUSS(1.)

      if(NAM .eq. 0) return                    !turn on/off nuclear interactions


!--------- Can nuclear interaction happen? -----

      if (avrrho == 0) then
        nuc_cl_l=1.0d6        !rescaled nuclear collision length
      else
        nuc_cl_l=collnt_cry(IS)/avrrho        !rescaled nuclear collision length
      end if

      zlm=-nuc_cl_l*log(dble(rndm4()))

!      zlm=0.0

      write(889,*) x_i,pv,Ueff,Et,Ec,N_am,avrrho,csref_tot_rsc,csref_inel_rsc,nuc_cl_l

      if(zlm.lt.DZ) then

!--------- Choose nuclear interaction --------


      aran=dble(rndm4())
      i=1
  10  if ( aran.gt.cprob_cry(i,IS) ) then
          i=i+1
          goto 10
      endif
      ichoix=i


!-------- Do the interaction ----------

!      ichoix=3

      if (ichoix.eq.1) then

        PROC = 'ch_absorbed'            !deep inelastic, impinging p disappeared

      endif

      if (ichoix.eq.2) then          ! p-n elastic
        PROC = 'ch_pne'
        write(*,*) "bn(IS)"
        bn(IS) = bnref_cry(IS) * cs(0,IS) / csref_tot_rsc
        t = -log(dble(rndm4()))/bn(IS)
      endif

      if ( ichoix .eq. 3 ) then   !p-p elastic
        PROC='ch_ppe'
        t = -log(dble(rndm4()))/bpp
      endif

      if ( ichoix .eq. 4 ) then   !single diffractive
        PROC='ch_diff'
        xm2 = exp( dble(rndm4()) * xln15s )
        PC = PC  *(1.d0 - xm2/ecmsq)
          if ( xm2 .lt. 2.d0 ) then
             bsd = 2.d0 * bpp
           elseif (( xm2 .ge. 2.d0 ).and. ( xm2 .le. 5.d0 )) then
                bsd = (106.d0-17.d0*xm2) *  bpp / 36.d0
           elseif ( xm2 .gt. 5.d0 ) then
                bsd = 7.d0 * bpp / 12.d0
          endif
        t = -log(dble(rndm4()))/bsd
      endif

      if ( ichoix .eq. 5 ) then
        PROC='ch_ruth'
           length_cry=1
           call funlux( cgen_cry(1,IS) ,xran_cry,length_cry)
           t=xran_cry(1)
      endif

!---------- calculate the related kick -----------


      if ( ichoix .eq. 4) then
        teta = sqrt(t)/PC_in_dan                !DIFF has changed PC!!!
      else
        teta = sqrt(t)/PC
      endif


! 100  va  =2d0*rndm4()-1d0
!      vb = dble(rndm4())
!      va2 = va*va
!      vb2 = vb*vb
!      r2 = va2 + vb2
!      if ( r2.gt.1.d0) go to 100
!      tx = teta * (2.d0*va*vb) / r2
!      tz = teta * (va2 - vb2) / r2

!c 100  va  =2d0*rndm4()-1d0
!c      vb  =2d0*rndm4()-1d0
!c      va2 = va*va
!c      vb2 = vb*vb
!c      r2 = va2 + vb2
!c      if ( r2.gt.1.d0) go to 100
!c      f = SQRT(-2d0*LOG(r2)/r2)
!c      tx = teta*f*va
!c      tz = teta*f*vb


      tx= teta * RAN_GAUSS(1.0d0)
      tz= teta * RAN_GAUSS(1.0d0)

      tx = tx * 1000
      tz = tz * 1000

!---------- change p angle --------

      XP = XP + tx
      YP = YP + tz

      end if                !close the if(zlm.lt.DZ)

      xp=xp/1000
      yp=yp/1000


!-----------------------------------------------
! print out the energy loss and process experienced
!-------------------

               if (PROC(1:2).eq.'AM')then
                 PROC_dan=1
               elseif (PROC(1:2).eq.'VR') then
                 PROC_dan=2
               elseif (PROC(1:2).eq.'CH')then
                 PROC_dan=3
               elseif (PROC(1:2).eq.'VC') then
                 PROC_dan=4
               elseif (PROC(1:3).eq.'out')then
                 PROC_dan=-1
               elseif (PROC(1:11).eq.'ch_absorbed') then
                 PROC_dan=15
               elseif (PROC(1:2).eq.'DC')then
                 PROC_dan=6
               elseif (PROC(1:6).eq.'ch_pne')then
                 PROC_dan=17
               elseif (PROC(1:6).eq.'ch_ppe')then
                 PROC_dan=18
               elseif (PROC(1:7).eq.'ch_diff')then
                 PROC_dan=19
               elseif (PROC(1:7).eq.'ch_ruth')then
                 PROC_dan=20

               endif

!      write(889,*) PC_in_dan, PC, PROC_dan, tx,tz,
!     & xp-xp_in-(kxmcs/1000), yp-yp_in-(kymcs/1000)

!-------- Block Data ----------

               write(*,*) "MOVE_CH done"

      RETURN
      END
!
      BLOCK DATA
      implicit none
      double precision DLAI(4),SAI(4)
      double precision DLRI(4),DLYI(4),AI(4),DES(4)
      double precision eUm(4)
      COMMON /ALAST/DLAI,SAI
      COMMON/CRYS/ DLRI,DLYI,AI,DES
      COMMON/eUc/  eUm
!-----4 substances: Si(110),W(110),C,Ge----------------------------
      DATA DLRI/0.0937,.0035,0.188,.023/         &! radiation  length(m), updated from pdg for Si
     &    ,DLYI/.455, .096, .400, .162/          &! nuclear length(m)
     &    ,AI /0.96E-7, 0.56E-7, 0.63E-7, 1.E-7/ & !Si110 1/2 interplan. dist. mm
     &    ,DES/0.56,  3.0,  0.6, 1./             & ! energy deposition in subst(GeV/m)
     &    ,DLAI/1.6,  0.57, 2.2, 1.0/            & ! elastic length(m)
     &    ,SAI /42.,  140.,  42., 50./           & ! elastic scat. r.m.s(mr)
     &    ,eUm/21.34,  21.,   21.,   21./         ! only for Si(110) potent. [eV]


!------------- daniele, more data needed---------

      double precision cprob_cry(0:5,1:4)
!      double precision cs(5,4)
      double precision csref_cry(0:5,1:4)
      double precision anuc_cry(4)
      double precision collnt_cry(4)
      double precision bnref_cry(4)
      double precision pptref_cry,pperef_cry,sdcoe_cry,pref_cry
      double precision pptco_cry,ppeco_cry,freeco_cry


      common/scat_cry/cprob_cry,csref_cry,collnt_cry,bnref_cry,anuc_cry
      common/scat_cry/pptref_cry,pperef_cry,sdcoe_cry,pref_cry
      common/scat_cry/pptco_cry,ppeco_cry,freeco_cry


      integer i

      data (cprob_cry(0,i),i=1,4)/4*0.0d0/
      data (cprob_cry(5,i),i=1,4)/4*1.0d0/

! pp cross-sections and parameters for energy dependence
      data pptref_cry,pperef_cry/0.04d0,0.007d0/
      data sdcoe_cry,pref_cry/0.00068d0,450.0d0/
      data pptco_cry,ppeco_cry,freeco_cry/0.05788d0,0.04792d0,1.618d0/

! Atomic mass [g/mole] from pdg

      data (anuc_cry(i),i=1,4)/28.08d0,0.0d0,0.0d0,0.0d0/     !implemented only Si

! Total and nuclear cross-sections [barn], implemeted only Si
      data csref_cry(0,1),csref_cry(1,1)/0.664d0, 0.430d0/  !from pdg
!      data csref_cry(0,1),csref_cry(1,1)/0.762d0, 0.504d0/  !with glauber's approx (NIMB 268 (2010) 2655-2659)
      data csref_cry(0,2),csref_cry(1,2)/0.0d0, 0.0d0/
      data csref_cry(0,3),csref_cry(1,3)/0.0d0, 0.0d0/
      data csref_cry(0,4),csref_cry(1,4)/0.0d0, 0.0d0/

      data csref_cry(5,1)/0.039d-2/
      data csref_cry(5,2)/0.0d0/
      data csref_cry(5,3)/0.0d0/
      data csref_cry(5,4)/0.0d0/



! Nuclear Collision length [m] from pdg (only for Si for the moment)

      data (collnt_cry(i),i=1,4)/0.3016d0,0.0d0,0.0d0,0.0d0/

! Nuclear elastic slope from Schiz et al.,PRD 21(3010)1980

      data (bnref_cry(i),i=1,4)/123.2d0,0.0d0,0.0d0,0.0d0/      !(only for Si for the moment)

! For calculation of dE/dX due to ionization

      double precision rho(4),z(4),Ime(4),anuc_cry2(4) !Daniele: material parameters for dE/dX calculation
      double precision k,re,me,mp !Daniele: parameters for dE/dX calculation (const,electron radius,el. mass, prot.mass)

      double precision emr_cry(4) !nuclear radius for Ruth scatt.

      common/ion/rho,z,Ime,k,re,me,mp,anuc_cry2,emr_cry

      data (rho(i),i=1,4)/2.33d0,0.0d0,0.0d0,0.0d0/  !material density [g/cm^3]
      data (z(i),i=1,4)/14.0d0,0.0d0,0.0d0,0.0d0/    !atomic number
      data (Ime(i),i=1,4)/173.0d-6,0.0d0,0.0d0,0.0d0/!mean exitation energy [MeV]
      data (anuc_cry2(i),i=1,4)/28.08d0,0.0d0,0.0d0,0.0d0/     !implemented only Si
      data (emr_cry(i),i=1,4)/0.306d0,0.0d0,0.0d0,0.0d0/      !nuclear radius take from R_Be*(A_Si/A_Be)^1/3, where R_Be=0.302


      data k/0.307075/ !constant in front bethe-bloch [MeV g^-1 cm^2]
      data re/2.818d-15/  !electron radius [m]
      data me/0.510998910/ !electron mass [MeV/c^2]
      data mp/938.272013/ !proton mass [MeV/c^2]

! For calculation of dechanneling length

      double precision aTF,dP,u1

      common/dech/aTF,dP,u1

      data aTF/0.194d-10/ !screening function [m]
      data dP/1.92d-10/   !distance between planes (110) [m]
      data u1/0.075d-10/  !thermal vibrations amplitude

      double precision tlcut_cry,hcut_cry(4)   !param. for Ruth scatt.
      double precision cgen_cry(200,4)
      integer mcurr_cry

      common/ruth_scat_cry/tlcut_cry,hcut_cry
      common/ruth_scat_cry/cgen_cry,mcurr_cry

      data tlcut_cry/0.0009982d0/
      data (hcut_cry(i),i=1,4)/0.02d0,0.0d0,0.0d0,0.0d0/


      END

!------------------daniele: definition of rutherford scattering formula--------!

      function ruth_cry(t_cry)
      implicit none
      integer nmat_cry,mcurr_cry
      parameter(nmat_cry=4)
!      double precision z(4),emr_cry(4),tlcut_cry,hcut_cry(4)
      double precision tlcut_cry,hcut_cry(4)

      double precision rho(4),z(4),Ime(4),anuc_cry2(4) !Daniele: material parameters for dE/dX calculation
      double precision k,re,me,mp !Daniele: parameters for dE/dX calculation (const,electron radius,el. mass, prot.mass)

      double precision emr_cry(4)

      double precision cgen_cry(200,4)

      common/ion/rho,z,Ime,k,re,me,mp,anuc_cry2,emr_cry

!      common/ion/z,emr_cry
      common/ruth_scat_cry/tlcut_cry,hcut_cry
      common/ruth_scat_cry/cgen_cry,mcurr_cry

      double precision ruth_cry,t_cry
      double precision cnorm,cnform
      parameter(cnorm=2.607d-4,cnform=0.8561d3)
!      parameter(cnorm=2.607d-5,cnform=0.8561d3) !daniele: corrected constant
!c      write(6,'('' t,exp'',2e15.8)')t,t*cnform*EMr(mcurr)**2
!      write(*,*)mcurr_cry,emr_cry(mcurr_cry),z(mcurr_cry),t_cry

       ruth_cry=cnorm*exp(-t_cry*cnform*emr_cry(mcurr_cry)**2)*(z(mcurr_cry)/t_cry)**2

      end



!~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STUFF BELOW IS NOT IN USE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~



!------------------------------past treatment commented (c always comment, cc old line)----------------------------------------------------------------------------------



!c      xp=xp*1000
!c      yp=yp*1000
!      write(*,*)'xp initial:', xp, 'yp initial', yp
!. DEI - dE/dx stoping energy
!c      PC    = PC - DEI*DZ    ! energy lost beacause of ionization process[GeV]
!. Coulomb scattering
!. DYA - rms of coloumb scattering
!c      DYA = (14.0/PC)*SQRT(DZ/DLr)             !MCS (mrad)
!      write(*,*)'dya=',dya
!c      XP    = XP+DYA*RAN_GAUSS(1.)
!c      YP    = YP+DYA*RAN_GAUSS(1.)

!c      if(NAM .eq. 0) return
!.  Elastic scattering
!c      IF (rndm4() .LE. DZ/DLAI(IS)) THEN
!c        PROC = 'mcs'
!         write(*,*)'case 1'
!        XP    = XP+SAI(IS)*RAN_GAUSS(1.)/PC
!        YP    = YP+SAI(IS)*RAN_GAUSS(1.)/PC
!c        xp    = xp+196.*RAN_GAUSS(1.)/PC ! angle elast. scattering in R plane
!c        yp    = yp+196.*RAN_GAUSS(1.)/PC
!c      ENDIF
!.  Diffraction interaction
!c      IF (rndm4() .LE. DZ/(DLY*6.143)) THEN
!c        PROC = 'diff'
!         write(*,*)'case 2'

!######################################################
! daniele: comment old treatement to implement the same as standard SixTrack
!######################################################

!        XP    = XP+ 257.0*RAN_GAUSS(1.)/PC ! angle elast. scattering in R plane
!        YP    = YP+ 257.0*RAN_GAUSS(1.)/PC
!        W_p = rndm4()
!        PC = PC -0.5*(0.3*PC)**W_p             ! m0*c = 1 GeV/c for particle

!############## end comment -> "new" treatment of diffractive interactions #########

!         ecmsq = 2 * 0.93828d0 * PC
!         xln15s=log(0.15*ecmsq)
!         bpp = 8.5d0 + 1.086d0 * log(sqrt(ecmsq))
!         xm2 = exp( dble(rndm4()) * xln15s )
!         PC = PC  *(1.d0 - xm2/ecmsq)

!         if ( xm2 .lt. 2.d0 ) then
!              bsd = 2.d0 * bpp
!            elseif (( xm2 .ge. 2.d0 ).and. ( xm2 .le. 5.d0 )) then
!              bsd = (106.d0-17.d0*xm2) *  bpp / 26.d0
!            elseif ( xm2 .gt. 5.d0 ) then
!              bsd = 7.d0 * bpp / 12.d0
!         endif

!         t = -log(dble(rndm4()))/bsd
!         teta = sqrt(t)/PC
!      10 va  =2d0*rndm4()-1d0
!         vb = dble(rndm4())
!         va2 = va*va
!         vb2 = vb*vb
!         r2 = va2 + vb2
!         if ( r2.gt.1.d0) go to 10
!         tx = teta * (2.d0*va*vb) / r2
!         tz = teta * (va2 - vb2) / r2

!         XP = XP + tx
!         YP = YP + tz

!#################### end of "new" treatment ###############

!c      ENDIF
!.  Inelastic interaction
!c      IF (rndm4() .LE. DZ/DLY) THEN
!        PC = 0.                                         !daniele: needed for debugged energy loss in the crystal, otherwise SixTrack get stuck if PC=0 (anyway the particle is "forgot" from now)
!c        PROC = 'absorbed'
!c      ENDIF
!      write(*,*)'xp final:', xp, 'yp final', yp
!c      xp=xp/1000
!c      yp=yp/1000


!---------------DANIELE--------------------------------
! print out the energy loss and process experienced
!-------------------

!c               if (PROC(1:2).eq.'AM')then
!c                 PROC_dan=1
!c               elseif (PROC(1:2).eq.'VR') then
!c                 PROC_dan=2
!c               elseif (PROC(1:2).eq.'CH')then
!c                 PROC_dan=3
!c               elseif (PROC(1:2).eq.'VC') then
!c                 PROC_dan=3
!c               elseif (PROC(1:3).eq.'out')then
!c                 PROC_dan=-1
!c              elseif (PROC(1:8).eq.'absorbed') then
!c                 PROC_dan=5
!c               elseif (PROC(1:2).eq.'DC')then
!c                 PROC_dan=6
!c               elseif (PROC(1:3).eq.'mcs')then
!c                 PROC_dan=7
!c               elseif (PROC(1:4).eq.'diff')then
!c                 PROC_dan=8
!c               endif

!c      write(889,*) PC_in_dan, PC, PROC_dan

!--------END DANIELE------------------------------------

!c      RETURN
!c      END
!
!c      BLOCK DATA
!c      implicit none
!c      double precision DLAI(4),SAI(4)
!c      double precision DLRI(4),DLYI(4),AI(4),DES(4)
!c      double precision eUm(4)
!c      COMMON /ALAST/DLAI,SAI
!c      COMMON/CRYS/ DLRI,DLYI,AI,DES
!c      COMMON/eUc/  eUm
!-----4 substances: Si(110),W(110),C,Ge----------------------------
!c      DATA DLRI/0.09336,.0035,0.188,.023/    ! radiation  length(m)
!c     +    ,DLYI/.455, .096, .400, .162/      ! nuclear length(m)
!c     +    ,AI /0.96E-7, 0.56E-7, 0.63E-7, 1.E-7/  !Si110 1/2 interplan. dist. mm
!c     +    ,DES/0.56,  3.0,  0.6, 1./         ! energy deposition in subst(GeV/m)
!c     +    ,DLAI/1.6,  0.57, 2.2, 1.0/        ! elastic length(m)
!c     +    ,SAI /42.,  140.,  42., 50./       ! elastic scat. r.m.s(mr)
!c     +    ,eUm/21.34,  21.,   21.,   21./    ! only for Si(110) potent. [eV]
!c      END







!                        'MATH.F'   10.07.01
!1  RANNOR - Gauss distribution simulation procedure.
!2  RF RNDM- ­  [0,1]        *
!3  RNDMST - .
!4  ...
!.**************************************************************************
!     subroutine for the generation of random numbers with a gaussian
!     distribution
!.**************************************************************************
!.**************************************************************************
!       SUBROUTINE RANNOR(A,B) !gaussian with mean 0 and sigma 1
! !
! !1    Gauss distribution simulation procedure
! !
!       Y = RNDM(1.)
!       IF (Y .EQ. 0.) Y = RNDM(1.)
!       Z = RNDM(1.)
!       X = 6.2831853*Z
!       A1= SQRT(-2.0*ALOG(Y))
!       A = A1*SIN(X)
!       B = A1*COS(X)
! !
!       RETURN
!       END
! !
!       REAL function RNDM(rdummy)
!       REAL          u,c,cd,cm, rrdummy
!       COMMON/RANDOM/  u(97),c,cd,cm,i,j
! !
!       rrdummy = rdummy
! !
!       RNDM = u(i)-u(j)
!       if (RNDM .LT. 0.) RNDM = RNDM+1.
!       U(i) = RNDM
!       i = i-1
!       if ( i .EQ. 0 ) i = 97
!       j = j-1
!       if ( j .EQ. 0 ) j = 97
!       c = c-cd
!       if ( C .LT. 0.) c = c+cm
!       RNDM = RNDM - c
!       if ( RNDM .LT. 0. ) RNDM = RNDM+1.
! !
!       return
!       END
! ! !.**************************************************************************
! !
! !.**************************************************************************
!       subroutine RNDMST(NA1,NA2,NA3,NB1)
!       REAL u,c,cd,cm
!       COMMON/RANDOM/ u(97),c,cd,cm,i,j
! !
!       ma1 = na1
!       ma2 = na2
!       ma3 = na3
!       mb1 = nb1
!       i   = 97
!       j   = 33
! !
!       do ii2 = 1,97
!         s = 0.0
!         t = 0.5
!         do ii1 = 1,24
!           mat  = MOD(MOD(ma1*ma2,179)*ma3,179)
!           ma1  = ma2
!           ma2  = ma3
!           ma3  = mat
!           mb1  = MOD(53*mb1+1,169)
!           if (MOD(MB1*MAT,64) .GE. 32 ) s = s+t
!           t = 0.5*t
!         end do
!         u(ii2) = s
!       end do
! !
!       c  =   362436./16777216.
!       cd =  7654321./16777216.
!       cm = 16777213./16777216.
! !
!       return
!       END
! !.**************************************************************************
! !
! !.**************************************************************************
!       subroutine RNDMIN(uin,cin,cdin,cmin,iin,jin)
!       REAL uin(97),cin,cdin,cmin
!       REAL u,c,cd,cm
!       COMMON/RANDOM/ u(97),c,cd,cm,i,j
! !
!       do kkk = 1,97
!         u(kkk) = uin(kkk)
!       end do
! !
!       c  = cin
!       cd = cdin
!       cm = cmin
!       i  = iin
!       j  = jin
! !
!       return
!       END
! !.**************************************************************************
! !
! !.**************************************************************************
!       subroutine RNDMOU(UOUT,COUT,CDOUT,CMOUT,IOUT,JOUT)
!       REAL uout(97),cout,cdout,cmout
!       REAL u,c,cd,cm
!       COMMON/RANDOM/ u(97),c,cd,cm,i,j
! !
!       do kkk = 1,97
!         uout(kkk) = u(kkk)
!       end do
! !
!       COUT  = C
!       CDOUT = CD
!       CMOUT = CM
!       IOUT  = I
!       JOUT  = J
! !
!       return
!       END
! !.**************************************************************************
! !
! !.**************************************************************************
!       subroutine RNDMTE(IO)
! !
! !.    *******************************************************************
! !.    *  SUBROUTINE RNDMTE(IO)
! !.    *******************************************************************
! !
!       REAL uu(97)
!       REAL u(6),x(6),d(6)
!       DATA u / 6533892.0 , 14220222.0 ,  7275067.0 ,  6172232.0 ,  8354498.0 , 10633180.0 /
! !
!       call RNDMOU(UU,CC,CCD,CCM,II,JJ)
!       call RNDMST(12,34,56,78)
! !
!       do ii1 = 1,20000
!         xx = rndm4()
!       end do
! !
!       sd = 0.0
!       do ii2 = 1,6
!         x(II2)  = 4096.*(4096.*rndm4())
!         d(II2)  = x(II2)-u(II2)
!         sd = sd+d(II2)
!       end do
! !
!       call RNDMIN(uu,cc,ccd,ccm,ii,jj)
!       if (IO.EQ.1 .OR. SD .NE. 0.) write(6,10) (u(I),x(I),d(I),i=1,6)
! !
!  10   format('  === TEST OF THE RANDOM-GENERATOR ===',/, &
!      &       '    EXPECTED VALUE    CALCULATED VALUE     DIFFERENCE',/, &
!      &       6(F17.1,F20.1,F15.3,/), &
!      &       '  === END OF TEST ;', &
!      &       '  GENERATOR HAS THE SAME STATUS AS BEFORE CALLING RNDMTE')
!       return
!       END


