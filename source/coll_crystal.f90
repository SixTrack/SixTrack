! ============================================================================ !
!
!  Crystal Collimation Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ============================================================================ !

module coll_crystal

  use floatPrecision
  use numerical_constants
  use coll_materials, only : nmat

  implicit none

  integer, parameter :: max_ncoll = 99

  logical,              save :: bool_create
  integer, allocatable, save :: bool_proc(:)
  integer, allocatable, save :: bool_proc_old(:)

  integer,           save :: iProc
  logical,           save :: changed_tilt1(max_ncoll)
  logical,           save :: changed_tilt2(max_ncoll)

  real(kind=fPrec),  save :: Rcurv    ! Crystal geometrical parameters [m]
  real(kind=fPrec),  save :: C_xmax   ! Crystal geometrical parameters [m]
  real(kind=fPrec),  save :: C_ymax   ! Crystal geometrical parameters [m]
  real(kind=fPrec),  save :: Alayer   ! Crystal amorphous layer [mm]
  real(kind=fPrec),  save :: miscut   ! Crystal miscut angle in rad
  integer,           save :: C_orient ! Crystal orientation [0-1]

  ! Rutherford Scatter
  real(kind=fPrec), parameter :: tlcut_cry   = 0.0009982d0
  real(kind=fPrec), save      :: cgen_cry(200,nmat)
  integer,          save      :: mcurr_cry

  real(kind=fPrec), save :: enr,mom,betar,gammar,bgr !Daniele: energy,momentum,beta relativistic, gamma relativistic
  real(kind=fPrec), save :: Tmax,plen !Daniele: maximum energy tranfer in single collision, plasma energy (see pdg)

  real(kind=fPrec), parameter :: re  = 2.818d-15   ! electron radius [m]
  real(kind=fPrec), parameter :: me  = 0.510998910 ! electron mass [MeV/c^2]
  real(kind=fPrec), parameter :: mp  = 938.272013  ! proton mass [MeV/c^2]

  real(kind=fPrec), parameter :: aTF = 0.194d-10   ! Screening function [m]
  real(kind=fPrec), parameter :: dP  = 1.92d-10    ! Distance between planes (110) [m]
  real(kind=fPrec), parameter :: u1  = 0.075d-10   ! Thermal vibrations amplitude

  ! pp cross-sections and parameters for energy dependence
  real(kind=fPrec), parameter :: pptref_cry = 0.04d0
  real(kind=fPrec), parameter :: pperef_cry = 0.007d0
  real(kind=fPrec), parameter :: sdcoe_cry  = 0.00068d0
  real(kind=fPrec), parameter :: pref_cry   = 450.0d0
  real(kind=fPrec), parameter :: pptco_cry  = 0.05788d0
  real(kind=fPrec), parameter :: ppeco_cry  = 0.04792d0
  real(kind=fPrec), parameter :: freeco_cry = 1.618d0

  ! Crystal Specific Material Arrays
  real(kind=fPrec), save :: dlri(nmat)   = zero
  real(kind=fPrec), save :: dlyi(nmat)   = zero
  real(kind=fPrec), save :: ai(nmat)     = zero
  real(kind=fPrec), save :: eUm(nmat)    = zero
  real(kind=fPrec), save :: collnt(nmat) = zero ! Nuclear Collision length [m]

  ! Processes
  integer, parameter :: proc_out         =  -1
  integer, parameter :: proc_AM          =   1
  integer, parameter :: proc_VR          =   2
  integer, parameter :: proc_CH          =   3
  integer, parameter :: proc_VC          =   4
  integer, parameter :: proc_absorbed    =   5
  integer, parameter :: proc_DC          =   6
  integer, parameter :: proc_pne         =   7
  integer, parameter :: proc_ppe         =   8
  integer, parameter :: proc_diff        =   9
  integer, parameter :: proc_ruth        =  10
  integer, parameter :: proc_ch_absorbed =  15
  integer, parameter :: proc_ch_pne      =  17
  integer, parameter :: proc_ch_ppe      =  18
  integer, parameter :: proc_ch_diff     =  19
  integer, parameter :: proc_ch_ruth     =  20
  integer, parameter :: proc_TRVR        = 100
  integer, parameter :: proc_TRAM        = 101

contains

subroutine cry_expandArrays(npart_new)

  use mod_alloc
  use numerical_constants

  integer, intent(in) :: npart_new

  call alloc(bool_proc,     npart_new, 0, "bool_proc")
  call alloc(bool_proc_old, npart_new, 0, "bool_proc_old")

end subroutine cry_expandArrays

subroutine cry_init

  use coll_materials

  integer m

  ! dlri : Radiation length(m), updated from PDG for Si
  ! dlyi : Nuclear length(m)
  ! ai   : Si110 1/2 interplan. dist. mm, Ge taken from A. Fomin, Si from initial implementation
  ! eUm  : Only for Si(110) and Ge(110) potent. [eV], Ge taken from A. Fomin, Si from initial implementation

  ! Si
  m = collmat_getCollMatID("Si")
  dlri(m)   = 0.0937
  dlyi(m)   = 0.4652
  ai(m)     = 0.96e-7
  eUm(m)    = 21.34
  collnt(m) = 0.3016d0

  ! W
  m = collmat_getCollMatID("W")
  dlri(m)   = 0.0035
  dlyi(m)   = 0.096
  ai(m)     = 0.56e-7
  eUm(m)    = 21.0
  collnt(m) = 0.0d0

  ! C
  m = collmat_getCollMatID("C")
  dlri(m)   = 0.188
  dlyi(m)   = 0.400
  ai(m)     = 0.63e-7
  eUm(m)    = 21.0
  collnt(m) = 0.0d0

  ! Ge
  m = collmat_getCollMatID("Ge")
  dlri(m)   = 0.02302
  dlyi(m)   = 0.2686
  ai(m)     = 1.0e-7
  eUm(m)    = 40.0
  collnt(m) = 0.1632d0

end subroutine cry_init

end module coll_crystal

!===============================================================
! MDA: added icoll and MAX_NPART parameters
!
      SUBROUTINE collimate_cry(icoll, iturn, ie, C_LENGTH, &
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
      use coll_common, only : cry_proc, xp_pencil0, yp_pencil0, x_pencil, y_pencil, pencil_dx, ipencil
      use floatPrecision
      use coll_crystal

      IMPLICIT NONE
!

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
         real(kind=fPrec)            p0
         real(kind=fPrec)            zlm
         real(kind=fPrec)            x,xp
         real(kind=fPrec)            shift
         real(kind=fPrec)            x_shift, xp_shift,s_shift !coordinates after shift/rotation
         real(kind=fPrec)            x_rot, xp_rot,s_rot
         real(kind=fPrec)            x_temp, xp_temp,s_temp !!all the _temp variables are used when you hit the cry from below
         real(kind=fPrec)            tilt_int, x_int,xp_int,s_int       !all the _int variables are used when you hit the cry from below (int=interaction point)
         real(kind=fPrec)            x00
         real(kind=fPrec)            z
         real(kind=fPrec)            z00
         real(kind=fPrec)            zp
         real(kind=fPrec)            p
         real(kind=fPrec)            dpop
         real(kind=fPrec)            s
         real(kind=fPrec)            a_eq,b_eq,c_eq,Delta
         real(kind=fPrec)            ENOM
!
         real(kind=fPrec)            x_PRINT,xp_PRINT,y_PRINT,yp_PRINT
         real(kind=fPrec)            AX,BX,AY,BY
         real(kind=fPrec)            X_NORM,XP_NORM,Y_NORM,YP_NORM
         real(kind=fPrec)            EMITX0
         real(kind=fPrec)            EMITY0
!
         integer         LHIT(MAX_NPART)
         integer         LHIT_TURN(MAX_NPART)
         integer         PART_ABS(MAX_NPART)
         integer         PART_ABS_TURN(MAX_NPART)
!
!
         real(kind=fPrec)            x_in(MAX_NPART)
         real(kind=fPrec)            xp_in(MAX_NPART)
         real(kind=fPrec)            y_in(MAX_NPART)
         real(kind=fPrec)            yp_in(MAX_NPART)
         real(kind=fPrec)            p_in(MAX_NPART)    !be careful: [Gev]
         real(kind=fPrec)            s_in(MAX_NPART)
!        adding variables for the pencil beam. Variables in the absolute reference frame.
         real(kind=fPrec)            x_in0(MAX_NPART)
         real(kind=fPrec)            xp_in0(MAX_NPART)
         real(kind=fPrec)            y_in0(MAX_NPART)
         real(kind=fPrec)            yp_in0(MAX_NPART)
         real(kind=fPrec)            p_in0(MAX_NPART)    !be careful: [Gev]
         real(kind=fPrec)            s_in0(MAX_NPART)
         real(kind=fPrec)            s_impact
         integer                     flagsec(MAX_NPART)
         logical                     dowrite_impact

!
         real(kind=fPrec)            IMPACT(MAX_NPART)
         real(kind=fPrec)            INDIV(MAX_NPART)
         real(kind=fPrec)            LINT(MAX_NPART)
         integer                     name(MAX_NPART)
!
         real(kind=fPrec)            x_out(MAX_NPART)
         real(kind=fPrec)            xp_out(MAX_NPART)
         real(kind=fPrec)            y_out(MAX_NPART)
         real(kind=fPrec)            yp_out(MAX_NPART)
         real(kind=fPrec)            p_out(MAX_NPART)
         real(kind=fPrec)            s_out(MAX_NPART)
!
         real(kind=fPrec)            fracab
         real(kind=fPrec)            drift_length
         real(kind=fPrec)            mirror
         real(kind=fPrec)            tiltangle
         real(kind=fPrec)            tiltangle2
!
         real(kind=fPrec)      C_LENGTH       !Length in m
         real(kind=fPrec)      C_ROTATION     !Rotation angle vs vertical in radian
         real(kind=fPrec)      C_APERTURE     !Aperture in m
         real(kind=fPrec)      C_OFFSET       !Offset in m
         real(kind=fPrec)      C_TILT(2)      !Tilt in radian
         real(kind=fPrec)      C_TILT0(2)      !Tilt in radian
         real(kind=fPrec)      cry_bend

!
!
!
        integer ie,iturn,nabs_total
!
!
!
        real(kind=fPrec)        AMPLZ
!
      real(kind=fPrec) XP_tangent
!
      real(kind=fPrec) Cry_tilt_part                    !crystal tilt [rad]
      real(kind=fPrec) Cry_tilt                         !crystal tilt [rad]
      real(kind=fPrec) Cry_tilt0                        !tilt of the crystal for having channeling (determined by divergence of the beam) [rad]
      real(kind=fPrec) Cry_length                       !original length (from the db) [m]

      integer n_chan
      integer n_VR
      integer n_amorphous

! MDA: parameters assignment from new CollDB format

      mirror = 1.0d0

      Rcurv = cdb_cryBend(icoll)
      Alayer = cdb_cryThick(icoll)
      C_xmax = cdb_cryXDim(icoll)
      C_ymax = cdb_cryYDim(icoll)
      C_orient = cdb_cryOrient(icoll)
      miscut = cdb_cryMisCut(icoll)
      mat = cdb_cMaterialID(icoll)

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

      do j = 1, nev
!
        impact(j) = -1.
        lint(j)   = -1.
        indiv(j)  = -1.

        if (ITURN .eq. 1) then                                                          !daniele
                bool_proc_old(j)=-1                                                     !daniele
        else                                                                            !daniele
               if (bool_proc(j).ne.-1)     &
     &                       bool_proc_old(j)=bool_proc(j)                  !daniele
        endif                                                                           !daniele
        iProc=proc_out
        bool_proc(j)=-1                                                                 !daniele
        cry_proc(j)=-1


        nabs = 0

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

!        write(*,*) "p at entrance", p




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
!           write(*,*) j,'- s=',s
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
!          write(*,*) "debug - S shift" ,  S_shift
!          write(*,*) "debug - X shift" ,  X_shift
!
!    2nd transformation: rotation
         S_rot =X_shift*sin(Cry_tilt)+S_shift*cos(Cry_tilt)
         X_rot = X_shift*cos(Cry_tilt)-S_shift*sin(Cry_tilt)
         XP_rot= XP - Cry_tilt
!          write(*,*) "debug - cry tilt", Cry_tilt
!          write(*,*) "debug - S rot" ,  S_rot
!          write(*,*) "debug - X rot" ,  X_rot
!          write(*,*) "debug - XP rot" ,  XP_rot
!    3rd transformation: drift to the new coordinate s=0
         XP=XP_rot
         X= X_rot - XP_rot*S_rot
         Z= Z - ZP*S_rot
         S=0
!          write(*,*) "debug - S cryRF" ,  S_rot
!          write(*,*) "debug - X cryRF" ,  X_rot
!          write(*,*) "debug - XP cryRF" ,  XP_rot
!
!  NOW CHECK IF THE PARTICLE HIT the crystal
!

!          write(*,*) "debug - checking if crystal hit"

         if (x .ge. 0 .and. x.lt.C_xmax) then !Daniele: check that part. hit cry
           s_impact=s_in0(j) !(for the first impact)
!           write(*,*) "HIT", x, C_xmax
!          write(*,*) "HIT"
!           stop
           write(*,*)'hit the cry entrance face', x, C_xmax
!           write(*,*)'impact at s,x = ', s_impact,x_in0(j)
!           write(*,*)'with angle xp = ',xp
!           write(*,*)'s before', s
           CALL CRYST(mat,X,XP,Z,ZP,p,cry_length,j)
!           write(*,*) "p after exit", p
           write(*,*) "xp after crystal coll routine exit", XP
           s=Rcurv*sin(cry_bend)
           zlm=Rcurv*sin(cry_bend)
           if(iProc /= proc_out) then
             NHIT = NHIT + 1
             LHIT(j) = ie
             LHIT_TURN(j) = ITURN
             IMPACT(j) = X_in0(j)
             INDIV(j) = XP_in0(j)
           endif
         else
           write(*,*) "NOT hit the cry entrance face", x, C_xmax
           if (x < 0) then ! Marco: crystal hit from below
            write(*,*) "Crystal hit from below"
            XP_tangent=sqrt((-2*X*Rcurv+X**2)/(Rcurv**2))
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
                          CALL CRYST(mat,X,XP,Z,ZP,p,(cry_length-(tilt_int*Rcurv)),j)
            !              write(*,*) "p after exit", p
                          s=Rcurv*sin(cry_bend-tilt_int)
                          zlm=Rcurv*sin(cry_bend-tilt_int)
                          if(iProc /= proc_out) then
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
           else ! Marco: crystal hit from above
            write(*,*) "Crystal hit from above"
            XP_tangent=asin((Rcurv*(1-cos(cry_bend))-x)/sqrt(2*Rcurv*(Rcurv-x)*(1-cos(cry_bend))+x**2))
            !           stop
            !           write(*,*)j,'-','tangent',xp_tangent,'angle',xp
            !           write(*,*) "Rcurv", Rcurv, "cry_bend", cry_bend
            !           write(*,*)'s tan',Rcurv*sin(XP_tangent)
            !           write(*,*) 's tot', c_length,Rcurv*sin(cry_bend)
                       if ( XP .le. XP_tangent  ) then

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
                          CALL CRYST(mat,X,XP,Z,ZP,p,(cry_length-(tilt_int*Rcurv)),j)
            !              write(*,*) "p after exit", p
                          s=Rcurv*sin(cry_bend-tilt_int)
                          zlm=Rcurv*sin(cry_bend-tilt_int)
                          if(iProc /= proc_out) then
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
!                       stop
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


!  ----------------------DANIELE----------
! debugged assignation of the process experienced in the crystal
! ----------------

              bool_proc(j)=iProc
              cry_proc(j)=iProc
              if(iProc == proc_AM) then
                n_amorphous = n_amorphous + 1
              else if(iProc == proc_VR) then
                n_VR = n_VR + 1
              else if(iProc == proc_CH) then
                n_chan = n_Chan + 1
              else if(iProc == proc_absorbed) then
                nabs = 1
              else if(iProc == proc_ch_absorbed) then
                nabs = 1
              end if

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

!        write(*,*) "XP", XP, mirror
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
           XP_IN(J) = XP_IN(J) + MIRROR*XP_PENCIL0
           YP_IN(J) = YP_IN(J) + MIRROR*YP_PENCIL0
           X_IN(J) = X_IN(J) + MIRROR*X_PENCIL(ICOLL)
           Y_IN(J) = Y_IN(J) + MIRROR*Y_PENCIL(ICOLL)

           IF (.NOT. CHANGED_TILT1(ICOLL) .AND. MIRROR.GT.0.) THEN
                   WRITE (*,*) 'NEVER!!!'
                   C_TILT(1) = XP_PENCIL0*COS(C_ROTATION)+SIN(C_ROTATION)*YP_PENCIL0
                   WRITE(*,*) 'INFO> Changed tilt1  ICOLL  to  ANGLE  ',ICOLL, C_TILT(1)
!
                   CHANGED_TILT1(ICOLL) = .true.
           ELSEIF (.NOT. CHANGED_TILT2(ICOLL) .AND. MIRROR.LT.0.) THEN
                   C_TILT(2) = -1.*(XP_PENCIL0*COS(C_ROTATION)+SIN(C_ROTATION)*YP_PENCIL0)
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
             if(flagsec(j).eq.0 .and. iProc /= proc_out) then
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
!      write(*,*) p

      endif  !collimator with length = 0
      return
      end


! ================================================================================================ !
!  Subroutine for the movements of the particles in the crystal
!  Simple tranport protons in crystal 2
!   J         - number of element
!   S         - longitudinal coordinate
!   IS        - number of substance 1-4: Si,W,C,Ge(110)
!   x,xp,y,yp - coordinates at input of crystal
!   PC        - momentum of particle*c [GeV]
!   W         - weigth of particle
! ================================================================================================ !
subroutine cryst(is,x,xp,y,yp,pc,length,j)

  use mod_ranlux
  use mod_funlux
  use mod_common_main
  use floatPrecision
  use coll_crystal
  use coll_materials, only : zatom, exenergy, rho, anuc

  implicit none

  integer,          intent(in)    :: is
  real(kind=fPrec), intent(inout) :: x
  real(kind=fPrec), intent(inout) :: xp
  real(kind=fPrec), intent(inout) :: y
  real(kind=fPrec), intent(inout) :: yp
  real(kind=fPrec), intent(inout) :: pc
  real(kind=fPrec), intent(in)    :: length
  integer,          intent(in)    :: j

  integer nam,zn                        ! Switch on/off the nuclear interaction (NAM) and the MCS (ZN)
  real(kind=fPrec) ymax,ymin            ! Crystal geometrical parameters
  real(kind=fPrec) s_length             ! Element length along s
  real(kind=fPrec) DESt                 ! Changed energy loss by ionization now calculated and not tabulated
  real(kind=fPrec) x0,y0                ! Coordinates of the particle [m]
  real(kind=fPrec) s                    ! Long coordinates of the particle [m]
  real(kind=fPrec) a_eq,b_eq,c_eq,Delta ! Second order equation param.
  real(kind=fPrec) Ang_rms, Ang_avr     ! Volume reflection mean angle [rad]
  real(kind=fPrec) c_v1                 ! Fitting coefficient
  real(kind=fPrec) c_v2                 ! Fitting coefficient
  real(kind=fPrec) Dechan               ! Probability for dechanneling
  real(kind=fPrec) Lrefl, Srefl         ! Distance of the reflection point [m]
  real(kind=fPrec) Vcapt                ! Volume capture probability
  real(kind=fPrec) Chann                ! Channeling probability
  real(kind=fPrec) N_atom               ! Probability for entering channeling near atomic planes
  real(kind=fPrec) Dxp                  ! Variation in angle
  real(kind=fPrec) xpcrit               ! Critical angle for curved crystal[rad]
  real(kind=fPrec) xpcrit0              ! Critical angle for str. crystal [rad]
  real(kind=fPrec) Rcrit                ! Critical curvature radius [m]
  real(kind=fPrec) ratio                ! X=Rcurv/Rcrit
  real(kind=fPrec) TLdech2              ! Tipical dechanneling length(1) [m]
  real(kind=fPrec) TLdech1              ! Tipical dechanneling length(2) [m]
  real(kind=fPrec) tdech, Ldech,Sdech   ! Angle, lenght, and S coordinate of dechanneling point
  real(kind=fPrec) Rlength, Red_S       ! Reduced length/s coordinate (in case of dechanneling)
  real(kind=fPrec) Am_length            ! Amorphous length
  real(kind=fPrec) Length_xs, Length_ys ! Amorphous length
  real(kind=fPrec) xp_rel               ! Xp-miscut angle in mrad
  real(kind=fPrec) alpha                ! Par for new chann prob
  real(kind=fPrec) Pvr                  ! Prob for VR->AM transition

  real(kind=fPrec) const_dech,xpin,ypin,tchan,L_chan

  nam = 1 ! Switch on/off the nuclear interaction (NAM) and the MCS (ZN)
  zn  = 1

  ! dE/dX and dechanneling length calculation
  mom    = pc*c1e3               ! [GeV]
  enr    = (mom*mom+mp*mp)**half ! [MeV]
  gammar = enr/mp
  betar  = mom/enr
  bgr    = betar*gammar

  tmax = (two*me*bgr**2)/(one + two*gammar*me/mp+(me/mp)**2)  ! [MeV]
  plen = ((rho(is)*zatom(is)/anuc(is))**half)*28.816e-6_fPrec ! [MeV]

  const_dech = (256.0/(9.0*pi**2))* &
    (one/(log(two*me*gammar/(exenergy(is)*c1e3))-one))*((aTF*dP)/(re*me)) ! [m/MeV]
  const_dech = const_dech*c1e3 ! [m/GeV]

  s        = 0
  s_length = Rcurv*(sin(length/Rcurv))
  L_chan   = length

  if(miscut < zero .and. x > zero .and. x < -length*tan(miscut)) then
    L_chan = -x/sin(miscut)
  end if

  tchan   = L_chan/Rcurv
  xp_rel = xp-miscut

  ymin = -C_ymax/two
  ymax =  C_ymax/two

  ! FIRST CASE: p don't interact with crystal
  if(y < ymin .or. y > ymax .or. x > c_xmax) then
    x     = x + xp*s_length
    y     = y + yp*s_length
    iProc = proc_out
    return

  ! SECOND CASE: p hits the amorphous layer
  else if(x < alayer .or. y-ymin < alayer .or. ymax-y < alayer) then
    x0    = x
    y0    = y
    a_eq  = one + xp**2
    b_eq  = two*x*xp - two*xp*Rcurv
    c_eq  = x**2 - two*x*Rcurv
    Delta = b_eq**2 - four*a_eq*c_eq
    s     = (-b_eq+sqrt(Delta))/(two*a_eq)
    if(s >= s_length) then
      s = s_length
    end if
    x         =  xp*s + x0
    Length_xs = sqrt((x-x0)**2 + s**2)
    if(yp >= zero .and. y + yp*s <= ymax) then
      Length_ys = yp*Length_xs
    else if(yp < zero .and. y + yp*s >= ymin) then
      Length_ys = yp*Length_xs
    else
      s         = (ymax-y)/yp
      Length_ys = sqrt((ymax-y)**2 + s**2)
      x         = x0 + xp*s
      Length_xs = sqrt((x-x0)**2 + s**2)
    end if
    Am_length = sqrt(Length_xs**2 + Length_ys**2)
    s     = s/two
    x     = x0 + xp*s
    y     = y0 + yp*s
    iProc = proc_AM
    call calc_ion_loss_cry(is,pc,am_length,dest)
    call move_am(is,nam,am_length,dest,dlyi(is),dlri(is),xp,yp,pc)
    x = x + xp*(s_length-s)
    y = y + yp*(s_length-s)
    return

  elseif ((x.gt.(c_xmax-alayer)) .and. x.lt.(c_xmax)  ) then
    iProc = proc_AM
    call calc_ion_loss_cry(is,pc,s_length,dest)
    call move_am(is,nam,s_length,dest,dlyi(is),dlri(is),xp,yp,pc)
    return

  end if

  ! THIRD CASE: the p interacts with the crystal.
  ! Define typical angles/probabilities for orientation 110
  xpcrit0 = (2.0e-9*eUm(is)/pc)**half  ! Critical angle (rad) for straight crystals
  Rcrit   = pc/(2.0e-6*eum(is))*ai(is) ! Critical curvature radius [m]

  ! If R>Rcritical=>no channeling is possible (ratio<1)
  ratio  = Rcurv/Rcrit
  xpcrit = xpcrit0*(Rcurv-Rcrit)/Rcurv ! Critical angle for curved crystal

  c_v1 =  1.7 ! Fitting coefficient ???
  c_v2 = -1.5 ! Fitting coefficient ???

  if(ratio <= one) then ! no possibile channeling
    Ang_rms = c_v1*0.42*xpcrit0*sin(1.4*ratio)  ! rms scattering
    Ang_avr = c_v2*xpcrit0*0.05*ratio           ! average angle reflection
    Vcapt   = zero                              ! probability of VC

  else if(ratio <= three) then ! Strongly bent xstal
    Ang_rms = c_v1*0.42*xpcrit0*sin(1.571*0.3*ratio+0.85) ! rms scattering
    Ang_avr = c_v2*xpcrit0*(0.1972*ratio-0.1472)          ! avg angle reflection
    Vcapt   = 0.0007*(ratio-0.7)/PC**0.2                  ! correction by sasha drozdin/armen
    ! K=0.00070 is taken based on simulations using CATCH.f (V.Biryukov)

  else ! Rcry >> Rcrit
    Ang_rms = c_v1*xpcrit0*(1./ratio)        ! rms scattering
    Ang_avr = c_v2*xpcrit0*(1.-1.6667/ratio) ! average angle for VR
    Vcapt   = 0.0007*(ratio-0.7)/PC**0.2     ! probability for VC correction by sasha drozdin/armen
    ! K=0.0007 is taken based on simulations using CATCH.f (V.Biryukov)

  end if

  if(C_orient == 2) then
    Ang_avr = Ang_avr * 0.93
    Ang_rms = Ang_rms * 1.05
    xpcrit  = xpcrit * 0.98
  end if

  if(abs(xp_rel) < xpcrit) then
    alpha  = xp_rel/xpcrit
    Chann  = ((0.9-alpha*alpha*0.9)**half)*(one-(one/ratio))**half ! Saturation at 95%
    N_atom = 0.1

    ! if they can channel: 2 options
    if(rndm4() <= chann) then ! option 1:channeling

      TLdech1 = const_dech*PC*(one-one/ratio)**2 ! Updated calculate typical dech. length(m)
      if(rndm4() <= n_atom) then
        TLdech1 = (const_dech/200.d0)*PC*(one-one/ratio)**2  ! Updated dechanneling length (m)
      end if

      Dechan = -log(rndm4())  ! probability of dechanneling
      Ldech  = TLdech1*Dechan ! actual dechan. length

      ! careful: the dechanneling lentgh is along the trajectory
      ! of the particle -not along the longitudinal coordinate...
      if(ldech < l_chan) then
        iProc = proc_DC
        Dxp   = Ldech/Rcurv ! change angle from channeling [mrad]
        Sdech = Ldech*cos(miscut + half*Dxp)
        x     = x + Ldech*(sin(half*Dxp+miscut)) ! trajectory at channeling exit
        xp    = xp + Dxp + two*(rndm4()-half)*xpcrit
        y     = y + yp * Sdech

        call calc_ion_loss_cry(is,pc,ldech,dest)
        pc = pc - half*dest*Ldech ! energy loss to ionization while in CH [GeV]
        x  = x  + half*(s_length-Sdech)*xp
        y  = y  + half*(s_length-Sdech)*yp

        call calc_ion_loss_cry(is,pc,s_length-sdech,dest)
        call move_am(is,nam,s_length-sdech,dest,dlyi(is),dlri(is),xp,yp,pc)
        x = x + half*(s_length-Sdech)*xp
        y = y + half*(s_length-Sdech)*yp
      else
        iProc = proc_CH
        xpin  = XP
        ypin  = YP

        call move_ch(is,nam,l_chan,x,xp,yp,pc,rcurv,rcrit) ! check if a nuclear interaction happen while in CH
        if(iProc /= proc_CH) then
          ! if an nuclear interaction happened, move until the middle with initial xp,yp then
          ! propagate until the "crystal exit" with the new xp,yp accordingly with the rest
          ! of the code in "thin lens approx"
          x = x + half*L_chan*xpin
          y = y + half*L_chan*ypin
          x = x + half*L_chan*XP
          y = y + half*L_chan*YP

          call calc_ion_loss_cry(is,pc,length,dest)
          pc = pc - dest*length ! energy loss to ionization [GeV]
        else
          Dxp = L_chan/Rcurv + half*ran_gauss(one)*xpcrit ! change angle[rad]
          xp  = Dxp
          x   = x+ L_chan*(sin(half*Dxp+miscut)) ! trajectory at channeling exit
          y   = y + s_length * yp

          call calc_ion_loss_cry(is,pc,length,dest)
          pc = pc - half*dest*length ! energy loss to ionization [GeV]
        end if
      end if

    else ! Option 2: VR

      ! good for channeling but don't channel (1-2)
      iProc = proc_VR

      xp = xp + 0.45*(xp/xpcrit+1)*Ang_avr
      x  = x  + half*s_length*xp
      y  = y  + half*s_length*yp

      call calc_ion_loss_cry(is,pc,s_length,dest)
      call move_am(is,nam,s_length,dest,dlyi(is),dlri(is),xp,yp,pc)

      x = x + half*s_length*xp
      y = y + half*s_length*yp

    end if

  else ! case 3-2: no good for channeling. check if the  can VR

    Lrefl = xp_rel*Rcurv ! distance of refl. point [m]
    Srefl = sin(xp_rel/two + miscut)*Lrefl

    if(Lrefl > zero .and. Lrefl < Length) then ! VR point inside

      ! 2 options: volume capture and volume reflection

      if(rndm4() > Vcapt .or. ZN == zero) then ! Option 1: VR

        iProc = proc_VR
        x     = x + xp*Srefl
        y     = y + yp*Srefl
        Dxp   = Ang_avr
        xp    = xp + Dxp + Ang_rms*ran_gauss(one)
        x     = x + half*xp*(s_length - Srefl)
        y     = y + half*yp*(s_length - Srefl)

        call calc_ion_loss_cry(is,pc,s_length-srefl,dest)
        call move_am(is,nam,s_length-srefl,dest,dlyi(is),dlri(is),xp,yp,pc)
        x = x + half*xp*(s_length - Srefl)
        y = y + half*yp*(s_length - Srefl)

      else ! Option 2: VC

        x = x + xp*Srefl
        y = y + yp*Srefl

        TLdech2 = (const_dech/c1e1)*pc*(one-one/ratio)**2    ! updated typical dechanneling length(m)
        Ldech   = TLdech2*(sqrt(0.01-log(rndm4())) - 0.1)**2 ! updated DC length
        tdech   = Ldech/Rcurv
        Sdech   = Ldech*cos(xp + half*tdech)

        if(Ldech < Length-Lrefl) then

          iProc = proc_DC
          Dxp   = Ldech/Rcurv + half*ran_gauss(one)*xpcrit
          x     = x + Ldech*(sin(half*Dxp+xp)) ! trajectory at channeling exit
          y     = y + Sdech*yp
          xp    =  Dxp
          Red_S = s_length - Srefl - Sdech
          x     = x + half*xp*Red_S
          y     = y + half*yp*Red_S

          call calc_ion_loss_cry(is,pc,srefl,dest)
          pc = pc - dest*Srefl ! "added" energy loss before capture

          call calc_ion_loss_cry(is,pc,sdech,dest)
          pc = pc - half*dest*Sdech ! "added" energy loss while captured

          call calc_ion_loss_cry(is,pc,red_s,dest)
          call move_am(is,nam,red_s,dest,dlyi(is),dlri(is),xp,yp,pc)
          x = x + half*xp*Red_S
          y = y + half*yp*Red_S

        else

          iProc   = proc_VC
          Rlength = Length-Lrefl
          tchan   = Rlength/Rcurv
          Red_S   = Rlength*cos(xp + half*tchan)

          call calc_ion_loss_cry(is,pc,lrefl,dest)
          pc   = pc - dest*Lrefl ! "added" energy loss before capture
          xpin = xp
          ypin = yp

          call move_ch(is,nam,rlength,x,xp,yp,pc,rcurv,rcrit) ! check if a nuclear interaction happen while in ch
          if(iProc /= proc_VC) then
            ! if an nuclear interaction happened, move until the middle with initial xp,yp then propagate until
            ! the "crystal exit" with the new xp,yp accordingly with the rest of the code in "thin lens approx"
            x = x + half*Rlength*xpin
            y = y + half*Rlength*ypin
            x = x + half*Rlength*XP
            y = y + half*Rlength*YP

            call calc_ion_loss_cry(is,pc,rlength,dest)
            pc = pc - dest*Rlength
          else
            Dxp = (Length-Lrefl)/Rcurv
            x   = x + sin(half*Dxp+xp)*Rlength ! trajectory at channeling exit
            y   = y + red_S * yp
            xp  = Length/Rcurv + half*ran_gauss(one)*xpcrit ! [mrad]

            call calc_ion_loss_cry(is,pc,rlength,dest)
            pc = pc - half*dest*Rlength  ! "added" energy loss once captured
          end if
        end if
      end if

    else

      ! Case 3-3: move in amorphous substance (big input angles)
      ! Modified for transition vram daniele
      if(xp_rel > L_chan/Rcurv + two*xpcrit .or. xp_rel < -xpcrit) then
        iProc = proc_AM
        x     = x + half*s_length*xp
        y     = y + half*s_length*yp
        if(zn > zero) then
          call calc_ion_loss_cry(is,pc,s_length,dest)
          call move_am(is,nam,s_length,dest,dlyi(is),dlri(is),xp,yp,pc)
        end if
        x = x + half*s_length*xp
        y = y + half*s_length*yp
      else
        Pvr = (xp_rel-(L_chan/Rcurv))/(two*xpcrit)
        if(rndm4() > Pvr) then
          iProc = proc_TRVR
          x     = x + xp*Srefl
          y     = y + yp*Srefl

          Dxp = -three*Ang_rms*xp_rel/(two*xpcrit) + Ang_avr + (three*Ang_rms*(L_chan/Rcurv)/(two*xpcrit))
          xp  = xp+Dxp
          x   = x + half*xp*(s_length-Srefl)
          y   = y + half*yp*(s_length-Srefl)

          call calc_ion_loss_cry(is,pc,s_length-srefl,dest)
          call move_am(is,nam,s_length-srefl,dest,dlyi(is),dlri(is),xp,yp,pc)
          x = x + half*xp*(s_length - Srefl)
          y = y + half*yp*(s_length - Srefl)
        else
          iProc = proc_TRAM
          x = x + xp*Srefl
          y = y + yp*Srefl
          Dxp = -one*(13.6/pc)*sqrt(s_length/dlri(is))*0.001*xp_rel/(two*xpcrit) + &
            (13.6/pc)*sqrt(s_length/DLRi(is))*0.001*(one+(L_chan/Rcurv)/(two*xpcrit))
          xp = xp+Dxp
          x  = x + half*xp*(s_length-Srefl)
          y  = y + half*yp*(s_length-Srefl)

          call calc_ion_loss_cry(is,pc,s_length-srefl,dest)
          call move_am(is,nam,s_length-srefl,dest,dlyi(is),dlri(is),xp,yp,pc)
          x = x + half*xp*(s_length - Srefl)
          y = y + half*yp*(s_length - Srefl)
        end if
      end if
    end if
  end if

end subroutine cryst

! ================================================================================================ !
!  Subroutine for the calculazion of the energy loss by ionisation
! ================================================================================================ !
subroutine calc_ion_loss_cry(is,pc,dz,EnLo)

  use mod_ranlux
  use mod_funlux
  use floatPrecision
  use coll_crystal
  use coll_materials, only : zatom, exenergy, rho, anuc

  implicit none

  integer,          intent(in)  :: is
  real(kind=fPrec), intent(in)  :: pc
  real(kind=fPrec), intent(in)  :: dz
  real(kind=fPrec), intent(out) :: EnLo

  real(kind=fPrec) thl,tt,cs_tail,prob_tail
  real(kind=fPrec), parameter :: k = 0.307075 ! Constant in front bethe-bloch [mev g^-1 cm^2]

  thl       = four*k*zatom(is)*dz*c1e2*rho(is)/(anuc(is)*betar**2) ! [MeV]
  EnLo      = ((k*zatom(is))/(anuc(is)*betar**2))*( &
    half*log((two*me*bgr*bgr*Tmax)/(c1e6*exenergy(is)**2)) - &
    betar**2-log(plen/(exenergy(is)*c1e3))-log(bgr)+half)
  EnLo      = EnLo*rho(is)*0.1*dz ! [GeV]
  Tt        = EnLo*c1e3+thl ! [MeV]
  cs_tail   = ((k*zatom(is))/(anuc(is)*betar**2))*((half*((one/Tt)-(one/Tmax))) - &
    (log(Tmax/Tt)*(betar**2)/(two*Tmax)) + ((Tmax-Tt)/(four*(gammar**2)*(mp**2))))
  prob_tail = cs_tail*rho(is)*dz*c1e2

  if(rndm4() < prob_tail) then
    EnLo = ((k*zatom(is))/(anuc(is)*betar**2))*( &
      half*log((two*me*bgr*bgr*Tmax)/(c1e6*exenergy(is)**2)) - &
      betar**2-log(plen/(exenergy(is)*c1e3))-log(bgr)+half   + &
      (TMax**2)/(eight*(gammar**2)*(mp**2)))
    EnLo = EnLo*rho(is)*0.1 ! [GeV/m]
  else
    EnLo = EnLo/dz ! [GeV/m]
  end if

end subroutine calc_ion_loss_cry

! ================================================================================================ !
!  Subroutine for the movement in the amorphous
! ================================================================================================ !
subroutine move_am(is,nam,dz,dei,dly,dlr,xp,yp,pc)

  use mod_ranlux
  use mod_funlux
  use floatPrecision
  use coll_crystal
  use coll_materials, only : anuc, hcut, bnref, csref

  implicit none

  integer,          intent(in)    :: is
  integer,          intent(in)    :: nam
  real(kind=fPrec), intent(in)    :: dz
  real(kind=fPrec), intent(in)    :: dei
  real(kind=fPrec), intent(in)    :: dly
  real(kind=fPrec), intent(in)    :: dlr
  real(kind=fPrec), intent(inout) :: xp
  real(kind=fPrec), intent(inout) :: yp
  real(kind=fPrec), intent(inout) :: pc

  integer i,length_cry,ichoix
  real(kind=fPrec) t,xran_cry(1),bn,cs(0:5),cprob(0:5),freep,zlm,xp_in,yp_in,xm2,xln15s,tz,tx,tlow, &
    thigh,teta,pptot,ppsd,ppel,pc_in,kymcs,kxmcs,ecmsq,dya,bsd,bpp,aran
  real(kind=fPrec), external :: ruth_cry

  xp_in = xp
  yp_in = yp
  pc_in = pc

  ! New treatment of scattering routine based on standard sixtrack routine
  ! useful calculations for cross-section and event topology calculation
  ecmsq  = two*0.93828_fPrec*pc
  xln15s = log(0.15*ecmsq)

  ! New models, see Claudia's thesis
  pptot = 0.041084d0-0.0023302d0*log(ecmsq)+0.00031514d0*log(ecmsq)**2
  ppel  = (11.7d0-1.59d0*log(ecmsq)+0.134d0*log(ecmsq)**2)/1000
  ppsd  = (4.3d0+0.3d0*log(ecmsq))/1000
  bpp   = 7.156d0+1.439d0*log(sqrt(ecmsq))

  ! Distribution for Ruth. scatt.
  tlow      = tlcut_cry
  mcurr_cry = is
  thigh     = hcut(is)
  call funlxp(ruth_cry,cgen_cry(1,is),tlow,thigh)

  ! Cross-section calculation
  ! freep: number of nucleons involved in single scattering
  freep = freeco_cry * anuc(is)**(one/three)

  ! Compute pp and pn el+single diff contributions to cross-section (both added : quasi-elastic or qel later)
  cs(3) = freep*ppel
  cs(4) = freep*ppsd

  ! Correct TOT-CSec for energy dependence of qel
  ! TOT CS is here without a Coulomb contribution
  cs(0) = csref(0,is) + freep*(pptot - pptref_cry)
  bn    = bnref(is)*cs(0)/csref(0,is)

  ! Also correct inel-CS
  cs(1) = csref(1,is)*cs(0)/csref(0,is)

  ! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
  cs(2) = cs(0) - cs(1) - cs(3) - cs(4)
  cs(5) = csref(5,is)

  ! Now add Coulomb
  cs(0) = cs(0) + cs(5)

  ! Calculate cumulative probability
  cprob(:) = zero
  cprob(5) = one
  do i=1,4
    cprob(i) = cprob(i-1)+cs(i)/cs(0)
  end do

  ! Multiple Coulomb Scattering
  xp  = xp*c1e3
  yp  = yp*c1e3
  pc  = pc-dei*dz ! Energy lost because of ionization process[GeV]

  dya   = (13.6/pc)*sqrt(dz/dlr) ! rms of coloumb scattering MCS (mrad)
  kxmcs = dya*ran_gauss(one)
  kymcs = dya*ran_gauss(one)

  xp = xp+kxmcs
  yp = yp+kymcs

  if(nam == 0) return ! Turn on/off nuclear interactions

  ! Can nuclear interaction happen?
  zlm = -collnt(is)*log(rndm4())

  if(zlm < dz) then
    ! Choose nuclear interaction
    aran = rndm4()
    i=1
10  if(aran > cprob(i)) then
      i = i+1
      goto 10
    end if
    ichoix = i

    ! Do the interaction
    t = 0 ! default value to cover ichoix=1
    select case(ichoix)
    case(1) ! Deep inelastic, impinging p disappeared
      iProc = proc_absorbed

    case(2) ! p-n elastic
      iProc = proc_pne
      t     = -log(rndm4())/bn

    case(3) ! p-p elastic
      iProc = proc_ppe
      t     = -log(rndm4())/bpp

    case(4) ! Single diffractive
      iProc = proc_diff
      xm2   = exp(rndm4()*xln15s)
      pc    = pc*(one - xm2/ecmsq)
      if(xm2 < two) then
        bsd = two*bpp
      else if(xm2 >= two .and. xm2 <= five) then
        bsd = (106.d0 - 17.d0*xm2)* bpp/36.d0
      else if(xm2 > five) then
        bsd = 7.d0*bpp/12.d0
      end if
      t = -log(rndm4())/bsd

    case(5)
      iProc      = proc_ruth
      length_cry = 1
      call funlux(cgen_cry(1,is),xran_cry,length_cry)
      t = xran_cry(1)

    end select

    ! Calculate the related kick
    if(ichoix == 4) then
      teta = sqrt(t)/pc_in ! DIFF has changed PC
    else
      teta = sqrt(t)/pc
    end if

    tx = teta*ran_gauss(one)*c1e3
    tz = teta*ran_gauss(one)*c1e3

    ! Change p angle
    xp = xp + tx
    yp = yp + tz
  end if

  xp = xp/c1e3
  yp = yp/c1e3

end subroutine move_am

! ================================================================================================ !
!  Subroutine for check if a nuclear interaction happen while in channeling
! ================================================================================================ !
subroutine move_ch(is,nam,dz,x,xp,yp,pc,r,rc)

  use crcoall
  use mod_ranlux
  use mod_funlux
  use floatPrecision
  use coll_common, only : coll_debug
  use coll_crystal
  use coll_materials, only : nmat, rho, anuc, hcut, bnref, csref, csect

  implicit none

  integer,          intent(in)    :: is
  integer,          intent(in)    :: nam
  real(kind=fPrec), intent(in)    :: dz
  real(kind=fPrec), intent(inout) :: x
  real(kind=fPrec), intent(inout) :: xp
  real(kind=fPrec), intent(inout) :: yp
  real(kind=fPrec), intent(inout) :: pc
  real(kind=fPrec), intent(in)    :: r
  real(kind=fPrec), intent(in)    :: rc

  integer i,np,length_cry,ichoix
  real(kind=fPrec) t,xran_cry(1),bn,cs(0:5),cprob(0:5),freep,zlm,xp_in,yp_in,xminU,xm2,xln15s,x_min,&
    x_max,x_i,Umin,Ueff,tz,tx,tlow,thigh,teta,rho_min,rho_max,pv,pptot,ppsd,ppel,PC_in,nuc_cl_l,    &
    N_am,Et,ecmsq,Ec,csref_inel_rsc,csref_tot_rsc,bsd,bpp,aran,avrrho
  real(kind=fPrec), external :: ruth_cry

  xp_in = xp
  yp_in = yp
  pc_in = pc

  ! New treatment of scattering routine based on standard sixtrack routine

  ! Useful calculations for cross-section and event topology calculation
  ecmsq = 2*0.93828d0*pc
  xln15s=log(0.15*ecmsq)

  ! New models, see Claudia's thesis
  pptot = 0.041084d0-0.0023302d0*log(ecmsq)+0.00031514d0*log(ecmsq)**2
  ppel  = (11.7d0-1.59d0*log(ecmsq)+0.134d0*log(ecmsq)**2)/1000
  ppsd  = (4.3d0+0.3d0*log(ecmsq))/1000
  bpp   = 7.156d0+1.439d0*log(sqrt(ecmsq))

  ! Distribution for Ruth. scatt.
  tlow      = tlcut_cry
  mcurr_cry = is
  thigh     = hcut(is)
  call funlxp(ruth_cry,cgen_cry(1,is),tlow,thigh)

  ! Rescale the total and inelastic cross-section accordigly to the average density seen
  x_i = x
  np  = int(x_i/dp)    ! Calculate in which crystalline plane the particle enters
  x_i = x_i - Np*dP    ! Rescale the incoming x at the left crystalline plane
  x_i = x_i - (dP/two) ! Rescale the incoming x in the middle of crystalline planes

  pv   = pc*c1e9*pc*c1e9/sqrt(pc*c1e9*pc*c1e9 + 93828.0e6_fPrec*93828.0e6_fPrec) ! Calculate pv=P/E
  Ueff = eUm(is)*(two*x_i/dp)*(two*x_i/dp)+pv*x_i/r                              ! Calculate effective potential
  Et   = pv*xp*xp/two+Ueff                                                       ! Calculate transverse energy
  Ec   = eUm(is)*(one-rc/r)*(one-rc/r)                                           ! Calculate critical energy in bent crystals

  ! To avoid negative Et
  xminU = -dp*dp*pc*c1e9/(eight*eUm(is)*r)
  Umin  = abs(eUm(is)*(two*xminU/dp)*(two*xminU/dP)+pv*xminU/R)
  Et    = Et+Umin
  Ec    = Ec+Umin

  ! Calculate min e max of the trajectory between crystalline planes
  x_min = -(dP/two)*Rc/R-(dP/two)*sqrt(Et/Ec)
  x_Max = -(dP/two)*Rc/R+(dP/two)*sqrt(Et/Ec)

  ! Change ref. frame and go back with 0 on the crystalline plane on the left
  x_min = x_min - dp/two
  x_max = x_max - dp/two

  ! Calculate the "normal density" in m^-3
  N_am  = rho(is)*6.022e23_fPrec*c1e6/anuc(is)

  ! Calculate atomic density at min and max of the trajectory oscillation
  rho_max = n_am*dp/two*(erf(x_max/sqrt(two*u1*u1)) - erf((dP-x_Max)/sqrt(two*u1*u1)))
  rho_min = N_am*dP/two*(erf(x_min/sqrt(two*u1*u1)) - erf((dP-x_min)/sqrt(two*u1*u1)))

  ! "zero-approximation" of average nuclear density seen along the trajectory
  avrrho  = (rho_max-rho_min)/(x_max-x_min)
  avrrho  = two*avrrho/N_am

  csref_tot_rsc  = csref(0,is)*avrrho ! Rescaled total ref cs
  csref_inel_rsc = csref(1,is)*avrrho ! Rescaled inelastic ref cs

  ! Cross-section calculation
  freep = freeco_cry * anuc(is)**(one/three)

  ! compute pp and pn el+single diff contributions to cross-section (both added : quasi-elastic or qel later)
  cs(3) = freep*ppel
  cs(4) = freep*ppsd

  ! correct TOT-CSec for energy dependence of qel
  ! TOT CS is here without a Coulomb contribution
  cs(0) = csref_tot_rsc + freep*(pptot - pptref_cry)

  ! Also correct inel-CS
  if(csref_tot_rsc == zero) then
    cs(1) = zero
  else
    cs(1) = csref_inel_rsc*cs(0)/csref_tot_rsc
  end if

  ! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
  cs(2) = cs(0) - cs(1) - cs(3) - cs(4)
  cs(5) = csref(5,is)

  ! Now add Coulomb
  cs(0) = cs(0) + cs(5)

  ! Calculate cumulative probability
  cprob(:) = zero
  cprob(5) = one
  if(cs(0) == zero) then
    do i=1,4
      cprob(i) = cprob(i-1)
    end do
  else
    do i=1,4
      cprob(i) = cprob(i-1) + cs(i)/cs(0)
    end do
  end if

  ! Multiple Coulomb Scattering
  xp = xp*c1e3
  yp = yp*c1e3

  ! Turn on/off nuclear interactions
  if(nam == 0) return

  ! Can nuclear interaction happen?
  ! Rescaled nuclear collision length
  if(avrrho == zero) then
    nuc_cl_l = c1e6
  else
    nuc_cl_l = collnt(is)/avrrho
  end if
  zlm = -nuc_cl_l*log(rndm4())

  ! write(889,*) x_i,pv,Ueff,Et,Ec,N_am,avrrho,csref_tot_rsc,csref_inel_rsc,nuc_cl_l

  if(zlm < dz) then
    ! Choose nuclear interaction
    aran = rndm4()
    i=1
10  if(aran > cprob(i)) then
      i=i+1
      goto 10
    end if
    ichoix = i

    ! Do the interaction
    select case(ichoix)
    case(1) ! deep inelastic, impinging p disappeared
      iProc = proc_ch_absorbed

    case(2) ! p-n elastic
      iProc = proc_ch_pne
      bn    = bnref(is)*cs(0)/csref_tot_rsc
      t     = -log(rndm4())/bn

    case(3) ! p-p elastic
      iProc = proc_ch_ppe
      t     = -log(rndm4())/bpp

    case(4) ! Single diffractive
      iProc = proc_ch_diff
      xm2   = exp(rndm4()*xln15s)
      pc    = pc*(one - xm2/ecmsq)
      if(xm2 < two) then
        bsd = two*bpp
      else if(xm2 >= two .and. xm2 <= five) then
        bsd = (106.d0 - 17.d0*xm2) * bpp / 36.d0
      else if(xm2 > five) then
        bsd = seven*bpp / 12.d0
      end if
      t = -log(rndm4())/bsd

    case(5)
      iProc      = proc_ch_ruth
      length_cry = 1
      call funlux(cgen_cry(1,is),xran_cry,length_cry)
      t = xran_cry(1)

    end select

    ! Calculate the related kick -----------
    if(ichoix == 4) then
      teta = sqrt(t)/PC_in ! DIFF has changed PC!!!
    else
      teta = sqrt(t)/PC
    end if

    tx = teta*ran_gauss(one)*c1e3
    tz = teta*ran_gauss(one)*c1e3

    ! Change p angle
    xp = xp + tx
    yp = yp + tz

  end if

  xp = xp/c1e3
  yp = yp/c1e3

end subroutine move_ch

! ================================================================================================ !
! Definition of rutherford scattering formula
! ================================================================================================ !
function ruth_cry(t_cry)

  use floatPrecision
  use coll_crystal
  use coll_materials

  implicit none

  real(kind=fPrec) ruth_cry,t_cry
  real(kind=fPrec), parameter :: cnorm  = 2.607e-4_fPrec
  real(kind=fPrec), parameter :: cnform = 0.8561e3_fPrec

  ruth_cry = cnorm*exp(-t_cry*cnform*emr(mcurr_cry)**2)*(zatom(mcurr_cry)/t_cry)**2

end function ruth_cry
