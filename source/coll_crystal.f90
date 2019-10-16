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
  real(kind=fPrec), save      :: cgen_cry(200,4)
  integer,          save      :: mcurr_cry

  real(kind=fPrec), save :: enr,mom,betar,gammar,bgr !Daniele: energy,momentum,beta relativistic, gamma relativistic
  real(kind=fPrec), save :: Tmax,plen !Daniele: maximum energy tranfer in single collision, plasma energy (see pdg)

  real(kind=fPrec), parameter :: re  = 2.818d-15   ! electron radius [m]
  real(kind=fPrec), parameter :: me  = 0.510998910 ! electron mass [MeV/c^2]
  real(kind=fPrec), parameter :: mp  = 938.272013  ! proton mass [MeV/c^2]

  real(kind=fPrec), parameter :: aTF = 0.194d-10   ! Screening function [m]
  real(kind=fPrec), parameter :: dP  = 1.92d-10    ! Distance between planes (110) [m]
  real(kind=fPrec), parameter :: u1  = 0.075d-10   ! Thermal vibrations amplitude

  ! Nuclear Collision length [m] from pdg (only for Si and Ge for the moment)
  real(kind=fPrec), parameter :: collnt_cry(4) = [0.3016d0,0.0d0,0.0d0,0.1632d0]

  real(kind=fPrec), public, save :: cprob_cry(0:5,1:4) = reshape([ &
    [zero, zero, zero, zero, zero, one], &
    [zero, zero, zero, zero, zero, one], &
    [zero, zero, zero, zero, zero, one], &
    [zero, zero, zero, zero, zero, one]  &
  ], shape=[6,4])

  ! pp cross-sections and parameters for energy dependence
  real(kind=fPrec), parameter :: pptref_cry = 0.04d0
  real(kind=fPrec), parameter :: pperef_cry = 0.007d0
  real(kind=fPrec), parameter :: sdcoe_cry  = 0.00068d0
  real(kind=fPrec), parameter :: pref_cry   = 450.0d0
  real(kind=fPrec), parameter :: pptco_cry  = 0.05788d0
  real(kind=fPrec), parameter :: ppeco_cry  = 0.04792d0
  real(kind=fPrec), parameter :: freeco_cry = 1.618d0

  ! Crystal Specific Material Arrays
  real(kind=fPrec), save :: dlri(nmat) = zero
  real(kind=fPrec), save :: dlyi(nmat) = zero
  real(kind=fPrec), save :: ai(nmat)   = zero
  real(kind=fPrec), save :: eUm(nmat)  = zero

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
  dlri(m) = 0.0937
  dlyi(m) = 0.4652
  ai(m)   = 0.96e-7
  eUm(m)  = 21.34

  ! W
  m = collmat_getCollMatID("W")
  dlri(m) = 0.0035
  dlyi(m) = 0.096
  ai(m)   = 0.56e-7
  eUm(m)  = 21.0

  ! C
  m = collmat_getCollMatID("C")
  dlri(m) = 0.188
  dlyi(m) = 0.400
  ai(m)   = 0.63e-7
  eUm(m)  = 21.0

  ! Ge
  m = collmat_getCollMatID("Ge")
  dlri(m) = 0.02302
  dlyi(m) = 0.2686
  ai(m)   = 1.0e-7
  eUm(m)  = 40.0

end subroutine cry_init

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
         CHARACTER(LEN=4) C_MATERIAL     !Material
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
        integer ie,iturn,nabs_total,new_mat
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
      character(len=*) name_coll             ! MDA: length to be defined by input
                           !string that contains the physical process
!                                            !=' ' if the particle does not pass by crystal, ='*' if there is interaction

! MDA: parameters assignment from new CollDB format

      mirror = 1.0d0

      Rcurv = cdb_cryBend(icoll)
      Alayer = cdb_cryThick(icoll)
      C_xmax = cdb_cryXDim(icoll)
      C_ymax = cdb_cryYDim(icoll)
      C_orient = cdb_cryOrient(icoll)
      miscut = cdb_cryMisCut(icoll)

      IF (C_MATERIAL(1:4).eq.'Si')THEN           ! MDA: changed the label for all materials to fit new format
           mat = 8
           new_mat = 8
      ELSEIF (C_MATERIAL(1:4).eq.'W')THEN
           mat = 9
           new_mat = 4
      ELSEIF (C_MATERIAL(1:4).eq.'C')THEN
           mat = 10
           new_mat = 6
      ELSEIF (C_MATERIAL(1:4).eq.'Ge')THEN
           mat = 11
           new_mat = 9
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
           CALL CRYST(mat-7,new_mat,X,XP,Z,ZP,p,cry_length,j)
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
                          CALL CRYST(mat-7,new_mat,X,XP,Z,ZP,p,(cry_length-(tilt_int*Rcurv)),j)
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
                          CALL CRYST(mat-7,new_mat,X,XP,Z,ZP,p,(cry_length-(tilt_int*Rcurv)),j)
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


!.**************************************************************************
!     SUBROUTINE FOR THE MOVEMENTS OF THE PARTICLES IN THE CRYSTAL
!.**************************************************************************
      SUBROUTINE CRYST(IS,is2,x,xp,y,yp,PC,Length,j)

      use mod_ranlux
      use mod_funlux
      use mod_common_main
      use floatPrecision
      use coll_crystal
      use coll_materials, only : zatom, exenergy, rho, anuc

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
      real(kind=fPrec) Length !crystal geometrical parameters
                                                  ! [m],[m],[m],[m],[rad]
      real(kind=fPrec) ymax,ymin       !crystal geometrical parameters
      real(kind=fPrec) s_length             !element length along s
      integer IS,j,is2                            !index of the material
!      integer counter
      real(kind=fPrec) DESt                  ! Daniele: changed energy loss by ionization now calculated and not tabulated
      integer NAM,ZN                        !switch on/off the nuclear
                                            !interaction (NAM) and the MCS (ZN)
      real(kind=fPrec) x,xp,y,yp,PC         !coordinates of the particle
                                            ![m],[rad],[m],[rad],[GeV]
      real(kind=fPrec) x0,y0                !coordinates of the particle [m]
      real(kind=fPrec) s                    !long coordinates of the particle [m]
      real(kind=fPrec) a_eq,b_eq,c_eq,Delta !second order equation param.
      real(kind=fPrec) Ang_rms, Ang_avr     !Volume reflection mean angle [rad]
      real(kind=fPrec) c_v1                 !fitting coefficient
      real(kind=fPrec) c_v2                 !fitting coefficient
      real(kind=fPrec) Dechan               !probability for dechanneling
      real(kind=fPrec) Lrefl, Srefl         !distance of the reflection point [m]
      real(kind=fPrec) Vcapt                !volume capture probability
      real(kind=fPrec) Chann                !channeling probability
      real(kind=fPrec) N_atom               !probability for entering
                                            !channeling near atomic planes
      real(kind=fPrec) Dxp                  !variation in angle
      real(kind=fPrec) xpcrit               !critical angle for curved crystal[rad]
      real(kind=fPrec) xpcrit0              !critical angle for str. crystal [rad]
      real(kind=fPrec) Rcrit                !critical curvature radius [m]
      real(kind=fPrec) ratio                !x=Rcurv/Rcrit
      real(kind=fPrec) Cry_length
      real(kind=fPrec) TLdech2              !tipical dechanneling length(1) [m]
      real(kind=fPrec) TLdech1              !tipical dechanneling length(2) [m]
      real(kind=fPrec) tdech, Ldech,Sdech   !angle, lenght, and S coordinate
                                            !of dechanneling point
      real(kind=fPrec) Rlength, Red_S       !reduced length/s coordinate
                                            !(in case of dechanneling)
      real(kind=fPrec) Am_length            !Amorphous length
      real(kind=fPrec) Length_xs, Length_ys !Amorphous length
      real(kind=fPrec) L_chan, tchan
      real(kind=fPrec) xp_rel               !xp-miscut angle in mrad

      real(kind=fPrec) const_dech,xpin,ypin
      real(kind=fPrec) alpha !Daniele: par for new chann prob
      real(kind=fPrec) Pvr !Daniele: prob for VR->AM transition

!
!
      NAM=1 !switch on/off the nuclear interaction (NAM) and the MCS (ZN)
      ZN=1

!-------------Daniele: dE/dX and dechanneling length calculation--------------------


!      write(*,*) "p entering routine", PC
      mom=PC*1.0d3 ! [GeV/c] -> [MeV/c]
      enr=(mom*mom+mp*mp)**0.5 ! [MeV]
      gammar=enr/mp
      betar=mom/enr
      bgr=betar*gammar

      Tmax=(2.0d0*me*bgr**2)/(1.0d0+2*gammar*me/mp+(me/mp)**2) ![MeV]

      plen=((rho(is2)*zatom(is2)/anuc(is2))**0.5)*28.816d-6 ![MeV]

      const_dech=(256.0/(9.0*(4.D0*DATAN(1.D0))**2))* &
     & (1.0/(log(2.0*me*gammar/(exenergy(is2)*c1e3))-1.0))*((aTF*dP)/(re*me))  ![m/MeV]

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
!      DESt = 0
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
        iProc = proc_out
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
        iProc = proc_AM
        CALL CALC_ION_LOSS_CRY(IS,is2,PC,AM_Length,DESt)
        CALL move_am(IS,is2,NAM,Am_Length,DESt,DLYi(IS2),dlri(is2),xp,yp,PC)
        x=x+xp*(s_length-s)
        y=y+yp*(s_length-s)
        GOTO 111
      ELSEIF ((x.GT.(C_xmax-Alayer)) .and. x.LT.(C_xmax)  ) THEN
        iProc = proc_AM
        CALL CALC_ION_LOSS_CRY(IS,is2,PC,s_length,DESt)
        CALL move_am(IS,is2,NAM,s_length,DESt,DLYi(IS2),DLRi(IS2), xp,yp,PC)
        WRITE(*,*)'Fix here!'
        GOTO 111
      END IF
!
! THIRD CASE: the p interacts with the crystal.
!. Define typical angles/probabilities for orientation 110
!
      xpcrit0 = (2.e-9*eUm(IS2)/PC)**0.5       ! critical angle (rad) for
                                              ! straight crystals
      Rcrit  = PC/(2.e-6*eUm(IS2))*AI(IS2)      ! critical curvature radius [m]
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
            iProc = proc_DC
            Dxp= Ldech/Rcurv             ! change angle from channeling [mrad]
            Sdech=Ldech*cos(miscut+0.5*Dxp)

            x  = x+ Ldech*(sin(0.5*Dxp+miscut))   ! trajectory at channeling exit
            xp = xp + Dxp + 2.0*(rndm4()-0.5)*xpcrit
            y= y + yp * Sdech
!            write(*,*) "Ldech", Ldech
!            write(*,*) "DESt", DESt
            CALL CALC_ION_LOSS_CRY(IS,is2,PC,Ldech,DESt)
            PC = PC - 0.5*DESt*Ldech          ! energy loss to ionization while in CH [GeV]

            x = x + 0.5*(s_length-Sdech)*xp
            y = y + 0.5*(s_length-Sdech)*yp

!            write(*,*) "s_length-Sdech", s_length-Sdech
            CALL CALC_ION_LOSS_CRY(IS,is2,PC,s_length-Sdech,DESt)
            CALL move_am(IS,is2,NAM,s_length-Sdech,DESt,DLYi(IS2),DLRi(IS2),xp,yp,PC)
           !next line new from sasha
            x = x + 0.5*(s_length-Sdech)*xp
            y = y + 0.5*(s_length-Sdech)*yp
          else
            iProc = proc_CH
            xpin=XP
            ypin=YP

!            write(*,*) "angles before entering CH subroutine (x,y):", xp, yp
            CALL MOVE_CH(IS,is2,NAM,L_chan,X,XP,YP,PC,Rcurv,Rcrit)  !daniele:check if a nuclear interaction happen while in CH
!            write(*,*) "angles after exiting CH subroutine (x,y):", xp, yp
!            stop

            if(iProc /= proc_CH) then             !daniele: if an nuclear interaction happened, move until the middle with initial xp,yp
            x = x + 0.5 * L_chan * xpin           !then propagate until the "crystal exit" with the new xp,yp
            y = y + 0.5 * L_chan * ypin           !accordingly with the rest of the code in "thin lens approx"
            x = x + 0.5 * L_chan * XP
            y = y + 0.5 * L_chan * YP
!            write(*,*) "Length", Length
            CALL CALC_ION_LOSS_CRY(IS,is2,PC,Length,DESt)
            PC = PC - DESt*Length       ! energy loss to ionization [GeV]
            else
            Dxp= L_chan/Rcurv + 0.5*RAN_GAUSS(1.0d0)*xpcrit ! change angle[rad]
            xp = Dxp
            !next line new from sasha
            x  = x+ L_chan*(sin(0.5*Dxp+miscut)) ! trajectory at channeling exit
!            xp = xp + Dxp + 2.0*(rndm4()-0.5)*xpcrit
            y = y + s_length * yp
            CALL CALC_ION_LOSS_CRY(IS,is2,PC,Length,DESt)
            PC = PC - 0.5*DESt*Length       ! energy loss to ionization [GeV]
            endif
          endif
        ELSE                                   !option 2: VR
                                               ! good for channeling
                                               ! but don't channel         (1-2)
          iProc = proc_VR
!          Dxp=0.5*(xp_rel/xpcrit+1)*Ang_avr
!          xp=xp+Dxp+Ang_rms*RAN_GAUSS(1.)
            !next line new from sasha
          xp=xp+0.45*(xp/xpcrit+1)*Ang_avr
          x = x + 0.5*s_length * xp
          y = y + 0.5*s_length * yp
          CALL CALC_ION_LOSS_CRY(IS,is2,PC,s_length,DESt)
          CALL move_am(IS,is2,NAM,s_length,DESt,DLYi(IS2),DLRi(IS2),xp ,yp,PC)
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
            iProc = proc_VR
            x = x + xp * Srefl
            y = y + yp * Srefl
            Dxp= Ang_avr
            xp = xp + Dxp + Ang_rms*RAN_GAUSS(1.0d0)
            x = x + 0.5* xp * (s_length - Srefl)
            y = y + 0.5* yp * (s_length - Srefl)
            CALL CALC_ION_LOSS_CRY(IS,is2,PC,s_length-Srefl,DESt)
            CALL move_am(IS,is2,NAM,s_length-Srefl,DESt,DLYi(IS2),DLRi(IS2),xp ,yp,PC)
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
              iProc = proc_DC
              Dxp= Ldech/Rcurv + 0.5*ran_gauss(1.0d0)*xpcrit
              x  = x+ Ldech*(sin(0.5*Dxp+xp))   ! trajectory at channeling exit
              y = y + Sdech * yp
              xp =  Dxp
              Red_S = s_length-Srefl -Sdech
              x = x + 0.5 * xp * Red_S
              y = y + 0.5 * yp * Red_S
!              write(*,*) "Srefl", Srefl
              CALL CALC_ION_LOSS_CRY(IS,is2,PC,Srefl,DESt)
              PC=PC - DESt * Srefl !Daniele: "added" energy loss before capture
!              write(*,*) "Sdech", Sdech
              CALL CALC_ION_LOSS_CRY(IS,is2,PC,Sdech,DESt)
              PC=PC - 0.5 * DESt * Sdech !Daniele: "added" energy loss while captured
              CALL CALC_ION_LOSS_CRY(IS,is2,PC,Red_S,DESt)
              CALL move_am(IS,is2,NAM,Red_S,DESt,DLYi(IS2),DLRi(IS2),xp,yp,PC)
              x = x + 0.5 * xp * Red_S
              y = y + 0.5 * yp * Red_S
            else
              iProc = proc_VC
              Rlength = Length-Lrefl
              tchan = Rlength / Rcurv
              Red_S=Rlength*cos(xp+0.5*tchan)
!              write(*,*) "Lrefl", Lrefl
              CALL CALC_ION_LOSS_CRY(IS,is2,PC,Lrefl,DESt)
              PC=PC - DESt*Lrefl  !Daniele: "added" energy loss before capture
              xpin=XP
              ypin=YP
              CALL MOVE_CH(IS,is2,NAM,Rlength,X,XP,YP,PC,Rcurv,Rcrit)  !daniele:check if a nuclear interaction happen while in CH

              if(iProc /= proc_VC) then             !daniele: if an nuclear interaction happened, move until the middle with initial xp,yp
              x = x + 0.5 * Rlength * xpin           !then propagate until the "crystal exit" with the new xp,yp
              y = y + 0.5 * Rlength * ypin           !accordingly with the rest of the code in "thin lens approx"
              x = x + 0.5 * Rlength * XP
              y = y + 0.5 * Rlength * YP
!              write(*,*) "Rlength", Rlength
              CALL CALC_ION_LOSS_CRY(IS,is2,PC,Rlength,DESt)
              PC=PC - DESt*Rlength
              else
              Dxp = (Length-Lrefl)/Rcurv
!              write(*,*) "Dxp", Dxp
              x  = x+ sin(0.5*Dxp+xp)*Rlength     ! trajectory at channeling exit
              y = y + red_S * yp
              xp =  Length/Rcurv + 0.5*ran_gauss(1.0d0)*xpcrit ! [mrad]
!              write(*,*) "xp at channeling exit", xp
!              write(*,*) "Rlength", Rlength
              CALL CALC_ION_LOSS_CRY(IS,is2,PC,Rlength,DESt)
              PC=PC - 0.5*DESt*Rlength  !Daniele: "added" energy loss once captured
              endif
            endif
          ENDIF
!.  case 3-3: move in amorphous substance (big input angles)---------------  MODIFIED FOR TRANSITION VRAM DANIELE
        else
           if(xp_rel .GT. L_chan/Rcurv+2.0*xpcrit .or. xp_rel .lt. -xpcrit) then
             iProc = proc_AM
             x = x + 0.5 * s_length * xp
             y = y + 0.5 * s_length * yp
            if(ZN .gt. 0) then
             CALL CALC_ION_LOSS_CRY(IS,is2,PC,s_length,DESt)
             CALL move_am(IS,is2,NAM,s_length,DESt,DLYi(IS2),DLRi(IS2), xp,yp,PC)
            endif
            x = x + 0.5 * s_length * xp
            y = y + 0.5 * s_length * yp
          else
!            Pvr=0.5*erf((xp_rel-(L_chan/Rcurv)-xpcrit)/(2.0*4.0*xpcrit*xpcrit)**0.5)+0.5
!            Pvr=0.2+(0.6/(2.0*xpcrit))*(xp_rel-(L_chan/Rcurv))
            Pvr=((xp_rel-(L_chan/Rcurv))/(2.0*xpcrit))
            if(rndm4() .gt. Pvr) then
            iProc = proc_TRVR
            x = x + xp * Srefl
            y = y + yp * Srefl
!            Dxp=(2.0*Ang_avr-Ang_rms)/(2.0*L_chan/Rcurv+2.0*xpcrit)*
!     +         xp_rel+L_chan/Rcurv*(2.0*Ang_avr-Ang_rms)/
!     +         (2.0*L_chan/Rcurv+2.0*xpcrit)-Ang_avr
            Dxp=-3.0*Ang_rms*xp_rel/(2.0*xpcrit)+Ang_avr+(3.0*Ang_rms*(L_chan/Rcurv)/(2.0*xpcrit))
!            write(*,*) xp_rel, Dxp, Ang_avr, Ang_rms
!            xp=xp+Dxp+Ang_rms*RAN_GAUSS(1.)
            xp=xp+Dxp
            x=x+0.5*xp*(s_length-Srefl)
            y=y+0.5*yp*(s_length-Srefl)
!            write(*,*) "s_length-Srefl", s_length-Srefl
            CALL CALC_ION_LOSS_CRY(IS,is2,PC,s_length-Srefl,DESt)
            CALL move_am(IS,is2,NAM,s_length-Srefl,DESt,DLYi(IS2),DLRi(IS2),xp ,yp,PC)
            x = x + 0.5 * xp * (s_length - Srefl)
            y = y + 0.5 * yp * (s_length - Srefl)
            else
            iProc = proc_TRAM
            x = x + xp * Srefl
            y = y + yp * Srefl
            Dxp=-1.0*(13.6/PC)*SQRT(s_length/DLRi(IS2))*0.001*xp_rel/    &
     &          (2.0*xpcrit)+(13.6/PC)*SQRT(s_length/DLRi(IS2))*0.001*   &
     &          (1.0+(L_chan/Rcurv)/(2.0*xpcrit))
            xp=xp+Dxp
            x=x+0.5*xp*(s_length-Srefl)
            y=y+0.5*yp*(s_length-Srefl)
!            write(*,*) "s_length-Srefl", s_length-Srefl
            CALL CALC_ION_LOSS_CRY(IS,is2,PC,s_length-Srefl,DESt)
            CALL move_am(IS,is2,NAM,s_length-Srefl,DESt,DLYi(IS2),DLRi(IS2),xp ,yp,PC)
            x = x + 0.5 * xp * (s_length - Srefl)
            y = y + 0.5 * yp * (s_length - Srefl)
            endif
          endif
        endif
       ENDIF

!       write(*,*) "p at the end of routine", PC

!      if (counter .eq. 0) then
111   write(833,*)'crystal parameters:\n Length:',Length, '\n Rcurv:'   &
     & , Rcurv ,'\n Critical Radius:', Rcrit, 'ratio',ratio             &
     &, '\n Critical angle for straight:',                              &
     & xpcrit0,'\n critical angle for curved crystal:', xpcrit,' \n Leng&
     &th:', Length, '\n xmax:', C_xmax, ' ymax:', ymax, '  C_orient: '  &
     &, C_orient                                                        &
     &, '\n Avg angle reflection:', Ang_avr, '\n full channeling angle: &
     &',(Length/Rcurv)
      END

! ================================================================================================ !
!  Subroutine for the calculazion of the energy loss by ionisation
! ================================================================================================ !
subroutine calc_ion_loss_cry(is,is2,pc,dz,enlo)

  use mod_ranlux
  use mod_funlux
  use floatPrecision
  use coll_crystal
  use coll_materials, only : zatom, exenergy, rho, anuc

  implicit none

  integer,          intent(in)  :: is
  integer,          intent(in)  :: is2
  real(kind=fPrec), intent(in)  :: pc
  real(kind=fPrec), intent(in)  :: dz
  real(kind=fPrec), intent(out) :: enlo

  real(kind=fPrec) thl,tt,cs_tail,prob_tail
  real(kind=fPrec), parameter :: k = 0.307075 ! Constant in front bethe-bloch [mev g^-1 cm^2]

  thl       = four*k*zatom(is2)*dz*c1e2*rho(is2)/(anuc(is2)*betar**2) ! [MeV]
  EnLo      = ((k*zatom(is2))/(anuc(is2)*betar**2))*( &
    half*log((two*me*bgr*bgr*Tmax)/(c1e6*exenergy(is2)**2)) - &
    betar**2-log(plen/(exenergy(is2)*c1e3))-log(bgr)+half)
  EnLo      = EnLo*rho(is2)*0.1*dz ! [GeV]
  Tt        = EnLo*c1e3+thl ! [MeV]
  cs_tail   = ((k*zatom(is2))/(anuc(is2)*betar**2))*((half*((one/Tt)-(one/Tmax))) - &
    (log(Tmax/Tt)*(betar**2)/(two*Tmax)) + ((Tmax-Tt)/(four*(gammar**2)*(mp**2))))
  prob_tail = cs_tail*rho(is2)*dz*c1e2

  if(rndm4() < prob_tail) then
    EnLo = ((k*zatom(is2))/(anuc(is2)*betar**2))*( &
      half*log((two*me*bgr*bgr*Tmax)/(c1e6*exenergy(is2)**2)) - &
      betar**2-log(plen/(exenergy(is2)*c1e3))-log(bgr)+half   + &
      (TMax**2)/(eight*(gammar**2)*(mp**2)))
    EnLo = EnLo*rho(is2)*0.1 ! [GeV/m]
  else
    EnLo = EnLo/dz ! [GeV/m]
  end if

end subroutine calc_ion_loss_cry

! ================================================================================================ !
!  Subroutine for the movement in the amorphous
! ================================================================================================ !
subroutine move_am(is,is2,nam,dz,dei,dly,dlr,xp,yp,pc)

  use mod_ranlux
  use mod_funlux
  use floatPrecision
  use coll_crystal
  use coll_materials, only : anuc, hcut, bnref, csref

  implicit none

  integer,          intent(in)    :: is
  integer,          intent(in)    :: is2
  integer,          intent(in)    :: nam
  real(kind=fPrec), intent(in)    :: dz
  real(kind=fPrec), intent(in)    :: dei
  real(kind=fPrec), intent(in)    :: dly
  real(kind=fPrec), intent(in)    :: dlr
  real(kind=fPrec), intent(inout) :: xp
  real(kind=fPrec), intent(inout) :: yp
  real(kind=fPrec), intent(inout) :: pc

  integer i,length_cry,ichoix
  real(kind=fPrec) t,xran_cry(1),bn(4),cs(0:5,1:4),freep(4),zlm,xp_in,yp_in,xm2,xln15s,tz,tx,tlow,  &
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
  mcurr_cry = is2
  thigh     = hcut(is2)
  call funlxp(ruth_cry,cgen_cry(1,is),tlow,thigh)

  ! Cross-section calculation
  ! freep: number of nucleons involved in single scattering
  freep(is) = freeco_cry * anuc(is2)**(one/three)

  ! Compute pp and pn el+single diff contributions to cross-section (both added : quasi-elastic or qel later)
  cs(3,is) = freep(is) * ppel
  cs(4,is) = freep(is) * ppsd

  ! Correct TOT-CSec for energy dependence of qel
  ! TOT CS is here without a Coulomb contribution
  cs(0,is) = csref(0,is2) + freep(is) * (pptot - pptref_cry)
  bn(is)   = bnref(is2) * cs(0,is) / csref(0,is2)

  ! Also correct inel-CS
  cs(1,is) = csref(1,is2) * cs(0,is) / csref(0,is2)

  ! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
  cs(2,is) = cs(0,is) - cs(1,is) - cs(3,is) - cs(4,is)
  cs(5,is) = csref(5,is2)
  ! Now add Coulomb
  cs(0,is) = cs(0,is) + cs(5,is)

  ! Calculate cumulative probability
  do i=1,4
    cprob_cry(i,is) = cprob_cry(i-1,is)+cs(i,is)/cs(0,is)
  end do

  ! Multiple Coulomb Scattering

  xp  = xp*c1e3
  yp  = yp*c1e3
  pc  = pc-dei*dz  ! Energy lost because of ionization process[GeV]

  dya   = (13.6/pc)*sqrt(dz/dlr) ! rms of coloumb scattering MCS (mrad)
  kxmcs = dya*ran_gauss(one)
  kymcs = dya*ran_gauss(one)

  xp = xp+kxmcs
  yp = yp+kymcs

  if(nam == 0) return ! Turn on/off nuclear interactions

  ! Can nuclear interaction happen?
  zlm = -collnt_cry(is)*log(rndm4())

  if(zlm < dz) then
    ! Choose nuclear interaction
    aran = rndm4()
    i=1
10  if(aran > cprob_cry(i,is)) then
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
      t     = -log(rndm4())/bn(is)

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
subroutine move_ch(is,is2,nam,dz,x,xp,yp,pc,r,rc)

  use crcoall
  use mod_ranlux
  use mod_funlux
  use floatPrecision
  use coll_common, only : coll_debug
  use coll_crystal
  use coll_materials, only : nmat, rho, anuc, hcut, bnref, csref

  implicit none

  integer,          intent(in)    :: is
  integer,          intent(in)    :: is2
  integer,          intent(in)    :: nam
  real(kind=fPrec), intent(in)    :: dz
  real(kind=fPrec), intent(inout) :: x
  real(kind=fPrec), intent(inout) :: xp
  real(kind=fPrec), intent(inout) :: yp
  real(kind=fPrec), intent(inout) :: pc
  real(kind=fPrec), intent(in)    :: r
  real(kind=fPrec), intent(in)    :: rc

  integer i,np,length_cry,ichoix
  real(kind=fPrec) t,xran_cry(1),bn(4),cs(0:5,1:4),freep(4),zlm,xp_in,yp_in,xminU,xm2,xln15s,x_min, &
    x_max,x_i,Umin,Ueff,tz,tx,tlow,thigh,teta,rho_min,rho_max,pv,pptot,ppsd,ppel,PC_in,nuc_cl_l,    &
    N_am,Et,ecmsq,Ec,csref_inel_rsc,csref_tot_rsc,bsd,bpp,aran,avrrho
  real(kind=fPrec), external :: ruth_cry

  if(coll_debug) then
    write(lout,"(2(a,1pe15.6))") "CRY> Angels before channeling: xp = ",xp," yp = ",yp
  end if

  xp_in = xp
  yp_in = yp
  pc_in = pc

  ! New treatment of scattering routine based on standard sixtrack routine

  ! Useful calculations for cross-section and event topology calculation
  ecmsq = 2 * 0.93828d0 * PC
  xln15s=log(0.15*ecmsq)

  ! New models, see Claudia's thesis
  pptot = 0.041084d0-0.0023302d0*log(ecmsq)+0.00031514d0*log(ecmsq)**2
  ppel  = (11.7d0-1.59d0*log(ecmsq)+0.134d0*log(ecmsq)**2)/1000
  ppsd  = (4.3d0+0.3d0*log(ecmsq))/1000
  bpp   = 7.156d0+1.439d0*log(sqrt(ecmsq))

  ! Distribution for Ruth. scatt.
  tlow      = tlcut_cry
  mcurr_cry = is2
  thigh     = hcut(is2)
  call funlxp(ruth_cry,cgen_cry(1,IS),tlow,thigh)

  ! Rescale the total and inelastic cross-section accordigly to the average density seen
  x_i = x
  np  = int(x_i/dp)    ! Calculate in which crystalline plane the particle enters
  x_i = x_i - Np*dP    ! Rescale the incoming x at the left crystalline plane
  x_i = x_i - (dP/two) ! Rescale the incoming x in the middle of crystalline planes

  pv   = pc*c1e9*pc*c1e9/sqrt(pc*c1e9*pc*c1e9 + 93828.0e6_fPrec*93828.0e6_fPrec) ! Calculate pv=P/E
  Ueff = eUm(is2)*(two*x_i/dp)*(two*x_i/dp)+pv*x_i/r                             ! Calculate effective potential
  Et   = pv*xp*xp/two+Ueff                                                       ! Calculate transverse energy
  Ec   = eUm(is2)*(one-rc/r)*(one-rc/r)                                          ! Calculate critical energy in bent crystals

  ! To avoid negative Et
  xminU = -dp*dp*pc*c1e9/(eight*eUm(is2)*r)
  Umin  = abs(eUm(is2)*(two*xminU/dp)*(two*xminU/dP)+pv*xminU/R)
  Et    = Et+Umin
  Ec    = Ec+Umin

  ! Calculate min e max of the trajectory between crystalline planes
  x_min = -(dP/two)*Rc/R-(dP/two)*sqrt(Et/Ec)
  x_Max = -(dP/two)*Rc/R+(dP/two)*sqrt(Et/Ec)

  ! Change ref. frame and go back with 0 on the crystalline plane on the left
  x_min = x_min - dp/two
  x_max = x_max - dp/two

  ! Calculate the "normal density" in m^-3
  N_am  = rho(is2)*6.022e23_fPrec*c1e6/anuc(is2)

  ! Calculate atomic density at min and max of the trajectory oscillation
  rho_max = n_am*dp/two*(erf(x_max/sqrt(two*u1*u1)) - erf((dP-x_Max)/sqrt(two*u1*u1)))
  rho_min = N_am*dP/two*(erf(x_min/sqrt(two*u1*u1)) - erf((dP-x_min)/sqrt(two*u1*u1)))

  ! "zero-approximation" of average nuclear density seen along the trajectory
  avrrho  = (rho_max-rho_min)/(x_max-x_min)
  avrrho  = two*avrrho/N_am

  csref_tot_rsc  = csref(0,is2)*avrrho ! Rescaled total ref cs
  csref_inel_rsc = csref(1,is2)*avrrho ! Rescaled inelastic ref cs

  ! Cross-section calculation
  freep(is) = freeco_cry * anuc(is2)**(one/three)

  ! compute pp and pn el+single diff contributions to cross-section (both added : quasi-elastic or qel later)
  cs(3,is) = freep(is) * ppel
  cs(4,is) = freep(is) * ppsd

  ! correct TOT-CSec for energy dependence of qel
  ! TOT CS is here without a Coulomb contribution
  cs(0,is) = csref_tot_rsc + freep(IS) * (pptot - pptref_cry)

  ! Also correct inel-CS
  if(csref_tot_rsc == zero) then
    cs(1,is) = zero
  else
    cs(1,is) = csref_inel_rsc * cs(0,is) / csref_tot_rsc
  end if

  ! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
  cs(2,is) = cs(0,is) - cs(1,is) - cs(3,is) - cs(4,is)
  cs(5,is) = csref(5,is2)
  ! Now add Coulomb
  cs(0,is) = cs(0,is) + cs(5,is)

  ! Calculate cumulative probability
  if(cs(0,is) == zero) then
    do i=1,4
      cprob_cry(i,is) = cprob_cry(i-1,is)
    end do
  else
    do i=1,4
      cprob_cry(i,is) = cprob_cry(i-1,is)+cs(i,is)/cs(0,is)
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
    nuc_cl_l = collnt_cry(is)/avrrho
  end if
  zlm = -nuc_cl_l*log(rndm4())

  ! write(889,*) x_i,pv,Ueff,Et,Ec,N_am,avrrho,csref_tot_rsc,csref_inel_rsc,nuc_cl_l

  if(zlm < dz) then
    ! Choose nuclear interaction
    aran = rndm4()
    i=1
10  if(aran > cprob_cry(i,is)) then
      i=i+1
      goto 10
    end if
    ichoix = i

    ! Do the interaction 
    select case(ichoix)
    case(1) ! deep inelastic, impinging p disappeared
      iProc = proc_ch_absorbed
    
    case(2) ! p-n elastic
      iProc  = proc_ch_pne
      bn(IS) = bnref(is2) * cs(0,IS) / csref_tot_rsc
      t      = -log(rndm4())/bn(IS)

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
        bsd = (106.d0 - 17.d0*xm2) *  bpp / 36.d0
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
