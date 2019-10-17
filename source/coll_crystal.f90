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
  use physical_constants

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
  real(kind=fPrec), parameter :: tlcut_cry = 0.0009982_fPrec
  real(kind=fPrec), save      :: cgen_cry(200,nmat)
  integer,          save      :: mcurr_cry

  real(kind=fPrec), save :: enr,mom,betar,gammar,bgr !Daniele: energy,momentum,beta relativistic, gamma relativistic
  real(kind=fPrec), save :: Tmax,plen !Daniele: maximum energy tranfer in single collision, plasma energy (see pdg)

  real(kind=fPrec), parameter :: re  = crade ! Electron radius [m]
  real(kind=fPrec), parameter :: me  = pmae  ! Electron mass [MeV/c^2]
  real(kind=fPrec), parameter :: mp  = pmap  ! Proton mass [MeV/c^2]

  real(kind=fPrec), parameter :: aTF = 0.194e-10_fPrec ! Screening function [m]
  real(kind=fPrec), parameter :: dP  = 1.920e-10_fPrec ! Distance between planes (110) [m]
  real(kind=fPrec), parameter :: u1  = 0.075e-10_fPrec ! Thermal vibrations amplitude

  ! pp cross-sections and parameters for energy dependence
  real(kind=fPrec), parameter :: pptref_cry = 0.040_fPrec
  real(kind=fPrec), parameter :: freeco_cry = 1.618_fPrec

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
  dlri(m)   =  0.0937_fPrec
  dlyi(m)   =  0.4652_fPrec
  ai(m)     =  0.96e-7_fPrec
  eUm(m)    = 21.34_fPrec
  collnt(m) =  0.3016_fPrec

  ! W
  m = collmat_getCollMatID("W")
  dlri(m)   =  0.0035_fPrec
  dlyi(m)   =  0.096_fPrec
  ai(m)     =  0.56e-7_fPrec
  eUm(m)    = 21.0_fPrec
  collnt(m) = zero

  ! C
  m = collmat_getCollMatID("C")
  dlri(m)   =  0.188_fPrec
  dlyi(m)   =  0.400_fPrec
  ai(m)     =  0.63e-7_fPrec
  eUm(m)    = 21.0_fPrec
  collnt(m) = zero

  ! Ge
  m = collmat_getCollMatID("Ge")
  dlri(m)   =  0.02302_fPrec
  dlyi(m)   =  0.2686_fPrec
  ai(m)     =  1.0e-7_fPrec
  eUm(m)    = 40.0_fPrec
  collnt(m) =  0.1632_fPrec

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
      use mod_common, only : napx
      use coll_common, only : cry_proc, xp_pencil0, yp_pencil0, x_pencil, y_pencil, pencil_dx, ipencil
      use floatPrecision
      use coll_crystal
      use mathlib_bouncer

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
         real(kind=fPrec)      cRot, sRot, cRRot, sRRot

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


  ! CRY ---------------------
  Rcurv = cdb_cryBend(icoll)
  Alayer = cdb_cryThick(icoll)
  C_xmax = cdb_cryXDim(icoll)
  C_ymax = cdb_cryYDim(icoll)
  C_orient = cdb_cryOrient(icoll)
  miscut = cdb_cryMisCut(icoll)
  Cry_bend =  Cry_length/Rcurv
  ! CRY ---------------------

  mat = cdb_cMaterialID(icoll)
  p0  = enom

  ! Initialise scattering processes
  call k2coll_scatin(p0)

  nhit   = 0
  fracab = zero
  mirror = one

  ! CRY ---------------------
  NEV = NP
  c_tilt0(1) = c_tilt(1)
  c_tilt0(2) = c_tilt(2)
  tiltangle  = c_tilt0(1)
  n_chan  = 0
  n_VR    = 0
  n_amorphous = 0
  ! CRY ---------------------

  cRot  = cos_mb(c_rotation)
  sRot  = sin_mb(c_rotation)
  cRRot = cos_mb(-c_rotation)
  sRRot = sin_mb(-c_rotation)

  do j=1,napx

    if(part_abs(j) /= 0 .and. part_abs_turn(j) /= 0) then
      ! Don't do scattering process for particles already absorbed
      cycle
    end if

    impact(j) = -one
    lint(j)   = -one
    indiv(j)  = -one

    ! CRY ---------------------
    if(iturn == 1) then
      bool_proc_old(j) = -1
    else
      if(bool_proc(j) /= -1) then
        bool_proc_old(j) = bool_proc(j)
      end if
    end if
    iProc = proc_out
    bool_proc(j) = -1
    cry_proc(j)  = -1
    ! CRY ---------------------

    ! CRY ---------------------
    nabs = 0
    s    = 0
    ! CRY ---------------------

    x   = x_in(j)
    xp  = xp_in(j)
    z   = y_in(j)
    zp  = yp_in(j)
    p   = p_in(j)

    ! CRY ---------------------
    x_temp   = zero
    x_int    = zero
    x_rot    = zero
    x_shift  = zero
    s_temp   = zero
    s_int    = zero
    s_rot    = zero
    s_shift  = zero
    xp_int   = zero
    xp_temp  = zero
    xp_rot   = zero
    xp_shift = zero
    shift    = zero
    tilt_int = zero
    dpop     = (p-p0)/p0

    s_in0(j)   = s_in(j)
  ! x_in0(j)   = x
    xp_in0(j)  = xp
    y_in0(j)   = z
    yp_in0(j)  = zp
    p_in0(j)   = p
    ! CRY ---------------------

    x  =  x_in(j)*cRot + sRot*y_in(j)
    z  =  y_in(j)*cRot - sRot*x_in(j)
    xp = xp_in(j)*cRot + sRot*yp_in(j)
    zp = yp_in(j)*cRot - sRot*xp_in(j)

    x  = (x - c_aperture/two) - c_offset

    ! Include collimator tilt
    if(tiltangle > zero) then
      xp = xp - tiltangle
    end if
    if(tiltangle < zero) then
      x  = x + sin_mb(tiltangle) * c_length
      xp = xp - tiltangle
    end if

    ! CRY ---------------------
    x_in0(j) = x
    ! CRY ---------------------

    ! CRY ---------------------
    ! Transform in the crystal rerference system
    ! 1st transformation: shift of the center of the reference frame
    if(Cry_tilt < zero) then
      S_shift = S
      shift   = Rcurv*(1-cos_mb(Cry_tilt))
      if(Cry_tilt < -Cry_bend) then
        shift = Rcurv*(cos_mb(-Cry_tilt) - cos_mb(Cry_bend-Cry_tilt))
      end if
      X_shift = X-shift
    else
      S_shift = S
      X_shift = X
    end if

    ! 2nd transformation: rotation
    S_rot  = X_shift*sin_mb(Cry_tilt)+S_shift*cos_mb(Cry_tilt)
    X_rot  = X_shift*cos_mb(Cry_tilt)-S_shift*sin_mb(Cry_tilt)
    XP_rot = XP - Cry_tilt

    ! 3rd transformation: drift to the new coordinate s=0
    XP = XP_rot
    X  = X_rot - XP_rot*S_rot
    Z  = Z - ZP*S_rot
    S  = zero

    ! Now check if the particle hit the crystal
    if(x >= zero .and. x < C_xmax) then ! Check that part. hit cry

      s_impact = s_in0(j) ! (for the first impact)
      call cryst(mat,x,xp,z,zp,p,cry_length,j)
      s   = Rcurv*sin_mb(cry_bend)
      zlm = Rcurv*sin_mb(cry_bend)
      if(iProc /= proc_out) then
        nhit = nhit + 1
        lhit(j) = ie
        lhit_turn(j) = iturn
        impact(j) = x_in0(j)
        indiv(j) = xp_in0(j)
      end if

    else

      if(x < zero) then ! Crystal hit from below
        XP_tangent = sqrt((-(two*x)*Rcurv + x**2)/(Rcurv**2))
        if(xp >= XP_tangent) then
          ! if it hits the crystal, calculate in which point and apply the
          ! transformation and drift to that point
          a_eq  = (one + xp**2)
          b_eq  = (two*xp)*(x-Rcurv)
          c_eq  = -(two*x)*Rcurv + x**2
          Delta = b_eq**2 - four*(a_eq*c_eq)
          S_int = (-b_eq - sqrt(Delta))/(two*a_eq)

          if(S_int < Rcurv*sin_mb(cry_bend)) then
            !  transform to a new ref system:shift and rotate
            x_int  = xp*s_int + x
            xp_int = xp
            z      = z + zp*s_int
            x      = zero
            s      = zero

            tilt_int = s_int/Rcurv
            xp       = xp-tilt_int

            call cryst(mat,x,xp,z,zp,p,(cry_length-(tilt_int*rcurv)),j)
            s   = Rcurv*sin_mb(cry_bend-tilt_int)
            zlm = Rcurv*sin_mb(cry_bend-tilt_int)
            if(iProc /= proc_out) then
              x_rot    = x_int
              s_rot    = s_int
              xp_rot   = xp_int
              s_shift  =  s_rot*cos_mb(-cry_tilt) + x_rot*sin_mb(-cry_tilt)
              x_shift  = -s_rot*sin_mb(-cry_tilt) + x_rot*cos_mb(-cry_tilt)
              xp_shift = xp_rot + cry_tilt

              if(Cry_tilt < zero) then
                s_impact  = s_shift
                x_in0(j)  = x_shift+shift
                xp_in0(j) = xp_shift
              else
                x_in0(j)  = x_shift
                s_impact  = s_shift
                xp_in0(j) = xp_shift
              end if

              nhit = nhit + 1
              lhit(j) = ie
              lhit_turn(j) = iturn
              impact(j) = x_in0(j)
              indiv(j) = xp_in0(j)
            end if

            ! un-rotate
            X_temp  = X
            S_temp  = S
            XP_temp = XP
            S       =  S_temp*cos_mb(-tilt_int) + X_temp*sin_mb(-tilt_int)
            X       = -S_temp*sin_mb(-tilt_int) + X_temp*cos_mb(-tilt_int)
            XP      = XP_temp + tilt_int

            ! 2nd: shift back the 2 axis
            X = X+X_int
            S = S+S_int

          else

            s = rcurv*sin_mb(cry_length/rcurv)
            x = x + s*xp
            z = z + s*zp

          end if

        else

          s = rcurv*sin_mb(cry_length/rcurv)
          x = x + s*xp
          z = z + s*zp

        end if

      else ! Crystal hit from above

        xp_tangent = asin_mb((rcurv*(one-cos_mb(cry_bend))-x)/sqrt(two*rcurv*(rcurv-x)*(one-cos_mb(cry_bend))+x**2))
        if(xp <= xp_tangent) then
          ! if it hits the crystal, calculate in which point and apply the
          ! transformation and drift to that point
          a_eq  = one + xp**2
          b_eq  = (two*xp)*(x - Rcurv)
          c_eq  = -(two*x)*Rcurv + x**2
          Delta = b_eq**2 - four*(a_eq*c_eq)
          S_int = (-b_eq-sqrt(Delta))/(two*a_eq)

          if(s_int < rcurv*sin_mb(cry_bend)) then
            !  transform to a new ref system:shift and rotate
            x_int  = xp*s_int + x
            xp_int = xp
            z      = z + zp*s_int
            x      = zero
            s      = zero

            tilt_int = s_int/rcurv
            xp       = xp-tilt_int

            call cryst(mat,x,xp,z,zp,p,(cry_length-(tilt_int*rcurv)),j)
            s   = Rcurv*sin_mb(cry_bend-tilt_int)
            zlm = Rcurv*sin_mb(cry_bend-tilt_int)
            if(iProc /= proc_out) then
              x_rot    = x_int
              s_rot    = s_int
              xp_rot   = xp_int
              s_shift  =  s_rot*cos_mb(-cry_tilt)+x_rot*sin_mb(-cry_tilt)
              x_shift  = -s_rot*sin_mb(-cry_tilt)+x_rot*cos_mb(-cry_tilt)
              xp_shift = xp_rot + cry_tilt

              if(Cry_tilt < zero) then
                s_impact  = s_shift
                x_in0(j)  = x_shift + shift
                xp_in0(j) = xp_shift
              else
                x_in0(j)  = x_shift
                s_impact  = s_shift
                xp_in0(j) = xp_shift
              end if

              nhit = nhit + 1
              lhit(j) = ie
              lhit_turn(j) = iturn
              impact(j) = x_in0(j)
              indiv(j) = xp_in0(j)
            end if

            ! un-rotate
            x_temp  = x
            s_temp  = s
            xp_temp = xp
            s  =  s_temp*cos_mb(-tilt_int)+x_temp*sin_mb(-tilt_int)
            x  = -s_temp*sin_mb(-tilt_int)+x_temp*cos_mb(-tilt_int)
            xp = xp_temp + tilt_int

            ! 2nd: shift back the 2 axis
            x = x+x_int
            s = s+s_int

          else

            s = rcurv*sin_mb(cry_length/rcurv)
            x = x + s*xp
            z = z + s*zp

          end if

        else

          s = rcurv*sin_mb(cry_length/rcurv)
          x = x+s*xp
          z = z+s*zp

        end if
      end if
    end if

    ! transform back from the crystal to the collimator reference system
    ! 1st: un-rotate the coordinates
    x_rot  = x
    s_rot  = s
    xp_rot = xp

    s_shift  =  s_rot*cos_mb(-cry_tilt)+x_rot*sin_mb(-cry_tilt)
    x_shift  = -s_rot*sin_mb(-cry_tilt)+x_rot*cos_mb(-cry_tilt)
    xp_shift = xp_rot + cry_tilt

    ! 2nd: shift back the reference frame
    if(cry_tilt < zero) then
      s  = s_shift
      x  = x_shift + shift
      xp = xp_shift
    else
      x  = x_shift
      s  = s_shift
      xp = xp_shift
    end if

    ! 3rd: shift to new S=Length position
    x = xp*(c_length-s) + x
    z = zp*(c_length-s) + z
    s = c_length
    ! CRY ---------------------

    ! CRY ---------------------
    nabs = 0
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
    ! CRY ---------------------

    ! Transform back to particle coordinates with opening and offset
    if(part_abs(j) == 0) then
      if(tiltangle > zero) then
        x  = x  + tiltangle*c_length
        xp = xp + tiltangle
      else if(tiltangle < zero) then
        x  = x + tiltangle*c_length
        xp = xp + tiltangle
        x  = x - sin_mb(tiltangle) * c_length
      end if

      ! Transform back to particle coordinates with opening and offset
      z00 = z
      x00 = x + mirror*c_offset
      x   = x + c_aperture/2 + mirror*c_offset

      ! Now mirror at the horizontal axis for negative X offset
      x   = mirror * x
      xp  = mirror * xp

      ! Last do rotation into collimator frame
      x_in(j)  =  x*cRRot +  z*sRRot
      y_in(j)  =  z*cRRot -  x*sRRot
      xp_in(j) = xp*cRRot + zp*sRRot
      yp_in(j) = zp*cRRot - xp*sRRot

      p_in(J) = P
      s_in(J) = s_in(J) + S

      if(nabs == 1) then
        fracab = fracab + 1
        part_abs(j) = ie
        part_abs_turn(j) = iturn
        lint(j) = zlm
        x = 99.99d-3
        z = 99.99d-3
      end if
    end if

    if(flagsec(j) == 0 .and. iProc /= proc_out) then
      flagsec(j) = 1
    end if

  end do

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
  use mathlib_bouncer

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

  real(kind=fPrec) const_dech,xpin,ypin,tchan,L_chan,mep

  real(kind=fPrec), parameter :: c_v1 =  1.7_fPrec ! Fitting coefficient
  real(kind=fPrec), parameter :: c_v2 = -1.5_fPrec ! Fitting coefficient

  nam = 1 ! Switch on/off the nuclear interaction (NAM) and the MCS (ZN)
  zn  = 1

  ! dE/dX and dechanneling length calculation
  mom    = pc*c1e3                ! [GeV]
  enr    = sqrt(mom**2 + pmap**2) ! [MeV]
  gammar = enr/pmap
  betar  = mom/enr
  bgr    = betar*gammar
  mep    = pmae/pmap  ! Electron/proton

  tmax = (two*pmae*bgr**2)/(one + two*gammar*mep + mep**2)  ! [MeV]
  plen = sqrt((rho(is)*zatom(is))/anuc(is))*28.816e-6_fPrec ! [MeV]

  const_dech = ((256.0_fPrec/(nine*pi**2)) * &
    (one/(log_mb(((two*pmae)*gammar)/(exenergy(is)*c1e3)) - one))) * ((aTF*dP)/(re*pmae)) ! [m/MeV]
  const_dech = const_dech*c1e3 ! [m/GeV]

  s        = zero
  s_length = Rcurv*(sin_mb(length/Rcurv))
  L_chan   = length

  if(miscut < zero .and. x > zero .and. x < -length*tan_mb(miscut)) then
    L_chan = -x/sin_mb(miscut)
  end if

  tchan  = L_chan/Rcurv
  xp_rel = xp - miscut

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
    b_eq  = (two*x)*xp - (two*xp)*Rcurv
    c_eq  = x**2 - (two*x)*Rcurv
    Delta = b_eq**2 - (four*a_eq)*c_eq
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

  else if(x > c_xmax-alayer .and. x < c_xmax) then
    iProc = proc_AM
    call calc_ion_loss_cry(is,pc,s_length,dest)
    call move_am(is,nam,s_length,dest,dlyi(is),dlri(is),xp,yp,pc)
    return

  end if

  ! THIRD CASE: the p interacts with the crystal.
  ! Define typical angles/probabilities for orientation 110
  xpcrit0 = sqrt((c2m9*eUm(is))/pc)    ! Critical angle (rad) for straight crystals
  Rcrit   = (pc/(c2m6*eUm(is)))*ai(is) ! Critical curvature radius [m]

  ! If R>Rcritical=>no channeling is possible (ratio<1)
  ratio  = Rcurv/Rcrit
  xpcrit = (xpcrit0*(Rcurv-Rcrit))/Rcurv ! Critical angle for curved crystal

  if(ratio <= one) then ! no possibile channeling
    Ang_rms = ((c_v1*0.42_fPrec)*xpcrit0)*sin_mb(1.4_fPrec*ratio) ! RMS scattering
    Ang_avr = ((c_v2*xpcrit0)*c5m2)*ratio                         ! Average angle reflection
    Vcapt   = zero                                                ! Probability of VC

  else if(ratio <= three) then ! Strongly bent crystal
    Ang_rms = ((c_v1*0.42_fPrec)*xpcrit0)*sin_mb(0.4713_fPrec*ratio + 0.85_fPrec) ! RMS scattering
    Ang_avr = (c_v2*xpcrit0)*(0.1972_fPrec*ratio - 0.1472_fPrec)                  ! Average angle reflection
    Vcapt   = 7.0e-4_fPrec*(ratio - 0.7_fPrec)/pc**c2m1                           ! Correction by sasha drozdin/armen
    ! K=0.0007 is taken based on simulations using CATCH.f (V.Biryukov)

  else ! Rcry >> Rcrit
    Ang_rms = (c_v1*xpcrit0)*(one/ratio)                ! RMS scattering
    Ang_avr = (c_v2*xpcrit0)*(one - 1.6667_fPrec/ratio) ! Average angle for VR
    Vcapt   = 7.0e-4_fPrec*(ratio - 0.7_fPrec)/pc**c2m1 ! Probability for VC correction by sasha drozdin/armen
    ! K=0.0007 is taken based on simulations using CATCH.f (V.Biryukov)

  end if

  if(C_orient == 2) then
    Ang_avr = Ang_avr*0.93_fPrec
    Ang_rms = Ang_rms*1.05_fPrec
    xpcrit  = xpcrit*0.98_fPrec
  end if

  if(abs(xp_rel) < xpcrit) then
    alpha  = xp_rel/xpcrit
    Chann  = sqrt(0.9_fPrec*(one - alpha**2))*sqrt(one-(one/ratio)) ! Saturation at 95%
    N_atom = c1m1

    ! if they can channel: 2 options
    if(rndm4() <= chann) then ! option 1:channeling

      TLdech1 = (const_dech*pc)*(one-one/ratio)**2 ! Updated calculate typical dech. length(m)
      if(rndm4() <= n_atom) then
        TLdech1 = ((const_dech/c2e2)*pc)*(one-one/ratio)**2  ! Updated dechanneling length (m)
      end if

      Dechan = -log_mb(rndm4()) ! Probability of dechanneling
      Ldech  = TLdech1*Dechan   ! Actual dechan. length

      ! careful: the dechanneling lentgh is along the trajectory
      ! of the particle -not along the longitudinal coordinate...
      if(ldech < l_chan) then
        iProc = proc_DC
        Dxp   = Ldech/Rcurv ! Change angle from channeling [mrad]
        Sdech = Ldech*cos_mb(miscut + half*Dxp)
        x     = x  + Ldech*(sin_mb(half*Dxp+miscut)) ! Trajectory at channeling exit
        xp    = xp + Dxp + (two*(rndm4()-half))*xpcrit
        y     = y  + yp * Sdech

        call calc_ion_loss_cry(is,pc,ldech,dest)
        pc = pc - half*dest*Ldech ! Energy loss to ionization while in CH [GeV]
        x  = x  + (half*(s_length-Sdech))*xp
        y  = y  + (half*(s_length-Sdech))*yp

        call calc_ion_loss_cry(is,pc,s_length-sdech,dest)
        call move_am(is,nam,s_length-sdech,dest,dlyi(is),dlri(is),xp,yp,pc)
        x = x + (half*(s_length-Sdech))*xp
        y = y + (half*(s_length-Sdech))*yp
      else
        iProc = proc_CH
        xpin  = XP
        ypin  = YP

        call move_ch(is,nam,l_chan,x,xp,yp,pc,rcurv,rcrit) ! check if a nuclear interaction happen while in CH
        if(iProc /= proc_CH) then
          ! if an nuclear interaction happened, move until the middle with initial xp,yp then
          ! propagate until the "crystal exit" with the new xp,yp accordingly with the rest
          ! of the code in "thin lens approx"
          x = x + (half*L_chan)*xpin
          y = y + (half*L_chan)*ypin
          x = x + (half*L_chan)*XP
          y = y + (half*L_chan)*YP

          call calc_ion_loss_cry(is,pc,length,dest)
          pc = pc - dest*length ! energy loss to ionization [GeV]
        else
          Dxp = L_chan/Rcurv + (half*ran_gauss(one))*xpcrit ! Change angle[rad]
          xp  = Dxp
          x   = x + L_chan*(sin(half*Dxp+miscut)) ! Trajectory at channeling exit
          y   = y + s_length * yp

          call calc_ion_loss_cry(is,pc,length,dest)
          pc = pc - (half*dest)*length ! energy loss to ionization [GeV]
        end if
      end if

    else ! Option 2: VR

      ! good for channeling but don't channel (1-2)
      iProc = proc_VR

      xp = xp + (0.45_fPrec*(xp/xpcrit + one))*Ang_avr
      x  = x  + (half*s_length)*xp
      y  = y  + (half*s_length)*yp

      call calc_ion_loss_cry(is,pc,s_length,dest)
      call move_am(is,nam,s_length,dest,dlyi(is),dlri(is),xp,yp,pc)

      x = x + (half*s_length)*xp
      y = y + (half*s_length)*yp

    end if

  else ! case 3-2: no good for channeling. check if the  can VR

    Lrefl = xp_rel*Rcurv ! Distance of refl. point [m]
    Srefl = sin_mb(xp_rel/two + miscut)*Lrefl

    if(Lrefl > zero .and. Lrefl < Length) then ! VR point inside

      ! 2 options: volume capture and volume reflection

      if(rndm4() > Vcapt .or. ZN == zero) then ! Option 1: VR

        iProc = proc_VR
        x     = x + xp*Srefl
        y     = y + yp*Srefl
        Dxp   = Ang_avr
        xp    = xp + Dxp + Ang_rms*ran_gauss(one)
        x     = x  + (half*xp)*(s_length - Srefl)
        y     = y  + (half*yp)*(s_length - Srefl)

        call calc_ion_loss_cry(is,pc,s_length-srefl,dest)
        call move_am(is,nam,s_length-srefl,dest,dlyi(is),dlri(is),xp,yp,pc)
        x = x + (half*xp)*(s_length - Srefl)
        y = y + (half*yp)*(s_length - Srefl)

      else ! Option 2: VC

        x = x + xp*Srefl
        y = y + yp*Srefl

        TLdech2 = (const_dech/c1e1)*pc*(one-one/ratio)**2          ! Updated typical dechanneling length(m)
        Ldech   = TLdech2*(sqrt(c1m2 - log_mb(rndm4())) - c1m1)**2 ! Updated DC length
        tdech   = Ldech/Rcurv
        Sdech   = Ldech*cos_mb(xp + half*tdech)

        if(Ldech < Length-Lrefl) then

          iProc = proc_DC
          Dxp   = Ldech/Rcurv + (half*ran_gauss(one))*xpcrit
          x     = x + Ldech*(sin_mb(half*Dxp+xp)) ! Trajectory at channeling exit
          y     = y + Sdech*yp
          xp    =  Dxp
          Red_S = (s_length - Srefl) - Sdech
          x     = x + (half*xp)*Red_S
          y     = y + (half*yp)*Red_S

          call calc_ion_loss_cry(is,pc,srefl,dest)
          pc = pc - dest*Srefl ! "added" energy loss before capture

          call calc_ion_loss_cry(is,pc,sdech,dest)
          pc = pc - (half*dest)*Sdech ! "added" energy loss while captured

          call calc_ion_loss_cry(is,pc,red_s,dest)
          call move_am(is,nam,red_s,dest,dlyi(is),dlri(is),xp,yp,pc)
          x = x + (half*xp)*Red_S
          y = y + (half*yp)*Red_S

        else

          iProc   = proc_VC
          Rlength = Length - Lrefl
          tchan   = Rlength/Rcurv
          Red_S   = Rlength*cos_mb(xp + half*tchan)

          call calc_ion_loss_cry(is,pc,lrefl,dest)
          pc   = pc - dest*Lrefl ! "added" energy loss before capture
          xpin = xp
          ypin = yp

          call move_ch(is,nam,rlength,x,xp,yp,pc,rcurv,rcrit) ! check if a nuclear interaction happen while in ch
          if(iProc /= proc_VC) then
            ! if an nuclear interaction happened, move until the middle with initial xp,yp then propagate until
            ! the "crystal exit" with the new xp,yp accordingly with the rest of the code in "thin lens approx"
            x = x + (half*Rlength)*xpin
            y = y + (half*Rlength)*ypin
            x = x + (half*Rlength)*XP
            y = y + (half*Rlength)*YP

            call calc_ion_loss_cry(is,pc,rlength,dest)
            pc = pc - dest*Rlength
          else
            Dxp = (Length-Lrefl)/Rcurv
            x   = x + sin_mb(half*Dxp+xp)*Rlength ! Trajectory at channeling exit
            y   = y + red_S*yp
            xp  = Length/Rcurv + (half*ran_gauss(one))*xpcrit ! [mrad]

            call calc_ion_loss_cry(is,pc,rlength,dest)
            pc = pc - (half*dest)*Rlength  ! "added" energy loss once captured
          end if
        end if
      end if

    else

      ! Case 3-3: move in amorphous substance (big input angles)
      ! Modified for transition vram daniele
      if(xp_rel > L_chan/Rcurv + two*xpcrit .or. xp_rel < -xpcrit) then
        iProc = proc_AM
        x     = x + (half*s_length)*xp
        y     = y + (half*s_length)*yp
        if(zn > zero) then
          call calc_ion_loss_cry(is,pc,s_length,dest)
          call move_am(is,nam,s_length,dest,dlyi(is),dlri(is),xp,yp,pc)
        end if
        x = x + (half*s_length)*xp
        y = y + (half*s_length)*yp
      else
        Pvr = (xp_rel-(L_chan/Rcurv))/(two*xpcrit)
        if(rndm4() > Pvr) then
          iProc = proc_TRVR
          x     = x + xp*Srefl
          y     = y + yp*Srefl

          Dxp = (((-three*Ang_rms)*xp_rel)/(two*xpcrit) + Ang_avr) + ((three*Ang_rms)*(L_chan/Rcurv))/(two*xpcrit)
          xp  = xp + Dxp
          x   = x + (half*xp)*(s_length-Srefl)
          y   = y + (half*yp)*(s_length-Srefl)

          call calc_ion_loss_cry(is,pc,s_length-srefl,dest)
          call move_am(is,nam,s_length-srefl,dest,dlyi(is),dlri(is),xp,yp,pc)
          x = x + (half*xp)*(s_length - Srefl)
          y = y + (half*yp)*(s_length - Srefl)
        else
          iProc = proc_TRAM
          x = x + xp*Srefl
          y = y + yp*Srefl
          Dxp = ((((-one*(13.6_fPrec/pc))*sqrt(s_length/dlri(is)))*c1m3)*xp_rel)/(two*xpcrit) + &
            (((13.6_fPrec/pc)*sqrt(s_length/DLRi(is)))*c1m3)*(one+(L_chan/Rcurv)/(two*xpcrit))
          xp = xp+Dxp
          x  = x + (half*xp)*(s_length-Srefl)
          y  = y + (half*yp)*(s_length-Srefl)

          call calc_ion_loss_cry(is,pc,s_length-srefl,dest)
          call move_am(is,nam,s_length-srefl,dest,dlyi(is),dlri(is),xp,yp,pc)
          x = x + (half*xp)*(s_length - Srefl)
          y = y + (half*yp)*(s_length - Srefl)
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
  use mathlib_bouncer

  implicit none

  integer,          intent(in)  :: is
  real(kind=fPrec), intent(in)  :: pc
  real(kind=fPrec), intent(in)  :: dz
  real(kind=fPrec), intent(out) :: EnLo

  real(kind=fPrec) thl,tt,cs_tail,prob_tail
  real(kind=fPrec), parameter :: k = 0.307075_fPrec ! Constant in front bethe-bloch [mev g^-1 cm^2]

  thl       = (((((four*k)*zatom(is))*dz)*c1e2)*rho(is))/(anuc(is)*betar**2) ! [MeV]
  EnLo      = ((k*zatom(is))/(anuc(is)*betar**2)) * ( &
    half*log_mb(((((two*pmae)*bgr)*bgr)*Tmax)/(c1e6*exenergy(is)**2)) - &
    betar**2 - log_mb(plen/(exenergy(is)*c1e3)) - log_mb(bgr) + half    &
  )
  EnLo      = ((EnLo*rho(is))*c1m1)*dz ! [GeV]
  Tt        = (EnLo*c1e3)+thl          ! [MeV]

  cs_tail   = ((k*zatom(is))/(anuc(is)*betar**2)) * ((half*((one/Tt)-(one/Tmax))) - &
    (log_mb(Tmax/Tt)*(betar**2)/(two*Tmax)) + ((Tmax-Tt)/((four*(gammar**2))*(pmap**2))))
  prob_tail = ((cs_tail*rho(is))*dz)*c1e2

  if(rndm4() < prob_tail) then
    EnLo = ((k*zatom(is))/(anuc(is)*betar**2)) * ( &
      half*log_mb((two*pmae*bgr*bgr*Tmax)/(c1e6*exenergy(is)**2)) -      &
      betar**2 - log_mb(plen/(exenergy(is)*c1e3)) - log_mb(bgr) + half + &
      TMax**2/(eight*(gammar**2)*(pmap**2)) &
    )
    EnLo = (EnLo*rho(is))*c1m1 ! [GeV/m]
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
  use mathlib_bouncer

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
  ecmsq  = ((two*pmap)*c1m3)*pc
  xln15s = log_mb(0.15_fPrec*ecmsq)

  ! New models, see Claudia's thesis
  pptot = (0.041084_fPrec - 0.0023302_fPrec*log_mb(ecmsq)) + 0.00031514_fPrec*log_mb(ecmsq)**2
  ppel  = (11.7_fPrec - 1.59_fPrec*log_mb(ecmsq) + 0.134_fPrec*log_mb(ecmsq)**2)/c1e3
  ppsd  = (4.3_fPrec + 0.3_fPrec*log_mb(ecmsq))/c1e3
  bpp   = 7.156_fPrec + 1.439_fPrec*log_mb(sqrt(ecmsq))

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
  bn    = (bnref(is)*cs(0))/csref(0,is)

  ! Also correct inel-CS
  cs(1) = (csref(1,is)*cs(0))/csref(0,is)

  ! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
  cs(2) = ((cs(0) - cs(1)) - cs(3)) - cs(4)
  cs(5) = csref(5,is)

  ! Now add Coulomb
  cs(0) = cs(0) + cs(5)

  ! Calculate cumulative probability
  cprob(:) = zero
  cprob(5) = one
  do i=1,4
    cprob(i) = cprob(i-1) + cs(i)/cs(0)
  end do

  ! Multiple Coulomb Scattering
  xp  = xp*c1e3
  yp  = yp*c1e3
  pc  = pc - dei*dz ! Energy lost because of ionization process[GeV]

  dya   = (13.6_fPrec/pc)*sqrt(dz/dlr) ! RMS of coloumb scattering MCS (mrad)
  kxmcs = dya*ran_gauss(one)
  kymcs = dya*ran_gauss(one)

  xp = xp+kxmcs
  yp = yp+kymcs

  if(nam == 0) return ! Turn on/off nuclear interactions

  ! Can nuclear interaction happen?
  zlm = -collnt(is)*log_mb(rndm4())

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
      t     = -log_mb(rndm4())/bn

    case(3) ! p-p elastic
      iProc = proc_ppe
      t     = -log_mb(rndm4())/bpp

    case(4) ! Single diffractive
      iProc = proc_diff
      xm2   = exp_mb(rndm4()*xln15s)
      pc    = pc*(one - xm2/ecmsq)
      if(xm2 < two) then
        bsd = two*bpp
      else if(xm2 >= two .and. xm2 <= five) then
        bsd = ((106.0_fPrec - 17.0_fPrec*xm2)*bpp)/36.0_fPrec
      else if(xm2 > five) then
        bsd = 7.0_fPrec*bpp/12.0_fPrec
      end if
      t = -log_mb(rndm4())/bsd

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

    tx = (teta*ran_gauss(one))*c1e3
    tz = (teta*ran_gauss(one))*c1e3

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
  use mathlib_bouncer

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
  ecmsq  = ((two*pmap)*c1m3)*pc
  xln15s = log_mb(0.15_fPrec*ecmsq)

  ! New models, see Claudia's thesis
  pptot = (0.041084_fPrec - 0.0023302_fPrec*log_mb(ecmsq)) + 0.00031514_fPrec*log_mb(ecmsq)**2
  ppel  = (11.7_fPrec - 1.59_fPrec*log_mb(ecmsq) + 0.134_fPrec*log_mb(ecmsq)**2)/c1e3
  ppsd  = (4.3_fPrec + 0.3_fPrec*log_mb(ecmsq))/c1e3
  bpp   = 7.156_fPrec + 1.439_fPrec*log_mb(sqrt(ecmsq))

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

  pv   = pc**2/sqrt(pc**2 + (pmap*c1m3)**2)*c1e9          ! Calculate pv=P/E
  Ueff = eUm(is)*((two*x_i)/dp)*((two*x_i)/dp) + pv*x_i/r ! Calculate effective potential
  Et   = (pv*xp**2)/two + Ueff                            ! Calculate transverse energy
  Ec   = (eUm(is)*(one-rc/r))*(one-rc/r)                  ! Calculate critical energy in bent crystals

  ! To avoid negative Et
  xminU = ((-dp**2*pc)*c1e9)/(eight*eUm(is)*r)
  Umin  = abs((eUm(is)*((two*xminU)/dp))*((two*xminU)/dP) + pv*xminU/R)
  Et    = Et + Umin
  Ec    = Ec + Umin

  ! Calculate min e max of the trajectory between crystalline planes
  x_min = (-(dP/two)*Rc)/R - (dP/two)*sqrt(Et/Ec)
  x_Max = (-(dP/two)*Rc)/R + (dP/two)*sqrt(Et/Ec)

  ! Change ref. frame and go back with 0 on the crystalline plane on the left
  x_min = x_min - dp/two
  x_max = x_max - dp/two

  ! Calculate the "normal density" in m^-3
  N_am  = ((rho(is)*6.022e23_fPrec)*c1e6)/anuc(is)

  ! Calculate atomic density at min and max of the trajectory oscillation
  rho_max = ((N_am*dp)/two)*(erf(x_max/sqrt(two*u1**2)) - erf((dP-x_Max)/sqrt(two*u1**2)))
  rho_min = ((N_am*dP)/two)*(erf(x_min/sqrt(two*u1**2)) - erf((dP-x_min)/sqrt(two*u1**2)))

  ! "zero-approximation" of average nuclear density seen along the trajectory
  avrrho  = (rho_max - rho_min)/(x_max - x_min)
  avrrho  = (two*avrrho)/N_am

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
    cs(1) = (csref_inel_rsc*cs(0))/csref_tot_rsc
  end if

  ! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
  cs(2) = ((cs(0) - cs(1)) - cs(3)) - cs(4)
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
  zlm = -nuc_cl_l*log_mb(rndm4())

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
      bn    = (bnref(is)*cs(0))/csref_tot_rsc
      t     = -log_mb(rndm4())/bn

    case(3) ! p-p elastic
      iProc = proc_ch_ppe
      t     = -log_mb(rndm4())/bpp

    case(4) ! Single diffractive
      iProc = proc_ch_diff
      xm2   = exp_mb(rndm4()*xln15s)
      pc    = pc*(one - xm2/ecmsq)
      if(xm2 < two) then
        bsd = two*bpp
      else if(xm2 >= two .and. xm2 <= five) then
        bsd = ((106.0_fPrec - 17.0_fPrec*xm2)*bpp)/36.0_fPrec
      else if(xm2 > five) then
        bsd = (seven*bpp)/12.0_fPrec
      end if
      t = -log_mb(rndm4())/bsd

    case(5)
      iProc      = proc_ch_ruth
      length_cry = 1
      call funlux(cgen_cry(1,is),xran_cry,length_cry)
      t = xran_cry(1)

    end select

    ! Calculate the related kick -----------
    if(ichoix == 4) then
      teta = sqrt(t)/pc_in ! DIFF has changed PC!!!
    else
      teta = sqrt(t)/pc
    end if

    tx = (teta*ran_gauss(one))*c1e3
    tz = (teta*ran_gauss(one))*c1e3

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
  use mathlib_bouncer

  implicit none

  real(kind=fPrec) ruth_cry,t_cry
  real(kind=fPrec), parameter :: cnorm  = 2.607e-4_fPrec
  real(kind=fPrec), parameter :: cnform = 0.8561e3_fPrec

  ruth_cry = (cnorm*exp_mb((-t_cry*cnform)*emr(mcurr_cry)**2))*(zatom(mcurr_cry)/t_cry)**2

end function ruth_cry
