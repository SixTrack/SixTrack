! ================================================================================================ !
!
!  Crystal Collimation Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  Written by: Igor Yazynin, Valentina Previtali and Daniele Mirarchi, BE-ABP-HSS
!  Re-written for SixTrack 5 by: Marco D'Andrea and Veronica K. Berglyd Olsen, BE-ABP-HSS (2019)
!
!  Last modified: 2021-12-03
!
! ================================================================================================ !
module coll_crystal

  use floatPrecision
  use numerical_constants
  use coll_materials, only : nmat

  implicit none

  integer,          private, save :: iProc       = 0
  integer,          private, save :: n_chan      = 0
  integer,          private, save :: n_VR        = 0
  integer,          private, save :: n_amorphous = 0

  ! Shared settings for the currently active crystal
  integer,          private, save :: c_orient   = 0    ! Crystal orientation [0-2]
  real(kind=fPrec), private, save :: c_rcurv    = zero ! Crystal geometrical parameters [m]
  real(kind=fPrec), private, save :: c_xmax     = zero ! Crystal geometrical parameters [m]
  real(kind=fPrec), private, save :: c_ymax     = zero ! Crystal geometrical parameters [m]
  real(kind=fPrec), private, save :: c_alayer   = zero ! Crystal amorphous layer [mm]
  real(kind=fPrec), private, save :: c_miscut   = zero ! Crystal miscut angle in rad
  real(kind=fPrec), private, save :: c_cpTilt   = zero ! Cosine of positive crystal tilt
  real(kind=fPrec), private, save :: c_spTilt   = zero ! Sine of positive crystal tilt
  real(kind=fPrec), private, save :: c_cnTilt   = zero ! Cosine of negative crystal tilt
  real(kind=fPrec), private, save :: c_snTilt   = zero ! Sine of negative crystal tilt
  real(kind=fPrec), private, save :: c_cBend    = zero ! Cosine of crystal bend
  real(kind=fPrec), private, save :: c_sBend    = zero ! Sine of crystal bend
  real(kind=fPrec), private, save :: cry_tilt   = zero ! Crystal tilt angle in rad
  real(kind=fPrec), private, save :: cry_length = zero ! Crystal length [m]
  real(kind=fPrec), private, save :: cry_bend   = zero ! Crystal bending angle in rad

  ! Rutherford Scatter
  real(kind=fPrec), parameter     :: tlcut_cry = 0.0009982_fPrec
  real(kind=fPrec), private, save :: cgen_cry(200,nmat)
  integer,          private, save :: mcurr_cry

  real(kind=fPrec), private, save :: enr
  real(kind=fPrec), private, save :: mom
  real(kind=fPrec), private, save :: betar
  real(kind=fPrec), private, save :: gammar
  real(kind=fPrec), private, save :: bgr
  real(kind=fPrec), private, save :: tmax
  real(kind=fPrec), private, save :: plen

  real(kind=fPrec), parameter :: aTF = 0.194e-10_fPrec ! Screening function [m]
  real(kind=fPrec), parameter :: dP  = 1.920e-10_fPrec ! Distance between planes (110) [m]
  real(kind=fPrec), parameter :: u1  = 0.075e-10_fPrec ! Thermal vibrations amplitude

  ! pp cross-sections and parameters for energy dependence
  real(kind=fPrec), parameter :: pptref_cry = 0.040_fPrec
  real(kind=fPrec), parameter :: freeco_cry = 1.618_fPrec

  ! Crystal Specific Material Arrays
  logical,          private, save :: validMat(nmat) = .false. ! True for materials the crystal module supports
  real(kind=fPrec), private, save :: dlri(nmat)     = zero
  real(kind=fPrec), private, save :: dlyi(nmat)     = zero
  real(kind=fPrec), private, save :: ai(nmat)       = zero
  real(kind=fPrec), private, save :: eUm(nmat)      = zero
  real(kind=fPrec), private, save :: collnt(nmat)   = zero    ! Nuclear Collision length [m]

  ! Processes
  integer, parameter :: proc_out         =  -1     ! Crystal not hit
  integer, parameter :: proc_AM          =   1     ! Amorphous
  integer, parameter :: proc_VR          =   2     ! Volume reflection
  integer, parameter :: proc_CH          =   3     ! Channeling
  integer, parameter :: proc_VC          =   4     ! Volume capture
  integer, parameter :: proc_absorbed    =   5     ! Absorption
  integer, parameter :: proc_DC          =   6     ! Dechanneling
  integer, parameter :: proc_pne         =   7     ! Proton-neutron elastic interaction
  integer, parameter :: proc_ppe         =   8     ! Proton-proton elastic interaction
  integer, parameter :: proc_diff        =   9     ! Single diffractive
  integer, parameter :: proc_ruth        =  10     ! Rutherford scattering
  integer, parameter :: proc_ch_absorbed =  15     ! Channeling followed by absorption
  integer, parameter :: proc_ch_pne      =  17     ! Channeling followed by proton-neutron elastic interaction
  integer, parameter :: proc_ch_ppe      =  18     ! Channeling followed by proton-proton elastic interaction
  integer, parameter :: proc_ch_diff     =  19     ! Channeling followed by single diffractive
  integer, parameter :: proc_ch_ruth     =  20     ! Channeling followed by Rutherford scattering
  integer, parameter :: proc_TRVR        = 100     ! Volume reflection in VR-AM transition region
  integer, parameter :: proc_TRAM        = 101     ! Amorphous in VR-AM transition region

contains

! ================================================================================================ !
!  Initialise the crystal module
! ================================================================================================ !
subroutine cry_init

  use coll_materials

  integer m

  ! dlri : Radiation length(m), updated from PDG for Si
  ! dlyi : Nuclear length(m)
  ! ai   : Si110 1/2 interplan. dist. mm, Ge taken from A. Fomin, Si from initial implementation
  ! eUm  : Only for Si(110) and Ge(110) potent. [eV], Ge taken from A. Fomin, Si from initial implementation

  ! Si
  m = collmat_getCollMatID("Si")
  validMat(m) = .true.
  dlri(m)     =  0.0937_fPrec
  dlyi(m)     =  0.4652_fPrec
  ai(m)       =  0.96e-7_fPrec
  eUm(m)      = 21.34_fPrec
  collnt(m)   =  0.3016_fPrec

  ! W
  m = collmat_getCollMatID("W")
  validMat(m) = .true.
  dlri(m)     =  0.0035_fPrec
  dlyi(m)     =  0.096_fPrec
  ai(m)       =  0.56e-7_fPrec
  eUm(m)      = 21.0_fPrec
  collnt(m)   = zero

  ! C
  m = collmat_getCollMatID("C")
  validMat(m) = .true.
  dlri(m)     =  0.188_fPrec
  dlyi(m)     =  0.400_fPrec
  ai(m)       =  0.63e-7_fPrec
  eUm(m)      = 21.0_fPrec
  collnt(m)   = zero

  ! Ge
  m = collmat_getCollMatID("Ge")
  validMat(m) = .true.
  dlri(m)     =  0.02302_fPrec
  dlyi(m)     =  0.2686_fPrec
  ai(m)       =  1.0e-7_fPrec
  eUm(m)      = 40.0_fPrec
  collnt(m)   =  0.1632_fPrec

end subroutine cry_init

! ================================================================================================ !
!  Initialise a crystal collimator
! ================================================================================================ !
subroutine cry_startElement(icoll, ie, emitX, emitY, o_tilt, o_length)

  use crcoall
  use coll_db
  use coll_common
  use coll_materials
  use mathlib_bouncer
  use mod_common_track

  integer,          intent(in)    :: icoll, ie
  real(kind=fPrec), intent(in)    :: emitX, emitY
  real(kind=fPrec), intent(inout) :: o_tilt
  real(kind=fPrec), intent(inout) :: o_length

  integer mat
  real(kind=fPrec) bendAng,cry_tilt0

  mat = cdb_cMaterialID(icoll)
  if(validMat(mat) .eqv. .false.) then
    write(lerr,"(a)") "COLL> ERROR Crystal collimation not supported for material '"//trim(colmats(mat))//"'"
    call prror
  end if

  if(modulo(cdb_cRotation(icoll),pi) < c1m9) then
    cry_tilt0 = -(sqrt(emitX/tbetax(ie))*talphax(ie))*cdb_cNSig(icoll)
  elseif (modulo(cdb_cRotation(icoll)-pi2,pi) < c1m9) then
    cry_tilt0 = -(sqrt(emitY/tbetay(ie))*talphay(ie))*cdb_cNSig(icoll)
  else
    write(lerr,"(a)") "COLL> ERROR Crystal collimator has to be horizontal or vertical"
    call prror
  end if

  cry_tilt = cdb_cryTilt(icoll) + cry_tilt0
  bendAng  = cdb_cLength(icoll)/cdb_cryBend(icoll)
  if(cry_tilt >= (-one)*bendAng) then
    cry_length = cdb_cryBend(icoll)*(sin_mb(bendAng + cry_tilt) - sin_mb(cry_tilt))
  else
    cry_length = cdb_cryBend(icoll)*(sin_mb(bendAng - cry_tilt) + sin_mb(cry_tilt))
  end if

  c_rcurv  = cdb_cryBend(icoll)
  c_alayer = cdb_cryThick(icoll)
  c_xmax   = cdb_cryXDim(icoll)
  c_ymax   = cdb_cryYDim(icoll)
  c_orient = cdb_cryOrient(icoll)
  c_miscut = cdb_cryMiscut(icoll)
  cry_bend = cry_length/c_rcurv
  c_cBend  = cos_mb(cry_bend)
  c_sBend  = sin_mb(cry_bend)
  c_cpTilt = cos_mb(cry_tilt)
  c_spTilt = sin_mb(cry_tilt)
  c_cnTilt = cos_mb(-cry_tilt)
  c_snTilt = sin_mb(-cry_tilt)

  n_chan      = 0
  n_VR        = 0
  n_amorphous = 0

  cry_proc(:) = proc_out

  o_tilt   = cry_tilt
  o_length = cry_length

end subroutine cry_startElement

! ================================================================================================ !
!  Subroutine for checking for interactions with crystal
! ================================================================================================ !
subroutine cry_doCrystal(ie,iturn,j,mat,x,xp,z,zp,s,p,x0,xp0,zlm,s_imp,isImp,nhit,nabs, &
  lhit,lhit_turn,part_abs,part_abs_turn,impact,indiv,c_length)

  use parpro
  use coll_common, only : cry_proc, cry_proc_prev, cry_proc_tmp
  use mathlib_bouncer

  integer,          intent(in)    :: ie
  integer,          intent(in)    :: iturn
  integer,          intent(in)    :: j
  integer,          intent(in)    :: mat

  real(kind=fPrec), intent(inout) :: x,xp
  real(kind=fPrec), intent(inout) :: z,zp
  real(kind=fPrec), intent(inout) :: s,p
  real(kind=fPrec), intent(inout) :: x0,xp0
  real(kind=fPrec), intent(inout) :: zlm,s_imp
  integer,          intent(inout) :: nhit,nabs
  integer,          intent(inout) :: lhit(npart)
  integer,          intent(inout) :: lhit_turn(npart)
  integer,          intent(inout) :: part_abs(npart)
  integer,          intent(inout) :: part_abs_turn(npart)
  real(kind=fPrec), intent(inout) :: impact(npart)
  real(kind=fPrec), intent(inout) :: indiv(npart)
  real(kind=fPrec), intent(in)    :: c_length

  logical,          intent(inout) :: isImp

  real(kind=fPrec) s_temp,s_shift,s_rot,s_int
  real(kind=fPrec) x_temp,x_shift,x_rot,x_int
  real(kind=fPrec) xp_temp,xp_shift,xp_rot,xp_int,xp_tangent
  real(kind=fPrec) tilt_int,shift,delta,a_eq,b_eq,c_eq
  real(kind=fPrec) s_P, x_P, s_P_tmp, x_P_tmp

  s_temp     = zero
  s_shift    = zero
  s_rot      = zero
  s_int      = zero
  x_temp     = zero
  x_shift    = zero
  x_rot      = zero
  x_int      = zero
  xp_temp    = zero
  xp_shift   = zero
  xp_rot     = zero
  xp_int     = zero
  xp_tangent = zero
  tilt_int   = zero
  shift      = zero
  delta      = zero
  a_eq       = zero
  b_eq       = zero
  c_eq       = zero
  s_imp      = zero


  ! Determining if particle previously interacted with a crystal and storing the process ID
  if(cry_proc_tmp(j) /= proc_out) then
    cry_proc_prev(j) = cry_proc_tmp(j)
  end if

  iProc       = proc_out
  cry_proc(j) = proc_out

  ! Transform in the crystal reference system
  ! 1st transformation: shift of the center of the reference frame
  if(cry_tilt < zero) then
    s_shift = s
    shift   = c_rcurv*(one - c_cpTilt)
    if(cry_tilt < -cry_bend) then
      shift = c_rcurv*(c_cnTilt - cos_mb(cry_bend - cry_tilt))
    end if
    x_shift = x - shift
  else
    s_shift = s
    x_shift = x
  end if

  ! 2nd transformation: rotation
  s_rot  = x_shift*c_spTilt + s_shift*c_cpTilt
  x_rot  = x_shift*c_cpTilt - s_shift*c_spTilt
  xp_rot = xp - cry_tilt

  ! 3rd transformation: drift to the new coordinate s=0
  xp = xp_rot
  x  = x_rot - xp_rot*s_rot
  z  = z - zp*s_rot
  s  = zero

  ! Check that particle hit the crystal
  if(x >= zero .and. x < c_xmax) then

    ! MISCUT first step: P coordinates (center of curvature of crystalline planes)
    s_P = (c_rcurv-c_xmax)*sin_mb(-c_miscut)
    x_P = c_xmax + (c_rcurv-c_xmax)*cos_mb(-c_miscut)

    call cry_interact(mat,x,xp,z,zp,p,cry_length,s_P,x_P)
    s   = c_rcurv*c_sBend
    zlm = c_rcurv*c_sBend
    if(iProc /= proc_out) then
      isImp        = .true.
      nhit         = nhit + 1
      lhit(j)      = ie
      lhit_turn(j) = iturn
      impact(j)    = x0
      indiv(j)     = xp0
    end if

  else

    if(x < zero) then ! Crystal can be hit from below
      xp_tangent = sqrt((-(two*x)*c_rcurv + x**2)/(c_rcurv**2))
    else              ! Crystal can be hit from above
      xp_tangent = asin_mb((c_rcurv*(one - c_cBend) - x)/sqrt(((two*c_rcurv)*(c_rcurv - x))*(one - c_cBend) + x**2))
    end if

    ! If the hit is below, the angle must be greater or equal than the tangent,
    ! or if the hit is above, the angle must be smaller or equal than the tangent
    if((x < zero .and. xp >= xp_tangent) .or. (x >= zero .and. xp <= xp_tangent)) then

      ! If it hits the crystal, calculate in which point and apply the transformation and drift to that point
      a_eq  = one + xp**2
      b_eq  = (two*xp)*(x - c_rcurv)
      c_eq  = -(two*x)*c_rcurv + x**2
      delta = b_eq**2 - four*(a_eq*c_eq)
      s_int = (-b_eq - sqrt(delta))/(two*a_eq)
      s_imp = s_int

      ! MISCUT first step: P coordinates (center of curvature of crystalline planes)
      s_P_tmp = (c_rcurv-c_xmax)*sin_mb(-c_miscut)
      x_P_tmp = c_xmax + (c_rcurv-c_xmax)*cos_mb(-c_miscut)

      if(s_int < c_rcurv*c_sBend) then
        ! Transform to a new reference system: shift and rotate
        x_int  = xp*s_int + x
        xp_int = xp
        z      = z + zp*s_int
        x      = zero
        s      = zero

        tilt_int = s_int/c_rcurv
        xp       = xp-tilt_int

        ! MISCUT first step (bis): transform P in new reference system
        ! Translation
        s_P_tmp = s_P_tmp - s_int
        x_P_tmp = x_P_tmp - x_int
        ! Rotation
        s_P = s_P_tmp*cos_mb(tilt_int) + x_P_tmp*sin_mb(tilt_int)
        x_P = -s_P_tmp*sin_mb(tilt_int) + x_P_tmp*cos_mb(tilt_int)

        call cry_interact(mat,x,xp,z,zp,p,cry_length-(tilt_int*c_rcurv),s_P,x_P)
        s   = c_rcurv*sin_mb(cry_bend - tilt_int)
        zlm = c_rcurv*sin_mb(cry_bend - tilt_int)
        if(iProc /= proc_out) then
          x_rot    = x_int
          s_rot    = s_int
          xp_rot   = xp_int
          s_shift  =  s_rot*c_cnTilt + x_rot*c_snTilt
          x_shift  = -s_rot*c_snTilt + x_rot*c_cnTilt
          xp_shift = xp_rot + cry_tilt

          if(cry_tilt < zero) then
            x0  = x_shift + shift
            xp0 = xp_shift
          else
            x0  = x_shift
            xp0 = xp_shift
          end if

          isImp        = .true.
          nhit         = nhit + 1
          lhit(j)      = ie
          lhit_turn(j) = iturn
          impact(j)    = x0
          indiv(j)     = xp0
        end if

        ! un-rotate
        x_temp  = x
        s_temp  = s
        xp_temp = xp
        s       =  s_temp*cos_mb(-tilt_int) + x_temp*sin_mb(-tilt_int)
        x       = -s_temp*sin_mb(-tilt_int) + x_temp*cos_mb(-tilt_int)
        xp      = xp_temp + tilt_int

        ! 2nd: shift back the 2 axis
        x = x + x_int
        s = s + s_int

      else

        s = c_rcurv*sin_mb(cry_length/c_rcurv)
        x = x + s*xp
        z = z + s*zp

      end if

    else

      s = c_rcurv*sin_mb(cry_length/c_rcurv)
      x = x + s*xp
      z = z + s*zp

    end if

  end if

  ! transform back from the crystal to the collimator reference system
  ! 1st: un-rotate the coordinates
  x_rot  = x
  s_rot  = s
  xp_rot = xp

  s_shift  =  s_rot*c_cnTilt + x_rot*c_snTilt
  x_shift  = -s_rot*c_snTilt + x_rot*c_cnTilt
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
  x = xp*(c_length - s) + x
  z = zp*(c_length - s) + z
  s = c_length

  nabs = 0
  cry_proc(j) = iProc
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

  ! Storing the process ID for the next interaction
  cry_proc_tmp(j) = cry_proc(j)

end subroutine cry_doCrystal

! ================================================================================================ !
!  Subroutine for the movements of the particles in the crystal
!  Simple tranport protons in crystal 2
! ================================================================================================ !
subroutine cry_interact(is,x,xp,y,yp,pc,length,s_P,x_P)

  use mod_ranlux
  use mod_funlux
  use mod_common_main
  use floatPrecision
  use coll_materials, only : zatom, exenergy, rho, anuc
  use mathlib_bouncer
  use physical_constants

  integer,          intent(in)    :: is  ! Material number
  real(kind=fPrec), intent(inout) :: x
  real(kind=fPrec), intent(inout) :: xp
  real(kind=fPrec), intent(inout) :: y
  real(kind=fPrec), intent(inout) :: yp
  real(kind=fPrec), intent(inout) :: pc
  real(kind=fPrec), intent(in)    :: length
  real(kind=fPrec), intent(in)    :: s_P
  real(kind=fPrec), intent(in)    :: x_P

  integer nam,zn                        ! Switch on/off the nuclear interaction (NAM) and the MCS (ZN)
  real(kind=fPrec) ymax,ymin            ! Crystal geometrical parameters
  real(kind=fPrec) s_length             ! Element length along s
  real(kind=fPrec) DESt                 ! Changed energy loss by ionization now calculated and not tabulated
  real(kind=fPrec) x0,y0                ! Coordinates of the particle [m]
  real(kind=fPrec) s                    ! Long coordinates of the particle [m]
  real(kind=fPrec) a_eq,b_eq,c_eq,delta ! Second order equation param.
  real(kind=fPrec) Ang_rms, Ang_avr     ! Volume reflection mean angle [rad]
  real(kind=fPrec) Dechan               ! Probability for dechanneling
  real(kind=fPrec) Lrefl,Srefl          ! Distance of the reflection point [m]
  real(kind=fPrec) Vcapt                ! Volume capture probability
  real(kind=fPrec) Chann                ! Channeling probability
  real(kind=fPrec) N_atom               ! Probability for entering channeling near atomic planes
  real(kind=fPrec) Dxp                  ! Variation in angle
  real(kind=fPrec) xpcrit               ! Critical angle for curved crystal[rad]
  real(kind=fPrec) xpcrit0              ! Critical angle for str. crystal [rad]
  real(kind=fPrec) Rcrit                ! Critical curvature radius [m]
  real(kind=fPrec) ratio                ! X=c_rcurv/Rcrit
  real(kind=fPrec) TLdech1,TLdech2      ! Typical dechanneling length(2) [m]
  real(kind=fPrec) tdech,Ldech,Sdech    ! Angle, lenght, and S coordinate of dechanneling point
  real(kind=fPrec) Rlength,Red_S        ! Reduced length/s coordinate (in case of dechanneling)
  real(kind=fPrec) am_len               ! Amorphous length
  real(kind=fPrec) len_xs,len_ys        ! Amorphous length
  real(kind=fPrec) xp_rel               ! Xp-c_miscut angle in mrad
  real(kind=fPrec) alpha                ! Par for new chann prob
  real(kind=fPrec) Pvr                  ! Prob for VR->AM transition

  ! Quantities for length and deflection calculation
  real(kind=fPrec) const_dech,xpin,ypin,tchan,tdefl,tP,L_chan,mep
  real(kind=fPrec) s_K,x_K,s_M,x_M,s_F,x_F,r,a
  real(kind=fPrec) A_F,B_F,C_F,alpha_F,beta_F

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
    (one/(log_mb(((two*pmae)*gammar)/(exenergy(is)*c1e3)) - one))) * ((aTF*dP)/(crade*pmae)) ! [m/MeV]
  const_dech = const_dech*c1e3 ! [m/GeV]

  s        = zero
  s_length = c_rcurv*sin_mb(length/c_rcurv)
  L_chan   = length

  ! MISCUT second step: fundamental coordinates (crystal edges and plane curvature radius)
  s_K = c_rcurv*sin_mb(length/c_rcurv)
  x_K = c_rcurv*(1-cos_mb(length/c_rcurv))
  s_M = (c_rcurv-c_xmax)*sin_mb(length/c_rcurv)
  x_M = c_xmax + (c_rcurv-c_xmax)*(1-cos_mb(length/c_rcurv))
  r   = sqrt(s_P**2 + (x-x_P)**2)

  ! MISCUT third step: F coordinates (exit point) on crystal exit face
  A_F = (tan_mb(length/c_rcurv))**2 + one
  B_F = ((-two)*(tan_mb(length/c_rcurv))**2)*c_rcurv + (two*tan_mb(length/c_rcurv))*s_P - two*x_P
  C_F = ((tan_mb(length/c_rcurv))**2)*(c_rcurv**2) - ((two*tan_mb(length/c_rcurv))*s_P)*c_rcurv + s_P**2 + x_P**2 - r**2

  x_F = (-B_F-sqrt(B_F**2-four*(A_F*C_F)))/(two*A_F)
  s_F = (-tan_mb(length/c_rcurv))*(x_F-c_rcurv)

  if(x_F >= x_K .and. x_F <= x_M .and. s_F >= s_M .and. s_F <= s_K) then
    ! No additional steps required for miscut
  else if (c_miscut == 0 .and. abs(x_F-x_K) <= c1m13 .and. abs(s_F-s_K) <= c1m13) then
    ! no miscut, entrance from below: exit point is K (lower edge)
    x_F = x_K
    s_F = s_K
  else if (c_miscut == 0 .and. abs(x_F-x_M) <= c1m13 .and. abs(s_F-s_M) <= c1m13) then
    ! no miscut, entrance from above: exit point is M (upper edge)
    x_F = x_M
    s_F = s_M
  else

    ! MISCUT Third step (bis): F coordinates (exit point)  on bent side
    if(c_miscut < 0) then
      ! Intersect with bottom side
      alpha_F = (c_rcurv-x_P)/x_P
      beta_F = -(r**2-s_P**2-x_P**2)/(two*s_P)
      A_F = alpha_F**2 + one
      B_F = two*(alpha_F*beta_F) - two*c_rcurv
      C_F = beta_F**2
    else
      ! Intersect with top side
      alpha_F = (c_rcurv-x_P)/s_P
      beta_F = -(r**2+c_xmax*(c_xmax-(two*c_rcurv))-s_P**2-x_P**2)/(two*s_P)
      A_F = alpha_F**2 + one
      B_F = two*(alpha_F*beta_F) - two*c_rcurv
      C_F = beta_F**2 - c_xmax*(c_xmax-two*c_rcurv)
    endif
    
    x_F = (-B_F-sqrt(B_F**2-four*(A_F*C_F)))/(two*A_F)
    s_F = alpha_F*x_F + beta_F
    
  endif

  ! MISCUT fourth step: deflection and length calculation
  a = sqrt(s_F**2+(x-x_F)**2)
  tP = acos_mb((2*(r**2)-a**2)/(2*(r**2)))
  tdefl = asin_mb((s_f-s_P)/r)
  L_chan = r*tP

  xp_rel = xp - c_miscut

  ymin = -c_ymax/two
  ymax =  c_ymax/two

  ! FIRST CASE: p don't interact with crystal
  if(y < ymin .or. y > ymax .or. x > c_xmax) then
    x     = x + xp*s_length
    y     = y + yp*s_length
    iProc = proc_out
    return

  ! SECOND CASE: p hits the amorphous layer
  else if(x < c_alayer .or. y-ymin < c_alayer .or. ymax-y < c_alayer) then
    x0    = x
    y0    = y
    a_eq  = one + xp**2
    b_eq  = (two*x)*xp - (two*xp)*c_rcurv
    c_eq  = x**2 - (two*x)*c_rcurv
    delta = b_eq**2 - (four*a_eq)*c_eq
    s     = (-b_eq+sqrt(delta))/(two*a_eq)
    if(s >= s_length) then
      s = s_length
    end if
    x         =  xp*s + x0
    len_xs = sqrt((x-x0)**2 + s**2)
    if(yp >= zero .and. y + yp*s <= ymax) then
      len_ys = yp*len_xs
    else if(yp < zero .and. y + yp*s >= ymin) then
      len_ys = yp*len_xs
    else
      s      = (ymax-y)/yp
      len_ys = sqrt((ymax-y)**2 + s**2)
      x      = x0 + xp*s
      len_xs = sqrt((x-x0)**2 + s**2)
    end if
    am_len = sqrt(len_xs**2 + len_ys**2)
    s     = s/two
    x     = x0 + xp*s
    y     = y0 + yp*s
    iProc = proc_AM
    call cry_calcIonLoss(is,pc,am_len,dest)
    call cry_moveAM(is,nam,am_len,dest,dlyi(is),dlri(is),xp,yp,pc)
    x = x + xp*(s_length-s)
    y = y + yp*(s_length-s)
    return

  else if(x > c_xmax-c_alayer .and. x < c_xmax) then
    iProc = proc_AM
    call cry_calcIonLoss(is,pc,s_length,dest)
    call cry_moveAM(is,nam,s_length,dest,dlyi(is),dlri(is),xp,yp,pc)
    return

  end if

  ! THIRD CASE: the p interacts with the crystal.
  ! Define typical angles/probabilities for orientation 110
  xpcrit0 = sqrt((c2m9*eUm(is))/pc)    ! Critical angle (rad) for straight crystals
  Rcrit   = (pc/(c2m6*eUm(is)))*ai(is) ! Critical curvature radius [m]

  ! If R>Rcritical=>no channeling is possible (ratio<1)
  ratio  = c_rcurv/Rcrit
  xpcrit = (xpcrit0*(c_rcurv-Rcrit))/c_rcurv ! Critical angle for curved crystal

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

  if(c_orient == 2) then
    Ang_avr = Ang_avr*0.93_fPrec
    Ang_rms = Ang_rms*1.05_fPrec
    xpcrit  = xpcrit*0.98_fPrec
  end if

  if(abs(xp_rel) < xpcrit) then
    alpha  = xp_rel/xpcrit
    Chann  = sqrt(0.9_fPrec*(one - alpha**2))*sqrt(one-(one/ratio)) ! Saturation at 95%
    N_atom = c1m1

    ! if they can channel: 2 options
    if(coll_rand() <= chann) then ! option 1:channeling

      TLdech1 = (const_dech*pc)*(one-one/ratio)**2 ! Updated calculate typical dech. length(m)
      if(coll_rand() <= n_atom) then
        TLdech1 = ((const_dech/c2e2)*pc)*(one-one/ratio)**2  ! Updated dechanneling length (m)
      end if

      Dechan = -log_mb(coll_rand()) ! Probability of dechanneling
      Ldech  = TLdech1*Dechan   ! Actual dechan. length

      ! careful: the dechanneling lentgh is along the trajectory
      ! of the particle -not along the longitudinal coordinate...
      if(ldech < l_chan) then
        iProc = proc_DC
        Dxp   = Ldech/r ! Change angle from channeling [mrad]
        Sdech = Ldech*cos_mb(c_miscut + half*Dxp)
        x     = x  + Ldech*(sin_mb(half*Dxp+c_miscut)) ! Trajectory at channeling exit
        xp    = xp + Dxp + (two*(coll_rand()-half))*xpcrit
        y     = y  + yp * Sdech

        call cry_calcIonLoss(is,pc,ldech,dest)
        pc = pc - half*dest*Ldech ! Energy loss to ionization while in CH [GeV]
        x  = x  + (half*(s_length-Sdech))*xp
        y  = y  + (half*(s_length-Sdech))*yp

        call cry_calcIonLoss(is,pc,s_length-sdech,dest)
        call cry_moveAM(is,nam,s_length-sdech,dest,dlyi(is),dlri(is),xp,yp,pc)
        x = x + (half*(s_length-Sdech))*xp
        y = y + (half*(s_length-Sdech))*yp
      else
        iProc = proc_CH
        xpin  = XP
        ypin  = YP

        call cry_moveCH(is,nam,l_chan,x,xp,yp,pc,c_rcurv,rcrit) ! check if a nuclear interaction happen while in CH
        if(iProc /= proc_CH) then
          ! if an nuclear interaction happened, move until the middle with initial xp,yp then
          ! propagate until the "crystal exit" with the new xp,yp accordingly with the rest
          ! of the code in "thin lens approx"
          x = x + (half*L_chan)*xpin
          y = y + (half*L_chan)*ypin
          x = x + (half*L_chan)*XP
          y = y + (half*L_chan)*YP

          call cry_calcIonLoss(is,pc,length,dest)
          pc = pc - dest*length ! energy loss to ionization [GeV]
        else
          Dxp = tdefl + (half*ran_gauss(zero))*xpcrit ! Change angle[rad]
        
          xp  = Dxp
          x   = x + L_chan*(sin_mb(half*Dxp)) ! Trajectory at channeling exit
          y   = y + s_length * yp

          call cry_calcIonLoss(is,pc,length,dest)
          pc = pc - (half*dest)*length ! energy loss to ionization [GeV]
        end if
      end if

    else ! Option 2: VR

      ! good for channeling but don't channel (1-2)
      iProc = proc_VR

      xp = xp + (0.45_fPrec*(xp_rel/xpcrit + one))*Ang_avr
      x  = x  + (half*s_length)*xp
      y  = y  + (half*s_length)*yp

      call cry_calcIonLoss(is,pc,s_length,dest)
      call cry_moveAM(is,nam,s_length,dest,dlyi(is),dlri(is),xp,yp,pc)

      x = x + (half*s_length)*xp
      y = y + (half*s_length)*yp

    end if

  else ! case 3-2: no good for channeling. check if the  can VR

    Lrefl = xp_rel*r ! Distance of refl. point [m]
    Srefl = sin_mb(xp_rel/two + c_miscut)*Lrefl

    if(Lrefl > zero .and. Lrefl < L_chan) then ! VR point inside

      ! 2 options: volume capture and volume reflection

      if(coll_rand() > Vcapt .or. ZN == zero) then ! Option 1: VR

        iProc = proc_VR
        x     = x + xp*Srefl
        y     = y + yp*Srefl
        Dxp   = Ang_avr
        xp    = xp + Dxp + Ang_rms*ran_gauss(zero)
        x     = x  + (half*xp)*(s_length - Srefl)
        y     = y  + (half*yp)*(s_length - Srefl)

        call cry_calcIonLoss(is,pc,s_length-srefl,dest)
        call cry_moveAM(is,nam,s_length-srefl,dest,dlyi(is),dlri(is),xp,yp,pc)
        x = x + (half*xp)*(s_length - Srefl)
        y = y + (half*yp)*(s_length - Srefl)

      else ! Option 2: VC

        x = x + xp*Srefl
        y = y + yp*Srefl

        TLdech2 = (const_dech/c1e1)*pc*(one-one/ratio)**2          ! Updated typical dechanneling length(m)
        Ldech   = TLdech2*(sqrt(c1m2 - log_mb(coll_rand())) - c1m1)**2 ! Updated DC length
        tdech   = Ldech/c_rcurv
        Sdech   = Ldech*cos_mb(xp + half*tdech)

        if(Ldech < Length-Lrefl) then

          iProc = proc_DC
          Dxp   = Ldech/c_rcurv + (half*ran_gauss(zero))*xpcrit
          x     = x + Ldech*(sin_mb(half*Dxp+xp)) ! Trajectory at channeling exit
          y     = y + Sdech*yp
          xp    =  Dxp
          Red_S = (s_length - Srefl) - Sdech
          x     = x + (half*xp)*Red_S
          y     = y + (half*yp)*Red_S

          call cry_calcIonLoss(is,pc,srefl,dest)
          pc = pc - dest*Srefl ! "added" energy loss before capture

          call cry_calcIonLoss(is,pc,sdech,dest)
          pc = pc - (half*dest)*Sdech ! "added" energy loss while captured

          call cry_calcIonLoss(is,pc,red_s,dest)
          call cry_moveAM(is,nam,red_s,dest,dlyi(is),dlri(is),xp,yp,pc)
          x = x + (half*xp)*Red_S
          y = y + (half*yp)*Red_S

        else

          iProc   = proc_VC
          Rlength = Length - Lrefl
          tchan   = Rlength/c_rcurv
          Red_S   = Rlength*cos_mb(xp + half*tchan)

          call cry_calcIonLoss(is,pc,lrefl,dest)
          pc   = pc - dest*Lrefl ! "added" energy loss before capture
          xpin = xp
          ypin = yp

          call cry_moveCH(is,nam,rlength,x,xp,yp,pc,c_rcurv,rcrit) ! Check if a nuclear interaction happen while in ch
          if(iProc /= proc_VC) then
            ! if an nuclear interaction happened, move until the middle with initial xp,yp then propagate until
            ! the "crystal exit" with the new xp,yp accordingly with the rest of the code in "thin lens approx"
            x = x + (half*Rlength)*xpin
            y = y + (half*Rlength)*ypin
            x = x + (half*Rlength)*XP
            y = y + (half*Rlength)*YP

            call cry_calcIonLoss(is,pc,rlength,dest)
            pc = pc - dest*Rlength
          else
            Dxp = (Length-Lrefl)/c_rcurv
            x   = x + sin_mb(half*Dxp+xp)*Rlength ! Trajectory at channeling exit
            y   = y + red_S*yp
            xp  = tdefl + (half*ran_gauss(zero))*xpcrit ! [mrad]

            call cry_calcIonLoss(is,pc,rlength,dest)
            pc = pc - (half*dest)*Rlength  ! "added" energy loss once captured
          end if
        end if
      end if

    else

      ! Case 3-3: move in amorphous substance (big input angles)
      ! Modified for transition vram daniele
      if(xp_rel > tdefl-c_miscut + two*xpcrit .or. xp_rel < -xpcrit) then
        iProc = proc_AM
        x     = x + (half*s_length)*xp
        y     = y + (half*s_length)*yp
        if(zn > zero) then
          call cry_calcIonLoss(is,pc,s_length,dest)
          call cry_moveAM(is,nam,s_length,dest,dlyi(is),dlri(is),xp,yp,pc)
        end if
        x = x + (half*s_length)*xp
        y = y + (half*s_length)*yp
      else
        Pvr = (xp_rel-(tdefl-c_miscut))/(two*xpcrit)
        if(coll_rand() > Pvr) then
          iProc = proc_TRVR
          x     = x + xp*Srefl
          y     = y + yp*Srefl

          Dxp = (((-three*Ang_rms)*xp_rel)/(two*xpcrit) + Ang_avr) + ((three*Ang_rms)*(tdefl-c_miscut))/(two*xpcrit)
          xp  = xp + Dxp
          x   = x + (half*xp)*(s_length-Srefl)
          y   = y + (half*yp)*(s_length-Srefl)

          call cry_calcIonLoss(is,pc,s_length-srefl,dest)
          call cry_moveAM(is,nam,s_length-srefl,dest,dlyi(is),dlri(is),xp,yp,pc)
          x = x + (half*xp)*(s_length - Srefl)
          y = y + (half*yp)*(s_length - Srefl)
        else
          iProc = proc_TRAM
          x = x + xp*Srefl
          y = y + yp*Srefl
          Dxp = ((((-one*(13.6_fPrec/pc))*sqrt(s_length/dlri(is)))*c1m3)*xp_rel)/(two*xpcrit) + &
            (((13.6_fPrec/pc)*sqrt(s_length/DLRi(is)))*c1m3)*(one+(tdefl-c_miscut)/(two*xpcrit))
          xp = xp+Dxp
          x  = x + (half*xp)*(s_length-Srefl)
          y  = y + (half*yp)*(s_length-Srefl)

          call cry_calcIonLoss(is,pc,s_length-srefl,dest)
          call cry_moveAM(is,nam,s_length-srefl,dest,dlyi(is),dlri(is),xp,yp,pc)
          x = x + (half*xp)*(s_length - Srefl)
          y = y + (half*yp)*(s_length - Srefl)
        end if
      end if
    end if
  end if

end subroutine cry_interact

! ================================================================================================ !
!  Subroutine for the calculazion of the energy loss by ionisation
! ================================================================================================ !
subroutine cry_calcIonLoss(is,pc,dz,EnLo)

  use mod_ranlux
  use mod_funlux
  use floatPrecision
  use coll_materials, only : zatom, exenergy, rho, anuc
  use mathlib_bouncer
  use physical_constants

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

  if(coll_rand() < prob_tail) then
    EnLo = ((k*zatom(is))/(anuc(is)*betar**2)) * ( &
      half*log_mb((two*pmae*bgr*bgr*Tmax)/(c1e6*exenergy(is)**2)) -      &
      betar**2 - log_mb(plen/(exenergy(is)*c1e3)) - log_mb(bgr) + half + &
      TMax**2/(eight*(gammar**2)*(pmap**2)) &
    )
    EnLo = (EnLo*rho(is))*c1m1 ! [GeV/m]
  else
    EnLo = EnLo/dz ! [GeV/m]
  end if

end subroutine cry_calcIonLoss

! ================================================================================================ !
!  Subroutine for the movement in the amorphous
! ================================================================================================ !
subroutine cry_moveAM(is,nam,dz,dei,dly,dlr,xp,yp,pc)

  use mod_ranlux
  use mod_funlux
  use floatPrecision
  use coll_materials, only : anuc, hcut, bnref, csref
  use mathlib_bouncer
  use physical_constants

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
  call funlxp(cry_ruth,cgen_cry(1,is),tlow,thigh)

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
  kxmcs = dya*ran_gauss(zero)
  kymcs = dya*ran_gauss(zero)

  xp = xp+kxmcs
  yp = yp+kymcs

  if(nam == 0) return ! Turn on/off nuclear interactions

  ! Can nuclear interaction happen?
  zlm = -collnt(is)*log_mb(coll_rand())

  if(zlm < dz) then
    ! Choose nuclear interaction
    aran = coll_rand()
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
      t     = -log_mb(coll_rand())/bn

    case(3) ! p-p elastic
      iProc = proc_ppe
      t     = -log_mb(coll_rand())/bpp

    case(4) ! Single diffractive
      iProc = proc_diff
      xm2   = exp_mb(coll_rand()*xln15s)
      pc    = pc*(one - xm2/ecmsq)
      if(xm2 < two) then
        bsd = two*bpp
      else if(xm2 >= two .and. xm2 <= five) then
        bsd = ((106.0_fPrec - 17.0_fPrec*xm2)*bpp)/36.0_fPrec
      else if(xm2 > five) then
        bsd = 7.0_fPrec*bpp/12.0_fPrec
      end if
      t = -log_mb(coll_rand())/bsd

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

    tx = (teta*ran_gauss(zero))*c1e3
    tz = (teta*ran_gauss(zero))*c1e3

    ! Change p angle
    xp = xp + tx
    yp = yp + tz
  end if

  xp = xp/c1e3
  yp = yp/c1e3

end subroutine cry_moveAM

! ================================================================================================ !
!  Subroutine for check if a nuclear interaction happen while in channeling
! ================================================================================================ !
subroutine cry_moveCH(is,nam,dz,x,xp,yp,pc,r,rc)

  use crcoall
  use mod_ranlux
  use mod_funlux
  use floatPrecision
  use coll_common, only : coll_debug
  use coll_materials, only : nmat, rho, anuc, hcut, bnref, csref, csect
  use mathlib_bouncer
  use physical_constants

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
  call funlxp(cry_ruth,cgen_cry(1,is),tlow,thigh)

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
  zlm = -nuc_cl_l*log_mb(coll_rand())

  ! write(889,*) x_i,pv,Ueff,Et,Ec,N_am,avrrho,csref_tot_rsc,csref_inel_rsc,nuc_cl_l

  if(zlm < dz) then
    ! Choose nuclear interaction
    aran = coll_rand()
    i=1
10  if(aran > cprob(i)) then
      i=i+1
      goto 10
    end if
    ichoix = i

    ! Do the interaction
    t = 0 ! default value to cover ichoix=1
    select case(ichoix)
    case(1) ! deep inelastic, impinging p disappeared
      iProc = proc_ch_absorbed

    case(2) ! p-n elastic
      iProc = proc_ch_pne
      bn    = (bnref(is)*cs(0))/csref_tot_rsc
      t     = -log_mb(coll_rand())/bn

    case(3) ! p-p elastic
      iProc = proc_ch_ppe
      t     = -log_mb(coll_rand())/bpp

    case(4) ! Single diffractive
      iProc = proc_ch_diff
      xm2   = exp_mb(coll_rand()*xln15s)
      pc    = pc*(one - xm2/ecmsq)
      if(xm2 < two) then
        bsd = two*bpp
      else if(xm2 >= two .and. xm2 <= five) then
        bsd = ((106.0_fPrec - 17.0_fPrec*xm2)*bpp)/36.0_fPrec
      else if(xm2 > five) then
        bsd = (seven*bpp)/12.0_fPrec
      end if
      t = -log_mb(coll_rand())/bsd

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

    tx = (teta*ran_gauss(zero))*c1e3
    tz = (teta*ran_gauss(zero))*c1e3

    ! Change p angle
    xp = xp + tx
    yp = yp + tz

  end if

  xp = xp/c1e3
  yp = yp/c1e3

end subroutine cry_moveCH

! ================================================================================================ !
! Definition of rutherford scattering formula
! ================================================================================================ !
function cry_ruth(t_cry)

  use floatPrecision
  use coll_materials
  use mathlib_bouncer

  real(kind=fPrec) cry_ruth,t_cry
  real(kind=fPrec), parameter :: cnorm  = 2.607e-4_fPrec
  real(kind=fPrec), parameter :: cnform = 0.8561e3_fPrec

  cry_ruth = (cnorm*exp_mb((-t_cry*cnform)*emr(mcurr_cry)**2))*(zatom(mcurr_cry)/t_cry)**2

end function cry_ruth

end module coll_crystal
