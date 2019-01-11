module aperture
  ! Aperture check module
  ! A.Mereghetti, P.Garcia Ortega and D.Sinuela Pastor, for the FLUKA Team
  ! J.Molson, BE/ABP-HSS
  ! K.Sjobak, BE/ABP-LAT

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

  use parpro !For nele, npart

  !contains pstop(npart) etc
  use mod_commonmn
  use crcoall
  use mod_common
  use mod_commons
  use mod_commont
  use mod_commond

  use mod_hions
  use mod_alloc
#ifdef HDF5
  use hdf5_output
#endif

  implicit none

  ! A.Mereghetti, P.Garcia Ortega and D.Sinuela Pastor, for the FLUKA Team
  ! last modified: 02-03-2018
  ! always in main code

  logical, save :: limifound                       ! limi block in fort.3

  logical, save :: aperture_debug                  ! Enable/disable debugging output for the aperture code

  integer, allocatable, save :: kape(:)            ! type of aperture (nele)
  ! aperture parameteres ape(9,nele)
  ! ape(1,:): hor dimension (RECT/RECTELLIPSE/OCT) [mm]
  ! ape(2,:): ver dimension (RECT/RECTELLIPSE/OCT) [mm]
  ! ape(3,:): hor dimension (CIRC/ELLI/RECTELLIPSE/RACETR) [mm]
  ! ape(4,:): ver dimension (CIRC/ELLI/RECTELLIPSE/RACETR) [mm]
  ! ape(5,:): m of sloped side (OCT) []
  ! ape(6,:): q of sloped side (OCT) [mm]
  ! ape(7,:): tilt angle of marker (all) [rad]
  ! ape(8,:): hor offset of marker (all) [mm]
  ! ape(9,:): ver offset of marker (all) [mm]
  real(kind=fPrec), allocatable, save ::  ape(:,:) !(9,nele)
  logical, allocatable, save :: lapeofftlt(:)      ! aperture is tilted/offcentred (nele)

  ! save (i.e. do not kill) lost particles
  logical, save :: apflag                          ! save or not
  integer, allocatable, save :: plost(:)           ! particle ID (npart)

  integer, save :: aperture_napxStart              ! initial napx

  ! dump aperture profile:
  logical, save :: ldmpaper                        ! dump or not
  integer, save :: aperunit                        ! fortran unit
  character(len=16), save :: aper_filename         ! file name
  logical, save :: ldmpaperMem                     ! dump aperture marker parameters as in memory
  ! File for aperture losses
  integer, save :: losses_unit                                             ! unit
  character(len=19), parameter :: losses_filename="aperture_losses.dat"    ! name

  ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
  ! last modified: 02-03-2018
  ! variables for back-tracking
  logical, save :: lbacktracking                      ! activate back-tracking
  real(kind=fPrec), allocatable, save :: xLast(:,:)   ! position after last thick element [mm] (2,npart)
  real(kind=fPrec), allocatable, save :: yLast(:,:)   ! angles after last thick element [mrad] (2,npart)
  real(kind=fPrec), allocatable, save :: ejfvLast(:)  ! linear momentum [MeV/c] (npart)
  real(kind=fPrec), allocatable, save :: ejvLast(:)   ! total energy [MeV] (npart)
  real(kind=fPrec), allocatable, save :: nucmLast(:)  ! nuclear mass [GeV/c2] (npart)
  real(kind=fPrec), allocatable, save :: sigmvLast(:) ! lag [mm] (npart)
  real(kind=fPrec), allocatable, save :: dpsvLast(:)  ! (npart)
  integer, allocatable, save :: naaLast(:)            ! nuclear mass [] (npart)
  integer, allocatable, save :: nzzLast(:)            ! atomic number [] (npart)
  real(kind=fPrec), save :: bktpre                 ! precision of back-tracking [m]
  integer, save :: iLast, ixLast                   ! indeces of last aperture marker
  integer, save :: iLastThick, ixLastThick         ! indeces of last thick element
  integer, save :: iBckTypeLast                    ! map of back-tracking - it follows kz values, eg:
                                                   ! -1: generic (eg after aperture check)
                                                   !  0: drift (the only one available)

  ! A.Mereghetti (CERN, BE/ABP-HSS), 2018-03-22
  ! x-sec at specific locations
  integer, save :: mxsec                           ! current number of requested x-secs
  integer, parameter :: nxsec=10                   ! max number of requested x-secs
  integer, save :: xsecunit(nxsec)                 ! fortran units
  character(len=16), save :: xsec_filename(nxsec)  ! file names
  real(kind=fPrec), save :: sLocMin(nxsec), sLocMax(nxsec), sLocDel(nxsec) ! locations
  integer, save :: nAzimuts(nxsec)                 ! number of points (azimuth angles)
  integer, parameter :: nAzimutDef=72              ! default number of points


  ! aperture types  -- kape
  ! no aperture     -- 0
  ! circle          -- 1
  ! rectangle       -- 2
  ! ellipse         -- 3
  ! rectellipse     -- 4
  ! octagon         -- 5
  ! racetrack       -- 6
  ! transition      -- 7
  character(len=2), parameter, dimension(-1:6) :: apeName=(/'TR','NA','CR','RE','EL','RL','OC','RT'/)

  ! precision parameters:
  real(kind=fPrec), parameter :: aPrec=c1m6 ! identify two ap. markers as identical [mm]
  real(kind=fPrec), parameter :: sPrec=c1m7 ! identify two ap. markers as at the same s-pos [m]

#ifdef HDF5
  integer, private, save :: aper_fmtLostPart
  integer, private, save :: aper_setLostPart
#endif

#ifdef CR
  ! For resetting file positions
  integer, private, save :: apefilepos=-1, apefilepos_cr
#endif

contains

subroutine aperture_expand_arrays(nele_new, npart_new)

  implicit none
  integer, intent(in) :: nele_new
  integer, intent(in) :: npart_new

  call alloc(kape,       nele_new, 0, 'kape')
  call alloc(lapeofftlt, nele_new, .false., 'lapeofftlt')
  call alloc(ape, 9,     nele_new, zero, 'ape')

  call alloc(plost,     npart_new, 0, "plost")        ! particle ID (npart)
  call alloc(xLast, 2,  npart_new, zero, "xLast")     ! position after last thick element [mm] (2,npart)
  call alloc(yLast, 2,  npart_new, zero, "yLast")     ! angles after last thick element [mrad] (2,npart)
  call alloc(ejfvLast,  npart_new, zero, "ejfvLast")  ! linear momentum [MeV/c] (npart)
  call alloc(ejvLast,   npart_new, zero, "ejvLast")   ! total energy [MeV] (npart)
  call alloc(nucmLast,  npart_new, zero, "nucmLast")  ! nuclear mass [GeV/c2] (npart)
  call alloc(sigmvLast, npart_new, zero, "sigmvLast") ! lag [mm] (npart)
  call alloc(dpsvLast,  npart_new, zero, "dpsvLast")  ! (npart)
  call alloc(naaLast,   npart_new, 0, "naaLast")      ! nuclear mass [] (npart)
  call alloc(nzzLast,   npart_new, 0, "nzzLast")      ! atomic number [] (npart)

end subroutine aperture_expand_arrays


subroutine aperture_comnul

  use numerical_constants
  implicit none
  integer ii, jj

  limifound=.false.

  aperture_debug=.false.

  apflag=.false.

  do ii=1,npart
    plost(ii)=0
  end do

  ldmpaper            = .false.
  aperunit            = -1
  aper_filename(1:16) = ' '
  ldmpaperMem         = .false.

  lbacktracking = .false. ! backtracking off by default
  ! do ii=1,npart
  !   do jj=1,2
  !     xLast(jj,ii) = zero
  !     yLast(jj,ii) = zero
  !   end do
  !   ejfvLast(jj)  = zero
  !   ejvLast(jj)   = zero
  !   nucmLast(jj)  = zero
  !   sigmvLast(jj) = zero
  !   dpsvLast(jj)  = zero
  !   naaLast(jj)   = 0
  !   nzzLast(jj)   = 0
  ! end do

  bktpre=c1m1 ! default precision: 0.1m
  iLast = 0
  ixLast = 0
  iLastThick = 0
  ixLastThick = 0
  iBckTypeLast = -1

  mxsec = 0
  do ii=1,nxsec
    xsecunit(ii) = 0
    xsec_filename(ii)(1:16) = ' '
    sLocMin(ii) = zero
    sLocMax(ii) = zero
    sLocDel(ii) = zero
    nAzimuts(ii) = nAzimutDef
  end do

  aperture_napxStart = 0

  return

end subroutine aperture_comnul

! ================================================================================================ !
!  Aperture module initialisation
!  V.K. Berglyd Olsen, BR-ABP-HSS
!  Last modified: 2018-05-15
! ================================================================================================ !
subroutine aperture_init

  use mod_units, only: f_open, f_requestUnit

  implicit none

#ifdef HDF5
  type(h5_dataField), allocatable :: setFields(:)
#endif
  logical isOpen,err

#ifdef HDF5
  if(h5_useForAPER) then
    call h5_initForAperture
#ifdef FLUKA
    allocate(setFields(17))
#else
    allocate(setFields(15))
#endif
    setFields(1)  = h5_dataField(name="TURN",         type=h5_typeInt)
    setFields(2)  = h5_dataField(name="BLOCK",        type=h5_typeInt)
    setFields(3)  = h5_dataField(name="BEZID",        type=h5_typeInt)
    setFields(4)  = h5_dataField(name="BEZ",          type=h5_typeChar, size=mNameLen)
    setFields(5)  = h5_dataField(name="SLOS",         type=h5_typeReal)
    setFields(6)  = h5_dataField(name="X",            type=h5_typeReal)
    setFields(7)  = h5_dataField(name="XP",           type=h5_typeReal)
    setFields(8)  = h5_dataField(name="Y",            type=h5_typeReal)
    setFields(9)  = h5_dataField(name="YP",           type=h5_typeReal)
    setFields(10) = h5_dataField(name="ETOT",         type=h5_typeReal)
    setFields(11) = h5_dataField(name="DE",           type=h5_typeReal)
    setFields(12) = h5_dataField(name="DT",           type=h5_typeReal)
    setFields(13) = h5_dataField(name="ATOMA",        type=h5_typeInt)
    setFields(14) = h5_dataField(name="ATOMZ",        type=h5_typeInt)
#ifdef FLUKA
    setFields(15) = h5_dataField(name="FLUKA_UID",    type=h5_typeInt)
    setFields(16) = h5_dataField(name="FLUKA_GEN",    type=h5_typeInt)
    setFields(17) = h5_dataField(name="FLUKA_WEIGHT", type=h5_typeReal)
#else
    setFields(15) = h5_dataField(name="PARTID",       type=h5_typeInt)
#endif
    call h5_createFormat("aperLostPart", setFields, aper_fmtLostPart)
    call h5_createDataSet("losses", h5_aperID, aper_fmtLostPart, aper_setLostPart)
  else
#endif

#ifdef CR
    if (apefilepos >= 0) then
      ! Expect the file to be opened already, in crcheck
      inquire( unit=losses_unit, opened=isOpen )
      if (.not.isOpen) then
        write(lout,"(2(a,i0),a)") "LIMI> ERROR The unit ",losses_unit," has apefilepos = ", apefilepos, " >= 0, "//&
          "but the file is NOT open. This is probably a bug."
        call prror(-1)
      end if
    else
#endif

    call f_requestUnit(losses_filename,losses_unit)
    inquire(unit=losses_unit, opened=isOpen) ! Was 999
      if(isOpen) then
        write(lout,"(a,i0,a)") "APER> ERROR Unit ",losses_unit," is already open."
        call prror(-1)
      end if

      call f_open(unit=losses_unit,file=losses_filename,formatted=.true.,mode='w',err=err)
#ifdef CR
      apefilepos=0
#endif

      write(losses_unit,"(a)") "# turn block bezid bez slos "// &
#ifdef FLUKA
        "fluka_uid fluka_gen fluka_weight "// &
#else
        "partid "// &
#endif
        "x xp y yp etot dE dT A_atom Z_atom "
      ! Flush file
      flush(losses_unit)
#ifdef CR
      apefilepos=apefilepos+1
    end if
#endif
#ifdef HDF5
  end if
#endif

end subroutine aperture_init

subroutine aperture_nul( ix )
  !-----------------------------------------------------------------------
  ! A.Mereghetti (CERN, BE-ABP-HSS), 2018-03-21
  ! initialise aperture marker to null
  !-----------------------------------------------------------------------
  implicit none
  integer ix, jj
  kape(ix)=0
  do jj=1,9
     ape(jj,ix)=zero
  end do
  lapeofftlt(ix)=.false.
end subroutine aperture_nul


subroutine aperture_initCR( ix, aper )
  !-----------------------------------------------------------------------
  ! A.Mereghetti (CERN, BE-ABP-HSS), 2018-03-21
  ! initialise aperture marker to circle
  !-----------------------------------------------------------------------
  implicit none
  integer ix
  real(kind=fPrec) aper
  call aperture_nul( ix )
  kape(ix)=1
  ape(1,ix)=aper
  ape(2,ix)=aper
  ape(3,ix)=aper
  ape(4,ix)=aper
  ape(5,ix)=-one
  ape(6,ix)=ape(1,ix)*sqrt(two)
end subroutine aperture_initCR


subroutine aperture_initRE( ix, aprx, apry )
  !-----------------------------------------------------------------------
  ! A.Mereghetti (CERN, BE-ABP-HSS), 2018-03-21
  ! initialise aperture marker to rectangle
  !-----------------------------------------------------------------------
  implicit none
  integer ix
  real(kind=fPrec) aprx, apry
  call aperture_nul( ix )
  kape(ix)=2
  ape(1,ix)=aprx
  ape(2,ix)=apry
  ape(3,ix)=aprx*sqrt(two)
  ape(4,ix)=apry*sqrt(two)
  ape(5,ix)=-one
  ape(6,ix)=ape(2,ix)-ape(5,ix)*ape(1,ix)
end subroutine aperture_initRE


subroutine aperture_initEL( ix, apex, apey )
  !-----------------------------------------------------------------------
  ! A.Mereghetti (CERN, BE-ABP-HSS), 2018-03-21
  ! initialise aperture marker to ellipse
  !-----------------------------------------------------------------------
  implicit none
  integer ix
  real(kind=fPrec) apex, apey
  call aperture_nul( ix )
  kape(ix)=3
  ape(1,ix)=apex
  ape(2,ix)=apey
  ape(3,ix)=apex
  ape(4,ix)=apey
  ape(5,ix)=-one
  ape(6,ix)=sqrt(ape(3,ix)**2+ape(4,ix)**2)
end subroutine aperture_initEL


subroutine aperture_initRL( ix, aprx, apry, apex, apey )
  !-----------------------------------------------------------------------
  ! A.Mereghetti (CERN, BE-ABP-HSS), 2018-03-21
  ! initialise aperture marker to rectellipse
  !-----------------------------------------------------------------------
  implicit none
  integer ix
  real(kind=fPrec) aprx, apry, apex, apey
  call aperture_nul( ix )
  kape(ix)=4
  ape(1,ix)=aprx
  ape(2,ix)=apry
  ape(3,ix)=apex
  ape(4,ix)=apey
  ape(5,ix)=-one
  ape(6,ix)=max(ape(2,ix)-ape(5,ix)*ape(1,ix),sqrt(ape(3,ix)**2+ape(4,ix)**2))
end subroutine aperture_initRL


subroutine aperture_initOC( ix, aprx, apry, theta1, theta2 )
  !-----------------------------------------------------------------------
  ! A.Mereghetti (CERN, BE-ABP-HSS), 2018-03-21
  ! initialise aperture marker to octagon
  !-----------------------------------------------------------------------
  use mod_common, only : rad
  implicit none
  integer ix
  real(kind=fPrec) aprx, apry, theta1, theta2, N, x1, x2, y1, y2
  call aperture_nul( ix )
  kape(ix)=5
  x1=aprx
  y1=aprx*tan_mb(theta1*rad) ! Convert degrees to radians
  x2=apry/tan_mb(theta2*rad) ! Convert degrees to radians
  y2=apry
  ape(1,ix)=aprx
  ape(2,ix)=apry
  ! ellipse circumscribed to octagon
  N=(x1*y2+y1*x2)*(x1*y2-y1*x2)        ! N=x1^2*y2^2-y1^2*x2^2
  ape(3,ix)=sqrt(N/((y2+y1)*(y2-y1)))  ! a=sqrt(N/(y2^2-y1^2))
  ape(4,ix)=sqrt(N/((x1+x2)*(x1-x2)))  ! b=sqrt(N/(x1^2-x2^2))
  ! m and q of sloped side
  ape(5,ix)=(y2-y1)/(x2-x1)  ! m = (y2-y1)/(x2-x1)
  ape(6,ix)=y1-ape(5,ix)*x1  ! q = y1 -m*x1
end subroutine aperture_initOC


subroutine aperture_initRT( ix, aprx, apry, radius )
  !-----------------------------------------------------------------------
  ! A.Mereghetti (CERN, BE-ABP-HSS), 2018-03-22
  ! initialise aperture marker to racetrack
  !-----------------------------------------------------------------------
  implicit none
  integer ix
  real(kind=fPrec) aprx, apry, radius
  call aperture_nul( ix )
  kape(ix)=6
  ape(1,ix)=aprx
  ape(2,ix)=apry
  ape(3,ix)=radius
  ape(4,ix)=radius
  ape(5,ix)=-one
  ape(6,ix)=(sqrt(ape(3,ix)**2+ape(4,ix)**2)+(ape(1,ix)-ape(3,ix)))+(ape(2,ix)-ape(4,ix))
end subroutine aperture_initRT


subroutine aperture_initTR( ix, aprx, apry, apex, apey, theta1, theta2 )
  !-----------------------------------------------------------------------
  ! A.Mereghetti (CERN, BE-ABP-HSS), 2018-03-21
  ! initialise aperture marker to transition
  !-----------------------------------------------------------------------
  implicit none
  integer ix
  real(kind=fPrec) aprx, apry, apex, apey, theta1, theta2
  call aperture_nul( ix )
  kape(ix)=5
  ape(1,ix)=aprx
  ape(2,ix)=apry
  ape(3,ix)=apex
  ape(4,ix)=apey
  ! x1=aprx=ape(1,ix)
  ! y1=ape(1,ix)*tan_mb(theta1)
  ! x2=ape(2,ix)/tan_mb(theta2)
  ! y2=apry=ape(2,ix)
  ! m and q of sloped side
  ! m = (y2-y1)/(x2-x1)
  ! q = y1 -m*x1
  ape(5,ix)=(ape(2,ix)-ape(1,ix)*tan_mb(theta1))/(ape(2,ix)/tan_mb(theta2)-ape(1,ix))
  ape(6,ix)=ape(1,ix)*tan_mb(theta1)-ape(5,ix)*ape(1,ix)
end subroutine aperture_initTR


subroutine aperture_initroffpos( ix, xoff, yoff, tilt )
  !-----------------------------------------------------------------------
  ! A.Mereghetti (CERN, BE-ABP-HSS), 2018-03-22
  ! initialise offset/tilt of aperture marker
  !-----------------------------------------------------------------------
  implicit none
  integer ix
  real(kind=fPrec) tilt, xoff, yoff
  ape(7,ix)=tilt
  ape(8,ix)=xoff
  ape(9,ix)=yoff
  lapeofftlt(ix)=ape(7,ix).ne.zero.or.ape(8,ix).ne.zero.or.ape(9,ix).ne.zero
end subroutine aperture_initroffpos


subroutine aperture_saveLastMarker( i, ix )
  !-----------------------------------------------------------------------
  ! A.Mereghetti (CERN, BE-ABP-HSS), 2018-03-07
  ! save last aperture marker
  !-----------------------------------------------------------------------

  implicit none
  ! interface variables
  integer i, ix
  iLast = i
  ixLast = ix
  return

end subroutine aperture_saveLastMarker


subroutine aperture_saveLastCoordinates( i, ix, iBack )
  !-----------------------------------------------------------------------
  ! A.Mereghetti (CERN, BE-ABP-HSS), 2018-03-21
  ! save particle coordinates:
  ! - at last aperture check (iBack=-1)
  ! - at last thick element (iBack>=0)
  !-----------------------------------------------------------------------

  use mod_commonmn ! for napx, xv and yv
  implicit none
  ! interface variables
  integer i, ix, iBack
  ! temporary variables
  integer j

  do j=1,napx
    xLast(1,j) = xv1(j)
    xLast(2,j) = xv2(j)
    yLast(1,j) = yv1(j)
    yLast(2,j) = yv2(j)
    ejfvLast(j) = ejfv(j)
    ejvLast(j) = ejv(j)
    nucmLast(j) = nucm(j)
    sigmvLast(j) = sigmv(j)
    dpsvLast(j) = dpsv(j)
    naaLast(j) = naa(j)
    nzzLast(j) = nzz(j)
  end do

  iLastThick = i
  ixLastThick = ix
  iBckTypeLast = iBack

  return

end subroutine aperture_saveLastCoordinates


subroutine aperture_backTrackingInit
  !-----------------------------------------------------------------------
  ! A.Mereghetti (CERN, BE-ABP-HSS), 2018-03-07
  ! initialise variables for back-tracking
  !-----------------------------------------------------------------------
  use parpro     ! for nblo
  use mod_common ! for ic(i)
  use crcoall    ! for lout
  implicit none
  ! temporary variables
  integer i, ix

  i=1
  ix=ic(i)-nblo
  if( ix.lt.0 ) then
    write(lout,"(a)") "APER> ERROR Impossible to properly initialise backtracking:"
    write(lout,"(a)") "APER>       first element of lattice structure is not a single element"
    call prror(-1)
  end if
  if( kape(ix).eq.0 ) then
    write(lout,"(a)") "APER> ERROR Impossible to properly initialise backtracking:"
    write(lout,"(a)") "APER>       first element of lattice structure is not assigned an aperture profile"
    call prror(-1)
  end if

  call aperture_saveLastCoordinates( i, ix, -1 )
  call aperture_saveLastMarker( i, ix )

  return

end subroutine aperture_backTrackingInit

subroutine aperture_checkApeMarker(turn, i, ix, llost)
!-----------------------------------------------------------------------
!     P.Garcia Ortega, A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified:  8-12-2014
!     aperture check and dump lost particles
!     always in main code
!-----------------------------------------------------------------------
!     7 April 2014
!-----------------------------------------------------------------------

  use physical_constants

#ifdef FLUKA
  use mod_fluka
#endif

#ifdef ROOT
  use iso_c_binding
  use root_output
#endif

  use collimation, only : do_coll, part_abs_turn, ipart

  implicit none

! parameters
  integer turn  ! turn number
  integer i     ! element entry in the lattice
  integer ix    ! single element type index
  logical llost ! at least a particle loss

  integer ib2,ib3,ilostch,j,jj,jj1,jjx

! temporary variables
  logical lparID
  real(kind=fPrec) apxx, apyy, apxy, aps, apc, radius2
  real(kind=fPrec) xchk(2)

#ifdef ROOT
  character(len=mNameLen+1) this_name
#endif

! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
! last modified: 12-06-2014
! additional variables for back-tracking, when computing locations of
! lost particles
! inserted in main code by the 'backtrk' compilation flag
  integer niter       ! number of iterations
  integer kapert      ! temporal integer for aperture type
  logical llos        ! temporal logic array for interpolation
  real(kind=fPrec) xlos(2), ylos(2), aprr(9), step, length, slos, ejfvlos, ejvlos, nucmlos, sigmvlos, dpsvlos
  integer naalos, nzzlos

  integer npart_tmp ! Temporary holder for number of particles,
                    ! used to switch between collimat/standard version at runtime

  save

  !-----------------------------------------------------------------------
  ! check against current aperture marker
  !-----------------------------------------------------------------------

  ! go through all possible types
  select case(kape(ix))

  case (-1) ! Transition
    apxx = ape(3,ix)**2.
    apyy = ape(4,ix)**2.
    apxy = apxx * apyy
    do j=1,napx
      if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll)) then
        if(lapeofftlt(ix)) then
          if(lbacktracking) then
            call roffpos(xLast(1,j),xLast(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
          else
            call roffpos(xv1(j),xv2(j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
          end if
          llostp(j)=checkTR(xchk(1),xchk(2),ape(1,ix),ape(2,ix),ape(3,ix),ape(4,ix),apxx,apyy,apxy,ape(5,ix),ape(6,ix)).or. &
            isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
        else
          if(lbacktracking) then
            llostp(j)= &
              checkTR(xLast(1,j),xLast(2,j),ape(1,ix),ape(2,ix),ape(3,ix),ape(4,ix),apxx,apyy,apxy,ape(5,ix),ape(6,ix)) .or. &
              isnan_mb(xLast(1,j)).or.isnan_mb(xLast(2,j))
          else
            llostp(j)= &
              checkTR(xv1(j),xv2(j),ape(1,ix),ape(2,ix),ape(3,ix),ape(4,ix),apxx,apyy,apxy,ape(5,ix),ape(6,ix))       .or. &
              isnan_mb(xv1(j)).or.isnan_mb(xv2(j))
          end if
        end if
        llost=llost.or.llostp(j)
      end if
    end do

  case (1) ! circle
    radius2 = ape(3,ix)**2
    do j=1,napx

      if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll)) then

        if(lapeofftlt(ix)) then
          if(lbacktracking) then
            call roffpos(xLast(1,j),xLast(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
          else
            call roffpos(xv1(j),xv2(j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
          end if
          llostp(j)=checkCR( xchk(1),xchk(2),radius2 ) .or. &
            isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
        else
          if(lbacktracking) then
            llostp(j)=checkCR( xLast(1,j),xLast(2,j),radius2 ) .or. &
              isnan_mb(xLast(1,j)).or.isnan_mb(xLast(2,j))
          else
            llostp(j)=checkCR( xv1(j),xv2(j),radius2 ) .or. &
              isnan_mb(xv1(j)).or.isnan_mb(xv2(j))
          end if
        end if
      end if
    end do

  case (2) ! Rectangle
    do j=1,napx
      if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll) ) then
        if(lapeofftlt(ix)) then
          if(lbacktracking) then
            call roffpos(xLast(1,j),xLast(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
          else
            call roffpos(xv1(j),xv2(j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
          end if
          llostp(j)=checkRE( xchk(1),xchk(2),ape(1,ix),ape(2,ix) ) .or. &
            isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
        else
          if(lbacktracking) then
            llostp(j)=checkRE( xLast(1,j),xLast(2,j),ape(1,ix),ape(2,ix) ) .or. &
              isnan_mb(xLast(1,j)).or.isnan_mb(xLast(2,j))
          else
            llostp(j)=checkRE( xv1(j),xv2(j),ape(1,ix),ape(2,ix) ) .or. &
              isnan_mb(xv1(j)).or.isnan_mb(xv2(j))
          end if
        end if
        llost=llost.or.llostp(j)
      end if
    end do

  case (3) ! Ellipse
    apxx = ape(3,ix)**2.
    apyy = ape(4,ix)**2.
    apxy = apxx * apyy
    do j=1,napx
      if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll)) then
        if(lapeofftlt(ix)) then
          if(lbacktracking) then
            call roffpos(xLast(1,j),xLast(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
          else
            call roffpos(xv1(j),xv2(j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
          end if
          llostp(j)=checkEL( xchk(1),xchk(2),apxx,apyy,apxy ) .or. &
            isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
        else
          if(lbacktracking) then
            llostp(j)=checkEL( xLast(1,j),xLast(2,j),apxx,apyy,apxy ) .or. &
              isnan_mb(xLast(1,j)).or.isnan_mb(xLast(2,j))
          else
            llostp(j)=checkEL( xv1(j),xv2(j),apxx,apyy,apxy ) .or. &
              isnan_mb(xv1(j)).or.isnan_mb(xv2(j))
          end if
        end if
        llost=llost.or.llostp(j)
      end if
    end do

  case (4) ! RectEllipse
    apxx = ape(3,ix)**2.
    apyy = ape(4,ix)**2.
    apxy = apxx * apyy
    do j=1,napx
      if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll)) then
        if(lapeofftlt(ix)) then
          if(lbacktracking) then
            call roffpos(xLast(1,j),xLast(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
          else
            call roffpos(xv1(j),xv2(j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
          end if
          llostp(j)=checkRL( xchk(1),xchk(2),ape(1,ix),ape(2,ix),apxx,apyy,apxy ) .or. &
            isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
        else
          if(lbacktracking) then
            llostp(j)=checkRL( xLast(1,j),xLast(2,j),ape(1,ix),ape(2,ix),apxx,apyy,apxy ) .or. &
              isnan_mb(xLast(1,j)).or.isnan_mb(xLast(2,j))
          else
            llostp(j)=checkRL( xv1(j),xv2(j),ape(1,ix),ape(2,ix),apxx,apyy,apxy ) .or. &
              isnan_mb(xv1(j)).or.isnan_mb(xv2(j))
          end if
        end if
        llost=llost.or.llostp(j)
      end if
    end do

  case (5) ! Octagon
    do j=1,napx
      if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll)) then
        if(lapeofftlt(ix)) then
          if(lbacktracking) then
            call roffpos(xLast(1,j),xLast(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
          else
            call roffpos(xv1(j),xv2(j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
          end if
          llostp(j)=checkOC(xchk(1),xchk(2),ape(1,ix),ape(2,ix),ape(5,ix),ape(6,ix)).or. &
            isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
        else
          if(lbacktracking) then
            llostp(j)=checkOC(xLast(1,j),xLast(2,j),ape(1,ix),ape(2,ix),ape(5,ix),ape(6,ix)).or. &
              isnan_mb(xLast(1,j)).or.isnan_mb(xLast(2,j))
          else
            llostp(j)=checkOC(xv1(j),xv2(j),ape(1,ix),ape(2,ix),ape(5,ix),ape(6,ix)).or. &
              isnan_mb(xv1(j)).or.isnan_mb(xv2(j))
          end if
        end if
        llost=llost.or.llostp(j)
      end if
    end do

  case (6) ! Racetrack
    !   NB: it follows the MadX definition
    apxy = ape(3,ix)**2.
    do j=1,napx
      if((do_coll .and. part_abs_turn(j).eq.0) .or. (.not.do_coll)) then
        if(lapeofftlt(ix)) then
          if(lbacktracking) then
            call roffpos(xLast(1,j),xLast(2,j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
          else
            call roffpos(xv1(j),xv2(j),xchk(1),xchk(2),ape(7,ix),ape(8,ix),ape(9,ix))
          end if
          llostp(j)=checkRT(xchk(1),xchk(2),ape(1,ix),ape(2,ix),ape(3,ix),apxy).or. &
            isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
        else
          if(lbacktracking) then
            llostp(j)=checkRT(xLast(1,j),xLast(2,j),ape(1,ix),ape(2,ix),ape(3,ix),apxy).or. &
              isnan_mb(xLast(1,j)).or.isnan_mb(xLast(2,j))
          else
            llostp(j)=checkRT(xv1(j),xv2(j),ape(1,ix),ape(2,ix),ape(3,ix),apxy).or. &
              isnan_mb(xv1(j)).or.isnan_mb(xv2(j))
          end if
        end if
        llost=llost.or.llostp(j)
      end if
    end do

  end select

end subroutine aperture_checkApeMarker


subroutine aperture_reportLoss(turn, i, ix)
!-----------------------------------------------------------------------
!     P.Garcia Ortega, A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified:  8-12-2014
!     aperture check and dump lost particles
!     always in main code
!-----------------------------------------------------------------------
!     7 April 2014
!-----------------------------------------------------------------------

  use physical_constants

#ifdef FLUKA
  use mod_fluka
#endif
#ifdef HDF5
  use hdf5_output
#endif
#ifdef ROOT
  use iso_c_binding
  use root_output
#endif

  use collimation, only : do_coll, part_abs_turn, ipart

  implicit none

! parameters
  integer turn  ! turn number
  integer i     ! element entry in the lattice
  integer ix    ! single element type index

  integer ib2,ib3,ilostch,j,jj,jj1,jjx

! temporary variables
  logical lparID
  real(kind=fPrec) apxx, apyy, apxy, aps, apc, radius2
  real(kind=fPrec) xchk(2)

#ifdef ROOT
  character(len=mNameLen+1) this_name
#endif

! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
! last modified: 12-06-2014
! additional variables for back-tracking, when computing locations of
! lost particles
! inserted in main code by the 'backtrk' compilation flag
  integer niter       ! number of iterations
  integer kapert      ! temporal integer for aperture type
  logical llos        ! temporal logic array for interpolation
  logical lback       ! actually perform backtracking
  real(kind=fPrec) xlos(2), ylos(2), aprr(9), step, length, slos, ejfvlos, ejvlos, nucmlos, sigmvlos, dpsvlos
  integer naalos, nzzlos

  integer npart_tmp ! Temporary holder for number of particles,
                    ! used to switch between collimat/standard version at runtime

  save

  lback=.false.

  !-----------------------------------------------------------------------
  ! dump coordinates in case of losses
  ! if back-tracking is requested, get more detailed point of loss
  ! for the moment, only bi-section method
  !-----------------------------------------------------------------------

  if(lbacktracking.and.kape(ix).ne.0.and.iBckTypeLast.ge.0) then
    lback=.true.

    ! Length between elements
    length = dcum(i) - dcum(iLast)

    ! - pay attention to overflow:
    if( length .lt. zero ) then
      length = length+dcum(iu)
    end if

    ! - pay attention to too short thick elements
    if( length .le. bktpre ) then
      lback=.false.
    end if

  end if

  ! Number of iterations for bisection method (ln(2x/precision)/ln(2)+1)
  if(lback) then
    niter=nint(inv_ln2*log_mb((two*length)/bktpre)+2)
  end if

  do j=1,napx
    if(llostp(j)) then
      ! treat a lost particle

      ! ==============================================================
      ! point of loss
      if(lback) then
        ! A. Mereghetti and P. Garcia Ortega, for the FLUKA Team
        ! last modified: 21-03-2018
        ! back-track particles, in order to better estimate actual loss point

        ylos(1)=yLast(1,j)
        ylos(2)=yLast(2,j)

        ! actual algorithm
        llos = llostp(j)
        step = one

        do jj=1,niter
          ! current step (bisection method):
          if( llos ) then
            step = step - one / (two**(jj))
          else
            step = step + one / (two**(jj))
          end if

          ! - step discretized if last iteration, to compare with BeamLossPattern
          if(jj.eq.niter) then
            slos = int((dcum(iLast)+length*step)/bktpre+one)*bktpre
            step = (slos-dcum(iLast))/length
          end if

          ! - particle coordinates at current step
          select case(iBckTypeLast)
          case (0)
            ! back-track along a drift
            xlos(1) = xLast(1,j)  - yLast(1,j)*((one-step)*length)
            xlos(2) = xLast(2,j)  - yLast(2,j)*((one-step)*length)
            slos    = dcum(iLast) + (step*length)
          end select

          ! - aperture at current step
          call interp_aperture( iLast, ixLast, i, ix, kapert, aprr, slos )

          ! Check aperture
          if( lapeofftlt(ix).or.lapeofftlt(ixLast) ) then
            call roffpos( xlos(1), xlos(2), xchk(1),xchk(2), aprr(7), aprr(8), aprr(9) )
          else
            xchk(1) = xlos(1)
            xchk(2) = xlos(2)
          end if

          select case(kapert)
          case(-1) ! Transition
            apxx = aprr(3)**2.
            apyy = aprr(4)**2.
            apxy = apxx * apyy
            llos=checkTR(xchk(1),xchk(2),aprr(1),aprr(2),aprr(3),aprr(4),apxx,apyy,apxy,aprr(5),aprr(6)).or. &
              isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
          case (1) ! Circle
            radius2 = aprr(3)**2
            llos=checkCR(xchk(1),xchk(2),radius2) .or. &
              isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
          case (2) ! Rectangle
            llos=checkRE(xchk(1),xchk(2),aprr(1),aprr(2)) .or. &
              isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
          case (3) ! Ellipse
            apxx = aprr(3)**2.
            apyy = aprr(4)**2.
            apxy = apxx * apyy
            llos=checkEL( xchk(1),xchk(2),apxx,apyy,apxy )  .or. &
              isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
          case (4) ! RectEllipse
            apxx = aprr(3)**2.
            apyy = aprr(4)**2.
            apxy = apxx * apyy
            llos = checkRL( xchk(1),xchk(2),aprr(1),aprr(2),apxx, apyy, apxy ) .or. &
              isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
          case (5) ! Octagon
            llos=checkOC(xchk(1), xchk(2), aprr(1), aprr(2), aprr(5), aprr(6) ) .or. &
              isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
          case (6) ! RaceTrack
            llos=checkRT( xchk(1), xchk(2), aprr(1), aprr(2), aprr(3), aprr(3)**2. ) .or. &
              isnan_mb(xchk(1)).or.isnan_mb(xchk(2))
          end select
        end do !do jj=1,niter

        ! pay attention to overflow
        if( slos.gt.dcum(iu) ) then
          slos=slos-dcum(iu)
        end if

      else !if(lback)
        if(lbacktracking) then
          xlos(1) = xLast(1,j)
          xlos(2) = xLast(2,j)
          ylos(1) = yLast(1,j)
          ylos(2) = yLast(2,j)
          slos    = dcum(iLastThick)
        else
          xlos(1) = xv1(j)
          xlos(2) = xv2(j)
          ylos(1) = yv1(j)
          ylos(2) = yv2(j)
          slos    = dcum(i)
        end if
      end if ! if(lback)

      ! get ready for dumping infos
      if(lbacktracking) then
        ejfvlos = ejfvLast(j)
        ejvlos = ejvLast(j)
        nucmlos = nucmLast(j)
        sigmvlos = sigmv(j)
        dpsvlos = dpsvLast(j)
        naalos = naaLast(j)
        nzzlos = nzzLast(j)
      else
        ejfvlos = ejfv(j)
        ejvlos = ejv(j)
        nucmlos = nucm(j)
        sigmvlos = sigmv(j)
        dpsvlos = dpsv(j)
        naalos = naa(j)
        nzzlos = nzz(j)
      end if

      ! ==============================================================
      ! If lost particles aren't killed, the lost info is dumped only
      ! the first time they hit the aperture. Their secondaries generated
      ! from a lost particles are considered lost as well
      if( apflag ) then
        lparID = .false.
        jjx=1

        !TODO is this really needed?
        if (do_coll) then
          npart_tmp = npart
        else
          npart_tmp = napx
        endif

        do jj=1,npart_tmp
          if(plost(jj).ne.0) then
#ifdef FLUKA
            if( fluka_uid(j).eq.plost(jj).or. fluka_gen(j).eq.plost(jj) ) then
              lparID=.true.
            end if
#else
            if ( (     do_coll .and. (  ipart(j) .eq. plost(jj) )) .or. &
                 (.not.do_coll .and. ( partID(j) .eq. plost(jj) ))       ) then
              lparID=.true.
            end if
#endif

            jjx=jj+1 !points to the last zero
          end if
        end do

        if(lparID) then
          !old lost particle or secondary, don't print it
          goto 1982
        else
          !new lost particle, store ID and print it
#ifdef FLUKA
          plost(jjx) = fluka_uid(j)
#else
          if (do_coll) then
            plost(jjx) = ipart(j)
          else
            plost(jjx) = j
          endif
#endif
        end if !if(lparID) then
      end if !if( apflag ) then

#ifdef HDF5
      if(h5_useForAPER) then
        call h5_prepareWrite(aper_setLostPart, 1)
        call h5_writeData(aper_setLostPart, 1,  1, turn)
        call h5_writeData(aper_setLostPart, 2,  1, i)
        call h5_writeData(aper_setLostPart, 3,  1, ix)
        call h5_writeData(aper_setLostPart, 4,  1, bez(ix))
        call h5_writeData(aper_setLostPart, 5,  1, slos)
        call h5_writeData(aper_setLostPart, 6,  1, xlos(1)*c1m3)
        call h5_writeData(aper_setLostPart, 7,  1, xlos(2)*c1m3)
        call h5_writeData(aper_setLostPart, 8,  1, ylos(1)*c1m3)
        call h5_writeData(aper_setLostPart, 9,  1, ylos(2)*c1m3)
        call h5_writeData(aper_setLostPart, 10, 1, ejfvlos*c1m3)
        call h5_writeData(aper_setLostPart, 11, 1, (ejvlos*(nucm0/nucmlos)-e0)*c1e6)
        call h5_writeData(aper_setLostPart, 12, 1, (-(c1m3 * (sigmvlos/clight) ))* (e0/e0f))
        call h5_writeData(aper_setLostPart, 13, 1, naalos)
        call h5_writeData(aper_setLostPart, 14, 1, nzzlos)
#ifdef FLUKA
        call h5_writeData(aper_setLostPart, 15, 1, fluka_uid(j))
        call h5_writeData(aper_setLostPart, 16, 1, fluka_gen(j))
        call h5_writeData(aper_setLostPart, 17, 1, fluka_weight(j))
#endif
        if (do_coll) then
          call h5_writeData(aper_setLostPart, 15, 1, ipart(j))
        endif
#ifndef FLUKA
        if (.not. do_coll) then
          call h5_writeData(aper_setLostPart, 15, 1, partID(j))
        endif
#endif
        call h5_finaliseWrite(aper_setLostPart)
      else
  ! END of #ifdef HDF5
#endif

        ! Print to unit 999 (fort.999)
#ifdef FLUKA
        write(losses_unit,'(3(1X,I8),1X,A48,1X,F12.5,2(1X,I8),8(1X,1PE14.7),2(1X,I8))')&
#else
        write(losses_unit,'(3(1X,I8),1X,A48,1X,F12.5,1X,I8,7(1X,1PE14.7),2(1X,I8))')   &
#endif

     &       turn, i, ix, bez(ix), slos,                                     &
#ifdef FLUKA
     &       fluka_uid(j), fluka_gen(j), fluka_weight(j),                    &
#else
     &       partID(j),                                                      &
#endif

     &       xlos(1)*c1m3, ylos(1)*c1m3, xlos(2)*c1m3, ylos(2)*c1m3,         &
     &       ejfvlos*c1m3, (ejvlos*(nucm0/nucmlos)-e0)*c1e6,                 &
     &       (-(c1m3 * (sigmvlos/clight) ))* (e0/e0f),                       &
     &       naalos, nzzlos
#ifdef CR
        apefilepos=apefilepos+1
#endif
#ifdef HDF5
      end if
#endif

#if defined(ROOT)
! root output
      if(root_flag .and. root_ApertureCheck.eq.1) then
        this_name = trim(adjustl(bez(ix))) // C_NULL_CHAR
#if defined(FLUKA)
        call ApertureCheckWriteLossParticleF(turn, i, ix, this_name, len_trim(this_name), slos, &
          fluka_uid(j), fluka_gen(j), fluka_weight(j), &
          xlos(1)*c1m3, ylos(1)*c1m3, xlos(2)*c1m3, ylos(2)*c1m3, ejfvlos*c1m3, (ejvlos-e0)*c1e6, &
          (-(c1m3 * (sigmvlos/clight))) * (e0/e0f), naalos, nzzlos)
#else
        call ApertureCheckWriteLossParticle(turn, i, ix, this_name, len_trim(this_name), slos, plost(j),&
          xlos(1)*c1m3, ylos(1)*c1m3, xlos(2)*c1m3, ylos(2)*c1m3, ejfvlos*c1m3, (ejvlos-e0)*c1e6, &
          (-(c1m3 * (sigmvlos/clight))) * (e0/e0f), naalos, nzzlos)
#endif
      end if
#endif

#ifdef FLUKA
      if(((partID(j).le.aperture_napxStart) .and. fluka_enable) .or. .not.fluka_enable) then
#else
      if(((partID(j).le.aperture_napxStart) .and. do_coll) .or. .not.do_coll) then
#endif
        pstop(partID(j))=.true.
        ! Record for postpr
        if(.not.limifound.or.kape(ix).eq.0) then
          aperv(partID(j),1) = aper(1)
          aperv(partID(j),2) = aper(2)
        else
          aperv(partID(j),1) = min(ape(1,ix),ape(3,ix))
          aperv(partID(j),2) = min(ape(2,ix),ape(4,ix))
        end if
        ixv(partID(j))     = ix
        xvl(1,partID(j))   = xlos(1)
        xvl(2,partID(j))   = xlos(2)
        yvl(1,partID(j))   = ylos(1)
        yvl(2,partID(j))   = ylos(2)
        dpsvl(partID(j))   = dpsvlos
        ejvl(partID(j))    = ejvlos
        sigmvl(partID(j))  = sigmvlos
        numxv(partID(j))   = numx
        nnumxv(partID(j))  = numx
#ifdef FLUKA
      end if ! partID(j).le.aperture_napxStart
#else
      end if ! (partID(j).le.aperture_napxStart .and. do_coll) .or. .not.do_coll
#endif

1982  continue

    end if ! if(llostp(j))
  end do ! do j=1,napx

  ! flush loss particle file
#ifdef HDF5
  if(.not. h5_useForAPER) then
#endif
     flush(losses_unit)
#ifdef HDF5
  end if
#endif

end subroutine aperture_reportLoss


logical function checkRE( x, y, apex, apey )
!-----------------------------------------------------------------------
!     A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
!     last modified: 16-05-2014
!     check particle position against REctangle aperture
!     always in main code
!-----------------------------------------------------------------------
  implicit none
! parameters
  real(kind=fPrec) x, y, apex, apey
  checkRE = ( abs(x).gt.apex ).or.( abs(y).gt.apey )
  return
end function

logical function checkEL( x, y, apxx, apyy, apxy )
!-----------------------------------------------------------------------
!     A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
!     last modified: 16-05-2014
!     check particle position against ELlipse aperture
!     always in main code
!-----------------------------------------------------------------------
  implicit none

! parameters
  real(kind=fPrec) x, y, apxx, apyy, apxy

  checkEL = x**two*apyy+y**two*apxx .gt. apxy
  return
end function checkEL

logical function checkRL( x, y, apex, apey, apxx, apyy, apxy )
!-----------------------------------------------------------------------
!     A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
!     last modified: 16-05-2014
!     check particle position against Rect-Ellipse aperture
!     always in main code
!-----------------------------------------------------------------------
  implicit none

! parameters
  real(kind=fPrec) x, y, apex, apey, apxx, apyy, apxy

  checkRL = checkRE( x, y, apex, apey ) .or. checkEL( x, y, apxx, apyy, apxy )
  return
end function checkRL

logical function checkOC( x, y, ap1, ap2, m, q )
!-----------------------------------------------------------------------
!     A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
!     last modified: 22-03-2018
!     check particle position against OCtagon aperture
!     always in main code
!-----------------------------------------------------------------------
  implicit none

! parameters
  real(kind=fPrec) x, y, ap1, ap2, m, q

  checkOC = checkRE(x,y,ap1,ap2).or.(abs(y).gt.m*abs(x)+q)
  return
end function checkOC

logical function checkRT( x, y, apex, apey, r, r2 )
!-----------------------------------------------------------------------
!     A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
!     last modified: 19-05-2014
!     check particle position against RaceTrack aperture
!     always in main code
!-----------------------------------------------------------------------
  implicit none

! parameters
  real(kind=fPrec) x, y, apex, apey, r, r2

  checkRT = checkRE( x, y, apex+r, apey+r ) .or. ( ( (abs(x)-apex)**2.+(abs(y)-apey)**2.).gt.r2 )
  return
end function checkRT

logical function checkCR( x, y, radius2 )
!-----------------------------------------------------------------------
!     check particle position against CiRcle aperture
!     always in main code
!-----------------------------------------------------------------------
  implicit none

! parameters
  real(kind=fPrec) x, y, radius2

  checkCR = (x**2 + y**2) .gt. radius2
  return
end function checkCR

logical function checkTR( x, y, aprx, apry, apex, apey, apxx, apyy, apxy, m, q )
!-----------------------------------------------------------------------
!     A.Mereghetti (CERN, BE/ABP-HSS)
!     last modified: 22-03-2018
!     check particle position against Transition aperture
!     always in main code
!-----------------------------------------------------------------------
  implicit none

! parameters
  real(kind=fPrec) x, y, aprx, apry, apex, apey, apxx, apyy, apxy, m, q

  checkTR = checkRL(x,y,aprx,apry,apxx,apyy,apxy).or.checkOC(x,y,aprx,apry,m,q)
  if(aprx-apex.gt.zero.and.apry-apey.gt.zero) checkTR=checkTR.or.checkRT(x,y,aprx,apry,apex,apxx)
  return
end function checkTR

subroutine roffpos( x, y, xnew, ynew, tlt, xoff, yoff )
!-----------------------------------------------------------------------
!     A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
!     last modified: 16-05-2014
!     centre/rotate position of particles in case of offcentered/tilted
!        aperture types
!     always in main code
!
!     input parameters:
!        x : horizontal particle position [mm]
!        y : vertical   particle position [mm]
!        tlt:  tilt angle of the aperture profile [rad]
!        xoff: horizontal aperture offset [mm]
!        yoff: vertical   aperture offset [mm]
!
!     output parameters:
!        xnew : offcentered/tilted horizontal particle position [mm]
!        ynew : offcentered/tilted vertical   particle position [mm]
!
!-----------------------------------------------------------------------
  implicit none

! parameters
  real(kind=fPrec) x, y, xnew, ynew, tlt, xoff, yoff

! temporary variables
  real(kind=fPrec) theta, radio, xtmp, ytmp, ttmp

  xtmp = x-xoff
  ytmp = y-yoff
  theta = atan2_mb(ytmp, xtmp)
  radio = sqrt(xtmp**two + ytmp**two)
  ttmp = theta-tlt
  xnew = radio * cos_mb(ttmp)
  ynew = radio * sin_mb(ttmp)
  return
end subroutine roffpos

subroutine roffpos_inv( x, y, xnew, ynew, tlt, xoff, yoff )
!-----------------------------------------------------------------------
!     A.Mereghetti (CERN, BE/ABP-HSS), 2018-03-24
!     inverse of roffpos - but same interface
!-----------------------------------------------------------------------
  implicit none

! parameters
  real(kind=fPrec) x, y, xnew, ynew, tlt, xoff, yoff

! temporary variables
  real(kind=fPrec) theta, radio, ttmp

  theta = atan2_mb(y, x)
  radio = sqrt(x**two + y**two)
  ttmp = theta+tlt
  xnew = radio * cos_mb(ttmp)
  ynew = radio * sin_mb(ttmp)
  xnew = xnew+xoff
  ynew = ynew+yoff
  return
end subroutine roffpos_inv

subroutine contour_aperture_markers( itElUp, itElDw, lInsUp )
!-----------------------------------------------------------------------
! by A.Mereghetti
! last modified: 20-12-2016
! always in main code
! check elements itElUp (upstream) and itElDw (downstream) and
!   assign them (or insert) an aperture marker, in case;
! lInsUp: force the insertion of an aperture marker upstream
!   of itElUp (.false. in case of initialisation of aperture profile
!   of entire lattice, as the upstream element is the actual entry
!   in lattice sequence);
!-----------------------------------------------------------------------
  implicit none

! interface variables
  integer itElUp, itElDw
  logical lInsUp
! run time variables
  integer iElUp, iElDw, ixApeUp, ixApeDw, jj, iuold
  logical lExtremes, lsame

! do not overwrite interface variables
  iElUp=itElUp
  iElDw=itElDw
! handling extremes of lattice structure?
  lExtremes=iElUp.eq.iu.and.iElDw.eq.1

! upstream marker
  iuold=iu
  call contour_aperture_marker( iElUp, lInsUp )
! the addition of the upstream aperture marker may have
!    shifted by one the downstream entries
! NB: if lExtremes, the upstream marker is the last entry
!     in the lattice structure! Hence, no other entry is shifted!
  if( .not.lExtremes ) then
    if( iu-iuold.ne.0 ) then
      iElDw=iElDw+(iu-iuold)
      write(lout,"(a,i0)") "APER> ...inserted upstream marker - downstream entries shifted by ",iu-iuold
    else
      write(lout,"(a)")    "APER> ...no need to insert an upstream marker - no shift of downstream entries required."
    end if
  end if

! downstream marker
  iuold=iu
  call contour_aperture_marker( iElDw, .false. )
! the addition of the downstream aperture marker may have shifted by one the downstream entries
  if( iu-iuold.ne.0 ) then
! NB: if lExtremes, the downstream entry is the first entry
! in the lattice structure! Hence, if a new entry has been inserted,
! the upstream entry (at the end of the lattice structure) is
! shifted by 1
    if( lExtremes ) then
      iElUp=iElUp+(iu-iuold)
    end if
    write(lout,"(a,i0)") "APER> ...inserted downstream marker - downstream entries shifted by ",iu-iuold
  else
    write(lout,"(a)")    "APER> ...no need to insert a downstream marker - no shift of downstream entries required."
  end if

  if( lExtremes ) then
! check that the aperture markers at the extremities of accelerator
! lattice structure are the same
    ixApeUp=ic(iElUp)-nblo
    ixApeDw=ic(iElDw)-nblo
    lsame = sameAperture(ixApeUp,ixApeDw)
    if( .not.lsame ) then
      write(lout,"(a)") "APER> ERROR Different aperture markers at extremeties of accelerator lattice strucure"
      call dump_aperture_header( lout )
      call dump_aperture_marker( lout, ixApeUp, iElUp )
      call dump_aperture_marker( lout, ixApeDw, iElDw )
      call prror(-1)
    end if
  end if

end subroutine contour_aperture_markers

subroutine contour_aperture_marker( iEl, lInsUp )
!-----------------------------------------------------------------------
!     by A.Mereghetti
!     last modified: 20-12-2016
!     put an aperture marker at iEl
!     NB: it can be either a brand new entry in lattice sequence or
!         updating an existing one
!     interface variables:
!     - iEl: entry in lattice sequence to be checked
!     - lInsUp: if true, the new aperture marker is inserted upstream of iEl
!     always in main code
!-----------------------------------------------------------------------
#ifdef FLUKA
! import mod_fluka
! inserted in main code by the 'fluka' compilation flag
  use mod_fluka
#endif

  implicit none

! interface variables
  integer iEl
  logical lInsUp
! temporary variables
  integer i,ix,iSrcUp,iSrcDw,iApeUp,ixApeUp,iApeDw,ixApeDw,jj,itmpape,iNew,ixNew,check_SE_unique,INEESE,INEELS,ixApeNewFrom,ixEl
  real(kind=fPrec) tmpape(9), ddcum
  logical lconst,lApeUp,lApeDw,lAupDcum,lAdwDcum,lApe,lAss,lfit

! echo of input parameters
  write(lout,"(a)") ""
  write(lout,"(a)") "APER> Call to contour_aperture_marker"

! check upstream element
  ixEl=ic(iEl)-nblo
  if( iEl.eq.iu ) then
! end of lattice sequence: a marker might be needed
    if( ixEl.le.0 ) then
      ix=INEESE()
      iu=INEELS( 0 )
      ic(iu)=ix+nblo
      iEl=iu
      ixEl=ix
      bez(ixEl)='e.latt.aper'
      write(lout,"(a)") "APER> -> Inserted empty marker at end of lattice"
    end if
  else if( iEl.eq.1 ) then
! beginning of lattice sequence: a marker might be needed
    if( ixEl.le.0 ) then
      ix=INEESE()
      iu=INEELS( 1 )
      ic(1)=ix+nblo
      iEl=1
      ixEl=ix
      bez(ixEl)='s.latt.aper'
      write(lout,"(a)") "APER> -> Inserted empty marker at start of lattice"
#ifdef FLUKA
    else if( fluka_type(ixEl).eq.FLUKA_ELEMENT.or.fluka_type(ixEl).eq.FLUKA_ENTRY   ) then
! A.Mereghetti
! last modified: 18-01-2017
! force aperture marker upstream of FLUKA_ENTRY
! inserted in main code by the 'fluka' compilation flag
      ix=INEESE()
      iu=INEELS( 1 )
      ic(1)=ix+nblo
      iEl=1
      ixEl=ix
      bez(ixEl)='s.latt.aper'
      write(lout,"(a)") "APER> -> Inserted empty marker at start of lattice since first entry is a FLUKA element"
#endif
    end if
  else if( ixEl.le.0 ) then
    write(lout,"(a,i0,a)") "APER> ERROR Lattice element at: i=",iEl," is NOT a SINGLE ELEMENT."
    call prror(-1)
  end if

! echo
  write(lout,"(a)")                 "APER> Look for aperture markers closest to:"
  write(lout,"(a,i0,a,i0,a,e15.7)") "APER> i=",iEl," - ix=",ixEl," - name: '"//bez(ixEl)//"' - s=",dcum(iEl)

! candidate aperture marker
  if( lInsUp ) then
    iNew=iEl-1
  else
    iNew=iEl
  end if

  ixNew=ic(iNew)-nblo
  if( iEl.eq.iu ) then
!   end of lattice sequence
    iSrcUp=iNew
    iSrcDw=1
  else if( iEl.eq.1 ) then
!   beginning of lattice sequence:
    iSrcUp=iu
    iSrcDw=iEl
  else
    iSrcUp=iNew
    iSrcDw=iEl
  end if

! - get closest upstream aperture marker
! NB: no risk of overflow, as first/last element in lattice
! sequence should be aperture markers (and the first
! call of this function is meant to verify this assumption)
  call find_closest_aperture(iSrcUp,.true.,iApeUp,ixApeUp,lApeUp)
  if( iApeUp.eq.-1 .and. ixApeUp.eq.-1 ) then
    write(lout,"(a)") "APER> ERROR Could not find upstream marker"
    call prror(-1)
  end if
! - get closest downstream aperture marker
! NB: no risk of overflow, as first/last element in lattice
! sequence should be aperture markers (and the first
! call of this function is meant to verify this assumption)
  call find_closest_aperture(iSrcDw,.false.,iApeDw,ixApeDw,lApeDw)
  if( iApeDw.eq.-1 .and. ixApeDw.eq.-1 ) then
    write(lout,"(a)") "APER> ERROR Could not find downstream marker"
    call prror(-1)
  end if
! - echo found apertures
  call dump_aperture_header( lout )
  call dump_aperture_marker( lout, ixApeUp, iApeUp )
  call dump_aperture_marker( lout, ixApeDw, iApeDw )

! - checks:
! . iNew is iApeUp
  lApeUp=iApeUp.eq.iNew.and.ixApeUp.eq.ixNew

! . iNew is at the same s as iApeUp (inlcuding ring overflow)
  lAupDcum=abs(dcum(iNew)-dcum(iApeUp)).lt.sPrec.or.abs(dcum(iNew)-dcum(iApeUp)-dcum(iu)).lt.sPrec

! . iNew is iApeDw
  lApeDw=iApeDw.eq.iNew.and.ixApeDw.eq.ixNew

! . iNew is at the same s as ApeDw (inlcuding ring overflow)
  lAdwDcum=abs(dcum(iNew)-dcum(iApeDw)).lt.sPrec.or.abs(dcum(iNew)-dcum(iApeDw)-dcum(iu)).lt.sPrec

! . constant aperture?
  lconst = sameAperture( ixApeUp, ixApeDw )

! . can iNew be assigned an aperture marker?
! ie is it a single element and is it used anywhere else?
  lApe=lApeUp.or.lApeDw
  lAss=ixNew.gt.0.and.check_SE_unique(iNew,ixNew).eq.-1

! some action is needed
  if( .not.lApe ) then
! . iNew must be assigned an aperture
    ixApeNewFrom=-1
    lfit=.false.
    itmpape=0
    do jj=1,9
      tmpape(jj)=zero
    end do

!   . aperture profile
    if( lconst.or.lAupDcum ) then
!   constant aperture or upstream aperture marker at the same s-location
!   -> it is wise to use the upstream aperture
      ixApeNewFrom=ixApeUp
    else if( lAdwDcum ) then
!   same s-location as the closest downstream aperture marker
!   -> it is wise to use it!
      ixApeNewFrom=ixApeDw
    else
!   varying aperture -> we need to interpolate
      call interp_aperture( iApeUp, ixApeUp, iApeDw, ixApeDw, itmpape, tmpape, dcum(iNew) )
      lfit=.true.
    end if

!   . aperture entry
    if( .not.lAss ) then
!     ixNew cannot be assigned an aperture marker: we have to insert
!     a new entry in the lattice sequence
      if( lfit ) then
        ixNew=INEESE()
        bez(ixNew)=CrtApeName()
      end if
      iNew=iNew+1
      iu=INEELS( iNew )
    end if

!   . assign aperture profile
    if( lAss.or.lfit ) then
!     aperture model must be copied
      call copy_aperture( ixNew,ixApeNewFrom,itmpape,tmpape )
      ic(iNew)=ixNew+nblo
    else if( ixApeNewFrom.gt.-1 ) then
!     an existing aperture model can be assigned
      ic(iNew)=ixApeNewFrom+nblo
    else
!     this should never happen
      write(lout,"(a)") "APER> ERROR in aperture auto assignment."
      call prror(-1)
    end if
  end if

! echo for checking
  write(lout,"(a)") "APER> Echo results of assignment:"
  call dump_aperture_header( lout )
  call dump_aperture_marker( lout, ic(iNew)-nblo, iNew )

! go home, man
  iEl=iNew
  return

 1982 format (a16,2(1x,a2),8(1x,f15.5))
end subroutine contour_aperture_marker

subroutine find_closest_aperture( iStart, lUp, iEl, ixEl, lfound )
!-----------------------------------------------------------------------
!     by A.Mereghetti (CERN, BE/ABP-HSS), 2018-03-24
!     find aperture marker closest to iStart
!-----------------------------------------------------------------------
#ifdef FLUKA
! import mod_fluka
! inserted in main code by the 'fluka' compilation flag
  use mod_fluka
#endif

  implicit none
  ! interface variables
  integer iStart, iEl, ixEl
  logical lUp, lfound
  ! temporary variables
  integer i, ix, iEnd, iStep

  iEl=-1
  ixEl=-1
  lfound=.false.

  ! search
  if(lUp) then
     iEnd=1
     iStep=-1
  else
     iEnd=iu
     iStep=1
  end if

  do i=iStart,iEnd,iStep
    ix=ic(i)-nblo
    if(ix.gt.0) then
!   SINGLE ELEMENT
#ifdef FLUKA
!     inserted in main code by the 'fluka' compilation flag
!     aperture markers should not coincide with a FLUKA element
      if( kape(ix).ne.0.and.fluka_type(ix).eq.FLUKA_NONE ) then
#else
      if( kape(ix).ne.0 ) then
#endif
        iEl=i
        ixEl=ix
        lfound=.true.
        exit
      end if
    end if
 end do
 return

end subroutine find_closest_aperture

function CrtApeName() result(retValue)
!-----------------------------------------------------------------------
!     by A.Mereghetti (CERN, BE/ABP-HSS)
!     last modified: 01-12-2016
!     Create Aperture Name
!     always in main code
!-----------------------------------------------------------------------
  implicit none

  character(len=mNameLen) retValue
  integer iApe, ii
  data iApe / 0 /
  save iApe

  iApe=iApe+1
  write(retValue, "(A10,I6)") "auto.aper.", iApe

  do ii=11,16
    if( retValue(ii:ii) .eq. ' ' ) retValue(ii:ii)='0'
  end do

end function CrtApeName

logical function sameAperture( ixApeUp, ixApeDw )
!-----------------------------------------------------------------------
!     by A.Mereghetti (CERN, BE/ABP-HSS)
!     last modified: 21-03-2018
!     Verify that two aperture markers actually describe the same aperture
!       restriction
!-----------------------------------------------------------------------
  implicit none
  integer ixApeUp, ixApeDw, jj
  sameAperture=ixApeDw.eq.ixApeUp.or.kape(ixApeDw).eq.kape(ixApeUp)
  if(sameAperture) then
     do jj=1,9
        sameAperture=sameAperture.and.abs(ape(jj,ixApeDw)-ape(jj,ixApeUp)).lt.aPrec
        if(.not.sameAperture) exit
     end do
  end if
end function sameAperture

subroutine interp_aperture( iUp,ixUp, iDw,ixDw, oKApe,oApe, spos )
!-----------------------------------------------------------------------
!     by A.Mereghetti
!     last modified: 21-03-2018
!     interpolate aperture
!     always in main code
!-----------------------------------------------------------------------
  implicit none

! interface variables
  integer iUp, ixUp, iDw, ixDw, oKApe
  real(kind=fPrec) oApe(9), spos
! temporary variables
  real(kind=fPrec) ddcum, mdcum
  integer jj

  if( sameAperture(ixUp,ixDw ) ) then
     ! constant aperture - no need to interpolate
     oKApe=kape(ixUp)
     do jj=1,9
        oApe(jj)=ape(jj,ixUp)
     end do
  else
     ! non-constant aperture - interpolate
     ! type: we may interpolate the same aperture type
     oKApe=-1 ! transition
     if( kape(ixUp).eq.kape(ixDw) ) oKApe=kape(ixUp)

     ! actual interpolation
     ddcum = spos-dcum(iUp)
     if( ddcum.lt.zero ) ddcum=dcum(iu)+ddcum
     mdcum = dcum(iDw)-dcum(iUp)
     if( mdcum.lt.zero ) mdcum=dcum(iu)+mdcum
     do jj=1,9
        oApe(jj)=((ape(jj,ixDw)-ape(jj,ixUp))/mdcum)*ddcum+ape(jj,ixUp)
     end do
  end if
  return
end subroutine interp_aperture

subroutine copy_aperture( ixApeTo, ixApeFrom, nKApe, nApe )
!-----------------------------------------------------------------------
!     by A.Mereghetti
!     last modified: 02-12-2016
!     copy aperture, either from an existing one or from the one
!       received on the fly
!     always in main code
!-----------------------------------------------------------------------
  implicit none

! interface variables
  integer ixApeTo, ixApeFrom, nKApe
  real(kind=fPrec) nApe(9)
! temporary variables
  integer jj

  if( ixApeFrom.gt.0 ) then
! copy aperture marker from existing SINGLE ELEMENT
    kape(ixApeTo)=kape(ixApeFrom)
    do jj=1,9
      ape(jj,ixApeTo)=ape(jj,ixApeFrom)
    end do
  else
! copy aperture marker from temporary one
    kape(ixApeTo)=nKApe
    do jj=1,9
      ape(jj,ixApeTo)=nApe(jj)
    end do
  end if

end subroutine copy_aperture

subroutine dump_aperture_model
!-----------------------------------------------------------------------
!     by P.Garcia Ortega, for the FLUKA Team, and A.Mereghetti
!     last modified: 08-12-2016
!     dump all apertures declared in machine
!     always in main code
!-----------------------------------------------------------------------
  use parpro
  use mod_units, only: f_open
  implicit none

! temporary variables
  integer i, ix
  logical lopen,err

  integer iOld, ixOld, niter, oKApe, jj
  real(kind=fPrec) aprr(9),slos
  character(len=mNameLen), parameter :: interpolated = 'interpolated'

  write(lout,"(a)") str_divLine
  write(lout,"(a)") ""
  write(lout,"(a)") " DUMP OF APERTURE MODEL"
  write(lout,"(a)") ""

  inquire( unit=aperunit, opened=lopen )
  if( .not.lopen ) then
    if( aperunit.ne.0 ) then
      call f_open(unit=aperunit,file=aper_filename,formatted=.true.,mode='w',err=err)
      write(lout,"(a)") "APER> Profile dumped in file: '"//trim(aper_filename)//"'"
    end if
  end if

! Header
  call dump_aperture_header( aperunit )

! First element of lattice
  i=1
  ix=ic(i)-nblo
  if( kape(ix).eq.0 ) then
    write(lout,"(a)") "APER> ERROR First element of lattice structure is not assigned any aperture type"
    call prror(-1)
  end if
  call dump_aperture_marker( aperunit, ix, i )
  iOld=i
  ixOld=ix

  do i=2,iu
    ix=ic(i)-nblo
    if(ix.gt.0) then
      ! SINGLE ELEMENT
      if( kape(ix) .ne. 0 ) then
        if(lbacktracking) then
          ! Number of iterations
          if( (dcum(i)-dcum(iOld)).gt.zero) then
            niter = nint((dcum(i)-dcum(iOld))/bktpre+one)
            do jj=1,niter
              slos = int(dcum(iOld)/bktpre+jj)*bktpre
              if( slos.lt.dcum(iOld) .or. slos.gt.dcum(i) ) exit
              call interp_aperture(iOld,ixOld,i,ix,oKApe,aprr,slos)
              call dump_aperture( aperunit, interpolated, oKApe, slos, aprr )
            end do
          end if
          iOld=i
          ixOld=ix
        end if
        call dump_aperture_marker( aperunit, ix, i )
      end if
    end if
  end do

  return

end subroutine dump_aperture_model

#ifdef HDF5
! ================================================================================================ !
!  DUMP APERTURE MODEL - HDF5 Version
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-10-02
! ================================================================================================ !
subroutine dump_aperture_model_hdf5

  use parpro
  use hdf5_output
  use string_tools

  implicit none

  integer i, ix
  integer iOld, ixOld, niter, oKApe, jj
  real(kind=fPrec) aprr(9),slos

  type(h5_dataField), allocatable :: setFields(:)
  character(len=:),   allocatable :: colNames(:)
  character(len=:),   allocatable :: colUnits(:)
  integer :: modelFmt, modelSet, nSplit, nTmp
  logical :: spErr

  allocate(setFields(12))

  setFields(1)  = h5_dataField(name="NAME",  type=h5_typeChar, size=mNameLen)
  setFields(2)  = h5_dataField(name="TYPE",  type=h5_typeChar, size=3)
  setFields(3)  = h5_dataField(name="S",     type=h5_typeReal)
  setFields(4)  = h5_dataField(name="APER1", type=h5_typeReal)
  setFields(5)  = h5_dataField(name="APER2", type=h5_typeReal)
  setFields(6)  = h5_dataField(name="APER3", type=h5_typeReal)
  setFields(7)  = h5_dataField(name="APER4", type=h5_typeReal)
  setFields(8)  = h5_dataField(name="APER5", type=h5_typeReal)
  setFields(9)  = h5_dataField(name="APER6", type=h5_typeReal)
  setFields(10) = h5_dataField(name="ANGLE", type=h5_typeReal)
  setFields(11) = h5_dataField(name="XOFF",  type=h5_typeReal)
  setFields(12) = h5_dataField(name="YOFF",  type=h5_typeReal)

  call h5_createFormat("aperModelFmt", setFields, modelFmt)
  call h5_createDataSet("model", h5_aperID, modelFmt, modelSet)
  call chr_split("name type s aper1 aper2 aper3 aper4 aper5 aper6 angle xoff yoff",colNames,nSplit,spErr)
  call chr_split("text text m mm mm mm*rad mm*rad mm*rad mm*rad rad mm mm",colUnits,nSplit,spErr)
  call h5_writeDataSetAttr(modelSet,"colNames",colNames)
  call h5_writeDataSetAttr(modelSet,"colUnits",colUnits)

  deallocate(setFields)

  ! First element of lattice
  i  = 1
  ix = ic(i)-nblo
  if(kape(ix) == 0) then
    write(lout,"(a)") "APER> ERROR First element of lattice structure is not assigned any aperture type"
    call prror(-1)
  end if
  call dump_aperture_hdf5(bez(ix), kape(ix), dcum(i), ape(1:9,ix), modelSet, .false.)
  iOld  = i
  ixOld = ix

  do i=2,iu
    ix = ic(i)-nblo
    if(ix > 0) then
      ! SINGLE ELEMENT
      if(kape(ix) /= 0) then
        if(lbacktracking) then
          ! Number of iterations
          if((dcum(i)-dcum(iOld)) > zero) then
            niter = nint((dcum(i)-dcum(iOld))/bktpre+1)
            do jj=1,niter
              slos = int(dcum(iOld)/bktpre+jj)*bktpre
              if(slos < dcum(iOld) .or. slos > dcum(i)) exit
              call interp_aperture(iOld,ixOld,i,ix,oKApe,aprr,slos)
              call dump_aperture_hdf5("Interpolated", oKApe, slos, aprr, modelSet, .false.)
            end do
          end if
          iOld  = i
          ixOld = ix
        end if
        call dump_aperture_hdf5(bez(ix), kape(ix), dcum(i), ape(1:9,ix), modelSet, .false.)
      end if
    end if
  end do

  ! Force a write of whatever is left in the buffer
  call dump_aperture_hdf5(" ", 0, zero, ape(1:9,1), modelSet, .true.)

  write(lout,"(a,i0)") "APER> Aperture model dumped to HDF5 dataset ",modelSet

end subroutine dump_aperture_model_hdf5

subroutine dump_aperture_hdf5(apName, apType, apSPos, apArr, dataSet, isEnd)

  character(len=*), intent(in) :: apName
  integer,          intent(in) :: apType
  real(kind=fPrec), intent(in) :: apSPos
  real(kind=fPrec), intent(in) :: apArr(9)
  integer,          intent(in) :: dataSet
  logical,          intent(in) :: isEnd

  ! Cache
  integer                 :: nRec = 0
  character(len=mNameLen) :: tmpName(1000)
  character(len=3)        :: tmpType(1000)
  real(kind=fPrec)        :: tmpSPos(1000)
  real(kind=fPrec)        :: tmpVals(1000,9)

  save nRec, tmpName, tmpType, tmpSPos, tmpVals

  if(nRec >= 1000 .or. isEnd) then
    call h5_prepareWrite(dataSet, nRec)
    call h5_writeData(dataSet, 1,  nRec, tmpName(1:nRec))
    call h5_writeData(dataSet, 2,  nRec, tmpType(1:nRec))
    call h5_writeData(dataSet, 3,  nRec, tmpSPos(1:nRec))
    call h5_writeData(dataSet, 4,  nRec, tmpVals(1:nRec,1))
    call h5_writeData(dataSet, 5,  nRec, tmpVals(1:nRec,2))
    call h5_writeData(dataSet, 6,  nRec, tmpVals(1:nRec,3))
    call h5_writeData(dataSet, 7,  nRec, tmpVals(1:nRec,4))
    call h5_writeData(dataSet, 8,  nRec, tmpVals(1:nRec,5))
    call h5_writeData(dataSet, 9,  nRec, tmpVals(1:nRec,6))
    call h5_writeData(dataSet, 10, nRec, tmpVals(1:nRec,7))
    call h5_writeData(dataSet, 11, nRec, tmpVals(1:nRec,8))
    call h5_writeData(dataSet, 12, nRec, tmpVals(1:nRec,9))
    call h5_finaliseWrite(dataSet)
    nRec = 0
    if(isEnd) return
  end if
  nRec = nRec + 1

  tmpName(nRec)      = apName
  tmpType(nRec)      = apeName(apType)
  tmpSPos(nRec)      = apSPos
  tmpVals(nRec, 1:9) = apArr(1:9)
  select case(apType)
  case(1) ! Circle
    tmpVals(nRec, 2:6) = zero
  case(2) ! Rectangle
    tmpVals(nRec, 3:6) = zero
  case(3) ! Ellipse
    tmpVals(nRec, 1:2) = apArr(3:4)
    tmpVals(nRec, 3:6) = zero
  case(4) ! Rectellipse
    tmpVals(nRec, 5:6) = zero
  case(5) ! Octagon
    tmpVals(nRec, 3)   = atan2_mb(apArr(1)*apArr(5) + apArr(6), apArr(1))
    tmpVals(nRec, 4)   = atan2_mb(apArr(2), (apArr(2) - apArr(6))/apArr(5))
    tmpVals(nRec, 5:6) = zero
  case(6) ! Racetrack
    tmpVals(nRec, 4:6) = zero
  end select

end subroutine dump_aperture_hdf5
#endif

subroutine dumpMe
  implicit none

! temporary variables
  integer i, ix

  write(lout,"(a)") "APER> dumpMe -----------------------------------------------------------------------------"
  do i=1,iu
    ix=ic(i)-nblo
    if( ix.gt.0 ) then
      write(lout,"(a,i8,1x,a48,1x,f15.6,1x,i8)") "APER> ",i,bez(ix),dcum(i),kape(ix)
    else
      write(lout,"(a,i8,1x,a48,1x,f15.6)") "APER> ",i,bezb(ic(i)),dcum(i)
    end if
  end do
  write(lout,"(a)") "APER> dumpMe -----------------------------------------------------------------------------"

end subroutine dumpMe

subroutine dump_aperture( iunit, name, aptype, spos, ape )
!-----------------------------------------------------------------------
!     by A.Mereghetti
!     last modified: 08-12-2016
!     dump any aperture marker
!     always in main code
!-----------------------------------------------------------------------
  use mod_settings
  implicit none

! interface variables
  integer iunit
  integer aptype
  character(len=mNameLen) name
  real(kind=fPrec) ape(9)
  real(kind=fPrec) spos

  ! Don't print to stdout if quiet flag is enabled.
  if(st_quiet > 0 .and. iunit == 6) return

  ! dump info
  if(ldmpaperMem) then
     write(iunit,1984) name, apeName(aptype), spos, ape(1), ape(2), ape(3), ape(4), ape(5), ape(6), ape(7), ape(8), ape(9)
  else
     select case(aptype)
     case(-1) ! transition
        write(iunit,1984) name, apeName(aptype), spos, ape(1), ape(2), ape(3), ape(4), ape(5), ape(6), ape(7), ape(8), ape(9)
     case(0) ! not an aperture marker
        write(iunit,1984) name, apeName(aptype), spos, ape(1), ape(2), ape(3), ape(4), ape(5), ape(6), ape(7), ape(8), ape(9)
     case(1) ! Circle
        write(iunit,1984) name, apeName(aptype), spos, ape(1),   zero,   zero,   zero,   zero,   zero, ape(7), ape(8), ape(9)
     case(2) ! Rectangle
        write(iunit,1984) name, apeName(aptype), spos, ape(1), ape(2),   zero,   zero,   zero,   zero, ape(7), ape(8), ape(9)
     case(3) ! Ellipse
        write(iunit,1984) name, apeName(aptype), spos, ape(3), ape(4),   zero,   zero,   zero,   zero, ape(7), ape(8), ape(9)
     case(4) ! Rectellipse
        write(iunit,1984) name, apeName(aptype), spos, ape(1), ape(2), ape(3), ape(4),   zero,   zero, ape(7), ape(8), ape(9)
     case(5) ! Octagon
        ! get angles from points passing through x1,y1 and x2,y2
        ! x1=ape(1)
        ! y1=ape(1)*tan(theta1)
        ! x2=ape(2)/tan(theta2)
        ! y2=ape(2)
        write(iunit,1984) name, apeName(aptype), spos, ape(1), ape(2), atan2_mb(ape(1)*ape(5)+ape(6),ape(1)), &
             &         atan2_mb(ape(2),(ape(2)-ape(6))/ape(5)),   zero,   zero, ape(7), ape(8), ape(9)
     case(6) ! Racetrack
        write(iunit,1984) name, apeName(aptype), spos, ape(1), ape(2), ape(3),   zero,   zero,   zero, ape(7), ape(8), ape(9)
     end select
  end if
  return
 1984 format (1x,a16,1x,a6,10(1x,f15.5))
end subroutine dump_aperture

subroutine dump_aperture_marker( iunit, ixEl, iEl )
!-----------------------------------------------------------------------
!     by A.Mereghetti
!     last modified: 08-12-2016
!     dump single aperture marker, existing in aperture DB
!     always in main code
!-----------------------------------------------------------------------
  implicit none

! interface variables
  integer iunit, iEl, ixEl

  call dump_aperture( iunit, bez(ixEl), kape(ixEl), dcum(iEl), ape(1:9,ixEl) )
  return
end subroutine dump_aperture_marker

subroutine dump_aperture_header( iunit )
!-----------------------------------------------------------------------
!     by A.Mereghetti
!     last modified: 22-03-2018
!     dump header of aperture marker
!     always in main code
!-----------------------------------------------------------------------
  use mod_settings
  implicit none
  integer iunit
  ! Don't print to stdout if quiet flag is enabled.
  if(st_quiet > 0 .and. iunit == 6) return
  write(iunit,1984) '#', 'name', 'aptype', 's[m]', 'aper1[mm]', 'aper2[mm]', &
 &                  'aper3[mm][rad]', 'aper4[mm][rad]', 'aper5[mm][rad]', 'aper6[mm][rad]', &
 &                  'angle[rad]', 'xoff[mm]', 'yoff[mm]'
  return
 1984 format (a1,a48,1x,a6,1x,10(1x,a15))
end subroutine dump_aperture_header

subroutine dump_aperture_xsecs
  !-----------------------------------------------------------------------
  ! A.Mereghetti (CERN, BE/ABP-HSS), 22-03-2018
  ! dump cross-sections of apertures at specific locations (loop)
  !-----------------------------------------------------------------------
  use mod_units, only: f_open, f_close
  implicit none
  ! temporary variables
  logical lfound, lopen, lApeUp, lApeDw, err
  integer ixsec, ierro, iEl, ixEl, iApeUp, ixApeUp, iApeDw, ixApeDw, itmpape
  real(kind=fPrec) sLoc, tmpape(9)

  ! loop over requested lines
  do ixsec=1,mxsec
     ! from print_lastlines_to_stderr
     inquire(unit=xsecunit(ixsec),opened=lopen)
     if(lopen) then
        write(lout,"(a,i0)")"APER> ERROR Dump_aperture_xsecs. Could not open file unit '"//trim(xsec_filename(ixsec))//&
          "' with unit ",xsecunit(ixsec)
        call prror(-1)
     end if
     call f_open(unit=xsecunit(ixsec),file=xsec_filename(ixsec),formatted=.true.,mode='w',err=err)
     if(ierro .ne. 0) then
        write(lout,"(2(a,i0))") "APER> ERROR Opening file '"//trim(xsec_filename(ixsec))//&
          "' on unit # ",xsecunit(ixsec),", iostat = ",ierro
        call prror(-1)
     end if

     ! loop over s-locations
     sLoc=sLocMin(ixsec)
     do while(sLoc.le.sLocMax(ixsec))
        call find_entry_at_s( sLoc, .true., iEl, ixEl, lfound )
        if(.not.lfound) call prror(-1)
        ! get upstream aperture marker
        call find_closest_aperture(iEl,.true.,iApeUp,ixApeUp,lApeUp)
        if( iApeUp.eq.-1 .and. ixApeUp.eq.-1 ) then
           write(lout,"(a)") "APER> ERROR Could not find upstream aperture marker"
           call prror(-1)
        end if
        ! get downstream aperture marker
        call find_closest_aperture(iEl,.false.,iApeDw,ixApeDw,lApeDw)
        if( iApeDw.eq.-1 .and. ixApeDw.eq.-1 ) then
           write(lout,"(a)") "APER> ERROR Could not find downstream aperture marker"
           call prror(-1)
        end if
        ! interpolate and get aperture at desired location
        call interp_aperture( iApeUp, ixApeUp, iApeDw, ixApeDw, itmpape, tmpape, sLoc )
        ! dump the x-sec of the aperture
        call dump_aperture_xsec(xsecunit(ixsec),itmpape,tmpape,nAzimuts(ixsec),sLoc)
        sLoc=sLoc+sLocDel(ixsec)
     end do

     call f_close(xsecunit(ixsec))
  end do

  return
end subroutine dump_aperture_xsecs


subroutine dump_aperture_xsec( iunit, itmpape, tmpape, nAzim, sLoc )
  !-----------------------------------------------------------------------
  ! A.Mereghetti (CERN, BE/ABP-HSS), 22-03-2018
  ! dump the cross-sections of the machine aperture at a specific location
  !-----------------------------------------------------------------------
  implicit none
  ! interface variables
  integer iunit, itmpape, nAzim
  real(kind=fPrec) tmpape(9), sLoc
  ! temporary variables
  logical tmpOffTlt
  integer i
  real(kind=fPrec) xChk, yChk, nChk, thetaRay, xRay, yRay

  write(iunit,*)'# aperture at s=',sLoc
  write(iunit,*)'# type:',itmpape
  write(iunit,*)'# specifiers:'
  do i=1,9
     write(iunit,*)'# - ape(',i,')=',tmpape(i)
  end do
  write(iunit,*)'# number of points:',nAzim
  write(iunit,1981) '# ang[deg]', 'rad [mm]', 'x [mm]', 'y [mm]'
  tmpOffTlt=tmpape(7).ne.zero.or.tmpape(8).ne.zero.or.tmpape(9).ne.zero

  ! origin of ray:
  xRay=zero
  yRay=zero
  if(tmpOffTlt) call roffpos(xRay,yRay,xRay,yRay,tmpape(7),tmpape(8),tmpape(9))

  ! loop over rays
  select case(itmpape)
  case(-1) ! transition
     do i=1,nAzim
        thetaRay=(i/real(nAzim))*(two*pi) ! radians
        ! call (angle to aperture ref sys)
        call intersectTR(xRay,yRay,thetaRay-tmpape(7),tmpape(1),tmpape(2),tmpape(3),tmpape(4),tmpape(5),tmpape(6),xChk,yChk,nChk)
        ! go back to machine reference system
        if(tmpOffTlt) call roffpos_inv(xChk,yChk,xChk,yChk,tmpape(7),tmpape(8),tmpape(9))
        ! echo result of scan
        write(iunit,1982) thetaRay/rad,sqrt(xChk**2+yChk**2),xChk,yChk
     end do
  case(1) ! circle
     do i=1,nAzim
        thetaRay=(i/real(nAzim))*(two*pi) ! radians
        ! call (angle to aperture ref sys)
        call intersectCR(xRay,yRay,thetaRay-tmpape(7),tmpape(3),zero,zero,xChk,yChk,nChk)
        ! go back to machine reference system
        if(tmpOffTlt) call roffpos_inv(xChk,yChk,xChk,yChk,tmpape(7),tmpape(8),tmpape(9))
        ! echo result of scan
        write(iunit,1982) thetaRay/rad,sqrt(xChk**2+yChk**2),xChk,yChk
     end do
  case(2) ! rectangle
     do i=1,nAzim
        thetaRay=(i/real(nAzim))*(two*pi) ! radians
        ! call (angle to aperture ref sys)
        call intersectRE(xRay,yRay,thetaRay-tmpape(7),tmpape(1),tmpape(2),xChk,yChk,nChk)
        ! go back to machine reference system
        if(tmpOffTlt) call roffpos_inv(xChk,yChk,xChk,yChk,tmpape(7),tmpape(8),tmpape(9))
        ! echo result of scan
        write(iunit,1982) thetaRay/rad,sqrt(xChk**2+yChk**2),xChk,yChk
     end do
  case(3) ! ellipse
     do i=1,nAzim
        thetaRay=(i/real(nAzim))*(two*pi) ! radians
        ! call (angle to aperture ref sys)
        call intersectEL(xRay,yRay,thetaRay-tmpape(7),tmpape(3),tmpape(4),zero,zero,xChk,yChk,nChk)
        ! go back to machine reference system
        if(tmpOffTlt) call roffpos_inv(xChk,yChk,xChk,yChk,tmpape(7),tmpape(8),tmpape(9))
        ! echo result of scan
        write(iunit,1982) thetaRay/rad,sqrt(xChk**2+yChk**2),xChk,yChk
     end do
  case(4) ! rectellipse
     do i=1,nAzim
        thetaRay=(i/real(nAzim))*(two*pi) ! radians
        ! call (angle to aperture ref sys)
        call intersectRL(xRay,yRay,thetaRay-tmpape(7),tmpape(1),tmpape(2),tmpape(3),tmpape(4),xChk,yChk,nChk)
        ! go back to machine reference system
        if(tmpOffTlt) call roffpos_inv(xChk,yChk,xChk,yChk,tmpape(7),tmpape(8),tmpape(9))
        ! echo result of scan
        write(iunit,1982) thetaRay/rad,sqrt(xChk**2+yChk**2),xChk,yChk
     end do
  case(5) ! octagon
     do i=1,nAzim
        thetaRay=(i/real(nAzim))*(two*pi) ! radians
        ! call (angle to aperture ref sys)
        call intersectOC(xRay,yRay,thetaRay-tmpape(7),tmpape(1),tmpape(2),tmpape(5),tmpape(6),xChk,yChk,nChk)
        ! go back to machine reference system
        if(tmpOffTlt) call roffpos_inv(xChk,yChk,xChk,yChk,tmpape(7),tmpape(8),tmpape(9))
        ! echo result of scan
        write(iunit,1982) thetaRay/rad,sqrt(xChk**2+yChk**2),xChk,yChk
     end do
  case(6) ! racetrack
     do i=1,nAzim
        thetaRay=(i/real(nAzim))*(two*pi) ! radians
        ! call (angle to aperture ref sys)
        call intersectRT(xRay,yRay,thetaRay-tmpape(7),tmpape(1),tmpape(2),tmpape(3),xChk,yChk,nChk)
        ! go back to machine reference system
        if(tmpOffTlt) call roffpos_inv(xChk,yChk,xChk,yChk,tmpape(7),tmpape(8),tmpape(9))
        ! echo result of scan
        write(iunit,1982) thetaRay/rad,sqrt(xChk**2+yChk**2),xChk,yChk
     end do
  end select
  write(iunit,*)''
  write(iunit,*)''
  return
1981 FORMAT(1(1X,a10))
1982 FORMAT(4(1x,f10.5))

end subroutine dump_aperture_xsec

subroutine intersectCR( xRay, yRay, thetaRay, radius, x0, y0, xChk, yChk, nChk )
  ! 0.0<=thetaRay<=2pi!!!!!
  implicit none
  ! interface variables
  real(kind=fPrec) xRay, yRay, thetaRay, radius, x0, y0, xChk, yChk, nChk
  ! temp variables
  real(kind=fPrec) mRay, qRay, delta, tmpX0, tmpY0

  ! scanning ray:
  xChk=c1e3
  yChk=c1e3
  nChk=c1e3
  mRay=zero
  qRay=zero
  if(abs(thetaRay).lt.c1m6.or.abs(thetaRay/(two*pi)-one).lt.c1m6) then ! thetaRay=0.0 or thetaRay=2pi
     yChk=zero
     xChk=radius
  else if(abs(thetaRay/(pi/two)-one).lt.c1m6) then ! thetaRay=0.5pi
     yChk=radius
     xChk=zero
  else if(abs(thetaRay/pi-one).lt.c1m6) then ! thetaRay=pi
     yChk=zero
     xChk=-radius
  else if(abs(thetaRay/(pi*(three/two))-one).lt.c1m6) then ! thetaRay=1.5pi
     yChk=-radius
     xChk=zero
  else
     mRay=tan_mb(thetaRay)
     qRay=yRay-mRay*xRay
     if(zero.lt.thetaRay.and.thetaRay.lt.pi/two) then ! first quadrant
        tmpX0=abs(x0)
        tmpY0=abs(y0)
     else if(pi/two.lt.thetaRay.and.thetaRay.lt.pi) then ! second quadrant
        tmpX0=-abs(x0)
        tmpY0=abs(y0)
     else if(pi.lt.thetaRay.and.thetaRay.lt.pi*(three/two)) then ! second quadrant
        tmpX0=-abs(x0)
        tmpY0=-abs(y0)
     else ! fourth quadrant
        tmpX0=abs(x0)
        tmpY0=-abs(y0)
     end if
     delta=-((mRay*tmpX0-tmpY0)+qRay)**2+radius**2*(one+mRay**2)
     if(delta.lt.zero) return
     if((zero.lt.thetaRay.and.thetaRay.lt.pi/two) .or. & ! first quadrant
 &       (pi*(three/two).lt.thetaRay.and.thetaRay.lt.two*pi)) then ! fourth quadrant
        xChk=((tmpX0+mRay*(tmpY0-qRay))+sqrt(delta))/(one+mRay**2)
     else
        xChk=((tmpX0+mRay*(tmpY0-qRay))-sqrt(delta))/(one+mRay**2)
     end if
     yChk=mRay*xChk+qRay
  end if
  nChk=radius
  return
end subroutine intersectCR

subroutine intersectRE( xRay, yRay, thetaRay, xRe, yRe, xChk, yChk, nChk )
  ! 0.0<=thetaRay<=2pi!!!!!
  implicit none
  ! interface variables
  real(kind=fPrec) xRay, yRay, thetaRay, xRe, yRe, xChk, yChk, nChk
  ! temp variables
  real(kind=fPrec) mRay, qRay, xTmp(2), yTmp(2), nTmp(2)

  ! scanning ray:
  xChk=zero
  yChk=zero
  nChk=zero
  mRay=zero
  qRay=zero
  if(abs(thetaRay).lt.c1m6.or.abs(thetaRay/(two*pi)-one).lt.c1m6) then ! thetaRay=0.0 or thetaRay=2pi
     yChk=zero
     xChk=xRe
     nChk=xRe
  else if(abs(thetaRay/(pi/two)-one).lt.c1m6) then ! thetaRay=0.5pi
     yChk=yRe
     xChk=zero
     nChk=yRe
  else if(abs(thetaRay/pi-one).lt.c1m6) then ! thetaRay=pi
     yChk=zero
     xChk=-xRe
     nChk=xRe
  else if(abs(thetaRay/(pi*(three/two))-one).lt.c1m6) then ! thetaRay=1.5pi
     yChk=-yRe
     xChk=zero
     nChk=yRe
  else
     mRay=tan_mb(thetaRay)
     qRay=yRay-mRay*xRay
     if(zero.lt.thetaRay.and.thetaRay.lt.pi/two) then ! first quadrant
        xTmp(1)=xRe
        yTmp(2)=yRe
     else if(pi/two.lt.thetaRay.and.thetaRay.lt.pi) then ! second quadrant
        xTmp(1)=-xRe
        yTmp(2)=yRe
     else if(pi.lt.thetaRay.and.thetaRay.lt.pi*(three/two)) then ! third quadrant
        xTmp(1)=-xRe
        yTmp(2)=-yRe
     else ! fourth quadrant
        xTmp(1)=xRe
        yTmp(2)=-yRe
     end if
     yTmp(1)=mRay*xTmp(1)+qRay
     xTmp(2)=(yTmp(2)-qRay)/mRay
     nTmp(1)=xTmp(1)**2+yTmp(1)**2
     nTmp(2)=xTmp(2)**2+yTmp(2)**2
     if(nTmp(1).lt.nTmp(2)) then
        xChk=xTmp(1)
        yChk=yTmp(1)
     else
        xChk=xTmp(2)
        yChk=yTmp(2)
     end if
     nChk=sqrt(xChk**2+yChk**2)
  end if
  return
end subroutine intersectRE

subroutine intersectEL( xRay, yRay, thetaRay, aa, bb, x0, y0, xChk, yChk, nChk )
  ! 0.0<=thetaRay<=2pi!!!!!
  implicit none
  ! interface variables
  real(kind=fPrec) xRay, yRay, thetaRay, aa, bb, x0, y0, xChk, yChk, nChk
  ! temp variables
  real(kind=fPrec) mRay, qRay, delta, tmpX0, tmpY0

  ! scanning ray:
  xChk=c1e3
  yChk=c1e3
  nChk=c1e3
  mRay=zero
  qRay=zero
  if(abs(thetaRay).lt.c1m6.or.abs(thetaRay/(two*pi)-one).lt.c1m6) then ! thetaRay=0.0 or thetaRay=2pi
     yChk=zero
     xChk=aa
     nChk=aa
  else if(abs(thetaRay/(pi/two)-one).lt.c1m6) then ! thetaRay=0.5pi
     yChk=bb
     xChk=zero
     nChk=bb
  else if(abs(thetaRay/pi-one).lt.c1m6) then ! thetaRay=pi
     yChk=zero
     xChk=-aa
     nChk=aa
  else if(abs(thetaRay/(pi*(three/two))-one).lt.c1m6) then ! thetaRay=1.5pi
     yChk=-bb
     xChk=zero
     nChk=bb
  else
     mRay=tan_mb(thetaRay)
     qRay=yRay-mRay*xRay
     if(zero.lt.thetaRay.and.thetaRay.lt.pi/two) then ! first quadrant
        tmpX0=abs(x0)
        tmpY0=abs(y0)
     else if(pi/two.lt.thetaRay.and.thetaRay.lt.pi) then ! second quadrant
        tmpX0=-abs(x0)
        tmpY0=abs(y0)
     else if(pi.lt.thetaRay.and.thetaRay.lt.pi*(three/two)) then ! second quadrant
        tmpX0=-abs(x0)
        tmpY0=-abs(y0)
     else ! fourth quadrant
        tmpX0=abs(x0)
        tmpY0=-abs(y0)
     end if
     delta=-((mRay*tmpX0-tmpY0)+qRay)**2+(bb**2+aa**2*mRay**2)
     if(delta.lt.zero) return
     if((zero.lt.thetaRay.and.thetaRay.lt.pi/two).or. & ! first quadrant
 &       (pi*(three/two).lt.thetaRay.and.thetaRay.lt.two*pi)) then ! fourth quadrant
        xChk=((aa**2*(mRay*(tmpY0-qRay))+bb**2*tmpX0)+(aa*bb)*sqrt(delta))/(bb**2+aa**2*mRay**2)
     else
        xChk=((aa**2*(mRay*(tmpY0-qRay))+bb**2*tmpX0)+(-(aa*bb))*sqrt(delta))/(bb**2+aa**2*mRay**2)
     end if
     yChk=mRay*xChk+qRay
     nChk=sqrt(xChk**2+yChk**2)
  end if
  return
end subroutine intersectEL

subroutine intersectRL( xRay, yRay, thetaRay, xRe, yRe, aa, bb, xChk, yChk, nChk )
  ! 0.0<=thetaRay<=2pi!!!!!
  implicit none
  ! interface variables
  real(kind=fPrec) xRay, yRay, thetaRay, xRe, yRe, aa, bb, xChk, yChk, nChk
  ! temp variables
  real(kind=fPrec) xTmp(2), yTmp(2), nTmp(2)
  call intersectRE( xRay, yRay, thetaRay, xRe, yRe, xTmp(1), yTmp(1), nTmp(1) )
  call intersectEL( xRay, yRay, thetaRay, aa, bb, zero, zero, xTmp(2), yTmp(2), nTmp(2) )
  if(nTmp(1).lt.nTmp(2)) then
     xChk=xTmp(1)
     yChk=yTmp(1)
     nChk=nTmp(1)
  else
     xChk=xTmp(2)
     yChk=yTmp(2)
     nChk=nTmp(2)
  end if
  return
end subroutine intersectRL

subroutine intersectLN( xRay, yRay, thetaRay, mLine, qLine, xChk, yChk, nChk )
  ! 0.0<=thetaRay<=2pi!!!!!
  implicit none
  ! interface variables
  real(kind=fPrec) xRay, yRay, thetaRay, mLine, qLine, xChk, yChk, nChk
  ! temp variables
  real(kind=fPrec) mRay, qRay, mTmp, qTmp

  ! scanning ray:
  xChk=zero
  yChk=zero
  nChk=zero
  mRay=zero
  qRay=zero
  mTmp=zero
  qTmp=zero
  if(abs(thetaRay).lt.c1m6.or.abs(thetaRay/(two*pi)-one).lt.c1m6) then ! thetaRay=0.0 or thetaRay=2pi
     yChk=zero
     xChk=-qLine/mLine
     nChk=abs(qLine/mLine)
  else if(abs(thetaRay/(pi/two)-one).lt.c1m6) then ! thetaRay=0.5pi
     yChk=qLine
     xChk=zero
     nChk=abs(qLine)
  else if(abs(thetaRay/pi-one).lt.c1m6) then ! thetaRay=pi
     yChk=zero
     xChk=qLine/mLine
     nChk=abs(qLine/mLine)
  else if(abs(thetaRay/(pi*(three/two))-one).lt.c1m6) then ! thetaRay=1.5pi
     yChk=-qLine
     xChk=zero
     nChk=abs(qLine)
  else
     mRay=tan_mb(thetaRay)
     qRay=yRay-mRay*xRay
     if(zero.lt.thetaRay.and.thetaRay.lt.pi/two) then ! first quadrant
        mTmp=mLine
        qTmp=qLine
     else if(pi/two.lt.thetaRay.and.thetaRay.lt.pi) then ! second quadrant
        mTmp=-mLine
        qTmp=qLine
     else if(pi.lt.thetaRay.and.thetaRay.lt.pi*(three/two)) then ! third quadrant
        mTmp=mLine
        qTmp=-qLine
     else ! fourth quadrant
        mTmp=-mLine
        qTmp=-qLine
     end if
     xChk=-(qRay-qLine)/(mRay-mLine)
     yChk=mRay*xChk+qRay
     nChk=sqrt(xChk**2+yChk**2)
  end if
  return
end subroutine intersectLN

subroutine intersectOC( xRay, yRay, thetaRay, xRe, yRe, mOct, qOct, xChk, yChk, nChk )
  ! 0.0<=thetaRay<=2pi!!!!!
  use numerical_constants
  implicit none
  ! interface variables
  real(kind=fPrec) xRay, yRay, thetaRay, xRe, yRe, mOct, qOct, xChk, yChk, nChk
  ! temp variables
  real(kind=fPrec) xTmp(2), yTmp(2), nTmp(2)
  call intersectRE( xRay, yRay, thetaRay,  xRe,  yRe, xTmp(1), yTmp(1), nTmp(1) )
  call intersectLN( xRay, yRay, thetaRay, mOct, qOct, xTmp(2), yTmp(2), nTmp(2) )
  ! octagon part
  if(nTmp(1).lt.nTmp(2)) then
     xChk=xTmp(1)
     yChk=yTmp(1)
     nChk=nTmp(1)
  else
     xChk=xTmp(2)
     yChk=yTmp(2)
     nChk=nTmp(2)
  end if
  return
end subroutine intersectOC

subroutine intersectRT( xRay, yRay, thetaRay, xRe, yRe, radius, xChk, yChk, nChk )
  ! 0.0<=thetaRay<=2pi!!!!!
  implicit none
  ! interface variables
  real(kind=fPrec) xRay, yRay, thetaRay, xRe, yRe, radius, xChk, yChk, nChk
  ! temp variables
  real(kind=fPrec) xTmp(2), yTmp(2), nTmp(2)
  call intersectRE( xRay, yRay, thetaRay, xRe, yRe, xTmp(1), yTmp(1), nTmp(1) )
  call intersectCR( xRay, yRay, thetaRay, radius, xRe-radius, yRe-radius, xTmp(2), yTmp(2), nTmp(2) )
  if(nTmp(1).lt.nTmp(2)) then
     xChk=xTmp(1)
     yChk=yTmp(1)
     nChk=nTmp(1)
  else
     xChk=xTmp(2)
     yChk=yTmp(2)
     nChk=nTmp(2)
  end if
  return
end subroutine intersectRT

subroutine intersectTR( xRay, yRay, thetaRay, xRe, yRe, aa, bb, mOct, qOct, xChk, yChk, nChk )
  ! 0.0<=thetaRay<=2pi!!!!!
  implicit none
  ! interface variables
  real(kind=fPrec) xRay, yRay, thetaRay, xRe, yRe, aa, bb, mOct, qOct, xChk, yChk, nChk
  ! temp variables
  real(kind=fPrec) xTmp(2), yTmp(2), nTmp(2)
  call intersectRE( xRay, yRay, thetaRay, xRe, yRe, xTmp(1), yTmp(1), nTmp(1) )
  call intersectEL( xRay, yRay, thetaRay, aa, bb, xRe-aa, yRe-bb, xTmp(2), yTmp(2), nTmp(2) )
  if(nTmp(1).gt.nTmp(2)) then
     xTmp(1)=xTmp(2)
     yTmp(1)=yTmp(2)
     nTmp(1)=nTmp(2)
  end if
  call intersectLN( xRay, yRay, thetaRay, mOct, qOct, xTmp(2), yTmp(2), nTmp(2) )
  if(nTmp(1).lt.nTmp(2)) then
     xChk=xTmp(1)
     yChk=yTmp(1)
     nChk=nTmp(1)
  else
     xChk=xTmp(2)
     yChk=yTmp(2)
     nChk=nTmp(2)
  end if
  return
end subroutine intersectTR

! ================================================================================================ !
!  APERTURE LIMITATIONS PARSING
!  A. Mereghetti, P. Garcia Ortega and D. Sinuela Pastor, for the FLUKA Team
!  J. Molson, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-12-20
!  Input parsing split up, updated and moved from DATEN by VKBO.
!  Original LIMI block extended to deal with RectEllipse, Octagon and RaceTrack aperture types,
!    and with offset/tilting of profile.
!  Possibility to read the apertures from external file with LOAD keyword
! ================================================================================================ !
subroutine aper_parseLoadFile(load_file, iLine, iErr)

  use parpro, only : mInputLn
  use mod_units

  implicit none

  character(len=64),intent(in)    :: load_file
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=mInputLn) unitLine
  integer                 iErro, lineNo, loadunit
  logical                 err, lExist

  lineNo=0
  
  call f_requestUnit(trim(load_file),loadunit)
  inquire(file=load_file, exist=lExist)
  if(.not.lexist) then
    write(lout,"(a)") "LIMI> ERROR LOAD file '"//trim(load_file)//"' not found in the running folder."
    iErr = .true.
    return
  end if
  call f_open(unit=loadunit,file=load_file,formatted=.true.,mode='r',err=err)

  ! iterate through LOAD file
10 continue
  read(loadunit,"(a)",end=90,iostat=iErro) unitLine
  if(iErro > 0) then
    write(lout,"(a,i0)") "LIMI> ERROR Could not read from unit ",loadunit
    call prror(-1)
  end if
  lineNo = lineNo + 1

  if(len_trim(unitLine) == 0)  goto 10 ! Empty line, ignore
  if(unitLine(1:1) == "/")     goto 10 ! Comment line, ignore
  if(unitLine(1:1) == "!")     goto 10 ! Comment line, ignore
  if(unitLine(1:4) == "LIMI")  goto 10 ! header from MADX, ignore
  if(unitLine(1:4) == "NEXT")  goto 10 ! closure by MADX, ignore

  call aper_parseInputLine(unitLine, iLine, iErr)
  if(iErr) then
    write(lout,"(a)")      "LIMI> ERROR in external LIMI file."
    write(lout,"(a,i0,a)") "LIMI> Line ",lineNo,": '"//trim(unitLine)//"'"
    return
  end if
  goto 10

90 continue
  write(lout,"(a,i0,a)") "LIMI> Read ",lineNo," lines from external file."
  call f_close(loadunit)
  return

end subroutine aper_parseLoadFile

recursive subroutine aper_parseInputLine(inLine, iLine, iErr)

  use string_tools
  use sixtrack_input
  use mod_units

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=64)               :: load_file
  real(kind=fPrec) tmplen,tmpflts(3)
  integer          nSplit, i
  logical          spErr, lExist, apeFound, err

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "LIMI> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  select case(lnSplit(1))

  case("LOAD")
    ! P.G.Ortega and A.Mereghetti, 02-03-2018
    ! Reading apertures from external file
    if(nSplit .ne. 2 ) then
      write(lout,"(a,i0)") "LIMI> ERROR Wrong number of input parameters for keyword LOAD. Expected 2, got ",nSplit
      iErr = .true.
      return
    end if

    load_file = trim(lnSplit(2))
    write(lout,"(a)") "LIMI> Apertures will be read from file '"//trim(load_file)//"'"
    call aper_parseLoadFile(load_file, iLine, iErr)
    if(iErr) return

  case("PRIN","PRINT")
    ! P.G.Ortega and A.Mereghetti, 02-03-2018
    ! flag for dumping the aperture model
    if(nSplit < 2 .and. nSplit > 3 ) then
      write(lout,"(a,i0)") "LIMI> ERROR Wrong number of input parameters for keyword PRIN. Expected 2 or 3, got ",nSplit
      iErr = .true.
      return
    end if

    aper_filename = trim(lnSplit(2))
    call f_requestUnit(trim(aper_filename),aperunit)

    ldmpaper = .true.
    if(nSplit .eq. 3) then
      if(lnSPlit(3) == "MEM") then
        ldmpaperMem=.true.
      else
        write(lout,"(a,a)") "LIMI> ERROR Unknown third argument to PRIN keyword: ",lnSPlit(3)
        iErr = .true.
        return
      end if
    end if

    if(aperture_debug) then
      call sixin_echoVal("aperunit",aperunit,     "LIMI",iLine)
      call sixin_echoVal("filename",aper_filename,"LIMI",iLine)
    end if

  case("DEBUG")
    aperture_debug = .true.
    write(lout,"(a)") "LIMI> DEBUG is enabled."

  case("SAVE")
    ! P.G.Ortega and A.Mereghetti, 02-03-2018
    ! flag for saving particles at aperture check
    apflag = .true.
    write(lout,"(a)") "LIMI> SAVE is enabled."

  case("BACKTRKOFF")
    ! A.Mereghetti, 07-03-2018
    ! switch off back tracking
    lbacktracking=.false.
    write(lout,"(a)") "LIMI> Backtracking is disabled."

  case("PREC")
    ! A.Mereghetti and P.Garcia Ortega, 02-03-2018
    ! set precision for back-tracking
    if(nSplit < 2) then
      write(lout,"(a,i0)") "LIMI> ERROR Wrong number of input parameters for keyword PREC. Expected 2, got ",nSplit
      iErr = .true.
      return
    end if

    call chr_cast(lnSplit(2),tmplen,iErr)
    if(tmplen <= zero) then
      write(lout,"(a,e22.15)") "LIMI> WARNING Wrong precision value: ",tmplen
      write(lout,"(a,e22.15)") "LIMI>         Ignoring. Using default [m]: ",bktpre
    else
      bktpre = tmplen
    endif
    lbacktracking=.true.
    write(lout,"(a)") "LIMI> Backtracking is on."

  case("XSEC")
    write(lout,"(a)") "LIMI> ERROR Dump of aperture cross sections at specific locations are not available yet"
    iErr = .true.
    return

    ! A.Mereghetti, 22-03-2018
    ! ask for xsec at specific locations
    ! example input line:        XSEC myCrossSec.dat 12355.78 12356.78 0.1 180
    if(nSplit < 3) then
      write(lout,"(a,i0)") "LIMI> ERROR Wrong number of input parameters for keyword XSEC. Expected at least 3, got ",nSplit
      iErr = .true.
      return
    end if

    mxsec = mxsec + 1
    if(mxsec > nxsec) then
      write(lout,"(2(a,i0))") "LIMI> ERROR Too many xsecs! Asked for ",mxsec,", but max is ",nxsec
      iErr = .true.
      return
    end if

    xsec_filename(mxsec) = lnSplit(2)
    call chr_cast(lnSplit(3),sLocMin(mxsec),iErr)
    call f_requestUnit(xsec_filename(mxsec),xsecunit(mxsec))

    if(sLocMin(mxsec) < zero) then
      write(lout,"(a)") "LIMI> ERROR Negative min s-value for xsecs!"
      iErr = .true.
      return
    end if

    if(nSplit > 3) then
      call chr_cast(lnSplit(4),sLocMax(mxsec),iErr)
      if (sLocMax(mxsec).lt.zero) then
        write(lout,"(a)") "LIMI> ERROR Negative max s-value for xsecs!"
        iErr = .true.
        return
      end if
    end if

    if(nSplit > 4) then
      call chr_cast(lnSplit(5),sLocDel(mxsec),iErr)
      if(sLocDel(mxsec) < zero) sLocDel(mxsec)=-sLocDel(mxsec) ! increasing s-val
    end if

    if(sLocMax(mxsec) == zero) sLocMax(mxsec)=sLocMin(mxsec)
    if(sLocMax(mxsec) <  sLocMin(mxsec)) then
      block
        ! swap sMin and sMax
        real(kind=fPrec) tmpflts(1)
        tmpflts(1)     = sLocMax(mxsec)
        sLocMax(mxsec) = sLocMin(mxsec)
        sLocMin(mxsec) = tmpflts(1)
      end block
    end if

    if(nSplit > 5) call chr_cast(lnSplit(6),nAzimuts(mxsec),iErr)

    if(aperture_debug) then
      call sixin_echoVal("xsecunit",xsecunit(mxsec),"LIMI",iLine)
      call sixin_echoVal("filename",xsec_filename(mxsec),"LIMI",iLine)
      call sixin_echoVal("sLocMin", sLocMin(mxsec),"LIMI",iLine)
      if(nSplit > 3) call sixin_echoVal("sLocMax", sLocMax(mxsec), "LIMI",iLine)
      if(nSplit > 4) call sixin_echoVal("sLocDel", sLocDel(mxsec), "LIMI",iLine)
      if(nSplit > 5) call sixin_echoVal("nAzimuts",nAzimuts(mxsec),"LIMI",iLine)
    end if

  case default

    apeFound = .false.
    do i=1,il
      if(bez(i) == lnSplit(1)) then
        call aper_parseElement(inLine, i, iErr)
        apeFound = .true.
        if(kape(i) == -1) then
          if(nSplit > 8)  call chr_cast(lnSplit(9), tmpflts(1),iErr)
          if(nSplit > 9)  call chr_cast(lnSplit(10),tmpflts(2),iErr)
          if(nSplit > 10) call chr_cast(lnSplit(11),tmpflts(3),iErr)
        else
          if(nSplit > 6)  call chr_cast(lnSplit(7), tmpflts(1),iErr)
          if(nSplit > 7)  call chr_cast(lnSplit(8), tmpflts(2),iErr)
          if(nSplit > 8)  call chr_cast(lnSplit(9), tmpflts(3),iErr)
        end if
        call aperture_initroffpos(i,tmpflts(1),tmpflts(2),tmpflts(3))
        limifound = .true.
      end if
    end do
    if(.not.apeFound) then
      write(lout,"(a)") "LIMI> WARNING Unidentified element '"//lnSplit(1)//"', ignoring ..."
    end if

  end select

end subroutine aper_parseInputLine

subroutine aper_parseElement(inLine, iElem, iErr)

  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iElem
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  real(kind=fPrec) tmpflts(6)
  integer          nSplit, i
  logical          spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "LIMI> ERROR Failed to parse element."
    iErr = .true.
    return
  end if

  if(nSplit < 2) then
    write(lout,"(a)") "LIMI> ERROR Invalid entry."
    iErr = .true.
    return
  end if

  select case(lnSplit(2))

  case(apeName(1)) ! Circle
    if(nSplit < 3) then
      write(lout,"(a,i0)") "LIMI> ERROR Wrong number of input parameters for the '"//apeName(1)//&
        "' aperture marker. Expected 3, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(3),tmpflts(1),iErr)
    call aperture_initCR(iElem,tmpflts(1))

  case(apeName(2)) ! Rectangle
    if(nSplit < 4) then
      write(lout,"(a,i0)") "LIMI> ERROR Wrong number of input parameters for the '"//apeName(2)//&
        "' aperture marker. Expected 4, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(3),tmpflts(1),iErr)
    call chr_cast(lnSplit(4),tmpflts(2),iErr)
    call aperture_initRE(iElem,tmpflts(1),tmpflts(2))

  case(apeName(3)) ! Ellipse
    if(nSplit < 4) then
      write(lout,"(a,i0)") "LIMI> ERROR Wrong number of input parameters for the '"//apeName(3)//&
        "' aperture marker. Expected 4, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(3),tmpflts(1),iErr)
    call chr_cast(lnSplit(4),tmpflts(2),iErr)
    call aperture_initEL(iElem,tmpflts(1),tmpflts(2))

  case(apeName(4)) ! Rectellipse
    if(nSplit < 6) then
      write(lout,"(a,i0)") "LIMI> ERROR Wrong number of input parameters for the '"//apeName(4)//&
        "' aperture marker. Expected 6, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(3),tmpflts(1),iErr)
    call chr_cast(lnSplit(4),tmpflts(2),iErr)
    call chr_cast(lnSplit(5),tmpflts(3),iErr)
    call chr_cast(lnSplit(6),tmpflts(4),iErr)
    call aperture_initRL(iElem,tmpflts(1),tmpflts(2),tmpflts(3),tmpflts(4))

  case(apeName(5)) ! Octagon
    if(nSplit < 6) then
      write(lout,"(a,i0)") "LIMI> ERROR Wrong number of input parameters for the '"//apeName(5)//&
        "' aperture marker. Expected 6, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(3),tmpflts(1),iErr)
    call chr_cast(lnSplit(4),tmpflts(2),iErr)
    call chr_cast(lnSplit(5),tmpflts(3),iErr)
    call chr_cast(lnSplit(6),tmpflts(4),iErr)
    call aperture_initOC(iElem,tmpflts(1),tmpflts(2),tmpflts(3),tmpflts(4))

  case(apeName(6)) ! Racetrack
    if(nSplit < 5) then
      write(lout,"(a,i0)") "LIMI> ERROR Wrong number of input parameters for the '"//apeName(6)//&
        "' aperture marker. Expected 5, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(3),tmpflts(1),iErr)
    call chr_cast(lnSplit(4),tmpflts(2),iErr)
    call chr_cast(lnSplit(5),tmpflts(3),iErr)
    call aperture_initRT(iElem,tmpflts(1),tmpflts(2),tmpflts(3))

  case(apeName(-1)) ! Transition
    if(nSplit < 8) then
      write(lout,"(a,i0)") "LIMI> ERROR Wrong number of input parameters for the '"//apeName(-1)//&
        "' aperture marker. Expected 8, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(3),tmpflts(1),iErr)
    call chr_cast(lnSplit(4),tmpflts(2),iErr)
    call chr_cast(lnSplit(5),tmpflts(3),iErr)
    call chr_cast(lnSplit(6),tmpflts(4),iErr)
    call chr_cast(lnSplit(7),tmpflts(5),iErr)
    call chr_cast(lnSplit(8),tmpflts(6),iErr)
    call aperture_initTR(iElem,tmpflts(1),tmpflts(2),tmpflts(3),tmpflts(4),tmpflts(5),tmpflts(6))

  case default
    write(lout,"(a)") "LIMI> ERROR Aperture profile not identified for element '"//lnSplit(1)//"' value '"//lnSplit(2)//"'"
    iErr = .true.
    return

  end select

end subroutine aper_parseElement

subroutine aper_inputParsingDone

  use parpro
  use mod_settings

  implicit none

  integer i,ii

  if(limifound) then
    write(lout,"(a)") str_divLine
    write(lout,"(a)") ""
    write(lout,"(a)") "    DATA BLOCK APERTURE LIMITATIONS"
    write(lout,"(a)") ""
    if(st_quiet < 2) then ! Only print if quiet flag is less than 2 (default is 0).
      call dump_aperture_header(lout)
      do ii=1,il
        if ( kape(ii).ne.0 ) call dump_aperture_marker( lout, ii, 1 )
      enddo
    end if
    ! A.Mereghetti and P.Garcia Ortega, 12-06-2014
    ! echo precision for back-tracking
    if (lbacktracking) then
      write(lout,"(a,e22.15)") "    Back-tracking at aperture LIMIs is on, with precision [m]: ",bktpre
    else
      write(lout,"(a)")        "    No back-tracking, only checks at aperture LIMIs"
    end if
    write(lout,"(a)") ""
    write(lout,"(a)") str_divLine
  else
    write(lout,"(a)") "LIMI> No single element is assigned an aperture model."
    lbacktracking = .false.
  endif
  if (mxsec > 0) then
    do i=1,mxsec
      if (sLocDel(i) == zero) sLocDel(i)=bktpre
    enddo
  endif

end subroutine aper_inputParsingDone
! ================================================================================================ !
!  END APERTURE LIMITATIONS PARSING
! ================================================================================================ !

! ================================================================================================================================ !
!  Begin Checkpoint Restart
! ================================================================================================================================ !
#ifdef CR

! ================================================================================================================================ !
subroutine aper_crcheck_readdata(fileunit, readerr)

  implicit none

  integer, intent(in) :: fileunit
  logical, intent(out) :: readerr

  integer j

  read(fileunit,err=100,end=100) apefilepos_cr

  readerr = .false.
  return

100 continue
  readerr = .true.

end subroutine aper_crcheck_readdata

! ================================================================================================================================ !
subroutine aper_crcheck_positionFiles

  use crcoall
  use string_tools
  use mod_common
  use mod_units, only: f_open, f_close, f_requestUnit

  implicit none

  integer i,j
  logical lerror,lopen,err
  character(len=1024) arecord

  call f_requestUnit(losses_filename,losses_unit)
  write(93,"(a,i0)") "SIXTRACR> CRCHECK REPOSITIONING file of APERTURE LOSSES to apefilepos_cr = ",apefilepos_cr
  flush(93)

  inquire(unit=losses_unit, opened=lopen)
  if (.not. lopen) call f_open(unit=losses_unit,file=losses_filename,status='old',formatted=.true.,mode='rw',err=err)

  apefilepos = 0
  do j=1,apefilepos_cr
    read(losses_unit,'(a1024)',end=111,err=111,iostat=ierro) arecord
    apefilepos = apefilepos +1
  end do

  ! Crop aperture losses file
  ! This is not a FLUSH!
  endfile (losses_unit,iostat=ierro)

  ! Change from 'readwrite' to 'write'
  call f_close(losses_unit)
  call f_open(unit=losses_unit,file=losses_filename,status='old',formatted=.true.,mode='w+',err=err)

  return

111 continue
  write(93,*) 'SIXTRACR> APER_CRCHECK_POSITIONFILE *** ERROR *** reading file of APERTURE LOSSES, iostat=',ierro
  write(93,*) 'apefilepos=',apefilepos,' apefilepos_cr=',apefilepos_cr,' losses_unit=',losses_unit
  flush(93)
  write(lout,"(a)") "SIXTRACR> ERROR APER_CRCHECK_POSITIONFILES failure positioning file of APERTURE LOSSES"
  call prror(-1)

end subroutine aper_crcheck_positionFiles

! ================================================================================================================================ !
subroutine aper_crpoint(fileunit,lerror,ierro)

  implicit none

  integer, intent(in)    :: fileunit
  logical, intent(inout) :: lerror
  integer, intent(inout) :: ierro

  write(fileUnit,err=100,iostat=ierro) apefilepos
  endfile (fileunit,iostat=ierro)
  backspace (fileunit,iostat=ierro)
  return

100 continue
  lerror = .true.
  return

end subroutine aper_crpoint
! ================================================================================================================================ !

#endif
! ================================================================================================================================ !
!  End Checkpoint Restart
! ================================================================================================================================ !
end module aperture
