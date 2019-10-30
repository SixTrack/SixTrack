!>
!! @brief Module containing constants for beamgas part
!!
!<
module beamgascommon

  use floatPrecision
!       common to beamGasInit and beamGas
  integer, parameter :: bgmaxx=40000,bamount=1000
  integer bgmax, bgid, bgiddb(bgmaxx), ibgloc, pressID, njobs, njobthis, dpmjetevents
  real pressARRAY(2,bgmaxx)

!       bgParameters are s_null, n_null and n_here
!       these values are needed to know when enough particles are scattered
!       at a given point
!       s_null tells you how far the scattering process has gone so far
!       required that s_now > s_null (move s_null to s_now+ small delta
!       afterwards)
!       n_null tells you how many particles are scattered in previous gas
!       elements
!       n_here is a counter telling you how many particles are scattered
!       at this location
  real(kind=fPrec) bgParameters(3)
  real bgxpdb(bgmaxx),bgypdb(bgmaxx),bgEdb(bgmaxx)

  ! File unit variables
  integer, public, save :: bg_dpmJetUnit
  integer, public, save :: bg_scatterLocUnit
  integer, public, save :: bg_configUnit
  integer, public, save :: bg_pressProUnit
  integer, public, save :: bg_locLossesUnit

end module beamgascommon

module lorentzcommon
  use floatPrecision
  ! Common to lorentzBoost and createLorentzMatrix
  real(kind=fPrec) lorentzmatrix(4,4),new4MomCoord(4)
end module lorentzcommon

!>
!! \brief YIL subroutine beam gas scattering process.
!!
!!
!! This is a part of the inclusion of beam gas simulation for sixtrack
!!
!! Any "pressure element" (i.e. an element starting with\n
!! press) should call this function.\n
!! It will cycle through all primary particles and\n
!! scatter according to rules generated.
!!
!! @author Yngve Inntjore Levinsen <yngve.inntjore.levinsen@cern.ch>
!!
!! @date Last modified 21. Jul. 2010
!!
!!
!! \param myix The block ID number
!! \param mysecondary This is the array that tells which of the 64 particles are secondaries
!! \param totals This is the position around the ring, calculated from the start flag in fort.2
!! \param myenom This is an array with the energy of the 64 particles
!!
!! \warning This is the one-turn version
!! \return The subroutine does not return anything
!! \see thin6d, beamGasInit and rotateMatrix
!<
subroutine beamGas(myix, mysecondary, totals, myenom, turn, el_idx)
  use floatPrecision
  use numerical_constants
  use physical_constants
  use mathlib_bouncer

  use beamgascommon
  use lorentzcommon

  use crcoall
  use parpro
  use parbeam
  use mod_common
  use mod_common_track
  use mod_common_main
  use collimation, only : numeff, numeffdpop, part_abs_pos, part_abs_turn, nhit_stage

  implicit none

  real(kind=fPrec) myenom

! KNS (22/08/2017): The variable "secondary" is now visible
! from a +cd block as a COMMON variable;
! Renaming it to avoid conflict and to get the BEAMGAS version to compile.
! However, it has not been tested...
  integer mysecondary(npart)

! These are local subroutine stuff
  real(kind=fPrec) totals, oldCoordinates(5), protonmass, totMomentum, doLorentz, tmpPX, tmpPY

  integer choice,myix
  real(kind=fPrec) rotm(3,3), z(3),ztmp(3) ! the variable used to store rotation matrix and coordinates
! CHECK: Is ichar('0')=48 and so on for all systems??

  integer i,j,k,i_tmp

  integer turn, el_idx          ! KNS: turn and structure element idx

  pressID=0
  j=1

  do while (pressID.eq.0.and.j.le.bgmaxx)
    if ((pressARRAY(1,j).gt.(totals-c1m2)).and.(pressARRAY(1,j).le.(totals+c1m2))) pressID=j
    j=j+1
  end do

  if(pressID.eq.0) then
    write(lerr,"(a,e15.7)") "BEAMGAS> ERROR Couldn't find pressure marker at ",totals
    call prror
  end if

  doLorentz=0
  if((abs(yv1(1)).gt.3e-3_fPrec).or.(abs(yv2(1)).gt.3e-3_fPrec)) then ! do a Lorentz boost of DPMJET events
 !YIL warning: hardcoded mass of protons:
    protonmass=pmap
    doLorentz=1
    tmpPX=yv1(1) !don't think I can send array elements to functions??
    tmpPY=yv2(1)
    call createLorentzMatrix(myenom,tmpPX,tmpPY,protonmass)
  end if

  do j = 2,napx
  choice=0
  if ((nhit_stage(j).eq.0) .and. (part_abs_pos (j).eq.0 .and. part_abs_turn(j).eq.0) .and. (bgParameters(1).le.totals)) then
668 continue
! Warning: We round DOWN to the nearest integer at each
! location. It is needed in order not to run out of particles
! In generate_pmarkers.py the normalized sum is accordingly changed to 1
  if ((pressARRAY(2,pressID)*njobs*dpmjetevents).gt.(bgParameters(3)+1)) then
  bgParameters(3)=bgParameters(3)+1
  if (((bgParameters(2)+bgParameters(3)).gt.(dpmjetevents*njobthis)).and.((bgParameters(2)+bgParameters(3)).le.  &
    (dpmjetevents*(njobthis+1)))) then
!       The scattering id is increased by one for each interaction
  bgid=bgid+1

  do while (bgid.gt.bgiddb(ibgloc))
!    get to the right place in the lists
     ibgloc=ibgloc+1
  end do

  if(bgid.lt.bgiddb(ibgloc)) then ! no proton for this scattering event
!   check that this works correctly!!!!!!
    write(bg_locLossesUnit,*) partID(j),turn,totals,xv1(j),yv1(j),xv2(j),yv2(j),sigmv(j),(0-myenom)/myenom, &
      bgid+njobthis*dpmjetevents

!   writing down the scattering location information
    write(bg_scatterLocUnit,*) partID(j),turn,totals,xv1(j),yv1(j),xv2(j),yv2(j),sigmv(j),ejv(j),bgid+njobthis*dpmjetevents,&
      bgid,ejv(j),xv1(j),xv2(j),yv1(j),yv2(j)
!    part_abs(j) = 1
    part_abs_pos(j)  = el_idx
    part_abs_turn(j) = turn
  goto 669
  endif
  if(bgid.eq.bgiddb(ibgloc)) then ! a proton was found for this scattering event

    choice=ibgloc
! If several protons, the one with max energy is used
! THIS IS NECESSARY SINCE ALL PROTONS ARE READ INTO LIST!
  do i_tmp = 1,10
    if(bgiddb(choice).ne.bgiddb(ibgloc+i_tmp)) then
      exit
    endif
    if(bgEdb(ibgloc+i_tmp).gt.bgEdb(choice)) then
      choice=ibgloc+i_tmp
    end if
  end do

     oldCoordinates(1)=yv1(j)
     oldCoordinates(2)=yv2(j)
     oldCoordinates(3)=ejv(j)
     oldCoordinates(4)=xv1(j)
     oldCoordinates(5)=xv2(j)

     if(doLorentz.eq.one) then ! we need to boost the dpmjet event first:
     totMomentum=sqrt(bgEdb(choice)**2-(protonmass*c1m3)**2)
     tmpPX=bgxpdb(choice)*totMomentum
     tmpPY=bgypdb(choice)*totMomentum
     protonmass=protonmass*c1m3
     call lorentzBoost(tmpPX,tmpPY,totMomentum,protonmass)
     protonmass=protonmass*c1e3
!    This returns E,px,py,pz, need xp,yp
     totMomentum=sqrt(new4MomCoord(2)**2+new4MomCoord(3)**2+new4MomCoord(4)**2)
     call rotateMatrix(yv1(1),yv2(1),rotm) !we also need to "rotate back" before we're in the "same state"
     z(1) = (new4MomCoord(2)/totMomentum)
     z(2) = (new4MomCoord(3)/totMomentum)
     z(3) = (new4MomCoord(4)/totMomentum)
!    rotating the vector into the orbit reference system:
     z = matmul(rotm,z)
      if (z(3).eq.0) then
       write(lerr,"(a)")       "BEAMGAS> ERROR There is something wrong with your dpmjet event "
       write(lerr,"(a,e15.7)") "BEAMGAS>  * bgiddb(choice) = ",bgiddb(choice)
       write(lerr,"(a,e15.7)") "BEAMGAS>  * totMomentum    = ",totMomentum
       write(lerr,"(a,e15.7)") "BEAMGAS>  * new4MomCoord   = ",new4MomCoord
        call prror
      else
!      boosted xp event
       bgxpdb(choice) = z(1)
!      boosted yp event
       bgypdb(choice) = z(2)
       bgEdb(choice) = new4MomCoord(1) ! boosted energy
      endif
     endif ! doLorentz
     call rotateMatrix(yv1(j),yv2(j),rotm)
!    creating resulting vector [x,y,z] from dpmjet:
     z(1) = (bgxpdb(choice)) ! this is correct, since dpmjet gives xp=px/p and so on...
     z(2) = (bgypdb(choice))
     z(3) = sqrt(1-z(1)**2-z(2)**2)

!     rotating the vector into the orbit reference system:
      ztmp=z
      z=matmul(rotm,z)
! adding the angles to the yv vector:
  if (z(3).eq.0) then
    yv1(j) = acos_mb(zero)*c1e3
    yv2(j) = one*yv1(j)
  else
!   xp
    yv1(j) = atan_mb(z(1)/z(3))*c1e3
!   yp
    yv2(j) = atan_mb(z(2)/z(3))*c1e3
  endif
!    energy WARNING: I DO NOT KNOW ALL THE PLACES I NEED TO INSERT THE ENERGY????
     ejv(j) = bgEdb(choice)*c1e3
!YIL Copied this here, think these are all variables in need of an update
!++  Energy update, as recommended by Frank [comment from collimat part]
!
    ejfv(j)=sqrt(ejv(j)**2-nucm(j)**2)
    rvv(j)=(ejv(j)*e0f)/(e0*ejfv(j))
    dpsv(j)=(ejfv(j)*(nucm0/nucm(j))-e0f)/e0f
    oidpsv(j)=one/(one+dpsv(j))
    moidpsv(j)=mtc(j)/(one+dpsv(j))
    omoidpsv(j)=c1e3*((one-mtc(j))*oidpsv(j))
    dpsv1(j)=(dpsv(j)*c1e3)*oidpsv(j)
! writing down the scattering location information
  write(bg_scatterLocUnit,*) partID(j),turn,totals,xv1(j),        &
    yv1(j),xv2(j),yv2(j),sigmv(j),oldCoordinates(3),            &
    bgid+njobthis*dpmjetevents,bgid,ejv(j),oldCoordinates(4),      &
    oldCoordinates(5),oldCoordinates(1),oldCoordinates(2)
    mysecondary(j)=1

! if bgid.eq.bgiddb(ibgloc) end statement
  endif

! if njob correct range statement
  else if((bgParameters(2)+bgParameters(3)).le.(dpmjetevents*njobthis)) then
  goto 668
  endif
! if (pressARRAY(2,pressID)*njobs*dpmjetevents).gt.(bgParameters(3)+1)
  else
  bgParameters(1) = totals+c1m3
  bgParameters(2) = bgParameters(2)+bgParameters(3)
  bgParameters(3) = 0
! if (pressARRAY(2,pressID)*njobs*dpmjetevents).gt.(bgParameters(3)+1)
  endif
! check secondary if statement
  endif
! end j=1,napx statement
669 continue
  end do
end subroutine beamGas

!>
!! \brief YIL subroutine beam gas initiation.
!!
!! This function must be called during the initialization of\n
!! the simulation, if beam gas should be included.
!!
!! @author Yngve Inntjore Levinsen <yngve.inntjore.levinsen@cern.ch>
!!
!! @date Last updated 25. July 2009
!!
!! \warning This is the version with scattering only in first turn
!! \param myenom Needs to know nom. energy to know which events to skip
!! \return The subroutine does not return anything
!! \see beamGas and maincr
!! \todo pressure marker ID not used anymore, should be removed
!<
subroutine beamGasInit(myenom)
  use floatPrecision
  use numerical_constants
  use beamgascommon
  use crcoall
  use mod_units
  use parpro, only : npart

  implicit none

  integer check,j,i
  real(kind=fPrec) myenom,minenergy

  character(len=11) bg_var
  integer filereaderror, previousEvent,numberOfEvents
  real bg_val,ecutoff,pPOS,pVAL

  write(lout,"(a)") "BEAMGAS> Initialising"

  ! Get file units
  call f_requestUnit("dpmjet.eve",          bg_dpmJetUnit)
  call f_requestUnit("scatterLOC.txt",      bg_scatterLocUnit)
  call f_requestUnit("beamgas_config.txt",  bg_configUnit)
  call f_requestUnit("pressure_profile.txt",bg_pressProUnit)
  call f_requestUnit("localLOSSES.txt",     bg_locLossesUnit)

  open(bg_dpmJetUnit,file='dpmjet.eve')
  open(bg_scatterLocUnit,file='scatterLOC.txt')
  write(bg_scatterLocUnit,*)'# 1=name 2=turn 3=s 4=x 5=xp 6=y 7=yp 8=z 9=E',      &
    ' 10=eventID 11=dpmjetID 12=newEnergy 13=oldX 14=oldY 15=oldXP 16=oldYP'
  write(bg_scatterLocUnit,*)'# These are original coordinates of proton after impact, and old xp,yp'
  write(bg_scatterLocUnit,*)

! initialize pressure markers array
! DO THIS BEFORE OTHER STUFF, AS YOU MIGHT DO STUPID THINGS TO VARIABLES
! (LEARNED THE HARD WAY!!!)
  open(bg_configUnit,file='beamgas_config.txt')
  open(bg_pressProUnit,file='pressure_profile.txt')
  filereaderror=0
  do
     read(bg_configUnit,*,IOSTAT=filereaderror) bg_var, bg_val
  if (filereaderror.lt.0) then
! end of file
     exit
  else if (filereaderror.eq.0.and.bg_var.eq.'thisjob') then
    njobthis = bg_val
  else if (filereaderror.eq.0.and.bg_var.eq.'njobs') then
    njobs = bg_val
  else if (filereaderror.eq.0.and.bg_var.eq.'dpmjetev') then
    dpmjetevents = bg_val
  else if (filereaderror.eq.0.and.bg_var.eq.'ecutoff') then
    ecutoff = bg_val
  end if
  end do
  j=1
  do
    read(bg_pressProUnit,*,IOSTAT=filereaderror) pPOS, pVAL
    if (filereaderror.eq.0) then
      pressARRAY(1,j)=pPOS
      pressARRAY(2,j)=pVAL
      j=j+1
      if (j>bgmaxx) then
        write(lerr,"(a)") "BEAMGAS> ERROR Too many pressure markers!"
        call prror
      endif
    else if (filereaderror.lt.0) then
      ! means that end of file is reached
      exit
    else if (filereaderror.gt.0) then
      ! means that this line did not correspond to normal input
      ! do not need to perform anything (probably a comment line)
    end if
  end do
  do 1328 i = j,bgmaxx
    pressARRAY(1,i)=-one
    pressARRAY(2,i)=zero
1328  continue
! count the number of lines in dpmjet
  j=1
  previousEvent=0
  numberOfEvents=0
! Here you can set the energy acceptance (0.95 means at least 95% of nominal energy)
! 0.001 is because minenergy must be in GeV whereas myenom is in MeV
! Note to self: Remember to update this in batchrun.sh immediately! :)
  minenergy=ecutoff*myenom*c1m3
  filereaderror=0
  do
!    2212 is the proton id. We do not load other particles.
!    The other particles will be used to generate a complete file
!    afterwards.
!    ONLY LOAD PROTONS WITH ENERGY OFFSET BELOW 5%!!
     read(bg_dpmJetUnit,*,IOSTAT=filereaderror) bgiddb(j), check, bgxpdb(j), bgypdb(j), bgEdb(j)
     if (check.eq.2212.and.bgEdb(j).gt.minenergy) then
        if (bgiddb(j).ne.previousEvent) then
           previousEvent=bgiddb(j)
           numberOfEvents=numberOfEvents+1
        endif
        j=j+1
     endif
     if (filereaderror.lt.0) exit
!    If we have more events in the dpmjet file than
!    what we are supposed to simulate, we stop here...
     if (previousEvent.gt.dpmjetevents) exit
     if (numberOfEvents.gt.(bgmaxx-1)) then
     write(lerr,"(a)") "BEAMGAS> ERROR Too many dpmjet events!"
     call prror
  endif
  enddo
! number of lines in dpmjet - 1
  bgmax=j
  call f_close(bg_dpmJetUnit)
  write(lout,"(a,i0)") "BREAMGAS> Trackable events in dpmjet.eve: ",bgmax-1

  if (numberOfEvents.gt.npart) then
     write(lerr,"(2(a,i0))") "BEAMGAS> ERROR There were too many trackable events. Maximum for this run is ",npart,&
      " you generated ",numberOfEvents
     call prror
  endif

  write(lout,"(a,i0)") "BEAMGAS> This is job number: ", njobthis
  write(lout,"(a,i0)") "BEAMGAS> Total number of jobs: ", njobs
  write(lout,"(a,i0)") "BEAMGAS> Total number of particles in simulation: ", njobs*dpmjetevents
  call f_close(bg_configUnit)

  open(bg_locLossesUnit,file='localLOSSES.txt')
  write(bg_locLossesUnit,*) '# 1=name 2=turn 3=s 4=x 5=xp 6=y 7=yp 8=z 9=DE/E 10=CollisionID'
  write(bg_locLossesUnit,*) '# Note that name is not unique, but CollisionID is'
  write(bg_locLossesUnit,*) '# Note that s is particle coordinate, not bunch coordinate'

! YOU HAVE TO PUT THESE INITIALIZATIONS AT THE END OF THE ROUTINE
! FOR SOME STRANGE FORTRAN-REASON
  bgParameters(1)=zero
  bgParameters(2)=zero
  bgParameters(3)=zero

  bgid=0
  check=0
  ibgloc=1

end subroutine beamGasInit

!> \brief The routine returns a 3x3 rotation matrix for cartesian coordinates
!!
!! The function rotates cartesian coordinates based on an angle of the old and new\n
!! z-axis in the xz-plane (ax) and yz-plane (ay), given in milliradians. Typically\n
!! a particle with a small offset from the closed orbit (z-direction)
!!
!! @author Yngve Inntjore Levinsen <yngve.inntjore.levinsen@cern.ch>
!!
!! @date Last modified: 26. Mar. 2010
!!
!! \warning The angles of the particle should be in rad even though ax and ay
!! are in millirad! Dpmjet uses rad while sixtrack stores xp,yp in millirad!
!! \warning edit Mar10: changed sign of entire matrix, think it was wrong?
!! \param ax Angle in x-direction [millirad]
!! \param ay Angle in y-direction [millirad]
!! \param matrix 3x3 array which will contain the returned rotation matrix
!!
!! \return The subroutine returns a 3x3 rotation matrix
!! \see beamGas
!!
!<
subroutine rotateMatrix(ax,ay,matrix)
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  real(kind=fPrec) matrix(3,3)
  real(kind=fPrec) sinax, sinay, cosax, cosay
  real(kind=fPrec) ax,ay

  sinax = sin_mb(ax*c1m3)
  cosax = cos_mb(ax*c1m3)
  sinay = sin_mb(ay*c1m3)
  cosay = cos_mb(ay*c1m3)

  matrix(1,1)=cosax
  matrix(1,2)=-sinax*sinay
  matrix(1,3)=sinax*cosay

  matrix(2,1)=-sinax*sinay
  matrix(2,2)=cosay
  matrix(2,3)=sinay*cosax

  matrix(3,1)=-sinax
  matrix(3,2)=-sinay
  matrix(3,3)=cosax*cosay
end subroutine rotateMatrix

!>
!! \brief Performs Lorentz boost on a given coordinate set
!!
!! This code performs a Lorentz boost on the coordinates px,py,E
!! The Lorentz transfer matrix must be initialized first, using subroutine
!! createLorentzMatrix.
!!
!! @author Yngve Inntjore Levinsen <yngve.inntjore.levinsen@cern.ch>
!!
!! @date Last modified 28. July 2010
!!
!!
!! \param px [GeV] Momentum in x-direction
!! \param py [GeV] Momentum in y-direction
!! \param p [GeV] Total particle momentum
!! \param mass [GeV] Particle mass
!!
!! \return new4MomCoord will contain the 4-momentum coordinates after boost
!! \see createLorentzMatrix
!<
subroutine lorentzBoost(px,py,ptot,mass)
  use floatPrecision
  use numerical_constants
  use lorentzcommon

  implicit none

  real(kind=fPrec) px,py,ptot,mass
  real(kind=fPrec) oldcoord(4)

  integer i,j

  oldcoord(1)=sqrt(ptot**2+mass**2)
  oldcoord(2)=px
  oldcoord(3)=py
  oldcoord(4)=sqrt(ptot**2-px**2-py**2) ! E/c, px,py,pz
  do j=1,4
    new4MomCoord(j)=zero
  end do

! Matrix multiplication:
  do i=1,4
    do j=1,4
      new4MomCoord(i)=new4MomCoord(i)+lorentzmatrix(i,j)*oldcoord(j)
    end do
  end do
end subroutine lorentzBoost

!>
!! \brief Creates Lorentz transform matrix
!!
!! This subroutine sets up (or updates) the Lorentz matrix
!! used for Lorentz boost. This Lorentz boost is used for
!! implementing the crossing angle in distributions that
!! are coming from head-on collisions. Used for IR cross talk.
!! Because the boost shouldn't increase the energy of the distribution,
!! the entire matrix is divided by the gamma factor!
!!
!! @author Yngve Inntjore Levinsen <yngve.inntjore.levinsen@cern.ch>
!!
!! @date Last modified 28. July 2010
!!
!!
!! \param E [MeV] energy of the BEAM
!! \param xp [mrad] forward cosine in x-direction of the ORBIT coordinates
!! \param yp [mrad] forward cosine in y-direction of the ORBIT coordinates
!! \param mass [MeV] mass of the particle type
!!
!! \return Nothing
!! \warning Matrix divided by gamma factor!
!! \see lorentzBoost
!<
subroutine createLorentzMatrix(E,xp,yp,mass)
  use lorentzcommon
  use numerical_constants

  implicit none

  real(kind=fPrec) E,xp,yp,mass
! local variables:
  real(kind=fPrec) v0,gpart,p0,b(3),b2,b2inv,g

  integer i,j

  gpart=E/mass ! relativistic gamma for the particles
  p0=sqrt(E**2-mass**2)
  v0=p0/(gpart*mass)
  ! !         print v0
  b(1)=xp*c1m3*v0 ! relativistic beta...
  b(2)=yp*c1m3*v0
  b(3)=zero ! Assumed no movement of CM in longitudinal direction...
  b2=zero

  do j=1,3
    b2=b2+b(j)*b(j)
  end do

  if(b2>0) then
    b2inv=one/b2
  else
    b2inv=0
  endif

  g=one/sqrt(one-b2) ! relativistic gamma for the boost

  lorentzmatrix(1,1)=g /g

  do j=2,4
    lorentzmatrix(1,j)=g /g
    lorentzmatrix(1,j)=-b(j-1)*g /g
    lorentzmatrix(j,1)=-b(j-1)*g /g
    lorentzmatrix(j,j) = (one + (g-one)* b(j-1)**2*b2inv) /g
  end do


  lorentzmatrix(2,3) = ((g-1)* b(1)*b(2)*b2inv) /g
  lorentzmatrix(3,2) = (lorentzmatrix(2,3)) /g

  lorentzmatrix(2,4) = ((g-1)* b(1)*b(3)*b2inv) /g
  lorentzmatrix(4,2) = (lorentzmatrix(2,4)) /g

  lorentzmatrix(3,4) = ((g-1)* b(2)*b(3)*b2inv) /g
  lorentzmatrix(4,3) = (lorentzmatrix(3,4)) /g

end subroutine createLorentzMatrix

