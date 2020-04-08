! ================================================================================================ !
!  SixTrack Particles Module
!  V.K. Berglyd Olsen, K.N. Sjobak, BE-ABP-HSS
!  Last modified: 2018-08-12
! ================================================================================================ !
module mod_particles

  use floatPrecision

  implicit none

  logical, public, save :: part_isTracking = .false.

contains

! ================================================================================================ !
!  Set Initial Particle IDs
! ~~~~~~~~~~~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-01-11
!  Updated: 2019-08-13
!      Sets the default particle ID to run from 1 to npart, and each particle being its own parent,
!  that is, it is a primary particle. The pairID signifies its original paired particle. Each pair
!  must consist of one particle with an even ID and one with an odd ID.
! ================================================================================================ !
subroutine part_setParticleID
  use parpro
  use mod_common_main
  integer j
  do j=1,npart
    partID(j)   = j
    parentID(j) = j
  end do
  call part_setPairID
end subroutine part_setParticleID

! ================================================================================================ !
!  Set Initial Pair IDs
! ~~~~~~~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-08-13
!  Updated: 2019-08-13
!      This routine is called if the partID and parentID are already set, and pairID needs to be
!  populated. Each particle has a pairID associated with it and whether it's the first or second
!  particle of the pair.
! ================================================================================================ !
subroutine part_setPairID
  use parpro
  use mod_common_main
  integer j
  do j=1,npart
#ifndef G4COLLIMATION
    pairID(1,j) = (j+1)/2    ! The pairID of particle j
    pairID(2,j) = 2-mod(j,2) ! Either particle 1 or 2 of the pair
#else
!with g4, treat all particles as secondaries
    pairID(1,j) = 0
    pairID(2,j) = 0
#endif
  end do
  call updatePairMap
end subroutine part_setPairID

! ================================================================================================ !
!  Get Original Particle Index
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-08-21
!  Updated: 2019-08-21
!      This function returns the original particle index based on its pairID. This is a substiture
!  for the old reverse lookup map that is now used as a particle ID array.
! ================================================================================================ !
pure integer function part_getOrigIndex(j)
  use mod_common_main, only : pairID
  integer, intent(in) :: j
  part_getOrigIndex = (pairID(1,j) - 1)*2 + pairID(2,j)
end function part_getOrigIndex

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-10-31
!  Moved from maincr. Applies the closed orbit correction to the particle coordinates.
! ================================================================================================ !
subroutine part_applyClosedOrbit

  use mod_common
  use mod_commons
  use mod_common_track
  use mod_common_main

  if(iclo6 == 2) then
    xv1(1:napx)   = xv1(1:napx)   +  clo6(1)
    yv1(1:napx)   = yv1(1:napx)   + clop6(1)
    xv2(1:napx)   = xv2(1:napx)   +  clo6(2)
    yv2(1:napx)   = yv2(1:napx)   + clop6(2)
    sigmv(1:napx) = sigmv(1:napx) +  clo6(3)
    dpsv(1:napx)  = dpsv(1:napx)  + clop6(3)
  else if(idfor == 0) then
    xv1(1:napx)   = xv1(1:napx)   +   clo(1)*real(idz(1),fPrec)
    yv1(1:napx)   = yv1(1:napx)   +  clop(1)*real(idz(1),fPrec)
    xv2(1:napx)   = xv2(1:napx)   +   clo(2)*real(idz(2),fPrec)
    yv2(1:napx)   = yv2(1:napx)   +  clop(2)*real(idz(2),fPrec)
  end if
  call part_updatePartEnergy(3,.false.)

end subroutine part_applyClosedOrbit

! ================================================================================================ !
!  K.N. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-02
!  Updates the reference energy and momentum
! ================================================================================================ !
subroutine part_updateRefEnergy(refEnergy)

  use crcoall
  use mod_common
  use mod_common_main
  use numerical_constants, only : one, c1m6
  use physical_constants, only: clight

  implicit none

  real(kind=fPrec), intent(in) :: refEnergy

  real(kind=fPrec) e0o, e0fo

  if(e0 == refEnergy) return

  ! Save previous values
  e0o    = e0
  e0fo   = e0f

  ! Modify the reference particle
  e0     = refEnergy
  e0f    = sqrt(e0**2 - nucm0**2)
  gammar = nucm0/e0
  gamma0 = e0/nucm0
  beta0  = sqrt((one+gammar)*(one-gammar))
  brho   = (e0f/(clight*c1m6))/qq0

  ! Also update sigmv with the new beta0 = e0f/e0
  sigmv(1:napx) = ((e0f*e0o)/(e0fo*e0))*sigmv(1:napx)

  if(e0 <= pieni) then
    write(lerr,"(a)") "PART> ERROR Reference energy ~= 0"
    call prror
  end if

  call part_updatePartEnergy(1,.false.)

end subroutine part_updateRefEnergy

! ================================================================================================ !
!  K.N. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2019-01-10
!  Updates the relevant particle arrays after the particle's energy, momentum or delta has changed.
! ================================================================================================ !
subroutine part_updatePartEnergy(refArray, updateAngle)

  use crcoall
  use mod_common
  use mod_common_track
  use mod_common_main
  use numerical_constants

  implicit none

  integer, intent(in) :: refArray
  logical, intent(in) :: updateAngle

  !if(part_isTracking .and. refArray /= 1) then
  !  write(lerr,"(a)") "PART> ERROR During tracking, only energy updates are allowed in part_updatePartEnergy."
  !  call prror
  !end if

  if(updateAngle .and. refArray /= 2) then
    ! If momentum is updated before the call, then ejf0v must be too
    ejf0v(1:napx) = ejfv(1:napx)
  end if

  select case(refArray)
  case(1) ! Update from energy array
    ejfv(1:napx) = sqrt(ejv(1:napx)**2 - nucm(1:napx)**2)        ! Momentum [MeV/c]
    dpsv(1:napx) = (ejfv(1:napx)*(nucm0/nucm(1:napx))-e0f)/e0f   ! Delta_p/p0 = delta
  case(2) ! Update from momentum array
    ejv(1:napx)  = sqrt(ejfv(1:napx)**2 + nucm(1:napx)**2)       ! Energy [MeV]
    dpsv(1:napx) = (ejfv(1:napx)*(nucm0/nucm(1:napx))-e0f)/e0f   ! Delta_p/p0 = delta
  case(3) ! Update from delta array
    ejfv(1:napx) = ((nucm(1:napx)/nucm0)*(dpsv(1:napx)+one))*e0f ! Momentum [MeV/c]
    ejv(1:napx)  = sqrt(ejfv(1:napx)**2 + nucm(1:napx)**2)       ! Energy [MeV]
  case default
    write(lerr,"(a)") "PART> ERROR Internal error in part_updatePartEnergy"
    call prror
  end select

  ! Modify the Energy Dependent Arrays
  ! Keep in sync with checpoint/restart crstart
  ! dpsv1(1:napx)    = (dpsv(1:napx)*c1e3)/(one + dpsv(1:napx))
  dpd(1:napx)      = one + dpsv(1:napx)                      ! For thick tracking
  dpsq(1:napx)     = sqrt(dpd(1:napx))                       ! For thick tracking
  oidpsv(1:napx)   = one/(one + dpsv(1:napx))
  moidpsv(1:napx)  = mtc(1:napx)/(one + dpsv(1:napx))        ! Relative rigidity offset [MV/c^2]
  omoidpsv(1:napx) = ((one-mtc(1:napx))*oidpsv(1:napx))*c1e3
  dpsv1(1:napx)    = (dpsv(1:napx)*c1e3)*oidpsv(1:napx)
  rvv(1:napx)      = (ejv(1:napx)*e0f)/(e0*ejfv(1:napx))     ! Beta_0 / beta(j)

  if(updateAngle) then ! Update particle angles
    yv1(1:napx)    = (ejf0v(1:napx)/ejfv(1:napx))*yv1(1:napx)
    yv2(1:napx)    = (ejf0v(1:napx)/ejfv(1:napx))*yv2(1:napx)
  end if

  if(ithick == 1) call synuthck

end subroutine part_updatePartEnergy

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created:  2018-11-28
!  Modified: 2019-01-31
!  Dumps the state of the particle arrays to a binary or text file.
!  Note: This subroutine assumes the default integer size is 32 bit for the binary files. This is
!        currently the case for the compilers supported by SixTrack. If this changes, or is
!        explicitly otherwise requested in the buold files, tests based in binary particle state
!        files will fail.
! ================================================================================================ !
subroutine part_writeState(fileName, isText, withIons)

  use parpro
  use mod_units
  use mod_common
  use mod_common_main
  use mod_common_track
  use mod_settings
  use string_tools
  use numerical_constants

  use, intrinsic :: iso_fortran_env, only : int16, int32

  character(len=*), intent(in) :: fileName
  logical,          intent(in) :: isText
  logical,          intent(in) :: withIons

  character(len=225)  :: roundBuf
  real(kind=fPrec)    :: tmpTas(6,6)
  integer(kind=int16) :: iIons
  integer             :: fileUnit, i, j, iPrim, iLost
  logical             :: rErr, isPrim

  call f_requestUnit(fileName, fileUnit)

  ! Copy the tas matrix and reverse the scaling from umlauda
  tmpTas(1:6,1:6) = tas(1:6,1:6)
  tmpTas(1:5,6)   = tmpTas(1:5,6)*c1m3
  tmpTas(6,1:5)   = tmpTas(6,1:5)*c1e3

  if(isText) then
    call f_open(unit=fileUnit,file=fileName,formatted=.true.,mode="w",status="replace")

    write(fileUnit,"(a)")    "# Tracking"
    write(fileUnit,"(a,i0)") "# NPart Start     = ",napxo
    write(fileUnit,"(a,i0)") "# NPart End       = ",napx
    write(fileUnit,"(a,i0)") "# NPart Allocated = ",npart
    write(fileUnit,"(a,i0)") "# NTurns          = ",numl

    call chr_fromReal(nucm0, roundBuf( 2:25),17,3,rErr)
    call chr_fromReal(e0,    roundBuf(27:50),17,3,rErr)
    call chr_fromReal(e0f,   roundBuf(52:75),17,3,rErr)

    write(fileUnit,"(a)")     "#"
    write(fileUnit,"(a)")     "# Reference Particle"
    write(fileUnit,"(a,a24)") "# Mass [MeV]      = ",roundBuf( 2:25)
    write(fileUnit,"(a,a24)") "# Energy [MeV]    = ",roundBuf(27:50)
    write(fileUnit,"(a,a24)") "# Momentum [MeV]  = ",roundBuf(52:75)
    write(fileUnit,"(a,i0)")  "# Atomic Mass     = ",aa0
    write(fileUnit,"(a,i0)")  "# Atomic Number   = ",zz0
    write(fileUnit,"(a,i0)")  "# Charge          = ",qq0

    write(fileUnit,"(a)") "#"
    write(fileUnit,"(a)") "# Closed Orbit [x, xp, y, yp, sigma, dp]"

    roundBuf = " "
    call chr_fromReal(clo(1),  roundBuf(  2:25 ),17,3,rErr)
    call chr_fromReal(clop(1), roundBuf( 27:50 ),17,3,rErr)
    call chr_fromReal(clo(2),  roundBuf( 52:75 ),17,3,rErr)
    call chr_fromReal(clop(2), roundBuf( 77:100),17,3,rErr)
    write(fileUnit,"(a,a100)") "# 4D Closed Orbit =",roundBuf(1:100)

    roundBuf = " "
    call chr_fromReal(clo6(1),  roundBuf(  2:25 ),17,3,rErr)
    call chr_fromReal(clop6(1), roundBuf( 27:50 ),17,3,rErr)
    call chr_fromReal(clo6(2),  roundBuf( 52:75 ),17,3,rErr)
    call chr_fromReal(clop6(2), roundBuf( 77:100),17,3,rErr)
    call chr_fromReal(clo6(3),  roundBuf(102:125),17,3,rErr)
    call chr_fromReal(clop6(3), roundBuf(127:150),17,3,rErr)
    write(fileUnit,"(a,a150)") "# 6D Closed Orbit =",roundBuf(1:150)

    write(fileUnit,"(a)") "#"
    roundBuf = " "
    call chr_fromReal(qwc(1), roundBuf(  2:25 ),17,3,rErr)
    call chr_fromReal(qwc(2), roundBuf( 27:50 ),17,3,rErr)
    call chr_fromReal(qwc(3), roundBuf( 52:75 ),17,3,rErr)
    write(fileUnit,"(a,a75)") "# Tune            =",roundBuf(1:75)
    do i=1,6
      roundBuf = " "
      call chr_fromReal(tmpTas(i,1), roundBuf(  2:25 ),17,3,rErr)
      call chr_fromReal(tmpTas(i,2), roundBuf( 27:50 ),17,3,rErr)
      call chr_fromReal(tmpTas(i,3), roundBuf( 52:75 ),17,3,rErr)
      call chr_fromReal(tmpTas(i,4), roundBuf( 77:100),17,3,rErr)
      call chr_fromReal(tmpTas(i,5), roundBuf(102:125),17,3,rErr)
      call chr_fromReal(tmpTas(i,6), roundBuf(127:150),17,3,rErr)
      write(fileUnit,"(a,i0,a,a150)") "# TAS(",i,",1:6)      =",roundBuf(1:150)
    end do

    write(fileUnit,"(a)") "#"
    if(withIons) then
      write(fileUnit,"(a1,a7,1x,a8,1x,a10,2(1x,a4),1x,a8,9(1x,a24),3(1x,a4),1x,a11)") &
        "#","partID","parentID","pairID","lost","prim","turns","x","y","xp","yp","sigma","dp","p","e","mass","A","Z","Q","PDGid"
    else
      write(fileUnit,"(a1,a7,1x,a8,1x,a10,2(1x,a4),1x,a8,8(1x,a24))") &
        "#","partID","parentID","pairID","lost","prim","turns","x","y","xp","yp","sigma","dp","p","e"
    end if
    do j=1,npart
      roundBuf = " "
      isPrim = partID(j) == parentID(j)
      call chr_fromReal(xv1(j),  roundBuf(  2:25 ),17,3,rErr)
      call chr_fromReal(xv2(j),  roundBuf( 27:50 ),17,3,rErr)
      call chr_fromReal(yv1(j),  roundBuf( 52:75 ),17,3,rErr)
      call chr_fromReal(yv2(j),  roundBuf( 77:100),17,3,rErr)
      call chr_fromReal(sigmv(j),roundBuf(102:125),17,3,rErr)
      call chr_fromReal(dpsv(j), roundBuf(127:150),17,3,rErr)
      call chr_fromReal(ejfv(j), roundBuf(152:175),17,3,rErr)
      call chr_fromReal(ejv(j),  roundBuf(177:200),17,3,rErr)
      if(withIons) then
        call chr_fromReal(nucm(j), roundBuf(202:225),17,3,rErr)
        write(fileUnit, "(i8,1x,i8,1x,i8,a1,i1,2(1x,l4),1x,i8,a225,3(1x,i4),1x,i11)") &
          partID(j),parentID(j),pairID(1,j),".",pairID(2,j),llostp(j),isPrim,min(numxv(j),numx),roundBuf(1:225),&
          naa(j),nzz(j),nqq(j),pdgid(j)
      else
        write(fileUnit, "(i8,1x,i8,1x,i8,a1,i1,2(1x,l4),1x,i8,a200)") &
          partID(j),parentID(j),pairID(1,j),".",pairID(2,j),llostp(j),isPrim,min(numxv(j),numx),roundBuf(1:200)
      end if
    end do
  else
    ! Format
    ! Header: 440 bytes
    ! Record:  96 bytes
    ! + Ions:  24 bytes
    call f_open(unit=fileUnit,file=fileName,formatted=.false.,mode="w",status="replace",access="stream")
    if(withIons) then
      iIons = 1
    else
      iIons = 0
    end if
    write(fileUnit) napxo,napx,npart,numl                               ! 4x32bit
    write(fileUnit) nucm0,e0,e0f,aa0,zz0,qq0,iIons                      ! 3x64bit + 4x16bit
    write(fileUnit) clo(1),clop(1),clo(2),clop(2)                       ! 4x64bit
    write(fileUnit) clo6(1),clop6(1),clo6(2),clop6(2),clo6(3),clop6(3)  ! 6x64bit
    write(fileUnit) qwc(1),qwc(2),qwc(3)                                ! 3x64bit
    write(fileUnit) tmpTas(1,1:6)                                       ! 6x64bit
    write(fileUnit) tmpTas(2,1:6)                                       ! 6x64bit
    write(fileUnit) tmpTas(3,1:6)                                       ! 6x64bit
    write(fileUnit) tmpTas(4,1:6)                                       ! 6x64bit
    write(fileUnit) tmpTas(5,1:6)                                       ! 6x64bit
    write(fileUnit) tmpTas(6,1:6)                                       ! 6x64bit
    do j=1,npart
      ! These have to be set explicitly as ifort converts logical to integer differently than gfortran and nagfor
      if(partID(j) == parentID(j)) then
        iPrim = 1
      else
        iPrim = 0
      end if
      if(llostp(j)) then
        iLost = 1
      else
        iLost = 0
      end if
      write(fileUnit) partID(j),parentID(j),pairID(1,j),pairID(2,j) ! 4x32 bit
      write(fileUnit) iLost,iPrim,min(numxv(j),numx),0_int32        ! 4x32 bit
      write(fileUnit) xv1(j),xv2(j),yv1(j),yv2(j)                   ! 4x64 bit
      write(fileUnit) sigmv(j),dpsv(j),ejfv(j),ejv(j)               ! 4x64 bit
      if(withIons) then
        write(fileUnit) nucm(j),naa(j),nzz(j),nqq(j),pdgid(j)       ! 64 bit + 3x16 + 32 bit
        write(fileUnit) 0_int16, 0_int32                            ! Pad to nearest 64 bit
      end if
    end do
  end if

  call f_freeUnit(fileUnit)

end subroutine part_writeState

end module mod_particles
