! ================================================================================================ !
!  SixTrack Particles Module
!  V.K. Berglyd Olsen, K.N. Sjobak, BE-ABP-HSS
!  Last modified: 2018-08-12
! ================================================================================================ !
module mod_particles

  use crcoall
  use floatPrecision

  implicit none

  logical, public, save :: part_isTracking = .false.

contains

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2019-01-11
!  Set the initial particle IDs (called before tracking, formerly in trauthin/trauthck)
! ================================================================================================ !
subroutine part_setParticleID
  use parpro
  use mod_common_main
  integer i
  do i=1,npart
    partID(i)   = i
    parentID(i) = i
  end do
end subroutine part_setParticleID

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

  implicit none

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
  call part_updatePartEnergy(3)

end subroutine part_applyClosedOrbit

! ================================================================================================ !
!  K.N. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-02
!  Updates the reference energy and momentum
! ================================================================================================ !
subroutine part_updateRefEnergy(refEnergy)

  use mod_hions
  use mod_common
  use mod_common_main
  use numerical_constants

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
  betrel = sqrt((one+gammar)*(one-gammar))

  ! Also update sigmv with the new beta0 = e0f/e0
  sigmv(1:napx) = ((e0f*e0o)/(e0fo*e0))*sigmv(1:napx)

  if(e0 <= pieni) then
    write(lout,"(a)") "PART> ERROR Reference energy ~= 0"
    call prror
  end if

  call part_updatePartEnergy(1)

end subroutine part_updateRefEnergy

! ================================================================================================ !
!  K.N. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2019-01-10
!  Updates the relevant particle arrays after the particle's energy, momentum or delta has changed.
! ================================================================================================ !
subroutine part_updatePartEnergy(refArray,updateAngle)

  use mod_hions
  use mod_common
  use mod_common_track
  use mod_common_main
  use numerical_constants

  implicit none

  integer,           intent(in) :: refArray
  logical, optional, intent(in) :: updateAngle

  logical :: doUpdateAngle = .false.

  if(part_isTracking .and. refArray /= 1) then
    write(lout,"(a)") "PART> ERROR During tracking, only energy updates are allowed in part_updatePartEnergy."
    call prror
  end if

  if(present(updateAngle)) then
    doUpdateAngle = updateAngle
  end if

  if(doUpdateAngle .and. refArray /= 2) then
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
    write(lout,"(a)") "PART> ERROR Internal error in part_updatePartEnergy"
    call prror
  end select

  ! Modify the Energy Dependent Arrays
  dpsv1(1:napx)    = (dpsv(1:napx)*c1e3)/(one + dpsv(1:napx))
  dpd(1:napx)      = one + dpsv(1:napx)                      ! For thick tracking
  dpsq(1:napx)     = sqrt(dpd(1:napx))                       ! For thick tracking
  oidpsv(1:napx)   = one/(one + dpsv(1:napx))
  moidpsv(1:napx)  = mtc(1:napx)/(one + dpsv(1:napx))        ! Relative rigidity offset (mod_hions) [MV/c^2]
  omoidpsv(1:napx) = ((one-mtc(1:napx))*oidpsv(1:napx))*c1e3
  rvv(1:napx)      = (ejv(1:napx)*e0f)/(e0*ejfv(1:napx))     ! Beta_0 / beta(j)

  if(doUpdateAngle) then ! Update particle angles
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
! ================================================================================================ !
subroutine part_writeState(theState)

  use, intrinsic :: iso_fortran_env, only : int16, int32, real64

  use mod_units
  use parpro
  use mod_hions
  use mod_common
  use mod_common_main
  use mod_settings
  use string_tools

  implicit none

  integer, intent(in) :: theState

  character(len=225) :: roundBuf
  character(len=17)  :: fileName
  integer            :: fileUnit, j, k, iDummy
  logical            :: rErr, isPrim, isBin, noIons

  if(theState == 0) then
    if(st_initialState == 0) return ! No dump was requested in fort.3
    isBin    = st_initialState == 1 .or. st_initialState == 3
    noIons   = st_initialState == 1 .or. st_initialState == 2
    fileName = "initial_state"
  elseif(theState == 1) then
    if(st_finalState == 0) return ! No dump was requested in fort.3
    isBin    = st_finalState == 1 .or. st_finalState == 3
    noIons   = st_finalState == 1 .or. st_finalState == 2
    fileName = "final_state"
  else
    ! Nothing to do
    return
  end if

  if(isBin) then

    fileName = trim(fileName)//".bin"
    call f_requestUnit(fileName, fileUnit)
    call f_open(unit=fileUnit,file=fileName,formatted=.false.,mode="w",status="replace",access="stream")

    write(fileUnit) int(imc,  kind=int32)
    write(fileUnit) int(napx, kind=int32)
    write(fileUnit) int(napxo,kind=int32)
    write(fileUnit) int(npart,kind=int32)

    iDummy = 0

    do j=1,npart
      isPrim = partID(j) <= napxo
      write(fileUnit)     int(  partID(j), kind=int32)
      write(fileUnit)     int(parentID(j), kind=int32)
      write(fileUnit) logical(  llostp(j), kind=int32)
      write(fileUnit) logical(  isPrim,    kind=int32)
      write(fileUnit)    real(     xv1(j), kind=real64)
      write(fileUnit)    real(     xv2(j), kind=real64)
      write(fileUnit)    real(     yv1(j), kind=real64)
      write(fileUnit)    real(     yv2(j), kind=real64)
      write(fileUnit)    real(   sigmv(j), kind=real64)
      write(fileUnit)    real(    dpsv(j), kind=real64)
      write(fileUnit)    real(    ejfv(j), kind=real64)
      write(fileUnit)    real(     ejv(j), kind=real64)
      if(noIons) cycle ! Skip the ion columns
      write(fileUnit)    real(    nucm(j), kind=real64)
      write(fileUnit)     int(     naa(j), kind=int16)
      write(fileUnit)     int(     nzz(j), kind=int16)
    ! write(fileUnit)     int(     nqq(j), kind=int16) ! Not implemented yet
      write(fileUnit)     int(     iDummy, kind=int32) ! Pad to n x 64 bit
    end do

    call f_close(fileUnit)

  else

    fileName = trim(fileName)//".dat"
    call f_requestUnit(fileName, fileUnit)
    call f_open(unit=fileUnit,file=fileName,formatted=.true.,mode="w",status="replace")

    write(fileUnit,"(a,i0)") "# imc   = ",imc
    write(fileUnit,"(a,i0)") "# napx  = ",napx
    write(fileUnit,"(a,i0)") "# napxo = ",napxo
    write(fileUnit,"(a,i0)") "# npart = ",npart
    if(noIons) then
      write(fileUnit,"(a1,a7,1x,a8,2(1x,a4),8(1x,a24))") &
        "#","partID","parentID","lost","prim","x","y","xp","yp","sigma","dp","p","e"
    else
      write(fileUnit,"(a1,a7,1x,a8,2(1x,a4),9(1x,a24),2(1x,a4))") &
        "#","partID","parentID","lost","prim","x","y","xp","yp","sigma","dp","p","e","mass","A","Z"
    end if

    do j=1,npart
      roundBuf = " "
      isPrim = partID(j) <= napxo
      call chr_fromReal(xv1(j),  roundBuf(  2:25 ),17,3,rErr)
      call chr_fromReal(xv2(j),  roundBuf( 27:50 ),17,3,rErr)
      call chr_fromReal(yv1(j),  roundBuf( 52:75 ),17,3,rErr)
      call chr_fromReal(yv2(j),  roundBuf( 77:100),17,3,rErr)
      call chr_fromReal(sigmv(j),roundBuf(102:125),17,3,rErr)
      call chr_fromReal(dpsv(j), roundBuf(127:150),17,3,rErr)
      call chr_fromReal(ejfv(j), roundBuf(152:175),17,3,rErr)
      call chr_fromReal(ejv(j),  roundBuf(177:200),17,3,rErr)
      if(noIons) then
        write(fileUnit, "(i8,1x,i8,2(1x,l4),a200)") &
          partID(j),parentID(j),llostp(j),isPrim,roundBuf(1:200)
      else
        call chr_fromReal(nucm(j), roundBuf(202:225),17,3,rErr)
        write(fileUnit, "(i8,1x,i8,2(1x,l4),a225,2(1x,i4))") &
          partID(j),parentID(j),llostp(j),isPrim,roundBuf(1:225),naa(j),nzz(j)
      end if
    end do

    call f_close(fileUnit)

  end if

end subroutine part_writeState

end module mod_particles
