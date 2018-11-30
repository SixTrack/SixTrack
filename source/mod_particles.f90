! ================================================================================================ !
!  SixTrack Particles Module
!  V.K. Berglyd Olsen, K.N. Sjobak, BE-ABP-HSS
!  Last modified: 2018-08-12
! ================================================================================================ !
module mod_particles

  use crcoall
  use floatPrecision

  implicit none

contains

subroutine part_allocate
end subroutine part_allocate

subroutine part_expand
end subroutine part_expand

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-10-31
!  Moved from maincr. Applies the closed orbit correction to the particle coordinates.
! ================================================================================================ !
subroutine part_applyClosedOrbit

  use mod_common
  use mod_commons
  use mod_commont
  use mod_commonmn

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
  use mod_commonmn
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
    call prror(-1)
  end if

  call part_updatePartEnergy(1)

end subroutine part_updateRefEnergy

! ================================================================================================ !
!  K.N. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-02
!  Updates the relevant particle arrays after the particle's energy, momentum or delta has changed.
! ================================================================================================ !
subroutine part_updatePartEnergy(refArray)

  use mod_hions
  use mod_common
  use mod_commont
  use mod_commonmn
  use numerical_constants

  implicit none

  integer, intent(in) :: refArray

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
    call prror(-1)
  end select

  ! Modify the Energy Dependent Arrays
  dpsv1(1:napx)    = (dpsv(1:napx)*c1e3)/(one + dpsv(1:napx))
  dpd(1:napx)      = one + dpsv(1:napx)                      ! For thick tracking
  dpsq(1:napx)     = sqrt(dpd(1:napx))                       ! For thick tracking
  oidpsv(1:napx)   = one/(one + dpsv(1:napx))
  moidpsv(1:napx)  = mtc(1:napx)/(one + dpsv(1:napx))        ! Relative rigidity offset (mod_hions) [MV/c^2]
  omoidpsv(1:napx) = ((one-mtc(1:napx))*oidpsv(1:napx))*c1e3
  rvv(1:napx)      = (ejv(1:napx)*e0f)/(e0*ejfv(1:napx))     ! Beta_0 / beta(j)

  if(ithick == 1) call synuthck

end subroutine part_updatePartEnergy

end module mod_particles
