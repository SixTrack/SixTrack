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
  use mod_commonmn

  implicit none

  integer j

  if(iclo6 == 2) then
    xv1(1:napx)   = xv1(1:napx)   +  clo6v(1,1:napx)
    yv1(1:napx)   = yv1(1:napx)   + clop6v(1,1:napx)
    xv2(1:napx)   = xv2(1:napx)   +  clo6v(2,1:napx)
    yv2(1:napx)   = yv2(1:napx)   + clop6v(2,1:napx)
    sigmv(1:napx) = sigmv(1:napx) +  clo6v(3,1:napx)
    dpsv(1:napx)  = dpsv(1:napx)  + clop6v(3,1:napx)
  else if(idfor == 0) then
    xv1(1:napx)   = xv1(1:napx)   +   clov(1,1:napx) * real(idz(1),fPrec)
    yv1(1:napx)   = yv1(1:napx)   +  clopv(1,1:napx) * real(idz(1),fPrec)
    xv2(1:napx)   = xv2(1:napx)   +   clov(2,1:napx) * real(idz(2),fPrec)
    yv2(1:napx)   = yv2(1:napx)   +  clopv(2,1:napx) * real(idz(2),fPrec)
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

  real(kind=fPrec) e0o, e0fo
  integer          j

  select case(refArray)
  case(1) ! Update from energy array
    do j=1, napx
      ejfv(j) = sqrt(ejv(j)**2 - nucm(j)**2)        ! Momentum [MeV/c]
      dpsv(j) = (ejfv(j)*(nucm0/nucm(j))-e0f)/e0f   ! Delta_p/p0 = delta
    end do
  case(2) ! Update from momentum array
    do j=1, napx
      ejv(j)  = sqrt(ejfv(j)**2 + nucm(j)**2)       ! Energy [MeV]
      dpsv(j) = (ejfv(j)*(nucm0/nucm(j))-e0f)/e0f   ! Delta_p/p0 = delta
    end do
  case(3) ! Update from delta array
    do j=1, napx
      ejfv(j) = ((nucm(j)/nucm0)*(dpsv(j)+one))*e0f ! Momentum [MeV/c]
      ejv(j)  = sqrt(ejfv(j)**2 + nucm(j)**2)       ! Energy [MeV]
    end do
  case default
    write(lout,"(a)") "PART> ERROR Internal error in part_updatePartEnergy"
    call prror(-1)
  end select

  ! Modify the Energy Dependent Arrays
  do j=1, napx
    dpsv1(j)    = (dpsv(j)*c1e3)/(one + dpsv(j))
    dpd(j)      = one + dpsv(j)                      ! For thick tracking
    dpsq(j)     = sqrt(dpd(j))                       ! For thick tracking
    oidpsv(j)   = one/(one + dpsv(j))
    moidpsv(j)  = mtc(j)/(one + dpsv(j))             ! Relative rigidity offset (mod_hions) [MV/c^2]
    omoidpsv(j) = ((one-mtc(j))*oidpsv(j))*c1e3
    rvv(j)      = (ejv(j)*e0f)/(e0*ejfv(j))          ! Beta_0 / beta(j)
  end do

  if(ithick == 1) call synuthck

end subroutine part_updatePartEnergy

end module mod_particles
