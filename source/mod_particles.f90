! ================================================================================================ !
!  SixTrack Particles Module
!  V.K. Berglyd Olsen, K.N. Sjobak, BE-ABP-HSS
!  Last modified: 2018-08-12
! ================================================================================================ !
module mod_particles

  use crcoall
  use floatPrecision
  use numerical_constants
  use physical_constants

  implicit none

contains

subroutine part_allocate
end subroutine part_allocate

subroutine part_expand
end subroutine part_expand

! ================================================================================================ !
!  K.N. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-08-13
!  Updates the relevant particle arrays after the particle energies or the reference energy changed.
!  If only the ejfv array has been changed, pass e0 to refEnergy
! ================================================================================================ !
subroutine part_updateEnergy(refEnergy)

  use mod_hions
  use mod_common
  use mod_commont
  use mod_commonmn

  implicit none

  real(kind=fPrec), intent(in) :: refEnergy

  real(kind=fPrec) e0o, e0fo
  integer          j

  ! Modify the reference particle
  e0o  = e0
  e0fo = e0f
  if (e0 /= refEnergy) then
    e0     = refEnergy
    e0f    = sqrt(e0**2 - nucm0**2)
    gammar = nucm0/e0
  end if

  ! Modify the Energy
  do j=1, napx
    ejfv(j)    = sqrt(ejv(j)**2 - nucm(j)**2)       ! Momentum [MeV/c]
    dpsv(j)    = (ejfv(j)*(nucm0/nucm(j))-e0f)/e0f  ! Delta_p/p0 = delta
    dpsv1(j)   = (dpsv(j)*c1e3)/(one + dpsv(j))
    dpd(j)     = one + dpsv(j)
    dpsq(j)    = sqrt(dpd(j))
    oidpsv(j)  = one/(one + dpsv(j))
    moidpsv(j) = mtc(j)/(one + dpsv(j))             ! Relative rigidity offset (mod_hions) [MV/c^2]
    rvv(j)     = (ejv(j)*e0f)/(e0*ejfv(j))          ! Beta_0 / beta(j)

    ! Also update sigmv with the new beta0 = e0f/e0
    if(e0 /= e0o) sigmv(j) = ((e0f*e0o)/(e0fo*e0))*sigmv(j)
  end do
  if(ithick == 1 .and. e0 /= e0o) call synuthck

end subroutine part_updateEnergy

end module mod_particles
