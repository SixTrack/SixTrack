! ================================================================================================ !
!  Dist Module for Collimation
!  Extracted from collimation module
! ================================================================================================ !
module coll_dist

  use parpro
  use floatPrecision
  use numerical_constants

  implicit none

  ! Other Settings
  real(kind=fPrec),         public,  save :: cdist_energy   = zero

  ! Twiss
  real(kind=fPrec),         public,  save :: cdist_alphaX   = zero
  real(kind=fPrec),         public,  save :: cdist_alphaY   = zero
  real(kind=fPrec),         public,  save :: cdist_betaX    = zero
  real(kind=fPrec),         public,  save :: cdist_betaY    = zero
  real(kind=fPrec),         private, save :: cdist_gammaX   = zero
  real(kind=fPrec),         private, save :: cdist_gammaY   = zero
  real(kind=fPrec),         public,  save :: cdist_emitX    = zero
  real(kind=fPrec),         public,  save :: cdist_emitY    = zero

  ! Distribution Settings
  real(kind=fPrec),         public,  save :: cdist_ampX     = zero
  real(kind=fPrec),         public,  save :: cdist_ampY     = zero
  real(kind=fPrec),         public,  save :: cdist_smearX   = zero
  real(kind=fPrec),         public,  save :: cdist_smearY   = zero
  real(kind=fPrec),         public,  save :: cdist_spreadE  = zero
  real(kind=fPrec),         public,  save :: cdist_bunchLen = zero

  character(len=mFileName), public,  save :: cdist_fileName = " "

contains

subroutine cdist_makeDist(distFormat)

  use crcoall

  integer, intent(in) :: distFormat

  cdist_gammaX = (one + cdist_alphaX**2)/cdist_betaX
  cdist_gammaY = (one + cdist_alphaY**2)/cdist_betaY

  select case(distFormat)
  case(0)
    return
  case(1)
    call cdist_makeDist_fmt1
  case default
    write(lout,"(a)") "COLLDIST> ERROR Unknown distribution format. Valid is 0 to 6, got ",distFormat
  end select

end subroutine cdist_makeDist

subroutine cdist_makeDist_fmt1

  use crcoall
  use mod_ranlux
  use mathlib_bouncer
  use mod_common, only : napx
  use mod_common_main, only : xv1, xv2, yv1, yv2, ejv, sigmv

  implicit none

  real(kind=fPrec) :: emitX, emitY, sigmaX, sigmaY
  integer :: j

! real(kind=fPrec) cdist_alphaX,cdist_betaX,cdist_emitX,myemitx,cdist_ampX,cdist_smearX, &
! &mygammax,cdist_alphaY,cdist_betaY,cdist_emitY,myemity,cdist_ampY,cdist_smearY,mygammay,   &
! &xsigmax,ysigmay,cdist_energy

  ! write(lout,"(a)") 'COLL> Generation of particle distribution Version 1:'
  ! write(lout,"(a)") 'COLL> This routine generates particles in phase space X/XP and Y/YP ellipses, as defined in the input '
  ! write(lout,"(a)") 'COLL> parameters. Distribution is flat in the band. X and Y are fully uncorrelated.'

  ! write(outlun,*)
  ! write(outlun,*) 'Generation of particle distribution Version 1'
  ! write(outlun,*)
  ! write(outlun,*) 'This routine generates particles in phase space'
  ! write(outlun,*) 'X/XP and Y/YP ellipses, as defined in the input'
  ! write(outlun,*) 'parameters. Distribution is flat in the band.'
  ! write(outlun,*) 'X and Y are fully uncorrelated.'
  ! write(outlun,*)
  ! write(outlun,*) 'INFO>  Number of particles   = ', napx
  ! write(outlun,*) 'INFO>  Av number of x sigmas = ', cdist_ampX
  ! write(outlun,*) 'INFO>  +- spread in x sigmas = ', cdist_smearX
  ! write(outlun,*) 'INFO>  Av number of y sigmas = ', cdist_ampY
  ! write(outlun,*) 'INFO>  +- spread in y sigmas = ', cdist_smearY
  ! write(outlun,*) 'INFO>  Nominal beam energy   = ', cdist_energy
  ! write(outlun,*) 'INFO>  Sigma_x0 = ', sqrt(cdist_betaX*cdist_emitX)
  ! write(outlun,*) 'INFO>  Sigma_y0 = ', sqrt(cdist_betaY*cdist_emitY)
  ! write(outlun,*) 'INFO>  Beta x   = ', cdist_betaX
  ! write(outlun,*) 'INFO>  Beta y   = ', cdist_betaY
  ! write(outlun,*) 'INFO>  Alpha x  = ', cdist_alphaX
  ! write(outlun,*) 'INFO>  Alpha y  = ', cdist_alphaY
  ! write(outlun,*)

  do j=1, napx
    emitX  = cdist_emitX*(cdist_ampX + ((two*real(rndm4()-half,fPrec))*cdist_smearX))**2
    sigmaX = sqrt(cdist_betaX*emitX)
    xv1(j) = sigmaX * sin_mb(twopi*real(rndm4(),fPrec))
    if(rndm4() > half) then
      yv1(j) = sqrt(emitX/cdist_betaX-xv1(j)**2/cdist_betaX**2)-(cdist_alphaX*xv1(j))/cdist_betaX
    else
      yv1(j) = -one*sqrt(emitX/cdist_betaX-xv1(j)**2/cdist_betaX**2)-(cdist_alphaX*xv1(j))/cdist_betaX
    end if

    emitY  = cdist_emitY*(cdist_ampY + ((two*real(rndm4()-half,fPrec))*cdist_smearY))**2
    sigmaY = sqrt(cdist_betaY*emitY)
    xv2(j) = sigmaY * sin_mb(twopi*real(rndm4(),fPrec))
    if(rndm4() > half) then
      yv2(j) = sqrt(emitY/cdist_betaY-xv2(j)**2/cdist_betaY**2)-(cdist_alphaY*xv2(j))/cdist_betaY
    else
      yv2(j) = -one*sqrt(emitY/cdist_betaY-xv2(j)**2/cdist_betaY**2)-(cdist_alphaY*xv2(j))/cdist_betaY
    end if

    ejv(j)   = cdist_energy
    sigmv(j) = zero

  ! !++  Dangerous stuff, just for the moment
  !  if (cut_input) then
  !    !0.1d-3 -> c1m4
  !    if((.not. (myy(j).lt.-0.008e-3_fPrec .and. myyp(j).lt. c1m4 .and.myyp(j).gt.zero) ) .and. &
  ! &        (.not. (myy(j).gt. 0.008e-3_fPrec .and. myyp(j).gt.-c1m4 .and.myyp(j).lt.zero) ) ) then
  !      j = j - 1
  !    end if
  !  end if
  end do

end subroutine cdist_makeDist_fmt1

end module coll_dist
