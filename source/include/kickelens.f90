! start include/kickelens.f90
! 1) apply offset of e-lens
!    xelens = x(proton) - elens_offset_x
!    yelens = y(proton) - elens_offset_y
xelens=xv1(j)-elens_offset_x(ielens(ix))
yelens=xv2(j)-elens_offset_y(ielens(ix))
! 2) calculate radius
!    radial position of main beam relative to center of elens beam
!    -> internal parameters to calculate kick:
!    rrelens = sqrt(xelens**2+yelens**2)
rrelens=sqrt((xelens)**2+(yelens**2))
! 3) calculate kick
!    shape function: spatial charge density depends on type:
!    0        if r < R1
!    frrelens if R1 < r < R2
!    1        if r > R2
if (rrelens.gt.elens_r1(ielens(ix))) then ! rrelens <= r1 -> no kick from elens
  if (rrelens.lt.elens_r2(ielens(ix))) then ! r1 < rrelens < r2
    select case (elens_type(ielens(ix)))
    case (1)
      ! UNIFORM: eLens with uniform profile
      ! formula: (r^2-r1^2)/(r2^2-r1^2)
      frrelens=( (rrelens+elens_r1(ielens(ix)))*(rrelens-elens_r1(ielens(ix))) )/elens_geo_norm(ielens(ix))
    case (2)
      ! GAUSSIAN: eLens with Gaussian profile
      ! formula: (exp(-r1^2/2sig^2)-exp(-r^2/2sig^2))/(exp(-r1^2/2sig^2)-exp(-r2^2/2sig^2))
      frrelens=( exp_mb(-0.5*(elens_r1(ielens(ix))/elens_sig(ielens(ix)))**2)    &
                -exp_mb(-0.5*(rrelens             /elens_sig(ielens(ix)))**2) )/ &
                              elens_geo_norm(ielens(ix))
    case (3)
      ! RADIAL PROFILE: eLens with radial profile as from file
      ! formula: (cumul_J(r)-cumul_J(r1))/(cumul_J(r2)-cumul_J(r1))
      frrelens=(lininterp( rrelens, &
            elens_radial_profile_R(0:elens_radial_profile_nPoints(elens_iRadial(ielens(ix))),elens_iRadial(ielens(ix))), &
            elens_radial_profile_J(0:elens_radial_profile_nPoints(elens_iRadial(ielens(ix))),elens_iRadial(ielens(ix))), &
            elens_radial_profile_nPoints(elens_iRadial(ielens(ix)))+1)-elens_radial_fr1(ielens(ix)) )/elens_geo_norm(ielens(ix))
    case default
      write(lout,"(a,i0)") "ELENS> ERROR in kickelens: elens_type=",elens_type(ielens(ix))," not recognized. "
      write(lout,"(a)")    "ELENS>       Possible values for type are: 1 and 2"
      call prror(-1)
    end select
  else ! r1 < r2 <= rrelens
    frrelens = one
  endif
  ! 'radial kick'
  frrelens = ((((elens_theta_r2(ielens(ix))*elens_r2(ielens(ix)))/rrelens)*frrelens)*oidpsv(j))*mtc(j)
  if(elens_lThetaR2(ielens(ix))) then
    if(elens_I(ielens(ix)) < zero) then
      frrelens = frrelens*((rvv(j)+elens_beta_e(ielens(ix))*betrel)/(one+elens_beta_e(ielens(ix))*betrel))
    else
      frrelens = frrelens*((rvv(j)-elens_beta_e(ielens(ix))*betrel)/(one-elens_beta_e(ielens(ix))*betrel))
    end if
  endif
  yv1(j)=yv1(j)-(frrelens*xelens)/rrelens
  yv2(j)=yv2(j)-(frrelens*yelens)/rrelens
endif
! end include/kickelens.f90
