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
    case default
      write(lout,"(a,i0)") "ELENS> ERROR in kickelens: elens_type=",elens_type(ielens(ix))," not recognized. "
      write(lout,"(a)")    "ELENS>       Possible values for type are: 1 and 2"
      call prror(-1)
    end select
  else ! r1 < r2 <= rrelens
    frrelens = one
  endif
  frrelens = elens_r2(ielens(ix))/rrelens * frrelens
  yv1(j)=yv1(j)-elens_theta_r2(ielens(ix))*frrelens*xelens/rrelens * oidpsv(j)
  yv2(j)=yv2(j)-elens_theta_r2(ielens(ix))*frrelens*yelens/rrelens * oidpsv(j)
endif
! end include/kickelens.f90
