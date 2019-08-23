! start include/thcklin.f90
  puxve1=xv1(j)
  puzve1=yv1(j)
  puxve2=xv2(j)
  puzve2=yv2(j)
  sigmv(j)=(((((sigmv(j)+as(1,1,j,jx))+puxve1*((as(2,1,j,jx)+ as &
  (4,1,j,jx)*puzve1)+as(5,1,j,jx)*puxve1))+ puzve1*(as           &
  (3,1,j,jx)+as(6,1,j,jx)*puzve1))                               &
  +as(1,2,j,jx))+puxve2*(as(2,2,j,jx)+ as                        &
  (4,2,j,jx)*puzve2+as(5,2,j,jx)*puxve2))+ puzve2*(as            &
  (3,2,j,jx)+as(6,2,j,jx)*puzve2)
    xv1(j)=(al(1,1,j,jx)*puxve1+ al(2,1,j,jx)*puzve1)+          &
  real(idz1,fPrec)*al(5,1,j,jx)
    xv2(j)=(al(1,2,j,jx)*puxve2+ al(2,2,j,jx)*puzve2)+          &
  real(idz2,fPrec)*al(5,2,j,jx)
    yv1(j)=(al(3,1,j,jx)*puxve1+ al(4,1,j,jx)*puzve1)+          &
  real(idz1,fPrec)*al(6,1,j,jx)
    yv2(j)=(al(3,2,j,jx)*puxve2+ al(4,2,j,jx)*puzve2)+          &
  real(idz2,fPrec)*al(6,2,j,jx)
  ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
  ! last modified: 07-03-2018
  ! store old particle coordinates
  if (lbacktracking) call aperture_saveLastCoordinates(i,ix,kz(jx))
! end include/thcklin.f90
