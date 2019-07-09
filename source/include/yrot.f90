! start include/yrot.f90
cos_t = cos_mb(temp_angle)
sin_t = sin_mb(temp_angle)
tan_t = tan_mb(temp_angle)

do j=1,napx
  temptr(1)=c1m3*xv1(j)
  temptr(2)=(c1m3*yv1(j))/moidpsv(j)
  temptr(3)=c1m3*xv2(j)
  temptr(4)=(c1m3*yv2(j))/moidpsv(j)
  temptr(5)=(c1m3*sigmv(j))/(e0f/e0)
  temptr(6)=c1m3*((mtc(j)*ejv(j)-e0)/e0f)

  z_temp = sqrt((one + dpsv(j))**2 - temptr(2)**2 - temptr(4)**2)
  pttemp = 1 - (tan_t*temptr(2))/z_temp

  xv1(j) = (c1e3*temptr(1))/(cos_t*pttemp)
  xv2(j) = xv2(j) + c1e3*(tan_t*(temptr(4)*temptr(1))/(z_temp*pttemp))
  yv1(j) = c1e3*(cos_t*temptr(2) + sin_t*z_temp)*moidpsv(j)
  sigmv(j) = sigmv(j) - c1e3*((tan_t*temptr(1)*(one/(e0f/e0)+temptr(6))/(z_temp*pttemp))*(e0f/e0))
enddo
! end include/yrot.f90
