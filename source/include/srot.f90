cos_t = cos_mb(temp_angle)
sin_t = -sin_mb(temp_angle)
do j=1,napx
temptr(1)=xv1(j)
temptr(2)=yv1(j)
temptr(3)=xv2(j)
temptr(4)=yv2(j)
     
xv1(j) = temptr(1)*cos_t - temptr(3)*sin_t
yv1(j) = temptr(2)*cos_t - temptr(4)*sin_t
xv2(j) = temptr(1)*sin_t + temptr(3)*cos_t
yv2(j) = temptr(2)*sin_t + temptr(4)*cos_t
enddo