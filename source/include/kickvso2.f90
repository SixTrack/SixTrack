! start include/kickvso2.f90
pxf  = temptr(2)*onedp + temptr(1)*q_temp
pyf  = temptr(4)*onedp +  temptr(3) *q_temp
r_temp = (fppsig / (onedp**2)) * (strackz(i) * strackx(i))
z_temp = (fppsig / (onedp**2))* (strackz(i))
sigf = sigmv(j) - (half*(temptr(1)**2 +  temptr(3) **2)*r_temp)*c1m3
sigmv(j) =  (sigf + (temptr(1)*pyf - temptr(3)*pxf)*(z_temp*c1m3))
! end include/kickvso2.f90
