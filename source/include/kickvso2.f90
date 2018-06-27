print *, sigmv(j), j
sigf = sigmv(j) - (half*(temptr(1)**2 +  temptr(3) **2)*r_temp)*c1m3
sigmv(j) =  (sigf + (temptr(1)*pyf - temptr(3)*pxf)*(z_temp*c1m3))
!print *, xv(1,j), yv(1,j), xv(2,j), yv(2,j), sigmv(j)
print *, sigmv(j), j