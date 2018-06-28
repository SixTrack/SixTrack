r_temp = (fppsig / (onedp**2)) * (strackz(i) * strackx(i))
z_temp = (fppsig / (onedp**2))* (strackz(i))
sigf = sigmv(j) - (half*(temptr(1)**2 +  temptr(3) **2)*r_temp)*c1m3
sigmv(j) =  (sigf + (temptr(1)*pyf - temptr(3)*pxf)*(z_temp*c1m3))
