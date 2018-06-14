!!print *,"strackxxxxxx", strackx(i) !This i ks/2 
!!print *,"strackzzzzz", strackz(i)  !ksi/2

! x'-> px; y'->py
!Transformation to some more reasonable coordinates
temptr(1)=xv(1,j)*c1m3
temptr(2)=(yv(1,j)*c1m3)/moidpsv(j)
temptr(3)=xv(2,j)*c1m3
temptr(4)=(yv(2,j)*c1m3)/moidpsv(j)

temptr(6)=(ejv(j)-e0)/(e0f*(e0f/e0))


onedp   = sqrt( one + two*temptr(6) + ((e0f/e0)**2)*(temptr(6)**2))
fppsig  = ( one + ((e0f/e0) **2)*temptr(6) ) / onedp

!     Set up C,S, q_temp,r_temp,Z
costh_temp = cos_mb(strackz(i)/onedp)
sinth_temp = sin_mb(strackz(i)/onedp)
q_temp = -strackz(i) * strackx(i) / onedp
r_temp = fppsig / (onedp**2) * strackz(i) * strackx(i)
z_temp = fppsig / (onedp**2) * strackz(i)

pxf  = temptr(2) + temptr(1)*q_temp
pyf  = temptr(4) +  temptr(3) *q_temp
sigf = sigmv(j)*c1m3 - half*(temptr(1)**2 +  temptr(3) **2)*r_temp


!       r_tempipken formulae p.29 (3.37)
xv(1,j) =  (temptr(1)  * costh_temp  +  temptr(3)  * sinth_temp)*c1e3
yv(1,j) =  (pxf * costh_temp  +  pyf * sinth_temp)*c1e3
xv(2,j) = (-temptr(1)  * sinth_temp  +  temptr(3)  * costh_temp)*c1e3
yv(2,j) = (-pxf * sinth_temp  +  pyf * costh_temp)*c1e3
sigmv(j) =  (sigf + (temptr(1)*pyf - temptr(3)*pxf)*z_temp)*c1e3

yv(j,1) = yv(j,1)*moidpsv(j)
yv(j,2) = yv(j,2)*moidpsv(j)

