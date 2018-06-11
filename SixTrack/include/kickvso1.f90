print *,"strackxxxxxx", strackx(i) !This i ks/2 
print *,"strackzzzzz", strackz(i)  !ksi/2

! x'-> px; y'->py
!Transformation to some more reasonable coordinates
temptr(1)=xv(1,j)
temptr(2)=yv(1,j)/moidpsv(j)
temptr(3)=xv(2,j)
temptr(4)=yv(2,j)/moidpsv(j)
!temptr(5)=sigmv(j)
temptr(6)=((mtc(j)*ejv(j)-e0)/e0f)*c1e3*(e0/e0f) 




!     We do not use a constant deltap!!!!! WE use full 6D formulae!
onedp   = sqrt( one + two*temptr(6) + ((e0f/e0)**2)*(temptr(6)**2) )
fpsig   = onedp - one
fppsig  = ( one + ((e0f/e0) **2)*temptr(6) ) / onedp

!     Set up C,S, q_temp,r_temp,Z
costh_temp = cos_mb(strackz(i)/onedp)
sinth_temp = sin_mb(strackz(i)/onedp)
q_temp = -strackz(i) * strackx(i) / onedp
r_temp = fppsig / (onedp**2) * strackz(i) * strackx(i)
z_temp = fppsig / (onedp**2) * strackz(i)

pxf  = temptr(2) + temptr(1)*q_temp
pyf  = temptr(4) +  temptr(3) *q_temp
sigf = sigmv(j)*(e0f/e0) - half*(temptr(1)**2 +  temptr(3) **2)*r_temp

!       r_tempipken formulae p.29 (3.37)
xv(1,j) =  temptr(1)  * costh_temp  +  temptr(3)  * sinth_temp
yv(1,j) =  pxf * costh_temp  +  pyf * sinth_temp
xv(2,j) = -temptr(1)  * sinth_temp  +  temptr(3)  * costh_temp
yv(2,j) = -pxf * sinth_temp  +  pyf * costh_temp
sigmv(j) =  (sigf + (temptr(1)*pyf - temptr(3)*pxf)*z_temp)


yv(j,1) = yv(j,1)*mtc(j)/(one+dpsv(j))
yv(j,2) = yv(j,2)*mtc(j)/(one+dpsv(j))








!yv(1,j)=yv(1,j)-xv(2,j)*strackx(i)
!yv(2,j)=yv(2,j)+xv(1,j)*strackx(i)
!
! TODO: Check if ejf0v should be e0f?? or oidpsv=ejf0v(j)/ejfv(j)=1/(1+delta)
!crkve, cikve
!crkve=yv(1,j)-(((xv(1,j)*strackx(i))*strackz(i))*ejf0v(j))/ejfv(j) !hr02
!cikve=yv(2,j)-(((xv(2,j)*strackx(i))*strackz(i))*ejf0v(j))/ejfv(j) !hr02
!yv(1,j)=crkve*cos_mb((strackz(i)*ejf0v(j))/ejfv(j))+cikve*sin_mb((strackz(i)*ejf0v(j))/ejfv(j))
!yv(2,j)=cikve*cos_mb((strackz(i)*ejf0v(j))/ejfv(j))-crkve*sin_mb((strackz(i)*ejf0v(j))/ejfv(j))
!crkve=xv(1,j)*cos_mb((strackz(i)*ejf0v(j))/ejfv(j))+xv(2,j)*sin_mb((strackz(i)*ejf0v(j))/ejfv(j))
!cikve=xv(2,j)*cos_mb((strackz(i)*ejf0v(j))/ejfv(j))-xv(1,j)*sin_mb((strackz(i)*ejf0v(j))/ejfv(j))
!xv(1,j)=crkve
!xv(2,j)=cikve
!yv(1,j)=yv(1,j)+xv(2,j)*strackx(i)
!yv(2,j)=yv(2,j)-xv(1,j)*strackx(i)
