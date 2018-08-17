!start kickvso1.f90
onedp   =  (one+dpsv(j))/mtc(j)
fppsig  = ( one + ((e0f/e0) **2)*temptr(6) ) / onedp
!
temptr(1)=xv(1,j)
temptr(2)=yv(1,j)
temptr(3)=xv(2,j)
temptr(4)=yv(2,j)

temptr(6)=(ejv(j)-e0)/(e0f*(e0f/e0))


!     Set up C,S, q_temp,r_temp,Z
costh_temp = cos_mb(strackz(i)/onedp)
sinth_temp = sin_mb(strackz(i)/onedp)


q_temp = -strackz(i) * strackx(i) / (onedp)

pxf  = temptr(2) + temptr(1)*q_temp / (onedp) 
pyf  = temptr(4) +  temptr(3) *q_temp / (onedp)


!       r_tempipken formulae p.29 (3.37)
xv(1,j) =  (temptr(1)  * costh_temp  +  temptr(3)  * sinth_temp)
yv(1,j) =  (pxf * costh_temp  +  pyf * sinth_temp)
xv(2,j) = (-temptr(1)  * sinth_temp  +  temptr(3)  * costh_temp)
yv(2,j) = (-pxf * sinth_temp  +  pyf * costh_temp)

!end kickvso1.f90
