!print *,"strackxxxxxx", strackx(i) !This i ks/2 
!print *,"strackzzzzz", strackz(i)  !ksi/2

! x'-> px; y'->py
!Transformation to some more reasonable coordinates

temptr(1)=xv(1,j)
temptr(2)=yv(1,j)/moidpsv(j)
temptr(3)=xv(2,j)
temptr(4)=yv(2,j)/moidpsv(j)

temptr(6)=(ejv(j)-e0)/(e0f*(e0f/e0))

!yv(1,j) = temptr(2)
!yv(2,j) = temptr(4)
!print *, "before", xv(1,j), yv(1,j), xv(2,j), yv(2,j), sigmv(j), moidpsv(j), j
onedp   = sqrt( one + two*temptr(6) + ((e0f/e0)**2)*(temptr(6)**2))
fppsig  = ( one + ((e0f/e0) **2)*temptr(6) ) / onedp

!     Set up C,S, q_temp,r_temp,Z
costh_temp = cos_mb(strackz(i)/onedp)
sinth_temp = sin_mb(strackz(i)/onedp)
q_temp = -strackz(i) * strackx(i) / onedp
r_temp = (fppsig / (onedp**2)) * (strackz(i) * strackx(i))
z_temp = (fppsig / (onedp**2))* (strackz(i))


!print *,"costh_temp", costh_temp !This i ks/2 
!print *,"sinth_temp", sinth_temp  !ksi/2
pxf  = temptr(2) + temptr(1)*q_temp
pyf  = temptr(4) +  temptr(3) *q_temp

!print *, "pxf", pxf , temptr(2)
!print *, "pyf", pyf , temptr(4)

!       r_tempipken formulae p.29 (3.37)
xv(1,j) =  (temptr(1)  * costh_temp  +  temptr(3)  * sinth_temp)
yv(1,j) =  (pxf * costh_temp  +  pyf * sinth_temp)*moidpsv(j)
xv(2,j) = (-temptr(1)  * sinth_temp  +  temptr(3)  * costh_temp)
yv(2,j) = (-pxf * sinth_temp  +  pyf * costh_temp)*moidpsv(j)



!print *, "after", xv(1,j), yv(1,j), xv(2,j), yv(2,j), sigmv(j), moidpsv(j), j
