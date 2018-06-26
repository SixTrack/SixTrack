yv(1,j)=yv(1,j)-xv(2,j)*strackx(i)
yv(2,j)=yv(2,j)+xv(1,j)*strackx(i)
!
! TODO: Check if ejf0v should be e0f?? or oidpsv=ejf0v(j)/ejfv(j)=1/(1+delta)
!
crkve=yv(1,j)-(((xv(1,j)*strackx(i))*strackz(i))*ejf0v(j))/ejfv(j) !hr02
cikve=yv(2,j)-(((xv(2,j)*strackx(i))*strackz(i))*ejf0v(j))/ejfv(j) !hr02
yv(1,j)=crkve*cos_mb((strackz(i)*ejf0v(j))/ejfv(j))+cikve*sin_mb((strackz(i)*ejf0v(j))/ejfv(j))
yv(2,j)=cikve*cos_mb((strackz(i)*ejf0v(j))/ejfv(j))-crkve*sin_mb((strackz(i)*ejf0v(j))/ejfv(j))
crkve=xv(1,j)*cos_mb((strackz(i)*ejf0v(j))/ejfv(j))+xv(2,j)*sin_mb((strackz(i)*ejf0v(j))/ejfv(j))
cikve=xv(2,j)*cos_mb((strackz(i)*ejf0v(j))/ejfv(j))-xv(1,j)*sin_mb((strackz(i)*ejf0v(j))/ejfv(j))
xv(1,j)=crkve
xv(2,j)=cikve
yv(1,j)=yv(1,j)+xv(2,j)*strackx(i)
yv(2,j)=yv(2,j)-xv(1,j)*strackx(i)
