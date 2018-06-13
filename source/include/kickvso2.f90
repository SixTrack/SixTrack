crkve=sigmv(j)-half*(((((((xv(1,j)**2+xv(2,j)**2)*strackx(i))*strackz(i))*rvv(j))*ejf0v(j))/ejfv(j))*ejf0v(j))/ejfv(j)
sigmv(j)=crkve
crkve=yv(1,j)-(((xv(1,j)*strackx(i))*strackz(i))*ejf0v(j))/ejfv(j)
cikve=yv(2,j)-(((xv(2,j)*strackx(i))*strackz(i))*ejf0v(j))/ejfv(j)
sigmv(j)=sigmv(j)+((((((xv(1,j)*cikve-xv(2,j)*crkve)*strackz(i))*rvv(j))*ejf0v(j))/ejfv(j))*ejf0v(j))/ejfv(j)
