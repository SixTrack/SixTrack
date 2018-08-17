!start mul4v06.f90
crkveuk=crkve*xlvj-cikve*zlvj
cikve=crkve*zlvj+cikve*xlvj
crkve=crkveuk
yv1j=(yv1j+bbiv(k,1,i)*crkve)+aaiv(k,1,i)*cikve        !hr03
yv2j=(yv2j-bbiv(k,1,i)*cikve)+aaiv(k,1,i)*crkve        !hr03
!end mul4v06.f90
