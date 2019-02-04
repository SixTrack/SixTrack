! start include/mul4v06.f90
crkveuk=crkve*xlvj-cikve*zlvj
cikve=crkve*zlvj+cikve*xlvj
crkve=crkveuk
yv1j=(yv1j+bbiv(k,i)*crkve)+aaiv(k,i)*cikve        !hr03
yv2j=(yv2j-bbiv(k,i)*cikve)+aaiv(k,i)*crkve        !hr03
! end include/mul4v06.f90
