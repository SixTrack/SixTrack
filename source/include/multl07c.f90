! start include/multl07c.f90
l1=l-1
qu=qu+real(l1,fPrec)*(bb(l)*crkve+aa(l)*cikve)
qv=qv+real(l1,fPrec)*(bb(l)*cikve-aa(l)*crkve)
crkveuk=crkve*xl-cikve*zl
cikve=crkve*zl+cikve*xl
crkve=crkveuk
dyy1=(dyy1+bb(l)*crkve)+aa(l)*cikve
dyy2=(dyy2-bb(l)*cikve)+aa(l)*crkve
! end include/multl07c.f90
