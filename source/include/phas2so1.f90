! start include/phas2so1.f90
crkve=t(i,2)-(t(i,1)*qu)*qv                                  !hr02
cikve=t(i,4)-(t(i,3)*qu)*qv                                  !hr02
t(i,2)=crkve*cos_mb(qv)+cikve*sin_mb(qv)                     !hr02
t(i,4)=cikve*cos_mb(qv)-crkve*sin_mb(qv)                     !hr02
crkve=t(i,1)*cos_mb(qv)+t(i,3)*sin_mb(qv)                    !hr02
cikve=t(i,3)*cos_mb(qv)-t(i,1)*sin_mb(qv)                    !hr02
t(i,1)=crkve
t(i,3)=cikve
! end include/phas2so1.f90
