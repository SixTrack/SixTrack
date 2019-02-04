! start include/kickaso1.f90
crkve=y(1,1)-((x(1,1)*ed(IX))*ek(IX))/(one+dpp)              !hr02
cikve=y(1,2)-((x(1,2)*ed(IX))*ek(IX))/(one+dpp)              !hr02
dyy1=(crkve*cos_mb(ek(IX))/(one+dpp))+(cikve*sin_mb(ek(IX))/(one+dpp))-y(1,1)  !hr02
dyy2=cikve*cos_mb(ek(IX)/(one+dpp))-crkve*sin_mb(ek(IX)/(one+dpp))-y(1,2)      !hr02
mpe=20
qu=ed(IX)
qv=ek(IX)
! end include/kickaso1.f90
