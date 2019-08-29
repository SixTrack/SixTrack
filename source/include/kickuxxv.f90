! start include/kickuxxv.f90
#ifndef TILT
dyy1=ekk*cikve
dyy2=ekk*crkve
y(1,1)=y(1,1)+dyy1
y(1,2)=y(1,2)+dyy2
#else
dyy1=ekk*cikve
dyy2=ekk*crkve
y(1,1)=(y(1,1)+tiltc(k)*dyy1)-tilts(k)*dyy2
y(1,2)=(y(1,2)+tiltc(k)*dyy2)+tilts(k)*dyy1
#endif
! end include/kickuxxv.f90
