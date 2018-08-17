! start include/kickudpe.f90
#ifndef TILT
dyy1=ed(IX)*crkve
dyy2=ek(IX)*cikve
y(1,1)=y(1,1)+dyy1/(one+dpp)
y(1,2)=y(1,2)+dyy2/(one+dpp)
#else
dyy1=(ed(IX)*crkve)/(one+dpp)                                    !hr02
dyy2=(ek(IX)*cikve)/(one+dpp)                                    !hr02
y(1,1)=(y(1,1)+tiltc(k)*dyy1)-tilts(k)*dyy2                      !hr02
y(1,2)=(y(1,2)+tiltc(k)*dyy2)+tilts(k)*dyy1                      !hr02
#endif
! end include/kickudpe.f90
