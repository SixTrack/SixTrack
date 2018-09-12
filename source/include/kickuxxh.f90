! start include/kickuxxh.f90
#ifndef TILT
dyy1=ekk*crkve
dyy2=(-one*ekk)*cikve                                            !hr08
y(1,1)=y(1,1)+dyy1
y(1,2)=y(1,2)+dyy2
#else
dyy1=ekk*crkve
dyy2=(-one*ekk)*cikve                                            !hr08
y(1,1)=(y(1,1)+tiltc(k)*dyy1)-tilts(k)*dyy2                      !hr02
y(1,2)=(y(1,2)+tiltc(k)*dyy2)+tilts(k)*dyy1                      !hr02
#endif
! end include/kickuxxh.f90
