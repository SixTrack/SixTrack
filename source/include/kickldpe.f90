! start include/kickldpe.f90
#ifndef TILT
dyy1=(ed(IX)*crkve)/(one+dpp)                                    !hr02
dyy2=(ek(IX)*cikve)/(one+dpp)                                    !hr02
#else
dyy1=((ed(IX)*tiltc(k))*crkve-(ek(IX)*tilts(k))*cikve)/(one+dpp) !hr02
dyy2=((ek(IX)*tiltc(k))*cikve+(ed(IX)*tilts(k))*crkve)/(one+dpp) !hr02
#endif
! end include/kickldpe.f90
