! start include/kickldpe.f90
#ifndef TILT
dyy1=(ed(IX)*crkve)/(one+dpp)
dyy2=(ek(IX)*cikve)/(one+dpp)
#else
dyy1=((ed(IX)*tiltc(k))*crkve-(ek(IX)*tilts(k))*cikve)/(one+dpp)
dyy2=((ek(IX)*tiltc(k))*cikve+(ed(IX)*tilts(k))*crkve)/(one+dpp)
#endif
! end include/kickldpe.f90
