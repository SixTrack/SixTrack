! start include/kicklxxv.f90
#ifndef TILT
dyy1=ekk*cikve
dyy2=ekk*crkve
#else
dyy1=ekk*(tiltc(k)*cikve-tilts(k)*crkve)
dyy2=ekk*(tiltc(k)*crkve+tilts(k)*cikve)
#endif
! end include/kicklxxv.f90
