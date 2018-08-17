! start include/kickq09h.f90
#ifndef TILT
qu=(eight*ekk)*crkve                                               !hr02
qv=(eight*ekk)*cikve                                               !hr02
#else
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=(eight*ekk)*(tiltck*crkve+tiltsk*cikve)                         !hr02
qv=(eight*ekk)*(tiltck*cikve-tiltsk*crkve)                         !hr02
#endif
! end include/kickq09h.f90
