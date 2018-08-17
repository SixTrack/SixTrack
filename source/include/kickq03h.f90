! start include/kickq03h.f90
#ifndef TILT
qu=(ekk*two)*crkve                                               !hr02
qv=(ekk*two)*cikve                                               !hr02
#else
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=(ekk*two)*(tiltck*crkve+tiltsk*cikve)                         !hr02
qv=(ekk*two)*(tiltck*cikve-tiltsk*crkve)                         !hr02
#endif
! end include/kickq03h.f90
