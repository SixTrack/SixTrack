! start include/kickq08h.f90
#ifndef TILT
qu=(seven*ekk)*crkve                                               !hr02
qv=(seven*ekk)*cikve                                               !hr02
#else
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=(seven*ekk)*(tiltck*crkve+tiltsk*cikve)                         !hr02
qv=(seven*ekk)*(tiltck*cikve-tiltsk*crkve)                         !hr02
#endif
! end include/kickq08h.f90
