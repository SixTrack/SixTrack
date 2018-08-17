! start include/kickq07h.f90
#ifndef TILT
qu=(six*ekk)*crkve                                               !hr02
qv=(six*ekk)*cikve                                               !hr02
#else
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=(six*ekk)*(tiltck*crkve+tiltsk*cikve)                         !hr02
qv=(six*ekk)*(tiltck*cikve-tiltsk*crkve)                         !hr02
#endif
! end include/kickq07h.f90
