! start include/kickq05h.f90
#ifndef TILT
qu=(four*ekk)*crkve                                              !hr02
qv=(four*ekk)*cikve                                              !hr02
#else
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=(four*ekk)*(tiltck*crkve+tiltsk*cikve)                        !hr02
qv=(four*ekk)*(tiltck*cikve-tiltsk*crkve)                        !hr02
#endif
! end include/kickq05h.f90
