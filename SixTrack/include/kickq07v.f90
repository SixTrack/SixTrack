#ifndef TILT
qu=(six*ekk)*cikve                                               !hr02
qv=(-six*ekk)*crkve                                              !hr02
#else
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=(six*ekk)*(tiltck*cikve-tiltsk*crkve)                         !hr02
qv=(-six*ekk)*(tiltck*crkve+tiltsk*cikve)                        !hr02
#endif
