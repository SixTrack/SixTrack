! start include/kickq06v.f90
#ifndef TILT
qu=(five*ekk)*cikve                                               !hr02
qv=(-five*ekk)*crkve                                              !hr02
#else
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=(five*ekk)*(tiltck*cikve-tiltsk*crkve)                         !hr02
qv=(-five*ekk)*(tiltck*crkve+tiltsk*cikve)                        !hr02
#endif
! end include/kickq06v.f90
